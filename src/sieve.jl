const wheel = (1, 7, 11, 13, 17, 19, 23, 29)              # residues coprime to 30, ascending

# INV_MOD_30[r + 1] = inverse of r modulo 30 (for r coprime to 30)
const INV_MOD_30 = let table = fill(0, 30)
    for r in 1:29
        gcd(r, 30) == 1 && (table[r + 1] = invmod(r, 30))
    end
    Tuple(table)
end

# Residue → lane index, the inverse of `wheel`: WHEEL_IDX[r + 1] = i with wheel[i] == r.
const WHEEL_IDX = let table = fill(0, 30)
    for (i, r) in enumerate(wheel)
        table[r + 1] = i
    end
    Tuple(table)
end

# Per (residue-class pri, lane li) stride constants, so the large-prime loop places a multiple
# without dividing by 30. For a prime p ≡ wheel[pri] (mod 30) with pq = p÷30, its smallest
# lane-li multiple sits in block  pq·STRIDE_KW[pri][li] + STRIDE_C[pri][li].  (STRIDE_KW[pri][li]
# is the least k with p·k ≡ wheel[li] (mod 30); STRIDE_C folds in the ÷30.)
const STRIDE_KW = ntuple(pri -> ntuple(li -> mod(wheel[li] * INV_MOD_30[wheel[pri] + 1], 30), Val(8)), Val(8))
const STRIDE_C  = ntuple(pri -> ntuple(li -> (wheel[pri] * STRIDE_KW[pri][li] - wheel[li]) ÷ 30, Val(8)), Val(8))
# p² depends only on p's residue class mod 30, so p²÷30 = 30·pq² + 2·pq·wheel[pri] + STRIDE_PSQ_DIV[pri]
# and p²%30 = STRIDE_PSQ_RES[pri]; together cld(p²-w, 30) = p²÷30 + (STRIDE_PSQ_RES[pri] > w).
const STRIDE_PSQ_DIV = ntuple(pri -> wheel[pri]^2 ÷ 30, Val(8))
const STRIDE_PSQ_RES = ntuple(pri -> wheel[pri]^2 % 30, Val(8))

# build_combtab(p)[lane][(abs_chunk % p) + 1] is the mask that clears, in residue `lane`, the
# multiples of p. The pattern has period p chunks, since 30·64·p ≡ 0 (mod p).
function build_combtab(p::Int)
    ntuple(Val(8)) do lane
        res = wheel[lane]
        mask = fill(~UInt64(0), p)
        for chunk in 0:(p - 1), bit in 0:63
            if (30 * (64 * chunk + bit) + res) % p == 0
                mask[chunk + 1] &= ~(UInt64(1) << bit)
            end
        end
        mask
    end
end

# Small primes are cheaper to sieve word at a time rather than prime at a time.
const COMB_THRESH = 256
const COMB_PRIMES = Int[p for p in 7:(COMB_THRESH - 1) if all(p % d != 0 for d in 2:isqrt(p))]
const COMB_TABLE   = NTuple{8,Vector{UInt64}}[build_combtab(p) for p in COMB_PRIMES]

# AND the lane against its period-p comb, realigned to the window's absolute phase so the
# inner loop is a contiguous (vectorizable) slice-AND — no `%p` per chunk.
@inline function comb_clear!(chunks::Vector{UInt64}, comb::Vector{UInt64}, p::Int, first_chunk::Int, nchunks::Int)
    chunk = 1
    phase = first_chunk % p
    @inbounds while chunk ≤ nchunks
        run = min(p - phase, nchunks - chunk + 1)
        @simd for j in 0:(run - 1)
            chunks[chunk + j] &= comb[phase + 1 + j]
        end
        chunk += run
        phase = 0
    end
end

# Clear every `stride`-th bit of the lane, from bit `startbit` up to `nbits`.
@inline function stride_clear!(chunks::Vector{UInt64}, startbit::Int, stride::Int, nbits::Int)
    bit = startbit
    @inbounds while bit < nbits
        chunks[(bit >>> 6) + 1] &= ~(UInt64(1) << (bit & 63))
        bit += stride
    end
end

"""
Segmented mod-30 lane sieve: one BitVector per residue coprime to 30,
rather than 8 bits packed per 30-block.
Within a lane the multiples of a prime p form a constant-stride-p progression.
As an optimization:
    p < COMB_THRESH : AND the lane with a periodic stride-p comb
    p ≥ COMB_THRESH : scalar stride-p bit-clear
Lanes are indexed in absolute block space, so a comb mask is a fixed function of
(abs_chunk % p) and can be built once per prime.
Output byte-transposes the 8 lane words into packed block-major words and scans those.

A `SegmentedSieve` is an iterator of windows: iterating fills the lanes for the next window
and yields the sieve itself. `each_lane_prime`, `primes`, `primesmask`, and `eachprime`
consume those windows.

2, 3, 5 are the caller's responsibility (the wheel skips them).
"""
mutable struct SegmentedSieve
    lo::Int
    hi::Int
    seg_blocks::Int          # lane size (blocks per window)
    # Large base primes (≥ COMB_THRESH) as (p÷30, residue-class index); p and p² reconstruct from
    # these, so the raw prime values are never stored. Small primes are the COMB_PRIMES const.
    stride_primes::Vector{Tuple{Int32,Int8}}
    lanes::Vector{BitVector} # the current window: 8 lanes, ≤ seg_blocks bits each
    win_block::Int           # current window: first absolute block
    nblocks::Int             # current window: number of blocks
end

# Optimal number of blocks per window. Low cap is 32 KiB/lane allowing for sieving within L1.
# As hi grows, we start sieving by bigger primes so we grow the window, but cap at MiB/lane to stay within L2.
_seg_blocks(hi::Integer) = clamp(isqrt(hi) << 3, 1 << 18, 1 << 24)

# Chunk-aligned first block at or below `lo`
_first_block(lo::Integer) = (fld(lo, 30) >>> 6) << 6

function SegmentedSieve(lo::Int, hi::Int; seg_blocks::Int = _seg_blocks(hi))
    lo = max(7, lo)
    # Decompose the stride primes (≥ COMB_THRESH; smaller ones are the comb) into (p÷30, residue
    # class) once. Stream them with eachprime so we never allocate the full base-prime vector.
    stride_primes = Tuple{Int32,Int8}[]
    sq = isqrt(hi)
    if sq ≥ COMB_THRESH
        for p in eachprime(COMB_THRESH, sq)
            pq, pr = divrem(p, 30)
            push!(stride_primes, (pq % Int32, WHEEL_IDX[pr + 1] % Int8))
        end
    end
    num_blocks = min(seg_blocks, fld(hi, 30) - _first_block(lo) + 1)
    lanes = [BitVector(undef, num_blocks) for _ in 1:8]
    return SegmentedSieve(lo, hi, seg_blocks, stride_primes, lanes, 0, 0)
end

Base.IteratorSize(::Type{SegmentedSieve}) = Base.SizeUnknown()
Base.eltype(::Type{SegmentedSieve}) = SegmentedSieve

# Fill the next window into the lanes and yield the sieve itself (with `win_block`/`nblocks` set
# to the current window). Boundary bits (< lo, > hi, padding) are masked, so consumers need
# no per-prime range checks.
function Base.iterate(ss::SegmentedSieve, win_block::Int = _first_block(ss.lo))
    last_block = fld(ss.hi, 30)
    win_block > last_block && return nothing
    lanes = ss.lanes
    nblocks = min(ss.seg_blocks, last_block - win_block + 1)
    nchunks = cld(nblocks, 64)
    first_chunk = win_block ÷ 64                       # absolute chunk of local chunk 1
    win_last_block = win_block + nblocks - 1
    win_max = win_last_block == last_block ? ss.hi : 30 * win_last_block + 29
    max_prime = isqrt(win_max)

    # Loop-invariant boundary conditions (applied per lane below to make consumers range-check-free).
    lo_block, lo_res = divrem(ss.lo, 30)               # values < lo occupy blocks ≤ lo_block
    tail_local = last_block - win_block                # local block of hi (if hi falls in this window)
    hi_res = ss.hi - 30 * last_block                   # hi % 30, overflow-safe threshold
    tail_bits = nblocks - (nchunks - 1) * 64           # valid bits in the final chunk (64 ⇒ no padding)

    # One lane fully before the next, so it stays cache-hot across all its primes: fill it, comb
    # the small primes (restoring each prime's own bit), stride the large ones, then mask boundaries.
    @inbounds for lane in 1:8
        chunks = lanes[lane].chunks
        fill!(chunks, ~UInt64(0))
        w = wheel[lane]
        for ci in eachindex(COMB_PRIMES)               # small primes: periodic comb AND
            p = COMB_PRIMES[ci]
            p > max_prime && break
            comb_clear!(chunks, COMB_TABLE[ci][lane], p, first_chunk, nchunks)
            if win_block == 0 && p % 30 == w           # comb cleared p's own bit (p·1) and this is p's lane
                b = fld(p, 30)
                chunks[b >>> 6 + 1] |= UInt64(1) << (b & 63)
            end
        end
        for (pq32, pri) in ss.stride_primes            # large primes: scalar stride clear
            # Reconstruct p and p²÷30 from (pq, residue class) — no division. res_block is the phase:
            # lane-w multiples sit at blocks ≡ res_block (mod p). The only division left is mod p.
            pq = Int(pq32)
            r = wheel[pri]
            p = 30 * pq + r
            p > max_prime && break                     # primes are ascending; rest are too big for this window
            res_block = pq * STRIDE_KW[pri][lane] + STRIDE_C[pri][lane]
            min_block = 30 * pq * pq + 2 * pq * r + STRIDE_PSQ_DIV[pri] + (STRIDE_PSQ_RES[pri] > w)
            from_block = max(min_block, win_block)      # min_block = cld(p²-w, 30): first block ≥ p²
            start_block = from_block + mod(res_block - from_block, p)
            local_start = start_block - win_block
            local_start < nblocks && stride_clear!(chunks, local_start, p, nblocks)
        end
        if win_block ≤ lo_block                        # first window: holds values < lo (incl. value 1)
            for block in win_block:(lo_block - 1)      # whole blocks below lo
                loc = block - win_block
                chunks[loc >>> 6 + 1] &= ~(UInt64(1) << (loc & 63))
            end
            if w < lo_res                              # partial block holding lo
                loc = lo_block - win_block
                chunks[loc >>> 6 + 1] &= ~(UInt64(1) << (loc & 63))
            end
        end
        if tail_local < nblocks && w > hi_res          # clear values > hi in the tail block
            chunks[tail_local >>> 6 + 1] &= ~(UInt64(1) << (tail_local & 63))
        end
        if tail_bits < 64                              # zero the final chunk's padding
            chunks[nchunks] &= (UInt64(1) << tail_bits) - 1
        end
    end

    ss.win_block = win_block
    ss.nblocks = nblocks
    return (ss, win_block + ss.seg_blocks)
end

# Exchange the bit-blocks of `a` and `b` selected by `mask` / `shift`
@inline function delta_swap(a::UInt64, b::UInt64, shift::Int, mask::UInt64)
    t = (a ⊻ (b << shift)) & mask
    return a ⊻ t, b ⊻ (t >> shift)
end

# transpose 8 Int64s viewed as a matrix of 8x8 bytes
# A 3-stage delta-swap network at byte distances 4, 2, 1.
@inline function transpose_lanes(w::NTuple{8,UInt64})
    a1, a2, a3, a4, a5, a6, a7, a8 = w
    a1, a5 = delta_swap(a1, a5, 32, 0xFFFFFFFF00000000)
    a2, a6 = delta_swap(a2, a6, 32, 0xFFFFFFFF00000000)
    a3, a7 = delta_swap(a3, a7, 32, 0xFFFFFFFF00000000)
    a4, a8 = delta_swap(a4, a8, 32, 0xFFFFFFFF00000000)
    a1, a3 = delta_swap(a1, a3, 16, 0xFFFF0000FFFF0000)
    a2, a4 = delta_swap(a2, a4, 16, 0xFFFF0000FFFF0000)
    a5, a7 = delta_swap(a5, a7, 16, 0xFFFF0000FFFF0000)
    a6, a8 = delta_swap(a6, a8, 16, 0xFFFF0000FFFF0000)
    a1, a2 = delta_swap(a1, a2, 8,  0xFF00FF00FF00FF00)
    a3, a4 = delta_swap(a3, a4, 8,  0xFF00FF00FF00FF00)
    a5, a6 = delta_swap(a5, a6, 8,  0xFF00FF00FF00FF00)
    a7, a8 = delta_swap(a7, a8, 8,  0xFF00FF00FF00FF00)
    return (a1, a2, a3, a4, a5, a6, a7, a8)
end

# transpose of Int64 viewed as a matrix of 8x8 bits: bit (8i+j) ↔ bit (8j+i).
@inline function transpose8(x::UInt64)
    t = (x ⊻ (x >> 7))  & 0x00AA00AA00AA00AA
    x ⊻= t ⊻ (t << 7)
    t = (x ⊻ (x >> 14)) & 0x0000CCCC0000CCCC
    x ⊻= t ⊻ (t << 14)
    t = (x ⊻ (x >> 28)) & 0x00000000F0F0F0F0
    x ⊻= t ⊻ (t << 28)
    return x
end

# Call f(prime) for each prime in the filled window, in ascending order.
function each_lane_prime(f, ss::SegmentedSieve)
    win_block = ss.win_block
    chunks = ntuple(i -> ss.lanes[i].chunks, Val(8))            # the 8 lanes' word arrays
    nchunks = cld(ss.nblocks, 64)
    @inbounds for chunk in 1:nchunks
        w = ntuple(i -> chunks[i][chunk], Val(8))               # one word from each lane
        reduce(|, w) == 0 && continue                           # no primes in these 64 blocks
        packed = transpose_lanes(w)                             # group g now holds 8 blocks × 8 residues
        chunk_base = 30 * (win_block + (chunk - 1) * 64)
        for g in 1:8                                            # each group is 8 consecutive blocks
            word = transpose8(packed[g])                        # within-word: bit 8·block + (lane-1)
            group_base = chunk_base + 30 * 8 * (g - 1)
            while word != 0
                pos = trailing_zeros(word)
                f(group_base + 30 * (pos >> 3) + wheel[(pos & 7) + 1])
                word &= word - 1
            end
        end
    end
end

"""
    primesmask([lo,] hi)

Returns a prime sieve, as a `BitArray`, of the positive integers (from `lo`, if specified)
up to `hi`. Useful when working with either primes or composite numbers.
"""
function primesmask(lo::Int, hi::Int)
    0 < lo ≤ hi || throw(ArgumentError("The condition 0 < lo ≤ hi must be met."))
    sieve = falses(hi - lo + 1)
    lo ≤ 2 ≤ hi && (sieve[3 - lo] = true)
    lo ≤ 3 ≤ hi && (sieve[4 - lo] = true)
    lo ≤ 5 ≤ hi && (sieve[6 - lo] = true)
    hi < 7 && return sieve
    offset = lo - 1
    for window in SegmentedSieve(max(7, lo), hi)
        each_lane_prime(window) do p
            sieve[p - offset] = true
        end
    end
    return sieve
end
primesmask(lo::Integer, hi::Integer) = lo ≤ hi ≤ typemax(Int) ? primesmask(Int(lo), Int(hi)) :
    throw(ArgumentError("Both endpoints of the interval to sieve must be ≤ $(typemax(Int)), got $lo and $hi."))

primesmask(limit::Int) = primesmask(1, limit)
primesmask(n::Integer) = n ≤ typemax(Int) ? primesmask(Int(n)) :
    throw(ArgumentError("Requested number of primes must be ≤ $(typemax(Int)), got $n."))

# Lazy prime stream over [lo, hi], backed by a SegmentedSieve.
# Buffers one window's primes at a time (refilling on exhaustion); 2, 3, 5 are yielded first.
# All iteration state lives in the (mutable) struct, so the iteration `state` is a dummy.
mutable struct EachPrime
    lo::Int
    hi::Int
    sieve::Union{SegmentedSieve,Nothing}  # nothing when hi < 7
    next_block::Int                       # first block of the next window to fetch
    buffer::Vector{Int}                   # the current window's primes
    pos::Int                              # next index into buffer
    phase::Int                            # 0,1,2 ⇒ emit 2,3,5; ≥3 ⇒ window mode
end

"""
    eachprime([lo,] hi)

Lazily iterate the primes in `[lo, hi]` (from 2 if `lo` is omitted) in increasing order.
"""
function eachprime(lo::Integer, hi::Integer)
    lo ≤ hi || throw(ArgumentError("The condition lo ≤ hi must be met."))
    lo, hi = Int(lo), Int(hi)
    sieve = hi < 7 ? nothing : SegmentedSieve(max(7, lo), hi)
    next_block = sieve === nothing ? 0 : _first_block(sieve.lo)
    return EachPrime(lo, hi, sieve, next_block, Int[], 1, 0)
end
eachprime(hi::Integer) = eachprime(1, hi)

Base.IteratorSize(::Type{EachPrime}) = Base.SizeUnknown()
Base.eltype(::Type{EachPrime}) = Int

function Base.iterate(ep::EachPrime, ::Any = nothing)
    while ep.phase < 3                          # 2, 3, 5 precede every wheel prime
        p = (2, 3, 5)[ep.phase + 1]
        ep.phase += 1
        ep.lo ≤ p ≤ ep.hi && return (p, nothing)
    end
    while ep.pos > length(ep.buffer)            # refill from the next window(s)
        ep.sieve === nothing && return nothing
        next = iterate(ep.sieve, ep.next_block)
        next === nothing && return nothing
        window, ep.next_block = next            # window === ep.sieve (it yields itself)
        empty!(ep.buffer)
        each_lane_prime(v -> push!(ep.buffer, v), window)
        ep.pos = 1
    end
    p = ep.buffer[ep.pos]
    ep.pos += 1
    return (p, nothing)
end
