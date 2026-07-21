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
@inline function comb_clear!(chunks::Vector{UInt64}, comb::Vector{UInt64}, p::Int, phase::Int, nchunks::Int)
    chunk = 1
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
mutable struct SegmentedSieve{T<:Integer}
    win_block::T             # current window: first absolute block (block = value ÷ 30)
    lo::T
    hi::T
    seg_chunks::Int          # lane size (64-block chunks per window), so windows stay chunk-aligned
    # Large base primes (≥ COMB_THRESH) packed as 256·(p÷30) + wheel-index; p and p² reconstruct
    # from these, so the raw prime values are never stored. Small primes are the COMB_PRIMES const.
    # Sieving primes are bounded by feasibility (p ≲ 10^18), so this stays Int64 regardless of T.
    stride_primes::Vector{Int64}
    ncomb::Int               # active comb primes: COMB_PRIMES[1:ncomb] (those ≤ sieve_bound)
    lanes::Vector{BitVector} # the current window: 8 lanes, ≤ 64·seg_chunks bits each
    nchunks::Int             # current window: number of valid chunks (padding bits are zeroed)
end

# Optimal number of chunks per window. Low cap is 32 KiB/lane allowing for sieving within L1.
# As hi grows, we start sieving by bigger primes so we grow the window, but cap at MiB/lane to stay within L2.
_seg_chunks(hi::Integer) = Int(clamp(isqrt(hi) >>> 3, 1 << 12, 1 << 18))  # clamp before narrowing

# Chunk holding `lo` (window starts are chunk-aligned)
_first_chunk(lo::Integer) = fld(lo, 30) >>> 6

# `sieve_bound` limits sieving to primes ≤ it: survivors are then the `sieve_bound`-rough numbers
# in [lo, hi] (no prime factor ≤ sieve_bound), a superset of the primes. The default isqrt(hi)
# gives an exact prime sieve. 2, 3, 5 are always the caller's responsibility (the wheel skips them).
function SegmentedSieve(lo::T, hi::T; seg_chunks::Int = _seg_chunks(hi),
                        sieve_bound::Integer = isqrt(hi)) where {T<:Integer}
    lo = max(T(7), lo)
    # Decompose the stride primes (≥ COMB_THRESH; smaller ones are the comb) into 256·(p÷30) +
    # wheel-index once. Stream them with eachprime so we never allocate the full base-prime vector.
    # The sieving-prime ceiling `sq` is small (a feasible sieve can't hold more base primes), so it
    # fits Int even when the endpoints do not, and the packing bound below is far past any real sieve.
    stride_primes = Int64[]
    sqT = min(isqrt(hi), sieve_bound)
    sqT ≤ 30 * (typemax(Int64) ÷ 256) ||
        throw(ArgumentError("sieve_bound $sieve_bound exceeds the supported ceiling ~10^18"))
    sq = Int(sqT)
    if sq ≥ COMB_THRESH
        for p in eachprime(COMB_THRESH, sq)
            pq, pr = divrem(p, 30)
            push!(stride_primes, 256 * pq + WHEEL_IDX[pr + 1])
        end
    end
    ncomb = count(≤(sq), COMB_PRIMES)
    num_chunks = Int(min(seg_chunks, cld(fld(hi, 30) + 1, 64) - _first_chunk(lo)))
    lanes = [BitVector(undef, num_chunks << 6) for _ in 1:8]
    return SegmentedSieve{T}(zero(T), lo, hi, seg_chunks, stride_primes, ncomb, lanes, 0)
end

Base.IteratorSize(::Type{<:SegmentedSieve}) = Base.SizeUnknown()
Base.eltype(::Type{S}) where {S<:SegmentedSieve} = S

# Fill the next window into the lanes and yield the sieve itself (with `win_block`/`nchunks` set
# to the current window). Boundary bits (< lo, > hi, padding) are masked, so consumers need
# no per-prime range checks.
function Base.iterate(ss::SegmentedSieve{T}, win_chunk::T = _first_chunk(ss.lo)) where {T}
    win_block = win_chunk << 6                         # absolute block (T); relative offsets below are Int
    last_block = fld(ss.hi, 30)
    win_block > last_block && return nothing
    lanes = ss.lanes
    nblocks = Int(min(ss.seg_chunks << 6, last_block - win_block + 1))
    nchunks = cld(nblocks, 64)
    win_last_block = win_block + nblocks - 1
    win_max = win_last_block == last_block ? ss.hi : 30 * win_last_block + 29
    mp = isqrt(win_max)                                # cap the prime scan at √win_max, but never
    max_prime = mp > typemax(Int) ? typemax(Int) : Int(mp)  # above Int (all stored primes fit Int)

    # Loop-invariant boundary conditions (applied per lane below to make consumers range-check-free).
    lo_block = fld(ss.lo, 30)                          # values < lo occupy blocks ≤ lo_block
    lo_res = Int(ss.lo - 30 * lo_block)                # lo % 30
    tail_local = last_block - win_block                # local block of hi (if hi falls in this window)
    hi_res = ss.hi - 30 * last_block                   # hi % 30, overflow-safe threshold
    tail_bits = nblocks - (nchunks - 1) * 64           # valid bits in the final chunk (64 ⇒ no padding)

    # One lane fully before the next, so it stays cache-hot across all its primes: fill it, comb
    # the small primes (restoring each prime's own bit), stride the large ones, then mask boundaries.
    @inbounds for lane in 1:8
        chunks = lanes[lane].chunks
        fill!(chunks, ~UInt64(0))
        w = wheel[lane]
        for ci in 1:ss.ncomb                           # small primes ≤ sieve_bound: periodic comb AND
            p = COMB_PRIMES[ci]
            p > max_prime && break
            comb_clear!(chunks, COMB_TABLE[ci][lane], p, Int(mod(win_chunk, p)), nchunks)
            if win_block == 0 && p % 30 == w           # comb cleared p's own bit (p·1) and this is p's lane
                b = fld(p, 30)
                chunks[b >>> 6 + 1] |= UInt64(1) << (b & 63)
            end
        end
        for code in ss.stride_primes                   # large primes: scalar stride clear
            # Reconstruct p and p²÷30 from (pq, residue class) — no division. res_block is the phase:
            # lane-w multiples sit at blocks ≡ res_block (mod p). The only division left is mod p.
            pq, pri = fldmod(code, 256)
            r = wheel[pri]
            p = 30 * pq + r
            p > max_prime && break                     # primes are ascending; rest are too big for this window
            res_block = pq * STRIDE_KW[pri][lane] + STRIDE_C[pri][lane]
            # min_block = cld(p²-w, 30): first block ≥ p². p² can exceed Int for large-T sieves, so
            # this absolute block is computed in T; the relative `local_start` below is back in Int.
            pqT = T(pq)
            min_block = 30 * pqT * pqT + 2 * pqT * r + STRIDE_PSQ_DIV[pri] + (STRIDE_PSQ_RES[pri] > w)
            from_block = max(min_block, win_block)
            start_block = from_block + mod(res_block - from_block, p)
            local_start = Int(start_block - win_block)
            local_start < nblocks && stride_clear!(chunks, local_start, p, nblocks)
        end
        if win_block ≤ lo_block                        # first window: holds values < lo (incl. value 1)
            for block in win_block:(lo_block - 1)      # whole blocks below lo
                loc = Int(block - win_block)
                chunks[loc >>> 6 + 1] &= ~(UInt64(1) << (loc & 63))
            end
            if w < lo_res                              # partial block holding lo
                loc = Int(lo_block - win_block)
                chunks[loc >>> 6 + 1] &= ~(UInt64(1) << (loc & 63))
            end
        end
        if tail_local < nblocks && w > hi_res          # clear values > hi in the tail block
            loc = Int(tail_local)
            chunks[loc >>> 6 + 1] &= ~(UInt64(1) << (loc & 63))
        end
        if tail_bits < 64                              # zero the final chunk's padding
            chunks[nchunks] &= (UInt64(1) << tail_bits) - 1
        end
    end

    ss.win_block = win_block
    ss.nchunks = nchunks
    return (ss, win_chunk + ss.seg_chunks)
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
function each_lane_prime(f, ss::SegmentedSieve{T}) where {T}
    win_block = ss.win_block                                    # T; yielded values are T
    chunks = ntuple(i -> ss.lanes[i].chunks, Val(8))            # the 8 lanes' word arrays
    nchunks = ss.nchunks
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
mutable struct EachPrime{T<:Integer}
    lo::T
    hi::T
    sieve::Union{SegmentedSieve{T},Nothing}  # nothing when hi < 7
    next_chunk::T                            # first chunk of the next window to fetch
    buffer::Vector{T}                        # the current window's primes
    pos::Int                                 # next index into buffer
    phase::Int                               # 0,1,2 ⇒ emit 2,3,5; ≥3 ⇒ window mode
    verify::Bool                             # narrow window: sieve is rough, so isprime-filter survivors
end

"""
    eachprime([lo,] hi)

Lazily iterate the primes in `[lo, hi]` (from 2 if `lo` is omitted) in increasing order.
"""
function eachprime(lo::Integer, hi::Integer)
    lo ≤ hi || throw(ArgumentError("The condition lo ≤ hi must be met."))
    lo, hi = promote(lo, hi)
    T = typeof(hi)
    # Narrow window at large height: sieving all the way to √hi costs Θ(√hi) regardless of window
    # width W. When W < √hi/150 it is cheaper to presieve only to k ≈ W (killing composites down to
    # a small-factor floor) and BPSW-test the survivors — see docs/design/narrow-window-primes.md.
    # In that regime the sieve yields W-rough numbers (a superset of primes), so we isprime-filter.
    # A bound of W is the empirical optimum: below it, surviving composites (needing isprime) grow
    # faster than the sieving saved; above it, sieving cost grows. The bowl is shallow near W.
    sq = isqrt(hi)
    W = hi - max(T(7), lo)
    bound = W < sq ÷ 150 ? clamp(W, T(1000), sq) : sq
    verify = bound < sq
    sieve = hi < 7 ? nothing : SegmentedSieve(max(T(7), lo), hi; sieve_bound = bound)
    next_chunk = sieve === nothing ? zero(T) : _first_chunk(sieve.lo)
    return EachPrime{T}(lo, hi, sieve, next_chunk, T[], 1, 0, verify)
end
eachprime(hi::Integer) = eachprime(1, hi)

Base.IteratorSize(::Type{<:EachPrime}) = Base.SizeUnknown()
Base.eltype(::Type{EachPrime{T}}) where {T} = T

function Base.iterate(ep::EachPrime{T}, ::Any = nothing) where {T}
    while ep.phase < 3                          # 2, 3, 5 precede every wheel prime
        p = T((2, 3, 5)[ep.phase + 1])
        ep.phase += 1
        ep.lo ≤ p ≤ ep.hi && return (p, nothing)
    end
    while ep.pos > length(ep.buffer)            # refill from the next window(s)
        ep.sieve === nothing && return nothing
        next = iterate(ep.sieve, ep.next_chunk)
        next === nothing && return nothing
        window, ep.next_chunk = next            # window === ep.sieve (it yields itself)
        empty!(ep.buffer)
        if ep.verify                            # rough sieve: keep only the true primes among survivors
            each_lane_prime(v -> (isprime(v) && push!(ep.buffer, v)), window)
        else
            each_lane_prime(v -> push!(ep.buffer, v), window)
        end
        ep.pos = 1
    end
    p = ep.buffer[ep.pos]
    ep.pos += 1
    return (p, nothing)
end
