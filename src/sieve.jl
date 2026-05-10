# Primes generating functions
#     https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
#     https://en.wikipedia.org/wiki/Wheel_factorization
#     http://primesieve.org
#     Jonathan Sorenson, "An analysis of two prime number sieves", Computer Science Technical Report Vol. 1028, 1991
const wheel         = (4,  2,  4,  2,  4,  6,  2,  6)
const wheel_primes  = (7, 11, 13, 17, 19, 23, 29, 31)
const wheel_indices = (0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7)

@inline function wheel_index(n)
    d, r = divrem(n - 1, 30)
    return 8d + wheel_indices[r + 2]
end
@inline function wheel_prime(n)
    d, r = (n - 1) >>> 3, (n - 1) & 7
    return 30d + wheel_primes[r + 1]
end

# `wheel_prime(x)%p` is periodic every `8p` (8 because `wheel` has 8 elements)
# Therefore, it is also periodic every `64p`.
# Since a `BitVector` stores `UInt64` chunks, each prime has the same effect every `p` chunks.
# Since prime `p` modifies `1/p` bits, small primes<64 modify more than 1 bit per UInt64
# which makes it worth it to sieve them 1 chunk at a time.
const PRESIEVE_PRIMES = (7,11,13,17,19)
const PRESIEVE_CHUNKS = Vector{Vector{UInt64}}(undef, length(PRESIEVE_PRIMES))
for (i,p) in enumerate(PRESIEVE_PRIMES)
   PRESIEVE_CHUNKS[i] = BitVector([wheel_prime(x)%p!=0 for x in 1:(64*p)]).chunks
end

function _primesmask(limit::Int)
    limit < 7 && throw(ArgumentError("The condition limit ≥ 7 must be met."))
    n = wheel_index(limit)
    m = wheel_prime(n)
    sieve = BitVector(undef,n)
    # First sieve small primes 64 bits at a time.
    # Loop order here is slightly more computation, but reduces passes over memory
    @inbounds for i in eachindex(sieve.chunks)
        presieve_mask = ~UInt64(0)
        for (j,p) in enumerate(PRESIEVE_PRIMES)
             presieve_mask &= PRESIEVE_CHUNKS[j][(i-1) % p + 1]
        end
        sieve.chunks[i] = presieve_mask
    end
    # undo the fact that presieving sets the `PRESIEVE_PRIMES` to `false`
    # 0x1f is the bitmask for setting sieve[1:5] to true (and doesn't need to check the length of sieve)
    sieve.chunks[1] |= 0x1f

    # Then sieve remaining primes 1 match at a time.
    @inbounds for i = (length(PRESIEVE_PRIMES)+1):wheel_index(isqrt(limit))
        if sieve[i]
            p = wheel_prime(i)
            q = p^2
            j = (i - 1) & 7 + 1
            while q ≤ m
                sieve[wheel_index(q)] = false
                q += wheel[j] * p
                j = j & 7 + 1
            end
        end
    end
    return sieve
end

function _primesmask(lo::Int, hi::Int)
    7 ≤ lo ≤ hi || throw(ArgumentError("The condition 7 ≤ lo ≤ hi must be met."))
    lo == 7 && return _primesmask(hi)
    wlo, whi = wheel_index(lo - 1), wheel_index(hi)
    m = wheel_prime(whi)
    sieve = ones(Bool, whi - wlo)
    hi < 49 && return sieve
    small_sieve = _primesmask(isqrt(hi))
    @inbounds for i = 1:length(small_sieve)  # don't use eachindex here
        if small_sieve[i]
            p = wheel_prime(i)
            j = wheel_index(2 * div(lo - p - 1, 2p) + 1)
            r = widemul(p, wheel_prime(j + 1))
            r > m && continue # use widemul to avoid r <= m caused by overflow
            j = j & 7 + 1
            q = Int(r)
            # q < 0 indicates overflow when incrementing q inside loop
            while 0 ≤ q ≤ m
                sieve[wheel_index(q) - wlo] = false
                q += wheel[j] * p
                j = j & 7 + 1
            end
        end
    end
    return sieve
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
    wheel_sieve = _primesmask(max(7, lo), hi)
    lsi = lo - 1
    lwi = wheel_index(lsi)
    @inbounds for i = 1:length(wheel_sieve)   # don't use eachindex here
        sieve[wheel_prime(i + lwi) - lsi] = wheel_sieve[i]
    end
    return sieve
end
primesmask(lo::Integer, hi::Integer) = lo ≤ hi ≤ typemax(Int) ? primesmask(Int(lo), Int(hi)) :
    throw(ArgumentError("Both endpoints of the interval to sieve must be ≤ $(typemax(Int)), got $lo and $hi."))

primesmask(limit::Int) = primesmask(1, limit)
primesmask(n::Integer) = n ≤ typemax(Int) ? primesmask(Int(n)) :
    throw(ArgumentError("Requested number of primes must be ≤ $(typemax(Int)), got $n."))
