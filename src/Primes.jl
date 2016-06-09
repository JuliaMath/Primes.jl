# This includes parts that were formerly a part of Julia. License is MIT: http://julialang.org/license
__precompile__()
module Primes

import DataStructures: SortedDict

if VERSION >= v"0.5.0-dev+4340"
    if isdefined(Base,:isprime)
        import Base: isprime, primes, primesmask, factor
    else
        export isprime, primes, primesmask, factor, factorvec
    end
end

# Primes generating functions
#     https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
#     https://en.wikipedia.org/wiki/Wheel_factorization
#     http://primesieve.org
#     Jonathan Sorenson, "An analysis of two prime number sieves", Computer Science Technical Report Vol. 1028, 1991
const wheel         = [4,  2,  4,  2,  4,  6,  2,  6]
const wheel_primes  = [7, 11, 13, 17, 19, 23, 29, 31]
const wheel_indices = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7]

@inline function wheel_index(n)
    (d, r) = divrem(n - 1, 30)
    8d + wheel_indices[r+2]
end
@inline function wheel_prime(n)
    (d, r) = ((n - 1) >>> 3, (n - 1) & 7)
    30d + wheel_primes[r+1]
end

function _primesmask(limit::Int)
    limit < 7 && throw(ArgumentError("limit must be at least 7, got $limit"))
    n = wheel_index(limit)
    m = wheel_prime(n)
    sieve = ones(Bool, n)
    @inbounds for i = 1:wheel_index(isqrt(limit))
        if sieve[i]; p = wheel_prime(i)
            q = p * p
            j = (i - 1) & 7 + 1
            while q ≤ m
                sieve[wheel_index(q)] = false
                q = q + wheel[j] * p
                j = j & 7 + 1
            end
        end
    end
    return sieve
end

function _primesmask(lo::Int, hi::Int)
    7 ≤ lo ≤ hi || throw(ArgumentError("the condition 7 ≤ lo ≤ hi must be met"))
    lo == 7 && return _primesmask(hi)
    wlo, whi = wheel_index(lo - 1), wheel_index(hi)
    m = wheel_prime(whi)
    sieve = ones(Bool, whi - wlo)
    hi < 49 && return sieve
    small_sieve = _primesmask(isqrt(hi))
    @inbounds for i = 1:length(small_sieve)  # don't use eachindex here
        if small_sieve[i]; p = wheel_prime(i)
            j = wheel_index(2 * div(lo - p - 1, 2p) + 1)
            q = p * wheel_prime(j + 1); j = j & 7 + 1
            while q ≤ m
                sieve[wheel_index(q)-wlo] = false
                q = q + wheel[j] * p
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
    0 < lo ≤ hi || throw(ArgumentError("the condition 0 < lo ≤ hi must be met"))
    sieve = falses(hi - lo + 1)
    lo ≤ 2 ≤ hi && (sieve[3-lo] = true)
    lo ≤ 3 ≤ hi && (sieve[4-lo] = true)
    lo ≤ 5 ≤ hi && (sieve[6-lo] = true)
    hi < 7 && return sieve
    wheel_sieve = _primesmask(max(7, lo), hi)
    lsi = lo - 1
    lwi = wheel_index(lsi)
    @inbounds for i = 1:length(wheel_sieve)   # don't use eachindex here
        sieve[wheel_prime(i + lwi) - lsi] = wheel_sieve[i]
    end
    return sieve
end
primesmask{T<:Integer}(lo::T, hi::T) = lo <= hi <= typemax(Int) ? primesmask(Int(lo), Int(hi)) :
    throw(ArgumentError("both endpoints of the interval to sieve must be ≤ $(typemax(Int)), got $lo and $hi"))

primesmask(limit::Int) = primesmask(1, limit)
primesmask(n::Integer) = n <= typemax(Int) ? primesmask(Int(n)) :
    throw(ArgumentError("requested number of primes must be ≤ $(typemax(Int)), got $n"))

"""
    primes([lo,] hi)

Returns a collection of the prime numbers (from `lo`, if specified) up to `hi`.
"""
function primes(lo::Int, hi::Int)
    lo ≤ hi || throw(ArgumentError("the condition lo ≤ hi must be met"))
    list = Int[]
    lo ≤ 2 ≤ hi && push!(list, 2)
    lo ≤ 3 ≤ hi && push!(list, 3)
    lo ≤ 5 ≤ hi && push!(list, 5)
    hi < 7 && return list
    sizehint!(list,  5 + floor(Int, hi / (log(hi) - 1.12) -  lo / (log(max(lo,2)) - 1.12*(lo > 7))) ) # http://projecteuclid.org/euclid.rmjm/1181070157
    sieve = _primesmask(max(7, lo), hi)
    lwi = wheel_index(lo - 1)
    @inbounds for i = 1:length(sieve)   # don't use eachindex here
        sieve[i] && push!(list, wheel_prime(i + lwi))
    end
    return list
end
primes(n::Int) = primes(1, n)

const PRIMES = primes(2^16)

"""
    isprime(x::Integer) -> Bool

Returns `true` if `x` is prime, and `false` otherwise.

```jldoctest
julia> isprime(3)
true
```
"""
function isprime(n::Integer)
    # Small precomputed primes + Miller-Rabin for primality testing:
    #     https://en.wikipedia.org/wiki/Miller–Rabin_primality_test
    #     https://github.com/JuliaLang/julia/issues/11594
    for m in (2,3,5,7,11,13,17,19,23)
        n % m == 0 && return n == m
    end
    n < 841 && return n > 1
    s = trailing_zeros(n-1)
    d = (n-1) >>> s
    for a in witnesses(n)::Tuple{Vararg{Int}}
        x = powermod(a,d,n)
        x == 1 && continue
        t = s
        while x != n-1
            (t-=1) <= 0 && return false
            x = oftype(n, widemul(x,x) % n)
            x == 1 && return false
        end
    end
    return true
end

"""
    isprime(x::BigInt, [reps = 25]) -> Bool

Probabilistic primality test. Returns `true` if `x` is prime with high probability (pseudoprime);
and `false` if `x` is composite (not prime). The false positive rate is about `0.25^reps`.
`reps = 25` is considered safe for cryptographic applications (Knuth, Seminumerical Algorithms).

```jldoctest
julia> isprime(big(3))
true
```
"""
isprime(x::BigInt, reps=25) = ccall((:__gmpz_probab_prime_p,:libgmp), Cint, (Ptr{BigInt}, Cint), &x, reps) > 0


# Miller-Rabin witness choices based on:
#     http://mathoverflow.net/questions/101922/smallest-collection-of-bases-for-prime-testing-of-64-bit-numbers
#     http://primes.utm.edu/prove/merged.html
#     http://miller-rabin.appspot.com
#     https://github.com/JuliaLang/julia/issues/11594
#     Forišek and Jančina, "Fast Primality Testing for Integers That Fit into a Machine Word", 2015
#         (in particular, see function FJ32_256, from which the hash and bases were taken)
const bases = UInt16[15591,2018,166,7429,8064,16045,10503,4399,1949,1295,2776,
3620,560,3128,5212,2657,2300,2021,4652,1471,9336,4018,2398,20462,10277,8028,
2213,6219,620,3763,4852,5012,3185,1333,6227,5298,1074,2391,5113,7061,803,1269,
3875,422,751,580,4729,10239,746,2951,556,2206,3778,481,1522,3476,481,2487,3266,
5633,488,3373,6441,3344,17,15105,1490,4154,2036,1882,1813,467,3307,14042,6371,
658,1005,903,737,1887,7447,1888,2848,1784,7559,3400,951,13969,4304,177,41,
19875,3110,13221,8726,571,7043,6943,1199,352,6435,165,1169,3315,978,233,3003,
2562,2994,10587,10030,2377,1902,5354,4447,1555,263,27027,2283,305,669,1912,601,
6186,429,1930,14873,1784,1661,524,3577,236,2360,6146,2850,55637,1753,4178,8466,
222,2579,2743,2031,2226,2276,374,2132,813,23788,1610,4422,5159,1725,3597,3366,
14336,579,165,1375,10018,12616,9816,1371,536,1867,10864,857,2206,5788,434,8085,
17618,727,3639,1595,4944,2129,2029,8195,8344,6232,9183,8126,1870,3296,7455,
8947,25017,541,19115,368,566,5674,411,522,1027,8215,2050,6544,10049,614,774,
2333,3007,35201,4706,1152,1785,1028,1540,3743,493,4474,2521,26845,8354,864,
18915,5465,2447,42,4511,1660,166,1249,6259,2553,304,272,7286,73,6554,899,2816,
5197,13330,7054,2818,3199,811,922,350,7514,4452,3449,2663,4708,418,1621,1171,
3471,88,11345,412,1559,194]

function _witnesses(n::UInt64)
    i = ((n >> 16) $ n) * 0x45d9f3b
    i = ((i >> 16) $ i) * 0x45d9f3b
    i = ((i >> 16) $ i) & 255 + 1
    @inbounds return (Int(bases[i]),)
end
witnesses(n::Integer) =
        n < 4294967296      ? _witnesses(UInt64(n)) :
        n < 2152302898747   ? (2,3,5,7,11) :
        n < 3474749660383   ? (2,3,5,7,11,13) :
                              (2,325,9375,28178,450775,9780504,1795265022)

isprime(n::UInt128) =
    n <= typemax(UInt64) ? isprime(UInt64(n)) : isprime(BigInt(n))
isprime(n::Int128) = n < 2 ? false :
    n <= typemax(Int64)  ? isprime(Int64(n))  : isprime(BigInt(n))


# Trial division of small (< 2^16) precomputed primes +
# Pollard rho's algorithm with Richard P. Brent optimizations
#     https://en.wikipedia.org/wiki/Trial_division
#     https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
#     http://maths-people.anu.edu.au/~brent/pub/pub051.html
#
"""
    factor(n) -> DataStructures.SortedDict

Compute the prime factorization of an integer `n`. Returns a sorted dictionary. The
keys of the dictionary correspond to the factors, and hence are of the same type as `n`.
The value associated with each key indicates the number of times the factor appears in the
factorization.

```jldoctest
julia> factor(100) # == 2*2*5*5
DataStructures.SortedDict{Int64,Int64,Base.Order.ForwardOrdering} with 2 entries:
  2 => 2
  5 => 2
```
"""
function factor{T<:Integer}(n::T)
    0 < n || throw(ArgumentError("number to be factored must be > 0, got $n"))
    h = SortedDict{T,Int,Base.Order.ForwardOrdering}()
    n == 1 && return h
    isprime(n) && (h[n] = 1; return h)
    local p::T
    for p in PRIMES
        if n % p == 0
            h[p] = get(h,p,0)+1
            n = div(n,p)
            while n % p == 0
                h[p] = get(h,p,0)+1
                n = div(n,p)
            end
            n == 1 && return h
            isprime(n) && (h[n] = 1; return h)
        end
    end
    T <: BigInt || widemul(n-1,n-1) <= typemax(n) ? pollardfactors!(n, h) : pollardfactors!(widen(n), h)
end

function pollardfactors!{T<:Integer,K<:Integer}(n::T, h::Associative{K,Int})
    while true
        local c::T = rand(1:(n-1)), G::T = 1, r::K = 1, y::T = rand(0:(n-1)), m::K = 1900
        local ys::T, q::T = 1, x::T
        while c == n - 2
            c = rand(1:(n-1))
        end
        while G == 1
            x = y
            for i in 1:r
                y = (y*y)%n
                y = (y+c)%n
            end
            local k::K = 0
            G = 1
            while k < r && G == 1
                for i in 1:(m>(r-k)?(r-k):m)
                    ys = y
                    y = (y*y)%n
                    y = (y+c)%n
                    q = (q*(x>y?x-y:y-x))%n
                end
                G = gcd(q,n)
                k = k + m
            end
            r = 2 * r
        end
        G == n && (G = 1)
        while G == 1
            ys = (ys*ys)%n
            ys = (ys+c)%n
            G = gcd(x>ys?x-ys:ys-x,n)
        end
        if G != n
            isprime(G) ? h[G] = get(h,G,0) + 1 : pollardfactors!(G,h)
            G2 = div(n,G)
            isprime(G2) ? h[G2] = get(h,G2,0) + 1 : pollardfactors!(G2,h)
            return h
        end
    end
end

"""
    factorvec(n::Integer) -> Vector

Compute the prime factorization of `n` with multiplicities. Returns a vector with
the same type as `n`. The product of the returned vector will equal `n`.

```jldoctest
julia> factorvec(100)
4-element Array{Int64,1}:
 2
 2
 5
 5
```
"""
function factorvec{T<:Integer}(n::T)
    if n < 1
        throw(ArgumentError("number to be factored must be > 0, got $n"))
    elseif n == 1
        return Vector{T}(0)  # For consistency with factor(1)
    else
        return mapreduce(collect, vcat, [repeated(k,v) for (k,v) in factor(n)])
    end
end

end # module
