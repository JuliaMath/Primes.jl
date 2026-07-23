# This file includes code that was formerly a part of Julia. License is MIT: http://julialang.org/license
module Primes

using Base.Iterators: repeated, rest

import Base: iterate, eltype, IteratorSize, IteratorEltype
using Base: BitSigned
using Base.Checked: checked_neg
using IntegerMathUtils

export isprime, primes, primesmask, eachprime, factor, eachfactor, divisors, ismersenneprime, isrieselprime,
       nextprime, nextprimes, prevprime, prevprimes, prime, prodfactors, radical, totient

include("factorization.jl")
include("sieve.jl")

"""
    primes([lo,] hi)

Returns a collection of the prime numbers (from `lo`, if specified) up to `hi`.
"""
function primes(lo::Int, hi::Int)
    lo ≤ hi || throw(ArgumentError("The condition lo ≤ hi must be met."))
    list = Int[]
    lo ≤ 2 ≤ hi && push!(list, 2)
    lo ≤ 3 ≤ hi && push!(list, 3)
    lo ≤ 5 ≤ hi && push!(list, 5)
    hi < 7 && return list
    lo = max(2, lo)
    sizehint!(list, 5 + floor(Int, hi / (log(hi) - 1.12) - lo / (log(lo) - 1.12 * (lo > 7))) ) # http://projecteuclid.org/euclid.rmjm/1181070157
    # collect window by window; the sieve masks values < max(7, lo) and > hi
    for s in SegmentedSieve(max(7, lo), hi)
        each_lane_prime(s) do p
            push!(list, p)
        end
    end
    return list
end
primes(n::Int) = primes(1, n)

function _generate_min_factors(limit)
    function min_factor(n)
        n < 4 && return n
        for i in 3:2:isqrt(n)
            n%i == 0 && return i
        end
        return n
    end
    res = Int[]
    for i in 3:2:limit
        m = min_factor(i)
        push!(res, m==i ? 1 : m)
    end
    return res
end

const N_SMALL_FACTORS = 2^16
const _MIN_FACTOR = UInt8.(_generate_min_factors(N_SMALL_FACTORS))
# _min_factor(n) = the minimum factor of n for odd n, if 1<n<N_SMALL_FACTORS
function _min_factor(n::T) where T<:Integer
    m = _MIN_FACTOR[n>>1]
    return m==1 ? n : T(m)
end

"""
    isprime(n::Integer) -> Bool

Returns for values in the range of an INT64 variable:  `true` if `n` is prime, and `false`
otherwise for bigger values: `true` if `n` is probably prime, and `false` otherwise
(false-positive rate = 0.25^reps with reps=25 --> considered safe)

More detailed:
- for even numbers: returns deterministic and correct results regardless of type
- for values in the range of an `Int64` variable: returns deterministic and correct results
  (by Lookup-tables, trial-division, Miller-Rabin, Lucas-Test)
- for bigger values: returns probabilistic results from GNU Multiple Precision Arithmetic Library

```julia
julia> isprime(3)
true

julia> isprime(4)
false
```
"""
function isprime(n::T, reps::Integer=25) where {T<:Integer}
    n ≤ typemax(Int64) && return isprime(Int64(n))
    # GMP-style probable-prime test in the operand's own arithmetic
    # BPSW, then reps - 24 extra Miller-Rabin rounds so
    # the error bound tracks GMP's for the same reps.
    iseven(n) && return false
    for m in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
        n % m == 0 && return false
    end
    miller_rabbin_test(2, n) || return false

    if Base.hastypemax(T) && n >= isqrt(typemax(T))
        lucas_test(widen(n)) || return false
    else
        lucas_test(n) || return false
    end
    for _ in 1:max(reps - 24, 1)
        miller_rabbin_test(rand(3:n-2), n) || return false
    end
    return true
end

# Uses a polyalgorithm depending on the size of n.
# n < 2^16: lookup table (we already have this table because it helps factor also)
# n < 2^32: trial division + Miller-Rabbin test with base chosen by
#           Forišek and Jančina, "Fast Primality Testing for Integers That Fit into a Machine Word", 2015
#           (in particular, see function FJ32_256, from which the hash and bases were taken)
# n < 2^64: Baillie–PSW for primality testing.
#           Specifically, this consists of a Miller-Rabbin test and a Lucas test
# For more background on fast primality testing, see:
#     http://ntheory.org/pseudoprimes.html
#     http://ntheory.org/pseudoprimes.html
function isprime(n::Int64)
    iseven(n) && return n == 2
    if n < N_SMALL_FACTORS
        n < 2 && return false
        return _min_factor(n) == n
    end
    for m in (3, 5, 7, 11, 13, 17, 19, 23)
        n % m == 0 && return false
    end
    if n < 2^32
        return miller_rabbin_test(_witnesses(UInt64(n)), n)
    end
    miller_rabbin_test(2, n) || return false
    return lucas_test(widen(n))
end

const bases = UInt16[
    15591,  2018,   166,  7429,  8064, 16045, 10503,  4399,  1949,  1295,  2776,  3620,
      560,  3128,  5212,  2657,  2300,  2021,  4652,  1471,  9336,  4018,  2398, 20462,
    10277,  8028,  2213,  6219,   620,  3763,  4852,  5012,  3185,  1333,  6227,  5298,
     1074,  2391,  5113,  7061,   803,  1269,  3875,   422,   751,   580,  4729, 10239,
      746,  2951,   556,  2206,  3778,   481,  1522,  3476,   481,  2487,  3266,  5633,
      488,  3373,  6441,  3344,    17, 15105,  1490,  4154,  2036,  1882,  1813,   467,
     3307, 14042,  6371,   658,  1005,   903,   737,  1887,  7447,  1888,  2848,  1784,
     7559,  3400,   951, 13969,  4304,   177,    41, 19875,  3110, 13221,  8726,   571,
     7043,  6943,  1199,   352,  6435,   165,  1169,  3315,   978,   233,  3003,  2562,
     2994, 10587, 10030,  2377,  1902,  5354,  4447,  1555,   263, 27027,  2283,   305,
      669,  1912,   601,  6186,   429,  1930, 14873,  1784,  1661,   524,  3577,   236,
     2360,  6146,  2850, 55637,  1753,  4178,  8466,   222,  2579,  2743,  2031,  2226,
     2276,   374,  2132,   813, 23788,  1610,  4422,  5159,  1725,  3597,  3366, 14336,
      579,   165,  1375, 10018, 12616,  9816,  1371,   536,  1867, 10864,   857,  2206,
     5788,   434,  8085, 17618,   727,  3639,  1595,  4944,  2129,  2029,  8195,  8344,
     6232,  9183,  8126,  1870,  3296,  7455,  8947, 25017,   541, 19115,   368,   566,
     5674,   411,   522,  1027,  8215,  2050,  6544, 10049,   614,   774,  2333,  3007,
    35201,  4706,  1152,  1785,  1028,  1540,  3743,   493,  4474,  2521, 26845,  8354,
      864, 18915,  5465,  2447,    42,  4511,  1660,   166,  1249,  6259,  2553,   304,
      272,  7286,    73,  6554,   899,  2816,  5197, 13330,  7054,  2818,  3199,   811,
      922,   350,  7514,  4452,  3449,  2663,  4708,   418,  1621,  1171,  3471,    88,
    11345,   412,  1559,   194
]

function _witnesses(n::UInt64)
    i = xor((n >> 16), n) * 0x45d9f3b
    i = xor((i >> 16), i) * 0x45d9f3b
    i = xor((i >> 16), i) & 255 + 1
    @inbounds return Int(bases[i])
end

function miller_rabbin_test(a, n::T) where T<:Integer
    s = trailing_zeros(n - 1)
    d = (n - 1) >>> s
    x::T = powermod(a, d, n)
    if x != 1
        t = s
        while x != n - 1
            (t -= 1) ≤ 0 && return false
            x = widemul(x, x) % n
            x == 1 && return false
        end
    end
    return true
end

# Requires n^2 + 2n < typemax(T)
function lucas_test(n::T) where T<:Integer
    isqrt(n)^2 == n && return false
    # Selfridge parameters: P = 1 and
    # D the first of 5, -7, 9, -11, ... with kronecker(D, n) == -1.
    local Dr::T, Q::T, kr
    for (s, d) in zip(Iterators.cycle((1, -1)), Iterators.countfrom(5, 2))
        Dr = s > 0 ? T(d) : n - T(d)               # D mod n
        Q  = s > 0 ? n - T((d - 1) >> 2) : T((d + 1) >> 2)  # (1-D)/4 mod n
        kr = kronecker(Dr, n)
        kr != 1 && break
    end
    kr == 0 && return false
    U::T, V::T, Qk::T = 1, 1, Q
    k::T = n + 1
    trail = trailing_zeros(k)
    k >>= trail
    # get digits 1 at a time since digits allocates
    for b in ndigits(k, base=2)-2:-1:0
        U = mod(U*V, n)
        # V*V+2(n-Qk) rather than V*V-2*Qk to prevent Unsigned underflow
        V = mod(V * V + 2*(n - Qk), n)
        Qk = mod(Qk * Qk, n)
        if isodd(k>>b)
            Qk = mod(Qk * Q, n)
            U, V = U + V, V + U*Dr
            # adding n makes them even
            # so we can divide by 2 without causing problems
            isodd(U) && (U += n)
            isodd(V) && (V += n)
            U = mod(U >> 1, n)
            V = mod(V >> 1, n)
        end
    end
    U == 0 && return true
    for _ in 1:trail
        V == 0 && return true
        V = mod(V*V + 2*(n - Qk), n)
        Qk = mod(Qk * Qk, n)
    end
    return false
end
"""
    isprime(x::BigInt, [reps = 25]) -> Bool
Probabilistic primality test. Returns `true` if `x` is prime with high probability (pseudoprime);
and `false` if `x` is composite (not prime). The false positive rate is about `0.25^reps`.
`reps = 25` is considered safe for cryptographic applications (Knuth, Seminumerical Algorithms).

is_probably_prime is inherited from the module IntegerMathUtils that provides a wrapper to access
functionality from the  GNU Multiple Precision Arithmetic Library (GMP) library
```julia
julia> isprime(big(3))
true
```
"""
isprime(x::BigInt, reps::Integer=25) = x>1 && is_probably_prime(x; reps=reps)

struct FactorIterator{T<:Integer}
    n::T
    FactorIterator(n::T) where {T} = new{T}(n)
end

IteratorSize(::Type{<:FactorIterator}) = Base.SizeUnknown()
IteratorEltype(::Type{<:FactorIterator}) = Base.HasEltype()
eltype(::Type{FactorIterator{T}}) where {T} = Tuple{T, Int}
Base.isempty(f::FactorIterator) = f.n == 1

# Iterates over the factors of n in an arbitrary order
# Uses a variety of algorithms depending on the size of n to find a factor.
#     https://en.algorithmica.org/hpc/algorithms/factorization
# Cache of small factors for small numbers (n < 2^16)
# Trial division of small (< 2^16) precomputed primes and
# Lenstra elliptic curve algorithm
#     https://en.wikipedia.org/wiki/Trial_division
#     https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization
#

"""
    eachfactor(n::Integer)->FactorIterator
Returns a lazy iterator of factors of `n` in `(factor, multiplicity)` pairs.
This can be very useful for computing multiplicative functions since for small numbers (e.g. numbers with no factor `>2^16`),
allocating the storage required for `factor(n)` can introduce significant overhead.
"""
eachfactor(n::Integer) = FactorIterator(n)

# state[1] is the current number to factor (this decreases when factors are found)
# state[2] is the prime to start trial division with.
function iterate(f::FactorIterator{T}, state=(f.n, T(3))) where T
    n, p::T = state
    if n <= p
        n == 1 && return nothing
        if n < 0
            # if n is typemin, we can't negate it properly
            # instead we set p=n which we can detect here.
            if isa(n, BitSigned) && n == typemin(T)
                if n != p
                    return (T(-1), 1), (n, n)
                end
                return (T(2), 8 * sizeof(T) - 1), T.((1, 1))
            end
            return  (T(-1), 1), (-n, p)
        end
        n == 0 && return (T(n), 1), (T(1), p)
    end
    tz = trailing_zeros(n)
    tz>0 && return (T(2), Int(tz)), (n >> tz, p)
    if n <= N_SMALL_FACTORS
        p = T(_min_factor(n))
        num_p = 1
        while true
            n = div(n, p)
            n == 1 && break
            _min_factor(n) == p || break
            num_p += 1
        end
        return (p, num_p), (n, p)
    elseif p == 3 && isprime(n)
        return (n, 1), (T(1), n)
    end
    for p in p:T(2):T(N_SMALL_FACTORS)
        _min_factor(p) == p || continue
        num_p = 0
        while true
            q, r = divrem(n, T(p)) # T(p) so julia <1.9 uses fast `divrem` for `BigInt`
            r == 0 || break
            num_p += 1
            n = q
        end
        if num_p > 0
            return (p, num_p), (n, p+2)
        end
        p*p > n && break
    end
    # if n < 2^32, then if it wasn't prime, we would have found the factors by trial division
    if n <= 2^32 || isprime(n)
        return (n, 1), (T(1), n)
    end
    if ispower(n) # check for n=root^r
        r = find_exponent(n)
        root = iroot(n, r)
        if isprime(root)
            return (root, r), (1, root+2)
        else
            f = first(eachfactor(root))
            return f, (n÷f, root+2)
        end
    end

    # lenstrafactor's modular arithmetic needs intermediates up to 2n^2
    fits = !Base.hastypemax(T) || widemul(n, n) ≤ typemax(n) ÷ 2
    p = fits ? lenstrafactor(n) : lenstrafactor(widen(n))
    num_p = 0
    while true
        q, r = divrem(n, p)
        r != 0 && return (p, num_p), (n, p)
        num_p += 1
        n = q
    end
end

# the curve arithmetic uses plain +/- with negative intermediates, so run
# unsigned inputs in the equally-sized signed domain (the widening check
# already caps n^2 ≤ typemax(n)÷2, so n fits the signed type)
lenstrafactor(n::Unsigned) = lenstrafactor(signed(typeof(n))(n))

function lenstrafactor(n::T; use_stage2::Bool=true) where{T<:Integer}
    # bounds and runs per bound taken from
    # https://www.rieselprime.de/ziki/Elliptic_curve_method
    # the 200 tier is a cheap pre-pass that catches most factors ≤ ~20 bits
    B1s = Int[200, 2e3, 11e3, 5e4, 25e4, 1e6, 3e6, 11e6,
                43e6, 11e7, 26e7, 85e7, 29e8, 76e8, 25e9]
    # published counts assume GMP-ECM's much larger B2; scaled 1.5x for B2 = 25*B1
    runs = Int[6, 38, 135, 450, 1050, 2700, 7650, 2700, 15900,
              29000, 73500, 186000, 315000, 510000, 15*10^5, 15*10^6]
    for (B1, run) in zip(B1s, runs)
        small_primes = primes(B1)
        for σ in 6:2*run+6
            res = lenstra_stage_1(n, σ, small_primes, B1; use_stage2)
            if res != 1
                return isprime(res) ? res : lenstrafactor(res; use_stage2)
            end
        end
    end
    throw(ArgumentError("This number is too big to be factored with this algorithm effectively"))
end

# x-only arithmetic on Montgomery curves B*y^2 = x^3 + A*x^2 + x, projective X:Z.
# mod-reduced values live in [0, N), but intermediate sums/differences may reach
# (-N, 2N), so this arithmetic needs signed T and |products| ≤ 2N^2 to fit in T —
# both guaranteed by the callers (widening check + unsigned-to-signed conversion).
# _addmod keeps sums that get squared in [0, N); squaring an unreduced sum would
# need 4N^2 of headroom, i.e. cost half a bit of the widening threshold.
@inline _addmod(a, b, N) = (c = a + b; c ≥ N ? c - N : c)

# doubling, with a24 = (A+2)/4
@inline function _xdbl(X, Z, a24, N)
    s = _addmod(X, Z, N)
    d = X - Z
    s2 = mod(s*s, N)
    d2 = mod(d*d, N)
    t = s2 - d2     # 4XZ
    return mod(s2*d2, N), mod(t*(d2 + mod(a24*t, N)), N)
end

# differential addition: P+Q given their difference D = P-Q
@inline function _xadd(XP, ZP, XQ, ZQ, XD, ZD, N)
    u = mod((XP - ZP)*(XQ + ZQ), N)
    v = mod((XP + ZP)*(XQ - ZQ), N)
    s = _addmod(u, v, N)
    d = u - v
    return mod(ZD*mod(s*s, N), N), mod(XD*mod(d*d, N), N)
end

# Montgomery ladder: k*(X:Z) for k ≥ 1
@inline function _xladder(k, X, Z, a24, N)
    k == 1 && return X, Z
    X0, Z0 = X, Z
    X1, Z1 = _xdbl(X, Z, a24, N)
    for i in 8*sizeof(k)-leading_zeros(k)-2:-1:0
        if (k >> i) & 1 == 0
            X1, Z1 = _xadd(X0, Z0, X1, Z1, X, Z, N)
            X0, Z0 = _xdbl(X0, Z0, a24, N)
        else
            X0, Z0 = _xadd(X1, Z1, X0, Z0, X, Z, N)
            X1, Z1 = _xdbl(X1, Z1, a24, N)
        end
    end
    return X0, Z0
end

# One ECM stage-1 attempt on the Suyama curve with parameter σ ≥ 6.
# Returns a nontrivial factor of N, or 1 if this curve fails.
function lenstra_stage_1(N::T, σ, small_primes, plimit; use_stage2::Bool=true) where T <: Integer
    # Suyama parametrization: guarantees torsion 12 and a known starting point
    u = mod(T(σ)*T(σ) - 5, N)
    v = mod(4*T(σ), N)
    X = mod(mod(u*u, N)*u, N)   # u^3
    Z = mod(mod(v*v, N)*v, N)   # v^3
    # a24 = (v-u)^3*(3u+v) / (16*u^3*v); a non-invertible denominator is itself a factor
    den = mod(mod(T(16)*X, N)*v, N)
    g, di, _ = gcdx(den, N)
    g > 1 && return g == N ? T(1) : g
    w = v - u
    num = mod(mod(mod(w*w, N)*w, N)*mod(3u + v, N), N)
    a24 = mod(num*mod(di, N), N)
    for B in small_primes
        # multiply in the largest power of B that stays ≤ plimit
        t = B*B > plimit ? B : B^trunc(Int64, log(B, plimit))
        X, Z = _xladder(t, X, Z, a24, N)
        Z == 0 && return T(1)  # hit the identity mod N: curve failed
    end
    g = gcd(Z, N)
    1 < g < N && return g
    g == N && return T(1)
    return use_stage2 ? _lenstra_stage2(N, X, Z, a24, plimit) : T(1)
end

# Stage 2 (standard continuation): look for a single prime q ∈ (B1, B2] with
# q*Q ≡ identity mod p. Writing q = m*D ± j, that means x(mD*Q) = x(j*Q), i.e.
# p | X_m*Z_j - X_j*Z_m; accumulate those cross products and gcd once at the end.
# Tests every mD ± j with gcd(j, D) = 1 rather than sieving primes out of (B1, B2] —
# a ~3x overhead that avoids enumerating primes up to B2.
# B2 = 25*B1 puts stage-2 cost at ~half of stage 1, which is throughput-optimal for
# this linear-cost continuation: marginal returns decay like 1/(q*ln q) while the
# cost per unit B2 is flat, so the optimum sits below cost balance.
function _lenstra_stage2(N::T, X, Z, a24, B1) where T <: Integer
    B2 = 25 * B1
    D = 210
    # baby steps: j*Q for odd j ≤ D/2, chained by Q_{j+2} = Q_j + 2Q (diff Q_{j-2})
    nb = (D ÷ 2 + 1) ÷ 2
    Xb = Vector{T}(undef, nb)
    Zb = Vector{T}(undef, nb)
    Xb[1], Zb[1] = X, Z
    X2, Z2 = _xdbl(X, Z, a24, N)
    Xp, Zp = X, Z
    Xj, Zj = _xadd(X2, Z2, X, Z, X, Z, N)  # 3Q
    Xb[2], Zb[2] = Xj, Zj
    for idx in 3:nb
        Xn, Zn = _xadd(Xj, Zj, X2, Z2, Xp, Zp, N)
        Xp, Zp = Xj, Zj
        Xj, Zj = Xn, Zn
        Xb[idx], Zb[idx] = Xj, Zj
    end
    # giant steps: R = m*D*Q, stepping by D*Q
    XD, ZD = _xladder(D, X, Z, a24, N)
    m0 = fld(B1, D)
    XRp, ZRp = _xladder(m0 - 1, XD, ZD, a24, N)
    XR, ZR = _xladder(m0, XD, ZD, a24, N)
    acc = T(1)
    for m in m0:cld(B2, D)
        for idx in 1:nb
            j = 2idx - 1
            gcd(j, D) == 1 || continue
            acc = mod(acc*(mod(XR*Zb[idx], N) - mod(Xb[idx]*ZR, N)), N)
        end
        acc == 0 && return T(1)  # N divides the product: factor unrecoverable, next curve
        XRn, ZRn = _xadd(XR, ZR, XD, ZD, XRp, ZRp, N)
        XRp, ZRp = XR, ZR
        XR, ZR = XRn, ZRn
    end
    g = gcd(acc, N)
    return 1 < g < N ? g : T(1)
end

function factor!(n::T, h::AbstractDict{K,Int}) where {T<:Integer,K<:Integer}
    for (p, num_p) in eachfactor(n)
        increment!(h, num_p, p)
    end
    return h
end

"""
    factor(n::Integer) -> Primes.Factorization

Compute the prime factorization of an integer `n`. The returned
object, of type `Factorization`, is an associative container whose
keys correspond to the factors, in sorted order. The value associated
with each key indicates the multiplicity (i.e. the number of times the
factor appears in the factorization).

```julia
julia> factor(100)
2^2 ⋅ 5^2
```

For convenience, a negative number `n` is factored as `-1*(-n)` (i.e. `-1` is considered
to be a factor), and `0` is factored as `0^1`:

```julia
julia> factor(-9)
-1 ⋅ 3^2

julia> factor(0)
0

julia> collect(factor(0))
1-element Array{Pair{Int64,Int64},1}:
 0=>1
```
"""
factor(n::T) where {T<:Integer} = factor!(n, Factorization{T}())


"""
    factor(ContainerType, n::Integer) -> ContainerType

Return the factorization of `n` stored in a `ContainerType`, which must be a
subtype of `AbstractDict` or `AbstractArray`, a `Set`, or an `BitSet`.

```julia
julia> factor(DataStructures.SortedDict, 100)
DataStructures.SortedDict{Int64,Int64,Base.Order.ForwardOrdering} with 2 entries:
  2 => 2
  5 => 2
```

When `ContainerType <: AbstractArray`, this returns the list
of all prime factors of `n` with multiplicities, in sorted order.

```julia
julia> factor(Vector, 100)
4-element Array{Int64,1}:
 2
 2
 5
 5

julia> prod(factor(Vector, 100)) == 100
true
```

When `ContainerType == Set`, this returns the distinct prime
factors as a set.

```julia
julia> factor(Set, 100)
Set([2,5])
```
"""
factor(::Type{D}, n::T) where {T<:Integer, D<:AbstractDict} = factor!(n, D(Dict{T,Int}()))
factor(::Type{A}, n::T) where {T<:Integer, A<:AbstractArray} = A(factor(Vector{T}, n))
factor(::Type{Vector{T}}, n::T) where {T<:Integer} =
    mapreduce(collect, vcat, [repeated(k, v) for (k, v) in factor(n)], init=Vector{T}())
factor(::Type{S}, n::T) where {T<:Integer, S<:Union{Set,BitSet}} = S(keys(factor(n)))
factor(::Type{T}, n) where {T<:Any} = throw(MethodError(factor, (T, n)))

"""
    prodfactors(factors)

Compute `n` (or the radical of `n` when `factors` is of type `Set` or
`BitSet`) where `factors` is interpreted as the result of
`factor(typeof(factors), n)`. Note that if `factors` is of type
`AbstractArray` or `Primes.Factorization`, then `prodfactors` is equivalent
to `Base.prod`.

```jldoctest
julia> prodfactors(factor(100))
100
```
"""
function prodfactors end

prodfactors(factors::AbstractDict{K}) where {K} = isempty(factors) ? one(K) : prod(p^i for (p, i) in factors)
prodfactors(factors::Union{AbstractArray, Set, BitSet}) = prod(factors)

"""
    Base.prod(factors::Primes.Factorization{T}) -> T

Compute `n` where `factors` is interpreted as the result of `factor(n)`.
"""
Base.prod(factors::Factorization) = prodfactors(factors)

"""
    radical(n::Integer)

Compute the radical of `n`, i.e. the largest square-free divisor of `n`.
This is equal to the product of the distinct prime numbers dividing `n`.

```jldoctest
julia> radical(2*2*3)
6
```
"""
radical(n) = n==1 ? one(n) : prod(p for (p, num_p) in eachfactor(n))

"""
    ismersenneprime(M::Integer; [check::Bool = true]) -> Bool

Lucas-Lehmer deterministic test for Mersenne primes. `M` must be a
Mersenne number, i.e. of the form `M = 2^p - 1`, where `p` is a prime
number. Use the keyword argument `check` to enable/disable checking
whether `M` is a valid Mersenne number; to be used with caution.
Return `true` if the given Mersenne number is prime, and `false`
otherwise.

```jldoctest
julia> ismersenneprime(2^11 - 1)
false

julia> ismersenneprime(2^13 - 1)
true
```
"""
function ismersenneprime(M::Integer; check::Bool = true)
    if check
        d = ndigits(M, base=2)
        (M >= 0 && d == count_ones(M) && isprime(d)) ||
            throw(ArgumentError("The argument given is not a valid Mersenne Number (`M = 2^p - 1`)."))
    end
    M < 7 && return M == 3
    return ll_primecheck(M)
end

"""
    isrieselprime(k::Integer, Q::Integer) -> Bool

Lucas-Lehmer-Riesel deterministic test for N of the form `N = k * Q`,
with `0 < k < Q`, `Q = 2^n - 1` and `n > 0`, also known as Riesel primes.
Returns `true` if `R` is prime, and `false` otherwise or
if the combination of k and n is not supported.

```jldoctest
julia> isrieselprime(1, 2^11 - 1)  # == ismersenneprime(2^11 - 1)
false

julia> isrieselprime(3, big(2)^607 - 1)
true
```
"""
function isrieselprime(k::Integer, Q::Integer)
    n = ndigits(Q, base=2)
    0 < k < Q || throw(ArgumentError("The condition 0 < k < Q must be met."))
    if k == 1 && isodd(n)
        return n % 4 == 3 ? ll_primecheck(Q, 3) : ll_primecheck(Q)
    elseif k == 3 && (n % 4) % 3 == 0
        return ll_primecheck(Q, 5778)
    else
        # TODO: Implement a case for (k % 6) % 4 == 1 && ((k % 3) * powermod(2, n, 3)) % 3 < 2
        error("The LLR test is not currently implemented for numbers of this form.")
    end
end

# LL backend -- not for export
function ll_primecheck(X::Integer, s::Integer = 4)
    S, N = BigInt(s), BigInt(ndigits(X, base=2))
    X < 7 && throw(ArgumentError("The condition X ≥ 7 must be met."))
    for i in 1:(N - 2)
        S = (S^2 - 2) % X
    end
    return S == 0
end

"""
    totient(f::Factorization{T}) -> T

Compute the Euler totient function of the number whose prime factorization is
given by `f`. This method may be preferable to [`totient(::Integer)`](@ref)
when the factorization can be reused for other purposes.
"""
function totient(f::Factorization{T}) where T <: Integer
    if iszero(sign(f))
        throw(ArgumentError("ϕ(0) is not defined"))
    end
    result = one(T)
    for (p, k) in f
        result *= p^(k-1) * (p - 1)
    end
    result
end

"""
    totient(n::Integer) -> Integer

Compute the Euler totient function ``ϕ(n)``, which counts the number of
positive integers less than or equal to ``n`` that are relatively prime to
``n`` (that is, the number of positive integers `m ≤ n` with `gcd(m, n) == 1`).
The totient function of `n` when `n` is negative is defined to be
`totient(abs(n))`.
"""
function totient(n::T) where T<:Integer
    n = abs(n)
    if n == 0
        throw(ArgumentError("ϕ(0) is not defined"))
    end
    result = one(T)
    for (p, k) in eachfactor(n)
        result *= p^(k-1) * (p - 1)
    end
    result
end

# add: checked add (when makes sense), result of same type as first argument

add(n::BigInt, x::Int) = n + x
add(n::Integer, x::Int) = Base.checked_add(n, oftype(n, x))
sub(n::BigInt, x::Int) = n - x
sub(n::Integer, x::Int) = Base.checked_sub(n, oftype(n, x))

# add_! : "may" mutate the Integer argument (only for BigInt currently)

# modify a BigInt in-place
function add_!(n::BigInt, x::Int)
    if x < 0
        ccall((:__gmpz_sub_ui, :libgmp), Cvoid, (Ref{BigInt}, Ref{BigInt}, Culong), n, n, -x)
    else
        ccall((:__gmpz_add_ui, :libgmp), Cvoid, (Ref{BigInt}, Ref{BigInt}, Culong), n, n, x)
    end
    n
end

# checked addition, without mutation
add_!(n::Integer, x::Int) = add(n, x)
sub_!(n::BigInt, x::Int) = add_!(n, -x)
sub_!(n::Integer, x::Int) = sub(n, x)

"""
    nextprime(n::Integer, i::Integer=1; interval::Integer=1)

The `i`-th smallest prime not less than `n` (in particular,
`nextprime(p) == p` if `p` is prime). If `i < 0`, this is equivalent to
prevprime(n, -i). Note that for `n::BigInt`, the returned number is
only a pseudo-prime (the function [`isprime`](@ref) is used
internally). See also [`prevprime`](@ref).

If `interval` is provided, primes are sought in increments of `interval`.
This can be useful to ensure the presence of certain divisors in `p-1`.
The range of possible values for `interval` is currently `1:typemax(Int)`.

```jldoctest
julia> nextprime(4)
5

julia> nextprime(5)
5

julia> nextprime(4, 2)
7

julia> nextprime(5, 2)
7

julia> nextprime(2^10+1; interval=1024)
12289

julia> gcd(12289 - 1, 1024) # 1024 | p - 1
1024
```
"""
function nextprime(n::Integer, i::Integer=1; interval::Integer=1)
    i == 0 && throw(DomainError(i))
    i < 0 && return prevprime(n, -i; interval=interval)
    interval < 1 && throw(DomainError(interval, "interval must be >= 1"))
    # TODO: lift the following condition
    interval > typemax(Int) && throw(DomainError(interval, "interval must be <= $(typemax(Int))"))
    interval = oftype(n, interval)
    if n < 2
        n = interval == 1 ?
            oftype(n, 2) :
            # smallest value >= 2 whose difference from n is a multiple of interval
            oftype(n, n + interval * (1 + (n-1)÷(-interval)))
    end
    if n == 2
        if i <= 1
            return n
        else
            n += interval
            i -= 1
        end
    else
        n += iseven(n) ? interval : zero(n)
    end
    # n can now be safely mutated
    # @assert n >= 3
    if iseven(n)
        @assert iseven(interval)
        throw(DomainError((n, interval),
            "`n` and `interval` should not be both even (there is then no correct answer)."))
    end
    # @assert isodd(n)
    interval = Int(interval)
    isodd(interval) && (interval = Base.checked_mul(interval, 2))

    while true
        while !isprime(n)
            n = add_!(n, interval)
        end
        i -= 1
        i <= 0 && break
        n = add_!(n, interval)
    end
    n
end


"""
    prevprime(n::Integer, i::Integer=1; interval::Integer=1)

The `i`-th largest prime not greater than `n` (in particular
`prevprime(p) == p` if `p` is prime). If `i < 0`, this is equivalent to
`nextprime(n, -i)`. Note that for `n::BigInt`, the returned number is
only a pseudo-prime (the function [`isprime`](@ref) is used internally). See
also [`nextprime`](@ref).

If `interval` is provided, primes are sought in increments of `interval`.
This can be useful to ensure the presence of certain divisors in `p-1`.
The range of possible values for `interval` is currently `1:typemax(Int)`.

```jldoctest
julia> prevprime(4)
3

julia> prevprime(5)
5

julia> prevprime(5, 2)
3
```
"""
function prevprime(n::Integer, i::Integer=1; interval::Integer=1)
    i <= 0 && return nextprime(n, -i; interval=interval)
    interval < 1 && throw(DomainError(interval, "interval must be >= 1"))
    interval > typemax(Int) && throw(DomainError(interval, "interval must be <= $(typemax(Int))"))
    interval = Int(interval)

    n += zero(n) # deep copy of n, which is mutated below

    # A bit ugly, but this lets us skip half of the isprime tests when isodd(interval)
    @inline function decrement(n)
        n = sub_!(n, interval)
        iseven(n) && n != 2 ? # n obviously not prime
            sub_!(n, interval) :
            n
    end

    while true
        while !isprime(n)
            n < 2 && throw(ArgumentError("There is no prime less than or equal to $n"))
            n = decrement(n)
        end
        i -= 1
        i <= 0 && break
        n = decrement(n)
    end
    n
end

"""
    prime(::Type{<:Integer}=Int, i::Integer)

The `i`-th prime number.

```jldoctest
julia> prime(1)
2

julia> prime(3)
5

```
"""
prime(::Type{T}, i::Integer) where {T<:Integer} = i < 0 ? throw(DomainError(i)) : nextprime(T(2), i)
prime(i::Integer) = prime(Int, i)


struct NextPrimes{T<:Integer}
    start::T
end

function iterate(np::NextPrimes, state=np.start)
    p = nextprime(state)
    (p, add(p, 1))
end

IteratorSize(::Type{<:NextPrimes}) = Base.IsInfinite()
IteratorEltype(::Type{<:NextPrimes}) = Base.HasEltype()

eltype(::Type{NextPrimes{T}}) where {T} = T

"""
    nextprimes(start::Integer)

Return an iterator over all primes greater than or equal to `start`,
in ascending order.
"""
nextprimes(start::Integer) = NextPrimes(start)

"""
    nextprimes(T::Type=Int)

Return an iterator over all primes, with type `T`.
Equivalent to `nextprimes(T(1))`.
"""
nextprimes(::Type{T}=Int) where {T<:Integer} = nextprimes(one(T))

"""
    nextprimes(start::Integer, n::Integer)

Return an array of the first `n` primes greater than or equal to `start`.

# Example

```
julia> nextprimes(10, 3)
3-element Array{Int64,1}:
 11
 13
 17
```
"""
nextprimes(start::T, n::Integer) where {T<:Integer} =
    collect(T, Iterators.take(nextprimes(start), n))

struct PrevPrimes{T<:Integer}
    start::T
end

function iterate(np::PrevPrimes, state=np.start)
    if isone(state)
        nothing
    else
        p = prevprime(state)
        (p, p-one(p))
    end
end

IteratorSize(::Type{<:PrevPrimes}) = Base.SizeUnknown()
IteratorEltype(::Type{<:PrevPrimes}) = Base.HasEltype()

eltype(::Type{PrevPrimes{T}}) where {T} = T

"""
    prevprimes(start::Integer)

Return an iterator over all primes less than or equal to `start`,
in descending order.

# Example

```
julia> collect(prevprimes(10))
4-element Array{Int64,1}:
 7
 5
 3
 2
```
"""
prevprimes(start::Integer) = PrevPrimes(max(one(start), start))

"""
    prevprimes(start::Integer, n::Integer)

Return an array of the first `n` primes less than or equal to `start`,
in descending order. When there are less than `n` primes less than or
equal to `start`, those primes are returned without an error.

# Example

```
julia> prevprimes(10, 3)
3-element Array{Int64,1}:
 7
 5
 3

julia> prevprimes(10, 10)
4-element Array{Int64,1}:
 7
 5
 3
 2
```
"""
prevprimes(start::T, n::Integer) where {T<:Integer} =
    collect(T, Iterators.take(prevprimes(start), n))

"""
    divisors(n::Integer) -> Vector
Return a vector with the positive divisors of `n`.
For a nonzero integer `n` with prime factorization `n = p₁^k₁ ⋯ pₘ^kₘ`, `divisors(n)`
returns a vector of length (k₁ + 1)⋯(kₘ + 1) containing the divisors of `n` in
lexicographic (rather than numerical) order.
`divisors(-n)` is equivalent to `divisors(n)`.
For convenience, `divisors(0)` returns `[]`.
# Example
```jldoctest; filter = r"(\\s+#.*)?"
julia> divisors(60)
12-element Vector{Int64}:
  1         # 1
  2         # 2
  4         # 2 * 2
  3         # 3
  6         # 3 * 2
 12         # 3 * 2 * 2
  5         # 5
 10         # 5 * 2
 20         # 5 * 2 * 2
 15         # 5 * 3
 30         # 5 * 3 * 2
 60         # 5 * 3 * 2 * 2
julia> divisors(-10)
4-element Vector{Int64}:
  1
  2
  5
 10
julia> divisors(0)
Int64[]
```
"""
function divisors(n::T) where {T<:Integer}
    n = abs(n)
    if iszero(n)
        return T[]
    elseif isone(n)
        return [n]
    else
        return divisors(factor(n))
    end
end

"""
    divisors(f::Factorization) -> Vector
Return a vector with the positive divisors of the number whose factorization is `f`.
Divisors are sorted lexicographically, rather than numerically.
"""
function divisors(f::Factorization{T}) where T
    sgn = sign(f)
    if iszero(sgn) # n == 0
        return T[]
    elseif isempty(f) || length(f) == 1 && sgn < 0 # n == 1 or n == -1
        return [one(T)]
    end

    i = m = 1
    fs = rest(f, 1 + (sgn < 0))
    divs = Vector{T}(undef, prod(x -> x.second + 1, fs))
    divs[i] = one(T)

    for (p, k) in fs
        i = 1
        for _ in 1:k
            for j in i:(i+m-1)
                divs[j+m] = divs[j] * p
            end
            i += m
        end
        m += i - 1
    end

    return divs
end

end # module
