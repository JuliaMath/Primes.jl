# This includes parts that were formerly a part of Julia. License is MIT: http://julialang.org/license
VERSION < v"0.7.0-beta2.199" && __precompile__()
module Primes

using Base.Iterators: repeated

using Base: BitSigned
using Base.Checked: checked_neg

export isprime, primes, primesmask, factor, ismersenneprime, isrieselprime,
       nextprime, prevprime, prime, prodfactors, radical, totient

include("factorization.jl")

# Primes generating functions
#     https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
#     https://en.wikipedia.org/wiki/Wheel_factorization
#     http://primesieve.org
#     Jonathan Sorenson, "An analysis of two prime number sieves", Computer Science Technical Report Vol. 1028, 1991
const wheel         = [4,  2,  4,  2,  4,  6,  2,  6]
const wheel_primes  = [7, 11, 13, 17, 19, 23, 29, 31]
const wheel_indices = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7]

@inline function wheel_index(n)
    d, r = divrem(n - 1, 30)
    return 8d + wheel_indices[r + 2]
end
@inline function wheel_prime(n)
    d, r = (n - 1) >>> 3, (n - 1) & 7
    return 30d + wheel_primes[r + 1]
end

function _primesmask(limit::Int)
    limit < 7 && throw(ArgumentError("The condition limit ≥ 7 must be met."))
    n = wheel_index(limit)
    m = wheel_prime(n)
    sieve = ones(Bool, n)
    @inbounds for i = 1:wheel_index(isqrt(limit))
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
            q = p * wheel_prime(j + 1)
            j = j & 7 + 1
            while q ≤ m
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
    isprime(n::Integer) -> Bool

Returns `true` if `n` is prime, and `false` otherwise.

```julia
julia> isprime(3)
true
```
"""
function isprime(n::Integer)
    # Small precomputed primes + Miller-Rabin for primality testing:
    #     https://en.wikipedia.org/wiki/Miller–Rabin_primality_test
    #     https://github.com/JuliaLang/julia/issues/11594
    for m in (2, 3, 5, 7, 11, 13, 17, 19, 23)
        n % m == 0 && return n == m
    end
    n < 841 && return n > 1
    s = trailing_zeros(n - 1)
    d = (n - 1) >>> s
    for a in witnesses(n)::Tuple{Vararg{Int}}
        x = powermod(a, d, n)
        x == 1 && continue
        t = s
        while x != n - 1
            (t -= 1) ≤ 0 && return false
            x = oftype(n, widemul(x, x) % n)
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

```julia
julia> isprime(big(3))
true
```
"""
isprime(x::BigInt, reps=25) = ccall((:__gmpz_probab_prime_p,:libgmp),
                                    Cint, (Any, Cint), x, reps) > 0
# TODO: Change `Any` to `Ref{BigInt}` when 0.6 support is dropped.
# The two have the same effect but `Ref{BigInt}` won't be optimized on 0.6.


# Miller-Rabin witness choices based on:
#     http://mathoverflow.net/questions/101922/smallest-collection-of-bases-for-prime-testing-of-64-bit-numbers
#     http://primes.utm.edu/prove/merged.html
#     http://miller-rabin.appspot.com
#     https://github.com/JuliaLang/julia/issues/11594
#     Forišek and Jančina, "Fast Primality Testing for Integers That Fit into a Machine Word", 2015
#         (in particular, see function FJ32_256, from which the hash and bases were taken)
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
    @inbounds return (Int(bases[i]),)
end
witnesses(n::Integer) =
        n < 4294967296      ? _witnesses(UInt64(n)) :
        n < 2152302898747   ? (2, 3, 5, 7, 11) :
        n < 3474749660383   ? (2, 3, 5, 7, 11, 13) :
                              (2, 325, 9375, 28178, 450775, 9780504, 1795265022)

isprime(n::UInt128) =
    n ≤ typemax(UInt64) ? isprime(UInt64(n)) : isprime(BigInt(n))
isprime(n::Int128) = n < 2 ? false :
    n ≤ typemax(Int64)  ? isprime(Int64(n))  : isprime(BigInt(n))


# Trial division of small (< 2^16) precomputed primes +
# Pollard rho's algorithm with Richard P. Brent optimizations
#     https://en.wikipedia.org/wiki/Trial_division
#     https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
#     http://maths-people.anu.edu.au/~brent/pub/pub051.html
#
function factor!(n::T, h::AbstractDict{K,Int}) where {T<:Integer,K<:Integer}
    # check for special cases
    if n < 0
        h[-1] = 1
        if isa(n, BitSigned) && n == typemin(T)
            h[2] = 8 * sizeof(T) - 1
            return h
        else
            return factor!(checked_neg(n), h)
        end
    elseif n == 1
        return h
    elseif n == 0 || isprime(n)
        h[n] = 1
        return h
    end

    local p::T
    for p in PRIMES
        if n % p == 0
            h[p] = get(h, p, 0) + 1
            n = div(n, p)
            while n % p == 0
                h[p] = get(h, p, 0) + 1
                n = div(n, p)
            end
            n == 1 && return h
            isprime(n) && (h[n] = 1; return h)
        end
    end
    T <: BigInt || widemul(n - 1, n - 1) ≤ typemax(n) ? pollardfactors!(n, h) : pollardfactors!(widen(n), h)
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
if VERSION >= v"0.7.0-beta.81"
    factor(::Type{Vector{T}}, n::T) where {T<:Integer} =
        mapreduce(collect, vcat, [repeated(k, v) for (k, v) in factor(n)], init=Vector{T}())
else
    factor(::Type{Vector{T}}, n::T) where {T<:Integer} =
        mapreduce(collect, vcat, Vector{T}(), [repeated(k, v) for (k, v) in factor(n)])
end
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
radical(n) = prod(factor(Set, n))

function pollardfactors!(n::T, h::AbstractDict{K,Int}) where {T<:Integer,K<:Integer}
    while true
        c::T = rand(1:(n - 1))
        G::T = 1
        r::K = 1
        y::T = rand(0:(n - 1))
        m::K = 1900
        ys::T = 0
        q::T = 1
        x::T = 0
        while c == n - 2
            c = rand(1:(n - 1))
        end
        while G == 1
            x = y
            for i in 1:r
                y = y^2 % n
                y = (y + c) % n
            end
            k::K = 0
            G = 1
            while k < r && G == 1
                for i in 1:(m > (r - k) ? (r - k) : m)
                    ys = y
                    y = y^2 % n
                    y = (y + c) % n
                    q = (q * (x > y ? x - y : y - x)) % n
                end
                G = gcd(q, n)
                k += m
            end
            r *= 2
        end
        G == n && (G = 1)
        while G == 1
            ys = ys^2 % n
            ys = (ys + c) % n
            G = gcd(x > ys ? x - ys : ys - x, n)
        end
        if G != n
            isprime(G) ? h[G] = get(h, G, 0) + 1 : pollardfactors!(G, h)
            G2 = div(n,G)
            isprime(G2) ? h[G2] = get(h, G2, 0) + 1 : pollardfactors!(G2, h)
            return h
        end
    end
end

if VERSION >= v"0.7.0-beta.183"
    ndigits2(n) = ndigits(n, base=2)
else
    ndigits2(n) = ndigits(n, 2)
end

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
        d = ndigits2(M)
        M >= 0 && isprime(d) && (M >> d == 0) ||
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

julia> isrieselprime(3, 2^607 - 1)
true
```
"""
function isrieselprime(k::Integer, Q::Integer)
    n = ndigits2(Q)
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
    S, N = BigInt(s), BigInt(ndigits2(X))
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
    if !isempty(f) && first(first(f)) == 0
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
function totient(n::Integer)
    totient(factor(abs(n)))
end

# add_! : "may" mutate the Integer argument (only for BigInt currently)

# modify a BigInt in-place
function add_!(n::BigInt, x::Int)
    # TODO: Change `Any` to `Ref{BigInt}` when 0.6 support is dropped.
    # The two have the same effect but `Ref{BigInt}` won't be optimized on 0.6.
    if x < 0
        ccall((:__gmpz_sub_ui, :libgmp), Cvoid, (Any, Any, Culong), n, n, -x)
    else
        ccall((:__gmpz_add_ui, :libgmp), Cvoid, (Any, Any, Culong), n, n, x)
    end
    n
end

# checked addition, without mutation
add_!(n::Integer, x::Int) = Base.checked_add(n, oftype(n, x))

"""
    nextprime(n::Integer, i::Integer=1)

The `i`-th smallest prime not less than `n` (in particular,
`nextprime(p) == p` if `p` is prime). If `i < 0`, this is equivalent to
prevprime(n, -i). Note that for `n::BigInt`, the returned number is
only a pseudo-prime (the function [`isprime`](@ref) is used
internally). See also [`prevprime`](@ref).

```jldoctest
julia> nextprime(4)
5

julia> nextprime(5)
5

julia> nextprime(4, 2)
7

julia> nextprime(5, 2)
7
```
"""
function nextprime(n::Integer, i::Integer=1)
    i < 0 && return prevprime(n, -i)
    i == 0 && throw(DomainError(i))
    n < 2 && (n = oftype(n, 2))
    if n == 2
        if i <= 1
            return n
        else
            n += one(n)
            i -= 1
        end
    else
        n += iseven(n)
    end
    # n can now be safely mutated
    # @assert isodd(n) && n >= 3
    while true
        while !isprime(n)
            n = add_!(n, 2)
        end
        i -= 1
        i <= 0 && break
        n = add_!(n, 2)
    end
    n
end


"""
    prevprime(n::Integer, i::Integer=1)

The `i`-th largest prime not greater than `n` (in particular
`prevprime(p) == p` if `p` is prime). If `i < 0`, this is equivalent to
`nextprime(n, -i)`. Note that for `n::BigInt`, the returned number is
only a pseudo-prime (the function [`isprime`](@ref) is used internally). See
also [`nextprime`](@ref).

```jldoctest
julia> prevprime(4)
3

julia> prevprime(5)
5

julia> prevprime(5, 2)
3
```
"""
function prevprime(n::Integer, i::Integer=1)
    i <= 0 && return nextprime(n, -i)
    n += zero(n) # deep copy of n, which is mutated below
    while true
        n < 2 && throw(ArgumentError("There is no prime less than or equal to $n"))
        while !isprime(n)
            n = add_!(n, -1)
        end
        i -= 1
        i <= 0 && break
        n = add_!(n, -1)
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

end # module
