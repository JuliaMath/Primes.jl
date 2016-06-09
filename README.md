# Primes.jl

[![Primes](http://pkg.julialang.org/badges/Primes_0.4.svg)](http://pkg.julialang.org/?pkg=Primes)
[![Primes](http://pkg.julialang.org/badges/Primes_0.5.svg)](http://pkg.julialang.org/?pkg=Primes)
[![Build Status](https://travis-ci.org/JuliaMath/Primes.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Primes.jl)
[![Windows Build](https://ci.appveyor.com/api/projects/status/ao64pk44lwo0092r/branch/master?svg=true)](https://ci.appveyor.com/project/ararslan/primes-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/JuliaMath/Primes.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/Primes.jl?branch=master)

Julia functions for computing prime numbers.

## Usage

At the moment, this repository just contains the following functions which have been duplicated from Base Julia:

    factor(n) -> Dict

> Compute the prime factorization of an integer `n`. Returns a dictionary. The keys of the
dictionary correspond to the factors, and hence are of the same type as `n`. The value
associated with each key indicates the number of times the factor appears in the
factorization.
> ```julia
julia> factor(100) # == 2*2*5*5
Dict{Int64,Int64} with 2 entries:
  2 => 2
  5 => 2
> ```

    factorvec(n::Integer) -> Vector

> Compute the prime factorization of `n` with multiplicities. Returns a vector with
the same type as `n`. The product of the returned vector will equal `n`.

> ```julia
julia> factorvec(100)
4-element Array{Int64,1}:
 2
 2
 5
 5
> ```

    isprime(x::Integer) -> Bool

> Returns `true` if `x` is prime, and `false` otherwise.
> ```julia
julia> isprime(3)
true
> ```

    isprime(x::BigInt, [reps = 25]) -> Bool

> Probabilistic primality test. Returns `true` if `x` is prime with high probability (pseudoprime);
and `false` if `x` is composite (not prime). The false positive rate is about `0.25^reps`.
`reps = 25` is considered safe for cryptographic applications (Knuth, Seminumerical Algorithms).

> ```julia
julia> isprime(big(3))
true
> ```

    primes([lo,] hi)

> Returns a collection of the prime numbers (from `lo`, if specified) up to `hi`.

    primesmask([lo,] hi)

> Returns a prime sieve, as a `BitArray`, of the positive integers (from `lo`, if specified)
up to `hi`. Useful when working with either primes or composite numbers.

To avoid naming conflicts with Base, these are not exported for Julia version 0.4. In this case you will need to explicitly import the symbols:

    using Primes
    import Primes: isprime, primes, primesmask, factor
