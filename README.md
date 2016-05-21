# Primes.jl

[![Build Status](https://travis-ci.org/JuliaMath/Primes.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Primes.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaMath/Primes.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/Primes.jl?branch=master)

Julia functions for computing prime numbers.

At the moment, this repository just contains the following functions which have been duplicated from Base Julia:
* `isprime`
* `primes`
* `primesmask`
* `factor`

None of these are currently exported, so need to be explicitly imported before use.
