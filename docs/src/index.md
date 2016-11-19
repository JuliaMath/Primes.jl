# Primes.jl

This package provides functions for computing prime numbers in Julia.

## Installation

The package is available for Julia versions 0.4 and up.
To install it, run
```julia
Pkg.add("Primes")
```
from the Julia REPL.

## Note

Prior to Julia 0.5, these (or similar) functions were available in Julia's Base module.
Because of this, the symbols from this package are not exported on Julia 0.4 to avoid
naming conflicts.
In this case, the symbols will need to be explicitly imported or called with the prefix
`Primes`.
For example,
```julia
using Primes
import Primes: isprime, primes, primesmask, factor
```
This is not necessary for Julia versions 0.5 or later.
