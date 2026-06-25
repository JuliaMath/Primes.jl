# Benchmarks

Primes.jl ships a benchmark suite in `benchmarks/` for measuring the performance of
its factorisation algorithms, and comparing the GMP-optimised BigInt code path against
the generic (allocating) implementation.

## Running Benchmarks

Install dependencies once:

```bash
julia --project=benchmarks -e 'using Pkg; Pkg.instantiate()'
```

Run all benchmarks:

```bash
julia --project=benchmarks benchmarks/run_all.jl
```

Or run individual suites:

```bash
julia --project=benchmarks benchmarks/ecm_gmp_vs_generic.jl
julia --project=benchmarks benchmarks/ecm_microbenchmarks.jl
julia --project=benchmarks benchmarks/factorization_endtoend.jl
```

## Benchmark Suites

### ECM: GMP In-Place vs Generic

`ecm_gmp_vs_generic.jl` compares the two ECM code paths for `BigInt` inputs:

- **GMP in-place** (`ecm_factor(n::BigInt, ...)`): Uses `Base.GMP.MPZ.*!` functions for
  zero-allocation arithmetic in the Montgomery ladder hot loop.
- **Generic** (via `_ecm_scalar_mul`): Functional style that allocates a new `BigInt` for
  every modular operation.

### ECM Micro-Benchmarks

`ecm_microbenchmarks.jl` measures individual ECM building blocks across four execution
modes:

| Operation | GMP in-place | Generic BigInt | UInt128 | UInt256 (LLVM) |
|-----------|-------------|----------------|---------|----------------|
| Modular multiply | `_mulmod!` | `_mulmod` | `_mulmod` | `_mulmod` |
| Point addition | `_ecm_add!` | `_ecm_add` | `_ecm_add` | `_ecm_add` |
| Point doubling | `_ecm_double!` | `_ecm_double` | `_ecm_double` | `_ecm_double` |
| Scalar multiply | `_ecm_scalar_mul!` | `_ecm_scalar_mul` | `_ecm_scalar_mul` | `_ecm_scalar_mul` |

!!! note
    `widen(UInt128) == BigInt` in base Julia, so UInt128 still allocates through BigInt.
    `widen(UInt256) == UInt512` via BitIntegers.jl, giving pure LLVM fixed-width arithmetic
    with zero allocations — but GMP's hand-tuned assembly is faster per-operation.

### End-to-End Factorisation

`factorization_endtoend.jl` benchmarks `factor(n)` for semiprimes ranging from 12 to 45+
digits, showing how the polyalgorithm (trial division → Pollard rho → ECM → MPQS) scales
with input size.
