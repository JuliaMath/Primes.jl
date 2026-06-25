# Benchmarks

Benchmark suite for Primes.jl, comparing GMP-optimised in-place arithmetic against
the generic (allocating) code path for BigInt operations, and measuring end-to-end
factorisation performance across different number sizes.

## Prerequisites

```bash
julia --project=benchmarks -e 'using Pkg; Pkg.instantiate()'
```

## Running

Run all benchmarks:

```bash
julia --project=benchmarks benchmarks/run_all.jl
```

Or run individual suites:

```bash
# GMP in-place vs generic (allocating) ECM — the core comparison
julia --project=benchmarks benchmarks/ecm_gmp_vs_generic.jl

# Micro-benchmarks for individual ECM operations
julia --project=benchmarks benchmarks/ecm_microbenchmarks.jl

# End-to-end factor(n) at various sizes
julia --project=benchmarks benchmarks/factorization_endtoend.jl
```

## Benchmark Suites

### `ecm_gmp_vs_generic.jl`

Compares the two ECM code paths for `BigInt` inputs:

| Path | Style | Allocations |
|------|-------|-------------|
| `ecm_factor(n::BigInt, ...)` | In-place GMP (`Base.GMP.MPZ.*!`) | Near-zero in hot loop |
| Generic via `_ecm_scalar_mul` | Functional (new BigInt per op) | O(k) per scalar multiply |

Tests at two scales: a ~40-digit and a ~55-digit semiprime.

### `ecm_microbenchmarks.jl`

Benchmarks individual ECM building blocks across three execution modes:

| Operation | GMP in-place | Generic BigInt | UInt128 (LLVM) |
|-----------|-------------|----------------|----------------|
| `_mulmod` | `_mulmod!` | `_mulmod` | `_mulmod` |
| `_ecm_add` | `_ecm_add!` | `_ecm_add` | `_ecm_add` |
| `_ecm_double` | `_ecm_double!` | `_ecm_double` | `_ecm_double` |
| `_ecm_scalar_mul` | `_ecm_scalar_mul!` | `_ecm_scalar_mul` | `_ecm_scalar_mul` |

### `factorization_endtoend.jl`

Benchmarks `factor(n)` for semiprimes ranging from 12 to 45+ digits,
showing how the polyalgorithm (trial division → Pollard rho → ECM → MPQS) scales.
