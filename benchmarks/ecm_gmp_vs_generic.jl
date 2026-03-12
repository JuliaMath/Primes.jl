# ECM Benchmark: GMP In-Place vs Generic (Allocating) BigInt Path
#
# Compares the two ECM code paths for BigInt inputs:
#   1. BigInt-specialised: uses in-place GMP arithmetic (zero-allocation hot loop)
#   2. Generic: functional style, allocates new BigInts per operation
#
# Usage:
#   julia --project=benchmarks benchmarks/ecm_gmp_vs_generic.jl

using Primes
using BenchmarkTools

# Reproduce the generic (allocating) ECM path for BigInt, bypassing the
# BigInt-specialised method that normally intercepts dispatch.
function ecm_factor_generic(n::BigInt, B1::Int, num_curves::Int)
    prime_powers = Primes._ecm_prime_powers(B1)
    T = BigInt
    for _ in 1:num_curves
        curve = Primes._ecm_suyama(n)
        curve === nothing && continue
        curve isa Tuple || return curve
        x0, z0, a24 = curve

        QX, QZ = x0, z0
        degenerate = false
        acc = one(T)
        batch_count = 0
        for pk in prime_powers
            QX, QZ = Primes._ecm_scalar_mul(pk, QX, QZ, n, a24)
            acc = Primes._mulmod(acc, QZ, n)
            batch_count += 1
            if batch_count >= 100
                g = gcd(acc, n)
                if 1 < g < n
                    return g
                end
                if g == n
                    degenerate = true
                    break
                end
                acc = one(T)
                batch_count = 0
            end
        end
        degenerate && continue
        if batch_count > 0
            g = gcd(acc, n)
            1 < g < n && return g
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Test cases at different scales
# ---------------------------------------------------------------------------
const CASES = [
    (
        name  = "small (~40-digit semiprime)",
        n     = big"824633720831" * big"1000000007",
        B1    = 2_000,
        curves = 50,
    ),
    (
        name  = "medium (~55-digit semiprime)",
        n     = big"780002082420426809" * big"810735269523504809437013569",
        B1    = 11_000,
        curves = 200,
    ),
]

# ---------------------------------------------------------------------------
# Run benchmarks
# ---------------------------------------------------------------------------
println("=" ^ 70)
println("  ECM Benchmark: GMP In-Place vs Generic (Allocating) BigInt")
println("=" ^ 70)

for case in CASES
    (; name, n, B1, curves) = case
    println("\n--- $name ---")
    println("  n      = $n  ($(ndigits(n)) digits)")
    println("  B1     = $B1")
    println("  curves = $curves")

    # Warm-up both paths
    Primes.ecm_factor(n, B1, curves)
    ecm_factor_generic(n, B1, curves)

    println("\n  [GMP in-place]")
    b_gmp = @benchmark Primes.ecm_factor($n, $B1, $curves) samples=20 evals=1
    display(b_gmp)

    println("\n  [Generic (allocating)]")
    b_gen = @benchmark ecm_factor_generic($n, $B1, $curves) samples=20 evals=1
    display(b_gen)

    med_gmp = median(b_gmp).time / 1e6
    med_gen = median(b_gen).time / 1e6
    alloc_gmp = median(b_gmp).memory / 1024
    alloc_gen = median(b_gen).memory / 1024

    println("\n  Summary (median):")
    println("    GMP in-place : $(round(med_gmp; digits=2)) ms,  $(round(alloc_gmp; digits=1)) KiB")
    println("    Generic alloc: $(round(med_gen; digits=2)) ms,  $(round(alloc_gen; digits=1)) KiB")
    if med_gmp > 0
        println("    Speedup      : $(round(med_gen / med_gmp; digits=2))×")
    end
    if alloc_gen > 0
        println("    Memory saved : $(round((1 - alloc_gmp / alloc_gen) * 100; digits=1))%")
    end
end

println("\n" * "=" ^ 70)
println("  Done.")
println("=" ^ 70)
