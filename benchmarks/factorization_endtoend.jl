# End-to-End Factorization Benchmarks
#
# Benchmarks `factor(n)` for semiprimes of increasing size to show the
# polyalgorithm's performance across different regimes:
#   - Small (< 20 digits): trial division + Pollard rho
#   - Medium (20-40 digits): ECM kicks in
#   - Large (40-60 digits): ECM + MPQS
#   - Very large (60+ digits): MPQS-dominated
#
# Usage:
#   julia --project=benchmarks benchmarks/factorization_endtoend.jl

using Primes
using BenchmarkTools

const CASES = [
    (name = "~12 digits",  n = big"824633720831"),
    (name = "~20 digits",  n = big"99999999999999999877" ),  # 20-digit semiprime-like
    (name = "~30 digits",  n = big"824633720831" * big"1000000007"),
    (name = "~40 digits",  n = big"2152302898747" * big"99999999999999999877"),
    (name = "~45 digits",  n = big"780002082420426809" * big"810735269523504809437013569"),
]

println("=" ^ 70)
println("  End-to-End factor(n) Benchmarks")
println("=" ^ 70)

for case in CASES
    (; name, n) = case
    println("\n--- $name  ($(ndigits(n)) digits) ---")
    println("  n = $n")

    # Warm-up and verify correctness
    f = factor(n)
    @assert prodfactors(f) == n "Factorization verification failed for $name"
    println("  factor(n) = $f")

    # Benchmark
    samples = ndigits(n) > 40 ? 5 : 20
    b = @benchmark factor($n) samples=samples evals=1
    display(b)

    med = median(b).time
    unit = med < 1e6 ? ("$(round(med / 1e3; digits=2)) μs") :
           med < 1e9 ? ("$(round(med / 1e6; digits=2)) ms") :
                       ("$(round(med / 1e9; digits=2)) s")
    println("  Median: $unit, $(round(median(b).memory / 1024; digits=1)) KiB allocated")
end

println("\n" * "=" ^ 70)
println("  Done.")
println("=" ^ 70)
