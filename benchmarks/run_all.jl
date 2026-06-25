# Run All Benchmarks
#
# Usage:
#   julia --project=benchmarks benchmarks/run_all.jl
#
# Or run individual benchmarks:
#   julia --project=benchmarks benchmarks/ecm_gmp_vs_generic.jl
#   julia --project=benchmarks benchmarks/ecm_microbenchmarks.jl
#   julia --project=benchmarks benchmarks/factorization_endtoend.jl

const BENCHDIR = @__DIR__

const SCRIPTS = [
    "ecm_microbenchmarks.jl",
    "ecm_gmp_vs_generic.jl",
    "factorization_endtoend.jl",
]

for script in SCRIPTS
    path = joinpath(BENCHDIR, script)
    println("\n" * "#" ^ 70)
    println("# Running: $script")
    println("#" ^ 70 * "\n")
    include(path)
end
