using Documenter, Primes

DocMeta.setdocmeta!(Primes, :DocTestSetup, :(using Primes); recursive = true)
makedocs(
    modules = [Primes],
    sitename = "Primes.jl",
    checkdocs = :none,
    pages = Any[
        "Home" => "index.md",
        "Functions" => "api.md",
        "Benchmarks" => "benchmarks.md"
    ],
)

deploydocs(
    repo = "github.com/JuliaMath/Primes.jl.git",
    target = "build",
)
