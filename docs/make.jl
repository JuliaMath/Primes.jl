using Documenter, Primes

makedocs(
    modules = [Primes],
    clean = false,
    format = :html,
    sitename = "Primes.jl",
    pages = Any[
        "Home" => "index.md",
        "Functions" => "api.md"
    ],
)

deploydocs(
    julia = "nightly",
    repo = "github.com/JuliaMath/Primes.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
