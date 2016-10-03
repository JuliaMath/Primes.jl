using Documenter, Primes

makedocs(
    modules = [Primes],
    format = Documenter.Formats.HTML,
    sitename = "Primes.jl",
    pages = Any[
        "Home" => "index.md",
        "Functions" => "api.md"
    ],
    doctest = false,
)

deploydocs(
    repo = "github.com/JuliaMath/Primes.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
