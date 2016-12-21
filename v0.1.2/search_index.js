var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Primes.jl-1",
    "page": "Home",
    "title": "Primes.jl",
    "category": "section",
    "text": "This package provides functions for computing prime numbers in Julia."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is available for Julia versions 0.4 and up. To install it, runPkg.add(\"Primes\")from the Julia REPL."
},

{
    "location": "index.html#Note-1",
    "page": "Home",
    "title": "Note",
    "category": "section",
    "text": "Prior to Julia 0.5, these (or similar) functions were available in Julia's Base module. Because of this, the symbols from this package are not exported on Julia 0.4 to avoid naming conflicts. In this case, the symbols will need to be explicitly imported or called with the prefix Primes. For example,using Primes\nimport Primes: isprime, primes, primesmask, factorThis is not necessary for Julia versions 0.5 or later."
},

{
    "location": "api.html#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#Prime-number-functions-1",
    "page": "Functions",
    "title": "Prime number functions",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#Base.factor",
    "page": "Functions",
    "title": "Base.factor",
    "category": "Function",
    "text": "factor(n::Integer) -> Dict\n\nCompute the prime factorization of an integer n. Returns a dictionary. The keys of the dictionary correspond to the factors, and hence are of the same type as n. The value associated with each key indicates the number of times the factor appears in the factorization.\n\njulia> factor(100) # == 2^2 * 5^2\nDict{Int64,Int64} with 2 entries:\n  2 => 2\n  5 => 2\n\nFor convenience, a negative number n is factored as -1*(-n) (i.e. -1 is considered to be a factor), and 0 is factored as 0^1:\n\njulia> factor(-9) # == -1 * 3^2\nDict{Int64,Int64} with 2 entries:\n  -1 => 1\n   3 => 2\n\njulia> factor(0)\nDict{Int64,Int64} with 1 entries:\n  0 => 1\n\n\n\nfactor(ContainerType, n::Integer) -> ContainerType\n\nReturn the factorization of n stored in a ContainerType, which must be a subtype of Associative or AbstractArray, a Set, or an IntSet.\n\njulia> factor(DataStructures.SortedDict, 100)\nDataStructures.SortedDict{Int64,Int64,Base.Order.ForwardOrdering} with 2 entries:\n  2 => 2\n  5 => 2\n\nWhen ContainerType <: AbstractArray, this returns the list of all prime factors of n with multiplicities, in sorted order.\n\njulia> factor(Vector, 100)\n4-element Array{Int64,1}:\n 2\n 2\n 5\n 5\n\njulia> prod(factor(Vector, 100)) == 100\ntrue\n\nWhen ContainerType == Set, this returns the distinct prime factors as a set.\n\njulia> factor(Set, 100)\nSet([2,5])\n\n\n\n"
},

{
    "location": "api.html#Prime-factorization-1",
    "page": "Functions",
    "title": "Prime factorization",
    "category": "section",
    "text": "Primes.factor"
},

{
    "location": "api.html#Base.primes",
    "page": "Functions",
    "title": "Base.primes",
    "category": "Function",
    "text": "primes([lo,] hi)\n\nReturns a collection of the prime numbers (from lo, if specified) up to hi.\n\n\n\n"
},

{
    "location": "api.html#Generating-prime-numbers-1",
    "page": "Functions",
    "title": "Generating prime numbers",
    "category": "section",
    "text": "Primes.primes"
},

{
    "location": "api.html#Base.isprime",
    "page": "Functions",
    "title": "Base.isprime",
    "category": "Function",
    "text": "isprime(x::Integer) -> Bool\n\nReturns true if x is prime, and false otherwise.\n\njulia> isprime(3)\ntrue\n\n\n\nisprime(x::BigInt, [reps = 25]) -> Bool\n\nProbabilistic primality test. Returns true if x is prime with high probability (pseudoprime); and false if x is composite (not prime). The false positive rate is about 0.25^reps. reps = 25 is considered safe for cryptographic applications (Knuth, Seminumerical Algorithms).\n\njulia> isprime(big(3))\ntrue\n\n\n\n"
},

{
    "location": "api.html#Base.primesmask",
    "page": "Functions",
    "title": "Base.primesmask",
    "category": "Function",
    "text": "primesmask([lo,] hi)\n\nReturns a prime sieve, as a BitArray, of the positive integers (from lo, if specified) up to hi. Useful when working with either primes or composite numbers.\n\n\n\n"
},

{
    "location": "api.html#Identifying-prime-numbers-1",
    "page": "Functions",
    "title": "Identifying prime numbers",
    "category": "section",
    "text": "Primes.isprime\nPrimes.primesmask"
},

]}
