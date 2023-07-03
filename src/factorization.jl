# implementation of a sorted dict (not optimized for speed) for storing
# the factorization of an integer

struct Factorization{T} <: AbstractDict{T, Int}
    pe::Vector{Pair{T, Int}} # Prime-Exponent

    function Factorization{T}() where T
        # preallocates enough space that numbers smaller than 2310 won't need to resize
        v = Vector{Pair{T, Int}}(undef, 4)
        empty!(v)
        new{T}(v)
    end
end

function Factorization{T}(d::AbstractDict) where T
    f = Factorization{T}()
    append!(f.pe, sort!(collect(d)))
    f
end

Factorization(d::AbstractDict{T}) where T = Factorization{T}(d)
Base.convert(::Type{Factorization}, d::AbstractDict) = Factorization(d)

Base.iterate(f::Factorization, state...) = iterate(f.pe, state...)

function Base.get(f::Factorization, p, default)
    found = searchsortedfirst(f.pe, p=>0, by=first)
    (found > length(f.pe) || first(f.pe[found])) != p  ? default : last(f.pe[found])
end

Base.getindex(f::Factorization, p) = get(f, p, 0)

function Base.setindex!(f::Factorization{T}, e::Int, p) where T
    found = searchsortedfirst(f.pe, p=>0, by=first)
    if found > length(f.pe)
        push!(f.pe, T(p)=>e)
    elseif first(f.pe[found]) != p
        insert!(f.pe, found, T(p)=>e)
    else
        f.pe[found] = T(p)=>e
    end
    f
end

"""
    impliments f[p] += e faster
"""
function increment!(f::Factorization{T}, e::Int, p) where T
    found = searchsortedfirst(f.pe, p=>0, by=first)
    if found > length(f.pe)
        push!(f.pe, T(p)=>e)
    elseif first(f.pe[found]) != p
        insert!(f.pe, found, T(p)=>e)
    else
        f.pe[found] = T(p)=>(last(f.pe[found])+e)
    end
    f
end
function increment!(f::AbstractDict, e::Int, p)
    f[p] = get(f, p, 0) + e
    return f
end

Base.length(f::Factorization) = length(f.pe)

Base.show(io::IO, ::MIME{Symbol("text/plain")}, f::Factorization) =
    join(io, isempty(f) ? "1" : [(e == 1 ? "$p" : "$p^$e") for (p,e) in f.pe], " * ")

Base.sign(f::Factorization) = isempty(f.pe) ? one(keytype(f)) : sign(first(f.pe[1]))
