# implementation of a sorted dict (not optimized for speed) for storing
# the factorization of an integer

struct Factorization{T<:Integer} <: AbstractDict{T, Int}
    pe::Vector{Pair{T, Int}} # Prime-Exponent

    function Factorization{T}() where {T<:Integer}
        # preallocates enough space that numbers smaller than 2310 won't need to resize
        v = Vector{Pair{T, Int}}(undef, 4)
        empty!(v)
        new{T}(v)
    end
end

function Factorization{T}(d::AbstractDict) where T<:Integer
    f = Factorization{T}()
    append!(f.pe, sort!(collect(d)))
    f
end

Factorization(d::AbstractDict{T}) where {T<:Integer} = Factorization{T}(d)
Base.convert(::Type{Factorization}, d::AbstractDict) = Factorization(d)

Base.iterate(f::Factorization, state...) = iterate(f.pe, state...)

function Base.get(f::Factorization, p, default)
    found = searchsortedfirst(f.pe, p, by=first)
    (found > length(f.pe) || first(f.pe[found])) != p  ? default : last(f.pe[found])
end

Base.getindex(f::Factorization, p::Integer) = get(f, p, 0)

function Base.setindex!(f::Factorization{T}, e::Int, p::Integer) where T
    found = searchsortedfirst(f.pe, p, by=first)
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
function increment!(f::Factorization{T}, e::Int, p::Integer) where T
    found = searchsortedfirst(f.pe, p, by=first)
    if found > length(f.pe)
        push!(f.pe, T(p)=>e)
    elseif first(f.pe[found]) != p
        insert!(f.pe, found, T(p)=>e)
    else
        f.pe[found] = T(p)=>(last(f.pe[found])+e)
    end
    f
end
increment!(f::AbstractDict, e::Int, p::Integer) = (h[p] = get(h, p, 0) + 1)

Base.length(f::Factorization) = length(f.pe)

Base.show(io::IO, ::MIME{Symbol("text/plain")}, f::Factorization) =
    join(io, isempty(f) ? "1" : [(e == 1 ? "$p" : "$p^$e") for (p,e) in f.pe], " * ")
