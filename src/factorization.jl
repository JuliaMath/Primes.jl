# implementation of a sorted dict (not optimized for speed) for storing
# the factorization of an integer

struct Factorization{T<:Integer} <: AbstractDict{T, Int}
    pe::Vector{Pair{T, Int}} # Prime-Exponent

    Factorization{T}() where {T<:Integer} = new{T}(Vector{Pair{T, Int}}())
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
    found = searchsorted(f.pe, p, by=first)
    isempty(found) ?
        default :
        last(f.pe[first(found)])
end

Base.getindex(f::Factorization, p::Integer) = get(f, p, 0)

function Base.setindex!(f::Factorization{T}, e::Int, p::Integer) where T
    found = searchsorted(f.pe, p, by=first)
    if isempty(found)
        insert!(f.pe, first(found), T(p)=>e)
    else
        f.pe[first(found)] = T(p)=>e
    end
    f
end

Base.length(f::Factorization) = length(f.pe)

Base.show(io::IO, ::MIME{Symbol("text/plain")}, f::Factorization) =
    join(io, isempty(f) ? "1" : [(e == 1 ? "$p" : "$p^$e") for (p,e) in f.pe], " ⋅ ")
