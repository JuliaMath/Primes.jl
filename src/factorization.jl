# implementation of a sorted dict (not optimized for speed) for storing
# the factorization of an integer

immutable Factorization{T<:Integer} <: Associative{T, Int}
    pe::Vector{Pair{T, Int}} # Prime-Exponent

    # Factorization{T}() where {T} = new(Vector{Pair{T, Int}}())
    (::Type{Factorization{T}}){T<:Integer}() = new{T}(Vector{Pair{T, Int}}())
end

function (::Type{Factorization{T}}){T<:Integer}(d::Associative)
    f = Factorization{T}()
    append!(f.pe, sort!(collect(d)))
    f
end

Base.convert{T}(::Type{Factorization}, d::Associative{T}) = Factorization{T}(d)

Base.start(f::Factorization) = start(f.pe)
Base.next(f::Factorization, i) = next(f.pe, i)
Base.done(f::Factorization, i) = done(f.pe, i)


function Base.get(f::Factorization, p, default)
    found = searchsorted(f.pe, p, by=first)
    isempty(found) ?
        default :
        last(f.pe[first(found)])
end

Base.getindex(f::Factorization, p::Integer) = get(f, p, 0)

function Base.setindex!{T}(f::Factorization{T}, e::Int, p::Integer)
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
