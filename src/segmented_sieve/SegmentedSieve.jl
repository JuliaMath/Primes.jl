module SegmentedSieve

const ps = (1, 7, 11, 13, 17, 19, 23, 29)

"""
Population count of a vector of UInt8s for counting prime numbers.
See https://github.com/JuliaLang/julia/issues/34059
"""
function vec_count_ones(xs::Union{Vector{UInt8}, Base.FastContiguousSubArray{UInt8}})
    n = length(xs)
    count = 0
    chunks = n ÷ sizeof(UInt)
    GC.@preserve xs begin
        ptr = Ptr{UInt}(pointer(xs))
        for i in 1:chunks
            count += count_ones(unsafe_load(ptr, i))
        end
    end

    @inbounds for i in 8chunks+1:n
        count += count_ones(xs[i])
    end

    count
end

function to_idx(x)
    x ==  1 && return 1
    x ==  7 && return 2
    x == 11 && return 3
    x == 13 && return 4
    x == 17 && return 5
    x == 19 && return 6
    x == 23 && return 7
    return 8
end

include("generate_sieving_loop.jl")
include("sieve_small.jl")
include("presieve.jl")
include("siever.jl")
include("sieve.jl")

end # module
