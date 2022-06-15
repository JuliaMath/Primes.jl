import Base: iterate

export SegmentedSieve

function generate_siever_primes(small_sieve::SmallSieve, segment_lo)
    xs = small_sieve.xs
    sievers = Vector{Siever}(undef, vec_count_ones(xs))
    j = 0
    @inbounds for i = eachindex(xs)
        x = xs[i]
        while x != 0x00
            sievers[j += 1] = Siever(compute_prime(x, i), segment_lo)
            x &= x - 0x01
        end
    end
    return sievers
end

struct SegmentIterator{T<:AbstractUnitRange}
    range::T
    segment_length::Int
    first_byte::Int
    last_byte::Int
    sievers::Vector{Siever}
    presieve_buffer::Vector{UInt8}
    segment::Vector{UInt8}
end

struct Segment{Tr,Ts}
    range::Tr
    segment::Ts
end

function Base.show(io::IO, s::Segment)
    # compute left padding
    padding = floor(Int, log10(last(s.range))) + 1

    padding_str = " " ^ padding

    print(io, padding_str, "  ")
    for p in ps
        print(io, lpad(p, 2, "0"), "  ")
    end

    println()

    for (start, byte) in zip(s.range, s.segment)
        mask = 0b00000001
        print(lpad(start, padding, "0"), "  ")
        for i = 1 : 8
            print(io, (byte & mask) == mask ? " x  " : " .  ")
            mask <<= 1
        end
        println()
    end

    io
end

function SegmentIterator(range::T, segment_length::Integer) where {T<:AbstractUnitRange}
    from, to = first(range), last(range)
    first_byte, last_byte = cld(first(range), 30), cld(last(range), 30)
    sievers = generate_siever_primes(SmallSieve(isqrt(to)), 30 * (first_byte - 1) + 1)
    presieve_buffer = create_presieve_buffer()
    xs = zeros(UInt8, segment_length)

    return SegmentIterator{T}(range, segment_length, first_byte, last_byte, sievers, presieve_buffer, xs)
end

function iterate(iter::SegmentIterator, segment_index_start = iter.first_byte)
    @inbounds begin
        if segment_index_start ≥ iter.last_byte
            return nothing
        end

        from, to = first(iter.range), last(iter.range)

        segment_index_next = min(segment_index_start + iter.segment_length, iter.last_byte + 1)
        segment_curr_len = segment_index_next - segment_index_start

        # Presieve
        apply_presieve_buffer!(iter.segment, iter.presieve_buffer, segment_index_start, segment_index_next - 1)

        # Set the preceding so many bits before `from` to 0
        if segment_index_start == iter.first_byte
            if iter.first_byte === 1
                iter.segment[1] = 0b11111110 # just make 1 not a prime.
            end
            for i = 1 : 8
                30 * (segment_index_start - 1) + ps[i] >= from && break
                iter.segment[1] &= wheel_mask(ps[i])
            end
        end

        # Set the remaining so many bits after `to` to 0
        if segment_index_next == iter.last_byte + 1
            for i = 8 : -1 : 1
                to ≥ 30 * (segment_index_next - 2) + ps[i] && break
                iter.segment[segment_curr_len] &= wheel_mask(ps[i])
            end
        end

        # Sieve the interval, but skip the pre-sieved primes
        xs = iter.segment

        for p_idx in 5:length(iter.sievers)
            p            = iter.sievers[p_idx]
            last_idx     = 0
            n_bytes      = segment_index_next - segment_index_start
            byte_idx     = p.byte_index - segment_index_start + 1
            wheel_idx    = p.wheel_index
            increment    = p.prime_div_30
            @sieve_loop :unroll :save_on_exit
            iter.sievers[p_idx] = Siever(increment, segment_index_start + byte_idx - 1, last_idx)
        end

        segment_start = 30 * (segment_index_start - 1)
        segment_stop = 30 * (segment_index_next - 1) - 1

        segment_index_start += iter.segment_length

        return Segment(segment_start:30:segment_stop, view(xs, Base.OneTo(segment_curr_len))), segment_index_start
    end
end
