function compute_offset(p::Integer, segment_lo::Integer)
    p² = p * p
    if p² >= segment_lo
        byte = cld(p², 30)

        # Wheel index stores both the index of the prime number mod 30
        # and the index of the active multiple. We start crossing off
        # p * p, so that would be to_idx(p % 30) twice. We combine those values
        # as a number between 1 ... 64.
        wheel = to_idx(p % 30)
        wheel_index = 8 * (wheel - 1) + wheel

        return byte, wheel_index
    else
        # p * q will be the first number to cross off
        q = cld(segment_lo, p)
        q_quot, q_rem = divrem(q, 30)

        remainders = (1, 7, 11, 13, 17, 19, 23, 29, 31)
        i = 1
        while remainders[i] < q_rem
            i += 1
        end

        # maybe wrap around
        q_rem = i == 9 ? 1 : remainders[i]
        q_quot = i == 9 ? q_quot + 1 : q_quot
        i = i == 9 ? 1 : i

        # Our actual first acceptable multiple
        q = 30q_quot + q_rem

        byte = cld(p * q, 30)
        wheel_index = 8 * (to_idx(p % 30) - 1) + i
        
        return byte, wheel_index
    end
end

struct Siever
    prime_div_30::Int

    # byte_index is the integer range 30(byte_index - 1) up to 30byte_index - 1
    byte_index::Int

    # Stores the next prime number to be crossed off. 
    # If `p` is the prime number and `q` the next multiple to be stored
    # Wheel index 8 * to_idx(p % 30) * to_idx(q % 30)
    wheel_index::Int

    function Siever(p::Int, segment_lo::Int)
        byte, wheel = compute_offset(p, segment_lo)
        return new(p ÷ 30, byte, wheel)
    end
    
    Siever(prime_div_30, byte_index, wheel_index) = new(prime_div_30, byte_index, wheel_index)
end

function Base.show(io::IO, p::Siever)
    print(io, 30p.prime_div_30 + ps[(p.wheel_index - 1) ÷ 8 + 1], " (", p.byte_index, ", ", p.wheel_index, ")")
end
