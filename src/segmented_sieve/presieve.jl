"""
The {2, 3, 5}-wheel is efficient because it compresses memory
perfectly (1 byte = 30 numbers) and removes 4/15th of all the
multiples already. We don't get that memory efficiency when 
extending the wheel to {2, 3, 5, 7}, since we need to store
48 bits per 210 numbers, which is could be done with one 64-bit
integer per 210 numbers, which in fact compresses worse 
(64 bits / 210 numbers) than the {2, 3, 5}-wheel 
(8 bits / 30 numbers).

What we can do however, is compute the repeating pattern that
the first `n` primes create, and copy that pattern over. That
is, we look at a the numbers modulo p₁ * p₂ * ⋯ * pₙ.

For instance, when presieving all multiples of {2, 3, ..., 19}
we allocate a buffer for the range 1 : 2 * 3 * ... * 19 = 
1:9_699_690. In a {2, 3, 5} wheel this means a buffer of 
9_699_690 ÷ 30 = 323_323 bytes.
"""
function create_presieve_buffer()
    n_bytes = 7 * 11 * 13 * 17
    xs = fill(0xFF, n_bytes)

    @inbounds for p in (7, 11, 13, 17)
        p²        = p * p
        byte_idx  = p² ÷ 30 + 1
        wheel     = to_idx(p)
        wheel_idx = 8 * (wheel - 1) + wheel
        increment = 0
        @sieve_loop :unroll
    end

    @inbounds xs[1] = 0b11100001 # remove 7, 11, 13 and 17
    return xs
end

"""
When applying the presieve buffer, we have to compute the offset in 
"""
function apply_presieve_buffer!(xs::Vector{UInt8}, buffer::Vector{UInt8}, byte_start, byte_stop)

    len = byte_stop - byte_start + 1

    # todo, clean this up a bit.
    from_idx = (byte_start - 1) % length(buffer) + 1
    to = min(len, length(buffer) - from_idx + 1)

    # First copy the remainder of buffer at the front
    copyto!(view(xs, Base.OneTo(to)), view(buffer, from_idx:from_idx + to - 1))
    from = to + 1

    # Then copy buffer multiple times
    while from + length(buffer) - 1 <= len
        copyto!(view(xs, from : from + length(buffer) - 1), buffer)
        from += length(buffer)
    end

    # And finally copy the remainder of buffer again
    copyto!(view(xs, from:len), view(buffer, Base.OneTo(length(from:len))))

    xs
end