"""
For a prime number p and a multiple q, the wheel index encodes the prime number index of 
p and q modulo 30, which can be encoded to a single number from 1 ... 64. This allows us
to jump immediately into the correct loop at the correct offset.
"""
create_jump(wheel_index, i) = :($wheel_index === $i && @goto $(Symbol(:x, i)))

create_label(wheel_index) = :(@label $(Symbol(:x, wheel_index)))

wheel_mask(prime_mod_30)::UInt8 = ~(0x01 << (to_idx(prime_mod_30) - 1))

"""
For any prime number `p` we compute its prime number index modulo 30 (here `wheel`) and we
generate the loop that crosses of the next 8 multiples that, modulo 30, are 
p * {1, 7, 11, 13, 17, 19, 23, 29}.
"""
function unrolled_loop(p_idx)
    p = ps[p_idx]

    # First push the stopping criterion
    unrolled_loop_body = Any[:(byte_idx > unrolled_max && break)]

    # Cross off the 8 next multiples
    for q in ps
        div, rem = divrem(p * q, 30)
        bit = wheel_mask(rem)
        push!(unrolled_loop_body, :(xs[byte_idx + increment * $(q - 1) + $div] &= $bit))
    end

    # Increment the byte index to where the next / 9th multiple is located
    push!(unrolled_loop_body, :(byte_idx += increment * 30 + $p))

    quote
        while true
            $(unrolled_loop_body...)
        end
    end
end

"""
The fan-in / fan-out phase that crosses off one multiple and then checks bounds; this is
before and after the unrolled loop starts and finishes respectively.
"""
function single_loop_item_not_unrolled(p_idx, q_idx, save_on_exit = true)
    # Our prime number modulo 30
    p = ps[p_idx]

    ps_next = (1, 7, 11, 13, 17, 19, 23, 29, 31)
    
    # Label name
    jump_idx = 8 * (p_idx - 1) + q_idx

    # Current and next multiplier modulo 30
    q_curr, q_next = ps_next[q_idx], ps_next[q_idx + 1]

    # Get the bit mask for crossing off p * q_curr
    div_curr, rem_curr = divrem(p * q_curr, 30)
    bit = wheel_mask(rem_curr)

    # Compute the increments for the byte index for the next multiple
    incr_bytes = p * q_next ÷ 30 - div_curr
    incr_multiple = q_next - q_curr
    
    quote
        # Todo: this generates an extra jump, maybe conditional moves are possible?
        if byte_idx > n_bytes

            # For a segmented sieve we store where we exit the loop, since that is the
            # entrypoint in the next loop; to avoid modulo computation to find the offset
            $(save_on_exit ? :(last_idx = $jump_idx) : nothing)
            @goto out
        end

        # Cross off the multiple
        xs[byte_idx] &= $bit

        # Increment the byte index to where the next multiple is located
        byte_idx += increment * $incr_multiple + $incr_bytes 
    end
end

"""
Full loop generates a potentially unrolled loop for a particular wheel
that may or may not save the exit point.
"""
function full_loop_for_wheel(wheel, unroll = true, save_on_exit = true)
    loop_statements = []

    for i = 1 : 8
        push!(loop_statements, create_label(8 * (wheel - 1) + i))
        unroll && i == 1 && push!(loop_statements, unrolled_loop(wheel))
        push!(loop_statements, single_loop_item_not_unrolled(wheel, i, save_on_exit))
    end

    quote
        while true
            $(loop_statements...)
        end
    end
end

"""
Generates a sieving loop that crosses off multiples of a given prime number.

    @sieve_loop :unroll :save_on_exit
    @sieve_loop
"""
macro sieve_loop(options...)
    unroll, save_on_exit = :(:unroll) in options, :(:save_on_exit) in options
    
    # When crossing off p * q where `p` is the siever prime and `q` the current multiplier
    # we have that p and q are {1, 7, 11, 13, 17, 19, 23, 29} mod 30.
    # For each of these 8 possibilities for `p` we create a loop, and per loop we
    # create 8 entrypoints to jump into. The first entrypoint is the unrolled loop for
    # whenever we can remove 8 multiples at the same time when all 8 fit in the interval
    # between byte_start:byte_next_start-1. Otherwise we can only remove one multiple at
    # a time. With 8 loops and 8 entrypoints per loop we have 64 different labels, numbered
    # x1 ... x64.

    # As an example, take p = 7 as a prime number and q = 23 as the first multiplier, and
    # assume our number line starts at 1 (so byte 1 represents 1:30, byte 2 represent 31:60). 
    # We have to cross off 7 * 23 = 161 first, which has byte index 6. Our prime number `p`
    # is in the 2nd spoke of the wheel and q is in the 7th spoke. This means we have to jump
    # to the 7th label in the 2nd loop; that is label 8 * (2 - 1) + 7 = 15. There we cross 
    # off the multiple (since 161 % 30 = 11 is the 3rd spoke, we "and" the byte with 0b11011111)
    # Then we move to 7 * 29 (increment the byte index accordingly), cross it off as well.
    # And now we enter the unrolled loop where 7 * {31, 37, ..., 59} are crossed off, then 
    # 7 * {61, 67, ..., 89} etc. Lastly we reach the end of the sieving interval, we cross
    # off the remaining multiples one by one, until the byte index is passed the end.
    # When that is the case, we save at which multiple / label we exited, so we can jump
    # there without computation when the next interval of the number line is sieved.

    esc(quote
        $(unroll ? :(unrolled_max = n_bytes - increment * 28 - 28) : nothing)

        # Create jumps inside loops
        $([create_jump(:wheel_idx, i) for i = 1 : 64]...)

        # # Create loops
        $([full_loop_for_wheel(wheel, unroll, save_on_exit) for wheel in 1 : 8]...)

        # Point of exit
        @label out
    end)
end