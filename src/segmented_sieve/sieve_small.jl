@inline compute_prime(x::UInt8, i) = 30 * (i - 1) + ps[trailing_zeros(x) + 0x01] 

"""
An iterator for finding small (<= 1_000_000) primes.

    for p in SmallSieve(1_000_000)
        println(p)
    end

Uses n ÷ 30 bytes of memory. Skips 2, 3, and 5; starts at 7.
"""
struct SmallSieve
    xs::Vector{UInt8}

    function SmallSieve(n::Integer)
        # Unrolled loop without segments
        n_bytes = cld(n, 30)
        xs = fill(0xFF, n_bytes)

        # Ensure `1` is not a prime number
        @inbounds xs[1] &= wheel_mask(1)

        # And ensure numbers > n are not prime since we are not considering them
        @inbounds for i = 8 : -1 : 1
            n >= 30 * (n_bytes - 1) + ps[i] && break
            xs[n_bytes] &= wheel_mask(ps[i])
        end

        hi = isqrt(n)

        @inbounds for i = eachindex(xs)
            x = xs[i]
            while x != 0x00
                # The next prime number
                p = compute_prime(x, i)

                # Are we done yet?
                p > hi && @goto done
                
                # Otherwise cross off multiples of p starting at p²
                p²        = p * p
                byte_idx  = p² ÷ 30 + 1
                wheel     = to_idx(p % 30)
                wheel_idx = 8 * (wheel - 1) + wheel
                increment = i - 1
                @sieve_loop :unroll # Just unroll -- no segmented business here
                x &= x - 0x01
            end
        end

        @label done

        return new(xs)
    end
end

# Yes, this is an O(n) computation, but it should be much faster than computing
# the prime numbers from the bitmask anways
Base.length(s::SmallSieve) = vec_count_ones(s.xs)
Base.eltype(s::SmallSieve) = Int

@inline function Base.iterate(s::SmallSieve)
    @inbounds for i = eachindex(s.xs)
        x = s.xs[i]
        x !== 0x00 && return compute_prime(x, i), (x & (x - 0x01), i)
    end
    
    return nothing
end

@inline function Base.iterate(s::SmallSieve, state::Tuple{UInt8, Int})
    x, i = state
    @inbounds while true
        x !== 0x00 && return compute_prime(x, i), (x & (x - 0x01), i)
        i === length(s.xs) && return nothing
        x = s.xs[i += 1]
    end
end

function small_primes(n)
    # This runs 20% faster than collect(SmallSieve(n)) :(
    xs = SmallSieve(n).xs
    primes = Vector{Int}(undef, vec_count_ones(xs))
    j = 0
    @inbounds for i = eachindex(xs)
        x = xs[i]
        while x != 0x00
            primes[j += 1] = compute_prime(x, i)
            x &= x - 0x01
        end
    end

    primes
end