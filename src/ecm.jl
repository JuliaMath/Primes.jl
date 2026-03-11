# Elliptic Curve Method (ECM) for integer factorization
# Ref: Lenstra (1987) "Factoring integers with elliptic curves"
# Ref: Montgomery (1987) "Speeding the Pollard and Elliptic Curve Methods of Factorization"

"""
Point on a Montgomery curve in projective coordinates (X:Z).
The point at infinity is represented by Z == 0.
"""
struct MontgomeryCurvePoint
    X::BigInt
    Z::BigInt
end

"""
In-place modular reduction: sets r = n mod d (non-negative remainder).
"""
function _mpz_fdiv_r!(r::BigInt, n::BigInt, d::BigInt)
    ccall((:__gmpz_fdiv_r, :libgmp), Cvoid, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), r, n, d)
end

"""
Preallocated scratch space for ECM point arithmetic.
Avoids BigInt allocation in the hot Montgomery ladder loop.
t1-t6: scratch for add!/double!; R0/R1/tmp: scratch for scalar_mul!
"""
struct ECMBuffers
    t1::BigInt
    t2::BigInt
    t3::BigInt
    t4::BigInt
    t5::BigInt
    t6::BigInt
    R0_X::BigInt
    R0_Z::BigInt
    R1_X::BigInt
    R1_Z::BigInt
    tmp_X::BigInt
    tmp_Z::BigInt
end

ECMBuffers() = ECMBuffers(BigInt(), BigInt(), BigInt(), BigInt(), BigInt(), BigInt(),
                          BigInt(), BigInt(), BigInt(), BigInt(), BigInt(), BigInt())

"""
In-place: mulmod!(dst, a, b, n, tmp) sets dst = (a * b) mod n using tmp as scratch.
"""
@inline function _mulmod!(dst::BigInt, a::BigInt, b::BigInt, n::BigInt, tmp::BigInt)
    Base.GMP.MPZ.mul!(tmp, a, b)
    _mpz_fdiv_r!(dst, tmp, n)
end

"""
Differential addition on Montgomery curve: given P, Q and P-Q, compute P+Q.
Uses projective coordinates and in-place arithmetic to avoid allocations.
"""
function _ecm_add!(res_X::BigInt, res_Z::BigInt,
                   P_X::BigInt, P_Z::BigInt, Q_X::BigInt, Q_Z::BigInt,
                   diff_X::BigInt, diff_Z::BigInt, n::BigInt, buf::ECMBuffers)
    t1, t2, t3, t4, t5, t6 = buf.t1, buf.t2, buf.t3, buf.t4, buf.t5, buf.t6
    # u = (P.X - P.Z) * (Q.X + Q.Z) mod n
    Base.GMP.MPZ.sub!(t1, P_X, P_Z)         # t1 = P.X - P.Z
    Base.GMP.MPZ.add!(t2, Q_X, Q_Z)         # t2 = Q.X + Q.Z
    _mulmod!(t5, t1, t2, n, t3)             # t5 = u

    # v = (P.X + P.Z) * (Q.X - Q.Z) mod n
    Base.GMP.MPZ.add!(t1, P_X, P_Z)         # t1 = P.X + P.Z
    Base.GMP.MPZ.sub!(t2, Q_X, Q_Z)         # t2 = Q.X - Q.Z
    _mulmod!(t6, t1, t2, n, t3)             # t6 = v

    # add = u + v, sub = u - v
    Base.GMP.MPZ.add!(t1, t5, t6)           # t1 = add = u + v
    Base.GMP.MPZ.sub!(t2, t5, t6)           # t2 = sub = u - v

    # X = diff.Z * add^2 mod n
    _mulmod!(t3, t1, t1, n, t4)             # t3 = add^2 mod n
    _mulmod!(res_X, diff_Z, t3, n, t4)      # res_X = diff.Z * add^2 mod n

    # Z = diff.X * sub^2 mod n
    _mulmod!(t3, t2, t2, n, t4)             # t3 = sub^2 mod n
    _mulmod!(res_Z, diff_X, t3, n, t4)      # res_Z = diff.X * sub^2 mod n
end

"""
In-place point doubling on Montgomery curve with parameter a24 = (a+2)/4.
"""
function _ecm_double!(res_X::BigInt, res_Z::BigInt,
                      P_X::BigInt, P_Z::BigInt,
                      n::BigInt, a24::BigInt, buf::ECMBuffers)
    t1, t2, t3, t4, t5, t6 = buf.t1, buf.t2, buf.t3, buf.t4, buf.t5, buf.t6
    # u = (P.X + P.Z)^2 mod n
    Base.GMP.MPZ.add!(t1, P_X, P_Z)         # t1 = P.X + P.Z
    _mulmod!(t5, t1, t1, n, t3)             # t5 = u = (P.X+P.Z)^2 mod n

    # v = (P.X - P.Z)^2 mod n
    Base.GMP.MPZ.sub!(t1, P_X, P_Z)         # t1 = P.X - P.Z
    _mulmod!(t6, t1, t1, n, t3)             # t6 = v = (P.X-P.Z)^2 mod n

    # diff = u - v
    Base.GMP.MPZ.sub!(t1, t5, t6)           # t1 = diff = u - v

    # X = u * v mod n
    _mulmod!(res_X, t5, t6, n, t3)          # res_X = u * v mod n

    # Z = diff * (v + a24 * diff) mod n
    _mulmod!(t2, a24, t1, n, t3)            # t2 = a24 * diff mod n
    Base.GMP.MPZ.add!(t2, t6)               # t2 = v + a24 * diff
    _mulmod!(res_Z, t1, t2, n, t3)          # res_Z = diff * (v + a24*diff) mod n
end

"""
Montgomery ladder scalar multiplication: compute [k]P on Montgomery curve.
Uses preallocated buffers to avoid allocation in the inner loop.
Returns the point [k]P as (res_X, res_Z).
"""
function _ecm_scalar_mul!(res_X::BigInt, res_Z::BigInt,
                          k::BigInt, P_X::BigInt, P_Z::BigInt,
                          n::BigInt, a24::BigInt, buf::ECMBuffers)
    R0_X, R0_Z = buf.R0_X, buf.R0_Z
    R1_X, R1_Z = buf.R1_X, buf.R1_Z
    tmp_X, tmp_Z = buf.tmp_X, buf.tmp_Z

    # R0 = P, R1 = 2P
    Base.GMP.MPZ.set!(R0_X, P_X)
    Base.GMP.MPZ.set!(R0_Z, P_Z)
    _ecm_double!(R1_X, R1_Z, P_X, P_Z, n, a24, buf)

    bits = ndigits(k, base=2)
    for i in (bits - 2):-1:0
        if isodd(k >> i)
            _ecm_add!(tmp_X, tmp_Z, R0_X, R0_Z, R1_X, R1_Z, P_X, P_Z, n, buf)
            Base.GMP.MPZ.set!(R0_X, tmp_X)
            Base.GMP.MPZ.set!(R0_Z, tmp_Z)
            _ecm_double!(tmp_X, tmp_Z, R1_X, R1_Z, n, a24, buf)
            Base.GMP.MPZ.set!(R1_X, tmp_X)
            Base.GMP.MPZ.set!(R1_Z, tmp_Z)
        else
            _ecm_add!(tmp_X, tmp_Z, R0_X, R0_Z, R1_X, R1_Z, P_X, P_Z, n, buf)
            Base.GMP.MPZ.set!(R1_X, tmp_X)
            Base.GMP.MPZ.set!(R1_Z, tmp_Z)
            _ecm_double!(tmp_X, tmp_Z, R0_X, R0_Z, n, a24, buf)
            Base.GMP.MPZ.set!(R0_X, tmp_X)
            Base.GMP.MPZ.set!(R0_Z, tmp_Z)
        end
    end
    Base.GMP.MPZ.set!(res_X, R0_X)
    Base.GMP.MPZ.set!(res_Z, R0_Z)
end

"""
    ecm_factor(n::BigInt, B1::Int, num_curves::Int) -> Union{BigInt, Nothing}

Attempt to find a non-trivial factor of `n` using the Elliptic Curve Method.
Computes [m]P where m = lcm(1..B1) = prod(p^floor(log_p(B1)) for p prime ≤ B1).
Uses batched gcd (accumulate Z coordinates, check periodically) to reduce gcd calls.
Returns a factor or `nothing` if none found within the curve budget.
"""
function ecm_factor(n::BigInt, B1::Int, num_curves::Int)::Union{BigInt, Nothing}
    # Precompute prime powers for Stage 1
    prime_powers = BigInt[]
    for p in primes(B1)
        pk = BigInt(p)
        while pk * p <= B1
            pk *= p
        end
        push!(prime_powers, pk)
    end

    buf = ECMBuffers()
    Q_X = BigInt()
    Q_Z = BigInt()
    tmp_mul = BigInt()  # scratch for acc * Q.Z

    for _ in 1:num_curves
        # Generate random curve via σ parameter (Suyama's parametrization)
        σ = BigInt(rand(6:10^9))
        u = mod(σ * σ - 5, n)
        v = mod(4 * σ, n)
        x0 = mod(u * u * u, n)
        z0 = mod(v * v * v, n)

        vu_diff = mod(v - u, n)
        a24_num = mod(vu_diff^3 * mod(3 * u + v, n), n)
        a24_den = mod(16 * x0 * v, n)

        g = gcd(a24_den, n)
        if g > 1 && g < n
            return g
        end
        if g == n
            continue
        end

        a24_den_inv = invmod(a24_den, n)
        a24 = mod(a24_num * a24_den_inv, n)

        Base.GMP.MPZ.set!(Q_X, x0)
        Base.GMP.MPZ.set!(Q_Z, z0)

        # Stage 1: multiply Q by each prime power, with batched gcd
        degenerate = false
        acc = BigInt(1)
        batch_count = 0
        for pk in prime_powers
            _ecm_scalar_mul!(Q_X, Q_Z, pk, Q_X, Q_Z, n, a24, buf)
            Base.GMP.MPZ.mul!(tmp_mul, acc, Q_Z)
            _mpz_fdiv_r!(acc, tmp_mul, n)
            batch_count += 1

            if batch_count >= 100
                g = gcd(acc, n)
                if g > 1 && g < n
                    return g
                end
                if g == n
                    degenerate = true
                    break
                end
                Base.GMP.MPZ.set_si!(acc, 1)
                batch_count = 0
            end
        end

        degenerate && continue

        if batch_count > 0
            g = gcd(acc, n)
            if g > 1 && g < n
                return g
            end
        end
    end
    return nothing
end
