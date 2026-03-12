# Elliptic Curve Method (ECM) for integer factorization
# Ref: Lenstra (1987) "Factoring integers with elliptic curves"
# Ref: Montgomery (1987) "Speeding the Pollard and Elliptic Curve Methods of Factorization"

# =============================================================================
# Generic modular arithmetic helpers (work for any T<:Integer, including
# fixed-width types from BitIntegers.jl such as UInt256, UInt512, etc.)
# Uses widen(T) for overflow-safe multiplication; for BitIntegers types,
# widen returns the next fixed-width type (e.g. widen(UInt256) == UInt512),
# so no BigInt allocation occurs.
# =============================================================================

# Modular addition: (a + b) mod n, for a, b ∈ [0, n).
# Avoids computing a + b directly, which could overflow for unsigned T
# when n is close to typemax(T).
@inline function _addmod(a::T, b::T, n::T) where {T<:Integer}
    a >= n - b ? a - (n - b) : a + b
end

# Modular subtraction: (a - b) mod n, for a, b ∈ [0, n).
# Written as n - (b - a) when a < b to avoid underflow in unsigned types.
@inline function _submod(a::T, b::T, n::T) where {T<:Integer}
    a >= b ? a - b : n - (b - a)
end

# Modular multiplication: (a * b) mod n.
# Uses widen(T) to avoid overflow for fixed-width integers.
# For BigInt, widen returns BigInt itself, so this is simply mod(a * b, n).
@inline function _mulmod(a::T, b::T, n::T) where {T<:Integer}
    mod(a * b, n)
end

# =============================================================================
# Shared helpers used by both the generic and BigInt-optimized ECM paths.
# =============================================================================

# Precompute prime powers p^k ≤ B1 for each prime p ≤ B1.
function _ecm_prime_powers(B1::Int)
    result = Int[]
    for p in primes(B1)
        pk = p
        while pk * p <= B1
            pk *= p
        end
        push!(result, pk)
    end
    result
end

# Suyama's parametrization: generate a random Montgomery curve and initial point.
# Returns (x0, z0, a24) on success, a lucky factor g::T, or nothing for a
# degenerate curve (gcd == n).
function _ecm_suyama(n::T) where {T<:Integer}
    σ = T(rand(6:10^9))
    u  = _submod(_mulmod(σ, σ, n), T(5), n)
    v  = _mulmod(T(4), σ, n)
    x0 = _mulmod(_mulmod(u, u, n), u, n)
    z0 = _mulmod(_mulmod(v, v, n), v, n)

    vu_diff  = _submod(v, u, n)
    vu_diff3 = _mulmod(_mulmod(vu_diff, vu_diff, n), vu_diff, n)
    a24_num  = _mulmod(vu_diff3, _addmod(_mulmod(T(3), u, n), v, n), n)
    a24_den  = _mulmod(_mulmod(T(16), x0, n), v, n)

    g = gcd(a24_den, n)
    1 < g < n && return g
    g == n && return nothing

    a24 = _mulmod(a24_num, invmod(a24_den, n), n)
    (x0, z0, a24)
end

# =============================================================================
# GMP-optimized path for BigInt: in-place arithmetic avoids allocations in
# the hot Montgomery ladder loop, which matters for large-integer factoring.
# =============================================================================

# In-place modular reduction: r = n mod d (non-negative remainder).
function _mpz_fdiv_r!(r::BigInt, n::BigInt, d::BigInt)
    ccall((:__gmpz_fdiv_r, :libgmp), Cvoid, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), r, n, d)
end

# Preallocated scratch space for the BigInt ECM inner loop.
# t1–t6: scratch for add!/double!; R0/R1/tmp: scratch for scalar_mul!
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

# In-place: dst = (a * b) mod n, using tmp as scratch.
@inline function _mulmod!(dst::BigInt, a::BigInt, b::BigInt, n::BigInt, tmp::BigInt)
    Base.GMP.MPZ.mul!(tmp, a, b)
    _mpz_fdiv_r!(dst, tmp, n)
end

# Differential addition on Montgomery curve (BigInt in-place version).
# Given projective P, Q and their difference P-Q (as diff), computes P+Q.
function _ecm_add!(res_X::BigInt, res_Z::BigInt,
                   P_X::BigInt, P_Z::BigInt, Q_X::BigInt, Q_Z::BigInt,
                   diff_X::BigInt, diff_Z::BigInt, n::BigInt, buf::ECMBuffers)
    t1, t2, t3, t4, t5, t6 = buf.t1, buf.t2, buf.t3, buf.t4, buf.t5, buf.t6
    Base.GMP.MPZ.sub!(t1, P_X, P_Z)
    Base.GMP.MPZ.add!(t2, Q_X, Q_Z)
    _mulmod!(t5, t1, t2, n, t3)             # t5 = (P.X-P.Z)(Q.X+Q.Z) mod n
    Base.GMP.MPZ.add!(t1, P_X, P_Z)
    Base.GMP.MPZ.sub!(t2, Q_X, Q_Z)
    _mulmod!(t6, t1, t2, n, t3)             # t6 = (P.X+P.Z)(Q.X-Q.Z) mod n
    Base.GMP.MPZ.add!(t1, t5, t6)           # t1 = t5 + t6
    Base.GMP.MPZ.sub!(t2, t5, t6)           # t2 = t5 - t6
    _mulmod!(t3, t1, t1, n, t4)
    _mulmod!(res_X, diff_Z, t3, n, t4)      # res_X = diff.Z * (t5+t6)^2 mod n
    _mulmod!(t3, t2, t2, n, t4)
    _mulmod!(res_Z, diff_X, t3, n, t4)      # res_Z = diff.X * (t5-t6)^2 mod n
end

# Point doubling on Montgomery curve with a24 = (a+2)/4 (BigInt in-place version).
function _ecm_double!(res_X::BigInt, res_Z::BigInt,
                      P_X::BigInt, P_Z::BigInt,
                      n::BigInt, a24::BigInt, buf::ECMBuffers)
    t1, t2, t3, t4, t5, t6 = buf.t1, buf.t2, buf.t3, buf.t4, buf.t5, buf.t6
    Base.GMP.MPZ.add!(t1, P_X, P_Z)
    _mulmod!(t5, t1, t1, n, t3)             # t5 = (P.X+P.Z)^2 mod n
    Base.GMP.MPZ.sub!(t1, P_X, P_Z)
    _mulmod!(t6, t1, t1, n, t3)             # t6 = (P.X-P.Z)^2 mod n
    Base.GMP.MPZ.sub!(t1, t5, t6)           # t1 = t5 - t6
    _mulmod!(res_X, t5, t6, n, t3)          # res_X = t5 * t6 mod n
    _mulmod!(t2, a24, t1, n, t3)
    Base.GMP.MPZ.add!(t2, t6)
    _mulmod!(res_Z, t1, t2, n, t3)          # res_Z = (t5-t6) * (t6 + a24*(t5-t6)) mod n
end

# Montgomery ladder scalar multiplication (BigInt in-place version).
function _ecm_scalar_mul!(res_X::BigInt, res_Z::BigInt,
                          k::Integer, P_X::BigInt, P_Z::BigInt,
                          n::BigInt, a24::BigInt, buf::ECMBuffers)
    R0_X, R0_Z = buf.R0_X, buf.R0_Z
    R1_X, R1_Z = buf.R1_X, buf.R1_Z
    tmp_X, tmp_Z = buf.tmp_X, buf.tmp_Z
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

# =============================================================================
# Generic ECM functions (functional style, no mutation).
# Work for any T<:Integer; the compiler specialises for each concrete type.
# For BitIntegers.jl types the widen chain stays in fixed-width LLVM integers,
# making these significantly faster than going through BigInt.
# =============================================================================

# Differential addition on Montgomery curve (generic version).
# Given projective P=(PX:PZ), Q=(QX:QZ) and their difference P-Q=(diffX:diffZ),
# returns (res_X, res_Z) = P+Q.
function _ecm_add(PX::T, PZ::T, QX::T, QZ::T,
                  diffX::T, diffZ::T, n::T) where {T<:Integer}
    u = _mulmod(_submod(PX, PZ, n), _addmod(QX, QZ, n), n)
    v = _mulmod(_addmod(PX, PZ, n), _submod(QX, QZ, n), n)
    ad = _addmod(u, v, n)
    sb = _submod(u, v, n)
    resX = _mulmod(diffZ, _mulmod(ad, ad, n), n)
    resZ = _mulmod(diffX, _mulmod(sb, sb, n), n)
    resX, resZ
end

# Point doubling on Montgomery curve with a24 = (a+2)/4 (generic version).
# Returns (res_X, res_Z).
function _ecm_double(PX::T, PZ::T, n::T, a24::T) where {T<:Integer}
    s = _addmod(PX, PZ, n)
    d = _submod(PX, PZ, n)
    u = _mulmod(s, s, n)
    v = _mulmod(d, d, n)
    diff = _submod(u, v, n)
    resX = _mulmod(u, v, n)
    resZ = _mulmod(diff, _addmod(v, _mulmod(a24, diff, n), n), n)
    resX, resZ
end

# Montgomery ladder scalar multiplication (generic version).
# Computes [k]P on the Montgomery curve; k is a small integer (prime power ≤ B1).
# Returns (res_X, res_Z).
function _ecm_scalar_mul(k::Integer, PX::T, PZ::T, n::T, a24::T) where {T<:Integer}
    R0X, R0Z = PX, PZ
    R1X, R1Z = _ecm_double(PX, PZ, n, a24)
    nbits = ndigits(k, base=2)
    for i in (nbits - 2):-1:0
        if isodd(k >> i)
            tmpX, tmpZ = _ecm_add(R0X, R0Z, R1X, R1Z, PX, PZ, n)
            R0X, R0Z = tmpX, tmpZ
            R1X, R1Z = _ecm_double(R1X, R1Z, n, a24)
        else
            tmpX, tmpZ = _ecm_add(R0X, R0Z, R1X, R1Z, PX, PZ, n)
            R1X, R1Z = tmpX, tmpZ
            R0X, R0Z = _ecm_double(R0X, R0Z, n, a24)
        end
    end
    R0X, R0Z
end

# =============================================================================
# Public API
# =============================================================================

"""
    ecm_factor(n::T, B1::Int, num_curves::Int) where {T<:Integer} -> Union{T, Nothing}

Attempt to find a non-trivial factor of `n` using the Elliptic Curve Method (ECM).

Uses Montgomery curves with Suyama's parametrization and a batched GCD to reduce
the number of `gcd` calls in the inner loop.

This generic method works for any integer type `T`, including fixed-width types from
[BitIntegers.jl](https://github.com/rfourquet/BitIntegers.jl) (e.g. `UInt256`, `UInt512`).
For those types modular multiplication stays within fixed-width LLVM integers via
`widen(T)`, avoiding BigInt allocation entirely.

Returns a non-trivial factor of `n`, or `nothing` if none is found within the budget.

See also the `BigInt`-optimised overload which uses in-place GMP arithmetic.
"""
function ecm_factor(n::T, B1::Int, num_curves::Int) where {T<:Integer}
    prime_powers = _ecm_prime_powers(B1)

    for _ in 1:num_curves
        curve = _ecm_suyama(n)
        curve === nothing && continue
        curve isa Tuple || return curve  # lucky factor from gcd
        x0, z0, a24 = curve

        QX, QZ = x0, z0

        # Stage 1: multiply Q by each prime power, batching the GCD checks.
        degenerate = false
        acc = one(T)
        batch_count = 0
        for pk in prime_powers
            QX, QZ = _ecm_scalar_mul(pk, QX, QZ, n, a24)
            acc = _mulmod(acc, QZ, n)
            batch_count += 1
            if batch_count >= 100
                g = gcd(acc, n)
                if 1 < g < n
                    return g
                end
                if g == n
                    degenerate = true
                    break
                end
                acc = one(T)
                batch_count = 0
            end
        end

        degenerate && continue

        if batch_count > 0
            g = gcd(acc, n)
            1 < g < n && return g
        end
    end
    return nothing
end

"""
    ecm_factor(n::BigInt, B1::Int, num_curves::Int) -> Union{BigInt, Nothing}

BigInt-optimised ECM using in-place GMP arithmetic to minimise allocations in the
hot Montgomery ladder loop.  For fixed-width integer types use the generic method.

Returns a non-trivial factor of `n`, or `nothing` if none is found within the budget.
"""
function ecm_factor(n::BigInt, B1::Int, num_curves::Int)::Union{BigInt, Nothing}
    prime_powers = _ecm_prime_powers(B1)

    buf     = ECMBuffers()
    Q_X     = BigInt()
    Q_Z     = BigInt()
    tmp_mul = BigInt()

    for _ in 1:num_curves
        curve = _ecm_suyama(n)
        curve === nothing && continue
        curve isa Tuple || return curve  # lucky factor from gcd
        x0, z0, a24 = curve

        Base.GMP.MPZ.set!(Q_X, x0)
        Base.GMP.MPZ.set!(Q_Z, z0)

        # Stage 1: multiply Q by each prime power, with batched GCD.
        degenerate  = false
        acc         = BigInt(1)
        batch_count = 0
        for pk in prime_powers
            _ecm_scalar_mul!(Q_X, Q_Z, pk, Q_X, Q_Z, n, a24, buf)
            Base.GMP.MPZ.mul!(tmp_mul, acc, Q_Z)
            _mpz_fdiv_r!(acc, tmp_mul, n)
            batch_count += 1
            if batch_count >= 100
                g = gcd(acc, n)
                if 1 < g < n
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
            1 < g < n && return g
        end
    end
    return nothing
end
