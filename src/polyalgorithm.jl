# Polyalgorithm dispatch for large integer factorization
# Ref: Cohen (1993) "A Course in Computational Algebraic Number Theory", §1.7

"""
    _check_perfect_power(n::T) -> Union{Tuple{T, Int}, Nothing}

Check if `n` is a perfect power `k^d` for `d ∈ {11, 7, 5, 3, 2}`.
Returns `(k, d)` if found, `nothing` otherwise.
"""
function _check_perfect_power(n::T)::Union{Tuple{T, Int}, Nothing} where {T<:Integer}
    for d in (11, 7, 5, 3, 2)
        if d == 2
            r = isqrt(n)
            if r * r == n
                return (r, 2)
            end
        else
            # Integer k-th root via floating point approximation + correction
            r = _integer_kth_root(n, d)
            if r > 0 && r^d == n
                return (T(r), d)
            end
        end
    end
    return nothing
end

"""
Compute the integer k-th root of n, i.e. floor(n^(1/k)).
Uses Newton's method for BigInt accuracy.
"""
function _integer_kth_root(n::T, k::Int)::T where {T<:Integer}
    n <= 0 && return T(0)
    if n < typemax(Int)
        # Use floating point as initial guess
        r = T(round(BigInt, BigFloat(n)^(BigFloat(1)/k)))
    else
        # Newton's method for large BigInt
        r = BigInt(round(BigInt, BigFloat(n)^(BigFloat(1)/k)))
        # Refine with Newton's method
        for _ in 1:100
            r_new = ((k - 1) * r + n ÷ r^(k - 1)) ÷ k
            if r_new >= r
                break
            end
            r = r_new
        end
        r = T(r)
    end
    # Check r-1, r, r+1 due to floating-point imprecision
    for candidate in (r - T(1), r, r + T(1))
        candidate > 0 || continue
        if candidate^k == n
            return candidate
        end
    end
    return T(0)
end

"""
    _find_factor(n::T) -> T

Find a non-trivial factor of composite `n` using a polyalgorithm:
1. Perfect power check
2. ECM (Elliptic Curve Method)
3. MPQS (Multiple Polynomial Quadratic Sieve)
"""
function _find_factor(n::T)::T where {T<:Integer}
    # 1. Perfect power check
    pp = _check_perfect_power(n)
    if pp !== nothing
        return pp[1]
    end

    # Convert to BigInt for ECM/MPQS
    nb = BigInt(n)

    # 2. Progressive ECM with increasing B1 bounds
    # Keep ECM budget moderate; MPQS is faster for balanced semiprimes.
    # ECM excels when one factor is much smaller than the other.
    d = ndigits(nb)
    ecm_schedule = if d >= 55
        # For large numbers, brief ECM then fall through to MPQS
        [(B1=2000, curves=10), (B1=11000, curves=20)]
    elseif d >= 40
        [(B1=2000, curves=25), (B1=11000, curves=90), (B1=50000, curves=200)]
    else
        [(B1=2000, curves=25), (B1=11000, curves=90)]
    end

    for (B1, curves) in ecm_schedule
        result = ecm_factor(nb, B1, curves)
        if result !== nothing
            return T(result)
        end
    end

    # 3. MPQS fallback
    result = mpqs_factor(nb)
    return T(result)
end
