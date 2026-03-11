# Polyalgorithm dispatch for large integer factorization
# Ref: Cohen (1993) "A Course in Computational Algebraic Number Theory", §1.7

"""
    _find_factor(n::T) -> T

Find a non-trivial factor of composite `n` using a polyalgorithm:
1. Perfect power check (via IntegerMathUtils.ispower/iroot)
2. ECM (Elliptic Curve Method)
3. MPQS (Multiple Polynomial Quadratic Sieve)
"""
function _find_factor(n::T)::T where {T<:Integer}
    # 1. Perfect power check using IntegerMathUtils (GMP-backed)
    if ispower(n)
        d = find_exponent(n)
        r = iroot(n, d)
        return T(r)
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
