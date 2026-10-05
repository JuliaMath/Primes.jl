# Multiple Polynomial Quadratic Sieve (MPQS) for integer factorization
# Ref: Silverman (1987) "The Multiple Polynomial Quadratic Sieve"
# Ref: Pomerance (1982) "Analysis and Comparison of Some Integer Factoring Algorithms"
# Ref: Knuth & Trabb Pardo (1976) "Analysis of a Simple Factorization Algorithm"
# Ref: Crandall & Pomerance (2005) "Prime Numbers: A Computational Perspective", Ch.6

"""
Context for MPQS factorization.
"""
struct MPQSContext
    n::BigInt           # kn (multiplied)
    n_orig::BigInt      # original n
    factor_base::Vector{Int}
    fb_size::Int
    sqrt_kn_mod::Vector{Int}
    sieve_interval::Int
    # ceil(2^64 / p): d < 2^32 is divisible by p iff d * magic ≤ magic - 1 (wrapping)
    fb_magic::Vector{UInt64}
end

"""
A smooth or partially-smooth relation found during sieving.
large_prime == 0: fully smooth relation
large_prime > 0, large_prime2 == 0: single large prime partial
large_prime > 0, large_prime2 > 0: double large prime partial (large_prime < large_prime2)
"""
struct SmoothRelation
    # ax + b, reduced mod n. BigInt regardless of the working type: combining partials
    # multiplies two of these mod n.
    ax_plus_b::BigInt
    # g(x) over the factor base as column indices, one entry per power: 1 is the sign,
    # j + 1 is factor_base[j]. Sparse because most relations are partials that never
    # pair, and a dense row per relation costs a factor-base-sized copy each.
    factors::Vector{Int32}
    a_indices::Vector{Int}  # FB indices of primes composing polynomial `a`
    large_prime::Int     # first unfactored prime (0 if fully smooth)
    large_prime2::Int    # second unfactored prime (0 if single LP or smooth)
    squared_lps::Vector{Int} # large primes that combining has squared; each adds one factor to y
end

# Partials waiting for a partner. `single` holds one relation per large prime, kept so
# that every later relation sharing that prime combines with it. `waiting` indexes
# double-large-prime relations under both primes until one of them lands in `single`.
struct PartialPool
    single::Dict{Int, SmoothRelation}
    waiting::Dict{Int, Vector{SmoothRelation}}
end
PartialPool() = PartialPool(Dict{Int, SmoothRelation}(), Dict{Int, Vector{SmoothRelation}}())

# A polynomial is the triple (a, b, c) in the working type `T`, where
# c = (b² - kn) / a, so that g(x) = Q(x)/a = a·x² + 2b·x + c can be evaluated
# without ever forming Q(x) = (ax + b)² - kn itself.

# Precomputed table: digit count → (fb_size, sieve_interval)
# Ref: Silverman (1987), Crandall & Pomerance (2005) Ch.6
# Standard MPQS parameters: factor base sizes from L-smoothness formula.
# Large prime variation significantly amplifies effective relation yield.
# Ref: L-smoothness formula calibrated against reference implementations
const _MPQS_PARAMS = [
    (digits=30, fb_size=200,   sieve_interval=55000),
    (digits=35, fb_size=400,   sieve_interval=55000),
    (digits=40, fb_size=700,   sieve_interval=55000),
    (digits=45, fb_size=1200,  sieve_interval=55000),
    (digits=50, fb_size=1900,  sieve_interval=66000),
    (digits=55, fb_size=2900,  sieve_interval=66000),
    (digits=58, fb_size=3700,  sieve_interval=95000),
    (digits=60, fb_size=5500,  sieve_interval=250000),
    (digits=62, fb_size=6200,  sieve_interval=250000),
    (digits=65, fb_size=7500,  sieve_interval=250000),
    (digits=68, fb_size=9000,  sieve_interval=250000),
    (digits=70, fb_size=10000, sieve_interval=250000),
    (digits=73, fb_size=13000, sieve_interval=250000),
    (digits=76, fb_size=16000, sieve_interval=250000),
]

"""
Select MPQS parameters (factor base size, sieve interval) based on digit count.
"""
function _mpqs_select_params(n)
    d = ndigits(n)
    # Find bracketing entries and interpolate
    if d <= _MPQS_PARAMS[1].digits
        return _MPQS_PARAMS[1].fb_size, _MPQS_PARAMS[1].sieve_interval
    end
    if d >= _MPQS_PARAMS[end].digits
        return _MPQS_PARAMS[end].fb_size, _MPQS_PARAMS[end].sieve_interval
    end
    for i in 1:(length(_MPQS_PARAMS) - 1)
        lo, hi = _MPQS_PARAMS[i], _MPQS_PARAMS[i + 1]
        if lo.digits <= d <= hi.digits
            t = (d - lo.digits) / (hi.digits - lo.digits)
            fb = round(Int, lo.fb_size + t * (hi.fb_size - lo.fb_size))
            si = round(Int, lo.sieve_interval + t * (hi.sieve_interval - lo.sieve_interval))
            return fb, si
        end
    end
    return _MPQS_PARAMS[end].fb_size, _MPQS_PARAMS[end].sieve_interval
end

"""
Select optimal Knuth multiplier for MPQS.
Scores k ∈ {1,3,5,...,47} using Silverman's formula.
Ref: Silverman (1987) §4
"""
function _select_knuth_multiplier(n)::Int
    best_k = 1
    best_score = -Inf
    small_primes = primes(200)

    for k in 1:2:47
        kn = BigInt(k) * n
        score = 0.0

        # Special handling for p=2
        kn_mod8 = Int(mod(kn, 8))
        if kn_mod8 == 1
            score += 2 * log(2.0)
        elseif kn_mod8 == 5
            score += log(2.0)
        end

        # Score odd primes
        for p in small_primes
            p == 2 && continue
            if mod(k, p) == 0
                score += log(Float64(p)) / p
            elseif powermod(kn, div(p - 1, 2), p) == 1  # Legendre(kn, p) == 1
                score += 2 * log(Float64(p)) / p
            end
        end

        # Penalize large k slightly
        score -= 0.5 * log(Float64(k))

        if score > best_score
            best_score = score
            best_k = k
        end
    end
    return best_k
end

"""
Build factor base: primes p where Legendre(kn, p) == 1, plus p=2.
Also computes sqrt(kn) mod p and log₂(p) for each prime.
"""
function _build_factor_base(kn::Integer, fb_size_target::Int)
    factor_base = Int[2]
    sqrt_kn_mod = Int[mod(kn, 2) == 0 ? 0 : 1]

    p = 3
    while length(factor_base) < fb_size_target
        if isprime(p) && powermod(kn, div(p - 1, 2), p) == 1
            push!(factor_base, p)
            # kn mod p fits in an Int, so the root is found in machine arithmetic
            push!(sqrt_kn_mod, _tonelli_shanks(Int(mod(kn, p)), p))
        end
        p += 2
    end

    log_primes = UInt8[floor(UInt8, log2(p)) for p in factor_base]

    return factor_base, sqrt_kn_mod, log_primes
end

"""
Tonelli-Shanks algorithm for computing square root mod p.
Returns r such that r² ≡ n (mod p).
"""
function _tonelli_shanks(n::Int, p::Int)::Int
    n = mod(n, p)
    n == 0 && return 0
    p == 2 && return n

    # Factor out powers of 2 from p-1
    q = p - 1
    s = 0
    while iseven(q)
        q >>= 1
        s += 1
    end

    if s == 1
        # p ≡ 3 (mod 4)
        return powermod(n, div(p + 1, 4), p)
    end

    # Find a non-residue z
    z = 2
    while powermod(z, div(p - 1, 2), p) != p - 1
        z += 1
    end

    M = s
    c = powermod(z, q, p)
    t = powermod(n, q, p)
    R = powermod(n, div(q + 1, 2), p)

    while true
        t == 1 && return R
        # Find least i such that t^(2^i) ≡ 1 (mod p)
        i = 0
        temp = t
        while temp != 1
            temp = Int(mod(widemul(temp, temp), p))
            i += 1
        end
        b = powermod(c, 1 << (M - i - 1), p)
        M = i
        c = Int(mod(widemul(b, b), p))
        t = Int(mod(widemul(t, c), p))
        R = Int(mod(widemul(R, b), p))
    end
end

"""
Divide the factor base out of `remainder`, from prime `j0` on, pushing each power's
column onto `factors`. Returns `(remainder, next_j, status)`: status 0 when the scan
finished or the remainder reached 1; `j > 0` when it stopped at prime j with
remainder < p_j², so what is left is 1 or a prime; -1 when a BigInt remainder first
fits Int128.
"""
@inline function _tf_scan!(remainder::R, j0::Int, ctx::MPQSContext, pos0::Int,
                           offset1::Vector{Int}, offset2::Vector{Int},
                           factors::Vector{Int32}) where {R<:Integer}
    fb = ctx.factor_base
    @inbounds for j in j0:ctx.fb_size
        p = fb[j]
        o1 = offset1[j]
        # o1 < 0 marks a prime with no sieve root, which is divided directly instead.
        if o1 >= 0
            # A hardware `rem` per prime per candidate dominates otherwise. Adding p keeps
            # the difference positive; it stays below 2M + 1 + p < 2^32.
            magic = ctx.fb_magic[j]
            hit = (pos0 - o1 + p) % UInt64 * magic <= magic - 1
            if !hit
                o2 = offset2[j]
                hit = (o2 != o1) && (pos0 - o2 + p) % UInt64 * magic <= magic - 1
            end
            hit || continue
        end

        pR = R(p)
        divided = false
        while true
            q, r = divrem(remainder, pR)
            iszero(r) || break
            push!(factors, Int32(j + 1))
            divided = true
            remainder = q
        end
        divided || continue
        isone(remainder) && return remainder, j + 1, 0
        remainder < pR * pR && return remainder, j + 1, j
        R === BigInt && ndigits(remainder, base=2) < 128 && return remainder, j + 1, -1
    end
    return remainder, ctx.fb_size + 1, 0
end

"""
Root-guided trial factoring of g(x), which is held in the working type `T`.
`factors` is a reused buffer, emptied and refilled per candidate.
Returns a SmoothRelation or nothing.
"""
@inline function _trial_factor_guided(ax_b::T, gx::T, n_orig::BigInt, ctx::MPQSContext,
                               large_prime_bound::Int, dlp_bound::Int, dlp_bound_sq::Int,
                               sieve_pos::Int,
                               offset1::Vector{Int}, offset2::Vector{Int},
                               a_indices::Vector{Int}, factors::Vector{Int32}
                               )::Union{SmoothRelation, Nothing} where {T<:Integer}
    fb = ctx.factor_base
    fb_size = ctx.fb_size

    empty!(factors)
    gx < 0 && push!(factors, Int32(1))

    # Roots are 0-based (sieve position offset+1), so compare against sieve_pos - 1.
    pos0 = sieve_pos - 1
    remainder, j, status = _tf_scan!(abs(gx), 1, ctx, pos0, offset1, offset2, factors)
    if status == -1
        # A BigInt g(x) usually fits Int128 after its first few small factors; finish there.
        remainder, j, status = _tf_scan!(Int128(remainder), j, ctx, pos0, offset1, offset2,
                                         factors)
    end
    early_exit_j = status > 0 ? status : fb_size + 1

    # Handle early exit: remainder has at most one prime factor
    if early_exit_j <= fb_size && !isone(remainder)
        rem_int = Int(remainder)
        idx = searchsortedfirst(fb, rem_int)
        if idx <= fb_size && fb[idx] == rem_int
            push!(factors, Int32(idx + 1))
            remainder = one(remainder)
        end
    end

    # A relation is kept rarely enough that promoting ax+b to a BigInt here costs nothing.
    if isone(remainder)
        return SmoothRelation(mod(BigInt(ax_b), n_orig), copy(factors), a_indices, 0, 0, Int[])
    elseif remainder <= large_prime_bound && remainder > 1
        lp = Int(remainder)
        if isprime(lp)
            return SmoothRelation(mod(BigInt(ax_b), n_orig), copy(factors), a_indices, lp, 0, Int[])
        end
    elseif remainder <= dlp_bound_sq && remainder > 1 && !isprime(Int(remainder))
        # Double large prime: split the composite remainder, which is ≤ dlp_bound_sq and
        # so fits an Int. `eachfactor` rather than `lenstrafactor`, which throws on the
        # prime powers and small factors that turn up here.
        r = Int(remainder)
        f = first(first(eachfactor(r)))
        p1 = Int(min(f, div(r, f)))
        p2 = Int(max(f, div(r, f)))
        if p1 > 1 && p2 > 1 && p1 <= dlp_bound && p2 <= dlp_bound && isprime(p1) && isprime(p2)
            axb = mod(BigInt(ax_b), n_orig)
            # p² is already a square, so the relation is full as it stands
            p1 == p2 && return SmoothRelation(axb, copy(factors), a_indices, 0, 0, [p1])
            return SmoothRelation(axb, copy(factors), a_indices, p1, p2, Int[])
        end
    end

    return nothing
end

"""
Combine relation `r` with `single`, whose only unfactored prime is `shared_lp`, also a
large prime of `r`. The shared prime squares; `r`'s other large prime, if any, remains.
"""
function _combine_partials(r::SmoothRelation, single::SmoothRelation,
                           shared_lp::Int, ctx::MPQSContext)::SmoothRelation
    rest = r.large_prime == shared_lp ? r.large_prime2 : r.large_prime
    return SmoothRelation(mod(r.ax_plus_b * single.ax_plus_b, ctx.n_orig),
                          vcat(r.factors, single.factors),
                          vcat(r.a_indices, single.a_indices),
                          rest, 0,
                          vcat(r.squared_lps, single.squared_lps, shared_lp))
end

"""
    _gf2_eliminate(relations; max_deps=64) -> Vector{Vector{Int}}

GF(2) dependencies among parity vectors: sets of relation indices whose vectors XOR
to zero, at most `max_deps` of them.
"""
function _gf2_eliminate(relations::Vector{BitVector}; max_deps::Int=64)::Vector{Vector{Int}}
    nr = length(relations)
    nr == 0 && return Vector{Int}[]
    m = length(relations[1])

    # Column j of `mat` is prime equation j: bit r is set when relation r has an odd
    # exponent of that prime.
    nw = cld(nr, 64)
    mat = zeros(UInt64, nw, m)
    for (r, rel) in enumerate(relations), j in findall(rel)
        mat[(r - 1) >> 6 + 1, j] |= UInt64(1) << ((r - 1) & 63)
    end

    # Gauss-Jordan over the equations, pivoting on relations. No elimination history
    # is kept: the null space is read off the reduced matrix afterwards.
    pivot_of = fill(-1, m)
    is_pivot = falses(nr)
    @inbounds for i in 1:m
        w = findfirst(!iszero, @view mat[:, i])
        w === nothing && continue
        piv = (w - 1) * 64 + trailing_zeros(mat[w, i])
        pivot_of[i] = piv
        is_pivot[piv + 1] = true
        pbit = UInt64(1) << (piv & 63)
        for j in 1:m
            if j != i && mat[w, j] & pbit != 0
                for k in w:nw   # words below w are zero in equation i
                    mat[k, j] ⊻= mat[k, i]
                end
            end
        end
    end

    # Each free relation f gives one dependency: f together with the pivot relation of
    # every equation that still contains f.
    deps = Vector{Int}[]
    for f in 0:nr-1
        is_pivot[f + 1] && continue
        fw, fbit = f >> 6 + 1, UInt64(1) << (f & 63)
        dep = [f + 1]
        for i in 1:m
            pivot_of[i] >= 0 && mat[fw, i] & fbit != 0 && push!(dep, pivot_of[i] + 1)
        end
        push!(deps, dep)
        length(deps) == max_deps && break
    end
    return deps
end

# Exponent parities of Q(x) = a·g(x) over the factor base columns (1 is the sign).
function _parity(r::SmoothRelation, ncols::Int)
    v = falses(ncols)
    for col in r.factors
        v[col] ⊻= true
    end
    for ai in r.a_indices
        v[ai + 1] ⊻= true
    end
    return v
end

"""
Extract a factor from a dependency set using the stored factorizations.
`factors` lists g(x)'s factor base columns, one per power (column 1 is the sign);
a_indices stores which FB primes compose the polynomial's `a` value.
Q(x) = a · g(x), so exp_Q(p) = exp_g(p) + exp_a(p).
"""
function _extract_factor(n_orig::BigInt, kn::BigInt, k::Int,
                         dependency::Vector{Int},
                         relations::Vector{SmoothRelation},
                         factor_base::Vector{Int})::Union{BigInt, Nothing}
    x = BigInt(1)
    fb_size = length(factor_base)
    total_exp = zeros(Int, fb_size)

    for idx in dependency
        r = relations[idx]
        x = mod(x * r.ax_plus_b, n_orig)

        for col in r.factors
            col > 1 && (total_exp[col - 1] += 1)
        end
        # Add a exponents: each prime in a_indices contributes exponent 1
        for ai in r.a_indices
            total_exp[ai] += 1
        end
    end

    # Compute y = ∏ p^(e/2) mod n (all exponents should be even by GF(2) dependency)
    y = BigInt(1)
    for j in 1:fb_size
        e = total_exp[j]
        if e > 0
            y = mod(y * powermod(BigInt(factor_base[j]), e ÷ 2, n_orig), n_orig)
        end
    end

    # Each squared large prime appears as lp², contributing one lp to y.
    for idx in dependency
        for lp in relations[idx].squared_lps
            y = mod(y * BigInt(lp), n_orig)
        end
    end

    # x ≡ ±y (mod n) is the useless case; either sign of the difference can be the
    # one that splits n, so try both.
    for d in (x - y, x + y)
        g = _split_off_multiplier(gcd(abs(d), n_orig), n_orig, k)
        g === nothing || return g
    end
    return nothing
end

# A gcd that is divisible by the Knuth multiplier is a factor of kn, not of n;
# divide k back out and keep the result only if it really divides n.
function _split_off_multiplier(g::BigInt, n_orig::BigInt, k::Int)
    1 < g < n_orig || return nothing
    g2 = gcd(g, BigInt(k))
    if 1 < g2 < g
        g = div(g, g2)
    end
    return (1 < g < n_orig && iszero(mod(n_orig, g))) ? g : nothing
end

"""
Pick an `a` for the next batch of polynomials, as a product of `s` factor base primes.
Returns (a, a_indices, B_components), or nothing if no unused `a` was found.
"""
function _generate_siqs_a(::Type{T}, ctx::MPQSContext, used_a_sets::Set{Vector{Int}}) where {T<:Integer}
    fb = ctx.factor_base
    fb_size = ctx.fb_size
    target_a = isqrt(2 * ctx.n) ÷ ctx.sieve_interval

    # Sieve yield falls off sharply either side of target_a, so `a` is built to hit it:
    # draw s-1 primes near target_a^(1/s), then let the factor base pick the last.
    # Each `a` costs a pass of invmods and s·fb reductions but serves 2^(s-1)
    # polynomials, so s is set to make its primes small (near 2000): that setup cost
    # dwarfs the lost sieve root of each prime dividing `a`. Small factor bases end
    # below 2000, so s is also kept large enough for its primes to fit well inside.
    log_target = log(Float64(target_a))
    s = max(2, round(Int, log_target / log(2000.0)),
            ceil(Int, log_target / log(fb[fb_size] / 2)))
    center = exp(log_target / s)
    lo = clamp(searchsortedfirst(fb, floor(Int, center / 1.5)), 2, fb_size - s + 1)
    hi = min(max(searchsortedlast(fb, ceil(Int, center * 1.5)), lo + 4s), fb_size)

    for _ in 1:20
        indices = unique(sort(rand(lo:hi, s - 1)))
        length(indices) == s - 1 || continue
        partial = prod(BigInt(fb[i]) for i in indices)
        last = searchsortedfirst(fb, target_a ÷ partial)
        # Clamping instead would put `a` far off target, so drop the draw.
        (2 <= last <= fb_size && last ∉ indices) || continue
        push!(indices, last)
        sort!(indices)
        indices in used_a_sets && continue

        push!(used_a_sets, indices)
        a = partial * fb[last]
        return (T(a), indices, _crt_components(T, a, indices, fb, ctx.sqrt_kn_mod))
    end
    return nothing
end

# B_j = √(kn) mod q_j · (a/q_j) · inv(a/q_j, q_j) mod a, one per prime q_j of `a`.
# Summing them under every choice of sign gives the 2^(s-1) roots of b² ≡ kn (mod a).
function _crt_components(::Type{T}, a::BigInt, indices::Vector{Int},
                         fb::Vector{Int}, sqrt_kn_mod::Vector{Int}) where {T<:Integer}
    return map(indices) do idx
        q = BigInt(fb[idx])
        a_div_q = a ÷ q
        T(mod(sqrt_kn_mod[idx] * a_div_q * invmod(mod(a_div_q, q), q), a))
    end
end

# x mod p in [0, p). `mod(::BigInt, ::Int)` promotes p to a fresh BigInt on every call,
# which dominated per-polynomial setup once the working type is BigInt.
_modp(x::Integer, p::Int) = Int(mod(x, p))
_modp(x::BigInt, p::Int) =
    Int(ccall((:__gmpz_fdiv_ui, Base.GMP.libgmp), Culong, (Ref{BigInt}, Culong), x, p))

"""
Compute initial sieve root offsets from b, once per a-value.
offset1[j], offset2[j] ∈ [0, p-1]: sieve positions are offset+1, offset+1+p, ...
Use -1 as sentinel for "no root" (p | a and 2b ≡ 0 mod p).
"""
function _compute_siqs_roots!(offset1::Vector{Int}, offset2::Vector{Int},
                               b::T, c::T,
                               inv_a::Vector{Int},
                               sqrt_kn_mod::Vector{Int},
                               factor_base::Vector{Int},
                               fb_size::Int, M::Int) where {T<:Integer}
    @inbounds for j in 1:fb_size
        p = factor_base[j]
        if inv_a[j] == 0
            _set_single_root!(offset1, offset2, j, p, b, c, M)
        else
            sqr = sqrt_kn_mod[j]
            b_mod_p = _modp(b, p)
            ai = inv_a[j]
            r1 = mod((sqr - b_mod_p) * ai, p)
            r2 = mod((-sqr - b_mod_p) * ai, p)
            offset1[j] = mod(r1 + M, p)
            offset2[j] = mod(r2 + M, p)
        end
    end
end

# For p | a, Q(x)/a is linear mod p — 2b·x + c — so it has one root, not two.
# -1 in both offsets marks "no root at all", when 2b vanishes mod p as well.
@inline function _set_single_root!(offset1::Vector{Int}, offset2::Vector{Int},
                                   j::Int, p::Int, b::Integer, c::Integer, M::Int)
    b2 = mod(2 * _modp(b, p), p)
    if b2 == 0
        offset1[j] = -1
        offset2[j] = -1
    else
        o = mod(mod(-_modp(c, p) * invmod(b2, p), p) + M, p)
        offset1[j] = o
        offset2[j] = o
    end
end

# inv(a) mod p for every factor base prime. 0 marks the primes that divide a.
function _precompute_inv_a!(inv_a::Vector{Int}, factor_base::Vector{Int},
                            fb_size::Int, a::Integer)
    @inbounds for j in 1:fb_size
        p = factor_base[j]
        a_mod_p = _modp(a, p)
        inv_a[j] = a_mod_p == 0 ? 0 : invmod(a_mod_p, p)
    end
end

# B_delta[v][j] = 2·B_v·inv(a) mod p: how far prime j's roots move when the sign of B_v
# flips. Precomputing these is what makes the Gray-code walk over b-values cheap.
function _precompute_b_deltas!(B_delta::Vector{Vector{Int}}, B_comps::Vector,
                               inv_a::Vector{Int}, factor_base::Vector{Int}, fb_size::Int)
    for (v, Bv) in enumerate(B_comps)
        delta = B_delta[v]
        @inbounds for j in 1:fb_size
            p = factor_base[j]
            delta[j] = mod(2 * _modp(Bv, p) * inv_a[j], p)
        end
    end
end

# Slide every two-root prime's offsets by ±delta, the Gray-code step between b-values.
function _shift_roots!(offset1::Vector{Int}, offset2::Vector{Int}, delta::Vector{Int},
                       inv_a::Vector{Int}, factor_base::Vector{Int}, fb_size::Int, forward::Bool)
    @inbounds for j in 1:fb_size
        inv_a[j] == 0 && continue
        p = factor_base[j]
        d = delta[j]
        if forward
            o1 = offset1[j] + d; o1 >= p && (o1 -= p); offset1[j] = o1
            o2 = offset2[j] + d; o2 >= p && (o2 -= p); offset2[j] = o2
        else
            o1 = offset1[j] - d; o1 < 0 && (o1 += p); offset1[j] = o1
            o2 = offset2[j] - d; o2 < 0 && (o2 += p); offset2[j] = o2
        end
    end
end

# Subtract logp from every sieve position in the arithmetic progression pos, pos+p, ...
@inline function _sieve_stride!(sieve::Vector{UInt8}, sieve_len::Int, pos::Int,
                                p::Int, logp::UInt8)
    @inbounds while pos <= sieve_len
        sieve[pos] -= logp
        pos += p
    end
end

# Primes below _BLOCK_THRESH are sieved one _SIEVE_BLOCK-byte block at a time, so their
# dense updates stay in cache; larger primes hit a block about once and sweep the whole
# interval. Tuned at 70 digits: 1.4× over an unblocked 500 KB sieve.
const _SIEVE_BLOCK = 65536
const _BLOCK_THRESH = 16384

"""
Fast SIQS sieve: constant initialization (fill!) + unclamped subtraction + small prime
skipping. `pos1`/`pos2` are scratch for the blocked primes' next positions.
"""
function _siqs_sieve!(sieve::Vector{UInt8}, sieve_len::Int,
                      offset1::Vector{Int}, offset2::Vector{Int},
                      factor_base::Vector{Int}, log_primes::Vector{UInt8},
                      fb_size::Int, sieve_start_idx::Int, log_init::UInt8,
                      pos1::Vector{Int}, pos2::Vector{Int})
    fill!(sieve, log_init)
    nsmall = max(sieve_start_idx - 1, searchsortedfirst(factor_base, _BLOCK_THRESH) - 1)

    # o1 < 0 marks p | a with 2b ≡ 0 (mod p): no root at all. o2 == o1 is exactly the
    # p | a single-root case, which must not be subtracted twice; those few primes are
    # sieved unblocked. pos1 == 0 marks either as done.
    @inbounds for j in sieve_start_idx:nsmall
        o1, o2 = offset1[j], offset2[j]
        if o1 < 0
            pos1[j] = 0
        elseif o2 == o1
            _sieve_stride!(sieve, sieve_len, o1 + 1, factor_base[j], log_primes[j])
            pos1[j] = 0
        else
            pos1[j], pos2[j] = minmax(o1, o2) .+ 1
        end
    end
    @inbounds for block_start in 1:_SIEVE_BLOCK:sieve_len
        block_end = min(block_start + _SIEVE_BLOCK - 1, sieve_len)
        for j in sieve_start_idx:nsmall
            r1 = pos1[j]
            r1 == 0 && continue
            r2 = pos2[j]
            p = factor_base[j]
            logp = log_primes[j]
            while r2 <= block_end
                sieve[r1] -= logp
                sieve[r2] -= logp
                r1 += p
                r2 += p
            end
            # r2 - r1 < p, so after one more step on r1 the pair swaps order.
            if r1 <= block_end
                sieve[r1] -= logp
                r1, r2 = r2, r1 + p
            end
            pos1[j] = r1
            pos2[j] = r2
        end
    end

    @inbounds for j in nsmall+1:fb_size
        p = factor_base[j]
        logp = log_primes[j]
        o1 = offset1[j]
        o1 < 0 && continue
        o2 = offset2[j]
        if o2 == o1
            _sieve_stride!(sieve, sieve_len, o1 + 1, p, logp)
            continue
        end
        # Both roots in one loop: two independent update chains overlap in the pipeline.
        r1, r2 = minmax(o1, o2) .+ 1
        while r2 <= sieve_len
            sieve[r1] -= logp
            sieve[r2] -= logp
            r1 += p
            r2 += p
        end
        r1 <= sieve_len && (sieve[r1] -= logp)
    end
end

"""
Collect smooth candidates from sieve. Candidates are positions where the sieve
value underflowed (>= 0x80), indicating sufficient factorization over the factor base.
"""
function _siqs_collect!(a::T, b::T, c::T, a_factors::Vector{Int}, ctx::MPQSContext,
                        sieve::Vector{UInt8}, offset1::Vector{Int}, offset2::Vector{Int},
                        relations::Vector{SmoothRelation},
                        partial_relations::PartialPool,
                        M::Int, large_prime_bound::Int, dlp_bound::Int, dlp_bound_sq::Int,
                        tf_factors::Vector{Int32}) where {T<:Integer}
    sieve_len = length(sieve)
    num_chunks = div(sieve_len, 8)
    body_len = num_chunks * 8
    chunks = reinterpret(UInt64, @view sieve[1:body_len])

    # Scan 8 bytes at a time; a chunk with no high bit set holds no candidate, which is
    # the common case. The final iteration covers the tail, with no chunk test.
    @inbounds for j in 1:(num_chunks + 1)
        if j <= num_chunks
            chunks[j] & 0x8080808080808080 == 0 && continue
            first_pos, last_pos = (j - 1) * 8 + 1, j * 8
        else
            first_pos, last_pos = body_len + 1, sieve_len
        end
        for i in first_pos:last_pos
            sieve[i] < 0x80 && continue
            _process_candidate!(i, a, b, c, M, ctx,
                                large_prime_bound, dlp_bound, dlp_bound_sq,
                                offset1, offset2, a_factors, tf_factors,
                                relations, partial_relations)
        end
    end
end

"""
Process a single sieve candidate at position `i`.
"""
@inline function _process_candidate!(i::Int, a::T, b::T, c::T, M::Int,
                             ctx::MPQSContext,
                             large_prime_bound::Int, dlp_bound::Int, dlp_bound_sq::Int,
                             offset1::Vector{Int}, offset2::Vector{Int},
                             a_factors::Vector{Int}, tf_factors::Vector{Int32},
                             relations::Vector{SmoothRelation},
                             partial_relations::PartialPool) where {T<:Integer}
    x = T(i - M - 1)
    # g(x) = Q(x)/a by Horner. Evaluating g directly keeps every intermediate at g's
    # width; forming Q(x) = (ax+b)² - kn would need twice that.
    gx = (a * x + 2 * b) * x + c
    iszero(gx) && return
    ax_b = a * x + b

    relation = _trial_factor_guided(ax_b, gx, ctx.n_orig, ctx,
                                    large_prime_bound, dlp_bound, dlp_bound_sq,
                                    i, offset1, offset2,
                                    a_factors, tf_factors)
    relation === nothing && return

    _store_relation!(relation, relations, partial_relations, ctx)
end

"""
File a relation: smooth ones go into the pool, partials are combined with a stored
single-large-prime relation sharing a large prime, or parked to wait for one.
"""
function _store_relation!(relation::SmoothRelation,
                          relations::Vector{SmoothRelation},
                          pool::PartialPool, ctx::MPQSContext)
    lp1, lp2 = relation.large_prime, relation.large_prime2
    if lp1 == 0
        push!(relations, relation)
    elseif lp2 == 0
        other = get(pool.single, lp1, nothing)
        if other !== nothing
            push!(relations, _combine_partials(relation, other, lp1, ctx))
            return
        end
        pool.single[lp1] = relation
        # Double-large-prime relations parked on lp1 now reduce to their other prime.
        # Drained one at a time: the recursion can consume entries of this list too.
        ws = get(pool.waiting, lp1, nothing)
        ws === nothing && return
        while !isempty(ws)
            w = pop!(ws)
            rest = w.large_prime == lp1 ? w.large_prime2 : w.large_prime
            haskey(pool.waiting, rest) && filter!(x -> x !== w, pool.waiting[rest])
            _store_relation!(_combine_partials(w, relation, lp1, ctx), relations, pool, ctx)
        end
        delete!(pool.waiting, lp1)
    else
        key = haskey(pool.single, lp1) ? lp1 : haskey(pool.single, lp2) ? lp2 : 0
        if key == 0
            push!(get!(Vector{SmoothRelation}, pool.waiting, lp1), relation)
            push!(get!(Vector{SmoothRelation}, pool.waiting, lp2), relation)
        else
            _store_relation!(_combine_partials(relation, pool.single[key], key, ctx),
                             relations, pool, ctx)
        end
    end
end

"""
    mpqs_factor(n::Integer) -> BigInt

Return a non-trivial factor of `n` using the Self-Initializing Quadratic Sieve
(SIQS variant of MPQS), which requires `n` composite and not a perfect power.

Polynomial arithmetic runs in the narrowest type that holds it. The largest
intermediate is Horner's `(a·x + 2b)·x` at roughly `M·√(2kn)`, so anything fitting
126 bits sieves in `Int128`; only the mod-`n` steps need `n`'s full width.
"""
function mpqs_factor(n::Integer)
    nb = BigInt(n)
    # Relation collection runs until it splits n, so anything MPQS cannot split must be
    # refused rather than looped on. A prime has no split, and x² ≡ y² (mod p^e) does
    # not reliably find one for a perfect power. `eachfactor` excludes both already.
    nb > 3 || throw(ArgumentError("mpqs_factor needs n > 3, got $n"))
    isprime(n) && throw(ArgumentError("mpqs_factor needs a composite n, got the prime $n"))
    ispower(n) && throw(ArgumentError("mpqs_factor cannot reliably split the perfect power $n; factor its root"))
    k = _select_knuth_multiplier(n)
    kn = BigInt(k) * n
    fb_size_target, sieve_interval = _mpqs_select_params(n)
    width = ndigits(isqrt(2 * kn), base=2) + ndigits(sieve_interval, base=2) + 4
    return width <= 126 ?
        _mpqs_factor(Int128, nb, k, kn, fb_size_target, sieve_interval) :
        _mpqs_factor(BigInt, nb, k, kn, fb_size_target, sieve_interval)
end

function _mpqs_factor(::Type{T}, n::BigInt, k::Int, kn::BigInt,
                      fb_size_target::Int, sieve_interval::Int)::BigInt where {T<:Integer}
    factor_base, sqrt_kn_mod, log_primes = _build_factor_base(kn, fb_size_target)
    actual_fb_size = length(factor_base)

    fb_magic = UInt64[div(typemax(UInt64), p) + 1 for p in factor_base]
    ctx = MPQSContext(kn, n, factor_base, actual_fb_size, sqrt_kn_mod, sieve_interval, fb_magic)

    relations = SmoothRelation[]
    partial_relations = PartialPool()
    used_a_sets = Set{Vector{Int}}()
    target_relations = actual_fb_size + 50

    M = sieve_interval
    sieve_len = 2 * M + 1

    # Preallocate buffers (reused across all polynomials)
    sieve = Vector{UInt8}(undef, sieve_len)
    offset1 = Vector{Int}(undef, actual_fb_size)
    offset2 = Vector{Int}(undef, actual_fb_size)
    sieve_pos1 = Vector{Int}(undef, actual_fb_size)
    sieve_pos2 = Vector{Int}(undef, actual_fb_size)

    # Trial factoring buffer, reused across all candidates
    tf_factors = Int32[]

    # Constant sieve init: candidates are detected by UInt8 underflow (>= 0x80).
    # log_init calibrated to match reference implementation threshold.
    # nbits ≈ log2(sqrt(2kn)), up_to = 1.5 tolerance factor
    nbits = ndigits(isqrt(2 * kn), base=2)
    logp_max = floor(Int, log2(Float64(factor_base[end])))
    up_to = 1.5
    log_init = UInt8(clamp(nbits - round(Int, up_to * logp_max), 1, 127))

    # Skip small primes in sieve (they hit too many positions relative to their
    # log contribution). They are still checked during trial factoring.
    skip_limit = nbits > 120 ? 512 : (nbits > 90 ? 256 : 128)
    sieve_start_idx = 1
    while sieve_start_idx <= actual_fb_size && factor_base[sieve_start_idx] < skip_limit
        sieve_start_idx += 1
    end

    # Large prime bounds
    p_max = factor_base[end]
    large_prime_bound = (2 + (p_max >> 16)) * p_max * floor(Int, log2(Float64(p_max)))
    # Double large prime: accept remainder = p1 * p2 where both < dlp_bound
    dlp_bound = p_max * 100
    dlp_bound_sq = dlp_bound * dlp_bound

    # Gray-code root deltas, one row per prime factor of `a`. s grows with n, so this
    # count is a starting size, not a bound.
    B_delta = [Vector{Int}(undef, actual_fb_size) for _ in 1:10]
    inv_a = Vector{Int}(undef, actual_fb_size)

    # Sieve until the relation matrix yields a factor. No cap on a-values: given the
    # preconditions above a dependency eventually splits n, so a cap could only turn a
    # slow factorization into a spurious failure.
    while true
        if length(relations) >= target_relations
            dependencies = _gf2_eliminate([_parity(r, actual_fb_size + 1) for r in relations])
            for dep in dependencies
                result = _extract_factor(n, kn, k, dep, relations, factor_base)
                result === nothing || return result
            end
            # Every dependency gave x ≡ ±y (mod n) and so a trivial gcd. Each further
            # relation adds another dependency, each independently splitting n with
            # probability ≥ 1/2, so collect more rather than giving up.
            target_relations = length(relations) + 50
        end

        result = _generate_siqs_a(T, ctx, used_a_sets)
        result === nothing && continue
        a, a_indices, B_comps = result
        s = length(a_indices)
        while length(B_delta) < s
            push!(B_delta, Vector{Int}(undef, actual_fb_size))
        end

        _precompute_inv_a!(inv_a, factor_base, actual_fb_size, a)
        _precompute_b_deltas!(B_delta, B_comps, inv_a, factor_base, actual_fb_size)

        # Initial b (all positive CRT signs)
        b = mod(sum(B_comps), a)
        if mod(widemul(b, b), a) != mod(kn, a)
            b = a - b
        end
        # c ≈ M·√(kn/2) fits T, but b² alone does not — hence widemul. Once per
        # polynomial against a whole sieve pass, so the promotion is free.
        c = T(div(widemul(b, b) - kn, a))

        # Compute initial roots (once per a)
        _compute_siqs_roots!(offset1, offset2, b, c, inv_a,
                              sqrt_kn_mod, factor_base, actual_fb_size, M)

        # First polynomial
        _siqs_sieve!(sieve, sieve_len, offset1, offset2,
                     factor_base, log_primes, actual_fb_size, sieve_start_idx, log_init,
                     sieve_pos1, sieve_pos2)
        _siqs_collect!(a, b, c, a_indices, ctx, sieve, offset1, offset2,
                       relations, partial_relations, M, large_prime_bound,
                       dlp_bound, dlp_bound_sq, tf_factors)

        # Enough relations: hand back to the top of the loop, which runs the elimination.
        length(relations) >= target_relations && continue

        # Remaining b-values via Gray code incremental root update
        num_b = 1 << (s - 1)
        for i in 2:num_b
            gray_prev = (i - 2) ⊻ ((i - 2) >> 1)
            gray_curr = (i - 1) ⊻ ((i - 1) >> 1)
            flip_bit = trailing_zeros(gray_prev ⊻ gray_curr) + 1

            forward = (gray_curr >> (flip_bit - 1)) & 1 == 1
            # No reduction mod a: _shift_roots! assumes b moved by exactly ∓2B, and a wrap
            # by ±a would turn g(x) into g(x ± 1), off by one from every shifted root.
            b = forward ? b - 2 * B_comps[flip_bit] : b + 2 * B_comps[flip_bit]
            _shift_roots!(offset1, offset2, B_delta[flip_bit], inv_a,
                          factor_base, actual_fb_size, forward)

            # The primes dividing a are not shifted, so re-derive their single roots.
            c = T(div(widemul(b, b) - kn, a))
            for idx in a_indices
                _set_single_root!(offset1, offset2, idx, factor_base[idx], b, c, M)
            end

            # Sieve and collect
            _siqs_sieve!(sieve, sieve_len, offset1, offset2,
                         factor_base, log_primes, actual_fb_size, sieve_start_idx, log_init,
                         sieve_pos1, sieve_pos2)
            _siqs_collect!(a, b, c, a_indices, ctx, sieve, offset1, offset2,
                           relations, partial_relations, M, large_prime_bound,
                           dlp_bound, dlp_bound_sq, tf_factors)

            length(relations) >= target_relations && break
        end
    end
end
