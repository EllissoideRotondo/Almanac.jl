# Numeric helpers of the SPK evaluators: elementary functions that `Taylor` does not
# have yet, and the two functions that carry derivatives through an iterative solution.
# The vector routines are the `_m_` kernels of `math/kernels.jl`.

# --- Taylor series of sinh, cosh and atan ---------------------------------------------

"""
    _sinhcosh(x::Taylor) -> (sinh(x), cosh(x))

Coupled recurrences `s_k = (1/k) Σ_{j=1}^{k} j x_j c_{k-j}` and
`c_k = (1/k) Σ_{j=1}^{k} j x_j s_{k-j}`, started from `sinh(x_0)` and `cosh(x_0)`.
"""
@generated function _sinhcosh(x::Taylor{N, T}) where {N, T <: Real}
    ex = _unpack(:a, :x, N)
    push!(ex, :(s0 = sinh(a0)), :(c0 = cosh(a0)))
    for k in 1:N
        ss = _sum([:($j * $(_sym(:a, j)) * $(_sym(:c, k - j))) for j in 1:k])
        cs = _sum([:($j * $(_sym(:a, j)) * $(_sym(:s, k - j))) for j in 1:k])
        push!(ex, :($(_sym(:s, k)) = $ss / $k), :($(_sym(:c, k)) = $cs / $k))
    end
    return Expr(:block, Expr(:meta, :inline), ex..., :(($(_pack(:s, N)), $(_pack(:c, N)))))
end

@inline Base.sinh(x::Taylor) = _sinhcosh(x)[1]
@inline Base.cosh(x::Taylor) = _sinhcosh(x)[2]

"""
    atan(y::Taylor, x::Taylor)

Series of the angle of the point `(x, y)`: `r_0 = atan(y_0, x_0)` and
`r' = (x y' - y x') / (x² + y²)`, by the quotient recurrence on the coefficients of
the derivative. Requires `x_0² + y_0² > 0` (not checked).
"""
@generated function Base.atan(y::Taylor{N, T}, x::Taylor{N, T}) where {N, T <: Real}
    ex = vcat(_unpack(:s, :y, N), _unpack(:c, :x, N))
    push!(ex, :(r0 = atan(s0, c0)))
    for m in 0:(N - 1)
        # q_m: coefficient m of x² + y²; d: coefficient m of x y' - y x';
        # e_m: coefficient m of the derivative of the angle.
        q = _sum(
            [
                :($(_sym(:c, j)) * $(_sym(:c, m - j)) + $(_sym(:s, j)) * $(_sym(:s, m - j)))
                    for j in 0:m
            ]
        )
        d = _sum(
            [
                :(
                    $(j + 1) * (
                        $(_sym(:c, m - j)) * $(_sym(:s, j + 1)) -
                            $(_sym(:s, m - j)) * $(_sym(:c, j + 1))
                    )
                ) for j in 0:m
            ]
        )
        push!(ex, :($(_sym(:q, m)) = $q))
        if m == 0
            push!(ex, :(e0 = $d / q0))
        else
            sub = _sum([:($(_sym(:q, j)) * $(_sym(:e, m - j))) for j in 1:m])
            push!(ex, :($(_sym(:e, m)) = ($d - $sub) / q0))
        end
        push!(ex, :($(_sym(:r, m + 1)) = $(_sym(:e, m)) / $(m + 1)))
    end
    return Expr(:block, Expr(:meta, :inline), ex..., _pack(:r, N))
end

# --- derivatives through an iterative solution ----------------------------------------
#
# CSPICE solves the two Kepler equations (`prop2b`, `kpsolv`) by bisection, which
# carries no derivative. The evaluators here make the solution with plain `Float64`
# numbers exactly as CSPICE does, and then make Newton steps in the epoch number type
# on the residual with its plain value removed. Such a step does not change the plain
# value of the unknown, and each step makes at least one more derivative order exact;
# `_ad_order` steps are sufficient.

"""
    _ad_order(::Type{T}) -> Int

Number of derivative orders that a value of the number type `T` carries: 0 for a
plain float, `N` plus the order of the coefficient type for `Taylor{N}`, one more
than the value type for a `ForwardDiff.Dual` (package extension). Another `Real`
type gives 8, which is sufficient for a type that carries not more than 8 orders.
"""
@inline _ad_order(::Type{<:AbstractFloat}) = 0
@inline _ad_order(::Type{Taylor{N, T}}) where {N, T} = N + _ad_order(T)
@inline _ad_order(::Type{<:Real}) = 8

"""
    _with_primal(x::T, p::Float64) -> T

`x` with its plain value (`primal(x)`) replaced by `p`; the derivative parts do not
change. Postcondition: `primal(result) == p`.
"""
@inline _with_primal(x::AbstractFloat, p::Float64) = oftype(x, p)
@inline _with_primal(x::Taylor, p::Float64) = _taylor(_with_primal(x.v, p), x.d)
@inline _with_primal(x::Real, p::Float64) = x + (p - primal(x))

"""
    _ad_pin(a, p::Float64)

`_with_primal(a, p)`: the value `p` of one formula with the derivatives of another, `a`.
For a `Float64` it is `p`, and the rules of the ChainRulesCore, Enzyme and Mooncake
extensions give it the derivative 1 with respect to `a` and 0 with respect to `p`, so
the reverse-mode packages take the derivatives of `a`. `p` is a constant or comes from a
function that the extensions treat as constant (`_AD_CONSTANT`), so that no derivative is
taken through it. Not inlined: Enzyme applies its rule to a call that is still there,
and would otherwise differentiate `a + (p - a)`, whose derivative is 0.
"""
@noinline _ad_pin(a::Real, p::Float64) = _with_primal(a, p)
# The select is always `p` (`a == p` is false for a NaN `p`), but LLVM cannot prove it.
# A body that returns `p` itself gets the `returned` attribute, and Enzyme on Julia 1.10
# then takes the result as constant and calls the rule for neither value nor derivative.
@noinline _ad_pin(a::Float64, p::Float64) = ifelse(isnan(p) & (a == p), a, p)

"""
`|a|` with the sign of `b`, as `d_sign` of f2c: a zero `b` of either sign gives `+|a|`,
and `|a|` is `a` unless `a < 0` (so `-0.0` stays). The branches are taken on `primal`.
"""
@inline _dsign(a::Real, b::Real) = (x = primal(a) < 0.0 ? -a : a; primal(b) >= 0.0 ? x : -x)
