module AlmanacForwardDiffExt

using Almanac
using ForwardDiff: ForwardDiff
using ForwardDiff: Dual, value, partials

# A query with a dual-number epoch does not run the evaluators on dual numbers. It takes
# the analytical time derivatives at the plain epoch (the code that the package image
# already holds) and assembles the dual result from them by the chain rule. The first
# call with a new `Dual` type then compiles only the small assembly functions below, and
# a query costs about as much as `spkder` of the same order.
#
# The partial of a result is therefore the next order of `spkder`/`pxder`. For an SPK
# type with a separate velocity series (type 3) the partial of the position is the
# stored velocity, as in `spkder`, and not the derivative of the position series.
#
# The generic evaluators stay reachable with `invoke` on the `et::Real` signature; the
# tests use this as an independent check of the analytical derivatives.

Almanac.primal(x::Dual) = Almanac.primal(value(x))

"`_ad_order(Dual)`: a dual number carries one derivative order more than its value type."
Almanac._ad_order(::Type{ForwardDiff.Dual{Tg, V, N}}) where {Tg, V, N} = 1 + Almanac._ad_order(V)

"`_with_primal(x::Dual, p)`: `x` with the plain value `p` and the same partials."
Almanac._with_primal(x::ForwardDiff.Dual{Tg}, p::Float64) where {Tg} =
    ForwardDiff.Dual{Tg}(Almanac._with_primal(ForwardDiff.value(x), p), ForwardDiff.partials(x))

"""
    lift(k, jet, x)

Return the derivative of order `k` of a function of time at the epoch `x`, where
`jet[m+1]` is the derivative of order `m` at `primal(x)`. Requires
`k + _ad_order(typeof(x)) < length(jet)`. For a dual `x` the value is order `k` at
`value(x)` and each partial is order `k+1` at `value(x)` times the partial of `x`.
This is the method of `_ad_lift` for a `Dual`: it builds the nested duals level by level,
not as the Taylor sum of the generic method.
"""
@inline lift(k::Int, jet::Tuple, ::Real) = jet[k + 1]
@inline function lift(k::Int, jet::Tuple, x::Dual{T}) where {T}
    v = lift(k, jet, value(x))
    d = lift(k + 1, jet, value(x))
    p = partials(x)
    return map((a, b) -> Dual{T}(a, b * p), v, d)
end

@inline Almanac._ad_lift(k::Int, jet::Tuple, x::Dual) = lift(k, jet, x)

# The methods of the derivative points (`@_ad_epoch_methods` of `src/ad.jl`) take the jet
# at the plain epoch and assemble the dual result from it with `lift`. A corrected
# position is assembled in the same way, from the jet of `spkder` with the correction at
# the plain epoch. The value is the position of the plain call, and the partials are the
# exact derivatives of the position function (see the docstring of `spkder`). For
# `spkezr` the value of the velocity is the velocity that CSPICE defines and its partials
# are the next orders of the same jet; for the `LT` forms that velocity is not the
# derivative of the position, so the partial of the position and the value of the
# velocity differ there, as they do between `spkder` and `spkezr`. The light time has the
# value of the plain call and the derivative of `|r| / c`, because the stellar aberration
# and the frame change are rotations. In the same way the value of order 1 of
# `frame_jet` (the rate block of `sxform`) is the rate of the plain call, and `pxder`
# keeps the exact derivative.
Almanac.@_ad_epoch_methods Dual

end # module
