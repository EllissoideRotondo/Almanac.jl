module AlmanacChainRulesCoreExt

using Almanac
using Almanac: KernelSet, Aberration, Orbit, element, elements, _orbit_property, _ORBIT_FIELDS,
    spk_jet, frame_rotation, frame_jet, abcorr_jet,
    abcorr_position, abcorr_state, _ma_to_ta, _ad_rule, _AD_CONSTANT, _ad_typeof
using ChainRulesCore: ChainRulesCore
using ChainRulesCore: AbstractZero, AbstractThunk, NoTangent, Tangent, unthunk, backing,
    RuleConfig, HasReverseMode, rrule_via_ad
using LinearAlgebra: dot
using StaticArrays: StaticArray

# Forward and reverse rules for the derivative points of `src/ad.jl`. They serve Zygote,
# Diffractor and every other package that reads ChainRulesCore rules, and Mooncake reads
# them through its extension. The epoch is any `Real`, so a `ForwardDiff.Dual` epoch
# (forward mode over a reverse-mode rule, as for a Hessian) gets its derivative from the
# ForwardDiff extension.

# Cotangent of a value times its derivative, summed over the whole structure.
_cdot(::AbstractZero, ::Any) = false
_cdot(a::AbstractThunk, d) = _cdot(unthunk(a), d)
_cdot(a::Tangent, d) = _cdot(backing(a), d)
_cdot(a::NamedTuple, d::StaticArray) = _cdot(a.data, Tuple(d))
_cdot(a::Tuple, d::StaticArray) = _cdot(a, Tuple(d))
_cdot(a::Tuple, d::Tuple) = sum(map(_cdot, a, d))
_cdot(a::AbstractArray, d::AbstractArray) = dot(a, d)
_cdot(a::Number, d::Number) = a * d

# Tangent of the result for the epoch tangent `Δ`.
_tangent(y::Tuple, dy::Tuple, Δ) = Tangent{typeof(y)}(map((a, b) -> _tangent(a, b, Δ), y, dy)...)
_tangent(::Any, dy, Δ) = dy * Δ

# The position of the epoch in the arguments of each derivative point.
const POINTS = (
    (spk_jet, :(Tuple{KernelSet, Int32, Real, Int32, Int32, Val}), 3),
    (frame_rotation, :(Tuple{KernelSet, Int32, Int32, Real}), 4),
    (frame_jet, :(Tuple{KernelSet, Int32, Int32, Real, Val}), 4),
    (abcorr_jet, :(Tuple{KernelSet, Int32, Real, Int32, Aberration, Int32, Val}), 3),
    (abcorr_position, :(Tuple{KernelSet, Int32, Real, Int32, Aberration, Int32}), 3),
    (abcorr_state, :(Tuple{KernelSet, Int32, Real, Int32, Aberration, Int32}), 3),
)

for (f, sig, iet) in POINTS
    types = eval(sig).parameters
    n = length(types)
    args = [Symbol(:a, i) for i in 1:n]
    typed = [:($(args[i])::$(types[i])) for i in 1:n]
    nots = [i == iet ? :d : :(NoTangent()) for i in 1:n]
    @eval begin
        function ChainRulesCore.frule(
                (_, $([Symbol(:Δ, i) for i in 1:n]...)),
                ::typeof($f), $(typed...)
            )
            y, dy = _ad_rule($f, $(args...))
            return y, _tangent(y, dy, $(Symbol(:Δ, iet)))
        end
        function ChainRulesCore.rrule(::typeof($f), $(typed...))
            y, dy = _ad_rule($f, $(args...))
            function pullback(ȳ)
                d = _cdot(unthunk(ȳ), dy)
                return (NoTangent(), $(nots...))
            end
            return y, pullback
        end
    end
end

# The Kepler solver of the Orbit layer: two arguments, the mean anomaly and the
# eccentricity, with the partials of the implicit function (see `src/analysis/brouwer.jl`).
function ChainRulesCore.frule((_, Δm, Δe), ::typeof(_ma_to_ta), m::Real, e::Real)
    y, (dm, de) = _ad_rule(_ma_to_ta, m, e)
    return y, dm * Δm + de * Δe
end
function ChainRulesCore.rrule(::typeof(_ma_to_ta), m::Real, e::Real)
    y, (dm, de) = _ad_rule(_ma_to_ta, m, e)
    function pullback(ȳ)
        g = unthunk(ȳ)
        return (NoTangent(), g * dm, g * de)
    end
    return y, pullback
end

# `element(o, name)` and the properties of an `Orbit` take the name at run time, so Zygote
# would compile the pullback of every branch of `element` for the first property it
# meets (`o.state` as well), which takes about 35 s. These rules differentiate only the
# branch of the name, with the AD system of the caller.
function ChainRulesCore.rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(element),
        o::Orbit, name::Symbol)
    name in elements() || element(o, name)    # throws SpiceError(:badarg)
    y, pb = rrule_via_ad(config, getfield(Almanac, Symbol(:_el_, name)), o)
    return y, ȳ -> (NoTangent(), pb(ȳ)[2], NoTangent())
end
function ChainRulesCore.rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(getproperty),
        o::O, s::Symbol) where {O<:Orbit}
    if s in _ORBIT_FIELDS
        field_pb(ȳ) = (NoTangent(), Tangent{O}(; NamedTuple{(s,)}((unthunk(ȳ),))...), NoTangent())
        return getfield(o, s), field_pb
    end
    y, pb = rrule_via_ad(config, _orbit_property, o, s)
    return y, ȳ -> (NoTangent(), pb(ȳ)[2], NoTangent())
end

# Lookups, the snapshot of a context and the parsing of strings have no derivative.
for f in _AD_CONSTANT
    @eval begin
        ChainRulesCore.frule(ṫ, ::$(_ad_typeof(f)), args...; kw...) = ($f(args...; kw...), NoTangent())
        ChainRulesCore.rrule(::$(_ad_typeof(f)), args...; kw...) =
            ($f(args...; kw...), _ -> (NoTangent(), map(_ -> NoTangent(), args)...))
    end
end

end # module
