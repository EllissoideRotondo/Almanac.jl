module AlmanacMooncakeExt

using Almanac
using Almanac: KernelSet, SpiceContext, Aberration, _ma_to_ta, _ad_pin, _ad_rule,
    _ad_scale, _AD_POINTS, _AD_CONSTANT, _AD_CONSTANT_ORBIT, _ad_typeof
using Mooncake: Mooncake
using Mooncake: DefaultCtx, CoDual, Dual, NoTangent, NoRData, RData, @is_primitive,
    @zero_derivative, zero_fcodual
using StaticArrays: StaticArray

# Forward and reverse rules of Mooncake for the derivative points of `src/ad.jl`. Mooncake
# does not differentiate the evaluators: a rule takes the value and the derivative from
# `_ad_rule`. Only the epoch has a tangent; the other arguments are lookups and codes.

# Mooncake tangent of a value (a static array is a struct with the field `data`).
_mt(x::StaticArray) = Mooncake.Tangent((data = Tuple(x),))
_mt(x::Tuple) = map(_mt, x)
_mt(x::Number) = x

# Reverse data of the result times its derivative, summed over the whole structure.
_mdot(::NoRData, ::Any) = 0.0
_mdot(a::RData, d::StaticArray) = _mdot(a.data.data, Tuple(d))
_mdot(a::Tuple, d::Tuple) = sum(map(_mdot, a, d))
_mdot(a::Number, d::Number) = a * d

# The rules of the derivative points (`_AD_POINTS`), for a `Float64` epoch.
for (f, types, iet) in _AD_POINTS
    sig = Tuple{typeof(f), [i == iet ? Float64 : types[i] for i in eachindex(types)]...}
    @eval begin
        @is_primitive DefaultCtx $sig
        function Mooncake.frule!!(::Dual{typeof($f)}, args::Vararg{Dual})
            y, dy = _ad_rule($f, map(Mooncake.primal, args)...)
            return Dual(y, _mt(_ad_scale(dy, Mooncake.tangent(args[$iet]))))
        end
        function Mooncake.rrule!!(::CoDual{typeof($f)}, args::Vararg{CoDual})
            y, dy = _ad_rule($f, map(Mooncake.primal, args)...)
            function pullback(ȳ)
                d = _mdot(ȳ, dy)
                return (NoRData(), ntuple(i -> i == $iet ? d : NoRData(), length(args))...)
            end
            return zero_fcodual(y), pullback
        end
    end
end

# The Kepler solver of the Orbit layer: both arguments, the mean anomaly and the
# eccentricity, have tangents (see `src/analysis/brouwer.jl`).
@is_primitive DefaultCtx Tuple{typeof(_ma_to_ta), Float64, Float64}
function Mooncake.frule!!(::Dual{typeof(_ma_to_ta)}, m::Dual, e::Dual)
    y, (dm, de) = _ad_rule(_ma_to_ta, Mooncake.primal(m), Mooncake.primal(e))
    return Dual(y, dm * Mooncake.tangent(m) + de * Mooncake.tangent(e))
end
function Mooncake.rrule!!(::CoDual{typeof(_ma_to_ta)}, m::CoDual, e::CoDual)
    y, (dm, de) = _ad_rule(_ma_to_ta, Mooncake.primal(m), Mooncake.primal(e))
    pullback(ȳ) = (NoRData(), ȳ * dm, ȳ * de)
    return zero_fcodual(y), pullback
end

# The value of one formula with the derivatives of another (see `src/ad.jl`).
@is_primitive DefaultCtx Tuple{typeof(_ad_pin), Float64, Float64}
Mooncake.frule!!(::Dual{typeof(_ad_pin)}, a::Dual, p::Dual) =
    Dual(_ad_pin(Mooncake.primal(a), Mooncake.primal(p)), Mooncake.tangent(a))
function Mooncake.rrule!!(::CoDual{typeof(_ad_pin)}, a::CoDual, p::CoDual)
    pullback(ȳ) = (NoRData(), ȳ, zero(ȳ))
    return zero_fcodual(_ad_pin(Mooncake.primal(a), Mooncake.primal(p))), pullback
end

# The snapshot, the context and the correction are data without derivatives, and the
# lookups and the parsing of strings have no derivative.
Mooncake.tangent_type(::Type{<:KernelSet}) = NoTangent
Mooncake.tangent_type(::Type{<:SpiceContext}) = NoTangent
Mooncake.tangent_type(::Type{Aberration}) = NoTangent
for f in (_AD_CONSTANT..., _AD_CONSTANT_ORBIT...)
    @eval @zero_derivative DefaultCtx Tuple{$(_ad_typeof(f)), Vararg}
end

end # module
