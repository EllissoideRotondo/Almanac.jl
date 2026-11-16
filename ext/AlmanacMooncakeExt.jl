module AlmanacMooncakeExt

using Almanac
using Almanac: KernelSet, SpiceContext, Aberration, spk_jet, frame_rotation, frame_jet,
    abcorr_jet, abcorr_position, abcorr_state, _ad_rule, _AD_CONSTANT, _ad_typeof
using Mooncake
using Mooncake: DefaultCtx, CoDual, Dual, NoTangent, NoRData, RData, @is_primitive,
    @zero_derivative, zero_fcodual
using StaticArrays: StaticArray

# Forward and reverse rules of Mooncake for the derivative points of `src/ad.jl`. Mooncake
# does not differentiate the evaluators: a rule takes the value and the derivative from
# `_ad_rule`. Only the epoch has a tangent; the other arguments are lookups and codes.

# Mooncake tangent of a value (a static array is a struct with the field `data`).
_mt(x::StaticArray) = Mooncake.Tangent((data=Tuple(x),))
_mt(x::Tuple) = map(_mt, x)
_mt(x::Number) = x

# The derivative `dy` times the epoch tangent `s`.
_scale(dy::Tuple, s) = map(d -> _scale(d, s), dy)
_scale(dy, s) = dy * s

# Reverse data of the result times its derivative, summed over the whole structure.
_mdot(::NoRData, ::Any) = 0.0
_mdot(a::RData, d::StaticArray) = _mdot(a.data.data, Tuple(d))
_mdot(a::Tuple, d::Tuple) = sum(map(_mdot, a, d))
_mdot(a::Number, d::Number) = a * d

const POINTS = (
    (spk_jet, :(Tuple{typeof(spk_jet),KernelSet,Int32,Float64,Int32,Int32,Val}), 3),
    (frame_rotation, :(Tuple{typeof(frame_rotation),KernelSet,Int32,Int32,Float64}), 4),
    (frame_jet, :(Tuple{typeof(frame_jet),KernelSet,Int32,Int32,Float64,Val}), 4),
    (abcorr_jet, :(Tuple{typeof(abcorr_jet),KernelSet,Int32,Float64,Int32,Aberration,Int32,Val}), 3),
    (abcorr_position, :(Tuple{typeof(abcorr_position),KernelSet,Int32,Float64,Int32,Aberration,Int32}), 3),
    (abcorr_state, :(Tuple{typeof(abcorr_state),KernelSet,Int32,Float64,Int32,Aberration,Int32}), 3),
)

for (f, sig, iet) in POINTS
    @eval begin
        @is_primitive DefaultCtx $sig
        function Mooncake.frule!!(::Dual{typeof($f)}, args::Vararg{Dual})
            y, dy = _ad_rule($f, map(Mooncake.primal, args)...)
            return Dual(y, _mt(_scale(dy, Mooncake.tangent(args[$iet]))))
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

# The snapshot, the context and the correction are data without derivatives, and the
# lookups and the parsing of strings have no derivative.
Mooncake.tangent_type(::Type{<:KernelSet}) = NoTangent
Mooncake.tangent_type(::Type{<:SpiceContext}) = NoTangent
Mooncake.tangent_type(::Type{Aberration}) = NoTangent
for f in _AD_CONSTANT
    @eval @zero_derivative DefaultCtx Tuple{$(_ad_typeof(f)),Vararg}
end

end # module
