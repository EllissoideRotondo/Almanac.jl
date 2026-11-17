module AlmanacEnzymeCoreExt

using Almanac
using Almanac: KernelSet, SpiceContext, Aberration, spk_jet, frame_rotation, frame_jet,
    abcorr_jet, abcorr_position, abcorr_state, _ma_to_ta, _ad_rule, _AD_CONSTANT, _ad_typeof
using EnzymeCore
using EnzymeCore: Annotation, Const, Active, Duplicated, BatchDuplicated
using EnzymeCore.EnzymeRules
using EnzymeCore.EnzymeRules: FwdConfig, RevConfig, AugmentedReturn, needs_primal,
    needs_shadow, width
using LinearAlgebra: dot

# Forward and reverse rules of Enzyme for the derivative points of `src/ad.jl`. Enzyme
# does not differentiate the evaluators: a rule takes the value and the derivative from
# `_ad_rule`. Only the epoch can be active; the other arguments are lookups and codes.
# In forward mode the epoch is any real number, so that ForwardDiff over Enzyme (a `Dual`
# epoch) takes its derivative from the ForwardDiff extension.

# The derivative `dy` times the epoch tangent `s`, and the dot product of a cotangent of
# the result with `dy`, over the whole structure (tuples of vectors, matrices, numbers).
_scale(dy::Tuple, s) = map(d -> _scale(d, s), dy)
_scale(dy, s) = dy * s
_edot(a::Tuple, d::Tuple) = sum(map(_edot, a, d))
_edot(a::AbstractArray, d::AbstractArray) = dot(a, d)
_edot(a::Number, d::Number) = a * d

_val(a::Annotation) = a.val

# The epoch tangents of each lane (an empty tuple for a constant epoch).
_lanes(et::Const, ::Val{W}) where {W} = ntuple(_ -> zero(et.val), Val(W))
_lanes(et::Duplicated, ::Val{1}) = (et.dval,)
_lanes(et::BatchDuplicated, ::Val{W}) where {W} = et.dval
_lanes(et::Annotation, ::Val{W}) where {W} = ntuple(i -> et.dval[i], Val(W))

function _forward(config::FwdConfig, f, iet::Int, args::Tuple)
    y, dy = _ad_rule(f, map(_val, args)...)
    W = width(config)
    return _fwd_return(config, y, map(s -> _scale(dy, s), _lanes(args[iet], Val(W))))
end

# The result of a forward rule with the value `y` and the shadows of each lane.
function _fwd_return(config::FwdConfig, y, shadows)
    W = width(config)
    if needs_primal(config) && needs_shadow(config)
        return W == 1 ? Duplicated(y, shadows[1]) : BatchDuplicated(y, shadows)
    elseif needs_shadow(config)
        return W == 1 ? shadows[1] : shadows
    elseif needs_primal(config)
        return y
    end
    return nothing
end

# The cotangent of the epoch for each lane, or `nothing` for a constant epoch.
_epoch_cotangent(::Const, dret, dy) = nothing
_epoch_cotangent(::Active, dret::Active, dy) = _edot(dret.val, dy)
_epoch_cotangent(::Active, dret::NTuple{W,Active}, dy) where {W} =
    map(d -> _edot(d.val, dy), dret)
_epoch_cotangent(et::Active, ::Type{<:Const}, dy) = zero(et.val)
_epoch_cotangent(et::Active, ::NTuple{W,Type{<:Const}}, dy) where {W} =
    ntuple(_ -> zero(et.val), Val(W))

@noinline _no_nested_reverse(et) = throw(ArgumentError("Almanac: Enzyme reverse mode " *
    "takes a Float64 epoch, not a $(typeof(et.val)); for a second derivative use " *
    "ForwardDiff over Enzyme forward mode, or over Zygote, ReverseDiff or Tracker"))

const POINTS = (
    (spk_jet, 3, 6), (frame_rotation, 4, 4), (frame_jet, 4, 5),
    (abcorr_jet, 3, 7), (abcorr_position, 3, 6), (abcorr_state, 3, 6),
)

for (f, iet, n) in POINTS
    args = [Symbol(:a, i) for i in 1:n]
    typed = [i == iet ? :($(args[i])::Annotation{<:Real}) : :($(args[i])::Const) for i in 1:n]
    plain = [i == iet ? :($(args[i])::Annotation{Float64}) : :($(args[i])::Const) for i in 1:n]
    @eval begin
        EnzymeRules.forward(config::FwdConfig, ::Const{typeof($f)}, ::Type, $(typed...)) =
            _forward(config, $f, $iet, ($(args...),))

        # The cotangent of a `Dual` is not a `Dual`, so ForwardDiff over reverse mode
        # would need a rule of its own; it throws instead of giving a wrong number.
        EnzymeRules.augmented_primal(::RevConfig, ::Const{typeof($f)}, ::Type, $(typed...)) =
            _no_nested_reverse($(args[iet]))

        function EnzymeRules.augmented_primal(config::RevConfig, ::Const{typeof($f)},
                ::Type, $(plain...))
            y, dy = _ad_rule($f, $([:($a.val) for a in args]...))
            return AugmentedReturn(needs_primal(config) ? y : nothing, nothing, dy)
        end

        function EnzymeRules.reverse(config::RevConfig, ::Const{typeof($f)}, dret, dy,
                $(plain...))
            d = _epoch_cotangent($(args[iet]), dret, dy)
            return ($([i == iet ? :d : :nothing for i in 1:n]...),)
        end
    end
end

# The Kepler solver of the Orbit layer: both arguments, the mean anomaly and the
# eccentricity, can be active (see `src/analysis/brouwer.jl`).
function EnzymeRules.forward(config::FwdConfig, ::Const{typeof(_ma_to_ta)}, ::Type,
        m::Annotation{<:Real}, e::Annotation{<:Real})
    y, (dm, de) = _ad_rule(_ma_to_ta, m.val, e.val)
    W = width(config)
    lm, le = _lanes(m, Val(W)), _lanes(e, Val(W))
    return _fwd_return(config, y, ntuple(i -> dm * lm[i] + de * le[i], Val(W)))
end
EnzymeRules.augmented_primal(::RevConfig, ::Const{typeof(_ma_to_ta)}, ::Type,
    m::Annotation{<:Real}, e::Annotation{<:Real}) = _no_nested_reverse(m)
function EnzymeRules.augmented_primal(config::RevConfig, ::Const{typeof(_ma_to_ta)},
        ::Type, m::Annotation{Float64}, e::Annotation{Float64})
    y, dy = _ad_rule(_ma_to_ta, m.val, e.val)
    return AugmentedReturn(needs_primal(config) ? y : nothing, nothing, dy)
end
EnzymeRules.reverse(config::RevConfig, ::Const{typeof(_ma_to_ta)}, dret, dy,
    m::Annotation{Float64}, e::Annotation{Float64}) =
    (_epoch_cotangent(m, dret, dy[1]), _epoch_cotangent(e, dret, dy[2]))

# The snapshot, the context and the correction are data without derivatives, and the
# lookups and the parsing of strings have no derivative.
EnzymeRules.inactive_type(::Type{<:KernelSet}) = true
EnzymeRules.inactive_type(::Type{<:SpiceContext}) = true
EnzymeRules.inactive_type(::Type{Aberration}) = true
for f in _AD_CONSTANT
    @eval EnzymeRules.inactive(::$(_ad_typeof(f)), args...) = nothing
end

end # module
