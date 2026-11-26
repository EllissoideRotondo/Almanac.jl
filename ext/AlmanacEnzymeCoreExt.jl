module AlmanacEnzymeCoreExt

using Almanac
using Almanac: KernelSet, SpiceContext, Aberration, _ad_pin, _ad_rule,
    _ad_scale, _AD_POINTS, _AD_CONSTANT, _ad_typeof
using EnzymeCore: EnzymeCore, EnzymeRules, Annotation, Const, Active, Duplicated,
    BatchDuplicated
using EnzymeCore.EnzymeRules: FwdConfig, RevConfig, AugmentedReturn, needs_primal,
    needs_shadow, width
using LinearAlgebra: dot

# Forward and reverse rules of Enzyme for the derivative points of `src/ad.jl`. Enzyme
# does not differentiate the evaluators: a rule takes the value and the derivative from
# `_ad_rule`. Only the epoch can be active; the other arguments are lookups and codes.
# In forward mode the epoch is any real number, so that ForwardDiff over Enzyme (a `Dual`
# epoch) takes its derivative from the ForwardDiff extension.

# The dot product of a cotangent of the result with the derivative `dy`, over the whole
# structure (tuples of vectors, matrices, numbers).
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
    return _fwd_return(config, y, map(s -> _ad_scale(dy, s), _lanes(args[iet], Val(W))))
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
_epoch_cotangent(::Active, dret::NTuple{W, Active}, dy) where {W} =
    map(d -> _edot(d.val, dy), dret)
_epoch_cotangent(et::Active, ::Type{<:Const}, dy) = zero(et.val)
_epoch_cotangent(et::Active, ::NTuple{W, Type{<:Const}}, dy) where {W} =
    ntuple(_ -> zero(et.val), Val(W))

@noinline _no_nested_reverse(et) = throw(
    ArgumentError(
        "Almanac: Enzyme reverse mode " *
            "takes a Float64 epoch, not a $(typeof(et.val)); for a second derivative use " *
            "ForwardDiff over Enzyme forward mode, or over Zygote, ReverseDiff or Tracker"
    )
)

# The rules of the derivative points (`_AD_POINTS`: the types of the arguments and the
# position of the epoch among them).
for (f, types, iet) in _AD_POINTS
    n = length(types)
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

        function EnzymeRules.augmented_primal(
                config::RevConfig, ::Const{typeof($f)},
                ::Type, $(plain...)
            )
            y, dy = _ad_rule($f, $([:($a.val) for a in args]...))
            return AugmentedReturn(needs_primal(config) ? y : nothing, nothing, dy)
        end

        function EnzymeRules.reverse(
                config::RevConfig, ::Const{typeof($f)}, dret, dy,
                $(plain...)
            )
            d = _epoch_cotangent($(args[iet]), dret, dy)
            return ($([i == iet ? :d : :nothing for i in 1:n]...),)
        end
    end
end

# The value of one formula with the derivatives of another (see `src/ad.jl`): the
# derivative is 1 with respect to `a` and 0 with respect to `p`.
function EnzymeRules.forward(
        config::FwdConfig, ::Const{typeof(_ad_pin)}, ::Type,
        a::Annotation{<:Real}, p::Annotation{Float64}
    )
    return _fwd_return(config, _ad_pin(a.val, p.val), _lanes(a, Val(width(config))))
end
EnzymeRules.augmented_primal(
    ::RevConfig, ::Const{typeof(_ad_pin)}, ::Type,
    a::Annotation{<:Real}, p::Annotation{Float64}
) = _no_nested_reverse(a)
EnzymeRules.augmented_primal(
    config::RevConfig, ::Const{typeof(_ad_pin)}, ::Type,
    a::Annotation{Float64}, p::Annotation{Float64}
) = AugmentedReturn(needs_primal(config) ? _ad_pin(a.val, p.val) : nothing, nothing, nothing)
EnzymeRules.reverse(
    config::RevConfig, ::Const{typeof(_ad_pin)}, dret, tape,
    a::Annotation{Float64}, p::Annotation{Float64}
) = (_epoch_cotangent(a, dret, 1.0), _epoch_cotangent(p, dret, 0.0))

# The snapshot, the context and the correction are data without derivatives, and the
# lookups and the parsing of strings have no derivative.
EnzymeRules.inactive_type(::Type{<:KernelSet}) = true
EnzymeRules.inactive_type(::Type{<:SpiceContext}) = true
EnzymeRules.inactive_type(::Type{Aberration}) = true
for f in _AD_CONSTANT
    @eval EnzymeRules.inactive(::$(_ad_typeof(f)), args...) = nothing
end

end # module
