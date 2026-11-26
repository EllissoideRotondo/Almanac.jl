module AlmanacChainRulesCoreExt

using Almanac
using Almanac: _ad_pin, _ad_rule, _AD_POINTS, _AD_CONSTANT, _ad_typeof
using ChainRulesCore: ChainRulesCore
using ChainRulesCore: AbstractZero, AbstractThunk, NoTangent, Tangent, unthunk, backing
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

# The rules of the derivative points (`_AD_POINTS`: the types of the arguments and the
# position of the epoch among them).
for (f, types, iet) in _AD_POINTS
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

# The value of one formula with the derivatives of another (see `src/ad.jl`).
ChainRulesCore.frule((_, Δa, _), ::typeof(_ad_pin), a::Real, p::Float64) = _ad_pin(a, p), Δa
ChainRulesCore.rrule(::typeof(_ad_pin), a::Real, p::Float64) =
    _ad_pin(a, p), ȳ -> (NoTangent(), unthunk(ȳ), NoTangent())

# Lookups, the snapshot of a context and the parsing of strings have no derivative.
for f in _AD_CONSTANT
    @eval begin
        ChainRulesCore.frule(ṫ, ::$(_ad_typeof(f)), args...; kw...) = ($f(args...; kw...), NoTangent())
        ChainRulesCore.rrule(::$(_ad_typeof(f)), args...; kw...) =
            ($f(args...; kw...), _ -> (NoTangent(), map(_ -> NoTangent(), args)...))
    end
end

end # module
