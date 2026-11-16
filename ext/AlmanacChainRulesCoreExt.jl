module AlmanacChainRulesCoreExt

using Almanac
using Almanac: KernelSet, Aberration, spk_jet, frame_rotation, frame_jet, abcorr_jet,
    abcorr_position, abcorr_state, _ad_rule, _AD_CONSTANT, _ad_typeof
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

# The position of the epoch in the arguments of each derivative point.
const POINTS = (
    (spk_jet, :(Tuple{KernelSet,Int32,Real,Int32,Int32,Val}), 3),
    (frame_rotation, :(Tuple{KernelSet,Int32,Int32,Real}), 4),
    (frame_jet, :(Tuple{KernelSet,Int32,Int32,Real,Val}), 4),
    (abcorr_jet, :(Tuple{KernelSet,Int32,Real,Int32,Aberration,Int32,Val}), 3),
    (abcorr_position, :(Tuple{KernelSet,Int32,Real,Int32,Aberration,Int32}), 3),
    (abcorr_state, :(Tuple{KernelSet,Int32,Real,Int32,Aberration,Int32}), 3),
)

for (f, sig, iet) in POINTS
    types = eval(sig).parameters
    n = length(types)
    args = [Symbol(:a, i) for i in 1:n]
    typed = [:($(args[i])::$(types[i])) for i in 1:n]
    nots = [i == iet ? :d : :(NoTangent()) for i in 1:n]
    @eval begin
        function ChainRulesCore.frule((_, $([Symbol(:Δ, i) for i in 1:n]...)),
                ::typeof($f), $(typed...))
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

# Lookups, the snapshot of a context and the parsing of strings have no derivative.
for f in _AD_CONSTANT
    @eval begin
        ChainRulesCore.frule(ṫ, ::$(_ad_typeof(f)), args...; kw...) = ($f(args...; kw...), NoTangent())
        ChainRulesCore.rrule(::$(_ad_typeof(f)), args...; kw...) =
            ($f(args...; kw...), _ -> (NoTangent(), map(_ -> NoTangent(), args)...))
    end
end

end # module
