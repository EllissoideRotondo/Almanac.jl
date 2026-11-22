module AlmanacReverseDiffExt

using Almanac
using ReverseDiff: ReverseDiff, TrackedReal

# A `TrackedReal` epoch does not run the evaluators on tracked numbers. Each derivative
# point takes the jet at the plain epoch and assembles the tracked result from it (see
# `src/ad.jl`), so the tape holds a few operations for each query.

Almanac.primal(x::TrackedReal) = Almanac.primal(ReverseDiff.value(x))
Almanac._ad_order(::Type{<:TrackedReal{V}}) where {V} = 1 + Almanac._ad_order(V)

Almanac.@_ad_epoch_methods TrackedReal

end # module
