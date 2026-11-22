module AlmanacTrackerExt

using Almanac
using Tracker: Tracker, TrackedReal, TrackedVector

# A `TrackedReal` epoch does not run the evaluators on tracked numbers. Each derivative
# point takes the jet at the plain epoch and assembles the tracked result from it (see
# `src/ad.jl`), so the tape holds a few operations for each query.

Almanac.primal(x::TrackedReal) = Almanac.primal(Tracker.data(x))
Almanac._ad_order(::Type{<:TrackedReal{V}}) where {V} = 1 + Almanac._ad_order(V)

Almanac.@_ad_epoch_methods TrackedReal

# An array of Tracker has the element type of its data (Float64) and returns tracked
# numbers, so an `Orbit` takes its elements one by one.
Almanac.Orbit(
    state::TrackedVector, et::Real, center::Integer, frame::Integer, gm::Real,
    radii::AbstractVector{<:Real} = Almanac._NO_RADII
) =
    Almanac.Orbit([state[i] for i in eachindex(state)], et, center, frame, gm, radii)

end # module
