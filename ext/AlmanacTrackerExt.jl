module AlmanacTrackerExt

using Almanac
using Almanac: KernelSet, Aberration
using Tracker: Tracker, TrackedReal, TrackedVector

# A `TrackedReal` epoch does not run the evaluators on tracked numbers. Each derivative
# point takes the jet at the plain epoch and assembles the tracked result from it (see
# `src/ad.jl`), so the tape holds a few operations for each query.

Almanac.primal(x::TrackedReal) = Almanac.primal(Tracker.data(x))
Almanac._ad_order(::Type{<:TrackedReal{V}}) where {V} = 1 + Almanac._ad_order(V)

Almanac.spk_jet(ks::KernelSet, targ::Int32, et::TrackedReal, ref::Int32, obs::Int32,
    v::Val) = Almanac._ad_spk_jet(ks, targ, et, ref, obs, v)
Almanac.frame_rotation(ks::KernelSet, from::Int32, to::Int32, et::TrackedReal) =
    Almanac._ad_frame_rotation(ks, from, to, et)
Almanac.frame_jet(ks::KernelSet, from::Int32, to::Int32, et::TrackedReal, v::Val) =
    Almanac._ad_frame_jet(ks, from, to, et, v)
Almanac.frame_jet(ks::KernelSet, from::Int32, to::Int32, et::TrackedReal, v::Val{0}) =
    Almanac._ad_frame_jet(ks, from, to, et, v)
Almanac.abcorr_jet(ks::KernelSet, targ::Int32, et::TrackedReal, ref::Int32,
    corr::Aberration, obs::Int32, v::Val) = Almanac._ad_abcorr_jet(ks, targ, et, ref, corr, obs, v)
Almanac.abcorr_position(ks::KernelSet, targ::Int32, et::TrackedReal, ref::Int32,
    corr::Aberration, obs::Int32) = Almanac._ad_abcorr_position(ks, targ, et, ref, corr, obs)
Almanac.abcorr_state(ks::KernelSet, targ::Int32, et::TrackedReal, ref::Int32,
    corr::Aberration, obs::Int32) = Almanac._ad_abcorr_state(ks, targ, et, ref, corr, obs)

# An array of Tracker has the element type of its data (Float64) and returns tracked
# numbers, so an `Orbit` takes its elements one by one.
Almanac.Orbit(state::TrackedVector, et::Real, center::Integer, frame::Integer, gm::Real,
    radii::AbstractVector{<:Real}=Almanac._NO_RADII) =
    Almanac.Orbit([state[i] for i in eachindex(state)], et, center, frame, gm, radii)

end # module
