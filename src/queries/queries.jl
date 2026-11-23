# Prepared queries: a state or a rotation whose names, frame and aberration correction
# are resolved once, on one snapshot of the kernels, and then evaluated at any number
# of epochs.

"""
    StateQuery

The state of one target relative to one observer in one frame with one aberration
correction, resolved once by [`state_query`](@ref). Evaluate it with
[`spkezr`](@ref)`(q, et)` or [`spkpos`](@ref)`(q, et)`, or call it: `q(et)` is
`spkezr(q, et)`.

The fields `targ`, `obs` and `ref` hold the NAIF codes of the target and the observer
and the frame ID, and `corr` the [`Aberration`](@ref).
"""
struct StateQuery
    ks::KernelSet
    targ::Int32
    obs::Int32
    ref::Int32
    corr::Aberration
end

"""
    state_query([ctx,] targ, ref, abcorr, obs) -> StateQuery

Resolve the names of the target `targ` and the observer `obs`, the frame `ref` and the
aberration correction `abcorr` once, for a state that a program needs at many epochs.
The arguments are those of [`spkezr`](@ref) without the epoch, in the same order.

The query holds the kernels of `ctx` as they are at the call: a later
[`furnsh`](@ref), [`unload`](@ref) or pool edit of `ctx` does not change it. The query
is immutable, so one query can be used from any number of threads.

[`spkezr`](@ref)`(q, et)` and [`spkpos`](@ref)`(q, et)` return what `spkezr` and
`spkpos` return for the same arguments and kernels, to the last bit, without the
lookup of the names at each epoch. The epoch can be a `ForwardDiff.Dual`, as for
`spkezr`. This call has no CSPICE counterpart.

Throws `SpiceError`:

- `:idcodenotfound` if `targ` or `obs` is not a known body name;
- `:unknownframe` if `ref` is not a known frame;
- the kinds of [`Aberration`](@ref) for a malformed `abcorr`.

Data that is missing for an epoch throws at the evaluation, as `spkezr` does.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> add_ephemeris!(et -> [7000.0, 0.001et, 0, 0, 0.001, 0], ctx, -1000, 399, "J2000", 0.0, 100.0)

julia> q = state_query(ctx, -1000, "J2000", "NONE", 399);

julia> q(50.0)[1] == spkezr(ctx, -1000, 50.0, "J2000", "NONE", 399)[1]
true
```

See also [`rotation_query`](@ref).
"""
function state_query(ks::KernelSet, targ::NameOrID, ref::NameOrID, abcorr, obs::NameOrID)
    corr = abcorr isa Aberration ? abcorr : Aberration(abcorr)
    t, o, f = body_code(ks, targ), body_code(ks, obs), frame_id(ks, ref)
    frame_def(ks, f)                         # the frame must be defined
    return StateQuery(ks, t, o, f, corr)
end
state_query(ctx::SpiceContext, targ::NameOrID, ref::NameOrID, abcorr, obs::NameOrID) =
    state_query(kernelset(ctx), targ, ref, abcorr, obs)
state_query(targ::NameOrID, ref::NameOrID, abcorr, obs::NameOrID) =
    state_query(default_context(), targ, ref, abcorr, obs)

spkezr(q::StateQuery, et::Real) = spkezr(q.ks, q.targ, et, q.ref, q.corr, q.obs)
spkpos(q::StateQuery, et::Real) = spkpos(q.ks, q.targ, et, q.ref, q.corr, q.obs)
(q::StateQuery)(et::Real) = spkezr(q, et)
Base.broadcastable(q::StateQuery) = Ref(q)

"""
    RotationQuery

The orientation of one frame relative to another, resolved once by
[`rotation_query`](@ref). Evaluate it with [`pxform`](@ref)`(q, et)` or
[`sxform`](@ref)`(q, et)`, or call it: `q(et)` is `pxform(q, et)`. The fields `from`
and `to` hold the frame IDs.
"""
struct RotationQuery
    ks::KernelSet
    from::Int32
    to::Int32
end

"""
    rotation_query([ctx,] from, to) -> RotationQuery

Resolve the frames `from` and `to` once, for a rotation that a program needs at many
epochs. [`pxform`](@ref)`(q, et)` and [`sxform`](@ref)`(q, et)` return what `pxform`
and `sxform` return for the same frames and kernels, to the last bit.

As a [`state_query`](@ref), the query holds the kernels of `ctx` as they are at the
call, does not change, and can be used from any number of threads. This call has no
CSPICE counterpart.

Throws `SpiceError(:unknownframe)` if a frame is not known.

# Examples
```jldoctest
julia> q = rotation_query(SpiceContext(), "J2000", "ECLIPJ2000");

julia> q(0.0) == pxform(SpiceContext(), "J2000", "ECLIPJ2000", 0.0)
true
```

See also [`state_query`](@ref).
"""
function rotation_query(ks::KernelSet, from::NameOrID, to::NameOrID)
    f, t = frame_id(ks, from), frame_id(ks, to)
    frame_def(ks, f)
    frame_def(ks, t)
    return RotationQuery(ks, f, t)
end
rotation_query(ctx::SpiceContext, from::NameOrID, to::NameOrID) =
    rotation_query(kernelset(ctx), from, to)
rotation_query(from::NameOrID, to::NameOrID) = rotation_query(default_context(), from, to)

pxform(q::RotationQuery, et::Real) = pxform(q.ks, q.from, q.to, et)
sxform(q::RotationQuery, et::Real) = sxform(q.ks, q.from, q.to, et)
(q::RotationQuery)(et::Real) = pxform(q, et)
Base.broadcastable(q::RotationQuery) = Ref(q)
