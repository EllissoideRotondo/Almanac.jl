# The analysis layer: an `Orbit` is a state about a central body, with the GM and the
# radii of that body, and the quantities of `elements.jl` read from it. The formulas
# follow ANISE (src/astro/orbit.rs) at commit b3aaadc.

"""
    Orbit(state, et, center, frame, gm[, radii])
    Orbit([ctx,] state, et, center, frame)

A state about a central body, in km and km/s, with the central body's GM (km³/s²) and
radii (km). No CSPICE call has this interface; the elements of the orbit are those of
`oscltx`, `recgeo` and `conics` for the same state, and the type matches the `Orbit` of
ANISE.

`state` is `[x, y, z, vx, vy, vz]` in the frame `frame`, relative to `center`, at the
epoch `et` (seconds past J2000 TDB). The orbit is taken in the frame it is given; no
rotation is applied. The number type `T` of the orbit is the promoted type of `state`,
`et` and `Float64`, so every element has the one type `T`, and a dual number in either
flows through every element. `gm` and `radii` stay `Float64`.

The form with `gm` needs no kernels, and `center` and `frame` are NAIF IDs. `radii`
defaults to NaN, which is what an unknown radius is; the altitudes, `:latitude` and
`:height` then throw `SpiceError(:badarg)`. The form without `gm` takes `center` and
`frame` as names or NAIF IDs and reads `BODYnnn_GM` and, when present, `BODYnnn_RADII`
of the kernel pool (of `ctx`, or of the default context); see [`orbit`](@ref) for the
state from the kernels.

The fields are `state`, `et`, `center`, `frame`, `gm` and `radii`. The properties
`position`, `velocity`, `hvec` and `evec` (the eccentricity vector) are `SVector{3,T}`;
every symbol of [`elements`](@ref) is a property too (see [`element`](@ref)).

Throws `SpiceError(:badarg)` if `state` does not have 6 elements or `radii` 3, and, in
the form without `gm`, `SpiceError(:kernelvarnotfound)` if the pool has no GM for
`center`.

# Examples
```jldoctest
julia> o = Orbit([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0], 0.0, 399, 1, 398600.4415);

julia> o.rmag, o.hz
(7000.0, 52500.0)
```
"""
struct Orbit{T<:Real}
    state::SVector{6,T}
    et::T
    center::Int32
    frame::Int32
    gm::Float64
    radii::SVector{3,Float64}
end

const _NO_RADII = SVector(NaN, NaN, NaN)

function Orbit(state::AbstractVector{S}, et::E, center::Integer, frame::Integer, gm::Real,
        radii::AbstractVector{<:Real}=_NO_RADII) where {S<:Real,E<:Real}
    length(state) == 6 || throw(SpiceError(:badarg,
        "Orbit: the state has $(length(state)) elements, not 6"))
    length(radii) == 3 || throw(SpiceError(:badarg,
        "Orbit: the radii have $(length(radii)) elements, not 3"))
    T = promote_type(S, E, Float64)
    return Orbit{T}(SVector{6,T}(state), T(et), Int32(center), Int32(frame), Float64(gm),
        SVector{3,Float64}(radii))
end

function Orbit(ctx::SpiceContext, state::AbstractVector{<:Real}, et::Real, center::NameOrID,
        frame::NameOrID)
    ks = kernelset(ctx)
    c, f = body_code(ks, center), frame_id(ks, frame)
    gm = bodvcd(ctx, c, "GM")[1]
    radii = bodfnd(ctx, c, "RADII") ? SVector{3,Float64}(bodvcd(ctx, c, "RADII")) : _NO_RADII
    return Orbit(state, et, c, f, gm, radii)
end
Orbit(state::AbstractVector{<:Real}, et::Real, center::NameOrID, frame::NameOrID) =
    Orbit(default_context(), state, et, center, frame)

"""
    orbit([ctx,] targ, et, ref, abcorr, obs) -> Orbit

The [`Orbit`](@ref) of `targ` about `obs` at `et` (seconds past J2000 TDB): the state of
[`spkezr`](@ref) (same arguments, same order, with the correction `abcorr`) in the frame
`ref`, with the center `obs` and its GM and radii from the kernel pool. The epoch of the
orbit is `et`, the epoch of the observer, and the state is not rotated: it stays in
`ref`. No CSPICE call has this interface; it matches the construction of an `Orbit` from
an ephemeris state in ANISE.

Throws the errors of [`spkezr`](@ref), and `SpiceError(:kernelvarnotfound)` if the pool
has no GM for `obs`.
"""
function orbit(ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID,
        abcorr::AbstractString, obs::NameOrID)
    state, _ = spkezr(ctx, targ, et, ref, abcorr, obs)
    return Orbit(ctx, state, et, obs, ref)
end
orbit(targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString, obs::NameOrID) =
    orbit(default_context(), targ, et, ref, abcorr, obs)

const _ORBIT_VECTORS = (:position, :velocity, :hvec, :evec)

function Base.getproperty(o::Orbit, s::Symbol)
    hasfield(Orbit, s) && return getfield(o, s)
    s === :position && return SVector(getfield(o, :state)[1], getfield(o, :state)[2],
        getfield(o, :state)[3])
    s === :velocity && return SVector(getfield(o, :state)[4], getfield(o, :state)[5],
        getfield(o, :state)[6])
    s === :hvec && return _hvec(o)
    s === :evec && return _evec_for(:evec, o)
    s in elements() && return element(o, s)
    throw(SpiceError(:badarg, "Orbit: no field, vector or element :$s; the elements are " *
        join(elements(), ", ")))
end

Base.propertynames(::Orbit) = (fieldnames(Orbit)..., _ORBIT_VECTORS..., elements()...)

function Base.show(io::IO, o::Orbit)
    print(io, "Orbit(", o.state, ", ", o.et, ", ", o.center, ", ", o.frame, ", ", o.gm, ")")
end
