# The analysis layer: an `Orbit` is a state about a central body, with the GM and the
# radii of that body, and the quantities of `elements.jl` read from it. The formulas
# follow ANISE (src/astro/orbit.rs) at commit b3aaadc.

"""
    Orbit(state, et, center, frame, gm[, radii])
    Orbit([ctx,] state, et, center, frame)

The state `[x, y, z, vx, vy, vz]` (km, km/s) relative to `center` in `frame` at `et`
(seconds past J2000 TDB), with the GM (km³/s²) and radii (km) of `center`; the `Orbit` of
[ANISE](https://github.com/nyx-space/anise). The state is stored as given, in `frame`.
This call has no CSPICE counterpart.

With `gm`, no kernels are read and `center` and `frame` are NAIF IDs; `radii` defaults to
NaN (unknown). Without `gm`, `center` and `frame` are names or IDs, and `BODYnnn_GM` and,
when present, `BODYnnn_RADII` come from the kernel pool. [`orbit`](@ref) takes the state
from the kernels too.

The number type `T` of every element is the promoted type of `state`, `et` and
`Float64`; `gm` and `radii` stay `Float64`. Fields: `state`, `et`,
`center`, `frame`, `gm`, `radii`. Properties: `position`, `velocity`, `hvec`, `evec` (the
eccentricity vector), each an `SVector{3,T}`, and every name of [`elements`](@ref).

Throws `SpiceError(:badarg)` if the vector `state` does not have 6 elements or `radii` 3,
and, in the form without `gm`, `SpiceError(:kernelvarnotfound)` if the pool has no GM for
`center`. A degenerate state (zero position or velocity) is accepted; its elements throw
when read.

# Examples
```jldoctest
julia> o = Orbit([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0], 0.0, 399, 1, 398600.4415);

julia> o.rmag, o.hz
(7000.0, 52500.0)
```
"""
struct Orbit{T <: Real}
    state::SVector{6, T}
    et::T
    center::Int32
    frame::Int32
    gm::Float64
    radii::SVector{3, Float64}
end

const _NO_RADII = SVector(NaN, NaN, NaN)

function Orbit(
        state::AbstractVector{S}, et::E, center::Integer, frame::Integer, gm::Real,
        radii::AbstractVector{<:Real} = _NO_RADII
    ) where {S <: Real, E <: Real}
    length(state) == 6 || throw(
        SpiceError(
            :badarg,
            "Orbit: the state has $(length(state)) elements, not 6"
        )
    )
    length(radii) == 3 || throw(
        SpiceError(
            :badarg,
            "Orbit: the radii have $(length(radii)) elements, not 3"
        )
    )
    T = promote_type(S, E, Float64)
    return Orbit{T}(
        SVector{6, T}(state), T(et), Int32(center), Int32(frame), Float64(gm),
        SVector{3, Float64}(radii)
    )
end

function Orbit(
        ks::KernelSet, state::AbstractVector{<:Real}, et::Real, center::NameOrID,
        frame::NameOrID
    )
    c, f = body_code(ks, center), frame_id(ks, frame)
    gm = bodvcd(ks, c, "GM")[1]
    radii = bodfnd(ks, c, "RADII") ? SVector{3, Float64}(bodvcd(ks, c, "RADII")) : _NO_RADII
    return Orbit(state, et, c, f, gm, radii)
end
Orbit(
    ctx::SpiceContext, state::AbstractVector{<:Real}, et::Real, center::NameOrID,
    frame::NameOrID
) = Orbit(kernelset(ctx), state, et, center, frame)
Orbit(state::AbstractVector{<:Real}, et::Real, center::NameOrID, frame::NameOrID) =
    Orbit(default_context(), state, et, center, frame)

"""
    orbit([ctx,] targ, et, ref, abcorr, obs) -> Orbit

Return the [`Orbit`](@ref) of `targ` about `obs`: the state that [`spkezr`](@ref) returns
for the same arguments, at the epoch `et` of the observer, with the GM and radii of `obs`
from the kernel pool. This call has no CSPICE counterpart.

Throws the errors of [`spkezr`](@ref), and `SpiceError(:kernelvarnotfound)` if the pool
has no GM for `obs`.
"""
function orbit(
        ks::KernelSet, targ::NameOrID, et::Real, ref::NameOrID,
        abcorr::AbstractString, obs::NameOrID
    )
    state, _ = spkezr(ks, targ, et, ref, abcorr, obs)
    return Orbit(ks, state, et, obs, ref)
end
orbit(
    ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString,
    obs::NameOrID
) = orbit(kernelset(ctx), targ, et, ref, abcorr, obs)
orbit(targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString, obs::NameOrID) =
    orbit(default_context(), targ, et, ref, abcorr, obs)

const _ORBIT_VECTORS = (:position, :velocity, :hvec, :evec)

# The field names as a tuple of symbols: `hasfield(Orbit, s)` with a name known only at run
# time calls into the runtime (a `ccall`), which Mooncake cannot differentiate.
const _ORBIT_FIELDS = fieldnames(Orbit)

# Inlined, so that a field or a vector named in the code compiles to its own branch: the AD
# systems that differentiate the compiled code (Mooncake, Enzyme) would otherwise compile
# every branch of `element` for each `o.position`.
@inline function Base.getproperty(o::Orbit, s::Symbol)
    s in _ORBIT_FIELDS && return getfield(o, s)
    s === :position && return _orbit_position(o)
    s === :velocity && return _orbit_velocity(o)
    return _orbit_property(o, s)
end

@inline _orbit_position(o::Orbit) =
    SVector(getfield(o, :state)[1], getfield(o, :state)[2], getfield(o, :state)[3])
@inline _orbit_velocity(o::Orbit) =
    SVector(getfield(o, :state)[4], getfield(o, :state)[5], getfield(o, :state)[6])

# The properties that are not fields: the vectors and the elements. (A function of its own,
# so that the rule of ChainRulesCore for `getproperty` can call it.)
function _orbit_property(o::Orbit, s::Symbol)
    s === :position && return _orbit_position(o)
    s === :velocity && return _orbit_velocity(o)
    s === :hvec && return _hvec(o)
    s === :evec && return _evec_for(:evec, o)
    s in elements() && return element(o, s)
    throw(
        SpiceError(
            :badarg, "Orbit: no field, vector or element :$s; the elements are " *
                join(elements(), ", ")
        )
    )
end

Base.propertynames(::Orbit) = (fieldnames(Orbit)..., _ORBIT_VECTORS..., elements()...)

function Base.show(io::IO, o::Orbit)
    return print(io, "Orbit(", o.state, ", ", o.et, ", ", o.center, ", ", o.frame, ", ", o.gm, ")")
end
