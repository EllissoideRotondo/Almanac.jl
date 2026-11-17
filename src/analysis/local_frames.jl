# Local orbit frames of an `Orbit` (RIC, VNC, RCN, TNW) and the states relative to them.
# The axes follow ANISE (src/astro/orbit.rs at b3aaadc, `dcm3x3_from_*_to_inertial`); the
# rate of the 6 by 6 matrix is the exact derivative along the two-body motion, where ANISE
# differences two states 1 ms apart.

const _LOCAL_FRAME_KINDS = (:ric, :vnc, :rcn, :tnw)

# The 3 by 3 matrix of the columns `a`, `b`, `c`. The element types can differ: the tracked
# numbers of ReverseDiff have a type for each origin.
@inline function _lf_cols(a::SVector{3}, b::SVector{3}, c::SVector{3})
    T = promote_type(eltype(a), eltype(b), eltype(c))
    return SMatrix{3, 3, T, 9}(a[1], a[2], a[3], b[1], b[2], b[3], c[1], c[2], c[3])
end

# The matrix of the local axes and its time derivative along the two-body motion. The
# angular momentum is constant there, so only the unit vectors of the radius and of the
# velocity turn: with `a = -gm r / |r|^3`, d(r/|r|)/dt = (v - r̂ (r̂.v)) / |r| and
# d(v/|v|)/dt = (a - v̂ (v̂.a)) / |v|. The checks are those of `hvec()` and `hmag()`
# in ANISE, which fail where h is zero (a rectilinear orbit, or a zero radius or velocity).
function _local_frame_pair(name::Symbol, o::Orbit, kind::Symbol)
    kind in _LOCAL_FRAME_KINDS || throw(
        SpiceError(
            :badarg,
            "$name: the kind is :$kind; it must be one of :ric, :vnc, :rcn, :tnw"
        )
    )
    # the checks of `hvec()` and `hmag()`, named for the call that was made
    _hvec_for(name, o)
    primal(norm(_hvec(o))) > 0 || _degenerate(name, "the angular momentum is zero")
    return _local_frame_axes(o, kind)
end

function _local_frame_axes(o::Orbit{T}, kind::Symbol) where {T}
    h = _hvec(o)
    hm = norm(h)
    r, v = o.position, o.velocity
    hh = h / hm
    z = zero(SVector{3, T})
    rm = norm(r)
    if kind === :ric || kind === :rcn
        rh = r / rm
        rd = (v - rh * dot(rh, v)) / rm
        return _lf_cols(rh, cross(hh, rh), hh), _lf_cols(rd, cross(hh, rd), z)
    end
    vm = norm(v)
    vh = v / vm
    a = (-getfield(o, :gm) / rm^3) * r
    vd = (a - vh * dot(vh, a)) / vm
    kind === :vnc && return _lf_cols(vh, hh, cross(vh, hh)), _lf_cols(vd, z, cross(vd, hh))
    return _lf_cols(vh, cross(hh, vh), hh), _lf_cols(vd, cross(hh, vd), z)
end

"""
    local_frame(o, kind) -> SMatrix{3,3}

Return the rotation from a local frame of the orbit `o` to the frame of `o`. No CSPICE
call has this interface; `Orbit::dcm3x3_from_ric_to_inertial` and its `vnc`, `rcn` and
`tnw` forms in ANISE.

`kind` is `:ric`, `:vnc`, `:rcn` or `:tnw`. With `r̂ = r/|r|`, `v̂ = v/|v|` and
`ĥ = r × v/|r × v|`, the columns are the axes of the local frame in the frame of `o`:

| `kind` | columns |
|:--|:--|
| `:ric` (radial, in-track, cross-track) | `r̂`, `ĥ × r̂`, `ĥ` |
| `:vnc` (velocity, normal, co-normal) | `v̂`, `ĥ`, `v̂ × ĥ` |
| `:rcn` (radial, cross, normal) | the same as `:ric`; only the axes are named differently |
| `:tnw` (tangential, normal, orbit normal) | `v̂`, `ĥ × v̂`, `ĥ` |

The frame of `o` is assumed inertial and is not checked. The element type is that of `o`.

Throws `SpiceError(:badarg)` for another `kind`, and `SpiceError(:degeneratecase)` where
the angular momentum is zero.

# Examples
```jldoctest
julia> o = Orbit([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0], 0.0, 399, 1, 398600.4415);

julia> local_frame(o, :ric) == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
true
```
"""
local_frame(o::Orbit, kind::Symbol) = _local_frame_pair(:local_frame, o, kind)[1]

"""
    local_frame_state(o, kind) -> SMatrix{6,6}

Return the rotation of states from a local frame of the orbit `o` to the frame of `o`,
`[R 0; Ṙ R]`, where `R` is [`local_frame`](@ref)`(o, kind)` and `Ṙ` (1/s) its exact rate
along the two-body motion of `o`. No CSPICE call has this interface;
`Orbit::dcm_from_ric_to_inertial` and its forms in ANISE, whose rate is a
[difference quotient](@ref "Beyond CSPICE").

Throws what [`local_frame`](@ref) throws.

# Examples
```jldoctest
julia> o = Orbit([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0], 0.0, 399, 1, 398600.4415);

julia> S = local_frame_state(o, :ric); S[5, 1], S[1, 4]
(0.0010714285714285715, 0.0)
```
"""
function local_frame_state(o::Orbit, kind::Symbol)
    R, Rd = _local_frame_pair(:local_frame_state, o, kind)
    T = eltype(R)
    z = zero(T)
    return SMatrix{6, 6, T, 36}(
        R[1, 1], R[2, 1], R[3, 1], Rd[1, 1], Rd[2, 1], Rd[3, 1],
        R[1, 2], R[2, 2], R[3, 2], Rd[1, 2], Rd[2, 2], Rd[3, 2],
        R[1, 3], R[2, 3], R[3, 3], Rd[1, 3], Rd[2, 3], Rd[3, 3],
        z, z, z, R[1, 1], R[2, 1], R[3, 1],
        z, z, z, R[1, 2], R[2, 2], R[3, 2],
        z, z, z, R[1, 3], R[2, 3], R[3, 3]
    )
end

"""
    relative_state(o, chief, kind) -> SVector{6}

Return the state (km, km/s) of the orbit `o` minus that of `chief`, in the local frame
`kind` of `chief`: `Rᵀ dr` and `Ṙᵀ dr + Rᵀ dv`, with `dr` and `dv` the differences of
position and velocity and `R` and `Ṙ` from [`local_frame_state`](@ref)`(chief, kind)`. The velocity is the rate of the offset seen in
the rotating frame, so it is not zero for equal inertial velocities. No CSPICE call has
this interface; `Orbit::ric_difference_unchecked` and its `vnc` and `tnw` forms in ANISE
(`:rcn` uses the `ric` one). The frame is assumed inertial and is not checked.

Throws `SpiceError(:badarg)` if the epoch, the center or the frame differ, or `kind` is
unknown, and `SpiceError(:degeneratecase)` where the angular momentum of `chief` is zero.

# Examples
```jldoctest
julia> chief = Orbit([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0], 0.0, 399, 1, 398600.4415);

julia> o = Orbit([7001.0, 0.0, 0.0, 0.0, 7.5, 0.0], 0.0, 399, 1, 398600.4415);

julia> round.(Tuple(relative_state(o, chief, :ric)); digits=6)
(1.0, 0.0, 0.0, 0.0, -0.001071, 0.0)
```
"""
function relative_state(o::Orbit{S}, chief::Orbit{U}, kind::Symbol) where {S, U}
    primal(o.et) == primal(chief.et) || throw(
        SpiceError(
            :badarg,
            "relative_state: the epochs differ ($(primal(o.et)) and $(primal(chief.et)))"
        )
    )
    o.center == chief.center || throw(
        SpiceError(
            :badarg,
            "relative_state: the centers differ ($(o.center) and $(chief.center))"
        )
    )
    o.frame == chief.frame || throw(
        SpiceError(
            :badarg,
            "relative_state: the frames differ ($(o.frame) and $(chief.frame))"
        )
    )
    R, Rd = _local_frame_pair(:relative_state, chief, kind)
    d = o.state - chief.state
    dr = SVector(d[1], d[2], d[3])
    dv = SVector(d[4], d[5], d[6])
    p = _m_mtxv(R, dr)
    v = _m_mtxv(Rd, dr) + _m_mtxv(R, dv)
    return SVector(p[1], p[2], p[3], v[1], v[2], v[3])
end
