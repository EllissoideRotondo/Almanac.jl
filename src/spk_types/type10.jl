# SPK type 10: two-line element sets in a generic segment, propagated with SGP4/SDP4
# and turned from the TEME frame to J2000, as `spkr10.f` and `spke10.f` of CSPICE N0067
# with the generic-segment readers (`sgmeta.f`, `sgfcon.f`, `sgfrvi.f`, `sgfpkt.f`) and
# the frame routine `zzteme.f` (`zzeprc76.f`, `zzenut80.f`, `zzwahr.f`, `zzmobliq.f`,
# `eul2xf`, `zztwovxf.f`).
#
# `spke10` initializes the propagator for each of the two element sets at every call,
# so there is no saved state to reproduce: `_t10_state` does the same. The element
# epoch is turned to UTC with the leap seconds (`ttrans` in `xxsgp4i`), so this type is
# the only one that needs a leap-second kernel; the evaluator takes the `LeapSeconds`
# of the kernel set as an argument and `spk_eval`, which has no kernel set, cannot
# evaluate it.

@noinline _t10_error(seg::Segment, msg::String) = throw(
    SpiceError(
        :baddaf,
        "SPK type 10 segment at $(seg.baddr):$(seg.eaddr): $msg"
    )
)

# --- the generic segment ------------------------------------------------------------------

"""
    _t10_layout(data, seg) -> NamedTuple

Addresses of the parts of the generic segment `seg`, from its meta data (`sgmeta`):
`con` is the address before the first constant, `pkt` the address before the first
packet record, `ref` the address before the first reference epoch, `n` the count of
packets, `pksz` the packet size (10 or 14) and `stride` the size of a packet record
(the packet and the `pktoff` doubles that precede it).

The last double of the segment is the count of meta items. A segment with 15 or 16
items has no packet offset (item 16), as `sgmeta` reads the old files; with 17 or more
items the first 17 are used. Throws `SpiceError(:baddaf)` unless the segment has
fixed-size packets of 10 or 14 doubles, one reference epoch for each packet with the
reference index type 4 (closest epoch), at least 8 constants, and all its parts
inside the segment.
"""
function _t10_layout(data::Vector{Float64}, seg::Segment)
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    (b >= 1 && e <= length(data) && len >= 16) || _t10_error(seg, "the segment is too short")
    nmeta = _seg_count(seg, @inbounds(data[e]), 15, _MAXCOUNT, "meta item count")
    m, hasoff = _sg_meta_base(nmeta, e)      # meta item k is data[m + k]
    span = e - m                             # doubles from the first meta item to the end
    span <= len || _t10_error(seg, "the meta data is longer than the segment")
    item(k, lo, what) = _seg_count(seg, @inbounds(data[m + k]), lo, _MAXCOUNT, what)
    conbas = item(1, 0, "constants base")
    ncon = item(2, 0, "constants count")
    rdrtyp = item(5, 0, "reference index type")
    refbas = item(6, 0, "reference base")
    nref = item(7, 0, "reference count")
    pdrtyp = item(10, 0, "packet directory type")
    pktbas = item(11, 0, "packet base")
    npkt = item(12, 1, "packet count")
    pksz = item(15, 1, "packet size")
    pktoff = hasoff ? item(16, 0, "packet offset") : 0
    ncon >= 8 || _t10_error(seg, "the segment has $ncon constants, 8 are needed")
    pdrtyp == 0 || _t10_error(seg, "the packets do not have a fixed size")
    (pksz == 10 || pksz == 14) || _t10_error(seg, "the packet size $pksz is not 10 or 14")
    rdrtyp == 4 || _t10_error(seg, "the reference index type $rdrtyp is not 4")
    nref == npkt || _t10_error(seg, "the counts of packets and reference epochs differ")
    stride = pksz + pktoff
    data_len = len - span
    (
        conbas + 8 <= data_len && pktoff <= data_len && npkt <= (data_len - pktbas) ÷ stride &&
            pktbas <= data_len && refbas <= data_len && nref <= data_len - refbas
    ) ||
        _t10_error(seg, "a part of the segment is outside the segment")
    return (;
        con = b - 1 + conbas, pkt = b - 1 + pktbas + pktoff, ref = b - 1 + refbas, n = npkt,
        pksz, stride,
    )
end

"The ten elements of packet `i` of the segment with the layout `lay`."
@inline _t10_elems(data::Vector{Float64}, lay, i::Int) =
    ntuple(k -> @inbounds(data[lay.pkt + (i - 1) * lay.stride + k]), Val(10))

"Index of the reference epoch closest to `x` (`sgfrvi` with the index type 4); of two at the same distance the later one."
@inline function _t10_closest(data::Vector{Float64}, lay, x::Float64)
    n = lay.n
    i = _lstled(data, lay.ref, n, x)
    if i == 0
        i = 1
    elseif i < n
        @inbounds (data[lay.ref + i + 1] - x <= x - data[lay.ref + i]) && (i += 1)
    end
    return i
end

"""
    _t10_record(data, seg, x::Float64) -> (geophs, elems1, elems2)

The eight geophysical constants and the two element sets that `spkr10` selects for
the epoch `x`: the reference epoch closest to `x` is found (`sgfrvi` with the index
type 4; of two epochs at the same distance the later one), and the sets are the
ones that bracket `x`. Before the first epoch and after the last one both sets are
the end set. Postcondition: `elems1[10] <= elems2[10]`.
"""
function _t10_record(data::Vector{Float64}, seg::Segment, x::Float64)
    lay = _t10_layout(data, seg)
    n = lay.n
    i = _t10_closest(data, lay, x)
    value = @inbounds data[lay.ref + i]
    from, to = x <= value ? (max(i - 1, 1), i) : (i, min(i + 1, n))
    geophs = ntuple(k -> @inbounds(data[lay.con + k]), Val(8))
    return geophs, _t10_elems(data, lay, from), _t10_elems(data, lay, to)
end

# --- TEME to J2000 --------------------------------------------------------------------------

"""
    _t10_eul2xf(ang::NTuple{6,T}, a::Int, b::Int, c::Int) -> (r::Mat3{T}, dr::Mat3{T})

Rotation `r = [ang[1]]_a [ang[2]]_b [ang[3]]_c` and its time derivative `dr` for the
angle rates `ang[4:6]`, as `eul2xf`: `dr` is the product of the skew matrix of the
angular velocity and `r`. Requires `b != a` and `b != c` (the two uses here are 3-2-3
and 1-3-1).
"""
@inline function _t10_eul2xf(ang::NTuple{6, T}, a::Int, b::Int, c::Int) where {T <: Real}
    r = _m_mxm(_m_rotate(ang[1], a), _m_mxm(_m_rotate(ang[2], b), _m_rotate(ang[3], c)))
    l = 6 - a - b
    d = b == a % 3 + 1 ? 1.0 : -1.0
    sa, ca = sincos(ang[1])
    sb, cb = sincos(ang[2])
    u, v = a == c ? (cb, d * sb) : (-d * sb, cb)
    o1 = -d * ang[4] + (-d * u) * ang[6]
    o2 = (-d * ca) * ang[5] + (-sa * v) * ang[6]
    o3 = sa * ang[5] + (-d * ca * v) * ang[6]
    w = Mat3{T}(
        ntuple(Val(9)) do k
            i, j = (k - 1) % 3 + 1, (k - 1) ÷ 3 + 1
            (i == l && j == b) ? o1 : (i == b && j == l) ? -o1 :
                (i == a && j == l) ? o2 : (i == l && j == a) ? -o2 :
                (i == b && j == a) ? o3 : (i == a && j == b) ? -o3 : zero(T)
        end
    )
    return r, _m_mxm(w, r)
end

"""
    _t10_precession(et::T) -> (r, dr)

Rotation from J2000 to the mean equator and equinox of date by the IAU 1976
precession, and its time derivative, as `zzeprc76`.
"""
@inline function _t10_precession(et::T) where {T <: Real}
    zeta, z, theta = _dyn_prec_angles(et)
    dzeta, dz, dtheta = _dyn_prec_rates(et)
    return _t10_eul2xf((-z, theta, -zeta, -dz, dtheta, -dzeta), 3, 2, 3)
end

"""
    _t10_nutation(et::T) -> (r, dr)

Rotation from the mean equator and equinox of date to the true ones by the IAU 1980
nutation, and its time derivative, as `zzenut80` with the mean obliquity of
`zzmobliq`.
"""
@inline function _t10_nutation(et::T) where {T <: Real}
    psi, eps, dpsi, deps = _dyn_wahr_rates(et)
    mob, dmob = _dyn_mobliq_state(et)
    return _t10_eul2xf((-mob - eps, -psi, mob, -dmob - deps, -dpsi, dmob), 1, 3, 1)
end

"""
    _t10_teme(et::T) -> (x, y, z, dx, dy, dz)

Axes of the TEME frame (true equator, mean equinox of date) in J2000 at `et`, and
their time derivatives, as the columns of the matrix `TM2J` of `zzteme`: the Z axis is
the pole of the true equator, the X axis is in the plane of that pole and the mean
equinox, and the frame is made from the two vectors as `zztwovxf` makes it.
"""
@inline function _t10_teme(et::T) where {T <: Real}
    r1, dr1 = _t10_precession(et)
    r2, dr2 = _t10_nutation(et)
    xj, dxj = _m_row(r1, 1), _m_row(dr1, 1)           # mean equinox in J2000
    zp, zv = _m_row(r2, 3), _m_row(dr2, 3)            # true pole in the mean frame
    # The pole in J2000: the product with the inverse state matrix of the precession.
    zj = _m_mtxv(r1, zp)
    dzj = Vec3{T}(
        ntuple(
            i -> dr1[1, i] * zp[1] + dr1[2, i] * zp[2] + dr1[3, i] * zp[3] +
                r1[1, i] * zv[1] + r1[2, i] * zv[2] + r1[3, i] * zv[3], Val(3)
        )
    )
    # `zztwovxf_`: the f2c `dvhat_` and `ducrss_` (`_m_dvhat` with a division).
    sz = _m_state(zj, dzj)
    z, dz = _m_pv(_m_dvhat(sz, false))
    sy = _m_ducrss(sz, _m_state(xj, dxj))
    y, dy = _m_pv(sy)
    x, dx = _m_pv(_m_ducrss(sy, sz))
    return x, y, z, dx, dy, dz
end

# --- evaluation -----------------------------------------------------------------------------

"""
    _t10_epoch(l::LeapSeconds, et::Float64) -> Float64

The element epoch `et` (TDB seconds past J2000) in days from 1949 December 31 00:00
UTC, as `xxsgp4i` computes it: the Julian date UTC of `ttrans` minus 2433281.5.
"""
@inline _t10_epoch(l::LeapSeconds, et::Float64) =
    _ttrans(l, _TT_TDB, _TT_JDUTC, _tvec(et))[1] - 2433281.5

"""
    _t10_propagate(l, geophs, elems, et::T) -> (pos, vel)

TEME state at `et` of the element set `elems` with the constants `geophs`: the
initialization and the propagation of `evsgp4`.
"""
@inline function _t10_propagate(
        l::LeapSeconds, geophs::NTuple{8, Float64},
        elems::NTuple{10, Float64}, et::T
    ) where {T <: Real}
    r = _sgp4_init(geophs, elems, _t10_epoch(l, elems[10]), _SGP4_AFSPC)
    return _sgp4_eval(r, (et - elems[10]) / 60.0)
end

"""
    _t10_state(data, seg, et::T, l::LeapSeconds) -> (pos::Vec3{T}, vel::Vec3{T})

Position (km) and velocity (km/s) in J2000 of a type 10 segment at `et`, as `spkr10`
and `spke10`. Between two element epochs `t1 < t2` the TEME state is the weighted
mean of the propagations of the two sets, with the weight
`(1 + cos(pi (et - t1)/(t2 - t1)))/2` for the earlier set, and the velocity has the
term of the weight rate; outside the element epochs it is the propagation of the end
set. The state is then turned to J2000 with the state matrix of `_t10_teme`, so the
velocity has the term of the frame rate.

Contract: `seg` is a type 10 descriptor and `data` the array of its file; `l` holds
the leap seconds of the kernel set. Discrete decisions use `primal(et)`. With a plain
epoch the result is the CSPICE state to rounding level. Throws `SpiceError(:baddaf)`
for a segment that is not a valid type 10 segment and the errors of `_sgp4_init` and
`_sgp4_eval`. No allocation.
"""
function _t10_state(data::Vector{Float64}, seg::Segment, et::T, l::LeapSeconds) where {T <: Real}
    geophs, e1, e2 = _t10_record(data, seg, primal(et))
    t1, t2 = e1[10], e2[10]
    p1, v1 = _t10_propagate(l, geophs, e1, et)
    t1 == t2 && return _t10_to_j2000(et, p1, v1)
    p2, v2 = _t10_propagate(l, geophs, e2, et)
    denom = t2 - t1
    arg = (et - t1) * _PI / denom
    dargdt = _PI / denom
    s, c = sincos(arg)
    w = c * 0.5 + 0.5
    dwdt = s * -0.5 * dargdt
    return _t10_to_j2000(
        et, w * p1 + (1.0 - w) * p2,
        (w * v1 + (1.0 - w) * v2) + (dwdt * p1 + (-dwdt) * p2)
    )
end

"""
    _t10_to_j2000(et::T, p::Vec3{T}, v::Vec3{T}) -> (pos, vel)

The TEME state `(p, v)` at `et` in J2000: the product with the state matrix `TM2J`
of `zzteme`, whose lower left block is the rate of the rotation, with the sums in the
order of `mxvg`.
"""
@inline function _t10_to_j2000(et::T, p::Vec3{T}, v::Vec3{T}) where {T <: Real}
    x, y, z, dx, dy, dz = _t10_teme(et)
    pos = Vec3{T}(
        x[1] * p[1] + y[1] * p[2] + z[1] * p[3],
        x[2] * p[1] + y[2] * p[2] + z[2] * p[3], x[3] * p[1] + y[3] * p[2] + z[3] * p[3]
    )
    vel = Vec3{T}(
        dx[1] * p[1] + dy[1] * p[2] + dz[1] * p[3] + x[1] * v[1] + y[1] * v[2] + z[1] * v[3],
        dx[2] * p[1] + dy[2] * p[2] + dz[2] * p[3] + x[2] * v[1] + y[2] * v[2] + z[2] * v[3],
        dx[3] * p[1] + dy[3] * p[2] + dz[3] * p[3] + x[3] * v[1] + y[3] * v[2] + z[3] * v[3]
    )
    return pos, vel
end

"""
    _t10_jet(l, data, seg, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}

Jet of a type 10 segment in the convention of `spk_eval`: order 0 is the CSPICE
position, order 1 the CSPICE velocity, and order k >= 2 is derivative k-1 of the
CSPICE velocity function, from a Taylor epoch. No allocation.
"""
@inline _t10_jet(
    l::LeapSeconds, data::Vector{Float64}, seg::Segment, et::T,
    ::Val{0}
) where {T <: Real} = (_t10_state(data, seg, et, l)[1],)
@inline _t10_jet(
    l::LeapSeconds, data::Vector{Float64}, seg::Segment, et::T,
    ::Val{1}
) where {T <: Real} = _t10_state(data, seg, et, l)
@inline function _t10_jet(
        l::LeapSeconds, data::Vector{Float64}, seg::Segment, et::T,
        ::Val{N}
    ) where {T <: Real, N}
    p, v = _t10_state(data, seg, taylor_var(et, Val(N - 1)), l)
    return (_values(p), _vec_derivatives(v)...)
end

"""
    _t10_needs_lsk(seg) -> Union{}

The type 10 branch of `spk_eval`: throws `SpiceError(:notsupported)`, because
`spk_eval` has no kernel set and a type 10 segment needs the leap seconds. A query
through a kernel set (`spkezr`, `spkgeo`, `spkder`, ...) evaluates the segment with
`_t10_jet`.
"""
@noinline _t10_needs_lsk(seg::Segment) = throw(
    SpiceError(
        :notsupported,
        "SPK data type 10 needs the leap seconds of the kernel set: spk_eval cannot " *
            "evaluate the segment at $(seg.baddr):$(seg.eaddr); use spkezr, spkgeo or spkder"
    )
)
