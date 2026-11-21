# Segment search and evaluation of the fixed-interval Chebyshev data types
# (SPK types 2 and 3, binary PCK type 2), with time derivatives of any order.

"""
    find_segment(segs::Vector, et::Float64) -> Int

Index of the first segment of `segs` (which is in priority order) with
`t0 <= et <= t1`, or 0 if no segment covers `et`.
"""
function find_segment(segs::Vector, et::Float64)
    @inbounds for i in eachindex(segs)   # invariant: no segment before `i` covers `et`
        s = segs[i]
        (s.t0 <= et <= s.t1) && return i
    end
    return 0
end

"""
    find_segment(list::SegmentList, et::Float64) -> Int

The index of the first segment of `list` with `t0 <= et <= t1` (0 for none), as
`find_segment` of the vector `list.segs`, from a binary search of the index of
`list` when it has one.
"""
function find_segment(list::SegmentList, et::Float64)
    b = list.bounds
    isempty(b) && return find_segment(list.segs, et)
    lo, hi = 0, length(b) + 1
    @inbounds while hi - lo > 1            # invariant: b[lo] <= et < b[hi] (b[0] = -∞, b[end+1] = +∞)
        mid = (lo + hi) >>> 1
        b[mid] <= et ? (lo = mid) : (hi = mid)
    end
    lo == 0 && return 0                    # before the first bound, or `et` is NaN
    @inbounds b[lo] == et && return Int(list.win[2lo - 1])
    lo == length(b) && return 0
    return @inbounds Int(list.win[2lo])
end

"""
    _cheb_record(data, seg, et) -> Int

Address of the first double (MID) of the record of `seg` that covers `primal(et)`, as
`spkr02` selects it: the record index comes from INIT and INTLEN and is clamped to
`1:nrec`, so the segment end epoch falls in the last record. Throws
`SpiceError(:baddaf)` if the record is not inside `data`.
"""
@inline function _cheb_record(data::Vector{Float64}, seg::Segment, et)
    k = unsafe_trunc(Int, (primal(et) - seg.init) / seg.intlen)
    k = clamp(k, 0, Int(seg.nrec) - 1)
    addr = seg.baddr + k * Int(seg.rsize)
    (addr >= 1 && addr + Int(seg.rsize) - 1 <= length(data)) ||
        throw(SpiceError(:baddaf, "segment record is outside the data array"))
    return addr
end

"""
    _clenshaw_step(c, s2::T, b1::NTuple{M,T}, b2::NTuple{M,T}) -> NTuple{M,T}

One step of the Clenshaw recurrence with derivatives: from `b1 = b_{k+1}` and
`b2 = b_{k+2}` (each with derivatives of orders 0 to M-1 with respect to `s`) return
`b_k` for the coefficient `c`, where `s2 = 2s`:
`b_k^(m) = 2m b_{k+1}^(m-1) + 2s b_{k+1}^(m) - b_{k+2}^(m)`, plus `c` for `m = 0`.
Orders 0 and 1 use the operation order of `chbint_` of CSPICE.
"""
@inline _clenshaw_step(c::Real, s2::T, b1::NTuple{M, T}, b2::NTuple{M, T}) where {M, T} =
    ntuple(
    m -> m == 1 ? c + (s2 * b1[1] - b2[1]) :
        (2 * (m - 1)) * b1[m - 1] + b1[m] * s2 - b2[m], Val(M)
)

"""
    _clenshaw_end(c0, s::T, b1::NTuple{M,T}, b2::NTuple{M,T}) -> NTuple{M,T}

Value and derivatives (orders 0 to M-1 with respect to `s`) of the Chebyshev series
from the last two Clenshaw terms `b1 = b_1` and `b2 = b_2` and the coefficient `c0`:
`f^(m) = m b_1^(m-1) + s b_1^(m) - b_2^(m)`, plus `c0` for `m = 0`.
"""
@inline _clenshaw_end(c0::Real, s::T, b1::NTuple{M, T}, b2::NTuple{M, T}) where {M, T} =
    ntuple(
    m -> m == 1 ? c0 + (s * b1[1] - b2[1]) :
        (m - 1) * b1[m - 1] + s * b1[m] - b2[m], Val(M)
)

"""
    _cheb3(data, a, ncoef, s::T, radius, ::Val{N}) -> NTuple{N+1,Vec3{T}}

Sum of the three Chebyshev series whose coefficients are `data[a:a+ncoef-1]`,
`data[a+ncoef:a+2ncoef-1]` and `data[a+2ncoef:a+3ncoef-1]` at `s`, with derivatives
to order N; derivative `m` is divided `m` times by `radius` (`d s / d t = 1 / radius`).

Requires `ncoef >= 1` and `a + 3ncoef - 1 <= length(data)` (not checked). One pass
over the coefficients, highest first, with the Clenshaw recurrence. For a plain epoch
the value and the first derivative come from the same floating-point operations as in
`chbval_` and `chbint_` of CSPICE, so they are equal to the last bit; this matters for
an angle of many turns, where a difference of one unit in the last place is large
relative to the reduced angle.
"""
@inline function _cheb3(
        data::AbstractVector{<:Real}, a::Int, ncoef::Int, s::T, radius::Real,
        ::Val{N}
    ) where {T, N}
    M = Val(N + 1)
    z = ntuple(_ -> zero(T), M)
    s2 = s * 2
    x1, x2 = z, z
    y1, y2 = z, z
    z1, z2 = z, z
    # invariant: (x1, x2), (y1, y2), (z1, z2) are the Clenshaw terms b_{k+1}, b_{k+2}
    # of the three series, with derivatives to order N
    for k in (ncoef - 1):-1:1
        @inbounds cx, cy, cz = data[a + k], data[a + ncoef + k], data[a + 2 * ncoef + k]
        x1, x2 = _clenshaw_step(cx, s2, x1, x2), x1
        y1, y2 = _clenshaw_step(cy, s2, y1, y2), y1
        z1, z2 = _clenshaw_step(cz, s2, z1, z2), z1
    end
    @inbounds cx, cy, cz = data[a], data[a + ncoef], data[a + 2 * ncoef]
    fx = _clenshaw_end(cx, s, x1, x2)
    fy = _clenshaw_end(cy, s, y1, y2)
    fz = _clenshaw_end(cz, s, z1, z2)
    return ntuple(M) do m
        v = Vec3{T}(fx[m], fy[m], fz[m])
        for _ in 2:m      # invariant: `v` has been divided by `radius` once per pass
            v = v / radius
        end
        v
    end
end

"Jet of order N of the first three coefficient sets (`nsets` sets in a record) of the record of `seg` at `et`; `skip` sets are passed over first."
@inline function _cheb_eval(
        data::Vector{Float64}, seg::Segment, et::T, nsets::Int,
        skip::Int, v::Val{N}
    ) where {T, N}
    addr = _cheb_record(data, seg, et)
    ncoef = (Int(seg.rsize) - 2) ÷ nsets
    @inbounds mid, radius = data[addr], data[addr + 1]
    s = (et - mid) / radius
    return _cheb3(data, addr + 2 + skip * ncoef, ncoef, s, radius, v)
end

@inline _spk3(data, seg, et, ::Val{0}) = _cheb_eval(data, seg, et, 6, 0, Val(0))
@inline function _spk3(data, seg, et, ::Val{N}) where {N}
    p = _cheb_eval(data, seg, et, 6, 0, Val(0))
    v = _cheb_eval(data, seg, et, 6, 3, Val(N - 1))
    return (p[1], v...)
end

# `spkpvn` of CSPICE signals SPICE(SPKTYPENOTSUPP); a binary PCK of another type has no
# message of its own there.
@noinline _unsupported(what, dtype) =
    throw(
    SpiceError(
        what == "SPK" ? :spktypenotsupp : :notsupported,
        "$what data type $dtype is not supported"
    )
)

"""
    spk_eval(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}

Position (km) of `seg.target` relative to `seg.center` in frame `seg.frame` at `et`
and its time derivatives to order N; order 1 is the velocity that CSPICE returns.

Contract: `seg` is a descriptor from `open_daf` and `data` is the array of the same
file; `seg.t0 <= primal(et) <= seg.t1` (not checked: outside it the nearest record is
extrapolated). Type 2: all orders are derivatives of the position series. Type 3:
order 0 comes from the position series and order `k >= 1` is derivative `k-1` of the
velocity series. Types 1, 5, 8, 9, 12 to 15 and 17 to 21 (`spk_more.jl`): order
`k >= 2` is derivative `k-1` of the velocity function that CSPICE defines. Type 10
and unknown types throw `SpiceError(:spktypenotsupp, ...)`. No allocation.
"""
@noinline _empty_segment(seg) = throw(
    SpiceError(
        :baddaf,
        "the segment of body or frame $(seg.target) has no data (addresses $(seg.baddr):$(seg.eaddr))"
    )
)

@inline function spk_eval(data::Vector{Float64}, seg::Segment, et::T, v::Val{N}) where {T <: Real, N}
    seg.eaddr < seg.baddr && _empty_segment(seg)
    if seg.dtype == 2
        return _cheb_eval(data, seg, et, 3, 0, v)
    elseif seg.dtype == 3
        return _spk3(data, seg, et, v)
    else
        return _spk_more(data, seg, et, v)      # the other types: spk_more.jl
    end
end

"""
    pck_eval(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}

Euler angles (radians) of a binary PCK type 2, 3 or 20 segment at `et` and their time
derivatives to order N, as `pcke02` of CSPICE gives them (for type 3 the derivatives of
the angle series, not its angular velocity series; see `_pck03_state`).

Contract: `seg` is a PCK descriptor from `open_daf` and `data` is the array of the
same file; `seg.t0 <= primal(et) <= seg.t1` (not checked). The three components are
the stored 3-1-3 angles `(PHI, DELTA, W)` with `PHI = RA + pi/2` and
`DELTA = pi/2 - DEC`, so the matrix from the inertial frame `seg.frame` to the body
frame is `R3(W) * R1(DELTA) * R3(PHI)`. Other types throw
`SpiceError(:notsupported, ...)`. No allocation.
"""
@inline function pck_eval(data::Vector{Float64}, seg::Segment, et::T, v::Val{N}) where {T <: Real, N}
    seg.eaddr < seg.baddr && _empty_segment(seg)
    seg.dtype == 2 && return _cheb_eval(data, seg, et, 3, 0, v)
    # Type 20 stores the rates of the angles and the angles at the middle of each
    # record, as SPK type 20 stores velocities and positions (`pcke20_`).
    seg.dtype == 20 && return _spk_more(data, seg, et, v)
    seg.dtype == 3 && return _pck03_angles(data, seg, et, v)
    return _unsupported("binary PCK", seg.dtype)
end

# Binary PCK type 3 (`pckr03_`, `pcke03_`): the layout of SPK type 14, with six series in
# a packet: RA, DEC and W (degrees), then the angular velocity (degrees/s, components in
# the body frame), which is not the derivative of the angles.

"Address of the packet (its MID) of the PCK type 3 segment `seg` for `et`, and the coefficient count of a series; throws `SpiceError(:pckkrectoolarge)` as `pckmat_`."
@inline function _pck03_packet(data::Vector{Float64}, seg::Segment, et::Float64)
    lay = _spk14_layout(data, seg)
    recsiz = 6 * lay.ncoef + 2
    recsiz > 130 && throw(
        SpiceError(
            :pckkrectoolarge, "storage for $recsiz double precision numbers is needed " *
                "for a PCK data record and only 130 locations were available"
        )
    )
    k = clamp(_lstled(data, seg.baddr - 1 + lay.refbas, lay.nref, et), 1, lay.npkt)
    return seg.baddr + lay.pktbas + (k - 1) * lay.stride + lay.pktoff, lay.ncoef
end

"""
    _pck03_angles(data, seg, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}

The 3-1-3 angles `(PHI, DELTA, W)` (radians) of the PCK type 3 segment `seg` at `et`
with the arithmetic of `pcke03_` (`PHI = pi/2 + RA`, `DELTA = pi/2 - DEC`, W not reduced),
and their time derivatives to order N, which are the derivatives of the angle series.
"""
@inline function _pck03_angles(data::Vector{Float64}, seg::Segment, et::T, v::Val{N}) where {T <: Real, N}
    a, ncoef = _pck03_packet(data, seg, primal(et))
    @inbounds mid, radius = data[a], data[a + 1]
    jet = _cheb3(data, a + 2, ncoef, (et - mid) / radius, radius, v)
    return ntuple(Val(N + 1)) do m
        ra, dec, w = _RPD * jet[m][1], _RPD * jet[m][2], _RPD * jet[m][3]
        m == 1 ? Vec3{T}(_HALFPI + ra, _HALFPI - dec, w) : Vec3{T}(ra, -dec, w)
    end
end

"""
    _pck03_state(data, seg, et::Float64) -> (rot::Mat3, drot::Mat3)

The rotation from the reference frame of the PCK type 3 segment `seg` to the body
frame at `et` and its derivative, as `pcke03_` gives them: `eul2m_` of the angles, and
column j of the derivative is `-av × rot[:, j]`, with the angular velocity `av` of the
segment.
"""
function _pck03_state(data::Vector{Float64}, seg::Segment, et::Float64)
    a, ncoef = _pck03_packet(data, seg, et)
    @inbounds mid, radius = data[a], data[a + 1]
    s = (et - mid) / radius
    (ang,) = _cheb3(data, a + 2, ncoef, s, radius, Val(0))
    (av,) = _cheb3(data, a + 2 + 3ncoef, ncoef, s, radius, Val(0))
    ra, dec, w = _RPD * ang[1], _RPD * ang[2], _RPD * ang[3]
    rot = _eul313(w, _HALFPI - dec, _HALFPI + ra)
    mav = Vec3{Float64}(-(_RPD * av[1]), -(_RPD * av[2]), -(_RPD * av[3]))
    c(j) = _m_cross(mav, Vec3{Float64}(rot[1, j], rot[2, j], rot[3, j]))
    c1, c2, c3 = c(1), c(2), c(3)
    return rot, Mat3{Float64}(c1[1], c1[2], c1[3], c2[1], c2[2], c2[3], c3[1], c3[2], c3[3])
end
