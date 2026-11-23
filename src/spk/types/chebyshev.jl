# Chebyshev types with their own record layout: type 14 (unequal intervals, a generic
# segment) and type 20 (velocity series only, position by integration).

"""
    _spk14(data, seg, et::T) -> (pos, vel)

State of a type 14 segment, as `spkr14` and `spke14`. The segment is a generic
segment: constants (one, the number of coefficients of a series), packets of fixed
size (`MID`, `RADIUS`, three position series, three velocity series), the start
epochs of the packets as reference values, a directory of them, and the meta data
that give the place and the count of each part. The packet is the last one whose
start epoch is not after the epoch.

Only the layout that `spk14b` of CSPICE writes is accepted (fixed packet size,
explicit reference values with the index type "last less than or equal"); another
one throws `SpiceError(:notsupported)`.
"""
function _spk14(data::Vector{Float64}, seg::Segment, et::T) where {T <: Real}
    lay = _spk14_layout(data, seg)
    ncoef = lay.ncoef
    _spk14_recsize(ncoef)
    a = _spk14_packet(data, seg, lay, primal(et))     # address of MID
    @inbounds mid, radius = data[a], data[a + 1]
    tau = (et - mid) / radius
    (pos,) = _cheb3(data, a + 2, ncoef, tau, radius, Val(0))
    (vel,) = _cheb3(data, a + 2 + 3 * ncoef, ncoef, tau, radius, Val(0))
    return pos, vel
end

"""
    _spk14_recsize(ncoef)

Throws `SpiceError(:spkrectoolarge)` if a type 14 record (`6 ncoef + 3` doubles) does
not fit the record buffer of 198 doubles of `spkpvn`, that is for series of more than
32 coefficients.
"""
@inline function _spk14_recsize(ncoef::Int)
    recsiz = 6 * ncoef + 3
    recsiz > 198 && _spk14_toolarge(recsiz)
    return nothing
end
@noinline _spk14_toolarge(recsiz::Int) = throw(
    SpiceError(
        :spkrectoolarge,
        "storage for $recsiz double precision numbers is needed for an SPK data " *
            "record and only 198 locations were available"
    )
)

"""
    _spk14_packet(data, seg, lay, et::Float64) -> Int

Address of the first double (`MID`) of the packet of the SPK type 14 or PCK type 3
segment `seg` with the layout `lay` for `et`: the last packet whose reference epoch is
not after `et`, and the first packet for an epoch before them all.
"""
@inline function _spk14_packet(data::Vector{Float64}, seg::Segment, lay, et::Float64)
    k = clamp(_lstled(data, seg.baddr - 1 + lay.refbas, lay.nref, et), 1, lay.npkt)
    return seg.baddr + lay.pktbas + (k - 1) * lay.stride + lay.pktoff
end

"""
    _chbigr(data, a, ncof, scale, mid, radius, x::T) -> (p::T, itgrl::T)

Value at `x` of the Chebyshev series with the `ncof` coefficients
`data[a:a+ncof-1] * scale` on the interval `mid ± radius`, and its integral from
`mid` to `x`, by the recurrences of `chbigr`; also the kernel of the public `chbigr`.
Requires `ncof >= 1`, `radius > 0` and the coefficients inside `data` (not checked).
"""
@inline function _chbigr(
        data::AbstractVector, a::Int, ncof::Int, scale::Float64, mid::Real,
        radius::Real, x::T
    ) where {T}
    cp(i) = @inbounds data[a + i - 1] * scale       # coefficient i, 1-based
    degp = ncof - 1
    s = (x - mid) / radius
    s2 = s * 2.0
    # Coefficients of the integral: a_i = (c_{i-1} - c_{i+1}) / (2 (i-1)), with the
    # special cases of the ends.
    a2 = ncof >= 3 ? cp(1) - cp(3) * 0.5 : cp(1)
    adegp1 = degp >= 2 ? cp(degp) * 0.5 / degp : zero(a2)
    adegp2 = degp >= 1 ? cp(degp + 1) * 0.5 / (degp + 1) : zero(a2)
    f0 = zero(T) + (degp == 0 ? a2 : adegp2)
    f1 = zero(T)
    z0, z1 = (degp == 0 ? a2 : adegp2), zero(a2)    # the integral series at s = 0
    w0 = zero(T)
    w1 = zero(T)
    i = ncof
    # invariant: f, z and w hold the Clenshaw sums of the terms above i-1 of the
    # integral series at s, of the integral series at 0, and of the series itself
    while i > 1
        ai = i == 2 ? a2 : i < ncof ? (cp(i - 1) - cp(i + 1)) * 0.5 / (i - 1) : adegp1
        f0, f1 = ai + (s2 * f0 - f1), f0
        z0, z1 = ai - z1, z0
        w0, w1 = cp(i) + (s2 * w0 - w1), w0
        i -= 1
    end
    itgrl = z1 + s * f0 - f1
    p = cp(1) + (s * w0 - w1)
    return p, radius * itgrl
end

"""
    _spk20(data, seg, et::T) -> (pos, vel)

State of a type 20 segment, as `spkr20` and `spke20`. The segment ends with DSCALE,
TSCALE, INITJD, INITFR, INTLEN (days), RSIZE and N; `open_daf` checks RSIZE and N and
puts them in the descriptor. A record holds, for each component, the coefficients
of the velocity series (in units of DSCALE km per TSCALE s) and then the position at
the middle of the record interval (in units of DSCALE km). The velocity is the
series; the position is the middle position plus the integral of the series.
"""
function _spk20(data::Vector{Float64}, seg::Segment, et::T) where {T <: Real}
    b, e = seg.baddr, seg.eaddr
    recsiz, nrec = Int(seg.rsize), Int(seg.nrec)
    (nrec >= 1 && recsiz >= 1 && e - b + 1 == Int128(recsiz) * nrec + 7) ||
        _badsegment(seg, "the record size and count do not agree with the length")
    # `spke20` takes NCOF = (RSIZE - 3) / 3 and refuses a count below 1.
    ncof = (recsiz - 3) ÷ 3
    ncof < 1 && _segerror(:invalidcount, seg, "the coefficient count $ncof is not positive")
    nterms = ncof + 1
    3 * nterms == recsiz || _badsegment(seg, "the record size $recsiz is not 3 (NCOF + 1)")
    @inbounds dscale, tscale, initjd, initfr, intlen = data[e - 6], data[e - 5], data[e - 4],
        data[e - 3], data[e - 2]
    (
        isfinite(dscale) && isfinite(tscale) && tscale != 0.0 && isfinite(initjd) &&
            isfinite(initfr) && isfinite(intlen) && intlen > 0.0
    ) ||
        _badsegment(seg, "the scales or the interval values are not valid")
    init = (initjd - _JD2000 + initfr) * _SPD
    intrvl = intlen * _SPD
    recno = clamp(_trunc_index((primal(et) - init) / intrvl) + 1, 1, nrec)
    recbeg = (initjd - _JD2000 + (recno - 1) * intlen) * _SPD
    radius = intrvl / 2.0
    mid = recbeg + initfr * _SPD + radius
    r = b + (recno - 1) * recsiz
    vscale = dscale / tscale
    vx, ix = _chbigr(data, r, ncof, vscale, mid, radius, et)
    vy, iy = _chbigr(data, r + nterms, ncof, vscale, mid, radius, et)
    vz, iz = _chbigr(data, r + 2 * nterms, ncof, vscale, mid, radius, et)
    @inbounds px, py, pz = dscale * data[r + nterms - 1], dscale * data[r + 2 * nterms - 1],
        dscale * data[r + 3 * nterms - 1]
    return Vec3{T}(px + ix, py + iy, pz + iz), Vec3{T}(vx, vy, vz)
end
