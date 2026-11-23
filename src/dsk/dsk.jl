# Digital shape kernels, type 2 (triangular plates): the geometry of the ray-surface
# intercept, the volume elements of segments, and the selection of segments.
#
# The functions here are ports of the CSPICE routines that the intercept needs. They
# follow the order of the floating-point operations of CSPICE, so that the result is
# the same plate and the same point to rounding level. All of them are pure functions
# of their arguments and none allocates.
#
# The vector arithmetic is that of the `_m_` kernels (`math/kernels.jl`), in the
# operation order of SPICELIB.

const _DSK_DPMAX = floatmax(Float64)
const _DSK_ZAXIS = Vec3{Float64}(0.0, 0.0, 1.0)
const _DSK_ORIGIN = Vec3{Float64}(0.0, 0.0, 0.0)

# Limits of CSPICE that are part of the file format or of the intercept.
const _DSK_MAXVRT = 16000002
const _DSK_MAXPLT = 2 * (_DSK_MAXVRT - 2)
const _DSK_MAXVOX = 100000000
const _DSK_MAXCGR = 100000
const _DSK_MAXVXL = 50000                 # voxels of one ray that `xdda` accepts

"The tolerances of CSPICE at start, for calls on a `KernelSet` without a context."
const _DSK_DEFAULT_TOLS = DSKTolerances()

"`reclat`: radius, longitude and latitude of the rectangular point `p`."
@inline function _dsk_reclat(p::Vec3{Float64})
    big = max(max(abs(p[1]), abs(p[2])), abs(p[3]))
    big > 0.0 || return 0.0, 0.0, 0.0
    x = p[1] / big
    y = p[2] / big
    z = p[3] / big
    radius = big * sqrt(x * x + y * y + z * z)
    lat = atan(z, sqrt(x * x + y * y))
    lon = (p[1] == 0.0 && p[2] == 0.0) ? 0.0 : atan(p[2], p[1])
    return radius, lon, lat
end

@inline function _dsk_latrec(r::Float64, lon::Float64, lat::Float64)
    slon, clon = sincos(lon)
    slat, clat = sincos(lat)
    return Vec3{Float64}(r * clon * clat, r * slon * clat, r * slat)
end

@inline function _dsk_cylrec(r::Float64, lon::Float64, z::Float64)
    slon, clon = sincos(lon)
    return Vec3{Float64}(r * clon, r * slon, z)
end

@noinline _dsk_throw(kind::Symbol, msg::String) = throw(SpiceError(kind, msg))

"`surfnm`: the outward unit normal of the ellipsoid with the semi-axes `a`, `b`, `c` at `p`; requires positive axes."
@inline function _dsk_surfnm(a::Float64, b::Float64, c::Float64, p::Vec3{Float64})
    (a > 0.0 && b > 0.0 && c > 0.0) || _dsk_throw(
        :badaxislength,
        "the ellipsoid axes must be positive"
    )
    m = min(min(a, b), c)
    a1 = m / a
    b1 = m / b
    c1 = m / c
    return _m_hat(Vec3{Float64}(p[1] * (a1 * a1), p[2] * (b1 * b1), p[3] * (c1 * c1)))
end

"`georec`: the rectangular point of planetodetic longitude, latitude and altitude on the spheroid `re`, `f`."
@inline function _dsk_georec(lon::Float64, lat::Float64, alt::Float64, re::Float64, f::Float64)
    re > 0.0 || _dsk_throw(:valueoutofrange, "the equatorial radius must be positive")
    f < 1.0 || _dsk_throw(:valueoutofrange, "the flattening coefficient must be less than 1")
    rp = re - f * re
    sphi, cphi = sincos(lat)
    slmbda, clmbda = sincos(lon)
    big = max(abs(re * cphi), abs(rp * sphi))
    x = re * cphi / big
    y = rp * sphi / big
    scale = 1.0 / (big * sqrt(x * x + y * y))
    base = Vec3{Float64}(
        scale * re * re * clmbda * cphi, scale * re * re * slmbda * cphi,
        scale * rp * rp * sphi
    )
    normal = _dsk_surfnm(re, re, rp, base)
    return 1.0 * base + alt * normal
end

"""
    _dsk_surfpt(positn, u, a, b, c) -> (point, found)

`surfpt`: the first point of the ray from `positn` along `u` on the ellipsoid with the
semi-axes `a`, `b`, `c`. Requires `u != 0` and positive axes.
"""
@inline function _dsk_surfpt(
        positn::Vec3{Float64}, u::Vec3{Float64}, a::Float64, b::Float64,
        c::Float64
    )
    _m_iszero(u) && _dsk_throw(:zerovector, "the ray direction is the zero vector")
    (a > 0.0 && b > 0.0 && c > 0.0) || _dsk_throw(
        :badaxislength,
        "the ellipsoid axes must be positive"
    )
    x = Vec3{Float64}(u[1] / a, u[2] / b, u[3] / c)
    y = Vec3{Float64}(positn[1] / a, positn[2] / b, positn[3] / c)
    p = _m_perp(y, x, false)
    yproj = y - p
    ymag = _m_norm(y)
    pmag = _m_norm(p)
    ux = _m_hat(x)
    sign = 1.0
    if ymag > 1.0
        pmag > 1.0 && return _DSK_ORIGIN, false
        _m_dot(yproj, x) > 0.0 && return _DSK_ORIGIN, false
        pmag == 1.0 && return Vec3{Float64}(p[1] * a, p[2] * b, p[3] * c), true
        sign = -1.0
    elseif ymag == 1.0
        return positn, true
    end
    scale = sqrt(max(0.0, 1 - pmag * pmag))
    pt = 1.0 * p + (sign * scale) * ux
    return Vec3{Float64}(pt[1] * a, pt[2] * b, pt[3] * c), true
end

"`nplnpt`: the distance of `point` from the line through `linpt` along `lindir`; requires `lindir != 0`."
@inline function _dsk_nplnpt_dist(linpt::Vec3{Float64}, lindir::Vec3{Float64}, point::Vec3{Float64})
    _m_iszero(lindir) && _dsk_throw(:zerovector, "the direction of the line is the zero vector")
    trans = point - linpt
    pnear = _m_proj(trans, lindir, false) + linpt
    return _m_norm(pnear - point)
end

# ---------------------------------------------------------------------------------------
# Planes, cones and quadratics (for the latitude boundaries of volume elements)
# ---------------------------------------------------------------------------------------

"`nvp2pl`: the plane (unit normal, constant >= 0) through `point` with the normal `normal`."
@inline function _dsk_nvp2pl(normal::Vec3{Float64}, point::Vec3{Float64})
    n = _m_hat(normal)
    c = _m_dot(point, n)
    c < 0.0 && return -n, -c
    return n, c
end

"`nvc2pl`: the plane (unit normal, constant >= 0) of `normal` and `konst`; requires `normal != 0`."
@inline function _dsk_nvc2pl(normal::Vec3{Float64}, konst::Float64)
    n, mag = _m_unorm(normal)
    mag == 0.0 && _dsk_throw(:zerovector, "the normal vector of the plane is the zero vector")
    c = konst / mag
    c < 0.0 && return -n, -c
    return n, c
end

"""
    _dsk_inrypl(vertex, dir, normal, konst) -> (nxpts, xpt)

`inrypl`: the intersection of the ray from `vertex` along `dir` with the plane
`normal . x = konst`. `nxpts` is 1, 0, or -1 if the ray lies in the plane.
"""
@inline function _dsk_inrypl(
        vertex::Vec3{Float64}, dir::Vec3{Float64}, normal::Vec3{Float64},
        konst::Float64
    )
    toobig = _DSK_DPMAX / 3.0
    vn = _m_norm(vertex)
    (vn >= toobig || konst >= toobig) && _dsk_throw(
        :vectortoobig,
        "the ray vertex or the plane constant is too large"
    )
    udir = _m_hat(dir)
    _m_iszero(udir) && _dsk_throw(:zerovector, "the ray direction is the zero vector")
    mscale = max(konst, vn)
    sclvtx = vertex
    sclcon = konst
    if mscale != 0.0
        sclvtx = (1.0 / mscale) * vertex
        sclcon = konst / mscale
    end
    mscale > 1.0 && (toobig /= mscale)
    prjvn = _m_dot(sclvtx, normal)
    prjdif = sclcon - prjvn
    if prjdif == 0.0
        return (_m_dot(normal, udir) == 0.0 ? -1 : 1), vertex
    end
    prjdir = _m_dot(udir, normal)
    same = (prjdir < 0.0 && prjdif < 0.0) || (prjdir > 0.0 && prjdif > 0.0)
    same || return 0, _DSK_ORIGIN
    abs(prjdif) >= abs(prjdir) * toobig && return 0, _DSK_ORIGIN
    scale = abs(prjdif) / abs(prjdir)
    return 1, mscale * (1.0 * sclvtx + scale * udir)
end

"""
    _dsk_frame(x) -> (x, y, z)

`frame`: a right-handed orthonormal frame whose first axis is the unit vector of `x`.
"""
@inline function _dsk_frame(xin::Vec3{Float64})
    x = _m_hat(xin)
    a = x[1] * x[1]
    b = x[2] * x[2]
    c = x[3] * x[3]
    if a + b + c == 0.0
        return Vec3{Float64}(1.0, 0.0, 0.0), Vec3{Float64}(0.0, 1.0, 0.0), _DSK_ZAXIS
    end
    if a <= b && a <= c
        f = sqrt(b + c)
        y2 = -x[3] / f
        y3 = x[2] / f
        return x, Vec3{Float64}(0.0, y2, y3), Vec3{Float64}(f, -x[1] * y3, x[1] * y2)
    elseif b <= a && b <= c
        f = sqrt(a + c)
        y3 = -x[1] / f                   # s1, s2, s3 = 2, 3, 1
        y1 = x[3] / f
        return x, Vec3{Float64}(y1, 0.0, y3), Vec3{Float64}(x[2] * y3, f, -x[2] * y1)
    end
    f = sqrt(a + b)
    y1 = -x[2] / f                       # s1, s2, s3 = 3, 1, 2
    y2 = x[1] / f
    return x, Vec3{Float64}(y1, y2, 0.0), Vec3{Float64}(-x[3] * y2, x[3] * y1, f)
end

"""
    _dsk_bquad(a, b, c, ub) -> (n, nx, r1, r2)

`zzbquad`: the real roots of `a x^2 + b x + c` with magnitude at most `ub`; `n` roots
are returned and `nx` are excluded. Requires `ub > 0` and coefficients below
`sqrt(dpmax)/100`.
"""
function _dsk_bquad(a::Float64, b::Float64, c::Float64, ub::Float64)
    big = sqrt(_DSK_DPMAX) / 100
    (abs(a) > big || abs(b) > big || abs(c) > big || abs(ub) > big || ub <= 0.0) &&
        _dsk_throw(:valueoutofrange, "coefficients or root bound of the quadratic are out of range")
    if a == 0.0
        b == 0.0 && return (c == 0.0 ? -1 : -2), 0, 0.0, 0.0
        if abs(c) <= abs(ub * b)
            r = -c / b
            return 1, 0, r, r
        end
        return 0, 1, 0.0, 0.0
    end
    dscrim = b * b - a * 4 * c
    dscrim < 0.0 && return 0, 0, 0.0, 0.0
    if dscrim == 0.0
        num1 = -b
        denom = a * 2
        abs(num1) >= abs(denom * ub) && return 0, 1, 0.0, 0.0
        r = num1 / a / 2
        return 1, 0, r, r
    end
    denom = a * 2
    sqdisc = sqrt(dscrim)
    num2, num1 = b > 0.0 ? (-b - sqdisc, -b + sqdisc) : (-b + sqdisc, -b - sqdisc)
    if abs(num2) <= abs(ub * denom)
        r2 = num2 / a / 2
        r1 = abs(r2) > 0.0 ? c / (a * r2) : 0.0
        return 2, 0, r1, r2
    end
    abs(num1) <= abs(ub * denom) && return 1, 1, num1 / a / 2, 0.0
    return 0, 2, 0.0, 0.0
end

"""
    _dsk_cnquad(a, b, c, ub) -> (n, r1, r2)

`zzcnquad`: the roots of the quadratic of the cone intersection with magnitude at most
`ub`, in a form that stays stable when `a` is small.
"""
function _dsk_cnquad(a::Float64, b::Float64, c::Float64, ub::Float64)
    invub = sqrt(_DSK_DPMAX) / 200.0
    if a == 0.0 && b == 0.0
        return (c == 0.0 ? -1 : -2), 0.0, 0.0
    end
    maxmag = max(max(abs(a), abs(b)), abs(c))
    c1 = a / maxmag
    c2 = b / maxmag
    c3 = c / maxmag
    maxix = 1
    abs(c2) > abs(c1) && (maxix = 2)
    abs(c3) > abs(maxix == 1 ? c1 : c2) && (maxix = 3)
    if maxix == 1
        c1 = c1 >= 0.0 ? 1.0 : -1.0
    elseif maxix == 2
        c2 = c2 >= 0.0 ? 1.0 : -1.0
    else
        c3 = c3 >= 0.0 ? 1.0 : -1.0
    end
    if abs(c1) >= 1.0e-8 || c1 == 0.0
        n, _, r1, r2 = _dsk_bquad(c1, c2, c3, ub)
        return n, r1, r2
    elseif abs(c3) >= 1.0e-8
        n, _, inv1, inv2 = _dsk_bquad(c3, c2, c1, invub)
        if n == 1
            abs(inv1 * ub) >= 1.0 && return 1, 1.0 / inv1, 0.0
            return 0, 0.0, 0.0
        elseif n == 2
            if abs(inv2 * ub) >= 1.0
                r1 = 1.0 / inv2
                abs(inv1 * ub) >= 1.0 && return 2, r1, 1.0 / inv1
                return 1, r1, 0.0
            end
            return 0, 0.0, 0.0
        end
        return 0, 0.0, 0.0
    end
    if b < 0.0
        c1 = -c1
        c2 = -c2
        c3 = -c3
    end
    if abs(c3) <= ub
        r1 = -c3
        abs(c1 * c3 - 1.0) < abs(c1 * ub) && return 2, r1, c3 - 1.0 / c1
        return 1, r1, 0.0
    end
    return 0, 0.0, 0.0
end

"""
    _dsk_cxbrut(apex, axis, angle, endpt1, endpt2) -> (xpt, isbrck)

`zzcxbrut`: if exactly one end of the segment is inside the cone, the crossing point
found by bisection. `xpt` has no meaning when `isbrck` is false. Requires `axis != 0`.
"""
function _dsk_cxbrut(
        apex::Vec3{Float64}, axis::Vec3{Float64}, angle::Float64,
        endpt1::Vec3{Float64}, endpt2::Vec3{Float64}
    )
    _m_iszero(axis) && _dsk_throw(:zerovector, "the cone axis is the zero vector")
    locang = angle
    locaxi = axis
    if angle > _HALFPI
        locang = _PI - angle
        locaxi = -axis
    end
    locaxi = _m_hat(locaxi)
    cosang = cos(locang)
    off1 = endpt1 - apex
    off2 = endpt2 - apex
    state1 = _m_dot(_m_hat(off1), locaxi) >= cosang
    state2 = _m_dot(_m_hat(off2), locaxi) >= cosang
    state1 != state2 || return _DSK_ORIGIN, false
    seg = off2 - off1
    low = 0.0
    high = 1.0
    delta = abs(high - low)
    prvdlt = 2.0
    nitr = 0
    x = off1
    while delta > 1.0e-15 && delta < prvdlt && nitr < 1000
        # invariant: the crossing is between the parameters low and high
        midpt = (low + high) / 2
        x = 1.0 * off1 + midpt * seg
        state = _m_dot(_m_hat(x), locaxi) >= cosang
        if state == state1
            low = midpt
        else
            high = midpt
        end
        prvdlt = delta
        delta = abs(high - low)
        nitr += 1
    end
    return apex + x, true
end

"""
    _dsk_sglatx(p1, p2) -> (minlat, minp, maxlat, maxp)

`zzsglatx`: the extreme planetocentric latitudes on the segment from `p1` to `p2` and
the points where they occur.
"""
function _dsk_sglatx(p1::Vec3{Float64}, p2::Vec3{Float64})
    lat1 = _dsk_reclat(p1)[3]
    lat2 = _dsk_reclat(p2)[3]
    minlat, maxlat, minp, maxp = lat1 <= lat2 ? (lat1, lat2, p1, p2) : (lat2, lat1, p2, p1)
    normal = _m_cross(p1, p2)
    _m_iszero(normal) && return minlat, minp, maxlat, maxp
    normal = _m_hat(normal)
    crease = _m_cross(_DSK_ZAXIS, normal)
    _m_iszero(crease) && return minlat, minp, maxlat, maxp
    crease = _m_hat(crease)
    dp1 = _m_dot(p1, crease)
    dp2 = _m_dot(p2, crease)
    if (dp1 < 0.0 && dp2 > 0.0) || (dp1 > 0.0 && dp2 < 0.0)
        pn, pc = _dsk_nvc2pl(crease, 0.0)
        nxpts, t = _dsk_inrypl(p1, p2 - p1, pn, pc)
        if nxpts == 1
            lat = _dsk_reclat(t)[3]
            if lat > maxlat
                maxlat = lat
                maxp = t
            elseif lat < minlat
                minlat = lat
                minp = t
            end
        end
    end
    return minlat, minp, maxlat, maxp
end

"""
    _dsk_incnsg(apex, axis, angle, endpt1, endpt2) -> (nxpts, xpt1, xpt2)

`incnsg`: the points where the segment from `endpt1` to `endpt2` meets the cone with
the given apex, axis and half angle. Requires `axis != 0`, `angle >= 0` and distinct
end points. `xpt1` and `xpt2` have a meaning only up to `nxpts`.
"""
function _dsk_incnsg(
        apex::Vec3{Float64}, axis::Vec3{Float64}, angle::Float64,
        endpt1::Vec3{Float64}, endpt2::Vec3{Float64}
    )
    xpt1 = _DSK_ORIGIN
    xpt2 = _DSK_ORIGIN
    uaxis, axmag = _m_unorm(axis)
    axmag == 0.0 && _dsk_throw(:zerovector, "the cone axis is the zero vector")
    angle < 0.0 && _dsk_throw(:invalidangle, "the cone angle must be non-negative")
    udir, dmag = _m_unorm(endpt2 - endpt1)
    dmag == 0.0 && _dsk_throw(:endpointsmatch, "the end points of the segment coincide")
    cosang = cos(angle)
    locang = angle
    if cosang < 0.0
        cosang = -cosang
        locang = _PI - angle
        uaxis = -uaxis
    end
    off1 = endpt1 - apex
    off2 = endpt2 - apex
    dp1 = _m_dot(_m_hat(off1), uaxis)
    dp2 = _m_dot(_m_hat(off2), uaxis)
    in1 = dp1 >= cosang
    in2 = dp2 >= cosang
    ((in1 && in2) || (dp1 < 0.0 && dp2 < 0.0)) && return 0, xpt1, xpt2
    if abs(cosang) < 1.0e-14
        pn, pc = _dsk_nvp2pl(uaxis, apex)
        nplnx, plnx = _dsk_inrypl(endpt1, udir, pn, pc)
        (nplnx == 1 && _m_norm(plnx - endpt1) <= dmag) && return 1, plnx, xpt2
        return 0, xpt1, xpt2
    end
    uuax = _m_dot(udir, uaxis)
    wuax = _m_dot(off1, uaxis)
    wu = _m_dot(off1, udir)
    w2 = _m_dot(off1, off1)
    ca2 = cosang * cosang
    a = uuax * uuax - ca2
    b = (wuax * uuax - ca2 * wu) * 2
    c = wuax * wuax - w2 * ca2
    n, s1, s2 = _dsk_cnquad(a, b, c, dmag)
    nxpts = 0
    if n > 0
        if s1 >= 0.0
            xpt1 = Vec3{Float64}(
                endpt1[1] + s1 * udir[1], endpt1[2] + s1 * udir[2],
                endpt1[3] + s1 * udir[3]
            )
            uv1, v1mag = _m_unorm(xpt1 - apex)
            coserr = v1mag > 0.0 ? abs(_m_dot(uv1, uaxis) - cosang) : 0.0
            (v1mag == 0.0 || coserr < 1.0e-10) && (nxpts = 1)
        end
        if n == 2 && s2 >= 0.0
            xpt2 = Vec3{Float64}(
                endpt1[1] + s2 * udir[1], endpt1[2] + s2 * udir[2],
                endpt1[3] + s2 * udir[3]
            )
            uv2, v2mag = _m_unorm(xpt2 - apex)
            coserr = v2mag > 0.0 ? abs(_m_dot(uv2, uaxis) - cosang) : 0.0
            if v2mag == 0.0 || coserr < 1.0e-10
                nxpts += 1
                nxpts == 1 && (xpt1 = xpt2)
            end
        end
    end
    nxpts < 2 || return nxpts, xpt1, xpt2
    if (in1 && !in2) || (in2 && !in1)
        if nxpts == 0
            x, isbrck = _dsk_cxbrut(apex, uaxis, locang, endpt1, endpt2)
            isbrck && return 1, x, xpt2
        end
        return nxpts, xpt1, xpt2
    end
    # Both end points are outside the cone: the segment can cross the cone twice.
    nxpts = 0
    uaxis, fx, fy = _dsk_frame(uaxis)
    xform = Mat3{Float64}(fx[1], fy[1], uaxis[1], fx[2], fy[2], uaxis[2], fx[3], fy[3], uaxis[3])
    xoff1 = _m_mxv(xform, off1)
    xoff2 = _m_mxv(xform, off2)
    _, _, maxlat, maxp = _dsk_sglatx(xoff1, xoff2)
    colat = _HALFPI - maxlat
    if colat < locang
        vtemp, isbrck = _dsk_cxbrut(_DSK_ORIGIN, _DSK_ZAXIS, locang, xoff1, maxp)
        if isbrck
            xpt1 = _m_mtxv(xform, vtemp) + apex
            nxpts = 1
        end
        vtemp, isbrck = _dsk_cxbrut(_DSK_ORIGIN, _DSK_ZAXIS, locang, maxp, xoff2)
        if isbrck
            xpt2 = _m_mtxv(xform, vtemp) + apex
            if nxpts == 1
                nxpts = 2
            else
                xpt1 = xpt2
                nxpts = 1
            end
        end
    elseif colat == locang
        xpt1 = maxp                      # as CSPICE: the point stays in the cone frame
        nxpts = 1
    end
    return nxpts, xpt1, xpt2
end

"`zzryxsph`: the first point of the ray (`udir` a unit vector) on the sphere of radius `r` about the origin."
@inline function _dsk_ryxsph(vertex::Vec3{Float64}, udir::Vec3{Float64}, r::Float64)
    cpar = _m_dot(vertex, udir)
    perp = 1.0 * vertex + (-cpar) * udir
    pmag2 = _m_dot(perp, perp)
    r2 = r * r
    pmag2 > r2 && return _DSK_ORIGIN, false
    s = sqrt(max(0.0, r2 - pmag2))
    vmag2 = _m_dot(vertex, vertex)
    if vmag2 > r2
        cpar > 0.0 && return _DSK_ORIGIN, false
        return Vec3{Float64}(perp[1] - s * udir[1], perp[2] - s * udir[2], perp[3] - s * udir[3]), true
    elseif vmag2 < r2
        return Vec3{Float64}(perp[1] + s * udir[1], perp[2] + s * udir[2], perp[3] + s * udir[3]), true
    end
    return perp, true
end

"`zzinrypl`: the point of the ray on the plane `uplnml . x = konst`, if it is within `maxd` of the vertex."
@inline function _dsk_zzinrypl(
        vertex::Vec3{Float64}, udir::Vec3{Float64}, uplnml::Vec3{Float64},
        konst::Float64, maxd::Float64
    )
    vtxcon = _m_dot(vertex, uplnml)
    dircon = _m_dot(udir, uplnml)
    (vtxcon > konst && dircon > 0.0) && return 0, _DSK_ORIGIN
    (vtxcon < konst && dircon < 0.0) && return 0, _DSK_ORIGIN
    if vtxcon == konst
        dircon != 0.0 && return 1, vertex
        return 0, _DSK_ORIGIN
    end
    h = abs(vtxcon - konst)
    lpar = abs(dircon)
    h > maxd * lpar && return 0, _DSK_ORIGIN
    lpar == 0.0 && return 0, _DSK_ORIGIN
    s = h / lpar
    return 1, Vec3{Float64}(vertex[1] + s * udir[1], vertex[2] + s * udir[2], vertex[3] + s * udir[3])
end

# ---------------------------------------------------------------------------------------
# Spheroids
# ---------------------------------------------------------------------------------------

"`ednmpt`: the point of the ellipsoid `a`, `b`, `c` whose outward normal is parallel to `normal`."
@inline function _dsk_ednmpt(a::Float64, b::Float64, c::Float64, normal::Vec3{Float64})
    (a > 0.0 && b > 0.0 && c > 0.0) || _dsk_throw(
        :badaxislength,
        "all ellipsoid semi-axis lengths must be strictly positive"
    )
    scale = max(max(a, b), c)
    sa = a / scale
    sb = b / scale
    sc = c / scale
    (sa > 0.0 && sb > 0.0 && sc > 0.0) || _dsk_throw(
        :axisunderflow,
        "a scaled semi-axis length underflowed"
    )
    _m_iszero(normal) && _dsk_throw(:zerovector, "the normal vector is the zero vector")
    na2 = normal[1] * sa * sa
    nb2 = normal[2] * sb * sb
    nc2 = normal[3] * sc * sc
    arg = na2 * normal[1] + nb2 * normal[2] + nc2 * normal[3]
    arg > 0.0 || _dsk_throw(:degeneratecase, "the scale factor of the normal is not positive")
    lambda = arg^(-0.5) * scale
    return Vec3{Float64}(lambda * na2, lambda * nb2, lambda * nc2)
end

"`zzelnaxx`: the x and y intercepts of the normal line of the ellipse `a`, `b` at planetodetic latitude `lat`."
@inline function _dsk_elnaxx(a::Float64, b::Float64, lat::Float64)
    (a > 0.0 && b > 0.0) || _dsk_throw(:nonpositiveaxis, "the semi-axis lengths must be positive")
    slat, clat = sincos(lat)
    ept = _dsk_ednmpt(a, b, b, Vec3{Float64}(clat, slat, 0.0))
    d1 = b / a
    d2 = a / b
    return (1.0 - d1 * d1) * ept[1], (1.0 - d2 * d2) * ept[2]
end

"`zzellbds`: the semi-axes of the spheroids that bound the altitude range `hmin:hmax` over the spheroid `a >= b`."
@inline function _dsk_ellbds(a::Float64, b::Float64, hmax::Float64, hmin::Float64)
    b > 0.0 || _dsk_throw(:nonpositiveradius, "the polar radius must be positive")
    b > a && _dsk_throw(:radiioutoforder, "the radii are out of order")
    b + hmin > 0.0 || _dsk_throw(:lowerboundtoolow, "the lower altitude bound is too low")
    (hmin < 0.0 && b + a / b * hmin <= 0.0) && _dsk_throw(
        :lowerboundtoolow,
        "the lower altitude bound is too low"
    )
    hmin > hmax && _dsk_throw(:boundsoutoforder, "the altitude bounds are out of order")
    amax, bmax = hmax >= 0.0 ? (a + hmax, b + hmax * (a / b)) : (a + hmax * (b / a), b + hmax)
    amin, bmin = hmin <= 0.0 ? (a + hmin, b + hmin * (a / b)) : (a + hmin * (b / a), b + hmin)
    return amax, bmax, amin, bmin
end

"""
    _dsk_pdcmpl(re, f, p, lat) -> Int

`zzpdcmpl`: 1, 0 or -1 as the planetodetic latitude of `p` on the spheroid `re`, `f` is
greater than, equal to or less than `lat`.
"""
function _dsk_pdcmpl(re::Float64, f::Float64, p::Vec3{Float64}, lat::Float64)
    if p[1] == 0.0 && p[2] == 0.0
        if p[3] > 0.0
            return lat == _HALFPI ? 0 : 1
        elseif p[3] == 0.0
            return lat > 0.0 ? -1 : lat == 0.0 ? 0 : 1
        end
        return lat == -_HALFPI ? 0 : -1
    end
    rp = re * (1.0 - f)
    _, yincpt = _dsk_elnaxx(re, rp, lat)
    if lat == 0.0 || yincpt == 0.0
        return p[3] > 0.0 ? 1 : p[3] == 0.0 ? 0 : -1
    end
    offpcl = _dsk_reclat(Vec3{Float64}(p[1], p[2], p[3] - yincpt))[3]
    cmp = offpcl > lat ? 1 : offpcl == lat ? 0 : -1
    if lat > 0.0
        yincpt > 0.0 && return cmp
        return p[3] > 0.0 ? cmp : -1
    end
    yincpt < 0.0 && return cmp
    return p[3] < 0.0 ? cmp : 1
end

"`zzpdpltc`: true if `p` is on the correct side of the equator for the latitude cone of `lat`."
@inline function _dsk_pdpltc(re::Float64, f::Float64, p::Vec3{Float64}, lat::Float64)
    re > 0.0 || _dsk_throw(:valueoutofrange, "the equatorial radius must be positive")
    f < 1.0 || _dsk_throw(:valueoutofrange, "the flattening coefficient must be less than 1")
    f <= 0.0 && return true
    p[3] * lat >= 0.0 && return true
    abs(lat) >= 0.01 && return false
    xxpt, _ = _dsk_elnaxx(re, re * (1.0 - f), lat)
    r = sqrt(max(p[1] * p[1] + p[2] * p[2], 0.0))
    return r >= xxpt
end

# ---------------------------------------------------------------------------------------
# Volume elements: inclusion of a point
# ---------------------------------------------------------------------------------------

"""
    _dsk_nrmlon(inmin, inmax, tol) -> (outmin, outmax)

`zznrmlon`: longitude bounds with `outmin < outmax <= outmin + 2pi`, both in
`[-2pi, 2pi]`. Throws `SpiceError` (`:valueoutofrange`, `:zeroboundsextent`).
"""
function _dsk_nrmlon(inmin::Float64, inmax::Float64, tol::Float64)
    pi2 = _TWOPI
    tol < 0.0 && _dsk_throw(:valueoutofrange, "the longitude tolerance must be non-negative")
    (inmin < -pi2 - tol || inmin > pi2 + tol) && _dsk_throw(
        :valueoutofrange,
        "the longitude lower bound is outside [-2*pi, 2*pi]"
    )
    inmin == inmax && _dsk_throw(
        :zeroboundsextent,
        "the longitude lower bound is equal to the upper bound"
    )
    outmin = max(-pi2, min(inmin, pi2))
    (inmax < -pi2 - tol || inmax > pi2 + tol) && _dsk_throw(
        :valueoutofrange,
        "the longitude upper bound is outside [-2*pi, 2*pi]"
    )
    outmax = max(-pi2, min(inmax, pi2))
    if outmax <= outmin + tol
        if outmax <= 0.0
            outmax = min(outmax + pi2, pi2)
            outmax < outmin && (outmin = max(outmin - pi2, -pi2))
        else
            outmin = max(outmin - pi2, -pi2)
        end
    end
    delta = outmax - outmin
    delta > pi2 + tol && (outmax = outmax - pi2)
    outmin == outmax && _dsk_throw(
        :zeroboundsextent,
        "after adjustment the longitude bounds are equal"
    )
    return outmin, outmax
end

"""
    _dsk_lon_out0(lon, b) -> Bool

The longitude test of `zzinlat0` and `zzinpdt0`: true if `lon`, moved by a turn
toward the longitude bounds `b[1]`, `b[2]` when it is outside them, is still outside
them by more than 1e-12. Throws as `_dsk_nrmlon` for invalid bounds.
"""
@inline function _dsk_lon_out0(lon::Float64, b::NTuple{6, Float64})
    minlon, maxlon = _dsk_nrmlon(b[1], b[2], 1.0e-12)
    loclon = lon
    if lon < minlon - 1.0e-12
        loclon = lon + _TWOPI
    elseif lon > maxlon + 1.0e-12
        loclon = lon - _TWOPI
    end
    return loclon < minlon - 1.0e-12 || loclon > maxlon + 1.0e-12
end

"`zzinlat0`: inclusion of (r, lon, lat) in the latitudinal element `b`, without margin; `exclud` names a coordinate to skip."
function _dsk_inlat0(r::Float64, lon::Float64, lat::Float64, b::NTuple{6, Float64}, exclud::Int)
    (r == 0.0 && b[5] == 0.0) && return true
    exclud != 3 && (r < b[5] || r > b[6]) && return false
    exclud != 2 && (lat < b[3] || lat > b[4]) && return false
    exclud != 1 && _dsk_lon_out0(lon, b) && return false
    return true
end

"The longitude test with margin that `zzinlat` and `zzinpdt` share; `coslat` is the cosine of the planetocentric latitude."
@inline function _dsk_lon_in(lon::Float64, b::NTuple{6, Float64}, margin::Float64, coslat::Float64)
    minlon, maxlon = _dsk_nrmlon(b[1], b[2], 1.0e-12)
    lonmrg = max(1.0e-12, abs(margin))
    dlon = lonmrg / max(abs(coslat), 1.0e-8)
    amin = minlon - dlon
    amax = maxlon + dlon
    if lon < amin
        lon = lon < amin - 1.0e-12 ? lon + _TWOPI : amin
    elseif lon > amax
        lon = lon > amax + 1.0e-12 ? lon - _TWOPI : amax
    end
    return !(lon < amin || lon > amax)
end

"""
    _dsk_inlat(p, b, margin, exclud) -> Bool

`zzinlat`: true if `p` is inside the latitudinal volume element with the bounds `b`
(longitude, latitude, radius) expanded by `margin`. `exclud` in 0:3 names a
coordinate that is not tested. Requires `margin >= 0`.
"""
function _dsk_inlat(p::Vec3{Float64}, b::NTuple{6, Float64}, margin::Float64, exclud::Int)
    r, lon, lat = _dsk_reclat(p)
    margin == 0.0 && return _dsk_inlat0(r, lon, lat, b, exclud)
    margin < 0.0 && _dsk_throw(:valueoutofrange, "the margin must be non-negative")
    (0 <= exclud <= 3) || _dsk_throw(:valueoutofrange, "EXCLUD must be in the range 0:3")
    (r <= margin && b[5] <= margin) && return true
    if exclud != 3
        aminr = max(0.0, b[5] * (1.0 - margin))
        amaxr = b[6] * (margin + 1.0)
        (r < aminr || r > amaxr) && return false
    end
    if exclud != 2
        aminlt = max(-_HALFPI, b[3] - margin)
        amaxlt = min(_HALFPI, b[4] + margin)
        (lat < aminlt || lat > amaxlt) && return false
    end
    if exclud != 1
        if lat <= _HALFPI - 1.0e-8 && lat >= -_HALFPI + 1.0e-8
            _dsk_lon_in(lon, b, margin, cos(lat)) || return false
        else
            _dsk_nrmlon(b[1], b[2], 1.0e-12)       # the bounds are still validated
        end
    end
    return true
end

"`zzinrec`: true if `p` is inside the rectangular volume element `b` expanded by `margin` times each edge length."
function _dsk_inrec(p::Vec3{Float64}, b::NTuple{6, Float64}, margin::Float64, exclud::Int)
    margin < 0.0 && _dsk_throw(:valueoutofrange, "the margin must be non-negative")
    (0 <= exclud <= 3) || _dsk_throw(:indexoutofrange, "EXCLUD must be in the range 0:3")
    for i in 1:3
        b[2i] - b[2i - 1] < 0.0 && _dsk_throw(
            :boundsoutoforder,
            "the bounds of a rectangular volume element are out of order"
        )
    end
    for i in 1:3                         # invariant: p passes coordinates 1 .. i-1
        exclud == i && continue
        delta = margin * abs(b[2i] - b[2i - 1])
        (p[i] < b[2i - 1] - delta || p[i] > b[2i] + delta) && return false
    end
    return true
end

"`zzinpdt0`: inclusion of `p` (with planetocentric longitude `lon`) in the planetodetic element `b`, without margin."
function _dsk_inpdt0(
        p::Vec3{Float64}, lon::Float64, b::NTuple{6, Float64}, re::Float64,
        f::Float64, exclud::Int
    )
    rp = re * (1.0 - f)
    if exclud != 2
        minlat = max(b[3] - 1.0e-12, -_HALFPI)
        maxlat = min(b[4] + 1.0e-12, _HALFPI)
        minrel = _dsk_pdcmpl(re, f, p, minlat)
        maxrel = _dsk_pdcmpl(re, f, p, maxlat)
        (minrel == -1 || maxrel == 1) && return false
    end
    exclud != 1 && _dsk_lon_out0(lon, b) && return false
    if exclud != 3
        emax, pmax, emin, pmin = if f >= 0.0
            _dsk_ellbds(re, rp, b[6], b[5])
        else
            pmx, emx, pmn, emn = _dsk_ellbds(rp, re, b[6], b[5])
            (emx, pmx, emn, pmn)
        end
        d1 = p[1] / emax
        d2 = p[2] / emax
        d3 = p[3] / pmax
        d1 * d1 + d2 * d2 + d3 * d3 > 1.0 && return false
        d1 = p[1] / emin
        d2 = p[2] / emin
        d3 = p[3] / pmin
        d1 * d1 + d2 * d2 + d3 * d3 < 1.0 && return false
    end
    return true
end

"""
    _dsk_inpdt(p, b, re, f, margin, exclud) -> Bool

`zzinpdt`: true if `p` is inside the planetodetic volume element with the bounds `b`
(longitude, latitude, altitude) on the spheroid `re`, `f`, expanded by `margin`.
"""
function _dsk_inpdt(
        p::Vec3{Float64}, b::NTuple{6, Float64}, re::Float64, f::Float64,
        margin::Float64, exclud::Int
    )
    (0 <= exclud <= 3) || _dsk_throw(:valueoutofrange, "EXCLUD must be in the range 0:3")
    _, lon, pcnlat = _dsk_reclat(p)
    margin == 0.0 && return _dsk_inpdt0(p, lon, b, re, f, exclud)
    margin < 0.0 && _dsk_throw(:valueoutofrange, "the margin must be non-negative")
    if exclud != 2
        amnlat = max(-_HALFPI - 1.0e-12, b[3] - margin)
        amxlat = min(_HALFPI + 1.0e-12, b[4] + margin)
        relmin = _dsk_pdcmpl(re, f, p, amnlat)
        relmax = _dsk_pdcmpl(re, f, p, amxlat)
        (relmin == -1 || relmax == 1) && return false
    end
    if exclud != 3
        amnalt = b[5] - margin * abs(b[5])
        amxalt = b[6] + margin * abs(b[6])
        altbds = (-_PI, _PI, -_HALFPI, _HALFPI, amnalt, amxalt)
        _dsk_inpdt0(p, lon, altbds, re, f, 1) || return false
    end
    if exclud != 1
        minlon, maxlon = _dsk_nrmlon(b[1], b[2], 1.0e-12)
        relmax = _dsk_pdcmpl(re, f, p, _HALFPI - 1.0e-8)
        relmin = _dsk_pdcmpl(re, f, p, -_HALFPI + 1.0e-8)
        if relmax != 1 && relmin != -1
            _dsk_lon_in(lon, b, margin, cos(pcnlat)) || return false
        end
    end
    return true
end

"""
    _dsk_invelt(p, corsys, corpar, bounds, margin, exclud) -> Bool

`zzinvelt`: true if `p` is inside the volume element `bounds` of the coordinate system
`corsys` (1 latitudinal, 3 rectangular, 4 planetodetic), expanded by `margin`.
"""
@inline function _dsk_invelt(
        p::Vec3{Float64}, corsys::Int, corpar::NTuple{10, Float64},
        bounds::NTuple{6, Float64}, margin::Float64, exclud::Int
    )
    margin < 0.0 && _dsk_throw(:valueoutofrange, "the margin must be non-negative")
    corsys == 1 && return _dsk_inlat(p, bounds, margin, exclud)
    corsys == 4 && return _dsk_inpdt(p, bounds, corpar[1], corpar[2], margin, exclud)
    corsys == 3 && return _dsk_inrec(p, bounds, margin, exclud)
    _dsk_throw(:notsupported, "the coordinate system code $corsys is not recognized")
end

# ---------------------------------------------------------------------------------------
# Volume elements: intersection with a ray
# ---------------------------------------------------------------------------------------

"""
    _dsk_raybox(vertex, raydir, boxori, extent) -> (xpt, found)

`zzraybox`: the point where the ray enters the box with the corner `boxori` and the
edge lengths `extent`; the vertex itself if it is inside. Requires `raydir != 0` and
positive extents.
"""
function _dsk_raybox(
        vertex::Vec3{Float64}, raydir::Vec3{Float64}, boxori::Vec3{Float64},
        extent::Vec3{Float64}
    )
    _m_iszero(raydir) && _dsk_throw(:zerovector, "the ray direction is the zero vector")
    udir = _m_hat(raydir)
    min(min(extent[1], extent[2]), extent[3]) <= 0.0 && _dsk_throw(
        :valueoutofrange,
        "all box extents must be strictly positive"
    )
    center = Vec3{Float64}(
        boxori[1] + extent[1] / 2, boxori[2] + extent[2] / 2,
        boxori[3] + extent[3] / 2
    )
    offset = vertex - center
    if abs(offset[1]) <= extent[1] / 2 && abs(offset[2]) <= extent[2] / 2 &&
            abs(offset[3]) <= extent[3] / 2
        return vertex, true
    end
    r = _m_norm(extent) * 0.50049999999999994
    sphxpt = offset
    if !(_m_norm(offset) < r)
        sphxpt, sphfnd = _dsk_surfpt(offset, udir, r, r, r)
        sphfnd || return _DSK_ORIGIN, false
    end
    sphvtx = Vec3{Float64}(
        sphxpt[1] + center[1] - boxori[1], sphxpt[2] + center[2] - boxori[2],
        sphxpt[3] + center[3] - boxori[3]
    )
    maxidx = 1
    maxt = -1.0
    for i in 1:3                         # invariant: maxt is the largest plane distance of axes 1 .. i-1
        s = sphvtx[i]
        (s < 0.0 || s > extent[i]) || continue
        near = s < 0.0 ? 0.0 : extent[i]
        plndst = near - s
        limit = r * 2 * abs(udir[i])
        abs(plndst) > limit && return _DSK_ORIGIN, false
        t = plndst / udir[i]
        t < 0.0 && return _DSK_ORIGIN, false
        if t > maxt
            maxidx = i
            maxt = t
        end
    end
    maxt < 0.0 && return _DSK_ORIGIN, false
    x = 1.0 * sphvtx + maxt * udir
    x1 = x[1]
    x2 = x[2]
    x3 = x[3]
    for i in 1:3
        s = sphvtx[i]
        xi = x[i]
        if i == maxidx
            xi = s > extent[i] ? extent[i] : 0.0
        elseif !(s < 0.0 || s > extent[i])
            (xi < 0.0 || xi > extent[i]) && return _DSK_ORIGIN, false
        elseif s < 0.0
            xi = max(xi, 0.0)
            xi > extent[i] && return _DSK_ORIGIN, false
        else
            xi = min(xi, extent[i])
            xi < 0.0 && return _DSK_ORIGIN, false
        end
        i == 1 ? (x1 = xi) : i == 2 ? (x2 = xi) : (x3 = xi)
    end
    return Vec3{Float64}(x1, x2, x3) + boxori, true
end

"""
    _dsk_ryt_lon(inside, vertex, udir, minlon, maxlon, maxr, margin, nxpts, xpt, mndist)
        -> (nxpts, xpt, mndist)

The step of `zzrytlat` and `zzrytpdt` for the two longitude boundary half-planes: an
intercept on the west or the east plane that `inside(p)` accepts (the inclusion test of
the element without the longitude) and that is on the side of the element replaces
`xpt` if it is nearer than `mndist`. Nothing happens if the longitude range is a whole
turn.
"""
@inline function _dsk_ryt_lon(
        inside::F, vertex::Vec3{Float64}, udir::Vec3{Float64}, minlon::Float64,
        maxlon::Float64, maxr::Float64, margin::Float64, nxpts::Int, xpt::Vec3{Float64},
        mndist::Float64
    ) where {F}
    if cos(maxlon - minlon) < 1.0
        smin, cmin = sincos(minlon)
        westb = Vec3{Float64}(smin, -cmin, 0.0)
        s = (_m_norm(vertex) + maxr) * 1.1
        nx, srfx = _dsk_zzinrypl(vertex, udir, westb, 0.0, s)
        if nx == 1 && inside(srfx)
            wback = _m_ucrss(westb, _DSK_ZAXIS)
            if _m_dot(srfx, wback) < margin * maxr
                dist = _m_norm(vertex - srfx)
                if dist < mndist
                    xpt = srfx
                    nxpts = 1
                    mndist = dist
                end
            end
        end
        smax, cmax = sincos(maxlon)
        eastb = Vec3{Float64}(-smax, cmax, 0.0)
        nx, srfx = _dsk_zzinrypl(vertex, udir, eastb, 0.0, s)
        if nx == 1 && inside(srfx)
            eback = _m_ucrss(_DSK_ZAXIS, eastb)
            if _m_dot(srfx, eback) < margin * maxr
                dist = _m_norm(vertex - srfx)
                if dist < mndist
                    xpt = srfx
                    nxpts = 1
                    mndist = dist
                end
            end
        end
    end
    return nxpts, xpt, mndist
end

"""
    _dsk_rytlat(vertex, raydir, b, margin) -> (nxpts, xpt)

`zzrytlat`: the nearest point of the ray on the surface of the latitudinal volume
element `b` expanded by `margin`; the vertex itself if it is inside.
"""
function _dsk_rytlat(
        vertex::Vec3{Float64}, raydir::Vec3{Float64}, b::NTuple{6, Float64},
        margin::Float64
    )
    _dsk_inlat(vertex, b, margin * 2, 0) && return 1, vertex
    minr = b[5]
    maxr = b[6]
    minlon, maxlon = _dsk_nrmlon(b[1], b[2], 1.0e-12)
    minlat = b[3]
    maxlat = b[4]
    nxpts = 0
    xpt = _DSK_ORIGIN
    udir = _m_hat(raydir)
    mndist = _DSK_DPMAX
    srfx, found = _dsk_ryxsph(vertex, udir, maxr)
    found || return 0, xpt
    vr, _, vlat = _dsk_reclat(vertex)
    if _dsk_inlat(srfx, b, margin, 3)
        xpt = srfx
        nxpts = 1
        vr > maxr && return nxpts, xpt
        mndist = _m_norm(vertex - xpt)
    end
    s = _m_norm(vertex) + maxr * 1.1
    endpt2 = 1.0 * vertex + s * udir
    if maxlat < _HALFPI
        angle = max(0.0, _HALFPI - maxlat)
        nx, srfx, xpt2 = _dsk_incnsg(_DSK_ORIGIN, _DSK_ZAXIS, angle, vertex, endpt2)
        if nx >= 1
            if _dsk_inlat(srfx, b, margin, 2)
                dist = _m_norm(vertex - srfx)
                if dist < mndist
                    xpt = srfx
                    nxpts = 1
                    vlat > maxlat && return nxpts, xpt
                    mndist = dist
                end
            end
            if nx == 2 && _dsk_inlat(xpt2, b, margin, 2)
                dist = _m_norm(vertex - xpt2)
                if dist < mndist
                    xpt = xpt2
                    nxpts = 1
                    mndist = dist
                end
            end
        end
    end
    if minlat > -_HALFPI
        angle = _HALFPI - minlat
        nx, srfx, xpt2 = _dsk_incnsg(_DSK_ORIGIN, _DSK_ZAXIS, angle, vertex, endpt2)
        if nx >= 1
            if _dsk_inlat(srfx, b, margin, 2)
                dist = _m_norm(vertex - srfx)
                if dist < mndist
                    xpt = srfx
                    nxpts = 1
                    vlat < minlat && return nxpts, xpt
                    mndist = dist
                end
            end
            if nx == 2 && _dsk_inlat(xpt2, b, margin, 2)
                dist = _m_norm(vertex - xpt2)
                if dist < mndist
                    return 1, xpt2
                end
            end
        end
    end
    nxpts, xpt, mndist = _dsk_ryt_lon(
        p -> _dsk_inlat(p, b, margin, 1), vertex, udir, minlon, maxlon, maxr, margin,
        nxpts, xpt, mndist
    )
    if minr > 0.0
        srfx, found = _dsk_ryxsph(vertex, udir, minr)
        if found && _dsk_inlat(srfx, b, margin, 3)
            dist = _m_norm(vertex - srfx)
            if dist < mndist
                xpt = srfx
                nxpts = 1
                mndist = dist
            end
        end
        srfx, found = _dsk_ryxsph(endpt2, -udir, minr)
        if found && _dsk_inlat(srfx, b, margin, 3)
            dist = _m_norm(vertex - srfx)
            if dist < mndist
                xpt = srfx
                nxpts = 1
                mndist = dist
            end
        end
    end
    return nxpts, xpt
end

"`zzrytrec`: the point where the ray enters the rectangular volume element `b` expanded by `margin`."
function _dsk_rytrec(
        vertex::Vec3{Float64}, raydir::Vec3{Float64}, b::NTuple{6, Float64},
        margin::Float64
    )
    l = Vec3{Float64}(b[2] - b[1], b[4] - b[3], b[6] - b[5])
    (l[1] > 0.0 && l[2] > 0.0 && l[3] > 0.0) || _dsk_throw(
        :badcoordbounds,
        "the bounds of a rectangular volume element are out of order"
    )
    _dsk_inrec(vertex, b, margin * 2, 0) && return 1, vertex
    d1 = margin * abs(l[1])
    d2 = margin * abs(l[2])
    d3 = margin * abs(l[3])
    boxori = Vec3{Float64}(b[1] - d1, b[3] - d2, b[5] - d3)
    extent = Vec3{Float64}(l[1] + d1 * 2, l[2] + d2 * 2, l[3] + d3 * 2)
    xpt, found = _dsk_raybox(vertex, raydir, boxori, extent)
    return (found ? 1 : 0), xpt
end

"The latitude cone step that `zzrytpdt` does for the upper (`upper`) and the lower latitude bound."
@inline function _dsk_pdt_cone(
        vertex::Vec3{Float64}, endpt2::Vec3{Float64}, b::NTuple{6, Float64},
        re::Float64, f::Float64, rp::Float64, margin::Float64, lat::Float64, angle::Float64,
        upper::Bool, nxpts::Int, xpt::Vec3{Float64}, mndist::Float64
    )
    _, yincpt = _dsk_elnaxx(re, rp, lat)
    apex = Vec3{Float64}(0.0, 0.0, yincpt)
    vtxang = _m_vsep(vertex - apex, _DSK_ZAXIS)
    nx, srfx, xpt2 = _dsk_incnsg(apex, _DSK_ZAXIS, angle, vertex, endpt2)
    if nx > 0
        xval1 = _dsk_pdpltc(re, f, srfx, lat)
        xval2 = nx == 2 ? _dsk_pdpltc(re, f, xpt2, lat) : false
        if xval1 && !xval2
            nx = 1
        elseif xval2 && !xval1
            nx = 1
            srfx = xpt2
        elseif !xval1 && !xval2
            nx = 0
        end
    end
    done = false
    if nx >= 1
        if _dsk_inpdt(srfx, b, re, f, margin, 2)
            dist = _m_norm(vertex - srfx)
            if dist < mndist
                xpt = srfx
                nxpts = 1
                if upper
                    (vtxang < angle && (lat < 0.0 || vertex[3] > 0.0)) && return nxpts, xpt, mndist, true
                else
                    (vtxang > angle && (lat > 0.0 || vertex[3] < 0.0)) && return nxpts, xpt, mndist, true
                end
                mndist = dist
            end
        end
        if nx == 2 && _dsk_inpdt(xpt2, b, re, f, margin, 2)
            dist = _m_norm(vertex - xpt2)
            if dist < mndist
                xpt = xpt2
                nxpts = 1
                mndist = dist
                done = !upper
            end
        end
    end
    return nxpts, xpt, mndist, done
end

"""
    _dsk_rytpdt(vertex, raydir, b, re, f, margin) -> (nxpts, xpt)

`zzrytpdt`: the nearest point of the ray on the surface of the planetodetic volume
element `b` (spheroid `re`, `f`) expanded by `margin`; the vertex itself if it is
inside. Requires `margin >= 0` and `raydir != 0`.
"""
function _dsk_rytpdt(
        vertex::Vec3{Float64}, raydir::Vec3{Float64}, b::NTuple{6, Float64},
        re::Float64, f::Float64, margin::Float64
    )
    margin < 0.0 && _dsk_throw(:valueoutofrange, "the margin must be non-negative")
    _m_iszero(raydir) && _dsk_throw(:zerovector, "the ray direction is the zero vector")
    _dsk_inpdt(vertex, b, re, f, margin, 0) && return 1, vertex
    rp = re * (1.0 - f)
    minlon, maxlon = _dsk_nrmlon(b[1], b[2], 1.0e-12)
    minlat = b[3]
    maxlat = b[4]
    amnalt = b[5] - margin * abs(b[5])
    amxalt = b[6] + margin * abs(b[6])
    emax, pmax, emin, pmin = if re >= rp
        _dsk_ellbds(re, rp, amxalt, amnalt)
    else
        pmx, emx, pmn, emn = _dsk_ellbds(rp, re, amxalt, amnalt)
        (emx, pmx, emn, pmn)
    end
    nxpts = 0
    xpt = _DSK_ORIGIN
    udir = _m_hat(raydir)
    mndist = _DSK_DPMAX
    srfx, found = _dsk_surfpt(vertex, udir, emax, emax, pmax)
    found || return 0, xpt
    if _dsk_inpdt(srfx, b, re, f, margin, 3)
        xpt = srfx
        nxpts = 1
        d1 = vertex[1] / emax
        d2 = vertex[2] / emax
        d3 = vertex[3] / pmax
        d1 * d1 + d2 * d2 + d3 * d3 > 1.0 && return nxpts, xpt
        mndist = _m_norm(vertex - xpt)
    end
    maxr = max(emax, pmax)
    s = _m_norm(vertex) + maxr * 1.1
    endpt2 = 1.0 * vertex + s * udir
    if maxlat < _HALFPI
        angle = max(0.0, _HALFPI - maxlat)
        nxpts, xpt, mndist, done = _dsk_pdt_cone(
            vertex, endpt2, b, re, f, rp, margin, maxlat,
            angle, true, nxpts, xpt, mndist
        )
        done && return nxpts, xpt
    end
    if minlat > -_HALFPI
        angle = _HALFPI - minlat
        nxpts, xpt, mndist, done = _dsk_pdt_cone(
            vertex, endpt2, b, re, f, rp, margin, minlat,
            angle, false, nxpts, xpt, mndist
        )
        done && return nxpts, xpt
    end
    nxpts, xpt, mndist = _dsk_ryt_lon(
        p -> _dsk_inpdt(p, b, re, f, margin, 1), vertex, udir, minlon, maxlon, maxr,
        margin, nxpts, xpt, mndist
    )
    srfx, found = _dsk_surfpt(vertex, udir, emin, emin, pmin)
    if found && _dsk_inpdt(srfx, b, re, f, margin, 3)
        dist = _m_norm(vertex - srfx)
        if dist < mndist
            xpt = srfx
            nxpts = 1
            mndist = dist
        end
    end
    srfx, found = _dsk_surfpt(endpt2, -udir, emin, emin, pmin)
    if found && _dsk_inpdt(srfx, b, re, f, margin, 3)
        dist = _m_norm(vertex - srfx)
        if dist < mndist
            xpt = srfx
            nxpts = 1
        end
    end
    return nxpts, xpt
end

@inline _dsk_bounds(d::DSKDescr) = (d.co1min, d.co1max, d.co2min, d.co2max, d.co3min, d.co3max)

"""
`zzrytelt`: the nearest point of the ray on the volume element of the segment descriptor
`d`, with `margin`. The system of `d` is 1, 3 or 4: `_dsk_body_segments` has thrown the
error of `zzsegbox` for any other, which comes before that of `zzrytelt`.
"""
@inline function _dsk_rytelt(
        vertex::Vec3{Float64}, raydir::Vec3{Float64}, d::DSKDescr,
        margin::Float64
    )
    b = _dsk_bounds(d)
    d.corsys == 1 && return _dsk_rytlat(vertex, raydir, b, margin)
    d.corsys == 3 && return _dsk_rytrec(vertex, raydir, b, margin)
    return _dsk_rytpdt(vertex, raydir, b, d.corpar[1], d.corpar[2], margin)
end

# ---------------------------------------------------------------------------------------
# Bounding spheres of volume elements
# ---------------------------------------------------------------------------------------

"The box center and radius that `zzlatbox` and `zzpdtbox` build from the cylindrical extents."
@inline function _dsk_cylbox(
        outrad::Float64, inrad::Float64, minlon::Float64, maxlon::Float64,
        maxz::Float64, minz::Float64
    )
    hdlon = (maxlon - minlon) / 2
    midlon = minlon + hdlon
    shd, chd = sincos(hdlon)
    lr = hdlon <= _HALFPI ? outrad - inrad * chd : (1.0 - chd) * outrad
    lt = hdlon > _HALFPI ? outrad * 2.0 : outrad * 2.0 * shd
    lz = maxz - minz
    lr = max(0.0, lr)
    lt = max(0.0, lt)
    lz = max(0.0, lz)
    midr = outrad - lr / 2
    center = _dsk_cylrec(midr, midlon, minz + lz / 2)
    return center, _m_norm(Vec3{Float64}(lr / 2, lt / 2, lz / 2))
end

"Longitude and latitude bounds of `b` as `zzlatbox` and `zzpdtbox` check and clamp them."
function _dsk_box_lonlat(b::NTuple{6, Float64})
    minlon = b[1]
    maxlon = b[2]
    maxlon <= minlon && (maxlon += _TWOPI)
    maxlon <= minlon && _dsk_throw(:badlongituderange, "the longitude bounds are out of order")
    b[3] > b[4] && _dsk_throw(:badlatitudebounds, "the latitude bounds are out of order")
    (b[3] < -_HALFPI - 1.0e-12 || b[4] > _HALFPI + 1.0e-12) && _dsk_throw(
        :badlatituderange, "a latitude bound is outside [-pi/2, pi/2]"
    )
    return minlon, maxlon, max(b[3], -_HALFPI), min(b[4], _HALFPI)
end

"`zzlatbox`: center and radius of a sphere that contains the latitudinal volume element `b`."
function _dsk_latbox(b::NTuple{6, Float64})
    minlon, maxlon, minlat, maxlat = _dsk_box_lonlat(b)
    minr = b[5]
    maxr = b[6]
    smin, cmin = sincos(minlat)
    smax, cmax = sincos(maxlat)
    if minlat >= 0.0
        outrad = cmin * maxr
        inrad = cmax * minr
        maxz = maxr * smax
        minz = minr * smin
    elseif maxlat <= 0.0
        outrad = cmax * maxr
        inrad = cmin * minr
        maxz = minr * smax
        minz = maxr * smin
    else
        outrad = maxr
        inrad = cos(max(abs(maxlat), abs(minlat))) * minr
        maxz = maxr * smax
        minz = maxr * smin
    end
    return _dsk_cylbox(outrad, inrad, minlon, maxlon, maxz, minz)
end

"`zzrecbox`: center and radius of the sphere around the rectangular volume element `b`."
function _dsk_recbox(b::NTuple{6, Float64})
    l1 = b[2] - b[1]
    l2 = b[4] - b[3]
    l3 = b[6] - b[5]
    (l1 > 0.0 && l2 > 0.0 && l3 > 0.0) || _dsk_throw(
        :boundsoutoforder,
        "the bounds of a rectangular volume element are out of order"
    )
    center = Vec3{Float64}(b[1] + l1 / 2, b[3] + l2 / 2, b[5] + l3 / 2)
    return center, _m_norm(Vec3{Float64}(l1 / 2, l2 / 2, l3 / 2))
end

"`zzpdtbox`: center and radius of a sphere that contains the planetodetic volume element `b` on the spheroid `re`, `f`."
function _dsk_pdtbox(b::NTuple{6, Float64}, re::Float64, f::Float64)
    re > 0.0 || _dsk_throw(:valueoutofrange, "the equatorial radius must be positive")
    f < 1.0 || _dsk_throw(:valueoutofrange, "the flattening coefficient must be less than 1")
    minlon, maxlon, minlat, maxlat = _dsk_box_lonlat(b)
    minalt = b[5]
    maxalt = b[6]
    xy(v) = _m_norm(Vec3{Float64}(v[1], v[2], 0.0))
    if minlat >= 0.0
        outrad = xy(_dsk_georec(0.0, minlat, maxalt, re, f))
        inrad = xy(_dsk_georec(0.0, maxlat, minalt, re, f))
    elseif maxlat <= 0.0
        outrad = xy(_dsk_georec(0.0, maxlat, maxalt, re, f))
        inrad = xy(_dsk_georec(0.0, minlat, minalt, re, f))
    else
        outrad = re + maxalt
        inrad = xy(_dsk_georec(0.0, max(abs(maxlat), abs(minlat)), minalt, re, f))
    end
    if minlat >= 0.0
        topv = _dsk_georec(0.0, maxlat, maxalt, re, f)
        botv = _dsk_georec(0.0, minlat, minalt, re, f)
    elseif maxlat < 0.0
        topv = _dsk_georec(0.0, maxlat, minalt, re, f)
        botv = _dsk_georec(0.0, minlat, maxalt, re, f)
    else
        topv = _dsk_georec(0.0, maxlat, maxalt, re, f)
        botv = _dsk_georec(0.0, minlat, maxalt, re, f)
    end
    return _dsk_cylbox(outrad, inrad, minlon, maxlon, topv[3], botv[3])
end

"`zzsegbox`: center and radius of the bounding sphere of the volume element of the descriptor `d`."
function _dsk_segbox(d::DSKDescr)
    b = _dsk_bounds(d)
    d.corsys == 1 && return _dsk_latbox(b)
    d.corsys == 3 && return _dsk_recbox(b)
    d.corsys == 4 && return _dsk_pdtbox(b, d.corpar[1], d.corpar[2])
    _dsk_throw(:notsupported, "the coordinate system $(d.corsys) is not supported")
end

# ---------------------------------------------------------------------------------------
# Plates
# ---------------------------------------------------------------------------------------

"`pltexp` on three vertices: each vertex moves away from the centroid by the fraction `delta`."
@inline function _dsk_pltexp(v1::Vec3{Float64}, v2::Vec3{Float64}, v3::Vec3{Float64}, delta::Float64)
    d = delta / 3.0
    c1 = d * (v1[1] + v2[1] + v3[1])
    c2 = d * (v1[2] + v2[2] + v3[2])
    c3 = d * (v1[3] + v2[3] + v3[3])
    s = delta + 1.0
    return Vec3{Float64}(s * v1[1] - c1, s * v1[2] - c2, s * v1[3] - c3),
        Vec3{Float64}(s * v2[1] - c1, s * v2[2] - c2, s * v2[3] - c3),
        Vec3{Float64}(s * v3[1] - c1, s * v3[2] - c2, s * v3[3] - c3)
end

@inline _dsk_pltnrm(v1::Vec3{Float64}, v2::Vec3{Float64}, v3::Vec3{Float64}) =
    _m_cross(v2 - v1, v3 - v2)

"""
    _dsk_insang(v, e1, e2, e3) -> (found, scale)

`insang`: true if `v` is inside the solid angle of the edge vectors `e1`, `e2`, `e3`;
`scale * v` is then the point of the ray along `v` in the plane of the three points.
"""
@inline function _dsk_insang(v::Vec3{Float64}, e1::Vec3{Float64}, e2::Vec3{Float64}, e3::Vec3{Float64})
    norm12 = _m_cross(e1, e2)
    vn12 = _m_dot(v, norm12)
    en = _m_dot(e3, norm12)
    en == 0.0 && return false, 0.0
    ((en > 0.0 && vn12 < 0.0) || (en < 0.0 && vn12 > 0.0)) && return false, 0.0
    vn23 = _m_dot(v, _m_cross(e2, e3))
    ((en > 0.0 && vn23 < 0.0) || (en < 0.0 && vn23 > 0.0)) && return false, 0.0
    vn31 = _m_dot(v, _m_cross(e3, e1))
    ((en > 0.0 && vn31 < 0.0) || (en < 0.0 && vn31 > 0.0)) && return false, 0.0
    denom = vn12 + vn23 + vn31
    denom == 0.0 && return false, 0.0
    return true, en / denom
end

"`npsgpt`: the point of the segment from `ep1` to `ep2` that is nearest to `point`, and its distance."
function _dsk_npsgpt(ep1::Vec3{Float64}, ep2::Vec3{Float64}, point::Vec3{Float64})
    seg = ep2 - ep1
    _m_iszero(seg) && return ep1, _m_norm(ep1 - point)
    trans = point - ep1
    lnear = _m_proj(trans, seg, false) + ep1
    dist = _m_norm(lnear - point)
    offdot = _m_dot(lnear - ep1, seg)
    offdot < 0.0 && return ep1, _m_norm(ep1 - point)
    offdot > _m_dot(seg, seg) && return ep2, _m_norm(ep2 - point)
    return lnear, dist
end

"`pltnp` on `Vec3` values: the point of the plate nearest to `point`, and its distance."
function _dsk_pltnp(point::Vec3{Float64}, v1::Vec3{Float64}, v2::Vec3{Float64}, v3::Vec3{Float64})
    e1 = v2 - v1
    e2 = v3 - v2
    e3 = v1 - v3
    normal = _m_cross(e1, e2)
    enorm1 = _m_cross(e1, normal)
    enorm2 = _m_cross(e2, normal)
    enorm3 = _m_cross(e3, normal)
    if _m_iszero(normal) || _m_iszero(enorm1) || _m_iszero(enorm2) || _m_iszero(enorm3)
        l1 = _m_norm(e1)
        l2 = _m_norm(e2)
        l3 = _m_norm(e3)
        (l1 == 0.0 && l2 == 0.0) && return v1, _m_norm(v1 - point)
        l1 > max(l2, l3) && return _dsk_npsgpt(v1, v2, point)
        l2 > max(l3, l1) && return _dsk_npsgpt(v2, v3, point)
        return _dsk_npsgpt(v3, v1, point)
    end
    perp = _m_perp(point - v1, normal, false)
    in1 = _m_dot(perp, enorm1) <= 0.0
    in2 = _m_dot(perp, enorm2) <= _m_dot(e1, enorm2)
    in3 = _m_dot(perp, enorm3) <= 0.0
    if in1 && in2 && in3
        pnear = v1 + perp
        return pnear, _m_norm(pnear - point)
    end
    (!in1 && in2 && in3) && return _dsk_npsgpt(v1, v2, point)
    (!in2 && in3 && in1) && return _dsk_npsgpt(v2, v3, point)
    (!in3 && in1 && in2) && return _dsk_npsgpt(v3, v1, point)
    np1, d1 = _dsk_npsgpt(v1, v2, point)
    np2, d2 = _dsk_npsgpt(v2, v3, point)
    np3, d3 = _dsk_npsgpt(v3, v1, point)
    d1 <= min(d2, d3) && return np1, d1
    d2 <= min(d3, d1) && return np2, d2
    return np3, d3
end

# ---------------------------------------------------------------------------------------
# Type 2 segments in the DAS arrays
# ---------------------------------------------------------------------------------------

"`dskgd` on the double array: the DSK descriptor of the segment `d`; throws `SpiceError(:invalidformat)` if the segment has fewer than 24 doubles."
function _dsk_descr(dps::Vector{Float64}, d::DLADescr, path::AbstractString)
    d.dsize >= 24 || throw(
        SpiceError(
            :invalidformat, "the size of the double precision " *
                "component of the DLA segment is $(d.dsize); a DSK segment has at least 24 " *
                "doubles ($path)"
        )
    )
    (d.dbase >= 0 && d.dbase + d.dsize <= length(dps)) || throw(
        SpiceError(
            :dasnosuchaddress,
            "the doubles $(d.dbase + 1):$(d.dbase + d.dsize) of a DLA segment are not in $path"
        )
    )
    v = view(dps, (d.dbase + 1):(d.dbase + 24))
    return DSKDescr(
        round(Int, v[1]), round(Int, v[2]), round(Int, v[3]), round(Int, v[4]),
        round(Int, v[5]), round(Int, v[6]), ntuple(k -> v[6 + k], 10), v[17], v[18], v[19],
        v[20], v[21], v[22], v[23], v[24]
    )
end

"""
    _dsk_header(ints, dps, d, path) -> DSK02Header

The fixed-size items of the type 2 segment `d`, as `dskb02` reads them. Throws
`SpiceError(:dasnosuchaddress)` if the fixed-size items are not inside the arrays.
"""
function _dsk_header(ints::Vector{Int32}, dps::Vector{Float64}, d::DLADescr, path::AbstractString)
    (d.ibase >= 0 && d.ibase + 10 <= length(ints) && d.dbase >= 0 && d.dbase + 34 <= length(dps)) ||
        throw(
        SpiceError(
            :dasnosuchaddress,
            "the fixed-size items of a type 2 DSK segment are not inside $path"
        )
    )
    i = d.ibase
    b = d.dbase
    return DSK02Header(
        ints[i + 1], ints[i + 2], ints[i + 3],
        (dps[b + 25], dps[b + 26], dps[b + 27], dps[b + 28], dps[b + 29], dps[b + 30]), dps[b + 34],
        Vec3{Float64}(dps[b + 31], dps[b + 32], dps[b + 33]),
        (Int(ints[i + 4]), Int(ints[i + 5]), Int(ints[i + 6])), ints[i + 7], ints[i + 10], ints[i + 8],
        ints[i + 9]
    )
end

const _DSK_NO_HEADER = DSK02Header(
    0, 0, 0, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0.0,
    Vec3{Float64}(0.0, 0.0, 0.0), (0, 0, 0), 0, 0, 0, 0
)

"Number of coarse voxels of a type 2 segment; requires `hdr.cgscal >= 1`."
@inline _dsk_ncgr(hdr::DSK02Header) = hdr.nvxtot ÷ (hdr.cgscal * (hdr.cgscal * hdr.cgscal))

"""
    _dsk_segment(ints, dps, d, file, path, scan = false) -> DSKSegment

The descriptor of the segment `d` of a DSK file with the arrays `ints`, `dps`. For a
type 2 segment the sizes of the header are checked against the sizes of the DLA
segment, so that only the pointers inside the data can be out of range later; throws
`SpiceError(:invalidformat)` if they do not fit. The bounding sphere has a NaN radius
if the coordinate system is not 1, 3 or 4 or its bounds are not valid; a query on
such a segment then throws the error that CSPICE throws. With `scan`, the voxel
pointers of a type 2 segment are read once for its `voxtop` (`_dsk_voxtop`); else its
elements are -1.
"""
function _dsk_segment(
        ints::Vector{Int32}, dps::Vector{Float64}, d::DLADescr, file::Integer,
        path::AbstractString, scan::Bool = false
    )
    dsc = _dsk_descr(dps, d, path)
    (
        d.ibase >= 0 && d.isize >= 0 && d.ibase + d.isize <= length(ints) &&
            d.dbase + d.dsize <= length(dps)
    ) || throw(
        SpiceError(
            :dasnosuchaddress,
            "a DLA segment is not inside the address spaces of $path"
        )
    )
    ctr, rad = try
        _dsk_segbox(dsc)
    catch e
        e isa SpiceError || rethrow()
        (_DSK_ORIGIN, NaN)
    end
    hdr = _DSK_NO_HEADER
    if dsc.dtype == 2
        hdr = _dsk_header(ints, dps, d, path)
        ok = hdr.nv >= 0 && hdr.np >= 0 && hdr.voxnpt >= 0 && hdr.voxnpl >= 0 &&
            hdr.vtxnpl >= 0 && hdr.nvxtot >= 0 && hdr.cgscal >= 1
        ok = ok && 34 + 3 * hdr.nv <= d.dsize &&
            10 + 3 * hdr.np + hdr.voxnpt + hdr.voxnpl + hdr.nv + hdr.vtxnpl +
            _dsk_ncgr(hdr) <= d.isize
        ok || throw(
            SpiceError(
                :invalidformat, "the item counts of a type 2 DSK segment " *
                    "(NV = $(hdr.nv), NP = $(hdr.np), NVXTOT = $(hdr.nvxtot), CGSCAL = $(hdr.cgscal), " *
                    "VOXNPT = $(hdr.voxnpt), VOXNPL = $(hdr.voxnpl), VTXNPL = $(hdr.vtxnpl)) do not " *
                    "fit the sizes of its DLA segment in $path"
            )
        )
    end
    seg = DSKSegment(Int32(file), d, dsc, ctr, rad, hdr, _DSK_NOTOP)
    (scan && dsc.dtype == 2) || return seg
    return DSKSegment(Int32(file), d, dsc, ctr, rad, hdr, _dsk_voxtop(ints, seg))
end

"""
    _dsk_voxtop(ints, seg) -> NTuple{17,Int}

Element `i`: the sum of the `min(2^(i-1), k)` largest plate counts of the voxels of
the type 2 segment `seg`, where `k` is the sum of its grid extents (more voxels than a
ray crosses); the sum of all counts if fewer voxels have a positive count. A voxel is
an entry of the voxel pointer array that points into the voxel-plate list, so voxels
that share a list count once each. See `_dsk_top_bound`.
"""
function _dsk_voxtop(ints::Vector{Int32}, seg::DSKSegment)
    offptr = _dsk_off_voxptr(seg)
    offlst = _dsk_off_voxlst(seg)
    voxnpl = seg.hdr.voxnpl
    voxnpt = seg.hdr.voxnpt
    # A damaged header gives no bound; `dskx02` reports it when the segment is used.
    (
        all(>=(1), seg.hdr.vgrext) && voxnpt >= 0 && 0 <= offptr && 0 <= offlst &&
            offptr + voxnpt <= length(ints) && offlst + voxnpl <= length(ints)
    ) || return _DSK_NOTOP
    k = min(sum(Int, seg.hdr.vgrext), voxnpt)     # no more counts than pointers
    heap = Int[]                         # a min-heap of the largest counts so far
    sizehint!(heap, k)
    for j in 1:voxnpt
        # invariant: heap holds the min(k, m) largest positive counts of the pointers
        # 1 .. j-1, of which m are positive
        p = Int(@inbounds ints[offptr + j])
        1 <= p <= voxnpl || continue
        c = Int(@inbounds ints[offlst + p])
        c > 0 || continue
        if length(heap) < k
            push!(heap, c)
            _dsk_sift_up!(heap, length(heap))
        elseif c > heap[1]
            heap[1] = c
            _dsk_sift_down!(heap, 1)
        end
    end
    top = cumsum(sort!(heap; rev = true))
    @assert length(top) <= k && issorted(top)
    isempty(top) && return ntuple(_ -> 0, 17)
    return ntuple(i -> top[min(1 << (i - 1), length(top))], 17)
end

"""
    _dsk_top_bound(top, left) -> Int

An upper bound of the plate counts of any `left` voxels of a segment with the
`voxtop` `top`. For `2^(i-1) = a <= left`, the counts of the ranks after `a` are at
most the smallest count of the ranks `a/2 + 1` to `a`, so at most their mean.
"""
function _dsk_top_bound(top::NTuple{17, Int}, left::Int)
    @assert left >= 0 && top[1] >= 0
    left == 0 && return 0
    i = min(8 * sizeof(Int) - leading_zeros(left), 17)     # 2^(i-1) <= left
    i == 1 && return top[1]
    a = 1 << (i - 1)
    per = (top[i] - top[i - 1]) ÷ (a ÷ 2)
    b = top[i] + (left - a) * per
    @assert b >= top[i]
    return b
end

function _dsk_sift_up!(h::Vector{Int}, i::Int)
    while i > 1 && h[i ÷ 2] > h[i]       # invariant: h is a min-heap but for the node i
        h[i ÷ 2], h[i] = h[i], h[i ÷ 2]
        i ÷= 2
    end
    return h
end

function _dsk_sift_down!(h::Vector{Int}, i::Int)
    n = length(h)
    while true                           # invariant: h is a min-heap but for the node i
        m = i
        2i <= n && h[2i] < h[m] && (m = 2i)
        2i + 1 <= n && h[2i + 1] < h[m] && (m = 2i + 1)
        m == i && return h
        h[i], h[m] = h[m], h[i]
        i = m
    end
    return
end

"""
    open_dsk(path) -> DSKFile

Open the DSK file at `path` and return its integer and double arrays and its segment
descriptors in file order.

Contract: `path` names a readable DAS file with the ID word `DAS/DSK`, of either byte
order, segregated or not. In the result the arrays are memory-mapped where the file
allows it, and every segment lies inside them. Throws `SpiceError` with the kinds of
`dasopr`, `:invalidfiletype` for a DAS file that is not a DSK, and `:invalidformat` or
`:dasnosuchaddress` for a segment that does not fit the file.
"""
function open_dsk(path::AbstractString)
    h = _das_load(path, "furnsh", false)
    segs = try
        rstrip(h.idword) == "DAS/DSK" || throw(
            SpiceError(
                :invalidfiletype,
                "the DAS file $path has the ID word $(repr(h.idword)), not DAS/DSK"
            )
        )
        DSKSegment[_dsk_segment(h.ints, h.dps, d, 0, path, true) for d in _dla_segments(h.ints, path)]
    catch
        _das_release!(h)                 # see `_unmap!`
        rethrow()
    end
    return DSKFile(String(path), :DSK, h.ints, h.dps, segs)
end

# Offsets of the variable-size integer items of a type 2 segment (the address before
# the first element), in the order of the file.
@inline _dsk_off_plates(s::DSKSegment) = s.dladsc.ibase + 10
@inline _dsk_off_voxptr(s::DSKSegment) = _dsk_off_plates(s) + 3 * s.hdr.np
@inline _dsk_off_voxlst(s::DSKSegment) = _dsk_off_voxptr(s) + s.hdr.voxnpt
@inline _dsk_off_vtxptr(s::DSKSegment) = _dsk_off_voxlst(s) + s.hdr.voxnpl
@inline _dsk_off_vtxlst(s::DSKSegment) = _dsk_off_vtxptr(s) + s.hdr.nv
@inline _dsk_off_cgrptr(s::DSKSegment) = _dsk_off_vtxlst(s) + s.hdr.vtxnpl
@inline _dsk_off_verts(s::DSKSegment) = s.dladsc.dbase + 34

"Vertex `vid` of the segment; requires `1 <= vid <= nv` (else `SpiceError(:indexoutofrange)`)."
@inline function _dsk_vertex(dps::Vector{Float64}, s::DSKSegment, vid::Integer)
    (1 <= vid <= s.hdr.nv) || _dsk_throw(
        :indexoutofrange,
        "the vertex number $vid is not in the range 1:$(s.hdr.nv)"
    )
    k = _dsk_off_verts(s) + 3 * (Int(vid) - 1)
    @inbounds return Vec3{Float64}(dps[k + 1], dps[k + 2], dps[k + 3])
end

"The three vertex numbers of plate `plid`; requires `1 <= plid <= np` (else `SpiceError(:indexoutofrange)`)."
@inline function _dsk_plate(ints::Vector{Int32}, s::DSKSegment, plid::Integer)
    (1 <= plid <= s.hdr.np) || _dsk_throw(
        :indexoutofrange,
        "the plate number $plid is not in the range 1:$(s.hdr.np)"
    )
    k = _dsk_off_plates(s) + 3 * (Int(plid) - 1)
    @inbounds return Int(ints[k + 1]), Int(ints[k + 2]), Int(ints[k + 3])
end

@noinline _dsk_bad_pointer(what::String, value::Integer, size::Integer) = throw(
    SpiceError(
        :indexoutofrange, "the $what $value of a type 2 DSK segment is not in the range 1:$size"
    )
)

# ---------------------------------------------------------------------------------------
# Ray-surface intercept in one type 2 segment
# ---------------------------------------------------------------------------------------

"""
    _dsk_x02(ints, dps, seg, tols, vertex, raydir) -> (plid, xpt, found)

`dskx02`: the plate and the point where the ray from `vertex` along `raydir` (both in
the frame of the segment) first meets the plates of the type 2 segment `seg`.

Contract: `ints`, `dps` are the arrays of the file of `seg`; `raydir != 0`
(`SpiceError(:rayiszerovector)`). `found` is false if the ray meets no plate inside
the volume element of the segment. The function allocates nothing and writes nothing.

The voxel grid is walked along the ray as `xdda` does it, and the plates of each
voxel are tested in file order; a hit counts if it lies in the voxel that is being
tested (within the tolerance `xfract`) and in the volume element of the segment
(within `sgreed`), and the nearest hit of the first voxel that has one is returned.
CSPICE first collects the voxels and the plates in buffers, splits the voxels in two
groups and removes the plates that occur twice. The result here is the same, because
a second test of a plate has the same outcome as the first. The limits of two buffers
are kept, and both throw `SpiceError(:arraytoosmall)`, whether the ray meets a plate or
not: `xdda` lists every voxel of the ray before a plate is tested, so a ray with 50000
voxels or more in the grid throws (`_dsk_xdda_rest`); and the plate lists of a group of
voxels hold 256000 plates at most (`_dsk_plate_buffer`).
"""
function _dsk_x02(
        ints::Vector{Int32}, dps::Vector{Float64}, seg::DSKSegment,
        tols::DSKTolerances, vertex::Vec3{Float64}, raydir::Vec3{Float64}
    )
    nohit = (0, _DSK_ORIGIN, false)
    _m_iszero(raydir) && _dsk_throw(:rayiszerovector, "the ray direction is the zero vector")
    udir = _m_hat(raydir)
    hdr = seg.hdr
    dsc = seg.dskdsc
    voxsiz = hdr.voxsiz
    voxsiz == 0.0 && _dsk_throw(:invalidvalue, "the voxel size of a DSK segment is zero")
    cgscal = hdr.cgscal
    cgscal < 1 && _dsk_throw(:valueoutofrange, "the coarse voxel scale of a DSK segment is less than 1")
    vg1, vg2, vg3 = hdr.vgrext
    grdext = Vec3{Float64}(vg1 * voxsiz, vg2 * voxsiz, vg3 * voxsiz)
    grdtol = max(max(grdext[1], grdext[2]), grdext[3]) * 1.0e-12
    cg1 = vg1 ÷ cgscal
    cg2 = vg2 ÷ cgscal
    cgscl2 = cgscal * cgscal
    ncgr = _dsk_ncgr(hdr)
    ncgr > _DSK_MAXCGR && _dsk_throw(:gridtoolarge, "the coarse voxel grid of a DSK segment is too large")
    xpdfrc = tols.xfract
    greedm = tols.sgreed
    xtol = xpdfrc * max(max(abs(grdext[1]), abs(grdext[2])), abs(grdext[3]))
    voxori = hdr.voxori
    vtx2, boxhit = _dsk_raybox(vertex, raydir, voxori, grdext)
    boxhit || return nohit
    coord = Vec3{Float64}(
        (vtx2[1] - voxori[1]) / voxsiz, (vtx2[2] - voxori[2]) / voxsiz,
        (vtx2[3] - voxori[3]) / voxsiz
    )

    # --- the set-up of `xdda` ------------------------------------------------------
    min(min(vg1, vg2), vg3) < 1 && _dsk_throw(
        :baddimensions,
        "the voxel grid dimensions of a DSK segment must be strictly positive"
    )
    vg = (vg1, vg2, vg3)
    for i in 1:3
        (coord[i] < vg[i] * -1.0e-12 || coord[i] > vg[i] * 1.0000000000010001) &&
            _dsk_throw(:vertexnotingrid, "the ray vertex is outside of the voxel grid")
    end
    ax1 = 1
    maxcmp = abs(udir[1])
    for i in 2:3
        if abs(udir[i]) > maxcmp
            ax1 = i
            maxcmp = abs(udir[i])
        end
    end
    ax2 = _DAS_NEXT[ax1]
    ax3 = _DAS_NEXT[ax2]
    ic1 = clamp(unsafe_trunc(Int, coord[ax1]) + 1, 1, vg[ax1])
    ic2 = clamp(unsafe_trunc(Int, coord[ax2]) + 1, 1, vg[ax2])
    ic3 = clamp(unsafe_trunc(Int, coord[ax3]) + 1, 1, vg[ax3])
    vo1 = _m_brcktd(coord[ax1] - (ic1 - 1), 0.0, 1.0)
    vo2 = _m_brcktd(coord[ax2] - (ic2 - 1), 0.0, 1.0)
    vo3 = _m_brcktd(coord[ax3] - (ic3 - 1), 0.0, 1.0)
    limit = 1.0e-20 / vg[ax1] * abs(udir[ax1])
    ax2err = _DSK_DPMAX
    ax3err = _DSK_DPMAX
    s12 = 0.0
    s13 = 0.0
    if abs(udir[ax2]) > limit
        s12 = abs(udir[ax1] / udir[ax2])
        if udir[ax1] > 0.0
            ax2err = udir[ax2] > 0.0 ? s12 * (1.0 - vo2) + vo1 - 1.0 : s12 * vo2 + vo1 - 1.0
        else
            ax2err = udir[ax2] > 0.0 ? s12 * (1.0 - vo2) - vo1 : s12 * vo2 - vo1
        end
    end
    if abs(udir[ax3]) > limit
        s13 = abs(udir[ax1] / udir[ax3])
        if udir[ax1] > 0.0
            ax3err = udir[ax3] > 0.0 ? s13 * (1.0 - vo3) + vo1 - 1.0 : s13 * vo3 + vo1 - 1.0
        else
            ax3err = udir[ax3] > 0.0 ? s13 * (1.0 - vo3) - vo1 : s13 * vo3 - vo1
        end
    end
    step1 = udir[ax1] > limit ? 1 : udir[ax1] < -limit ? -1 : 0
    step2 = udir[ax2] > limit ? 1 : udir[ax2] < -limit ? -1 : 0
    step3 = udir[ax3] > limit ? 1 : udir[ax3] < -limit ? -1 : 0

    # --- the vertex from which the plates are seen -----------------------------------
    vtxoff = vertex - voxori
    if vtxoff[1] < -grdtol || vtxoff[1] > grdtol + grdext[1] || vtxoff[2] < -grdtol ||
            vtxoff[2] > grdtol + grdext[2] || vtxoff[3] < -grdtol || vtxoff[3] > grdtol + grdext[3]
        vtx2 = 1.0 * voxori + voxsiz * coord + (-1.0) * udir
    end

    offplt = _dsk_off_plates(seg)
    offptr = _dsk_off_voxptr(seg)
    offlst = _dsk_off_voxlst(seg)
    offcgr = _dsk_off_cgrptr(seg)
    corsys = dsc.corsys
    bounds = _dsk_bounds(dsc)

    start = (ic1, ic2, ic3, ax2err, ax3err)
    have = false
    near = 0.0
    plid = 0
    xpt = _DSK_ORIGIN
    nvx = 0
    nc = 0                           # voxels of the walk that are in a non-empty coarse voxel
    nref = 0                         # entries of their plate lists
    while true
        # invariant: (ic1, ic2, ic3) is voxel number nvx + 1 of the ray in the axis order
        # (ax1, ax2, ax3), and it is inside the grid; no voxel before it has a valid hit
        nvx += 1
        vx1 = ax1 == 1 ? ic1 : ax2 == 1 ? ic2 : ic3
        vx2 = ax1 == 2 ? ic1 : ax2 == 2 ? ic2 : ic3
        vx3 = ax1 == 3 ? ic1 : ax2 == 3 ? ic2 : ic3
        g1 = (vx1 - 1) ÷ cgscal + 1
        g2 = (vx2 - 1) ÷ cgscal + 1
        g3 = (vx3 - 1) ÷ cgscal + 1
        cvid = g1 + cg1 * (g2 + g3 * cg2 - (cg2 + 1))
        (1 <= cvid <= ncgr) || _dsk_bad_pointer("coarse voxel number", cvid, ncgr)
        cgrptr = Int(@inbounds ints[offcgr + cvid])
        if cgrptr > 0
            nc += 1
            fx = vx1 - cgscal * (g1 - 1)
            fy = vx2 - cgscal * (g2 - 1)
            fz = vx3 - cgscal * (g3 - 1)
            j = cgrptr + (fx + cgscal * (fy - 1) + cgscl2 * (fz - 1)) - 1
            (1 <= j <= hdr.voxnpt) || _dsk_bad_pointer("voxel pointer index", j, hdr.voxnpt)
            voxptr = Int(@inbounds ints[offptr + j])
            if voxptr != -1
                (1 <= voxptr <= hdr.voxnpl) || _dsk_bad_pointer("voxel-plate pointer", voxptr, hdr.voxnpl)
                nplate = Int(@inbounds ints[offlst + voxptr])
                (nplate <= 0 || voxptr + nplate <= hdr.voxnpl) ||
                    _dsk_bad_pointer("voxel-plate list end", voxptr + nplate, hdr.voxnpl)
                nref += max(nplate, 0)
                for k in 1:nplate
                    # invariant: (have, near, plid, xpt) is the nearest valid hit among
                    # the plates 1 .. k-1 of this voxel
                    pj = Int(@inbounds ints[offlst + voxptr + k])
                    (1 <= pj <= hdr.np) || _dsk_bad_pointer("plate number", pj, hdr.np)
                    q = offplt + 3 * (pj - 1)
                    p1 = _dsk_vertex(dps, seg, @inbounds ints[q + 1])
                    p2 = _dsk_vertex(dps, seg, @inbounds ints[q + 2])
                    p3 = _dsk_vertex(dps, seg, @inbounds ints[q + 3])
                    e1, e2, e3 = _dsk_pltexp(
                        p1 - vtx2, p2 - vtx2,
                        p3 - vtx2, xpdfrc
                    )
                    hits, scale = _dsk_insang(udir, e1, e2, e3)
                    hits || continue
                    _m_dot(udir, _dsk_pltnrm(p1, p2, p3)) <= 0.0 || continue
                    (!have || scale < near) || continue
                    xpt2 = 1.0 * vtx2 + scale * udir
                    h1 = (xpt2[1] - voxori[1]) / voxsiz
                    h2 = (xpt2[2] - voxori[2]) / voxsiz
                    h3 = (xpt2[3] - voxori[3]) / voxsiz
                    invox = h1 > vx1 - xtol - 1 && h1 < vx1 + xtol && h2 > vx2 - xtol - 1 &&
                        h2 < vx2 + xtol && h3 > vx3 - xtol - 1 && h3 < vx3 + xtol
                    invox || continue
                    _dsk_invelt(xpt2, corsys, dsc.corpar, bounds, greedm, 0) || continue
                    xpt = xpt2
                    have = true
                    near = scale
                    plid = pj
                end
                if have
                    _dsk_xdda_rest(
                        nvx, ic1, ic2, ic3, ax2err, ax3err, s12, s13, step1, step2, step3,
                        vg[ax1], vg[ax2], vg[ax3]
                    )
                    _dsk_plate_buffer(
                        ints, seg, (ax1, ax2, ax3), (s12, s13), (step1, step2, step3), start,
                        nc, nref, (ic1, ic2, ic3, ax2err, ax3err)
                    )
                    return plid, xpt, true
                end
            end
        end
        # --- the step of `xdda` to the next voxel, counted before the grid test ------
        ic1, ic2, ic3, ax2err, ax3err = _dsk_xdda_step(
            ic1, ic2, ic3, ax2err, ax3err, s12, s13, step1, step2, step3
        )
        nvx + 1 > _DSK_MAXVXL && _dsk_xdda_full()
        (1 <= ic1 <= vg[ax1] && 1 <= ic2 <= vg[ax2] && 1 <= ic3 <= vg[ax3]) || break
    end
    _dsk_plate_buffer(
        ints, seg, (ax1, ax2, ax3), (s12, s13), (step1, step2, step3), start, nc, nref, nothing
    )
    return nohit
end

"The room of the plate buffer of `dskx02` (`PLATID`, `SOURCE`, `ORDVEC`)."
const _DSK_PLTBUF = 256000

"The fixed quantities of the voxel walk of `dskx02` on one ray: the axis order, the
steps and slopes of `xdda`, the grid extents in that order, the coarse grid, the
offsets of the arrays of the segment and its `voxtop`."
struct _DSKWalk
    ax1::Int
    ax2::Int
    ax3::Int
    s12::Float64
    s13::Float64
    step1::Int
    step2::Int
    step3::Int
    n1::Int
    n2::Int
    n3::Int
    cgscal::Int
    cg1::Int
    cg2::Int
    ncgr::Int
    offcgr::Int
    offptr::Int
    offlst::Int
    voxnpt::Int
    voxnpl::Int
    voxtop::NTuple{17, Int}
end

"""
    _dsk_walk_plates(ints, w, state, first, maxnc) -> (nc, nref)

Walk the voxels of the ray from the voxel `state` (indices and errors of `xdda` in the
axis order of `w`), that voxel itself if `first`, to the end of the grid or until `nc`
reaches `maxnc`: `nc` voxels in a non-empty coarse voxel, whose plate lists have
`nref` entries in all. A pointer out of range counts as an empty list; the walk does
not throw.
"""
function _dsk_walk_plates(ints::Vector{Int32}, w::_DSKWalk, state, first::Bool, maxnc::Int)
    ic1, ic2, ic3, ax2err, ax3err = state
    nc = 0
    nref = 0
    if !first
        ic1, ic2, ic3, ax2err, ax3err = _dsk_xdda_step(
            ic1, ic2, ic3, ax2err, ax3err, w.s12, w.s13, w.step1, w.step2, w.step3
        )
    end
    cgscl2 = w.cgscal * w.cgscal
    while nc < maxnc && 1 <= ic1 <= w.n1 && 1 <= ic2 <= w.n2 && 1 <= ic3 <= w.n3
        # invariant: (nc, nref) count the voxels of the walk before (ic1, ic2, ic3)
        vx1 = w.ax1 == 1 ? ic1 : w.ax2 == 1 ? ic2 : ic3
        vx2 = w.ax1 == 2 ? ic1 : w.ax2 == 2 ? ic2 : ic3
        vx3 = w.ax1 == 3 ? ic1 : w.ax2 == 3 ? ic2 : ic3
        g1 = (vx1 - 1) ÷ w.cgscal + 1
        g2 = (vx2 - 1) ÷ w.cgscal + 1
        g3 = (vx3 - 1) ÷ w.cgscal + 1
        cvid = g1 + w.cg1 * (g2 + g3 * w.cg2 - (w.cg2 + 1))
        cgrptr = 1 <= cvid <= w.ncgr ? Int(@inbounds ints[w.offcgr + cvid]) : 0
        if cgrptr > 0
            nc += 1
            j = cgrptr + (
                vx1 - w.cgscal * (g1 - 1) + w.cgscal * (vx2 - w.cgscal * (g2 - 1) - 1) +
                    cgscl2 * (vx3 - w.cgscal * (g3 - 1) - 1)
            ) - 1
            voxptr = 1 <= j <= w.voxnpt ? Int(@inbounds ints[w.offptr + j]) : -1
            if 1 <= voxptr <= w.voxnpl
                nref += max(Int(@inbounds ints[w.offlst + voxptr]), 0)
            end
        end
        ic1, ic2, ic3, ax2err, ax3err = _dsk_xdda_step(
            ic1, ic2, ic3, ax2err, ax3err, w.s12, w.s13, w.step1, w.step2, w.step3
        )
    end
    return nc, nref
end

@noinline _dsk_plate_buffer_full(n::Int) = _dsk_throw(
    :arraytoosmall,
    "the plate lists of a group of voxels of the ray have $n entries; the plate " *
        "buffer of dskx02 has room for $_DSK_PLTBUF"
)

"""
    _dsk_plate_buffer(ints, seg, axes, slopes, steps, start, nc, nref, hit)

The plate buffer of `dskx02` for the ray whose walk in the segment `seg` (axis order,
slopes and steps of `xdda`) begins at the voxel `start`: `nc`
voxels of the walk, to the voxel `hit` of the first valid hit (`nothing`: to the end of
the grid), are in a non-empty coarse voxel and have `nref` plate list entries. CSPICE
splits the `n` such voxels of the whole ray into a first group of `max(1, (n + 1) ÷ 2)`
and a second group of the rest, collects the plate lists of the first group, then, if
the first group has no hit, those of the second, and throws
`SpiceError(:arraytoosmall)` for a group with more than `_DSK_PLTBUF` entries. After a
hit, the rest of the ray is walked unless the voxels left (at most the steps left in
the three axes) are bound by `voxtop` to keep the whole ray within `_DSK_PLTBUF`; the
first group is walked again only if the ray has more than `_DSK_PLTBUF` entries in all.
"""
@noinline function _dsk_plate_buffer(
        ints::Vector{Int32}, seg::DSKSegment, axes::NTuple{3, Int},
        slopes::NTuple{2, Float64}, steps::NTuple{3, Int}, start, nc::Int, nref::Int, hit
    )
    @assert nc >= 0 && nref >= 0
    hdr = seg.hdr
    vg = hdr.vgrext
    cgscal = hdr.cgscal
    w = _DSKWalk(
        axes..., slopes..., steps..., vg[axes[1]], vg[axes[2]], vg[axes[3]], cgscal,
        vg[1] ÷ cgscal, vg[2] ÷ cgscal, _dsk_ncgr(hdr), _dsk_off_cgrptr(seg),
        _dsk_off_voxptr(seg), _dsk_off_voxlst(seg), hdr.voxnpt, hdr.voxnpl, seg.voxtop
    )
    n, total = nc, nref
    if hit !== nothing && w.voxtop[1] >= 0
        left = _dsk_moves_left(hit[1], w.step1, w.n1) + _dsk_moves_left(hit[2], w.step2, w.n2) +
            _dsk_moves_left(hit[3], w.step3, w.n3)
        nref + _dsk_top_bound(w.voxtop, left) <= _DSK_PLTBUF && return nothing
    end
    if hit !== nothing
        rc, rref = _dsk_walk_plates(ints, w, hit, false, typemax(Int))
        n += rc
        total += rref
    end
    total <= _DSK_PLTBUF && return nothing
    ngrp1 = max(1, (n + 1) ÷ 2)
    _, ref1 = _dsk_walk_plates(ints, w, start, true, ngrp1)
    ref1 > _DSK_PLTBUF && _dsk_plate_buffer_full(ref1)
    # The second group is collected if it exists and the first group has no hit.
    (n > ngrp1 && (hit === nothing || nc > ngrp1) && total - ref1 > _DSK_PLTBUF) &&
        _dsk_plate_buffer_full(total - ref1)
    return nothing
end

"One step of `xdda`: the indices (axis order of the walk) and errors of the next voxel."
@inline function _dsk_xdda_step(
        ic1::Int, ic2::Int, ic3::Int, ax2err::Float64, ax3err::Float64, s12::Float64,
        s13::Float64, step1::Int, step2::Int, step3::Int
    )
    if ax2err < 0.0 || ax3err < 0.0
        if ax2err < ax3err
            ic2 += step2
            ax2err += s12
        else
            ic3 += step3
            ax3err += s13
        end
    else
        ic1 += step1
        step2 != 0 && (ax2err += -1.0)
        step3 != 0 && (ax3err += -1.0)
    end
    return ic1, ic2, ic3, ax2err, ax3err
end

@noinline _dsk_xdda_full() = _dsk_throw(
    :arraytoosmall,
    "the ray has $_DSK_MAXVXL voxels or more in the voxel grid of a DSK segment; " *
        "the voxel list of xdda has room for $_DSK_MAXVXL"
)

# Steps left in one axis before the walk leaves the grid of extent `n` there.
@inline _dsk_moves_left(ic::Int, step::Int, n::Int) = step > 0 ? n - ic : step < 0 ? ic - 1 : 0

"""
    _dsk_xdda_rest(nvx, ic1, ic2, ic3, ax2err, ax3err, s12, s13, step1, step2, step3, n1, n2, n3)

The rest of the walk of `xdda` after a hit in the voxel number `nvx` (`n1`, `n2`, `n3`
are the extents of the grid in the axis order of the walk): `SpiceError(:arraytoosmall)`
if the ray has `_DSK_MAXVXL` voxels or more in the grid. Each step moves one index by
one in a fixed direction, so the walk is needed only when the steps left in the three
axes could reach the limit.
"""
function _dsk_xdda_rest(
        nvx::Int, ic1::Int, ic2::Int, ic3::Int, ax2err::Float64, ax3err::Float64,
        s12::Float64, s13::Float64, step1::Int, step2::Int, step3::Int,
        n1::Int, n2::Int, n3::Int
    )
    nvx + _dsk_moves_left(ic1, step1, n1) + _dsk_moves_left(ic2, step2, n2) +
        _dsk_moves_left(ic3, step3, n3) < _DSK_MAXVXL && return nothing
    while true
        # invariant: (ic1, ic2, ic3) is voxel number nvx of the ray, inside the grid
        ic1, ic2, ic3, ax2err, ax3err = _dsk_xdda_step(
            ic1, ic2, ic3, ax2err, ax3err, s12, s13, step1, step2, step3
        )
        nvx + 1 > _DSK_MAXVXL && _dsk_xdda_full()
        (1 <= ic1 <= n1 && 1 <= ic2 <= n2 && 1 <= ic3 <= n3) || return nothing
        nvx += 1
    end
    return
end

# ---------------------------------------------------------------------------------------
# Selection of segments
# ---------------------------------------------------------------------------------------

"True if `id` is in `srflst`, or if `srflst` is empty (every surface is accepted)."
@inline function _dsk_surfok(srflst, id::Int)
    isempty(srflst) && return true
    for s in srflst
        s == id && return true
    end
    return false
end

@inline _dsk_eligible(seg::DSKSegment, bodyid::Int32, srflst, et::Float64) =
    seg.dskdsc.center == bodyid && _dsk_surfok(srflst, seg.dskdsc.surfce) &&
    et >= seg.dskdsc.start && et <= seg.dskdsc.stop

"""
    _dsk_offset(ks, dsc) -> Vec3

Position of the center of the frame of the segment relative to the center body of the
segment, in the frame of the segment at the middle of its coverage (the vector
`STOFF` of `zzdsksba`): zero when the frame is centered on the body. Throws
`SpiceError(:noframeinfo)` if the frame is unknown.
"""
function _dsk_offset(ks::KernelSet, dsc::DSKDescr)
    fd = _frame_def_or_nothing(ks, dsc.frmcde)
    fd === nothing && _dsk_no_frame_info(dsc.frmcde)
    fd.center == dsc.center && return _DSK_ORIGIN
    return _dsk_offset_spk(ks, dsc, fd.center)
end

@noinline function _dsk_offset_spk(ks::KernelSet, dsc::DSKDescr, frmctr::Int32)
    et = (dsc.start + dsc.stop) / 2
    return spk_jet(ks, frmctr, et, Int32(dsc.frmcde), Int32(dsc.center), Val(0))[1]
end

# `zzdsksba` signals SPICE(NOFRAMEINFO) and `zzdsksph` SPICE(NOFRAMEDATA) for a segment
# frame that `frinfo` does not know.
@noinline _dsk_no_frame_info(id::Integer) = _dsk_throw(
    :noframeinfo, "could not look up frame info for segment frame having ID $id"
)
@noinline _dsk_no_frame_data(id::Integer) = _dsk_throw(
    :noframedata, "no frame specification was found for frame ID $id"
)

# `zzdskbss`, which begins every search for DSK segments, signals SPICE(NOLOADEDDSKFILES)
# when no DSK file is loaded; with a DSK file loaded for another body the search finds
# no segment, and the caller gives its own error.
@noinline _dsk_no_files() = throw(
    SpiceError(
        :noloadeddskfiles, "at least one DSK file " *
            "must be loaded before a search for DSK segments can start"
    )
)
@inline _dsk_check_files(ks::KernelSet) = (isempty(ks.dskfiles) && _dsk_no_files(); nothing)

@noinline function _dsk_no_segments(bodyid::Integer, srflst)
    isempty(srflst) && throw(
        SpiceError(
            :dskdatanotfound,
            "no segments were found matching the body ID $bodyid"
        )
    )
    throw(
        SpiceError(
            :dskdatanotfound, "no segments were found matching the body ID " *
                "$bodyid and the surface list <$(join(srflst, ", "))>"
        )
    )
end

@noinline _dsk_no_body_segments(bodyid::Integer) = _dsk_throw(
    :nodsksegments, "the input segment " *
        "list was empty. This may be due to no DSKs containing data for body $bodyid " *
        "having been loaded."
)

"""
    _dsk_body_segments(ks, bodyid) -> Vector{DSKSegment}

The loaded segments of `bodyid` in the search order, with the checks that `zzdsksba`
makes on every one of them when it builds the segment table of the body, before a
search selects segments by surface and time: the error of `zzsegbox` and
`SpiceError(:noframeinfo)`. Throws `SpiceError(:noloadeddskfiles)` if no DSK file is
loaded and `:nodsksegments` (from `zzdskbux` and `zzdskbun`) if the body has no
segment. The SPK offset of a frame centered on another body is looked up only for the
segments that a search uses (see the page of CSPICE defects).
"""
function _dsk_body_segments(ks::KernelSet, bodyid::Int32)
    _dsk_check_files(ks)
    segs = get(ks.dsk, bodyid, nothing)
    (segs === nothing || isempty(segs)) && _dsk_no_body_segments(bodyid)
    for seg in segs
        dsc = seg.dskdsc
        isnan(seg.boxrad) && _dsk_segbox(dsc)        # throws the error of `zzsegbox`
        _frame_def_or_nothing(ks, dsc.frmcde) === nothing && _dsk_no_frame_info(dsc.frmcde)
    end
    return segs
end

"""
    _dsk_segments_ok(ks, segs) -> Bool

Whether the checks that `_dsk_body_segments` makes on the loaded segments `segs` of a
body pass: `false` if one of them would throw. Does not throw.
"""
function _dsk_segments_ok(ks::KernelSet, segs::Vector{DSKSegment})
    isempty(segs) && return false
    for seg in segs                      # invariant: the segments before `seg` pass
        isnan(seg.boxrad) && return false  # `_dsk_segment` set NaN where `zzsegbox` throws
        _frame_def_or_nothing(ks, seg.dskdsc.frmcde) === nothing && return false
    end
    return true
end

"""
    _dsk_sph(ks, segs, bodyid, srflst) -> (minrad, maxrad)

`zzdsksph`: radii of an inner and an outer bounding sphere, centered on the body, of
all segments of `bodyid` with a surface in `srflst` (any time). Throws
`SpiceError(:noframedata)` if the frame of such a segment is unknown, and
`:dskdatanotfound` if there is no such segment.
"""
function _dsk_sph(ks::KernelSet, segs::Vector{DSKSegment}, bodyid::Int32, srflst)
    maxr = -1.0
    minr = _DSK_DPMAX
    for seg in segs                      # invariant: minr, maxr bound the segments before `seg`
        dsc = seg.dskdsc
        _dsk_surfok(srflst, dsc.surfce) || continue
        fd = _frame_def_or_nothing(ks, dsc.frmcde)
        fd === nothing && _dsk_no_frame_data(dsc.frmcde)
        offmag = fd.center == bodyid ? 0.0 : _m_norm(_dsk_offset_spk(ks, dsc, fd.center))
        if dsc.corsys == 1
            sgminr = dsc.co3min
            sgmaxr = dsc.co3max
        elseif dsc.corsys == 4
            re = dsc.corpar[1]
            f = dsc.corpar[2]
            rp = re * (1.0 - f)
            sgmaxr, sgminr = f >= 0.0 ? (re + dsc.co3max, rp + dsc.co3min) :
                (rp + dsc.co3max, re + dsc.co3min)
        elseif dsc.corsys == 3
            boxctr, boxrad = _dsk_recbox(_dsk_bounds(dsc))
            sgminr = max(_m_norm(boxctr) - boxrad, 0.0)
            sgmaxr = _m_norm(boxctr) + boxrad
        else
            _dsk_throw(:notsupported, "the coordinate system $(dsc.corsys) is not supported")
        end
        ctrmnr = offmag <= sgminr ? sgminr - offmag : offmag >= sgmaxr ? offmag - sgmaxr : 0.0
        sgmaxr += offmag
        minr = min(minr, ctrmnr)
        maxr = max(maxr, sgmaxr)
    end
    maxr < 0.0 && _dsk_no_segments(bodyid, srflst)
    return minr, maxr
end

"`zzdsksph` of the body `bodyid` (`_dsk_sph` on its loaded segments), as `zzsudski` calls it."
function _dsk_body_sph(ks::KernelSet, bodyid::Int32, srflst)
    _dsk_check_files(ks)
    segs = get(ks.dsk, bodyid, nothing)
    segs === nothing && _dsk_no_segments(bodyid, srflst)
    return _dsk_sph(ks, segs, bodyid, srflst)
end

"The ray in the frame of the segment: rotation from the fixed frame, vertex and direction."
@inline function _dsk_to_segment(
        ks::KernelSet, dsc::DSKDescr, fixfid::Int32, et::Float64,
        vertex::Vec3{Float64}, raydir::Vec3{Float64}
    )
    xform = Mat3{Float64}(frame_rotation(ks, fixfid, Int32(dsc.frmcde), et))
    segvtx = _m_mxv(xform, vertex) - _dsk_offset(ks, dsc)
    return xform, segvtx, _m_mxv(xform, raydir)
end

const _DSK_MAXHIT = 1000             # the hit list of `zzdskbux`
const _DSK_HITBUF = 64               # hits kept on the stack; more go to `_dsk_bux_heap`

"""
    _dsk_collect!(hitseg, hitdst, ks, segs, tols, bodyid, srflst, et, fixfid, vertex,
                  sphvtx, raydir) -> (nhit, multfr)

The first loop of `zzdskbux`: the eligible segments whose volume element the ray meets,
in list order. Hit `k` is the segment `segs[hitseg[k]]`, whose volume element the ray
enters at the distance `hitdst[k]` from its vertex. The first `length(hitseg)` hits
are stored and `nhit` counts all of them; hit 1001 is `SpiceError(:buffertoosmall)`.
`multfr` tells whether an eligible segment has a frame other than `fixfid`. `sphvtx` is
the vertex moved onto the outer bounding sphere; it is used for the segments in the
frame `fixfid` that come before the first segment of another frame, as in CSPICE.
"""
@inline function _dsk_collect!(
        hitseg::AbstractVector{Int}, hitdst::AbstractVector{Float64}, ks::KernelSet,
        segs::Vector{DSKSegment}, tols::DSKTolerances, bodyid::Int32, srflst,
        et::Float64, fixfid::Int32, vertex::Vec3{Float64}, sphvtx::Vec3{Float64},
        raydir::Vec3{Float64}
    )
    sgmarg = tols.sgreed
    room = length(hitseg)
    nhit = 0
    multfr = false
    for i in eachindex(segs)
        # invariant: nhit counts the hit segments before i, and the first min(nhit, room)
        # of them are stored
        seg = @inbounds segs[i]
        _dsk_eligible(seg, bodyid, srflst, et) || continue
        dsc = seg.dskdsc
        segvtx = sphvtx
        segdir = raydir
        if dsc.frmcde != fixfid
            multfr = true
            _, segvtx, segdir = _dsk_to_segment(ks, dsc, fixfid, et, vertex, raydir)
        elseif multfr
            segvtx = vertex
        end
        _dsk_nplnpt_dist(segvtx, segdir, seg.boxctr) <= seg.boxrad || continue
        nxpts, xpt = _dsk_rytelt(segvtx, segdir, dsc, sgmarg)
        nxpts > 0 || continue
        nhit == _DSK_MAXHIT && _dsk_throw(
            :buffertoosmall,
            "too many segments were hit by the input ray; the buffer size is $_DSK_MAXHIT"
        )
        nhit += 1
        if nhit <= room
            @inbounds hitseg[nhit] = i
            @inbounds hitdst[nhit] = _m_norm(segvtx - xpt)
        end
    end
    return nhit, multfr
end

"""
    _dsk_orderd!(iorder, a, n) -> iorder

`orderd`: the order vector of `a[1:n]` in `iorder[1:n]`, by the Shell sort of SPICELIB.
The sort is not stable; equal values come in the order that `orderd` gives them, which
is the order in which `zzdskbux` searches segments with equal entry distances.
"""
@inline function _dsk_orderd!(iorder::AbstractVector{Int}, a::AbstractVector{Float64}, n::Int)
    for i in 1:n
        @inbounds iorder[i] = i
    end
    gap = n ÷ 2
    while gap > 0
        for i in (gap + 1):n
            j = i - gap
            while j > 0
                jg = j + gap
                oj = @inbounds iorder[j]
                ojg = @inbounds iorder[jg]
                if @inbounds(a[oj]) <= @inbounds(a[ojg])
                    j = 0
                else
                    @inbounds iorder[j] = ojg
                    @inbounds iorder[jg] = oj
                end
                j -= gap
            end
        end
        gap ÷= 2
    end
    return iorder
end

"""
    _dsk_search(ks, segs, tols, hitseg, hitdst, iorder, nhit, multfr, et, fixfid,
                vertex, sphvtx, raydir) -> (xpt, segidx, plid, found)

The second loop of `zzdskbux` and its end: the hit segments in the order `iorder`, until
the nearest intercept found so far is no farther from the vertex than the entry point
of the next segment. `xpt` is in the frame `fixfid`, and `segidx` indexes `segs`.
"""
@inline function _dsk_search(
        ks::KernelSet, segs::Vector{DSKSegment}, tols::DSKTolerances,
        hitseg::AbstractVector{Int}, hitdst::AbstractVector{Float64},
        iorder::AbstractVector{Int}, nhit::Int, multfr::Bool, et::Float64, fixfid::Int32,
        vertex::Vec3{Float64}, sphvtx::Vec3{Float64}, raydir::Vec3{Float64}
    )
    found = false
    dmin = 0.0
    xpt = _DSK_ORIGIN
    winner = 0
    plid = 0
    for k in 1:nhit
        # invariant: (found, dmin, xpt, winner, plid) is the nearest intercept in the
        # segments of the hits iorder[1:k-1]
        idx = @inbounds hitseg[iorder[k]]
        seg = @inbounds segs[idx]
        dsc = seg.dskdsc
        segvtx = multfr ? vertex : sphvtx
        segdir = raydir
        if dsc.frmcde != fixfid
            _, segvtx, segdir = _dsk_to_segment(ks, dsc, fixfid, et, vertex, raydir)
        end
        dsc.dtype == 2 || _dsk_throw(
            :typenotsupported,
            "the DSK data type $(dsc.dtype) is not supported"
        )
        file = ks.dskfiles[seg.file]
        p, srfx, xfnd = _dsk_x02(file.ints, file.dps, seg, tols, segvtx, segdir)
        if xfnd
            d = _m_norm(srfx - segvtx)
            if !found || d < dmin
                found = true
                dmin = d
                xpt = srfx
                winner = idx
                plid = p
            end
        end
        (found && k < nhit && dmin <= @inbounds(hitdst[iorder[k + 1]])) && break
    end
    found || return _DSK_ORIGIN, 0, 0, false
    dsc = segs[winner].dskdsc
    if dsc.frmcde != fixfid
        xform = Mat3{Float64}(frame_rotation(ks, fixfid, Int32(dsc.frmcde), et))
        off = _dsk_offset(ks, dsc)
        _m_iszero(off) || (xpt = xpt + off)
        xpt = _m_mtxv(xform, xpt)
    end
    return xpt, winner, plid, true
end

"""
    _dsk_bux_ray(ks, tols, segs, maxrad, bodyid, srflst, et, fixfid, vertex, raydir)
        -> (xpt, segidx, plid, found)

`zzdskbux` after its set-up, for one ray: `segs` is the list of `_dsk_body_segments`
and `maxrad` the outer radius of `_dsk_sph`, which a caller with many rays computes
once. Allocates nothing when the ray meets at most `_DSK_HITBUF` volume elements.
"""
function _dsk_bux_ray(
        ks::KernelSet, tols::DSKTolerances, segs::Vector{DSKSegment}, maxrad::Float64,
        bodyid::Int32, srflst, et::Float64, fixfid::Int32, vertex::Vec3{Float64},
        raydir::Vec3{Float64}
    )
    nohit = (_DSK_ORIGIN, 0, 0, false)
    maxrad *= 1.01
    sphvtx = vertex
    if _m_norm(vertex) > maxrad
        sphvtx, xfnd = _dsk_surfpt(vertex, raydir, maxrad, maxrad, maxrad)
        xfnd || return nohit
    end
    hitseg = MVector{_DSK_HITBUF, Int}(undef)
    hitdst = MVector{_DSK_HITBUF, Float64}(undef)
    nhit, multfr = _dsk_collect!(
        hitseg, hitdst, ks, segs, tols, bodyid, srflst, et, fixfid, vertex, sphvtx, raydir
    )
    nhit == 0 && return nohit
    nhit > _DSK_HITBUF && return _dsk_bux_heap(
        ks, tols, segs, nhit, bodyid, srflst, et, fixfid, vertex, sphvtx, raydir
    )
    iorder = MVector{_DSK_HITBUF, Int}(undef)
    _dsk_orderd!(iorder, hitdst, nhit)
    return _dsk_search(
        ks, segs, tols, hitseg, hitdst, iorder, nhit, multfr, et, fixfid, vertex,
        sphvtx, raydir
    )
end

# `_dsk_bux_ray` for a ray that meets more than `_DSK_HITBUF` volume elements: the hits
# are collected again into buffers of length `nhit`.
@noinline function _dsk_bux_heap(
        ks::KernelSet, tols::DSKTolerances, segs::Vector{DSKSegment}, nhit::Int,
        bodyid::Int32, srflst, et::Float64, fixfid::Int32, vertex::Vec3{Float64},
        sphvtx::Vec3{Float64}, raydir::Vec3{Float64}
    )
    hitseg = Vector{Int}(undef, nhit)
    hitdst = Vector{Float64}(undef, nhit)
    iorder = Vector{Int}(undef, nhit)
    _, multfr = _dsk_collect!(
        hitseg, hitdst, ks, segs, tols, bodyid, srflst, et, fixfid, vertex, sphvtx, raydir
    )
    _dsk_orderd!(iorder, hitdst, nhit)
    return _dsk_search(
        ks, segs, tols, hitseg, hitdst, iorder, nhit, multfr, et, fixfid, vertex,
        sphvtx, raydir
    )
end

"""
    _dsk_bux(ks, tols, bodyid, srflst, et, fixfid, vertex, raydir) -> (xpt, segidx, plid, found)

`zzdskbux` with `zzsbfxri`: the first point where the ray from `vertex` along `raydir`
(in the body-fixed frame `fixfid`, at the epoch `et`) meets the surface that the
loaded DSK segments of `bodyid` with a surface in `srflst` define (an empty `srflst`
accepts every surface). All matching segments have the same priority and the nearest
intercept wins. `segidx` indexes `ks.dsk[bodyid]`, and `plid` is the plate.

Throws the errors of `_dsk_body_segments`, then `SpiceError(:dskdatanotfound)` if no
segment has a surface of `srflst`. Allocates nothing (see `_dsk_bux_ray`).
"""
function _dsk_bux(
        ks::KernelSet, tols::DSKTolerances, bodyid::Int32, srflst, et::Float64,
        fixfid::Int32, vertex::Vec3{Float64}, raydir::Vec3{Float64}
    )
    segs = _dsk_body_segments(ks, bodyid)
    _, maxrad = _dsk_sph(ks, segs, bodyid, srflst)
    return _dsk_bux_ray(ks, tols, segs, maxrad, bodyid, srflst, et, fixfid, vertex, raydir)
end

# ---------------------------------------------------------------------------------------
# Plate that holds a surface point, and the outward normal there
# ---------------------------------------------------------------------------------------

"""
    _dsk_ptpl02(ints, dps, seg, tols, point) -> (plid, v1, v2, v3, found)

`zzptpl02`: the plate of the type 2 segment `seg` that is nearest to `point` (frame of
the segment) among the plates of the voxel of `point`, if `point` is on it within the
tolerance `ptmemm` times the radius of the bounding sphere.
"""
function _dsk_ptpl02(
        ints::Vector{Int32}, dps::Vector{Float64}, seg::DSKSegment,
        tols::DSKTolerances, point::Vec3{Float64}
    )
    none = (0, _DSK_ORIGIN, _DSK_ORIGIN, _DSK_ORIGIN, false)
    hdr = seg.hdr
    dsc = seg.dskdsc
    hdr.voxsiz == 0.0 && _dsk_throw(:valueoutofrange, "the voxel size of a DSK segment is zero")
    hdr.cgscal == 0 && _dsk_throw(:valueoutofrange, "the coarse voxel scale of a DSK segment is zero")
    isnan(seg.boxrad) && _dsk_segbox(dsc)
    ptsrfm = tols.ptmemm
    limit = ptsrfm * seg.boxrad
    xpdfrc = tols.xfract
    _dsk_invelt(point, dsc.corsys, dsc.corpar, _dsk_bounds(dsc), ptsrfm, 0) || return none
    offset = point - hdr.voxori
    v1 = clamp(unsafe_trunc(Int, offset[1] / hdr.voxsiz) + 1, 1, hdr.vgrext[1])
    v2 = clamp(unsafe_trunc(Int, offset[2] / hdr.voxsiz) + 1, 1, hdr.vgrext[2])
    v3 = clamp(unsafe_trunc(Int, offset[3] / hdr.voxsiz) + 1, 1, hdr.vgrext[3])
    cs = hdr.cgscal
    g1 = (v1 - 1) ÷ cs + 1
    g2 = (v2 - 1) ÷ cs + 1
    g3 = (v3 - 1) ÷ cs + 1
    ptroff = (v3 - cs * (g3 - 1) - 1) * cs * cs + (v2 - cs * (g2 - 1) - 1) * cs + (v1 - cs * (g1 - 1))
    cgrvid = g1 + (hdr.vgrext[1] ÷ cs) * (g2 - 1 + (g3 - 1) * (hdr.vgrext[2] ÷ cs))
    ncgr = _dsk_ncgr(hdr)
    (1 <= cgrvid <= ncgr) || _dsk_bad_pointer("coarse voxel number", cgrvid, ncgr)
    cgrptr = Int(ints[_dsk_off_cgrptr(seg) + cgrvid])
    cgrptr < 1 && return none
    ptrloc = cgrptr - 1 + ptroff
    (1 <= ptrloc <= hdr.voxnpt) || _dsk_bad_pointer("voxel pointer index", ptrloc, hdr.voxnpt)
    pltptr = Int(ints[_dsk_off_voxptr(seg) + ptrloc])
    pltptr < 1 && return none
    pltptr <= hdr.voxnpl || _dsk_bad_pointer("voxel-plate pointer", pltptr, hdr.voxnpl)
    offlst = _dsk_off_voxlst(seg)
    nplate = Int(ints[offlst + pltptr])
    nplate < 1 && return none
    pltptr + nplate <= hdr.voxnpl || _dsk_bad_pointer("voxel-plate list end", pltptr + nplate, hdr.voxnpl)
    dmin = _DSK_DPMAX
    best = none
    for j in 1:nplate                    # invariant: `best` is the nearest valid plate of 1 .. j-1
        pid = Int(ints[offlst + pltptr + j])
        i1, i2, i3 = _dsk_plate(ints, seg, pid)
        p1 = _dsk_vertex(dps, seg, i1)
        p2 = _dsk_vertex(dps, seg, i2)
        p3 = _dsk_vertex(dps, seg, i3)
        x1, x2, x3 = _dsk_pltexp(p1, p2, p3, xpdfrc)
        normal = _m_hat(_dsk_pltnrm(x1, x2, x3))
        dist = _DSK_DPMAX
        if abs(_m_dot(point - x1, normal)) <= limit
            _, dist = _dsk_pltnp(point, x1, x2, x3)
        end
        if dist <= limit && dist < dmin
            dmin = dist
            best = (pid, p1, p2, p3, true)
        end
    end
    return best
end

"""
    _dsk_bun(ks, tols, bodyid, srflst, et, fixfid, point) -> Vec3

`zzdskbun` with `zzsbfnrm`: the outward unit normal, in the frame `fixfid`, of the
plate that holds the surface point `point` (frame `fixfid`). The first segment in the
search order whose volume element holds the point and that has a plate there decides.
As in `zzdskbun`, a first pass finds all the segments whose volume element holds the
point, and a second pass looks for the plate in them.

Throws the errors of `_dsk_body_segments`, then `SpiceError(:toomanyhits)`,
`:pointnotinsegment`, `:notsupported` (a data type other than 2) or `:pointoffsurface`.
"""
function _dsk_bun(
        ks::KernelSet, tols::DSKTolerances, bodyid::Int32, srflst, et::Float64,
        fixfid::Int32, point::Vec3{Float64}
    )
    segs = _dsk_body_segments(ks, bodyid)
    nhit = 0
    for seg in segs                      # invariant: nhit counts the hit segments before `seg`
        _dsk_bun_inside(ks, tols, seg, bodyid, srflst, et, fixfid, point)[1] || continue
        nhit == _DSK_MAXHIT && _dsk_throw(:toomanyhits, "too many segments contain the input point")
        nhit += 1
    end
    nhit == 0 && _dsk_throw(
        :pointnotinsegment,
        "the input point is not contained in any DSK segment for body $bodyid"
    )
    for seg in segs                      # invariant: no hit segment before `seg` has a plate at the point
        inside, xform, segpt = _dsk_bun_inside(ks, tols, seg, bodyid, srflst, et, fixfid, point)
        inside || continue
        dsc = seg.dskdsc
        dsc.dtype == 2 || _dsk_throw(:notsupported, "the DSK data type $(dsc.dtype) is not supported")
        file = ks.dskfiles[seg.file]
        _, p1, p2, p3, xfnd = _dsk_ptpl02(file.ints, file.dps, seg, tols, segpt)
        xfnd || continue
        normal = _m_hat(_dsk_pltnrm(p1, p2, p3))
        return dsc.frmcde != fixfid ? _m_mtxv(xform, normal) : normal
    end
    _dsk_throw(
        :pointoffsurface, "the input point is not on the surface defined by the " *
            "DSK segments for body $bodyid"
    )
end

# The test of the first loop of `zzdskbun`: whether `seg` is eligible and its volume
# element holds the point, with the rotation to the frame of `seg` and the point there.
@inline function _dsk_bun_inside(
        ks::KernelSet, tols::DSKTolerances, seg::DSKSegment, bodyid::Int32, srflst,
        et::Float64, fixfid::Int32, point::Vec3{Float64}
    )
    xform = one(Mat3{Float64})
    _dsk_eligible(seg, bodyid, srflst, et) || return false, xform, point
    dsc = seg.dskdsc
    segpt = point
    if dsc.frmcde != fixfid
        xform = Mat3{Float64}(frame_rotation(ks, fixfid, Int32(dsc.frmcde), et))
        segpt = _m_mxv(xform, point) - _dsk_offset(ks, dsc)
    end
    _m_norm(seg.boxctr - segpt) <= seg.boxrad || return false, xform, segpt
    inside = _dsk_invelt(segpt, dsc.corsys, dsc.corpar, _dsk_bounds(dsc), tols.sgreed, 0)
    return inside, xform, segpt
end
