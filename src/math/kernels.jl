# Private vector and matrix kernels of the files in `src/math/` (prefix `_m_`): input
# conversion, 3-vector and 3 by 3 matrix arithmetic in the operation order of SPICELIB,
# rotations and the derivatives of unit vectors. They are included early, so that every
# file can use them; the public calls are in `math/vectors.jl` and the files after it.

# --- constants of SPICELIB ---------------------------------------------------------

"`pi_()`, `twopi_()` and `halfpi_()` of CSPICE as plain doubles."
const _PI = Float64(pi)
const _TWOPI = 2 * Float64(pi)
const _HALFPI = Float64(pi) / 2
const _RPD = pi / 180                # radians per degree, `rpd_`

"The ID of the frame J2000."
const _J2000_FRAME = Int32(1)

"Seconds per day (`spd_()`) and the Julian date of J2000 (`j2000_()`)."
const _SPD = 86400.0
const _JD2000 = 2451545.0

# --- input conversion ---------------------------------------------------------------

@inline _m_f(x::Integer) = Float64(x)
@inline _m_f(x::AbstractIrrational) = Float64(x)
@inline _m_f(x::Real) = x

@noinline _m_baddim(what::String, n) =
    throw(SpiceError(:badarg, "$what must have $n elements"))

@inline function _m_v3(a::AbstractVector)
    length(a) == 3 || _m_baddim("a vector", 3)
    @inbounds return SVector(_m_f(a[1]), _m_f(a[2]), _m_f(a[3]))
end
@inline _m_v3(a::SVector{3, <:AbstractFloat}) = a

@inline function _m_v6(a::AbstractVector)
    length(a) == 6 || _m_baddim("a state vector", 6)
    @inbounds return SVector(
        _m_f(a[1]), _m_f(a[2]), _m_f(a[3]), _m_f(a[4]), _m_f(a[5]),
        _m_f(a[6])
    )
end

@inline function _m_v4(a::AbstractVector)
    length(a) == 4 || _m_baddim("a quaternion", 4)
    @inbounds return SVector(_m_f(a[1]), _m_f(a[2]), _m_f(a[3]), _m_f(a[4]))
end

@inline function _m_m3(m::AbstractMatrix)
    size(m) == (3, 3) || _m_baddim("a matrix", "3 by 3")
    @inbounds return SMatrix{3, 3}(
        _m_f(m[1, 1]), _m_f(m[2, 1]), _m_f(m[3, 1]), _m_f(m[1, 2]),
        _m_f(m[2, 2]), _m_f(m[3, 2]), _m_f(m[1, 3]), _m_f(m[2, 3]), _m_f(m[3, 3])
    )
end
@inline _m_m3(m::SMatrix{3, 3, <:AbstractFloat, 9}) = m

@inline function _m_m6(m::AbstractMatrix)
    size(m) == (6, 6) || _m_baddim("a state transformation", "6 by 6")
    return SMatrix{6, 6}(ntuple(k -> @inbounds(_m_f(m[k])), Val(36)))
end
@inline _m_m6(m::SMatrix{6, 6, <:AbstractFloat, 36}) = m

@inline _m_pv(s::SVector{6}) = SVector(s[1], s[2], s[3]), SVector(s[4], s[5], s[6])

@inline _m_state(p::SVector{3}, v::SVector{3}) = SVector(p[1], p[2], p[3], v[1], v[2], v[3])

@inline _m_maxabs(v::SVector{3}) =
    max(abs(primal(v[1])), abs(primal(v[2])), abs(primal(v[3])))

@inline _m_iszero(v::SVector{3}) =
    primal(v[1]) == 0.0 && primal(v[2]) == 0.0 && primal(v[3]) == 0.0

# --- private kernels (arguments are static vectors and matrices) -----------------------

@inline _m_dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

@inline _m_cross(a, b) = SVector(
    a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
)

"`sqrt(a ⋅ a)`, without the scaling of `vnorm`, for every element type (no comparison)."
@inline _norm3(a) = sqrt(_m_dot(a, a))

"Norm with the scaling by the largest component of `vnorm`."
@inline function _m_norm(v::SVector{3, T}) where {T}
    vmax = _m_maxabs(v)
    vmax == 0.0 && return zero(T)
    a = v[1] / vmax
    b = v[2] / vmax
    c = v[3] / vmax
    return vmax * sqrt(a * a + b * b + c * c)
end

"`(unit vector, norm)` as `unorm`: the zero vector for a zero input."
@inline function _m_unorm(v::SVector{3})
    vmag = _m_norm(v)
    primal(vmag) > 0.0 || return zero(v), vmag
    return SVector(v[1] / vmag, v[2] / vmag, v[3] / vmag), vmag
end

@inline _m_hat(v::SVector{3}) = _m_unorm(v)[1]

"Projection of `a` on `b`: with `recip` the scaling of `vproj_c` (a product with the reciprocal of the largest component), without it the division of `vproj_`."
@inline function _m_proj(a::SVector{3}, b::SVector{3}, recip::Bool)
    biga = _m_maxabs(a)
    bigb = _m_maxabs(b)
    (biga == 0.0 || bigb == 0.0) && return zero(a) - zero(b)
    t = recip ? (1.0 / biga) * a : a / biga
    r = recip ? (1.0 / bigb) * b : b / bigb
    scale = _m_dot(t, r) * biga / _m_dot(r, r)
    return scale * r
end

"Component of `a` perpendicular to `b`, with the scaling of `_m_proj`; a zero `b` gives `a`, as `vperp_` does."
@inline function _m_perp(a::SVector{3}, b::SVector{3}, recip::Bool)
    biga = _m_maxabs(a)
    bigb = _m_maxabs(b)
    biga == 0.0 && return zero(a) - zero(b)
    bigb == 0.0 && return a - zero(b)
    t = recip ? (1.0 / biga) * a : a / biga
    r = recip ? (1.0 / bigb) * b : b / bigb
    # `vproj_c` and `vproj_` scale their arguments again; t and r have largest
    # component 1, so that scaling changes nothing.
    v = _m_proj(t, r, recip)
    return biga * (t - v)
end

"Unit cross product as `ucrss`: each vector is scaled by its largest component first."
@inline function _m_ucrss(v1::SVector{3}, v2::SVector{3})
    m1 = _m_maxabs(v1)
    m2 = _m_maxabs(v2)
    t1 = m1 != 0.0 ? v1 / m1 : zero(v1)
    t2 = m2 != 0.0 ? v2 / m2 : zero(v2)
    return _m_hat(_m_cross(t1, t2))
end

"Rotation of `v` about `axis` by `theta`, as `vrotv_`."
@inline function _m_vrotv(v::SVector{3}, axis::SVector{3}, theta)
    primal(_m_norm(axis)) == 0.0 && return v - zero(axis) * zero(theta)
    x = _m_hat(axis)
    p = _m_proj(v, x, false)
    v1 = v - p
    v2 = _m_cross(x, v1)
    s, c = sincos(theta)
    return (c * v1 + s * v2) + p
end

"""
    _m_mxv(m, v) -> SVector{3}

Return the product of the 3 by 3 matrix `m` and the vector `v`, each component the sum
of its three products from left to right, as `mxv_c`. The StaticArrays product `m * v`
lets the compiler fuse a multiplication and an addition, and compilers differ in where
they do so; this result does not depend on the compiler. `_m_mtxv`, `_m_mxm`, `_m_mtxm`
and `_m_mxmt` follow the same rule, and every 3 by 3 product of the package uses them in
place of `*`.
"""
@inline _m_mxv(m::SMatrix{3, 3}, v::SVector{3}) = SVector(
    m[1, 1] * v[1] + m[1, 2] * v[2] + m[1, 3] * v[3],
    m[2, 1] * v[1] + m[2, 2] * v[2] + m[2, 3] * v[3],
    m[3, 1] * v[1] + m[3, 2] * v[2] + m[3, 3] * v[3]
)

@inline _m_mtxv(m::SMatrix{3, 3}, v::SVector{3}) = SVector(
    m[1, 1] * v[1] + m[2, 1] * v[2] + m[3, 1] * v[3],
    m[1, 2] * v[1] + m[2, 2] * v[2] + m[3, 2] * v[3],
    m[1, 3] * v[1] + m[2, 3] * v[2] + m[3, 3] * v[3]
)

@inline _m_rowcol(k::Int) = (k - 1) % 3 + 1, (k - 1) ÷ 3 + 1

@inline _m_mxm(a::SMatrix{3, 3}, b::SMatrix{3, 3}) = SMatrix{3, 3}(
    ntuple(Val(9)) do k
        i, j = _m_rowcol(k)
        a[i, 1] * b[1, j] + a[i, 2] * b[2, j] + a[i, 3] * b[3, j]
    end
)

@inline _m_mtxm(a::SMatrix{3, 3}, b::SMatrix{3, 3}) = SMatrix{3, 3}(
    ntuple(Val(9)) do k
        i, j = _m_rowcol(k)
        a[1, i] * b[1, j] + a[2, i] * b[2, j] + a[3, i] * b[3, j]
    end
)

@inline _m_mxmt(a::SMatrix{3, 3}, b::SMatrix{3, 3}) = SMatrix{3, 3}(
    ntuple(Val(9)) do k
        i, j = _m_rowcol(k)
        a[i, 1] * b[j, 1] + a[i, 2] * b[j, 2] + a[i, 3] * b[j, 3]
    end
)

@inline _m_xpose(m::SMatrix{3, 3}) = SMatrix{3, 3}(
    m[1, 1], m[1, 2], m[1, 3], m[2, 1], m[2, 2],
    m[2, 3], m[3, 1], m[3, 2], m[3, 3]
)

@inline _m_col(m::SMatrix{3, 3}, j::Int) = SVector(m[1, j], m[2, j], m[3, j])
@inline _m_row(m::SMatrix{3, 3}, i::Int) = SVector(m[i, 1], m[i, 2], m[i, 3])

@inline _m_cols(a::SVector{3}, b::SVector{3}, c::SVector{3}) =
    SMatrix{3, 3}(a[1], a[2], a[3], b[1], b[2], b[3], c[1], c[2], c[3])

@inline _m_rows(a::SVector{3}, b::SVector{3}, c::SVector{3}) =
    SMatrix{3, 3}(a[1], b[1], c[1], a[2], b[2], c[2], a[3], b[3], c[3])

"Determinant with the operation order of `det_c`."
@inline _m_det(m::SMatrix{3, 3}) =
    (m[1, 1] * (m[2, 2] * m[3, 3] - m[3, 2] * m[2, 3])) -
    (m[1, 2] * (m[2, 1] * m[3, 3] - m[3, 1] * m[2, 3])) +
    (m[1, 3] * (m[2, 1] * m[3, 2] - m[3, 1] * m[2, 2]))

"The index triple `(i1, i2, i3)` of `rotate_` for any integer axis (taken modulo 3)."
@inline function _m_axes(iaxis::Integer)
    t = mod(Int(iaxis), 3)
    return t == 1 ? (1, 2, 3) : t == 2 ? (2, 3, 1) : (3, 1, 2)
end

@inline function _m_rotate(angle::T, iaxis::Integer) where {T <: Real}
    s, c = sincos(angle)
    i1, i2, i3 = _m_axes(iaxis)
    z = zero(s)
    o = one(s)
    return SMatrix{3, 3}(
        ntuple(Val(9)) do k
            i, j = _m_rowcol(k)
            (i == i1 || j == i1) ? (i == j ? o : z) :
                i == j ? c : (i == i2 ? s : -s)
        end
    )
end

@inline function _m_rotmat(m::SMatrix{3, 3}, angle::Real, iaxis::Integer)
    s, c = sincos(angle)
    i1, i2, i3 = _m_axes(iaxis)
    return SMatrix{3, 3}(
        ntuple(Val(9)) do k
            i, j = _m_rowcol(k)
            i == i1 ? m[i1, j] - zero(c) :
                i == i2 ? c * m[i2, j] + s * m[i3, j] : -s * m[i2, j] + c * m[i3, j]
        end
    )
end

"`[a3]_x3 * [a2]_x2 * [a1]_x1` with the arithmetic of `eul2m_`, in the number type of the angles; the axes are not checked and are taken modulo 3."
@inline _m_eul2m_unchecked(a3, a2, a1, x3::Integer, x2::Integer, x1::Integer) =
    _m_rotmat(_m_rotmat(_m_rotate(a1, x1), a2, x2), a3, x3)

"Inverse of a matrix with orthogonal columns, with the checks of `invort_`: the transpose with each row divided by the squared length of its column."
@inline function _m_invort(m::SMatrix{3, 3})
    unit(i) = begin
        u, len = _m_unorm(_m_col(m, i))
        l0 = primal(len)
        l0 == 0.0 && throw(SpiceError(:zerolengthcolumn, "column $i of the matrix has norm zero"))
        (l0 < 1.0 && l0 * floatmax(Float64) < 1.0) && throw(
            SpiceError(
                :columntoosmall,
                "the length $l0 of column $i cannot be inverted"
            )
        )
        (1.0 / len) * u
    end
    return _m_rows(unit(1), unit(2), unit(3))
end

"True if `m` is a rotation within the tolerances, as `isrot_c` (no argument checks)."
@inline function _m_isrot(m::SMatrix{3, 3}, ntol::Real, dtol::Real)
    u1, n1 = _m_unorm(_m_col(m, 1))
    u2, n2 = _m_unorm(_m_col(m, 2))
    u3, n3 = _m_unorm(_m_col(m, 3))
    inside(x, tol) = (p = primal(x); 1.0 - tol <= p <= 1.0 + tol)
    normok = inside(n1, ntol) && inside(n2, ntol) && inside(n3, ntol)
    # `isrot_c` takes the determinant of the matrix whose rows are the unit columns.
    return normok && inside(_m_det(_m_rows(u1, u2, u3)), dtol)
end

"""
Unit state of `dvhat`: the unit position and its derivative; a zero position gives the
zero position and the velocity as it is. With `recip` the velocity is scaled as
`dvhat_c` does (through `vperp_c`); without it by the division of the f2c `dvhat_`, which
`ducrss_`, `dvsep_`, `zztwovxf_` and `surfpv_` call. The two differ by rounding, and only
the division keeps a velocity below about 5.6e-309 finite.
"""
@inline function _m_dvhat(s::SVector{6}, recip::Bool = true)
    p, v = _m_pv(s)
    u, len = _m_unorm(p)
    primal(len) == 0.0 && return _m_state(u, v - zero(u))
    return _m_state(u, (1.0 / len) * _m_perp(v, u, recip))
end

"Cross product of the positions of two states and its derivative, as `dvcrss_c`."
@inline function _m_dvcrss(s1::SVector{6}, s2::SVector{6})
    p1, v1 = _m_pv(s1)
    p2, v2 = _m_pv(s2)
    return _m_state(_m_cross(p1, p2), _m_cross(v1, p2) + _m_cross(p1, v2))
end

"Unit cross product of two states and its derivative, as `ducrss_`."
@inline function _m_ducrss(s1::SVector{6}, s2::SVector{6})
    f1 = _m_maxabs(_m_pv(s1)[1])
    f2 = _m_maxabs(_m_pv(s2)[1])
    a = f1 > 0.0 ? (1.0 / f1) * s1 : s1
    b = f2 > 0.0 ? (1.0 / f2) * s2 : s2
    return _m_dvhat(_m_dvcrss(a, b), false)
end

"`brcktd_`: `x` bracketed between the two ends `e1` and `e2`, in either order."
@inline _m_brcktd(x::Float64, e1::Float64, e2::Float64) =
    e1 < e2 ? max(e1, min(e2, x)) : max(e2, min(e1, x))

"Distance between two vectors of equal length, as `vdistg_`; the largest difference is taken with the `max` of f2c, which keeps the previous value at a NaN."
function _m_vdistg(v1::AbstractVector, v2::AbstractVector)
    scale = 0.0
    @inbounds for i in eachindex(v1, v2)       # invariant: scale bounds the differences seen
        d = abs(primal(_m_f(v1[i]) - _m_f(v2[i])))
        scale = scale >= d ? scale : d
    end
    dist = zero(_m_f(zero(eltype(v1))) - _m_f(zero(eltype(v2)))) / 1.0
    scale == 0.0 && return dist
    @inbounds for i in eachindex(v1, v2)       # invariant: dist is the sum of squares seen
        t = (_m_f(v1[i]) - _m_f(v2[i])) / scale
        dist += t * t
    end
    return scale * sqrt(dist)
end

"Angle between two unit vectors `u1` and `u2`, with the two-branch arcsine of `vsep`."
@inline function _m_usep(u1::SVector{3}, u2::SVector{3})
    d = primal(_m_dot(u1, u2))
    if d > 0.0
        return 2.0 * asin(0.5 * _m_norm(u1 - u2))
    elseif d < 0.0
        return Float64(pi) - 2.0 * asin(0.5 * _m_norm(u1 + u2))
    end
    return Float64(pi) / 2 - zero(_m_dot(u1, u2))
end

"Angle between `v1` and `v2` as `vsep_`: zero if one of them is the zero vector."
@inline function _m_vsep(v1::SVector{3}, v2::SVector{3})
    u1, m1 = _m_unorm(v1)
    u2, m2 = _m_unorm(v2)
    (primal(m1) == 0.0 || primal(m2) == 0.0) && return zero(asin(_m_dot(u1, u2)))
    return _m_usep(u1, u2)
end

# --- quaternions and axis-angle rotations ---------------------------------------------

"""
    _m_q2m(q0, q1, q2, q3, always::Bool) -> SMatrix{3,3}

Rotation matrix of the SPICE quaternion `(q0, q1, q2, q3)`, as `q2m_`: the products are
divided by the squared length unless it is 0 or 1. With `always` a number type other
than `Float64` is divided also when the plain part is 1, so that the derivatives of the
length cancel.
"""
@inline function _m_q2m(q0, q1, q2, q3, always::Bool)
    q01, q02, q03 = q0 * q1, q0 * q2, q0 * q3
    q12, q13, q23 = q1 * q2, q1 * q3, q2 * q3
    q1s, q2s, q3s = q1 * q1, q2 * q2, q3 * q3
    l2 = q0 * q0 + q1s + q2s + q3s
    pl = primal(l2)
    if pl != 0.0 && (pl != 1.0 || (always && !(l2 isa Float64)))
        f = 1.0 / l2
        q01 *= f; q02 *= f; q03 *= f
        q12 *= f; q13 *= f; q23 *= f
        q1s *= f; q2s *= f; q3s *= f
    end
    return SMatrix{3, 3}(
        1.0 - (q2s + q3s) * 2.0, (q12 + q03) * 2.0, (q13 - q02) * 2.0,
        (q12 - q03) * 2.0, 1.0 - (q1s + q3s) * 2.0, (q23 + q01) * 2.0,
        (q13 + q02) * 2.0, (q23 - q01) * 2.0, 1.0 - (q1s + q2s) * 2.0
    )
end

"Unit quaternion of the rotation matrix `r`, scalar part first and not negative, as `m2q_` without its rotation check."
@inline function _m_m2q_unchecked(r::SMatrix{3, 3})
    trace = r[1] + r[5] + r[9]
    mtrace = 1.0 - trace
    cc4 = 1.0 + trace
    s114 = mtrace + 2.0 * r[1]
    s224 = mtrace + 2.0 * r[5]
    s334 = mtrace + 2.0 * r[9]
    # The first of the four squares that is at least 1/4 gives the divisor; the order
    # of the tests is the one of `m2q_`.
    if 1.0 <= primal(cc4)
        c = sqrt(cc4 * 0.25)
        factor = 1.0 / (c * 4.0)
        s1 = (r[6] - r[8]) * factor
        s2 = (r[7] - r[3]) * factor
        s3 = (r[2] - r[4]) * factor
    elseif 1.0 <= primal(s114)
        s1 = sqrt(s114 * 0.25)
        factor = 1.0 / (s1 * 4.0)
        c = (r[6] - r[8]) * factor
        s2 = (r[4] + r[2]) * factor
        s3 = (r[7] + r[3]) * factor
    elseif 1.0 <= primal(s224)
        s2 = sqrt(s224 * 0.25)
        factor = 1.0 / (s2 * 4.0)
        c = (r[7] - r[3]) * factor
        s1 = (r[4] + r[2]) * factor
        s3 = (r[8] + r[6]) * factor
    else
        s3 = sqrt(s334 * 0.25)
        factor = 1.0 / (s3 * 4.0)
        c = (r[2] - r[4]) * factor
        s1 = (r[7] + r[3]) * factor
        s2 = (r[8] + r[6]) * factor
    end
    l2 = c * c + s1 * s1 + s2 * s2 + s3 * s3
    if primal(l2) != 1.0
        polish = 1.0 / sqrt(l2)
        c *= polish
        s1 *= polish
        s2 *= polish
        s3 *= polish
    end
    return primal(c) > 0.0 ? SVector(c, s1, s2, s3) : SVector(-c, -s1, -s2, -s3)
end

"Throw `:notarotation` unless `isrot(r, 0.1, 0.1)`, the check of `m2q_` and `m2eul_`."
@inline function _m_checkrot(r::SMatrix{3, 3})
    _m_isrot(r, 0.1, 0.1) || throw(SpiceError(:notarotation, "the input matrix is not a rotation"))
    return nothing
end

@inline function _m_m2q(r::SMatrix{3, 3})
    _m_checkrot(r)
    return _m_m2q_unchecked(r)
end

"Axis and angle in [0, pi] of the rotation with the quaternion `q`, as `raxisa_`; the identity gives the axis `(0, 0, 1)`."
@inline function _m_raxisa_q(q::SVector{4})
    v = SVector(q[2], q[3], q[4])
    if _m_iszero(v)
        return SVector(zero(q[1]), zero(q[1]), one(q[1])), zero(atan(q[1]))
    elseif primal(q[1]) == 0.0
        return v, Float64(pi) - zero(atan(q[1]))
    end
    return _m_hat(v), 2.0 * atan(_m_norm(v), q[1])
end

"Matrix that rotates vectors by `angle` about `axis`, as `axisar_`; a zero axis gives the identity."
@inline _m_axisar(axis::SVector{3}, angle) = _m_cols(
    _m_vrotv(SVector(1.0, 0.0, 0.0), axis, angle), _m_vrotv(SVector(0.0, 1.0, 0.0), axis, angle),
    _m_vrotv(SVector(0.0, 0.0, 1.0), axis, angle)
)

"Quotient `a / b` with the checks of `zzdiv_`: `:dividebyzero`, `:numericoverflow`, and zero for a quotient that underflows."
@inline function _m_div(a::Real, b::Real)
    a0 = primal(a)
    b0 = primal(b)
    b0 == 0.0 && throw(SpiceError(:dividebyzero, "numerical divide by zero"))
    a0 == 0.0 && return zero(a) * zero(b)
    expnt = trunc(log10(floatmax(Float64)))
    logn = log10(abs(a0))
    logd = log10(abs(b0))
    logn - logd > expnt && throw(SpiceError(:numericoverflow, "numerical overflow in a quotient"))
    logn - logd < -(expnt - 1.0) && return zero(a) * zero(b)
    return a / b
end

# --- coordinates, ellipsoids, lines and planes -----------------------------------------
#
# The kernels of `math/coordinates.jl` and `math/ellipsoid.jl` that the DSK files, which
# come before them, use too.

"""
`(radius, lon, lat)` of `v`. With `scaled` the longitude comes from the components
divided by the largest one, as in `reclat_c`; without it from the components themselves,
as in `reclat_`, which also decide whether the point is on the z axis (`[0, 5e-324, 1e10]`
has the longitude pi/2 in `reclat_` and 0 in `reclat_c`).
"""
@inline function _m_reclat(v::SVector{3}, scaled::Bool)
    vmax = _m_maxabs(v)
    z0 = zero(atan(v[1] / 1.0))
    vmax > 0.0 || return z0, z0, z0
    x1 = v[1] / vmax
    y1 = v[2] / vmax
    z1 = v[3] / vmax
    radius = vmax * sqrt(x1 * x1 + y1 * y1 + z1 * z1)
    lat = atan(z1, sqrt(x1 * x1 + y1 * y1))
    if scaled
        (primal(x1) == 0.0 && primal(y1) == 0.0) && return radius + z0, z0, lat
        return radius + z0, atan(y1, x1), lat
    end
    (primal(v[1]) == 0.0 && primal(v[2]) == 0.0) && return radius + z0, z0, lat
    return radius + z0, atan(v[2], v[1]), lat
end

@inline function _m_latrec(radius, lon, lat)
    slon, clon = sincos(lon)
    slat, clat = sincos(lat)
    return SVector(radius * clon * clat, radius * slon * clat, radius * slat)
end

@inline function _m_cylrec(r, clon, z)
    s, c = sincos(clon)
    return SVector(r * c, r * s, z - zero(r * c))
end

@noinline _m_badre(re, kind::Symbol = :valueoutofrange) = throw(
    SpiceError(
        kind,
        "the equatorial radius must be positive; it is $(primal(re))"
    )
)
@noinline _m_badflat(f) = throw(
    SpiceError(
        :valueoutofrange,
        "the flattening coefficient must be less than 1; it is $(primal(f))"
    )
)

@inline function _m_georec(lon, lat, alt, re, f)
    primal(re) > 0.0 || _m_badre(re)
    primal(f) < 1.0 || _m_badflat(f)
    rp = re - f * re
    sphi, cphi = sincos(lat)
    slmbda, clmbda = sincos(lon)
    big = max(abs(primal(re * cphi)), abs(primal(rp * sphi)))
    x = re * cphi / big
    y = rp * sphi / big
    scale = 1.0 / (big * sqrt(x * x + y * y))
    base = SVector(
        scale * re * re * clmbda * cphi, scale * re * re * slmbda * cphi,
        scale * rp * rp * sphi
    )
    normal = _m_surfnm(re, re, rp, base)
    return 1.0 * base + alt * normal
end

@inline _m_abs(x::Real) = primal(x) < 0.0 ? -x : x

@inline _m_max(x, y) = primal(y) > primal(x) ? y : x

@noinline _m_badaxis(a, b, c) = throw(
    SpiceError(
        :badaxislength,
        "the ellipsoid semi-axis lengths must be positive; they are $(primal(a)), $(primal(b)), $(primal(c))"
    )
)

@inline _m_axes_ok(a, b, c) = primal(a) > 0.0 && primal(b) > 0.0 && primal(c) > 0.0

@inline function _m_surfnm(a, b, c, point::SVector{3})
    _m_axes_ok(a, b, c) || _m_badaxis(a, b, c)
    m = primal(a) <= primal(b) ? a : b
    m = primal(m) <= primal(c) ? m : c
    a1 = m / a
    b1 = m / b
    c1 = m / c
    return _m_hat(SVector(point[1] * (a1 * a1), point[2] * (b1 * b1), point[3] * (c1 * c1)))
end

"`(found, point)` of `surfpt_`, with the semi-axes `a`, `b`, `c`; the point is zero when nothing is found."
@inline function _m_surfpt(positn::SVector{3}, u::SVector{3}, a, b, c)
    _m_iszero(u) && throw(SpiceError(:zerovector, "the ray direction is the zero vector"))
    _m_axes_ok(a, b, c) || _m_badaxis(a, b, c)
    x = SVector(u[1] / a, u[2] / b, u[3] / c)
    y = SVector(positn[1] / a, positn[2] / b, positn[3] / c)
    # In the scaled space the ellipsoid is the unit sphere; p is the point of the
    # line that is nearest to the center.
    p = _m_perp(y, x, false)
    yproj = y - p
    ymag = _m_norm(y)
    pmag = _m_norm(p)
    ux = _m_hat(x)
    none = zero(p)
    sign = 1.0
    if primal(ymag) > 1.0
        primal(pmag) > 1.0 && return false, none
        primal(_m_dot(yproj, x)) > 0.0 && return false, none
        primal(pmag) == 1.0 && return true, SVector(p[1] * a, p[2] * b, p[3] * c)
        sign = -1.0
    elseif primal(ymag) == 1.0
        return true, convert(typeof(none), positn)  # `vequ_`: keeps the sign of a zero
    end
    t = 1.0 - pmag * pmag
    scale = primal(t) <= 0.0 ? zero(t) : sqrt(t)  # `sqrt(max(0, t))` of f2c, which keeps a NaN
    q = 1.0 * p + (sign * scale) * ux
    return true, SVector(q[1] * a, q[2] * b, q[3] * c)
end

"`edpnt_` after its checks of the axes and of a zero `v`: `v` scaled to the ellipsoid."
@inline function _m_edpnt(v::SVector{3}, a, b, c)
    x = v[1] / a
    y = v[2] / b
    z = v[3] / c
    level = x * x + y * y + z * z
    primal(level) <= 0.0 && throw(
        SpiceError(
            :pointtoosmall,
            "the input point is too close to the origin to be scaled to the ellipsoid"
        )
    )
    sq = sqrt(level)
    return SVector(v[1] / sq, v[2] / sq, v[3] / sq)
end

"""
Nearest point of the line through `lp` along a nonzero `ld` to `p`, and its distance;
`recip` true projects as `nplnpt_c`, false as the f2c `nplnpt_`.
"""
@inline function _m_nplnpt(lp::SVector{3}, ld::SVector{3}, p::SVector{3}, recip::Bool)
    pnear = _m_proj(p - lp, ld, recip) + lp
    return pnear, _m_norm(pnear - p)
end

"`npsgpt`: the point of the segment from `e1` to `e2` that is nearest to `p`, and its distance."
@inline function _m_npsgpt(e1::SVector{3}, e2::SVector{3}, p::SVector{3})
    seg = e2 - e1
    _m_iszero(seg) && return e1, _m_norm(e1 - p)
    # `npsgpt_` calls the Fortran `nplnpt_`, whose projection is `vproj_`.
    lnear, dist = _m_nplnpt(e1, seg, p, false)
    offdot = _m_dot(lnear - e1, seg)
    primal(offdot) < 0.0 && return e1, _m_norm(e1 - p)
    primal(offdot) > primal(_m_dot(seg, seg)) && return e2, _m_norm(e2 - p)
    return lnear, dist
end


# --- planes as (unit normal, constant) ------------------------------------------------

"The plane `(normal, constant)` with the signs of both changed if `constant < 0`."
@inline _m_nvc(normal::SVector{3}, constant) =
    primal(constant) < 0.0 ? (-normal, -constant) : (normal, constant)

"`nvc2pl_` as a tuple `(normal, constant)`; throws `:zerovector` for a zero `normal`."
@inline function _m_nvc2pl(normal::SVector{3}, konst)
    u, mag = _m_unorm(normal)
    primal(mag) == 0.0 && throw(SpiceError(:zerovector, "the normal of a plane must not be zero"))
    return _m_nvc(u, konst / mag)
end

"`nvp2pl_` as a tuple `(normal, constant)`, without its check: a zero `normal` gives `(0, 0)`."
@inline function _m_nvp2pl(normal::SVector{3}, point::SVector{3})
    u = _m_hat(normal)
    # The zero (which gives the number type of `point` to the normal) is subtracted, not
    # added, because `-0.0 + 0.0` is `0.0`: a zero component keeps its sign, as in CSPICE.
    return _m_nvc(u - zero(point), _m_dot(point, u))
end

"""
`(nxpts, xpt)` of `inrypl_` for the ray from `vtx` along `d` and the plane of the unit
normal `normal` and the constant `constant >= 0`.
"""
@inline function _m_inrypl(vtx::SVector{3}, d::SVector{3}, normal::SVector{3}, constant)
    toobig = floatmax(Float64) / 3.0
    vnrm = _m_norm(vtx)
    primal(vnrm) >= toobig && throw(
        SpiceError(
            :vectortoobig,
            "the vertex of the ray is too far from the origin"
        )
    )
    primal(constant) >= toobig && throw(
        SpiceError(
            :vectortoobig,
            "the plane is too far from the origin"
        )
    )
    udir = _m_hat(d)
    _m_iszero(udir) && throw(SpiceError(:zerovector, "the direction of the ray is the zero vector"))
    mscale = _m_max(_m_abs(constant) - zero(vnrm), _m_abs(vnrm) - zero(constant))
    if primal(mscale) != 0.0
        sclvtx = (1.0 / mscale) * vtx
        sclcon = constant / mscale
    else
        sclvtx = vtx - zero(SVector{3, typeof(mscale)})
        sclcon = constant - zero(mscale)
    end
    primal(mscale) > 1.0 && (toobig = toobig / primal(mscale))
    prjdif = sclcon - _m_dot(sclvtx, normal)
    prjdir = _m_dot(udir, normal)
    none = zero(1.0 * sclvtx + prjdif * udir)
    if primal(prjdif) == 0.0
        return (primal(_m_dot(normal, udir)) == 0.0 ? -1 : 1), convert(typeof(none), vtx)
    end
    pd, pf = primal(prjdir), primal(prjdif)
    samesign = (pd > 0.0 && pf > 0.0) || (pd < 0.0 && pf < 0.0)
    (!samesign || abs(pf) >= abs(pd) * toobig) && return 0, none
    scale = _m_abs(prjdif) / _m_abs(prjdir)
    return 1, mscale * (1.0 * sclvtx + scale * udir)
end

# ---------------------------------------------------------------------------------
# Integer conversions of C
# ---------------------------------------------------------------------------------

"`x` cut toward zero as the C cast `(integer) x` on x86; the smallest 32-bit integer for NaN or out of range."
@inline _cint(x::Float64)::Int =
    (x > -2147483649.0 && x < 2147483648.0) ? Int(trunc(Int32, x)) : Int(typemin(Int32))

"`n` reduced to the range of a 32-bit integer, as the integer arithmetic of the C code."
@inline _wrap32(n::Int)::Int = Int(n % Int32)

"Integer nearest to `x`, a tie away from zero, as f2c `i_dnnt`."
@inline _nint(x::Float64)::Int = _cint(x >= 0.0 ? floor(x + 0.5) : -floor(0.5 - x))

"`x` rounded to a whole number, a tie away from zero, as f2c `d_nint`."
@inline _dnint(x::Float64)::Float64 = x >= 0.0 ? floor(x + 0.5) : -floor(0.5 - x)

# ---------------------------------------------------------------------------------
# Euler angles and orthonormal columns
# ---------------------------------------------------------------------------------

@noinline _m_badaxes(a, b, c) =
    throw(SpiceError(:badaxisnumbers, "axis numbers are $a, $b, $c; each must be 1, 2 or 3"))

"Throw `:badaxisnumbers` unless each of the axes `a`, `b`, `c` is 1, 2 or 3."
@inline function _m_checkaxes(a::Integer, b::Integer, c::Integer)
    (1 <= a <= 3 && 1 <= b <= 3 && 1 <= c <= 3) || _m_badaxes(a, b, c)
    return nothing
end

"`[angle3]_axis3 * [angle2]_axis2 * [angle1]_axis1` as `eul2m_`, with its axis check."
@inline function _m_eul2m(angle3, angle2, angle1, axis3::Integer, axis2::Integer, axis1::Integer)
    _m_checkaxes(axis3, axis2, axis1)
    return _m_eul2m_unchecked(angle3, angle2, angle1, axis3, axis2, axis1)
end

"`m` with orthonormal columns as `sharpr_`: column 1 scaled, column 3 the unit cross product of columns 1 and 2, column 2 that of columns 3 and 1."
@inline function _m_sharpr(m::SMatrix{3, 3})
    c1 = _m_hat(_m_col(m, 1))
    c3 = _m_ucrss(c1, _m_col(m, 2))
    c2 = _m_ucrss(c3, c1)
    return _m_cols(c1, c2, c3)
end

"The upper left and lower left 3 by 3 blocks of a state transformation: the rotation and its derivative."
@inline function _m_blocks(x::SMatrix{6, 6})
    r = SMatrix{3, 3}(ntuple(k -> x[(k - 1) % 3 + 1, (k - 1) ÷ 3 + 1], Val(9)))
    d = SMatrix{3, 3}(ntuple(k -> x[(k - 1) % 3 + 4, (k - 1) ÷ 3 + 1], Val(9)))
    return r, d
end
