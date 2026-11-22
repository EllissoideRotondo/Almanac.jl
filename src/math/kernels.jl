# Private vector and matrix kernels of the files in `src/math/` (prefix `_m_`): input
# conversion, 3-vector and 3 by 3 matrix arithmetic in the operation order of SPICELIB,
# rotations and the derivatives of unit vectors. They are included early, so that every
# file can use them; the public calls are in `math/vectors.jl` and the files after it.

# --- constants of SPICELIB ---------------------------------------------------------

"`pi_()`, `twopi_()` and `halfpi_()` of CSPICE as plain doubles."
const _PI = Float64(pi)
const _TWOPI = 2 * Float64(pi)
const _HALFPI = Float64(pi) / 2

"The ID of the frame J2000."
const _J2000_FRAME = Int32(1)

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

"""
Projection of `a` on `b`. With `recip` the scaling is a product with the reciprocal of
the largest component (
[`vproj_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vproj_c.html)); without it the scaling is a division (`vproj_`).
"""
@inline function _m_proj(a::SVector{3}, b::SVector{3}, recip::Bool)
    biga = _m_maxabs(a)
    bigb = _m_maxabs(b)
    (biga == 0.0 || bigb == 0.0) && return zero(a) - zero(b)
    t = recip ? (1.0 / biga) * a : a / biga
    r = recip ? (1.0 / bigb) * b : b / bigb
    scale = _m_dot(t, r) * biga / _m_dot(r, r)
    return scale * r
end

"""
Component of `a` perpendicular to `b`, with the scaling choice of `_m_proj`. A zero `b`
gives `a`, as `vperp_` of SPICELIB and the text of
[`vperp_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vperp_c.html) say.
"""
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

Contract: the product of the 3 by 3 matrix `m` and the vector `v`, each component the
sum of its three products from left to right, as
[`mxv_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mxv_c.html). The product `m * v` of
StaticArrays lets the compiler join a multiplication and an addition into one fused
operation, and compilers differ in where they do so; the result of this function does
not depend on the compiler, which matters where a component is a small difference of
large products.
The forms `_m_mtxv`, `_m_mxm`, `_m_mtxm` and `_m_mxmt` obey the same rule; the package uses
them in place of `*` for each 3 by 3 product.
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

"""
Contract: `[a3]_x3 * [a2]_x2 * [a1]_x1` with the arithmetic of `eul2m_` (`rotate_` of
`a1`, then `rotmat_` of `a2` and of `a3`), in the number type of the angles. The axes
are not checked; an axis outside 1:3 is taken modulo 3, as `rotate_` takes it.
"""
@inline _m_eul2m_unchecked(a3, a2, a1, x3::Integer, x2::Integer, x1::Integer) =
    _m_rotmat(_m_rotmat(_m_rotate(a1, x1), a2, x2), a3, x3)

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
Unit state of `dvhat`: the unit position and its derivative. With `recip` the
perpendicular part of the velocity is scaled by a product with a reciprocal, as the
C function
[`dvhat_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/dvhat_c.html) does (through `vperp_c`). Without it the scaling is a division,
as in the f2c routine `dvhat_` (through `vperp_`), which is the routine that `ducrss_`,
`dvsep_`, `zztwovxf_` and `surfpv_` call. The two forms differ by rounding, and the
reciprocal of a velocity below about 5.6e-309 is infinite, so only the division keeps
such a velocity finite. Postcondition: a zero position gives the zero position and the
velocity as it is.
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

"""
Quotient `a / b` with the checks of `zzdiv_`: `SpiceError(:dividebyzero)` for a zero
`b`, `SpiceError(:numericoverflow)` if the quotient is above the largest number, and
zero if it is below the smallest one.
"""
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

# ---------------------------------------------------------------------------------
# Integer conversions of C
# ---------------------------------------------------------------------------------

"Contract: returns `x` cut to an integer as the C cast `(integer) x` does on x86: toward zero, and the smallest 32-bit integer when `x` is NaN or out of range."
@inline _cint(x::Float64)::Int =
    (x > -2147483649.0 && x < 2147483648.0) ? Int(trunc(Int32, x)) : Int(typemin(Int32))

"Contract: returns `n` reduced to the range of a 32-bit integer, as integer arithmetic of the C code does."
@inline _wrap32(n::Int)::Int = Int(n % Int32)

"Contract: returns the integer nearest to `x`, a tie away from zero, as f2c `i_dnnt`."
@inline _nint(x::Float64)::Int = _cint(x >= 0.0 ? floor(x + 0.5) : -floor(0.5 - x))

"Contract: returns `x` rounded to a whole number, a tie away from zero, as f2c `d_nint`."
@inline _dnint(x::Float64)::Float64 = x >= 0.0 ? floor(x + 0.5) : -floor(0.5 - x)
