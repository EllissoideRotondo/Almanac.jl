# Brouwer-Lyddane mean elements (short-period terms only) of an `Orbit`: the port of
# src/astro/orbit_mean_elements.rs of ANISE (b3aaadc), itself a port of GMAT's
# StateConversionUtil::CartesianToBrouwerMeanShort. Angles are in rad here (ANISE: degrees),
# so the values differ from ANISE by the rounding of the conversions only. Everything that
# decides (a branch, a convergence test) does so on `primal`, as the rest of the package.

# The Earth's J2. The mean elements are those of the Earth only: a center other than
# body 399 is an error, as in ANISE.
const J2_EARTH = 1.082_626_925_638_815e-3

_brouwer_fail(msg::AbstractString) = throw(
    SpiceError(
        :degeneratecase,
        "the Brouwer mean elements are undefined: $msg"
    )
)

# The conversions below are also those of `keplerian_orbit`, `with_element` and `two_body`
# (build.jl), so their errors do not name the Brouwer elements.
_conv_fail(msg::AbstractString) = throw(SpiceError(:degeneratecase, "undefined: $msg"))

# ANISE `mean_anomaly_to_true_anomaly_rad` (src/astro/utils.rs): GTDS MathSpec 3-180, 3-181,
# 3-186 by Newton iteration. Its `MathError`s are `:degeneratecase`.
const _MA_EPSILON = 1.0e-12

# The iteration runs on the primal values, so the Float64 result is that of ANISE bit for bit.
# The derivative is carried by one more Newton step in the number type T at the converged value
# `x0`: `step` is the correction of that step, and `x0 + (step - primal(step))` has the primal
# value of `x0` and the exact first derivative of the solution (implicit function theorem). A
# Newton step doubles the number of exact orders, so the step is taken twice: orders 1 to 3, which
# nested ForwardDiff (Hessians) needs up to the second. The residual `r0` of the converged value
# (up to 1e-12) is taken out of the steps: they solve the equation for M - r0, whose root is
# exactly `x0`, so the derivative carries no residual times the derivative of the denominator
# (that was the error of d nu/de, 4.7e-12). The old loop returned an iterate whose dual part lagged the solution by one step, a derivative
# error up to 5e-8. A plain Float64 has no derivative, and `_carry` returns `x0`.
_carry(x0::Float64, step::Float64) = x0
_carry(x0, step) = x0 + (step - primal(step))

function _ma_to_ta(rm::T, ecc::T) where {T}
    return _ea_to_ta(_ma_to_ea(rm, ecc), ecc)
end

"""
    _ma_to_ea(rm, ecc) -> anomaly

The eccentric anomaly of the mean anomaly `rm` (ecc <= 1), or the hyperbolic anomaly
(ecc > 1): the Newton iteration of `_ma_to_ta`, with the derivative carried. The eccentric
anomaly is near `rm`, with 2pi added when it is negative; it is not reduced to [0, 2pi).
"""
function _ma_to_ea(rm::T, ecc::T) where {T}
    fail(what) = _conv_fail("the mean to true anomaly conversion $what")
    rp, ep = primal(rm), primal(ecc)
    if ep <= 1
        e2 = rp + ep * sin(rp)
        iter = 0
        while true
            iter += 1
            iter > 1000 && fail("did not converge in 1000 iterations")
            normalized = 1 - ep * cos(e2)
            abs(normalized) < _MA_EPSILON && fail("has a normalized anomaly below 1e-12")
            e1 = e2 - (e2 - ep * sin(e2) - rp) / normalized
            abs(e2 - e1) < _MA_EPSILON && break
            e2 = e1
        end
        r0 = rp - (e2 - ep * sin(e2))
        e2 = _carry(e2, (rm - (e2 - ecc * sin(e2)) - r0) / (1 - ecc * cos(e2)))
        e2 = _carry(e2, (rm - (e2 - ecc * sin(e2)) - r0) / (1 - ecc * cos(e2)))
        return primal(e2) < 0 ? e2 + 2 * T(π) : e2
    end
    f2 = 0.0
    iter = 0
    while true
        iter += 1
        iter > 1000 && fail("did not converge in 1000 iterations")
        normalizer = ep * cosh(f2) - 1
        abs(normalizer) < _MA_EPSILON && fail("has a normalized anomaly below 1e-12")
        f1 = f2 - (ep * sinh(f2) - f2 - rp) / normalizer
        abs(f2 - f1) < _MA_EPSILON && break
        f2 = f1
    end
    r0 = rp - (ep * sinh(f2) - f2)
    f2 = _carry(f2, (rm - (ecc * sinh(f2) - f2) - r0) / (ecc * cosh(f2) - 1))
    return _carry(f2, (rm - (ecc * sinh(f2) - f2) - r0) / (ecc * cosh(f2) - 1))
end

"""
    _ea_to_ta(e, ecc) -> anomaly

The true anomaly of the eccentric (ecc <= 1) or hyperbolic anomaly `e`, as ANISE takes it.

On an ellipse the plain value is that of ANISE, `2 atan(sqrt((1 + e)/(1 - e)) tan(E/2))`
wrapped to [0, 2pi), or `E` itself where |E - pi| < 1e-8. Its derivatives are taken from
`2 atan(sqrt(1 + e) sin(E/2), sqrt(1 - e) cos(E/2))`, the same function without the pole
of `tan` at every odd multiple of pi: near it the derivatives of ANISE's formula lose up
to 1e-8 of d2 nu/dM2 relative to its bound, and at an odd multiple other than pi (`E` is
not reduced) they are wrong by O(1).
"""
function _ea_to_ta(e::S, ecc::U) where {S, U}
    T = promote_type(S, U)
    fail(what) = _conv_fail("the mean to true anomaly conversion $what")
    if primal(ecc) <= 1
        ep, eccp = primal(e), primal(ecc)
        ta = if abs(ep - π) >= 1.0e-8
            normalized = 1 - eccp
            abs(normalized) < _MA_EPSILON && fail("has a normalized anomaly below 1e-12")
            ratio = (1 + eccp) / normalized
            ratio < 0 && fail("has a negative eccentricity ratio")
            2 * atan(sqrt(ratio) * tan(ep / 2))
        else
            ep      # ANISE: nu = E at E = pi
        end
        ta = ta < 0 ? ta + 2 * Float64(π) : ta
        S === Float64 && U === Float64 && return ta
        return _with_primal(2 * atan(sqrt(1 + ecc) * sin(e / 2), sqrt(1 - ecc) * cos(e / 2)), ta)
    end
    normalized = ecc - 1
    abs(primal(normalized)) < _MA_EPSILON && fail("has a normalized anomaly below 1e-12")
    ratio = (ecc + 1) / normalized
    primal(ratio) < 0 && fail("has a negative eccentricity ratio")
    ta = 2 * atan(sqrt(ratio) * tanh(e / 2))
    return primal(ta) < 0 ? ta + 2 * T(π) : ta
end

# `_ma_to_ta` is a derivative point (src/ad.jl): Enzyme, Mooncake and the users of
# ChainRulesCore see Float64 numbers, for which `primal` is the identity, so `_carry` drops
# the step and the branch at E = pi keeps dE/dM. They would differentiate the Newton loop
# (a derivative one step behind the solution) and take dE for dnu at E = pi (wrong by a
# factor up to 6 at ecc = 0.95). Their rules take the derivatives below, those of the
# implicit function M = E - e sin E (M = e sinh F - F) at the converged anomaly, with
# nu(E, e) of the conic: q = 1 - e cos E (e cosh F - 1), s = sqrt(|1 - e^2|),
#   dnu/dM = s / q^2,   dnu/de = sin E (2 - e^2 - e cos E) / (s q^2)
#                       (hyperbola: -sinh F (e^2 - 2 + e cosh F) / (s q^2)).
# They hold at E = pi as well, where ANISE returns E for nu. The tracked numbers of
# ReverseDiff and Tracker strip with `primal`, so `_carry` gives them the same derivative
# without a rule. With a `Dual` (ForwardDiff over a backend) the anomaly carries its exact
# derivatives, and so do the partials.
function _ad_rule(::typeof(_ma_to_ta), rm::Real, ecc::Real)
    m, e = promote(rm, ecc)
    a = _ma_to_ea(m, e)
    ta = _ea_to_ta(a, e)
    if primal(e) <= 1
        sa, ca = sincos(a)
        q = 1 - e * ca
        s = sqrt(1 - e^2)
        return ta, (s / q^2, sa * (2 - e^2 - e * ca) / (s * q^2))
    end
    q = e * cosh(a) - 1
    s = sqrt(e^2 - 1)
    return ta, (s / q^2, -sinh(a) * (e^2 - 2 + e * cosh(a)) / (s * q^2))
end

# ANISE `try_keplerian` (src/astro/orbit.rs): GMAT's KeplerianToCartesian, with its errors
# (parabola, true anomaly beyond the asymptotes, infinite radius, zero parameter).
function _kep_to_cart(sma::T, ecc::T, inc::T, raan::T, aop::T, ta::T, gm::Real) where {T}
    fail(what) = _conv_fail("the Keplerian to Cartesian conversion $what")
    ecc = primal(ecc) < 0 ? -ecc : ecc
    sma = if primal(ecc) > 1 && primal(sma) > 0
        -sma
    elseif primal(ecc) < 1 && primal(sma) < 0
        -sma
    else
        sma
    end
    abs(1 - primal(ecc)) >= ECC_EPSILON || fail("is parabolic (|ecc - 1| < 1e-11)")
    if primal(ecc) > 1
        w = mod(primal(ta), 2π)
        w = w >= π ? w - 2π : w
        abs(w) <= π - acos(1 / primal(ecc)) || fail("has a true anomaly beyond the asymptotes")
    end
    isfinite(primal(1 + ecc * cos(ta))) || fail("has an infinite radius")
    p = sma * (1 - ecc^2)
    abs(primal(p)) >= eps(Float64) || fail("has a zero semi-parameter")
    radius = p / (1 + ecc * cos(ta))
    sin_aop_ta, cos_aop_ta = sincos(aop + ta)
    sin_inc, cos_inc = sincos(inc)
    sin_raan, cos_raan = sincos(raan)
    sin_aop, cos_aop = sincos(aop)
    x = radius * (cos_aop_ta * cos_raan - cos_inc * sin_aop_ta * sin_raan)
    y = radius * (cos_aop_ta * sin_raan + cos_inc * sin_aop_ta * cos_raan)
    z = radius * sin_aop_ta * sin_inc
    sqrt_gm_p = sqrt(gm / p)
    cos_ta_ecc = cos(ta) + ecc
    sin_ta = sin(ta)
    vx = sqrt_gm_p * cos_ta_ecc * (-sin_aop * cos_raan - cos_inc * sin_raan * cos_aop) -
        sqrt_gm_p * sin_ta * (cos_aop * cos_raan - cos_inc * sin_raan * sin_aop)
    vy = sqrt_gm_p * cos_ta_ecc * (-sin_aop * sin_raan + cos_inc * cos_raan * cos_aop) -
        sqrt_gm_p * sin_ta * (cos_aop * sin_raan + cos_inc * cos_raan * sin_aop)
    vz = sqrt_gm_p * (cos_ta_ecc * sin_inc * cos_aop - sin_ta * sin_inc * sin_aop)
    return SVector{6, T}(x, y, z, vx, vy, vz)
end

# ANISE `try_keplerian_mean_anomaly`: the true anomaly from the mean one, then `try_keplerian`.
function _kep_ma_to_cart(sma, ecc, inc, raan, aop, ma, gm::Real)
    T = promote_type(typeof(sma), typeof(ecc), typeof(inc), typeof(raan), typeof(aop), typeof(ma))
    ta = _ma_to_ta(T(ma), T(ecc))
    return _kep_to_cart(T(sma), T(ecc), T(inc), T(raan), T(aop), ta, gm)
end

# The equinoctial elements `(a, h, k, p, q, lambda)` of Keplerian ones, with the mean anomaly.
function _kep_to_aeq(sma, ecc, inc, raan, aop, ma)
    s = sin(inc / 2)
    sr, cr = sincos(raan)
    sh, ch = sincos(aop + raan)
    return (sma, ecc * sh, ecc * ch, s * sr, s * cr, raan + aop + ma)
end

# ANISE `brouwer_mean_short_to_osculating_kep`: mean Brouwer-Lyddane elements (short
# period only) to osculating Keplerian ones, `(sma, ecc, inc, raan, aop, ma)` in km and rad,
# `re` the equatorial radius used to normalize. It checks only that the orbit is not
# hyperbolic, as ANISE does.
function _brouwer_mean_to_osc(sma, ecc, inc, raan, aop, ma, re::Real, j2::Real)
    T = promote_type(
        typeof(sma), typeof(ecc), typeof(inc), typeof(raan), typeof(aop),
        typeof(ma), Float64
    )
    k2 = j2 / 2
    smap = T(sma) / re
    eccp, incp, raanp, aopp, mean_anom = T(ecc), T(inc), T(raan), T(aop), T(ma)
    if primal(eccp) < 0
        eccp = -eccp
        mean_anom -= π
        aopp += π
    end
    primal(eccp) <= 1 || _brouwer_fail("the mean elements are those of a hyperbola")
    pseudostate = false
    if primal(incp) > deg2rad(175.0)
        incp = π - incp
        raanp = -raanp
        aopp = -aopp
        pseudostate = true
    end
    raanp = _wrap(raanp)
    aopp = _wrap(aopp)
    mean_anom = _wrap(mean_anom)

    eta = sqrt(1 - eccp^2)
    theta = cos(incp)
    p = smap * eta^2
    gm2 = k2 / smap^2
    gm2p = gm2 / eta^4
    tap = _ma_to_ta(mean_anom, eccp)
    rp = p / (1 + eccp * cos(tap))
    adr = smap / rp
    sin_incp = sin(incp)
    theta_sq = theta^2

    sma1 = smap + smap * gm2 * (
        (adr^3 - 1 / eta^3) * (-1 + 3 * theta_sq) +
            3 * (1 - theta_sq) * adr^3 * cos(2 * aopp + 2 * tap)
    )

    decc = eta^2 / 2 * (
        (
            3 * (1 / eta^6) * gm2 * (1 - theta_sq) * cos(2 * aopp + 2 * tap) *
                (3 * eccp * cos(tap)^2 + 3 * cos(tap) + eccp^2 * cos(tap)^3 + eccp)
        ) -
            (gm2p * (1 - theta_sq) * (3 * cos(2 * aopp + tap) + cos(3 * tap + 2 * aopp))) +
            (3 * theta_sq - 1) * gm2 / eta^6 *
            (
            eccp * eta + eccp / (1 + eta) + 3 * eccp * cos(tap)^2 + 3 * cos(tap) +
                eccp^2 * cos(tap)^3
        )
    )

    dinc = gm2p / 2 * theta * sin_incp * (
        3 * cos(2 * aopp + 2 * tap) +
            3 * eccp * cos(2 * aopp + tap) + eccp * cos(2 * aopp + 3 * tap)
    )

    draan = -gm2p / 2 * theta * (
        6 * (tap - mean_anom + eccp * sin(tap)) -
            3 * sin(2 * aopp + 2 * tap) - 3 * eccp * sin(2 * aopp + tap) -
            eccp * sin(2 * aopp + 3 * tap)
    )

    lgh = raanp + aopp + mean_anom +
        gm2p / 4 * (
        6 * (-1 - 2 * theta + 5 * theta_sq) * (tap - mean_anom + eccp * sin(tap)) +
            (3 + 2 * theta - 5 * theta_sq) * (
            3 * sin(2 * aopp + 2 * tap) +
                3 * eccp * sin(2 * aopp + tap) + eccp * sin(2 * aopp + 3 * tap)
        )
    ) +
        gm2p / 4 * eta^2 / (eta + 1) * eccp *
        (
        3 * (1 - theta_sq) * (
            sin(3 * tap + 2 * aopp) * (1 / 3 + adr^2 * eta^2 + adr) +
                sin(2 * aopp + tap) * (1 - adr^2 * eta^2 - adr)
        ) +
            2 * sin(tap) * (3 * theta_sq - 1) * (1 + adr^2 * eta^2 + adr)
    )

    eccpdl = -eta^3 / 4 * gm2p * (
        2 * (-1 + 3 * theta_sq) * (adr^2 * eta^2 + adr + 1) * sin(tap) +
            3 * (1 - theta_sq) * (
            (-adr^2 * eta^2 - adr + 1) * sin(2 * aopp + tap) +
                (adr^2 * eta^2 + adr + 1 / 3) * sin(2 * aopp + 3 * tap)
        )
    )

    ecosl = (eccp + decc) * cos(mean_anom) - eccpdl * sin(mean_anom)
    esinl = (eccp + decc) * sin(mean_anom) + eccpdl * cos(mean_anom)
    ecc1 = sqrt(ecosl^2 + esinl^2)
    ma1 = primal(ecc1) < 1.0e-11 ? zero(T) : _wrap(atan(esinl, ecosl))

    sin_half_i, cos_half_i = sincos(incp / 2)
    sinhalfisinh = (sin_half_i + cos_half_i * dinc / 2) * sin(raanp) + sin_half_i * draan * cos(raanp)
    sinhalficosh = (sin_half_i + cos_half_i * dinc / 2) * cos(raanp) - sin_half_i * draan * sin(raanp)
    sin_half_i_new = sqrt(sinhalfisinh^2 + sinhalficosh^2)
    sin_half_i_new = primal(sin_half_i_new) > 1 ? one(T) : sin_half_i_new
    inc1 = 2 * asin(sin_half_i_new)
    if abs(primal(inc1)) < 1.0e-9 || abs(primal(inc1) - π) < 1.0e-9
        raan1 = zero(T)
    else
        raan1 = _wrap(atan(sinhalfisinh, sinhalficosh))
    end
    aop1 = _wrap(lgh - ma1 - raan1)

    inc_out, raan_out = inc1, raan1
    if pseudostate
        inc_out = π - inc_out
        raan_out = 2 * T(π) - raan_out
    end
    return (sma1 * re, ecc1, inc_out, _wrap(raan_out), _wrap(aop1), _wrap(ma1))
end

# ANISE `calculate_brouwer_mean_short_elements`: the mean elements `(sma, ecc, inc, raan,
# aop, ma)` of the osculating orbit `o`, by the fixed-point iteration of GMAT. Its
# validity errors are `:degeneratecase`; the radii of the center are needed for the
# normalization.
function _brouwer_mean(name::Symbol, o::Orbit{T}) where {T}
    o.center == 399 || _brouwer_fail("the center is body $(o.center), not the Earth (399)")
    osc_sma = _el_sma(o)
    osc_ecc = _el_ecc(o)
    osc_inc = _el_inc(o)
    primal(osc_inc) <= π || _brouwer_fail("the orbit is retrograde beyond 180 degrees")
    primal(_el_periapsis(o)) > 3000 || _brouwer_fail("the periapsis is not above 3000 km")
    0 <= primal(osc_ecc) < 1 || _brouwer_fail("the orbit is not elliptic")
    re = _mean_eq_radius(name, o)

    gm = getfield(o, :gm)
    osc = (osc_sma, osc_ecc, osc_inc, _el_raan(o), _el_aop(o), _el_ma(o))
    pseudostate = false
    cart = getfield(o, :state)
    aeq = (
        _el_sma(o), _el_equinoctial_h(o), _el_equinoctial_k(o), _el_equinoctial_p(o),
        _el_equinoctial_q(o), _el_equinoctial_lambda(o),
    )
    if primal(osc[3]) > deg2rad(175.0)
        pseudostate = true
        osc = (osc[1], osc[2], π - osc[3], -osc[4], -osc[5], osc[6])
        cart = _kep_ma_to_cart(osc..., gm)
        flipped = Orbit{T}(
            cart, getfield(o, :et), getfield(o, :center), getfield(o, :frame),
            gm, getfield(o, :radii)
        )
        aeq = (
            _el_sma(flipped), _el_equinoctial_h(flipped), _el_equinoctial_k(flipped),
            _el_equinoctial_p(flipped), _el_equinoctial_q(flipped), _el_equinoctial_lambda(flipped),
        )
    end

    aeqmean = _kep_to_aeq(osc...)
    kep2 = _brouwer_mean_to_osc(osc..., re, J2_EARTH)
    aeq2 = _kep_to_aeq(kep2...)
    aeqmean2 = _shift(aeqmean, aeq, aeq2)

    emag = T(0.9)
    emag_old = T(1.0)
    ii = 0
    while primal(emag) > 1.0e-8
        sma, ecc, inc, raan, aop, ma = _equinoctial_to_kep(aeqmean2...)
        kep2 = _brouwer_mean_to_osc(sma, ecc, inc, raan, aop, ma, re, J2_EARTH)
        cart2 = _kep_ma_to_cart(kep2..., gm)
        emag = norm(cart - cart2) / norm(cart)
        aeq2 = _kep_to_aeq(kep2...)
        primal(emag_old) > primal(emag) || break      # not improving
        emag_old = emag
        aeqmean2 = _shift(aeqmean2, aeq, aeq2)
        ii > 100 && break
        ii += 1
    end

    sma, ecc, inc, raan, aop, ma = _equinoctial_to_kep(aeqmean2...)
    if primal(ecc) < 0
        ecc = -ecc
        aop += π
        ma -= π
    end
    if pseudostate
        inc = π - inc
        raan = -raan
    end
    return (sma, ecc, inc, _wrap(raan), _wrap(aop), _wrap(ma))
end

# `a .+ (b .- c)` of three tuples, by `map`: Zygote cannot take the cotangent of a broadcast
# over tuples.
_shift(a::Tuple, b::Tuple, c::Tuple) = map((x, y, z) -> x + (y - z), a, b, c)

_el_brouwer_sma(o::Orbit) = _brouwer_mean(:brouwer_sma, o)[1]
_el_brouwer_ecc(o::Orbit) = _brouwer_mean(:brouwer_ecc, o)[2]
_el_brouwer_inc(o::Orbit) = _brouwer_mean(:brouwer_inc, o)[3]
_el_brouwer_raan(o::Orbit) = _brouwer_mean(:brouwer_raan, o)[4]
_el_brouwer_aop(o::Orbit) = _brouwer_mean(:brouwer_aop, o)[5]
_el_brouwer_ma(o::Orbit) = _brouwer_mean(:brouwer_ma, o)[6]
