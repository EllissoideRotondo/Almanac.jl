# The scalar quantities of an `Orbit`. `element` is one `if` chain over the private
# `_el_<name>(o)`; `ELEMENT_NAMES` is the one ordered list of the names, in the order of
# the table of the design (later groups of elements are inserted in that order).

const ELEMENT_NAMES = (:x, :y, :z, :vx, :vy, :vz, :rmag, :vmag,
    :sma, :ecc, :inc, :raan, :aop, :ta, :ma, :ea, :ha,
    :aol, :tlong, :period, :mean_motion, :ta_dot,
    :energy, :c3, :fpa, :semi_parameter, :semi_minor_axis,
    :hmag, :hx, :hy, :hz, :periapsis, :apoapsis,
    :periapsis_altitude, :apoapsis_altitude, :sma_altitude, :altitude,
    :right_ascension, :declination, :velocity_declination,
    :latitude, :longitude, :height, :light_time,
    :bplane_bt, :bplane_br, :bplane_mag, :bplane_angle, :bplane_ltof,
    :equinoctial_a, :equinoctial_h, :equinoctial_k, :equinoctial_p, :equinoctial_q,
    :equinoctial_lambda,
    :brouwer_sma, :brouwer_ecc, :brouwer_inc, :brouwer_raan, :brouwer_aop, :brouwer_ma)

"""
    elements() -> NTuple{N,Symbol}

The names that [`element`](@ref) accepts, in the order of the table in
[Orbits](@ref "Element names"). No CSPICE call has this interface.
"""
elements() = ELEMENT_NAMES

"""
    element(o::Orbit, name::Symbol) -> T

One scalar quantity of the orbit `o`, in km, s and rad and their combinations; the
property `o.name` returns the same. `name` is one of [`elements`](@ref); the table in
[Orbits](@ref "Element names") gives the unit and the name in ANISE of each, and says
how each group behaves. The elements are those of `oscltx` (classical elements, period)
and `recgeo` (geodetic coordinates) where CSPICE has them; no CSPICE call has this
interface, and the other names follow `OrbitalElement` of ANISE.

The Keplerian elements are valid in an inertial frame, the geodetic ones in the
body-fixed frame of the center, and the B-plane ones need an inertial frame; no frame is
checked, and in another frame the numbers are those of the wrong quantity. Where a
quantity is undefined for the orbit (a circular orbit has no true anomaly; a zero
position or velocity has no angular momentum), `element` throws
`SpiceError(:degeneratecase)`. It never returns NaN. Two quantities are not finite or
not meaningful on a parabola: `:sma` is `±Inf` there, and `:period` is 0 for an orbit
that is not elliptic, so that `:ta_dot` is then `Inf`.

Throws `SpiceError(:badarg)` if `name` is not in [`elements`](@ref); also for the vectors
`position`, `velocity`, `hvec` and `evec`, which are properties only, so that a name
known only at run time never gives a union of scalar and vector. The altitudes,
`:latitude`, `:height` and the Brouwer elements throw `SpiceError(:badarg)` too if the
radii of the center are NaN (`o.radii`).

# Examples
```jldoctest
julia> o = Orbit([3.0, 4.0, 0.0, 0.0, 0.0, 2.0], 0.0, 399, 1, 4.0e5);

julia> element(o, :rmag), element(o, :hmag)
(5.0, 10.0)
```
"""
function element(o::Orbit, name::Symbol)
    name === :x && return _el_x(o)
    name === :y && return _el_y(o)
    name === :z && return _el_z(o)
    name === :vx && return _el_vx(o)
    name === :vy && return _el_vy(o)
    name === :vz && return _el_vz(o)
    name === :rmag && return _el_rmag(o)
    name === :vmag && return _el_vmag(o)
    name === :sma && return _el_sma(o)
    name === :ecc && return _el_ecc(o)
    name === :inc && return _el_inc(o)
    name === :raan && return _el_raan(o)
    name === :aop && return _el_aop(o)
    name === :ta && return _el_ta(o)
    name === :ma && return _el_ma(o)
    name === :ea && return _el_ea(o)
    name === :ha && return _el_ha(o)
    name === :aol && return _el_aol(o)
    name === :tlong && return _el_tlong(o)
    name === :period && return _el_period(o)
    name === :mean_motion && return _el_mean_motion(o)
    name === :ta_dot && return _el_ta_dot(o)
    name === :c3 && return _el_c3(o)
    name === :fpa && return _el_fpa(o)
    name === :semi_parameter && return _el_semi_parameter(o)
    name === :semi_minor_axis && return _el_semi_minor_axis(o)
    name === :periapsis && return _el_periapsis(o)
    name === :apoapsis && return _el_apoapsis(o)
    name === :right_ascension && return _el_right_ascension(o)
    name === :declination && return _el_declination(o)
    name === :velocity_declination && return _el_velocity_declination(o)
    name === :periapsis_altitude && return _el_periapsis_altitude(o)
    name === :apoapsis_altitude && return _el_apoapsis_altitude(o)
    name === :sma_altitude && return _el_sma_altitude(o)
    name === :altitude && return _el_altitude(o)
    name === :latitude && return _el_latitude(o)
    name === :longitude && return _el_longitude(o)
    name === :height && return _el_height(o)
    name === :energy && return _el_energy(o)
    name === :hmag && return _el_hmag(o)
    name === :hx && return _el_hx(o)
    name === :hy && return _el_hy(o)
    name === :hz && return _el_hz(o)
    name === :light_time && return _el_light_time(o)
    name === :bplane_bt && return _el_bplane_bt(o)
    name === :bplane_br && return _el_bplane_br(o)
    name === :bplane_mag && return _el_bplane_mag(o)
    name === :bplane_angle && return _el_bplane_angle(o)
    name === :bplane_ltof && return _el_bplane_ltof(o)
    name === :equinoctial_a && return _el_equinoctial_a(o)
    name === :equinoctial_h && return _el_equinoctial_h(o)
    name === :equinoctial_k && return _el_equinoctial_k(o)
    name === :equinoctial_p && return _el_equinoctial_p(o)
    name === :equinoctial_q && return _el_equinoctial_q(o)
    name === :equinoctial_lambda && return _el_equinoctial_lambda(o)
    name === :brouwer_sma && return _el_brouwer_sma(o)
    name === :brouwer_ecc && return _el_brouwer_ecc(o)
    name === :brouwer_inc && return _el_brouwer_inc(o)
    name === :brouwer_raan && return _el_brouwer_raan(o)
    name === :brouwer_aop && return _el_brouwer_aop(o)
    name === :brouwer_ma && return _el_brouwer_ma(o)
    throw(SpiceError(:badarg, "element: unknown name :$name; the names are " *
        join(elements(), ", ")))
end

"Throw the `:degeneratecase` error of `name`, with `what` for the reason."
_degenerate(name::Symbol, what::AbstractString) =
    throw(SpiceError(:degeneratecase,
        (name in ELEMENT_NAMES ? "element :$name" : "$name") * " is undefined: $what"))

_hvec(o::Orbit) = cross(o.position, o.velocity)

# ANISE's thresholds: `ECC_EPSILON` and `INC_EPSILON_DEG` (compared as radians) of
# src/astro/orbit.rs; `f64::EPSILON` guards the zero radius and velocity of `hvec`/`evec`.
const ECC_EPSILON = 1e-11
const INC_EPSILON = deg2rad(1e-6)

# hvec and evec as ANISE's `hvec()` and `evec()`: an error for a zero radius or velocity.
function _hvec_for(name::Symbol, o::Orbit)
    primal(_el_rmag(o)) > eps(Float64) || _degenerate(name, "the radius is zero")
    primal(_el_vmag(o)) > eps(Float64) || _degenerate(name, "the velocity is zero")
    return _hvec(o)
end

function _evec_for(name::Symbol, o::Orbit)
    r, v = o.position, o.velocity
    rm = norm(r)
    primal(rm) > eps(Float64) || _degenerate(name, "the radius is zero")
    return ((norm(v)^2 - getfield(o, :gm) / rm) * r - dot(r, v) * v) / getfield(o, :gm)
end

# The node vector `ẑ × h` and its two uses.
_node(h) = SVector(-h[2], h[1], zero(eltype(h)))

_el_x(o::Orbit) = getfield(o, :state)[1]
_el_y(o::Orbit) = getfield(o, :state)[2]
_el_z(o::Orbit) = getfield(o, :state)[3]
_el_vx(o::Orbit) = getfield(o, :state)[4]
_el_vy(o::Orbit) = getfield(o, :state)[5]
_el_vz(o::Orbit) = getfield(o, :state)[6]
_el_rmag(o::Orbit) = norm(o.position)
_el_vmag(o::Orbit) = norm(o.velocity)
# ANISE: energy_km2_s2 errors for a zero radius; hmag, hx, hy, hz go through `hvec()`,
# which errors for a zero radius or velocity.
function _el_energy(o::Orbit)
    primal(_el_rmag(o)) > eps(Float64) || _degenerate(:energy, "the radius is zero")
    return _el_vmag(o)^2 / 2 - getfield(o, :gm) / _el_rmag(o)
end
_el_hmag(o::Orbit) = norm(_hvec_for(:hmag, o))
_el_hx(o::Orbit) = _hvec_for(:hx, o)[1]
_el_hy(o::Orbit) = _hvec_for(:hy, o)[2]
_el_hz(o::Orbit) = _hvec_for(:hz, o)[3]
_el_light_time(o::Orbit) = _el_rmag(o) / clight()

# Keplerian elements: ANISE src/astro/orbit.rs (and orbit_hyperbolic.rs for :ha), angles in rad.

_el_ecc(o::Orbit) = norm(_evec_for(:ecc, o))
_el_sma(o::Orbit) = -getfield(o, :gm) / (2 * _el_energy(o))

function _el_inc(o::Orbit)
    h = _hvec_for(:inc, o)
    # ANISE: NaN for a rectilinear orbit (h = 0, 0/0), and when rounding leaves [-1, 1]
    return _nonan(:inc, acos(clamp(h[3] / norm(h), -1, 1)))
end

# The arc cosine of ANISE's `aop_deg` and `raan_deg`: where it has no value (NaN or
# |c| > 1), 180 degrees if c > 1 and 0 otherwise, so a NaN or c < -1 gives 0.
function _acos_or_fallback(c::T) where {T}
    (isnan(c) || abs(primal(c)) > 1) && return primal(c) > 1 ? T(π) : zero(T)
    return acos(c)
end

function _el_aop(o::Orbit{T}) where {T}
    h = _hvec_for(:aop, o)
    e = _evec_for(:aop, o)
    n = _node(h)
    c = dot(n, e) / (norm(n) * norm(e))
    (isnan(c) || abs(primal(c)) > 1) && return _acos_or_fallback(c)
    aop = acos(c)
    return primal(e[3]) < 0 ? 2 * T(π) - aop : aop
end

function _el_raan(o::Orbit{T}) where {T}
    n = _node(_hvec_for(:raan, o))
    c = n[1] / norm(n)
    (isnan(c) || abs(primal(c)) > 1) && return _acos_or_fallback(c)
    raan = acos(c)
    return primal(n[2]) < 0 ? 2 * T(π) - raan : raan
end

function _el_ta(o::Orbit{T}) where {T}
    ecc = _el_ecc_for(:ta, o)
    e = _evec_for(:ta, o)
    ta = acos(clamp(dot(e, o.position) / (ecc * _el_rmag(o)), -1, 1))
    return primal(dot(o.position, o.velocity)) < 0 ? 2 * T(π) - ta : ta
end

# The eccentricity, and the circular-orbit error of the anomalies built on the true one.
function _el_ecc_for(name::Symbol, o::Orbit)
    ecc = norm(_evec_for(name, o))
    primal(ecc) < ECC_EPSILON && _degenerate(name, "the orbit is circular (ecc = $(ecc) < 1e-11), " *
        "the line of apsides is undefined")
    return ecc
end

_wrap(x) = mod(x, 2 * oftype(x, π))

# ANISE returns NaN where a parabola or a zero velocity makes a quotient 0/0 or Inf*0;
# Almanac throws instead, so that no NaN escapes from `element`.
_nonan(name::Symbol, x) = isnan(x) ? _degenerate(name,
    "the formula is 0/0 or Inf*0 here (parabolic orbit or zero velocity); ANISE returns NaN") : x

function _el_ea(o::Orbit)
    ecc = _el_ecc_for(:ea, o)
    primal(ecc) > 1 && _degenerate(:ea, "the orbit is hyperbolic; ANISE returns NaN here")
    s, c = sincos(_el_ta(o))
    d = 1 + ecc * c
    return _nonan(:ea, atan(sqrt(1 - ecc^2) * s / d, (ecc + c) / d))
end

function _el_ma(o::Orbit)
    ecc = _el_ecc_for(:ma, o)
    abs(primal(ecc) - 1) < ECC_EPSILON && _degenerate(:ma, "the orbit is parabolic (|ecc - 1| < 1e-11)")
    if primal(ecc) < 1
        ea = _el_ea(o)
        return _wrap(ea - ecc * sin(ea))
    end
    s, c = sincos(_el_ta(o))
    sinh_h = sqrt(ecc^2 - 1) * s / (1 + ecc * c)
    return ecc * sinh_h - asinh(sinh_h)
end

function _el_ha(o::Orbit)
    ecc = _el_ecc(o)
    primal(ecc) <= 1 && _degenerate(:ha, "the orbit is not hyperbolic (ecc = $(ecc) <= 1)")
    s, c = sincos(_el_ta(o))
    return _nonan(:ha, _wrap(asinh(s * sqrt(ecc^2 - 1) / (1 + ecc * c))))
end

function _el_aol(o::Orbit)
    h = _hvec_for(:aol, o)
    h /= norm(h)
    n = _node(h)
    nn = norm(n)
    primal(nn) < INC_EPSILON && _degenerate(:aol, "the orbit is equatorial (|ẑ×ĥ| = $(nn) < 1e-6 deg)")
    n /= nn
    m = cross(h, n)
    r = o.position / _el_rmag(o)
    return _nonan(:aol, _wrap(atan(dot(r, m), dot(r, n))))   # NaN: rectilinear orbit, h = 0
end

_el_tlong(o::Orbit) = _wrap(_el_aop(o) + _el_raan(o) + _el_ta(o))

function _el_period(o::Orbit{T}) where {T}
    p2 = _el_sma(o)^3 / getfield(o, :gm)
    return primal(p2) > 0 && isfinite(primal(p2)) ? 2 * T(π) * sqrt(p2) : zero(T)
end

_el_mean_motion(o::Orbit) = sqrt(getfield(o, :gm) / abs(_el_sma(o))^3)
_el_ta_dot(o::Orbit{T}) where {T} = 2 * T(π) / _el_period(o)
_el_c3(o::Orbit) = -getfield(o, :gm) / _el_sma(o)

function _el_fpa(o::Orbit)
    ecc = _el_ecc(o)
    s, c = sincos(_el_ta(o))
    d = sqrt(1 + 2 * ecc * c + ecc^2)
    return _nonan(:fpa, atan(ecc * s / d, (1 + ecc * c) / d))
end

_el_semi_parameter(o::Orbit) = _nonan(:semi_parameter, _el_sma(o) * (1 - _el_ecc(o)^2))

function _el_semi_minor_axis(o::Orbit)
    ecc = _el_ecc(o)
    primal(ecc) <= 1 && return _nonan(:semi_minor_axis, _el_sma(o) * sqrt(1 - ecc^2))
    h = _hvec_for(:semi_minor_axis, o)
    return dot(h, h) / (getfield(o, :gm) * sqrt(ecc^2 - 1))
end

_el_periapsis(o::Orbit) = _nonan(:periapsis, _el_sma(o) * (1 - _el_ecc(o)))
_el_apoapsis(o::Orbit) = _nonan(:apoapsis, _el_sma(o) * (1 + _el_ecc(o)))

# Right ascension: atan2 of the position, 0 for a zero radius, so it does not throw.
function _el_right_ascension(o::Orbit)
    s = getfield(o, :state)
    return _wrap(atan(s[2], s[1]))
end

function _el_declination(o::Orbit)
    primal(_el_rmag(o)) > eps(Float64) || _degenerate(:declination, "the radius is zero")
    return asin(clamp(getfield(o, :state)[3] / _el_rmag(o), -1, 1))
end

function _el_velocity_declination(o::Orbit)
    primal(_el_vmag(o)) > eps(Float64) || _degenerate(:velocity_declination, "the velocity is zero")
    return asin(clamp(getfield(o, :state)[6] / _el_vmag(o), -1, 1))
end

"""
    keplerian(o::Orbit) -> NamedTuple{(:sma, :ecc, :inc, :raan, :aop, :ta)}

The classical elements of `o`: semi-major axis (km, negative for a hyperbola),
eccentricity, inclination, right ascension of the ascending node, argument of periapsis
and true anomaly (rad). `oscltx` returns these elements of a state with more; this call
matches `Orbit.to_keplerian_vec` of ANISE, with angles in rad.

The values are those of the properties of the same names (see [`element`](@ref)), so the
degenerate cases are the same: `:ta` throws `SpiceError(:degeneratecase)` for a circular
orbit (`ecc < 1e-11`), while `:raan` and `:aop` have their ANISE fallback values for an
equatorial orbit.

# Examples
```jldoctest
julia> o = Orbit([7000.0, 0.0, 0.0, 0.0, 8.0, 0.0], 0.0, 399, 1, 398600.4415);

julia> k = keplerian(o); (round(k.sma), round(k.ecc; digits=3), k.inc, k.ta)
(7990.0, 0.124, 0.0, 0.0)
```
"""
keplerian(o::Orbit) = (; sma=_el_sma(o), ecc=_el_ecc(o), inc=_el_inc(o), raan=_el_raan(o),
    aop=_el_aop(o), ta=_el_ta(o))

# Altitudes and geodetic coordinates: src/astro/orbit_geodetic.rs of ANISE. The radii of the
# center come from the pool; where it has none they are NaN and these elements are errors.
function _radii_missing(name::Symbol, o::Orbit)
    throw(SpiceError(:badarg, "element :$name needs the radii of body $(o.center), " *
        "which the pool does not hold"))
end

# ANISE `Frame::mean_equatorial_radius_km`: the mean of the two equatorial radii.
function _mean_eq_radius(name::Symbol, o::Orbit)
    r = getfield(o, :radii)
    (isnan(r[1]) || isnan(r[2])) && _radii_missing(name, o)
    return (r[1] + r[2]) / 2
end

_polar_radius(name::Symbol, o::Orbit) =
    isnan(getfield(o, :radii)[3]) ? _radii_missing(name, o) : getfield(o, :radii)[3]

_el_altitude(o::Orbit) = _el_rmag(o) - _mean_eq_radius(:altitude, o)
_el_sma_altitude(o::Orbit) = _el_sma(o) - _mean_eq_radius(:sma_altitude, o)
_el_periapsis_altitude(o::Orbit) = _el_periapsis(o) - _mean_eq_radius(:periapsis_altitude, o)
_el_apoapsis_altitude(o::Orbit) = _el_apoapsis(o) - _mean_eq_radius(:apoapsis_altitude, o)

# ANISE `latlongalt`: the procedure of Heikkinen, which is not iterative, on the ellipsoid
# of semi-axes (a, a, b) with a the mean equatorial radius. Returns the geodetic latitude
# (rad, ANISE: degrees) and the height (km). A NaN (the center of the body) is an error.
function _latlongalt(name::Symbol, o::Orbit)
    a_km = _mean_eq_radius(name, o)
    b_km = _polar_radius(name, o)
    s = getfield(o, :state)
    x, y, z = s[1], s[2], s[3]
    e2 = (a_km^2 - b_km^2) / a_km^2
    e_prime2 = (a_km^2 - b_km^2) / b_km^2
    p = sqrt(x^2 + y^2)
    big_f = 54 * b_km^2 * z^2
    big_g = p^2 + (1 - e2) * z^2 - e2 * (a_km^2 - b_km^2)
    c = (e2^2 * big_f * p^2) / big_g^3
    sc = (1 + c + sqrt(c^2 + 2 * c))^(1.0 / 3.0)
    k = sc + 1 + 1 / sc
    big_p = big_f / (3 * k^2 * big_g^2)
    big_q = sqrt(1 + 2 * e2^2 * big_p)
    r0 = (-big_p * e2 * p) / (1 + big_q) +
        sqrt(0.5 * a_km^2 * (1 + 1 / big_q) - (big_p * (1 - e2) * z^2) / (big_q * (1 + big_q)) -
            0.5 * big_p * p^2)
    big_u = sqrt((p - e2 * r0)^2 + z^2)
    big_v = sqrt((p - e2 * r0)^2 + (1 - e2) * z^2)
    z0 = b_km^2 * z / (a_km * big_v)
    alt = big_u * (1 - b_km^2 / (a_km * big_v))
    lat = atan((z + e_prime2 * z0) / p)
    return _nonan(name, lat), _nonan(name, alt)
end

_el_latitude(o::Orbit) = _latlongalt(:latitude, o)[1]
_el_height(o::Orbit) = _latlongalt(:height, o)[2]

# ANISE `longitude_deg`: the longitude in [-180, 180) degrees, 0 for a NaN position; the
# state is in the body-fixed frame.
function _el_longitude(o::Orbit{T}) where {T}
    s = getfield(o, :state)
    (isnan(s[1]) || isnan(s[2])) && return zero(T)
    a = mod(atan(s[2], s[1]), 2 * T(π))
    return primal(a) >= Float64(π) ? a - 2 * T(π) : a
end

# Equinoctial elements: src/astro/orbit_equinoctial.rs of ANISE (rad here). The mean
# longitude is raan + aop + ma and is not wrapped.
_el_equinoctial_a(o::Orbit) = _el_sma(o)
_el_equinoctial_h(o::Orbit) = _el_ecc(o) * sin(_el_aop(o) + _el_raan(o))
_el_equinoctial_k(o::Orbit) = _el_ecc(o) * cos(_el_aop(o) + _el_raan(o))
_el_equinoctial_p(o::Orbit) = sin(_el_inc(o) / 2) * sin(_el_raan(o))
_el_equinoctial_q(o::Orbit) = sin(_el_inc(o) / 2) * cos(_el_raan(o))
_el_equinoctial_lambda(o::Orbit) = _el_raan(o) + _el_aop(o) + _el_ma(o)

# ANISE `equinoctial_to_keplerian`: (sma, ecc, inc, raan, aop, ma), the angles not wrapped.
function _equinoctial_to_kep(sma, h, k, p, q, lambda)
    ecc = sqrt(h * h + k * k)
    s_sq = p * p + q * q                       # sin^2(i/2)
    inc = primal(s_sq) <= 1 ? acos(1 - 2 * s_sq) : oftype(s_sq, π)   # rounding beyond 1: ANISE gives acos(-1)
    raan = atan(p, q)
    aop_plus_raan = atan(h, k)
    return (sma, ecc, inc, raan, aop_plus_raan - raan, lambda - aop_plus_raan)
end

# B-plane: ANISE `BPlane::new` (src/astro/orbit_hyperbolic.rs). The triad is built on the
# +z axis of the frame of the orbit. Returns (b_t, b_r, s_hat) of the hyperbola.
function _bplane(name::Symbol, o::Orbit)
    ecc = _el_ecc(o)
    primal(ecc) > 1 || _degenerate(name, "the orbit is not hyperbolic (ecc = $(ecc) <= 1)")
    e_hat = _evec_for(name, o) / ecc
    h = _hvec_for(name, o)
    n_hat = cross(h / norm(h), e_hat)
    cos_nu = 1 / ecc
    sin_nu = sqrt(1 - cos_nu^2)
    s = cos_nu * e_hat + sin_nu * n_hat
    s_hat = s / norm(s)
    b_vec = _el_semi_minor_axis(o) * (sin_nu * e_hat - cos_nu * n_hat)
    t = SVector(s_hat[2], -s_hat[1], zero(eltype(s_hat)))
    t_norm = norm(t)
    primal(t_norm) > 1e-12 || _degenerate(name, "the incoming asymptote is parallel to the z axis " *
        "of the frame (|s x z| = $(t_norm) <= 1e-12)")
    t_hat = t / t_norm
    r_hat = cross(s_hat, t_hat)
    return dot(b_vec, t_hat), dot(b_vec, r_hat), s_hat
end

_el_bplane_bt(o::Orbit) = _bplane(:bplane_bt, o)[1]
_el_bplane_br(o::Orbit) = _bplane(:bplane_br, o)[2]
function _el_bplane_mag(o::Orbit)
    bt, br, _ = _bplane(:bplane_mag, o)
    return sqrt(bt^2 + br^2)
end
# ANISE `between_pm_180` of the arc tangent: [-π, π), so an angle of exactly π is -π.
function _el_bplane_angle(o::Orbit{T}) where {T}
    bt, br, _ = _bplane(:bplane_angle, o)
    a = atan(br, bt)
    return primal(a) >= Float64(π) ? a - 2 * T(π) : a
end
function _el_bplane_ltof(o::Orbit)
    _, _, s_hat = _bplane(:bplane_ltof, o)
    c3 = _el_c3(o)
    primal(c3) > 0 || _degenerate(:bplane_ltof, "c3 = $(c3) <= 0, the excess speed is undefined")
    return -dot(o.position, s_hat) / sqrt(c3)
end

# ANISE `vinf_periapsis_km` and `vinf_turn_angle_deg` (src/astro/orbit.rs, c3 form).
function _vinf_hyperbolic(name::Symbol, o::Orbit)
    ecc = _el_ecc(o)
    primal(ecc) > 1 || _degenerate(name, "the orbit is not hyperbolic (ecc = $(ecc) <= 1)")
    return _el_c3(o)
end

"""
    vinf_periapsis(o::Orbit, turn_angle::Real)

Return the periapsis radius in km (a number; a `Dual` where `o` or `turn_angle` is one)
of a flyby of the center of `o` that turns the excess velocity by `turn_angle` rad, for
the excess speed of the hyperbola `o`:
`(1 / cos((π - turn_angle) / 2) - 1) μ / C3`, where `C3 = o.c3` (km²/s²) is the square
of the excess speed and `μ = o.gm`. No CSPICE call has this interface; it matches
`Orbit.vinf_periapsis_km` of ANISE, with the angle in rad. It is the inverse of
[`vinf_turn_angle`](@ref).

`turn_angle` is in [0, π]: π gives 0 and 0 gives a very large radius; outside [0, π] the
formula is evaluated without a check. Throws `SpiceError(:degeneratecase)` if the orbit
is not hyperbolic (`ecc <= 1`, ANISE: `NotHyperbolic`) or the result is NaN.

# Examples
```jldoctest
julia> o = Orbit([7000.0, 0.0, 0.0, 0.0, 12.0, 0.0], 0.0, 399, 1, 398600.4415);

julia> rp = vinf_periapsis(o, vinf_turn_angle(o, 8000.0)); round(rp; digits=6)
8000.0
```
"""
function vinf_periapsis(o::Orbit, turn_angle::Real)
    c3 = _vinf_hyperbolic(:vinf_periapsis, o)
    cos_rho = cos((π - turn_angle) / 2)
    return _nonan(:vinf_periapsis, (1 / cos_rho - 1) * getfield(o, :gm) / c3)
end

"""
    vinf_turn_angle(o::Orbit, periapsis::Real)

Return the angle in rad (a number; a `Dual` where `o` or `periapsis` is one) by which a
flyby of the center of `o` with periapsis radius `periapsis` km turns the
excess velocity of the hyperbola `o`:
`π - 2 acos(1 / (1 + C3 periapsis / μ))`, where `C3 = o.c3` (km²/s²) is the square of the
excess speed and `μ = o.gm`. No CSPICE call has this interface; it matches
`Orbit.vinf_turn_angle_deg` of ANISE, with the angle in rad. It depends on the excess
speed only, not on where along the trajectory `o` is. It is the inverse of
[`vinf_periapsis`](@ref).

The result is in [0, π] for `periapsis >= 0` (π for 0), and is wrapped into [0, 2π) for a
negative one. Throws `SpiceError(:degeneratecase)` if the orbit is not hyperbolic
(`ecc <= 1`, ANISE: `NotHyperbolic`) or the arc cosine has no value (ANISE: NaN), as for
a negative `periapsis` below `-μ / C3`.

# Examples
```jldoctest
julia> o = Orbit([7000.0, 0.0, 0.0, 0.0, 12.0, 0.0], 0.0, 399, 1, 398600.4415);

julia> round(rad2deg(vinf_turn_angle(o, 8000.0)); digits=3)
77.113
```
"""
function vinf_turn_angle(o::Orbit, periapsis::Real)
    c3 = _vinf_hyperbolic(:vinf_turn_angle, o)
    x = 1 / (1 + c3 * (periapsis / getfield(o, :gm)))
    (isnan(x) || abs(primal(x)) > 1) && _degenerate(:vinf_turn_angle, "the arc cosine of $(x) has no value; ANISE returns NaN")
    return _wrap(π - 2 * acos(x))
end
