# The scalar quantities of an `Orbit`. `element` is one `if` chain over the private
# `_el_<name>(o)`; `ELEMENT_NAMES` is the one ordered list of the names, in the order of
# the table of the design (later groups of elements are inserted in that order).

const ELEMENT_NAMES = (:x, :y, :z, :vx, :vy, :vz, :rmag, :vmag,
    :sma, :ecc, :inc, :raan, :aop, :ta, :ma, :ea, :ha,
    :aol, :tlong, :period, :mean_motion, :ta_dot,
    :energy, :c3, :fpa, :semi_parameter, :semi_minor_axis,
    :hmag, :hx, :hy, :hz, :periapsis, :apoapsis,
    :right_ascension, :declination, :velocity_declination, :light_time)

"""
    elements() -> NTuple{N,Symbol}

The names that [`element`](@ref) accepts, in the order of the manual. This call has no
CSPICE counterpart.
"""
elements() = ELEMENT_NAMES

"""
    element(o::Orbit, name::Symbol) -> T

One scalar quantity of the orbit `o`, in km, s and rad and their combinations; the
property `o.name` returns the same. This call has no CSPICE counterpart; the names are
those of `StateElement` of ANISE (`element(o, :rmag)` is `rmag_km`, `:hmag` is `hmag`,
`:energy` is `energy_km2_s2`, `:light_time` is `light_time`). The light time is the
distance divided by [`clight`](@ref), in s.

The Keplerian elements (`:sma`, `:ecc`, `:inc`, `:raan`, `:aop`, `:ta`, `:ma`, `:ea`,
`:ha`, `:aol`, `:tlong`, `:period`, `:mean_motion`, `:ta_dot`, `:c3`, `:fpa`,
`:semi_parameter`, `:semi_minor_axis`, `:periapsis`, `:apoapsis`, `:right_ascension`,
`:declination`, `:velocity_declination`) follow the functions of `Orbit` of ANISE
(`sma_km`, `ecc`, `inc_deg`, ...) with angles in rad (the angles `:inc`, `:ta`, `:aop`,
`:raan`, ... are in [0, 2π) as in ANISE, `:ea` in (-π, π], `:declination` in [-π/2, π/2])
and `:mean_motion` and `:ta_dot` in rad/s. Where ANISE returns an error for an undefined
quantity, `element` throws `SpiceError(:degeneratecase)`: `:ta` and everything built on
it (`:ea`, `:ma`, `:fpa`, `:tlong`, `:ha`) for `ecc < 1e-11`, `:aol` for an orbit within
`1e-6` degrees of the equator, `:ma` for `|ecc - 1| < 1e-11`, `:ea` for `ecc > 1`, `:ha`
for `ecc <= 1`, and every element that needs the angular momentum or the eccentricity
vector for a zero position or velocity. Where ANISE substitutes a value, so does
`element`. Where ANISE returns NaN (a parabola, whose `:sma` is infinite, or a zero velocity:
`:ea`, `:fpa`, `:semi_parameter`, `:semi_minor_axis`, `:periapsis`, `:apoapsis`), `element`
throws `SpiceError(:degeneratecase)`, so no NaN is returned. `:raan` and `:aop` give 0 where the arc cosine has no value (an equatorial
orbit has `:raan == 0`), π where its argument rounds above 1. `:period` is 0 for an orbit
that is not elliptic, so `:ta_dot` is `Inf` there.

The vectors `position`, `velocity`, `hvec` and `evec` are properties only: `element` throws for
them, so that a name known only at run time never gives a union of scalar and vector.

Throws `SpiceError(:badarg)` if `name` is not in [`elements`](@ref).

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
    name === :energy && return _el_energy(o)
    name === :hmag && return _el_hmag(o)
    name === :hx && return _el_hx(o)
    name === :hy && return _el_hy(o)
    name === :hz && return _el_hz(o)
    name === :light_time && return _el_light_time(o)
    throw(SpiceError(:badarg, "element: unknown name :$name; the names are " *
        join(elements(), ", ")))
end

"Throw the `:degeneratecase` error of `name`, with `what` for the reason."
_degenerate(name::Symbol, what::AbstractString) =
    throw(SpiceError(:degeneratecase, "element :$name is undefined: $what"))

_hvec(o::Orbit) = cross(o.position, o.velocity)

# ANISE's thresholds: `ECC_EPSILON` and `INC_EPSILON_DEG` (compared as radians) of
# src/astro/orbit.rs; `f64::EPSILON` guards the zero radius and velocity of `hvec`/`evec`.
const ECC_EPSILON = 1e-11
const INC_EPSILON = deg2rad(1e-6)

# hvec and evec as ANISE's `hvec()` and `evec()`: an error for a zero radius or velocity.
function _hvec_for(name::Symbol, o::Orbit)
    _el_rmag(o) > eps(Float64) || _degenerate(name, "the radius is zero")
    _el_vmag(o) > eps(Float64) || _degenerate(name, "the velocity is zero")
    return _hvec(o)
end

function _evec_for(name::Symbol, o::Orbit)
    r, v = o.position, o.velocity
    rm = norm(r)
    rm > eps(Float64) || _degenerate(name, "the radius is zero")
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
_el_energy(o::Orbit) = _el_vmag(o)^2 / 2 - getfield(o, :gm) / _el_rmag(o)
_el_hmag(o::Orbit) = norm(_hvec(o))
_el_hx(o::Orbit) = _hvec(o)[1]
_el_hy(o::Orbit) = _hvec(o)[2]
_el_hz(o::Orbit) = _hvec(o)[3]
_el_light_time(o::Orbit) = _el_rmag(o) / clight()

# Keplerian elements: ANISE src/astro/orbit.rs (and orbit_hyperbolic.rs for :ha), angles in rad.

_el_ecc(o::Orbit) = norm(_evec_for(:ecc, o))
_el_sma(o::Orbit) = -getfield(o, :gm) / (2 * _el_energy_for(:sma, o))
_el_energy_for(name::Symbol, o::Orbit) = (_el_rmag(o) > eps(Float64) ||
    _degenerate(name, "the radius is zero"); _el_energy(o))

function _el_inc(o::Orbit)
    h = _hvec_for(:inc, o)
    return acos(clamp(h[3] / norm(h), -1, 1))     # ANISE: NaN when rounding leaves [-1, 1]
end

# The arc cosine of ANISE's `aop_deg` and `raan_deg`: where it has no value (NaN or
# |c| > 1), 180 degrees if c > 1 and 0 otherwise, so a NaN or c < -1 gives 0.
function _acos_or_fallback(c::T) where {T}
    (isnan(c) || abs(c) > 1) && return c > 1 ? T(π) : zero(T)
    return acos(c)
end

function _el_aop(o::Orbit{T}) where {T}
    h = _hvec_for(:aop, o)
    e = _evec_for(:aop, o)
    n = _node(h)
    c = dot(n, e) / (norm(n) * norm(e))
    (isnan(c) || abs(c) > 1) && return _acos_or_fallback(c)
    aop = acos(c)
    return e[3] < 0 ? 2 * T(π) - aop : aop
end

function _el_raan(o::Orbit{T}) where {T}
    n = _node(_hvec_for(:raan, o))
    c = n[1] / norm(n)
    (isnan(c) || abs(c) > 1) && return _acos_or_fallback(c)
    raan = acos(c)
    return n[2] < 0 ? 2 * T(π) - raan : raan
end

function _el_ta(o::Orbit{T}) where {T}
    ecc = _el_ecc_for(:ta, o)
    e = _evec_for(:ta, o)
    ta = acos(clamp(dot(e, o.position) / (ecc * _el_rmag(o)), -1, 1))
    return dot(o.position, o.velocity) < 0 ? 2 * T(π) - ta : ta
end

# The eccentricity, and the circular-orbit error of the anomalies built on the true one.
function _el_ecc_for(name::Symbol, o::Orbit)
    ecc = norm(_evec_for(name, o))
    ecc < ECC_EPSILON && _degenerate(name, "the orbit is circular (ecc = $(ecc) < 1e-11), " *
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
    ecc > 1 && _degenerate(:ea, "the orbit is hyperbolic; ANISE returns NaN here")
    s, c = sincos(_el_ta(o))
    d = 1 + ecc * c
    return _nonan(:ea, atan(sqrt(1 - ecc^2) * s / d, (ecc + c) / d))
end

function _el_ma(o::Orbit)
    ecc = _el_ecc_for(:ma, o)
    abs(ecc - 1) < ECC_EPSILON && _degenerate(:ma, "the orbit is parabolic (|ecc - 1| < 1e-11)")
    if ecc < 1
        ea = _el_ea(o)
        return _wrap(ea - ecc * sin(ea))
    end
    s, c = sincos(_el_ta(o))
    sinh_h = sqrt(ecc^2 - 1) * s / (1 + ecc * c)
    return ecc * sinh_h - asinh(sinh_h)
end

function _el_ha(o::Orbit)
    ecc = _el_ecc(o)
    ecc <= 1 && _degenerate(:ha, "the orbit is not hyperbolic (ecc = $(ecc) <= 1)")
    s, c = sincos(_el_ta(o))
    return _wrap(asinh(s * sqrt(ecc^2 - 1) / (1 + ecc * c)))
end

function _el_aol(o::Orbit)
    h = _hvec_for(:aol, o)
    h /= norm(h)
    n = _node(h)
    nn = norm(n)
    nn < INC_EPSILON && _degenerate(:aol, "the orbit is equatorial (|ẑ×ĥ| = $(nn) < 1e-6 deg)")
    n /= nn
    m = cross(h, n)
    r = o.position / _el_rmag(o)
    return _wrap(atan(dot(r, m), dot(r, n)))
end

_el_tlong(o::Orbit) = _wrap(_el_aop(o) + _el_raan(o) + _el_ta(o))

function _el_period(o::Orbit{T}) where {T}
    p2 = _el_sma(o)^3 / getfield(o, :gm)
    return p2 > 0 && isfinite(p2) ? 2 * T(π) * sqrt(p2) : zero(T)
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
    ecc <= 1 && return _nonan(:semi_minor_axis, _el_sma(o) * sqrt(1 - ecc^2))
    h = _hvec_for(:semi_minor_axis, o)
    return dot(h, h) / (getfield(o, :gm) * sqrt(ecc^2 - 1))
end

_el_periapsis(o::Orbit) = _nonan(:periapsis, _el_sma(o) * (1 - _el_ecc(o)))
_el_apoapsis(o::Orbit) = _nonan(:apoapsis, _el_sma(o) * (1 + _el_ecc(o)))

# ANISE returns NaN for a zero radius or velocity here; Almanac throws, so no NaN escapes.
function _el_right_ascension(o::Orbit)
    s = getfield(o, :state)
    return _wrap(atan(s[2], s[1]))
end

function _el_declination(o::Orbit)
    _el_rmag(o) > eps(Float64) || _degenerate(:declination, "the radius is zero")
    return asin(clamp(getfield(o, :state)[3] / _el_rmag(o), -1, 1))
end

function _el_velocity_declination(o::Orbit)
    _el_vmag(o) > eps(Float64) || _degenerate(:velocity_declination, "the velocity is zero")
    return asin(clamp(getfield(o, :state)[6] / _el_vmag(o), -1, 1))
end

"""
    keplerian(o::Orbit) -> NamedTuple{(:sma, :ecc, :inc, :raan, :aop, :ta)}

The classical elements of `o`: semi-major axis (km, negative for a hyperbola),
eccentricity, inclination, right ascension of the ascending node, argument of periapsis
and true anomaly (rad). This call has no CSPICE counterpart; it matches
`Orbit.to_keplerian_vec` of ANISE, with angles in rad.

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
