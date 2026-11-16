# Orbits

An orbit is a state about a central body. No CSPICE call has these interfaces; a docstring names the CSPICE routine that computes the same quantity (`oscltx`, `recgeo`, `conics`, `prop2b`) where there is one, and the function of ANISE that it matches. The deviations from ANISE are in [Beyond CSPICE](@ref "Beyond CSPICE"). Lengths are in km, times in seconds, angles in rad. [Orbits and relative motion](@ref) shows them in use.

## Orbit objects

```@docs
Orbit
orbit
```

## Elements

```@docs
element
elements
keplerian
vinf_periapsis
vinf_turn_angle
```

## Element names

The table lists every name of [`elements`](@ref), in the order that `elements()` returns, with its unit and the name in ANISE. The groups behave as follows.

- **Position and velocity** (`:x` to `:vmag`) are the state. The vectors `position`, `velocity`, `hvec` and `evec` are properties only.
- **Keplerian** (`:sma` to `:velocity_declination` in the order of `elements()`, which takes in `:energy`, `:hmag` to `:hz`, the apsides and the altitudes): angles in rad. `:inc` is in [0, π]; `:ta`, `:aop`, `:raan`, `:aol`, `:tlong`, `:ha` and `:right_ascension` are in [0, 2π], where `:raan` and `:aop` can be exactly 2π; `:ea` is in (-π, π]; the elliptic `:ma` is in [0, 2π] and the hyperbolic `:ma` is unbounded; `:declination` is in [-π/2, π/2]. `:raan` and `:aop` are 0 where the arc cosine has no value (an equatorial orbit has `:raan == 0`) and π where its argument rounds above 1, as in ANISE. `:mean_motion` is `sqrt(gm / |sma|³)`; for an orbit that is not elliptic `:period` is 0 and `:ta_dot` is `Inf`.
- **Undefined quantities** throw `SpiceError(:degeneratecase)`: `:ta` and what is built on it (`:ea`, `:ma`, `:fpa`, `:tlong`, `:ha`) for `ecc < 1e-11`; `:aol` for an orbit within 1e-6° of the equator; `:ma` for `|ecc - 1| < 1e-11`; `:ea` for `ecc > 1`; `:ha` for `ecc <= 1`; every quantity that needs the angular momentum (a zero position or velocity) and every one that needs the eccentricity vector or `:energy` (a zero position; a zero velocity gives `ecc = 1`).
- **B-plane** (`:bplane_bt` to `:bplane_ltof`): for a hyperbola in an inertial frame; the triad is built on the `+z` axis of the frame. `:bplane_angle` is in [-π, π). They throw for `ecc <= 1`, for an incoming asymptote within 1e-12 of the `z` axis, and `:bplane_ltof` also for `c3 <= 0`. See also [`vinf_periapsis`](@ref) and [`vinf_turn_angle`](@ref).
- **Altitudes** are the radius, `:sma` or an apsis radius minus the mean equatorial radius `(radii[1] + radii[2]) / 2` of the center. **Geodetic** `:latitude` (rad, in [-π/2, π/2]) and `:height` (km) are on the spheroid of semi-axes `a`, `a`, `b` with that `a` and `b = radii[3]`, in the body-fixed frame of the center; `:longitude` is `atan(y, x)` in [-π, π). They need the radii of the center.
- **Equinoctial**: `:equinoctial_h` and `:equinoctial_k` are `ecc * sin(aop + raan)` and `ecc * cos(aop + raan)`, `:equinoctial_p` and `:equinoctial_q` are `sin(inc / 2) * sin(raan)` and `sin(inc / 2) * cos(raan)`, and `:equinoctial_lambda` is `raan + aop + ma`, not wrapped. They throw where the Keplerian elements they use throw.
- **Brouwer mean elements** (short-period terms only, from the iteration of GMAT that ANISE ports) are valid for orbits about the Earth only, use the J2 of the Earth and need the radii of the center. `:brouwer_inc` is in [0, π] and the other angles in [0, 2π). They throw `SpiceError(:degeneratecase)` if the center is not body 399, the periapsis is not above 3000 km, the orbit is not elliptic or a Kepler equation has no solution. For an inclination above 175° `:brouwer_aop` has a known error; see [Beyond CSPICE](@ref "Beyond CSPICE").

```@example orbits_ref
using Almanac, Markdown

# name => (unit, name in ANISE)
info = Dict(
    :x => ("km", "X"), :y => ("km", "Y"), :z => ("km", "Z"),
    :vx => ("km/s", "VX"), :vy => ("km/s", "VY"), :vz => ("km/s", "VZ"),
    :rmag => ("km", "Rmag"), :vmag => ("km/s", "Vmag"),
    :sma => ("km", "SemiMajorAxis"), :ecc => ("1", "Eccentricity"),
    :inc => ("rad", "Inclination"), :raan => ("rad", "RAAN"), :aop => ("rad", "AoP"),
    :ta => ("rad", "TrueAnomaly"), :ma => ("rad", "MeanAnomaly"),
    :ea => ("rad", "EccentricAnomaly"), :ha => ("rad", "HyperbolicAnomaly"),
    :aol => ("rad", "AoL"), :tlong => ("rad", "TrueLongitude"),
    :period => ("s", "Period"), :mean_motion => ("rad/s", "mean_motion_deg_s"),
    :ta_dot => ("rad/s", "ta_dot_deg_s"),
    :energy => ("km²/s²", "Energy"), :c3 => ("km²/s²", "C3"),
    :fpa => ("rad", "FlightPathAngle"), :semi_parameter => ("km", "SemiParameter"),
    :semi_minor_axis => ("km", "SemiMinorAxis"),
    :hmag => ("km²/s", "Hmag"), :hx => ("km²/s", "HX"), :hy => ("km²/s", "HY"),
    :hz => ("km²/s", "HZ"),
    :periapsis => ("km", "PeriapsisRadius"), :apoapsis => ("km", "ApoapsisRadius"),
    :periapsis_altitude => ("km", "PeriapsisAltitude"),
    :apoapsis_altitude => ("km", "ApoapsisAltitude"),
    :sma_altitude => ("km", "SemiMajorAxisAltitude"), :altitude => ("km", "Altitude"),
    :right_ascension => ("rad", "RightAscension"), :declination => ("rad", "Declination"),
    :velocity_declination => ("rad", "VelocityDeclination"),
    :latitude => ("rad", "Latitude"), :longitude => ("rad", "Longitude"),
    :height => ("km", "Height"), :light_time => ("s", "LightTime"),
    :bplane_bt => ("km", "BPlaneBT"), :bplane_br => ("km", "BPlaneBR"),
    :bplane_mag => ("km", "BPlaneMag"), :bplane_angle => ("rad", "BPlaneAngle"),
    :bplane_ltof => ("s", "BPlaneLTOF"),
    :equinoctial_a => ("km", "equinoctial_a_km"), :equinoctial_h => ("1", "EquinoctialH"),
    :equinoctial_k => ("1", "EquinoctialK"), :equinoctial_p => ("1", "EquinoctialP"),
    :equinoctial_q => ("1", "EquinoctialQ"), :equinoctial_lambda => ("rad", "EquinoctialLambda"),
    :brouwer_sma => ("km", "BrouwerMeanShortSemiMajorAxis"),
    :brouwer_ecc => ("1", "BrouwerMeanShortEccentricity"),
    :brouwer_inc => ("rad", "BrouwerMeanShortInclination"),
    :brouwer_raan => ("rad", "BrouwerMeanShortRAAN"),
    :brouwer_aop => ("rad", "BrouwerMeanShortAoP"),
    :brouwer_ma => ("rad", "BrouwerMeanShortMeanAnomaly"))
@assert Set(keys(info)) == Set(elements())

rows = ["| `:$n` | $(info[n][1]) | `$(info[n][2])` |" for n in elements()]
Markdown.parse(join(["| Name | Unit | ANISE |", "|:--|:--|:--|", rows...], "\n"))
```

## Building and editing

```@docs
keplerian_orbit
with_element
with_apsides
```

## Propagation

```@docs
two_body
```

## Local frames

```@docs
local_frame
local_frame_state
relative_state
```
