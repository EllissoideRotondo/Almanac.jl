# Orbits

These calls have no CSPICE counterpart; they follow [ANISE](https://github.com/nyx-space/anise), and [Validation of the orbit calls](@ref) gives their validation and the [deviations from ANISE](@ref "Deviations from ANISE"). Lengths are in km, times in s, angles in rad. [Orbits and relative motion](@ref) shows them in use.

## Orbit objects

```@docs
Orbit
orbit
```

## Reading elements

```@docs
element
elements
keplerian
vinf_periapsis
vinf_turn_angle
```

## Element names

- **Angles.** `:inc` is in [0, π]; `:ta`, `:aop`, `:raan`, `:aol`, `:tlong`, `:ha` and `:right_ascension` in [0, 2π], where `:raan` and `:aop` can be exactly 2π; `:ea` in (-π, π]; the elliptic `:ma` in [0, 2π], the hyperbolic one unbounded; `:declination` in [-π/2, π/2]. Near a circular orbit (`ecc` down to 1e-9) the arc cosine amplifies rounding, and `:ta`, `:aop`, `:ma` and `:ea` agree with ANISE to 3e-8 rad.
- **Undefined** (`SpiceError(:degeneratecase)`): `:ta`, `:ea`, `:ma`, `:fpa`, `:tlong` and `:ha` for `ecc < 1e-11`; `:aol` within 1e-6° of the equator; `:ma` for `|ecc - 1| < 1e-11`; `:ea` for `ecc > 1`; `:ha` for `ecc <= 1`; whatever needs the angular momentum (zero position or velocity), the eccentricity vector or `:energy` (zero position; a zero velocity gives `ecc = 1`).
- **B-plane** (`:bplane_bt` to `:bplane_ltof`): for a hyperbola in an inertial frame, with the triad built on the `+z` axis of the frame; `:bplane_angle` is in [-π, π). They throw for `ecc <= 1` and for an incoming asymptote within 1e-12 of the `z` axis; `:bplane_ltof` also for `c3 <= 0`.
- **Rates.** `:mean_motion` is `sqrt(gm / |sma|³)`, `:ta_dot` is `2π / period` (ANISE's definition: a mean rate).
- **Altitudes** are a radius (`:rmag`, `:sma`, an apsis) minus the mean equatorial radius `(radii[1] + radii[2]) / 2`. `:latitude` (geodetic, in [-π/2, π/2]) and `:height` are on the spheroid with that equatorial radius and polar radius `radii[3]`, in the body-fixed frame of the center; `:longitude` is `atan(y, x)` in [-π, π).
- **Equinoctial**: `:equinoctial_h`, `:equinoctial_k` are `ecc * sin(aop + raan)`, `ecc * cos(aop + raan)`; `:equinoctial_p`, `:equinoctial_q` are `sin(inc / 2) * sin(raan)`, `sin(inc / 2) * cos(raan)`; `:equinoctial_lambda` is `raan + aop + ma`, not wrapped; they throw where the elements they use throw.
- **Brouwer mean elements** (short-period terms, the GMAT iteration that ANISE ports) use the J2 of the Earth. `:brouwer_inc` is in [0, π], the other angles in [0, 2π). They throw `SpiceError(:degeneratecase)` unless the center is body 399, the periapsis radius is above 3000 km and the orbit is elliptic, or if a Kepler equation has no solution. Above 175° of inclination `:brouwer_aop` has a [known error](@ref "Deviations from ANISE").

### Fallback values

These are returned instead of an error:

- `:raan` and `:aop` are 0 where their arc cosine has no value (`:raan` of an equatorial orbit) and π where its argument rounds above 1, as in ANISE.
- `:sma` is `±Inf` on a parabola.
- `:period` is 0 and `:ta_dot` is `Inf` on an orbit that is not elliptic.
- [`keplerian_orbit`](@ref) changes the sign of a negative `ecc`, and of an `sma` of the wrong sign for `ecc`.

### Units and ANISE names

The names in the order of `elements()`, with the unit in Almanac.jl and the name in ANISE, whose angles are in degrees:

```@eval
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

## Building and changing orbits

```@docs
keplerian_orbit
with_element
with_apsides
```

## Two-body propagation

```@docs
two_body
```

## Local frames

```@docs
local_frame
local_frame_state
relative_state
```
