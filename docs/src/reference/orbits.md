# Orbits

An orbit is a state about a central body. These calls have no CSPICE counterpart; each docstring names the function of ANISE that it matches. Lengths are in km, times in seconds, angles in rad. [Orbits and relative motion](@ref) shows them in use.

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

The table lists every name of [`elements`](@ref) with its unit and the name of ANISE. It is built from `elements()`, so a name that is added without a row here stops the build.

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
