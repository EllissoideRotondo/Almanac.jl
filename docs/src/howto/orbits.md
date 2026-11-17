# Orbits and relative motion

An [`Orbit`](@ref) is a state about a central body, with the GM and the radii of that body. It gives the orbital elements, is built from elements, advances by two-body motion and gives the state of one orbit relative to another. CSPICE has the conversions (`oscltx`, `conics`, `prop2b`) but no such object; the definitions follow [ANISE](https://github.com/nyx-space/anise), with the deviations listed in [Beyond CSPICE](@ref). Lengths are in km, times in s, angles in rad.

```@setup orbits
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Read the elements

[`orbit`](@ref) takes the arguments of [`spkezr`](@ref), and the GM and radii of the observer from the kernel pool:

```@example orbits
using Almanac

et = str2et("2025 MAR 18 12:00:00 TDB")
o = orbit("PROBE", et, "J2000", "NONE", "EARTH")
(o.sma, o.ecc)
```

Each name of [`elements`](@ref) is a property, and [`element`](@ref) takes it as a symbol:

```@example orbits
rad2deg.((o.inc, o.raan, o.aop, o.ta))
```

Altitudes are above the mean equatorial radius; [Element names](@ref) gives the unit and range of each name:

```@example orbits
(o.periapsis_altitude, o.apoapsis_altitude, element(o, :period) / 3600)
```

[`keplerian`](@ref) returns the six classical elements as a named tuple. An element that is undefined throws `SpiceError(:degeneratecase)` rather than returning NaN. A circular orbit (`ecc < 1e-11`) has no true anomaly, so `:ta`, `keplerian`, [`with_element`](@ref) (except for `:ta`) and [`two_body`](@ref) throw for it; a nearly circular orbit works.

## Change one element

[`with_element`](@ref) replaces one classical element. Raise the inclination to 51.6°, keeping the shape and the position along the orbit:

```@example orbits
o2 = with_element(o, :inc, deg2rad(51.6))
(rad2deg(o2.inc), o2.sma == o.sma, o2.ecc == o.ecc)
```

[`keplerian_orbit`](@ref) builds an orbit from elements; here about the Earth (399) in `J2000` (1). With `gm` given it reads no kernels, so the orbit has no radii and its altitudes throw. A nearly circular orbit 500 km above the equator:

```@example orbits
r = 6378.1366 + 500.0
leo = keplerian_orbit(et, 399, 1; sma=r, ecc=0.001, inc=deg2rad(98.0), raan=0.0,
    aop=0.0, ta=0.0, gm=398600.4418)
round(leo.period / 60; digits=2)
```

## Propagate

[`two_body`](@ref) advances an orbit, forward or back, by Kepler motion. The kernels of the manual hold two-body motion, so the result matches the SPK six hours later (km, km/s):

```@example orbits
later = et + 6 * 3600.0
p = two_body(o, later)
q = orbit("PROBE", later, "J2000", "NONE", "EARTH")
maximum(abs.(p.state - q.state))
```

## Relative motion

[`relative_state`](@ref) gives the state of one orbit minus a chief orbit, in a local frame of the chief: `:ric` (radial, in-track, cross-track), `:vnc`, `:rcn` or `:tnw`, whose axes [`local_frame`](@ref) lists. Both orbits must have the same epoch, center and frame:

```@example orbits
chief = two_body(o, later)
deputy = two_body(o2, later)
rel = relative_state(deputy, chief, :ric)
(radial = rel[1], in_track = rel[2], cross_track = rel[3])
```

## More calls

[`with_apsides`](@ref) sets the apsis radii, [`vinf_periapsis`](@ref) and [`vinf_turn_angle`](@ref) give the flyby geometry of a hyperbola, and `element` also returns B-plane, equinoctial and Brouwer mean elements.

## Dual numbers

A `ForwardDiff.Dual` epoch or state passes through these calls. The derivative of the position of `two_body` with respect to the epoch is the velocity:

```@example orbits
using ForwardDiff

v = ForwardDiff.derivative(t -> two_body(o, t).position, later)
maximum(abs.(v - p.velocity))
```

The other AD packages differentiate the same calls; see [Other AD packages](@ref).
