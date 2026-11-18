# Orbits and relative motion

An [`Orbit`](@ref) is a state about a central body, with the GM and the radii of that body. It gives the orbital elements, is built from elements, advances by two-body motion and gives the state of one orbit relative to another. CSPICE has the conversions (`oscltx`, `conics`, `prop2b`) but no such object; the definitions follow [ANISE](https://github.com/nyx-space/anise), with the [deviations from ANISE](@ref "Deviations from ANISE"). Lengths are in km, times in s, angles in rad.

```@setup orbits
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Read the elements

[`orbit`](@ref) takes the arguments of [`spkezr`](@ref), any aberration correction included, and the GM and radii of the observer from the kernel pool. The epoch of the orbit is `et`, the epoch of the observer:

```@example orbits
using Almanac

et = str2et("2025 MAR 18 12:00:00 TDB")
o = orbit("PROBE", et, "J2000", "NONE", "EARTH")
(o.sma, o.ecc)
```

`elements()` lists the names; each is a property, and [`element`](@ref) takes it as a symbol:

```@example orbits
rad2deg.((o.inc, o.raan, o.aop, o.ta))
```

Altitudes are above the mean equatorial radius; [Element names](@ref) gives the unit and range of each name:

```@example orbits
(o.periapsis_altitude, o.apoapsis_altitude, element(o, :period) / 3600)
```

[`keplerian`](@ref) returns the six classical elements as a named tuple. An element that is undefined throws `SpiceError(:degeneratecase)` rather than returning NaN; a few [fallback values](@ref "Fallback values") are returned instead.

An [`Orbit`](@ref) can also be made from a state vector and a GM, without kernels; the center and the frame are then NAIF IDs (399 is the Earth, 1 is `J2000`):

```@example orbits
leo0 = Orbit([7000.0, 0.0, 0.0, 0.0, 7.5, 1.0], et, 399, 1, 398600.4418)
(leo0.sma, leo0.ecc)
```

## Change one element

[`with_element`](@ref) replaces one classical element. Raise the inclination to 51.6°, keeping `sma`, `ecc`, `raan`, `aop` and `ta`; the kept elements change by rounding only:

```@example orbits
o2 = with_element(o, :inc, deg2rad(51.6))
(rad2deg(o2.inc), o2.sma - o.sma, o2.ta - o.ta)
```

[`keplerian_orbit`](@ref) builds an orbit from elements. With `gm` given it reads no kernels, so the orbit has no radii and its altitudes throw `SpiceError(:badarg)`. A nearly circular orbit with its semi-major axis 500 km above the equatorial radius:

```@example orbits
r = 6378.1366 + 500.0
leo = keplerian_orbit(et, 399, 1; sma=r, ecc=0.001, inc=deg2rad(98.0), raan=0.0,
    aop=0.0, ta=0.0, gm=398600.4418)
round(leo.period / 60; digits=2)
```

## Propagate by two-body motion

[`two_body`](@ref) moves an orbit to a new epoch (not by a step), forward or back, by Kepler motion without perturbations. The kernels of the manual hold two-body motion, so the result matches the SPK six hours later, here to a few µm (km, km/s):

```@example orbits
later = et + 6 * 3600.0
p = two_body(o, later)
q = orbit("PROBE", later, "J2000", "NONE", "EARTH")
maximum(abs.(p.state - q.state))
```

## Compute a relative state

[`relative_state`](@ref) gives the state of one orbit minus a chief orbit, in a local frame of the chief: `:ric` (radial, in-track, cross-track), `:vnc`, `:rcn` or `:tnw`, whose axes [`local_frame`](@ref) lists. Both orbits must have the same epoch, center and frame, so advance both:

```@example orbits
chief = two_body(o, later)
deputy = two_body(o2, later)
rel = relative_state(deputy, chief, :ric)
(radial = rel[1], in_track = rel[2], cross_track = rel[3])
```

## Limits of the orbit calls

- The frame is not checked. Keplerian, B-plane and local-frame results assume an inertial frame and are wrong, without an error, in a rotating one such as `IAU_EARTH`.
- The velocity of `relative_state` (`rel[4:6]`) is the rate of the offset seen in the rotating local frame, not the difference of the inertial velocities; it is not zero when those are equal. Orbits of different epoch, center or frame throw `SpiceError(:badarg)`.
- A circular orbit (`ecc < 1e-11`) has no true anomaly: `:ta`, `:ma`, `:ea`, `:fpa`, `:tlong`, `:equinoctial_lambda`, the Brouwer elements, `keplerian`, [`with_element`](@ref) (except for `:ta`), [`with_apsides`](@ref) and `two_body` throw. The state, `:sma`, `:ecc`, `:inc`, `:raan` and `:aop` (0 for an exactly circular state), `:aol` (if not equatorial), `:period`, the altitudes, the local frames and `relative_state` still work. A `ForwardDiff.Dual` state on an exactly circular orbit throws for `:ecc` too.
- Brouwer mean elements need an elliptic orbit about the Earth (399) with a periapsis radius above 3000 km; the B-plane elements need a hyperbola.

## Apsides, flybys and other elements

[`with_apsides`](@ref) sets the apsis radii, [`vinf_periapsis`](@ref) and [`vinf_turn_angle`](@ref) give the flyby geometry of a hyperbola, and `element` also returns B-plane, equinoctial and Brouwer mean elements.

## Differentiate an orbit

A `ForwardDiff.Dual` epoch, state or element passes through [these calls](@ref "Differentiate with ForwardDiff"). The derivative of the position of `two_body` with respect to the epoch is the velocity:

```@example orbits
using ForwardDiff

v = ForwardDiff.derivative(t -> two_body(o, t).position, later)
maximum(abs.(v - p.velocity))
```

The derivative throws where the value throws, as for the anomalies at `ecc < 1e-11`. A dual state of exactly zero radius, velocity, angular momentum or eccentricity vector throws `SpiceError(:degeneratecase)`, where a `Float64` state gives the value 0. The derivative of `:ta`, and of the calls built on it, is NaN at an exact apsis, where the arc cosine has an infinite slope; the value is not NaN.

The other AD packages differentiate the same calls; see [Differentiate with other AD packages](@ref).
