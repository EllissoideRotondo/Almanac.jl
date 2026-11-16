# Orbits and relative motion

An [`Orbit`](@ref) is a state about a central body with the GM and the radii of that body. It gives the elements of the orbit, builds an orbit from elements, changes one element, advances the orbit by two-body motion, and expresses the state of one orbit relative to another. None of this is in CSPICE; the definitions follow ANISE, and [Beyond CSPICE](@ref "Coverage and validation") says what was checked. Lengths are in km and angles in rad throughout.

```@setup orbits
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Read the elements

[`orbit`](@ref) takes the arguments of [`spkezr`](@ref) and returns the state of `PROBE` relative to the Earth, with the GM and the radii of the Earth from the kernels:

```@example orbits
using Almanac

et = str2et("2025 MAR 18 12:00:00 TDB")
o = orbit("PROBE", et, "J2000", "NONE", "EARTH")
(o.sma, o.ecc)
```

Each name in [`elements`](@ref) is a property of the orbit, and [`element`](@ref) takes the name as a symbol. The angles are in rad:

```@example orbits
rad2deg.((o.inc, o.raan, o.aop, o.ta))
```

The altitudes use the mean equatorial radius of the Earth, so the altitudes of the two apsides are those that the fixture was written with:

```@example orbits
(o.periapsis_altitude, o.apoapsis_altitude, element(o, :period) / 3600)
```

[`keplerian`](@ref) returns the six classical elements as a named tuple. An element that is undefined for the orbit, such as the true anomaly of a circular orbit, throws `SpiceError(:degeneratecase)`; ANISE returns NaN or an error there.

## Change one element

[`with_element`](@ref) returns a copy of the orbit with one of `:sma`, `:ecc`, `:inc`, `:raan`, `:aop` or `:ta` replaced. Raise the inclination of `PROBE` to 51.6°, keeping the shape and the position along the orbit:

```@example orbits
o2 = with_element(o, :inc, deg2rad(51.6))
(rad2deg(o2.inc), o2.sma == o.sma, o2.ecc == o.ecc)
```

[`keplerian_orbit`](@ref) builds an orbit from elements without a state. A circular orbit about the Earth, 500 km above the mean equatorial radius, needs the GM only:

```@example orbits
r = 6378.1366 + 500.0
leo = keplerian_orbit(et, 399, 1; sma=r, ecc=0.001, inc=deg2rad(98.0), raan=0.0,
    aop=0.0, ta=0.0, gm=398600.4418)
round(leo.period / 60; digits=2)
```

## Propagate

[`two_body`](@ref) advances an orbit to a later epoch by Kepler motion about the center. The kernels of the manual hold two-body motion, so the result can be compared with the state from the SPK six hours later:

```@example orbits
later = et + 6 * 3600.0
p = two_body(o, later)
q = orbit("PROBE", later, "J2000", "NONE", "EARTH")
maximum(abs.(p.state - q.state))
```

The two states agree to about 1e-8 (km and km/s). With a force model in the data, `two_body` and the SPK would differ by the perturbations.

## Relative motion

[`relative_state`](@ref) gives the state of one orbit relative to another, the chief, in a local frame of the chief. The kinds are `:ric` (radial, in-track, cross-track), `:vnc`, `:rcn` and `:tnw`; [`local_frame`](@ref) and [`local_frame_state`](@ref) return the rotation itself. The two orbits must have the same epoch, center and frame, so advance both:

```@example orbits
chief = two_body(o, later)
deputy = two_body(o2, later)
rel = relative_state(deputy, chief, :ric)
(radial = rel[1], in_track = rel[2], cross_track = rel[3])
```

The velocity part `rel[4:6]` is the rate of change of that offset as seen in the rotating frame, because the rate of the frame is included.

## Numbers other than `Float64`

A `ForwardDiff.Dual` epoch or state passes through every call, so the derivative of an element or of a propagated state is one `ForwardDiff` call. The derivative of the position of `two_body` with respect to the epoch is the velocity:

```@example orbits
using ForwardDiff

v = ForwardDiff.derivative(t -> two_body(o, t).position, later)
maximum(abs.(v - p.velocity))
```

See [Time derivatives](@ref) for the other derivative calls.
