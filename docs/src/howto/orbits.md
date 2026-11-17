# Orbits and relative motion

An [`Orbit`](@ref) is a state about a central body with the GM and the radii of that body. It gives the elements of the orbit, builds an orbit from elements, changes one element, advances the orbit by two-body motion, and expresses the state of one orbit relative to another. CSPICE has the conversions (`oscltx`, `conics`, `prop2b`) but no object that carries them; the definitions here follow ANISE, and [Beyond CSPICE](@ref) says what was checked and where Almanac.jl deviates. Epochs are TDB seconds past J2000, lengths are in km, times in s and angles in rad throughout.

```@setup orbits
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Read the elements

[`orbit`](@ref) takes the arguments of [`spkezr`](@ref) and returns the state of `PROBE` relative to the Earth in `J2000`, as `spkezr` gives it; the state is not rotated into another frame. The epoch of the orbit is `et`, the epoch of the observer. The GM comes from `BODY399_GM` of the kernels, and the radii from `BODY399_RADII` when it is there; without a GM the call throws `SpiceError(:kernelvarnotfound)`.

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

The altitudes are measured from the mean equatorial radius of the Earth, `(radii[1] + radii[2]) / 2`; [Element names](@ref) lists the units and the behaviour of every name:

```@example orbits
(o.periapsis_altitude, o.apoapsis_altitude, element(o, :period) / 3600)
```

[`keplerian`](@ref) returns the six classical elements as a named tuple. An element that is undefined for the orbit throws `SpiceError(:degeneratecase)`. This matters for circular orbits: with `ecc < 1e-11` the true anomaly is undefined, so `:ta`, [`keplerian`](@ref), [`with_element`](@ref) (for every name but `:ta`) and [`two_body`](@ref) all throw. A nearly circular orbit, as in the example below, works.

## Change one element

[`with_element`](@ref) returns a copy of the orbit with one of `:sma`, `:ecc`, `:inc`, `:raan`, `:aop` or `:ta` replaced. Raise the inclination of `PROBE` to 51.6°, keeping the shape and the position along the orbit:

```@example orbits
o2 = with_element(o, :inc, deg2rad(51.6))
(rad2deg(o2.inc), o2.sma == o.sma, o2.ecc == o.ecc)
```

[`keplerian_orbit`](@ref) builds an orbit from elements without a state. Its arguments after the epoch are the center and the frame as NAIF IDs: `399` is the Earth and `1` is `J2000`. With `gm` given it reads no kernels, so the orbit has no radii: the altitudes, `:latitude` and `:height` throw `SpiceError(:badarg)`. This orbit is nearly circular, 500 km above the mean equatorial radius:

```@example orbits
r = 6378.1366 + 500.0
leo = keplerian_orbit(et, 399, 1; sma=r, ecc=0.001, inc=deg2rad(98.0), raan=0.0,
    aop=0.0, ta=0.0, gm=398600.4418)
round(leo.period / 60; digits=2)
```

## Propagate

[`two_body`](@ref) advances an orbit to any other epoch, earlier or later, by Kepler motion about the center. The kernels of the manual hold two-body motion, so the result can be compared with the state from the SPK six hours later:

```@example orbits
later = et + 6 * 3600.0
p = two_body(o, later)
q = orbit("PROBE", later, "J2000", "NONE", "EARTH")
maximum(abs.(p.state - q.state))
```

The two states agree to this difference (km and km/s). With a force model in the data, `two_body` and the SPK would differ by the perturbations.

## Relative motion

[`relative_state`](@ref) gives the state of one orbit minus the other, the chief, in a local frame of the chief. The kinds are `:ric` (radial, in-track, cross-track), `:vnc`, `:rcn` and `:tnw`; the axes of each are in the table of [`local_frame`](@ref), which returns the rotation itself, and [`local_frame_state`](@ref) returns the matrix for states. The two orbits must have the same epoch, center and frame, or the call throws `SpiceError(:badarg)`, so advance both:

```@example orbits
chief = two_body(o, later)
deputy = two_body(o2, later)
rel = relative_state(deputy, chief, :ric)
(radial = rel[1], in_track = rel[2], cross_track = rel[3])
```

The velocity part `rel[4:6]` is the rate of change of the offset as seen in the rotating frame, so it is not zero even when the two inertial velocities are equal.

## More calls

[`with_apsides`](@ref) sets the apsis radii, [`vinf_periapsis`](@ref) and [`vinf_turn_angle`](@ref) give the flyby geometry of a hyperbola, and `element` also returns the B-plane, the equinoctial and the Brouwer mean elements (the last for orbits about the Earth only); [Element names](@ref) describes them.

## Dual numbers

With `using ForwardDiff`, a `ForwardDiff.Dual` epoch or state passes through these calls (the GM and the radii stay `Float64`), so the derivative of an element or of a propagated state is one `ForwardDiff` call. The derivative of the position of `two_body` with respect to the epoch is the velocity:

```@example orbits
using ForwardDiff

v = ForwardDiff.derivative(t -> two_body(o, t).position, later)
maximum(abs.(v - p.velocity))
```

Enzyme, Mooncake, Zygote, ReverseDiff and Tracker differentiate the same calls; see [Other AD packages](@ref) in [Time derivatives](@ref), which also lists the other derivative calls.
