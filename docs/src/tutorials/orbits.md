# Orbits and relative motion

An [`Orbit`](@ref) is a state about a central body, with the GM (km³/s²) and the radii (km) of that body. It gives the orbital elements, is built from elements, advances by two-body motion and gives the state of one orbit relative to another. CSPICE has the conversions (`oscltx`, `conics`, `prop2b`) but no such object; the definitions follow [ANISE](https://github.com/nyx-space/anise), with the [differences from ANISE](@ref "Differences from ANISE"). Lengths are in km, times in s, angles in rad, where ANISE uses degrees.

## Read the elements

[`orbit`](@ref) takes the arguments of [`spkezr`](@ref), any aberration correction included, and the GM and radii of the observer from the kernel pool, here from `gm_de440.tpc` and `pck00011.tpc`. The epoch of the orbit is `et`, the epoch of the observer:

```@example orbits
using Almanac, SPICEKernels                 # SPICEKernels downloads the NAIF generic kernels

furnsh(naif0012_tls())
furnsh(de440s())
furnsh(pck00011())
furnsh(gm_de440())

et = str2et("2025 MAR 18 12:00:00 TDB")
o = orbit("MOON", et, "J2000", "NONE", "EARTH")
(o.sma, o.ecc)
```

`elements()` lists the names; each is a property, and [`element`](@ref) takes it as a symbol:

```@example orbits
rad2deg.((o.inc, o.raan, o.aop, o.ta))
```

Altitudes are above the mean equatorial radius; [Element names](@ref) gives the unit and range of each name:

```@example orbits
(o.periapsis_altitude, o.apoapsis_altitude, element(o, :period) / 86400)
```

[`keplerian`](@ref) returns the six classical elements as a named tuple. An element that is undefined throws `SpiceError(:degeneratecase)` rather than returning NaN; a few [fallback values](@ref "Fallback values") are returned instead.

An [`Orbit`](@ref) can also be made from a state vector and a GM, without kernels; the center and the frame are then NAIF IDs (399 is the Earth, 1 is `J2000`):

```@example orbits
leo0 = Orbit([7000.0, 0.0, 0.0, 0.0, 7.5, 1.0], et, 399, 1, 398600.4418)
(leo0.sma, leo0.ecc)
```

## Propagate by two-body motion

[`two_body`](@ref)`(o, t)` returns the orbit at the epoch `t`, before or after the epoch of `o`, by Kepler motion without perturbations. `t` is an epoch, not a time step. The Sun perturbs the orbit of the Moon, so six hours later the result departs from the ephemeris (km, km/s):

```@example orbits
later = et + 6 * 3600.0
p = two_body(o, later)
q = orbit("MOON", later, "J2000", "NONE", "EARTH")
maximum(abs.(p.state - q.state))
```

## Change one element

[`keplerian_orbit`](@ref) builds an orbit from elements. With `gm` given it reads no kernels, so the orbit has no radii and its altitudes throw `SpiceError(:badarg)`. A nearly circular orbit with its semi-major axis 500 km above the equatorial radius of the Earth:

```@example orbits
r = 6378.1366 + 500.0
leo = keplerian_orbit(et, 399, 1; sma=r, ecc=0.001, inc=deg2rad(98.0), raan=0.0,
    aop=0.0, ta=0.0, gm=398600.4418)
round(leo.period / 60; digits=2)
```

[`with_element`](@ref) replaces one classical element. Raise the inclination by 0.1°, keeping `sma`, `ecc`, `raan`, `aop` and `ta`; the kept elements change by rounding only:

```@example orbits
leo2 = with_element(leo, :inc, leo.inc + deg2rad(0.1))
(rad2deg(leo2.inc), leo2.sma - leo.sma, leo2.ta - leo.ta)
```

## Compute a relative state

[`relative_state`](@ref) gives the state of one orbit, the deputy, minus that of a reference orbit, the chief, in a local frame of the chief: `:ric` (radial, in-track, cross-track), `:vnc`, `:rcn` or `:tnw`, whose axes [`local_frame`](@ref) lists. Both orbits must have the same epoch, center and frame, so advance both:

```@example orbits
chief = two_body(leo, later)
deputy = two_body(leo2, later)
rel = relative_state(deputy, chief, :ric)
(radial = rel[1], in_track = rel[2], cross_track = rel[3])
```

## Limits of the orbit calls

- The frame is not checked. Keplerian, B-plane and local-frame results assume an inertial frame and are wrong, without an error, in a rotating one such as `IAU_EARTH`.
- The velocity of `relative_state` (`rel[4:6]`) is the rate of the offset seen in the rotating local frame, not the difference of the inertial velocities; it is not zero when those are equal. Orbits of different epoch, center or frame throw `SpiceError(:badarg)`.
- A circular orbit (`ecc < 1e-11`) has no true anomaly, so the elements that need one, `keplerian`, [`with_element`](@ref) (except for `:ta`), [`with_apsides`](@ref) and `two_body` throw; the state, the other elements, the local frames and `relative_state` work. [Element names](@ref) lists the undefined cases.
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

The derivative throws where the value throws, as for the anomalies at `ecc < 1e-11`. A dual state of exactly zero radius, velocity, angular momentum or eccentricity vector throws `SpiceError(:degeneratecase)`, where a `Float64` state gives the value 0. `:ta`, `:raan` and `:aop` have the value of the arc cosine of ANISE and the derivative of the same angle as an arc tangent, so the derivative is finite at an apsis and at a node, where the arc cosine has an infinite slope.

The other AD packages differentiate the same calls; see [Differentiate with other AD packages](@ref).
