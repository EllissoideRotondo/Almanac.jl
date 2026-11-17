# Validation of the orbit calls

The calls of [Orbits and relative motion](@ref) are not CSPICE routines, so the
figures of [Known limits of accuracy](@ref) do not apply to them. They
follow the source of ANISE main (commit `b3aaadc`) function by function, with its error
and fallback cases. The released ANISE 0.10.6 and 1.0.0b0 differ from that source in
the flight path angle, inbound hyperbolas, the `TNW` frame and v-infinity.

Where CSPICE computes the same quantity (`oscltx`, `recgeo`, `conics`, `prop2b`),
CSPICE is the reference for the numbers. For everything else, and for what is an error
or a fallback, the ANISE source is the reference. Only the check [against
ANISE](@ref "Against ANISE") compares Almanac.jl with the output of ANISE itself.

## How to read the figures

- Differences between two numbers are relative to the larger one, unless stated
  otherwise.
- Each check is a test of the group `analysis_orbit`; sampled inputs come from a fixed
  seed.
- The figures are those of the default run. Figures in brackets are the larger of two
  exhaustive runs (`ALMANAC_TEST_FULL=1`, seeds 20261011 and 1, which sample 20 to 35
  times more); for the check against ANISE, they are one exhaustive run.
- A difference of an angle or of the shape is given in units of the error model of the
  sampled orbit, because the geometry of the sample amplifies the rounding of an arc
  cosine or of `1 - e`. The models are defined in `test/test_analysis_orbit.jl`:
  - `cond_shape` is `eps (1 + 1/|1 - e|)`;
  - `cond_angle` is `eps` times the amplifier of the angle (`1/ecc`, `1/sin inc`,
    `1/|sin angle|`);
  - `cond_kepler` is `1e-12 max(1, 1/sqrt(|1 - e²|))`, the stopping rule of the
    Kepler iteration;
  - `cond_state` is the larger of the amplifiers above, for a state rebuilt from its
    elements.

## Against CSPICE

- **`oscltx`.** `sma`, `ecc`, `inc`, `raan`, `aop`, `ma`, `ta` and `period` of 200
  random orbits. In units of `cond_shape`: `periapsis` 3.2 [8.4], `ecc` 9.5 [16],
  `sma` and `period` 3.2 [5.3]. In units of `cond_angle`: `inc` 1.0 [1.8], `raan` 1.9
  [2.0], `aop` 2.2 [3.3], `ta` 4.0 [6.0], `ma` 1.9 [4.8], `aol` 2.9 [3.0]. The larger
  raw differences of the anomalies are this amplification at the sampled orbits, not
  an error of a formula.
- **`recgeo`.** `latitude`, `longitude` and `height` at 900 random points about three
  bodies, one of them flattened by 0.1: latitude 3.9e-16 rad, longitude 4.4e-16 rad,
  height 1.1e-15 of the equatorial radius.
- **`conics`.** [`keplerian_orbit`](@ref) against the state of `conics` for 200 random
  orbits, and the round trip of the state. Position 1.2 [1.6] and velocity 0.96 [1.0]
  in units of `cond_kepler`; round trip of the state 8.8 [8.8] (`cond_state`), of `sma`
  and `ecc` 4.4e-15 [8.7e-15], of the angles 5.2 [9.1] (`cond_angle`).
- **`prop2b`.** [`two_body`](@ref) over ±10 periods. Position 1.3 [3.2] and velocity
  0.92 [5.5] in units of `cond_kepler`; `sma` 5.8e-15 [9.9e-15].

## Against published values

- **GMAT cases.** The Keplerian elements of three pairs of Cartesian and Keplerian
  states (circular inclined, elliptical, circular equatorial), the equatorial case, the
  two hyperbolic cases and 24 flight path angles, with the values that GMAT gave to the
  tests of ANISE; angles modulo 2π. Difference scaled by `max(|value|, 1)`: at most
  4.0e-13 (the limit is 1e-12); flight path angles 2.2e-16.
- **GMAT B-plane.** `bplane_bt`, `bplane_br`, `bplane_mag` and `bplane_angle` of one
  hyperbolic case. Absolute difference from GMAT: components 0.65 km, magnitude 0.23
  km, angle 2.0e-5 rad. The test of ANISE allows 1 km and 0.01°. The cause is not
  isolated.
- **Vallado.** Longitude, geodetic latitude, height and the altitudes of Examples 3-3
  and 3-4 of Vallado, with the radii that ANISE uses. Absolute: 2.2e-16 rad (latitude),
  1.1e-16 rad (longitude), 2.1e-12 km (height), 4.2e-12 km (altitudes).
- **Davis.** The hyperbolic anomaly of the state of the lecture of Davis, as ANISE
  tests it, and [`vinf_periapsis`](@ref) as the inverse of
  [`vinf_turn_angle`](@ref): hyperbolic anomaly 1.1e-11 rad, inverse 1.1e-13 km.
- **Test cases of ANISE.** The round trips of Brouwer mean elements (osculating, mean,
  osculating), [`with_element`](@ref) and [`with_apsides`](@ref) on the loops of
  `verif_with_init`, and `two_body` on `verif_orbit_at_epoch`, with the inputs and
  limits of the tests of ANISE:
  - Brouwer, low orbit: `sma` 3.1e-11 km, angles up to 3e-11 rad;
  - Brouwer, transfer orbit: `sma` 2.4e-9 km;
  - Brouwer, retrograde orbit: `sma` 3.0e-3 km, `ecc` 4.1e-7, `inc` 3.8e-6, `raan`
    4.3e-4, `aop` 1.09e-2, `ma` 3.6e-5 rad; the `aop` figure is a defect of ANISE,
    reproduced (see [Deviations from ANISE](@ref));
  - `with_element` and `with_apsides`: at most 1.5e-11;
  - `verif_orbit_at_epoch`: `aop` 3.3e-8 rad (ANISE allows 5.2e-8).

## Against ANISE

Every element, local frame and relative state of 12,000 random states about the Earth
and Mars (2000 each of elliptic, near-circular, near-equatorial, retrograde, outbound
hyperbolic and inbound hyperbolic), and [`keplerian_orbit`](@ref),
[`with_element`](@ref), [`with_apsides`](@ref) and [`two_body`](@ref) on the same
states, with the same state, epoch, GM and radii given to both libraries. In brackets,
one exhaustive run of 120,000 states. Differences are scaled by `max(|value|, 1)`,
angles are taken modulo 2π, and the range convention of every angle is checked.

- Most elements: at most 1.1e-12, reached by `period` [`tlong`, 8.4e-12].
  `light_time`: 1.0e-9 s, because ANISE rounds a duration to 1 ns.
- `ta`, `aop`, `ma`, `ea`, and the same angles after `with_element` or
  `with_apsides`: 2.2e-8 [2.7e-8], at the near-circular states (`ecc` down to 1e-9),
  where the arc cosine amplifies rounding.
- Brouwer elements: 2.6e-10 [3.6e-10], worst at a near-equatorial state.
- Local frames ([`local_frame`](@ref)) and the position of
  [`relative_state`](@ref): identical.
- `keplerian_orbit` 4e-14 [1.3e-12], `with_element` 1.1e-11 [4.8e-11],
  `with_apsides` 1.4e-8 [1.2e-8].
- `two_body`: 3.5e-9 [7.2e-9] for orbits that are not hyperbolic, 3.2e-11 [8.3e-11]
  for hyperbolas with `raan + aop` inside (-π, π].
- Rates of the local frames and velocity of `relative_state`: 1.3e-6 [8.0e-6]. This
  is at most 1.5 times 1e-12/(2 ms), the error that the stopping rule of the Kepler
  iteration of ANISE (1e-12 rad) gives to its 1 ms difference quotient, and 0.02 times
  that for outbound hyperbolas.

Not compared, because 1.0.0b0 defines them differently from the source that Almanac.jl
follows: `fpa`, `vinf_periapsis`, `vinf_turn_angle`, the inbound hyperbolic results of
the constructors and of the rates, the `TNW` frame, and the B-plane calls (absent in
1.0.0b0). Reported but not tested: `two_body` and the rates of hyperbolas with
`raan + aop > π`, a deliberate deviation described below.

This check is the opt-in group `test/anise/`, run on 2026-11-15 against ANISE
1.0.0b0. It installs Python and the `anise` wheel through CondaPkg:

```sh
julia --startup-file=no --project=test/anise test/anise/runtests.jl
```

## Internal checks

- **`with_element`.** The element that is set, and the five that must not change, over
  random orbits. The set `sma` (relative) or `ecc`: 13 [26] units of `cond_shape`; the
  set angle: 2.8 [6.6] units of `cond_angle`. The others: `sma` or `ecc` 1.4e-14
  [3.3e-14], angles 6.0 [14].
- **Identities.** `aol = aop + ta`, `c3 = 2 energy`, `mean_motion`, `ta_dot`, right
  ascension and declination, the hyperbolic functions of `ha`, the B-plane and
  v-infinity formulas, the conversion back from equinoctial elements,
  `relative_state(o, o) == 0`, orthonormal local frames, and the inverse of
  [`local_frame_state`](@ref):
  - scalar identities: at most 1.6e-15 [4.4e-15];
  - B-plane and v-infinity: 1.2e-15 [5.5e-15];
  - `aol`: 2.9 [3.0] units of `cond_angle`;
  - equinoctial conversion: 0.79 [1.5] units of `eps (16 + 1/sin inc)`;
  - orthonormality: 4.4e-16 [8.9e-16].
- **Derivatives with ForwardDiff.**
  - The exact rate of `local_frame_state` against `ForwardDiff.derivative` of the frame
    of `two_body`: 6.9 [51] in units of `eps |h|/r²`, the angular rate of the frame.
  - The derivative of the state of `two_body` with respect to a `Dual` epoch, against
    the velocity and acceleration of two-body motion: 3.7e-15 [8.6e-15] for the
    position, 3.8e-15 [9.2e-15] for the velocity.
  - The derivatives of the true anomaly through the Kepler solver (`Dual` mean anomaly
    and eccentricity, nested for the second derivative) against their closed forms:
    `dν/dM` 9.0e-15 [1.5e-14], `dν/de` 1.7e-15 [3.9e-15], `d²ν/dM²` 2.7e-12 [7.7e-12]
    relative to its bound.
- **Other AD packages against ForwardDiff.** Enzyme (forward and reverse), Mooncake
  (forward and reverse), Zygote, ReverseDiff and Tracker, in the groups `ad_enzyme`,
  `ad_mooncake`, `ad_zygote` and `ad_tracked`. The default run compares the gradients
  of `keplerian_orbit` with respect to the elements (ellipse, hyperbola, and mean
  anomaly π, the branch of the Kepler solver where ANISE returns E for ν) and the
  derivative of the local frame of `orbit` with respect to the epoch. The exhaustive
  tier adds `two_body`; the elements `sma`, `ecc`, `inc`, `raan`, `aop`, `ta`, `ma`,
  `latitude`, `bplane_bt` and `brouwer_sma`; `local_frame`; `relative_state` with
  respect to the state; and `:sma` of `orbit` with respect to the epoch. Results: at
  most 1.6e-15 for the gradients of the default run and 2.1e-15 over all cases (`:ta`
  with Tracker), and 1.3e-14 for `:sma` of `orbit` with respect to the epoch. The
  default limits are 5e-14 for the gradients and 1e-13 for the epoch.
- **Allocation.** AllocCheck.jl proves, for `Float64` inputs, that
  [`element`](@ref), [`keplerian`](@ref), `with_element`, `two_body`, `local_frame`,
  `local_frame_state` and `relative_state` cannot allocate. `keplerian_orbit`,
  `with_apsides`, [`Orbit`](@ref) and [`orbit`](@ref) are not in the proof.

## Deviations from ANISE

One defect of ANISE is reproduced, and three deviations from ANISE are deliberate.

- **Retrograde Brouwer elements** (the defect reproduced). As in ANISE, Brouwer
  elements are those of the Earth, and a center other than body 399 throws
  `SpiceError(:degeneratecase)`. For an inclination above 175° the
  iteration works on a pseudo-state, and ANISE maps back the inclination and `raan`
  but not the argument of periapsis. Almanac.jl does the same, so `:brouwer_aop` then
  carries an error of about twice the short-period correction of that angle (0.6° for
  an orbit of 7200 km). This behaviour of ANISE is reproduced, not corrected.
- **Hyperbolic `two_body`.** After a step, ANISE recovers the mean anomaly as the mean
  longitude minus `atan2(h, k)`. For a hyperbola, whose mean anomaly has no period,
  this is wrong by a multiple of 2π when `raan + aop` lies outside (-π, π].
  `two_body` advances the mean anomaly of a hyperbola directly, so there its states
  differ from those of ANISE and agree with `prop2b`. Elliptical orbits use the formula
  of ANISE.
- **Frame rate.** ANISE takes the rate of a local frame from the difference of two
  frames 1 ms apart. `local_frame_state` uses the exact rate of two-body motion, so it
  differs from ANISE by the error of that difference quotient.
- **No NaN values.** Where ANISE returns NaN for an undefined quantity, the calls
  throw `SpiceError(:degeneratecase)`. (A derivative can still be NaN; see
  [Dual numbers](@ref "Dual numbers").) The few values that are returned instead (not finite,
  or a fallback angle) are listed in [`element`](@ref) and [Element names](@ref).
