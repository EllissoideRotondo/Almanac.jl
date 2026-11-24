# Validation of the orbit calls

The calls of [Orbits and relative motion](@ref) are not CSPICE routines. They follow
the source of ANISE main (commit `b3aaadc`) function by function, including its error
and fallback cases. The released ANISE 1.0.0b0 differs from that source in the flight
path angle, inbound hyperbolas, the `TNW` frame and v-infinity.

Where CSPICE computes the same quantity, CSPICE is the reference for the numbers;
otherwise ANISE is. Differences below are relative unless a unit is given.

## What is checked

- **Against CSPICE** (`oscltx`, `recgeo`, `conics`, `prop2b`): elements, geodetic
  coordinates, [`keplerian_orbit`](@ref) and [`two_body`](@ref) agree to a few
  multiples of the rounding error that the sampled orbit amplifies (through `1 - e`
  near parabolic orbits, or through the arc cosine of a nearly circular or equatorial
  one).
- **Against published values:** the GMAT, Vallado and Davis cases and the Brouwer,
  `with_element` and `two_body` cases of the ANISE tests pass within the limits of
  those tests. The largest difference is the B-plane of the GMAT hyperbolic case,
  0.65 km (ANISE allows 1 km).
- **Against ANISE itself:** every element, local frame and relative state of random
  elliptic, near-circular, near-equatorial, retrograde and hyperbolic states, with the
  same inputs given to both libraries. Most elements agree to 1e-11, angles of
  near-circular orbits to 3e-8, and frame rates to 1e-5, the error of the difference
  quotient that ANISE uses (see below).
- **Identities and derivatives:** element identities, orthonormal local frames, the
  exact frame rate and the derivatives through the Kepler solver against ForwardDiff
  and closed forms; Enzyme, Mooncake, Zygote, ReverseDiff and Tracker against
  ForwardDiff.
- **Allocation:** AllocCheck.jl proves that [`element`](@ref), [`keplerian`](@ref),
  [`with_element`](@ref), `two_body`, [`local_frame`](@ref),
  [`local_frame_state`](@ref) and [`relative_state`](@ref) cannot allocate for
  `Float64` inputs.

The comparison with ANISE is not part of the default test run.

## Deviations from ANISE

- **Retrograde Brouwer elements** (an ANISE defect, reproduced). Above 175° of
  inclination ANISE computes on a pseudo-state, the prograde orbit with the inclination
  π - `inc` and the signs of `raan` and `aop` changed. It maps back the inclination and
  `raan` but not the argument of periapsis, so `:brouwer_aop` carries an error of about twice its
  short-period correction (0.6° for a 7200 km orbit). Brouwer elements are those of
  the Earth; another center throws `SpiceError(:degeneratecase)`.
- **Hyperbolic `two_body`.** ANISE recovers the mean anomaly as the mean longitude
  minus `atan2(h, k)`, which is off by a multiple of 2π for a hyperbola when
  `raan + aop` lies outside (-π, π]. `two_body` advances the mean anomaly of a
  hyperbola directly and agrees with `prop2b`.
- **Frame rate.** ANISE differences two frames 1 ms apart; `local_frame_state` uses the
  exact rate of two-body motion.
- **No NaN values.** Where ANISE returns NaN, the calls throw
  `SpiceError(:degeneratecase)`, also for a dual-number state whose value has no
  derivative ([Differentiate an orbit](@ref)). [Element names](@ref) lists the few
  non-finite and fallback values that are returned. The rule holds for values only:
  the derivative of a fallback value, such as the infinite `:ta_dot` of a hyperbola,
  can be NaN.
