# Tests and comparisons

CSPICE is the reference. The tests call each routine and its CSPICE counterpart
(through SPICE.jl) on the same input and compare the results, on NAIF kernels and on
kernels written by the CSPICE writers so that every data type is exercised; the
measured differences are in [Known limits of accuracy](@ref). The orbit calls, which
have no CSPICE counterpart, are compared with ANISE as well.

## Comparisons with CSPICE

- **Files:** a kernel written by Almanac.jl holds the data, word for word, of the one
  CSPICE writes from the same input. The bytes differ only where CSPICE leaves them
  unset.
- **Derivatives:** the first derivative of a state is compared with the velocity of
  CSPICE, and each higher derivative with automatic differentiation of the one below.
- **Errors:** an input that CSPICE refuses gives the kind of the CSPICE error. Checks
  that CSPICE does not make give kinds of Almanac.jl, such as `:badarg`.

## NAIF's test suite

The calls of NAIF's test suite of CSPICE N0067 (tspice_c, the SPICELIB test suite
translated to C) are recorded and replayed against Almanac.jl, and a replay fails on a
difference that is not in the baseline, the list of known differences from CSPICE.

## The orbit calls

The calls of [Orbits and relative motion](@ref) follow the source of ANISE function by
function, including its error and fallback cases, with the [differences from
ANISE](@ref "Differences from ANISE") of the API reference. Where CSPICE computes the
same quantity, CSPICE is the reference for the numbers; otherwise ANISE is. Differences
below are relative unless a unit is given.

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
  quotient that ANISE uses. This comparison is not part of the default test run.
- **Identities and derivatives:** element identities, orthonormal local frames, the
  exact frame rate and the derivatives through the Kepler solver against ForwardDiff
  and closed forms; Enzyme, Mooncake, Zygote, ReverseDiff and Tracker against
  ForwardDiff.

## Checks of every call

- **Threads:** answers from many threads, while kernels are loaded and unloaded, equal
  those of one thread.
- **Allocation:** AllocCheck.jl proves that the vector, matrix, coordinate and table
  calls, and the orbit calls [`element`](@ref), [`keplerian`](@ref),
  [`with_element`](@ref), `two_body`, [`local_frame`](@ref),
  [`local_frame_state`](@ref) and [`relative_state`](@ref), cannot allocate for
  `Float64` inputs. The tests measure the allocations of [`spkezr`](@ref), [`spkpos`](@ref),
  [`spkder`](@ref), [`pxform`](@ref), [`sxform`](@ref) and [`pxder`](@ref) call by call.
- **Package hygiene:** Aqua.jl and ExplicitImports.jl.
