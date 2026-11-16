# Analysis layer, part 1: `Orbit`, elements, local frames, and the checks against ANISE

Status: design, for approval. Branch: `analysis` (worktree `../Almanac-analysis`).

## Intent

Almanac.jl is a port of CSPICE N0067. ANISE (nyx-space, Rust and Python) has a set of
higher-level functions that CSPICE lacks: orbital elements of a state, local orbit frames,
eclipse percentages, ground stations, event finding over arbitrary quantities, instruments,
modern Earth frames, trajectory files. The user wants these in Almanac itself, as ANISE is
one package, developed on the `analysis` branch, and checked against ANISE so the results
are known to be right.

The work is about ten independent features. This document designs the first one and the
test machinery that every later feature reuses:

1. the `Orbit` type, its elements, element-based constructors and editing, two-body
   propagation;
2. the local orbit frames (RIC, VNC, RCN, TNW) and relative states;
3. the `anise` test group, which compares Almanac with the ANISE Python package;
4. a derivative test set that runs the new functions through every automatic
   differentiation system Almanac has rules for.

Later features, each with its own design: eclipse percentage, beta angle, Sun angle and
local solar times; ground stations with terrain masks and AER; events and arcs; instruments;
IAU 2000/2006 and TEME Earth frames; trajectories and CCSDS OEM/OPM/STK files; attitude
mathematics (MRP, TRIAD with rates).

## Goals

- Every quantity of ANISE's `Orbit` and `OrbitalElement` available on an Almanac `Orbit`,
  with the same definitions, so that the two agree to rounding on the same input.
- A small API: one type, one element function, property access, a few constructors.
- Every new function generic in the number type, so `ForwardDiff.Dual`, Enzyme, Mooncake,
  Zygote (ChainRules), ReverseDiff and Tracker can differentiate it with respect to the
  epoch and to the state.
- Element calls that do not allocate (Float64 input).

## Non-goals

- Matching ANISE's API names or its units. Almanac keeps km, s and rad; the comparison
  converts.
- ANISE's expression language (`ScalarExpr`, `VectorExpr`, `DcmExpr`) and its S-expression
  serialization. Julia closures replace it (events design, later).
- Numerical propagation, orbit determination: that is Nyx, not ANISE.

## Anti-goals

- A second, ANISE-style naming scheme next to the CSPICE one, or about fifty exported
  element functions (`sma`, `ecc`, `inc`, ...) that clash with user variables.
- Weakening the CSPICE claims of the manual: "differences from CSPICE are rounding" stays
  true of the CSPICE part, and the new part has its own section with its own references.
- Python as a dependency of Almanac or of the default test run.

## Constraints

- Julia 1.10 is the oldest supported version; no new runtime dependencies of Almanac.
- The conventions of the repository: `[ctx,]` first optional argument, `NameOrID` for
  bodies and frames, `SpiceError` with a kind symbol for errors, static arrays for small
  results. The `almanac-docs` and `almanac-tests` skills govern docstrings, manual and
  tests.
- The default test run (`test/parallel.jl`) stays Python-free and keeps its time budget.

## Approaches considered

**Where the code lives.** (a) A separate package on top of Almanac; (b) inside Almanac.
Picked: (b), by the user, as ANISE is one crate. The `add_frame!`/`add_ephemeris!` hooks
would have allowed (a); (b) costs a section of the manual kept apart from the CSPICE parity
claims.

**Shape of the API.** (A) Plain functions on a state vector, CSPICE style:
`sma(state, gm)`, `beta_angle(ctx, state, et, frame, center)`; no new type, but
frame-dependent quantities need many arguments. (B) An `Orbit` type, ANISE style, that
carries the state, epoch, center, frame and the center's constants. Picked: (B), by the
user, on condition that the API stays simple; the sketch below was accepted.

**Element access.** (1) One function with a symbol, `element(o, :sma)`; (2) one exported
function per element. ANISE has both forms: one method per element on `Orbit`
(`sma_km()`), and the `OrbitalElement` enum with `evaluate` for events and reports.
Picked: both ANISE forms in Julia terms, by the user: property access `o.sma` (the Julia
counterpart of a method; exports nothing, cannot clash, tab-completes) and
`element(o, :sma)` (the counterpart of the enum, for code that holds the element as a
value). Both run the same code.

**Checking against ANISE.** (i) A script outside the repository; (ii) a test group inside
it, with its own environment, not in the default run. Picked: (ii), by the user.

## The design

### `Orbit`

```julia
struct Orbit{T<:Real}
    state::SVector{6,T}      # km, km/s, in `frame`, relative to `center`
    et::T                    # TDB seconds past J2000
    center::Int32            # NAIF ID of the central body
    frame::Int32             # frame ID
    gm::Float64              # km^3/s^2 of the center
    radii::SVector{3,Float64}  # km, the center's triaxial radii; NaN when unknown
end
```

Constructors:

```julia
orbit([ctx,] targ, et, ref, abcorr, obs)        # from the kernels; argument order of spkezr
Orbit([ctx,] state, et, center, frame)          # GM and radii from the pool of ctx
Orbit(state, et, center, frame, gm[, radii])    # no kernels at all
keplerian_orbit([ctx,] et, center, frame; sma, ecc, inc, raan, aop, ta)
    # alternatives: periapsis= and apoapsis= in place of sma and ecc; ma= in place of ta
```

GM is `BODYnnn_GM` and the radii `BODYnnn_RADII` of the pool. A missing GM is an error at
construction; missing radii give NaN, and only the elements that need them fail. `T` is
the promoted type of the state and the epoch, so a `Dual` epoch through `orbit` or a
`Dual` state through `Orbit` both work.

The orbit is taken in the frame it is given; Keplerian elements are meaningful in an
inertial frame, geodetic ones in the center's body-fixed frame, as in ANISE. Almanac does
not rotate behind the user's back.

### Elements

`element(o, name::Symbol)` and the property `o.name` return one quantity. The symbols, with
ANISE's names (units are km, s, rad and their combinations):

| Symbols | ANISE |
|---|---|
| `:x :y :z :vx :vy :vz`, `:rmag`, `:vmag`; properties only: `position`, `velocity` | `X`..`VZ`, `Rmag`, `Vmag` |
| `:sma`, `:ecc`, `:inc`, `:raan`, `:aop`, `:ta`, `:ma`, `:ea`, `:ha` | `SemiMajorAxis`, `Eccentricity`, ..., `HyperbolicAnomaly` |
| `:aol`, `:tlong`, `:period`, `:mean_motion`, `:ta_dot` | `AoL`, `TrueLongitude`, `Period`, `mean_motion_deg_s`, `ta_dot_deg_s` |
| `:energy`, `:c3`, `:fpa`, `:semi_parameter`, `:semi_minor_axis` | `Energy`, `C3`, `FlightPathAngle`, ... |
| `:hmag`, `:hx`, `:hy`, `:hz`; properties only: `hvec`, `evec` | `hvec`, `Hmag`, `HX`, `HY`, `HZ`, `evec` |
| `:periapsis`, `:apoapsis`, `:periapsis_altitude`, `:apoapsis_altitude`, `:sma_altitude`, `:altitude` | `PeriapsisRadius`, ..., `Altitude` |
| `:right_ascension`, `:declination`, `:velocity_declination` | same |
| `:latitude`, `:longitude`, `:height` | geodetic `Latitude`, `Longitude`, `Height` |
| `:light_time` | `LightTime` |
| `:bplane_bt`, `:bplane_br`, `:bplane_mag`, `:bplane_angle`, `:bplane_ltof` | `BPlaneBT`, ..., `BPlaneLTOF` |
| `:equinoctial_a`, `:equinoctial_h`, `:equinoctial_k`, `:equinoctial_p`, `:equinoctial_q`, `:equinoctial_lambda` | `equinoctial_a_km`, `EquinoctialH`, ... |
| `:brouwer_sma`, `:brouwer_ecc`, `:brouwer_inc`, `:brouwer_raan`, `:brouwer_aop`, `:brouwer_ma` | `BrouwerMeanShort*` |

ANISE's `vinf_periapsis_km(turn_angle)` and `vinf_turn_angle_deg(periapsis)` take an
argument, so they are two functions, `vinf_periapsis(o, turn_angle)` and
`vinf_turn_angle(o, periapsis)`, not elements (builder's call, settled in review).

`propertynames(o)` lists the fields and every symbol.

Allocation (builder's call, settled in review): property access and `element` with a
constant symbol are type-stable and do not allocate; AllocCheck covers those. With a
symbol known only at run time, `element` returns a scalar for the scalar elements and is
type-stable when the caller asserts `::Float64`; the vector quantities (`:position`,
`:velocity`, `:hvec`, `:evec`) are properties only and not accepted by `element`, so a
run-time symbol never yields a union of scalar and vector. `keplerian(o)` returns
`(; sma, ecc, inc, raan, aop, ta)`. `elements()` returns the list of symbols, for the
manual and for loops.

Definitions follow ANISE's source exactly, including its edge cases, with these points
recorded because they are not obvious:

- Degenerate cases (builder's call, settled in review): Almanac does what the reference
  ANISE does, case by case. Where ANISE fails (for example `:ta` of a circular orbit,
  `ecc < 1e-11`; `:aol` of an equatorial one, `|ẑ×ĥ| < 1e-6·π/180`; the B-plane of an
  elliptic orbit; the anomalies of a parabolic one), Almanac throws
  `SpiceError(:degeneratecase)` with the same threshold. Where ANISE returns a fallback
  value (for example `:aop` and `:raan`, which turn the NaN of `acos` into 0 or π), Almanac
  returns the same value. The tests list each case with ANISE's behavior.
- `:altitude` uses the mean equatorial radius (mean of the first two radii); `:height` is
  geodetic, on the ellipsoid of the center with ANISE's closed-form algorithm.
- The Brouwer mean elements use the Earth's J2 (`1.082626925638815e-3`) whatever the
  center, as ANISE does, and are documented as valid for the Earth only; ANISE's validity
  test (`inc <= 180` deg, `ecc < 1`, periapsis >= 3000 km) becomes a `SpiceError`.

### Building and editing

```julia
with_element(o, :sma, 7000.0)          # a new Orbit, other elements kept
with_element(o, :inc, deg2rad(98.0))   # :sma :ecc :inc :raan :aop :ta, and
with_apsides(o, periapsis, apoapsis)   # both apsides at once
two_body(o, et2)                       # Kepler propagation to et2, ANISE's at_epoch
```

`with_element` covers ANISE's `with_*`, `set_*` and `add_*` (an `add` is
`with_element(o, :sma, o.sma + d)`). `two_body` follows ANISE's equinoctial method, so it
behaves for circular equatorial orbits.

### Local frames

```julia
local_frame(o, :ric)        # SMatrix{3,3}: local to the orbit's frame; also :vnc :rcn :tnw
local_frame_state(o, :ric)  # SMatrix{6,6}: the same for states, with its exact time derivative
relative_state(o, chief, :ric)   # SVector{6}: o minus chief, in chief's local frame
```

Axes as in ANISE: RIC = (r̂, ĥ×r̂, ĥ); VNC = (v̂, ĥ, v̂×ĥ); TNW = (v̂, ĥ×v̂, ĥ); RCN is the
same matrix as RIC in ANISE and is kept so. ANISE obtains the rate of the 6×6 matrix by
finite differences of two two-body states 1 ms apart; Almanac computes it exactly, so the
comparison of rates has a looser tolerance, stated in the test. `relative_state` requires
both orbits at the same epoch, center and frame (`SpiceError(:badarg)` otherwise), as
`ric_difference` does.

### Exports and manual

New exports: `Orbit`, `orbit`, `element`, `elements`, `keplerian`, `keplerian_orbit`,
`with_element`, `with_apsides`, `two_body`, `vinf_periapsis`, `vinf_turn_angle`,
`local_frame`, `local_frame_state`,
`relative_state`. Code in `src/analysis/`. Manual: a reference page "Orbits" and a how-to,
under a new "Beyond CSPICE" part of the navigation; `coverage.md` gets a separate section
"Beyond CSPICE" with its own validation table (ANISE, literature, derivatives).

## Checks

### Default test run (inside the repository, no Python)

A test group `analysis_orbit`:

- agreement with CSPICE where an equivalent exists: `oscltx` (sma, ecc, inc, raan, aop,
  ma, ta, period), `recgeo` (latitude, longitude, height), `prop2b` (two-body), `eqncpv`
  (equinoctial);
- published values: Vallado's worked examples for elements and for Kepler propagation;
- identities: element round trips (`keplerian_orbit(keplerian(o))`), orthonormal local
  frames, `relative_state(o, o) == 0`, the time derivative of `local_frame_state` against
  ForwardDiff;
- each degenerate case raises `SpiceError(:degeneratecase)`;
- AllocCheck on the element calls;
- regression values: a small table of ANISE results produced by the `anise` group, pinned
  as numbers.

A test group `analysis_derivatives`: every new function through ForwardDiff, Enzyme,
Mooncake, Zygote, ReverseDiff and Tracker by way of DifferentiationInterface.jl. The rules
for the last five live in extensions now being committed to `main` (the user has asked for
it); `main` is merged into `analysis` once they land, and the group is built after that
merge. It differentiates with
respect to the epoch (through `orbit`) and to the state; each result against ForwardDiff
and a central finite difference. The test environment gains those packages; a backend
that does not install or fails on a Julia version is skipped with a message, not silently.

### The `anise` group (inside the repository, opt-in)

`test/anise/` with its own `Project.toml` (PythonCall, CondaPkg, Almanac by path) and a
`CondaPkg.toml` that installs `anise` from PyPI. It is not in `GROUPS` of `runtests.jl`;
it runs with `julia --project=test/anise test/anise/runtests.jl`, like the tspice replay.

The reference is ANISE's main source, read at commit `b3aaadc` (1.0.0-beta, 2026-11-15),
decided by the user. The released versions differ from it in definitions: 0.10.6 has a
flight-path-angle precedence bug, a different hyperbolic mean anomaly, rejects inbound
hyperbolic states and has no TNW; 1.0.0b0 (2026-10-19) still lacks TNW and the hyperbolic
fixes. Almanac's code follows main now. The group is written now and pinned to the first
PyPI release that contains main's definitions (1.0.0 or a later beta); until that release
exists the group is not run, and the in-repository checks (CSPICE, literature, identities)
carry the feature. Building ANISE from git is not used.

For a few thousand random states (elliptic, near-circular, near-equatorial, hyperbolic,
retrograde) around the Earth and Mars, both libraries get the same state, epoch and GM
(and radii, given to ANISE's `Frame` explicitly, so no ANISE planetary-constant file is
needed). The group compares every element and every local frame and reports the largest
relative difference of each quantity. Each comparison has a tolerance; failures name the
state.

## Left to the builder

- File layout inside `src/analysis/`, internal helpers, the exact docstrings.
- The choice of random states and the tolerance of each quantity, from the measured
  differences, stated in the test.
- Whether `orbit` gets a `Dual` method directly or through the existing state rules
  (after the merge of the derivative rules from `main`).
- How `CondaPkg` is pointed at a pre-installed environment on a machine without internet.
