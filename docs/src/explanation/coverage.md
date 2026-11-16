# Coverage and validation

## What is available

Almanac.jl covers the library of CSPICE N0067; [What is absent](@ref) lists the routines replaced by Julia facilities.

| Area | Available |
|---|---|
| Kernel loading | SPK, CK, binary PCK, DSK, EK, text kernels (leap seconds, PCK, frames, spacecraft clock), meta-kernels |
| SPK data types | all types of CSPICE N0067: 1, 2, 3, 5, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21; `getelm`, `evsgp4` |
| States | `spkezr`, `spkpos`, `spkgeo`, `spkssb`, the integer-ID forms `spkez`, `spkezp` and `spkgps`, `spkltc`, `spkaps`, `spkacs`, `spkapo`, with all aberration corrections of CSPICE; constant-position and constant-velocity observers and targets (`spkcpo`, `spkcpt`, `spkcvo`, `spkcvt`) |
| Frames | inertial, body-fixed (text and binary PCK), fixed-offset, CK-based, dynamic (all families) and switch frames; `pxform`, `sxform`; `frinfo`, `cidfrm`, `cnmfrm`, `ccifrm`, `tkfram`, `tipbod`, `tisbod`, `pckfrm`, `pckcov`, `bltfrm`, `kplfrm` |
| Derivatives | `spkder` and `pxder` to any order; a `ForwardDiff.Dual` epoch in the calls that [Dual-number epochs](@ref) lists |
| Orientation and clocks | CK types 1 to 6, `ckgp`, `ckgpav`, `ckcov`, `ckobj`, `ckfrot`, `ckfxfm`; records of CK types 2 and 3 (`cknr02`, `ckgr02`, `cknr03`, `ckgr03`); spacecraft clock type 1 |
| Time | `str2et`, `utc2et`, `timout`, `et2utc`, `etcal`, `unitim`, `tparse`, `tpictr`, `timdef`, `tsetyr`, `tparch`, `deltet`, `et2lst`, `lspcn` |
| File access and writers | DAF read and write calls, SPK writers for all 16 types (type 14 through `spk14b`, `spk14a`, `spk14e`), CK writers for 5 types, binary PCK writers of types 2 and 20 (`pckopn`, `pckw02`, `pckw20`, `pckcls`), comment area, windows, `spkobj`, `spkcov`, `kinfo`, `kdata`; SPK segments by descriptor (`spksfs`, `spkpvn`, `spkuds`, and `spksub` for every type) |
| Kernel pool | `gdpool`, `gcpool`, `gipool`, `bodvrd`, `bodvcd`, `bodfnd`, body and frame name calls; `pdpool`, `pcpool`, `pipool`, `dvpool`, `lmpool`, `ldpool`, `clpool`, `boddef`; `expool`, `dtpool`, `gnpool`, `szpool`, `stpool`, `badkpv`; watchers `swpool`, `cvpool` for each context |
| Mathematics | vectors, matrices, rotations and quaternions, coordinate conversions with their Jacobians, ellipsoid, plane and ellipse calls, conic elements and two-body propagation, equinoctial elements (`eqncpv`), Chebyshev, Lagrange and Hermite interpolation, base 16 strings of doubles (`dp2hx`, `hx2dp`) |
| Shape kernels | DAS file access, DSK type 2 read and write, plate access calls, ray-surface intercept (`dskxv`, `dskxsi`, `dskx02`), `latsrf`, `srfnrm`, surface name calls (type 2 is the only DSK type of CSPICE N0067) |
| Event kernels | EK read access (`ekopr`, `ekssum`, `ekrcec`, `ekrced`, `ekrcei`), load and unload of one file (`eklef`, `ekuef`), table catalog (`ekntab`, `ektnam`, `ekccnt`, `ekcii`), queries (`ekpsel`, `ekfind`, `eknelt`, `ekgc`, `ekgd`, `ekgi`) with a result object for each query, writers by record (`ekopn`, `ekops`, `ekbseg`, `ekappr`, `ekinsr`, `ekdelr`, `ekace*`, `ekuce*`) and the fast writers, which write whole columns (`ekifld`, `ekacl*`, `ekffld`) |
| Geometry finder | `gfdist`, `gfsep`, `gfposc`, `gfrr`, `gfpa`, `gfsubc`, `gfsntc`, `gfilum`, `gfoclt`, `gftfov`, `gfrfov`, `gfuds`, `gfudb`, `gfevnt` for the six quantities of CSPICE, `gfocce`, `gffove`, `gfstol`, `gfrefn`, `uddf`, `uddc` |
| Geometry | `subpnt`, `subslr`, `sincpt`, `ilumin`, `illumg`, `illumf`, `phaseq`, `occult`, `limbpt`, `termpt`, `tangpt`, `getfov`, `getfvn`, `fovtrg`, `fovray`, `azlcpo`, `azlrec`, `recazl`, `dazldr`, `drdazl`, `srfrec`, on ellipsoids and on DSK data |
| Ephemeris and rotation tables | `ephemeris_table`, `table_position`, `table_state`, `table_jet`: the position of one body relative to another on a fixed interval, as one array of coefficients for the CPU or a GPU; `table_to_spk`: the table as an SPK segment of type 2; `rotation_table`, `table_rotation`, `table_rotation_jet`: the same for the rotation between two frames |

## How the results were checked

CSPICE itself is the reference: the tests call each routine and its CSPICE counterpart (through SPICE.jl) on the same input and compare, rather than comparing with stored numbers. Inputs are NAIF kernels (planetary ephemerides, lunar orientation, leap seconds, constants) and kernels the tests write with the CSPICE writers, so that every data type is exercised.

| Kind of check | What it compares |
|---|---|
| Differential | A result with the CSPICE result of the same call |
| File | A kernel written by Almanac.jl with the one CSPICE writes, byte for byte |
| Derivative | The first derivative with the velocity of CSPICE, and each higher derivative with automatic differentiation of the one below it |
| Error | The kind of the error, and the text of the message where a program can depend on it |
| Concurrency | Answers from many threads, while kernels are loaded and unloaded, with the answers of one thread |
| NAIF's test suite | The calls of tspice, SPICELIB's test suite (400 families), recorded through [PEPPER](https://github.com/jacobwilliams/pepper) and replayed against Almanac.jl; a difference from SPICELIB is checked against CSPICE (`test/tspice/`, not in the default run) |
| Allocation | A static proof (AllocCheck.jl) that 61 numeric routines cannot allocate: vector, matrix, rotation, coordinate and ellipsoid calls, conic elements, `unitim`, `deltet` and the table queries ([list](https://github.com/Orbitology/Almanac.jl/blob/main/test/test_alloccheck.jl)) |

| Measure | Value |
|---|---|
| Comparisons in the default run | about 40,000 |
| Comparisons in the exhaustive run (`ALMANAC_TEST_FULL=1`), with larger samples, longer search spans and larger files | about 100,000 |
| Line coverage of `src/` | 96.6 % |
| Julia versions tested | 1.10 and 1.13 |

## Known limits of accuracy

Apart from the deliberate differences in [Differences from CSPICE](@ref), results differ from CSPICE by rounding. For states of the Moon and the planets from de440s.bsp, in inertial and body-fixed frames, the default run measures:

| Quantity | Aberration corrections | Largest relative difference from CSPICE |
|---|---|---|
| Position, velocity and light time | all nine (`NONE` to `XCN+S`) | 7.3e-16 |

Three limits hold in CSPICE as well, because they come from the size of a number or from the CSPICE method:

- The rotation angle of the Earth in `IAU_EARTH` grows to about 50,000 rad this century; its rounding alone moves the Sun by some meters in that frame.
- An aberration-corrected position has millimeter steps at the boundaries between SPK [records](@ref "Exact derivatives and the rates of CSPICE"), because target and observer are evaluated at different epochs and cross a boundary at different times.
- `spkezr` and `spkpos` differ by up to 2e-11 relative with the converged corrections (`CN`, `CN+S`), because they stop iterating by different rules.

## Beyond CSPICE

The calls of [Orbits](@ref) (`Orbit`, `element`, `keplerian_orbit`, `two_body`, the local frames and the rest) have no CSPICE counterpart, so the claims above about agreement with CSPICE and about rounding do not cover them. They follow the source of ANISE (main, commit `b3aaadc`, 1.0.0-beta of 2026-11-15) function by function, including its error and fallback cases. A released ANISE differs from that commit in some definitions. The calls are checked in the test group `analysis_orbit`, with the results below measured in the default run; where CSPICE has an equivalent routine, it is the reference.

| Check | What it compares | Measured result |
|---|---|---|
| GMAT cases | The Keplerian elements of three Cartesian and Keplerian state pairs (circular inclined, elliptical, circular equatorial), the equatorial and the two hyperbolic cases, and 24 flight path angles, with the values GMAT gave to the tests of ANISE; angles modulo 2π | largest scaled difference 4.0e-13 (limit 1e-12); flight path angles 2.2e-16 |
| `oscltx` | `sma`, `ecc`, `inc`, `raan`, `aop`, `ma`, `ta` and `period` of 200 random orbits with the output of `oscltx` | relative 9.6e-15 (periapsis), 2.7e-15 (`ecc`), 1.3e-15 (`inc`), 4.0e-13 (`raan`), 1.1e-14 (`aop`), 2.2e-12 (`ta`), 1.1e-11 (`ma`), 2.8e-15 (`sma`, `period`) |
| `recgeo` | `latitude`, `longitude` and `height` at 900 random points about three bodies, one of them flattened by 0.1, with `recgeo` | latitude 3.9e-16, longitude 4.4e-16, height 1.1e-15 (in radians and equatorial radii) |
| Vallado | Longitude, geodetic latitude, height and the altitudes of Examples 3-3 and 3-4 of Vallado, with the radii that ANISE uses | 2.2e-16 rad (latitude), 1.1e-16 rad (longitude), 2.1e-12 km (height), 4.2e-12 km (altitudes) |
| GMAT B-plane | `bplane_bt`, `bplane_br`, `bplane_mag` and `bplane_angle` of one hyperbolic case with the values of GMAT | 0.117 km, 0.649 km, 0.230 km and 2.0e-5 rad (the limits of ANISE are 1 km and 0.01°; the cause of the 0.65 km is not isolated) |
| Davis | The hyperbolic anomaly of the state of the lecture of Davis, as ANISE tests it, and `vinf_periapsis` as the inverse of `vinf_turn_angle` | hyperbolic anomaly 1.1e-11 rad; inverse 1.1e-13 km |
| `conics` | `keplerian_orbit` against the state of `conics` for 200 random orbits | position 1.2e-12, velocity 2.5e-12 (relative) |
| `prop2b` | `two_body` over ±10 periods against `prop2b` | position 1.3e-12, velocity 9.2e-13, `sma` 5.8e-15 (relative) |
| ANISE test cases | The round trips of Brouwer mean elements, `with_element` and `with_apsides` on the loops of `verif_with_init`, and `two_body` on `verif_orbit_at_epoch`, as ANISE tests them | Brouwer: `sma` 3.1e-11 km and angles up to 3e-11 rad for a low orbit, 2.4e-9 km for a transfer orbit; `with_element` and `with_apsides` at most 1.5e-11; `verif_orbit_at_epoch` `aop` 3.3e-8 rad (the limit of ANISE is 5.2e-8) |
| Identities | `aol = aop + ta`, `c3 = 2 energy`, `mean_motion`, the hyperbolic functions of `ha`, the B-plane and the v-infinity formulas, the conversion back from equinoctial elements, `relative_state(o, o) == 0`, orthonormal local frames, the inverse of `local_frame_state` | at most 2.0e-14, except `aol` 2.2e-12; equinoctial conversion 3.1e-15; orthonormality 5.6e-16 |
| Derivatives | The exact rate of `local_frame_state` against `ForwardDiff.derivative` of the frame of `two_body`; the derivative of the state of `two_body` with respect to the epoch, a `Dual` epoch, against the velocity and the acceleration of two-body motion | 1.2e-15 1/s for the frame rate; 2.9e-12 relative for the state |
| Allocation | A static proof (AllocCheck.jl) that `element`, `keplerian`, `with_element`, `two_body`, `local_frame`, `local_frame_state` and `relative_state` cannot allocate ([list](https://github.com/Orbitology/Almanac.jl/blob/main/test/test_alloccheck.jl)) | no possible allocation in any of the seven |
| ANISE | Every element and local frame of a few thousand random states about the Earth and Mars with the same state, epoch and GM | planned (group `anise`, opt-in); not yet written, so not yet run |

Four points are deliberate and are not defects:

- **Brouwer elements.** They are those of the Earth. As in ANISE, a center other than body 399 is an error. For an inclination above 175° the argument of periapsis carries an error of ANISE's pseudo-state, about twice the short-period correction of that angle (0.6° for an orbit of 7200 km); it is reproduced, not corrected. The round trip of the retrograde test case differs by 1.09e-2 rad in `aop` for that reason.
- **Hyperbolic `two_body`.** ANISE recovers the mean anomaly after a step as the mean longitude minus `atan2(h, k)`; for a hyperbola, whose mean anomaly has no period, that is wrong by a multiple of 2π when `raan + aop` lies outside (-π, π]. `two_body` advances the mean anomaly of a hyperbola directly, so its states can differ from those of ANISE there; elliptical orbits use the formula of ANISE.
- **Frame rate.** ANISE takes the rate of a local frame from the difference of two frames 1 ms apart. `local_frame_state` uses the exact rate of two-body motion, so it differs by the error of that difference quotient.
- **No NaN.** Where ANISE returns NaN for an undefined quantity, the calls throw `SpiceError(:degeneratecase)`.

The test group is in the default run, which has a budget of 120 s for all groups; `analysis_orbit` takes 16 to 30 s of it.
