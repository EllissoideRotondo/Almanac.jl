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
| File | A kernel written by Almanac.jl with the one CSPICE writes: SPK, CK and binary PCK files byte for byte, DSK and EK files word for word |
| Derivative | The first derivative with the velocity of CSPICE, and each higher derivative with automatic differentiation of the one below it |
| Error | The kind of the error, and the text of the message where a program can depend on it |
| Concurrency | Answers from many threads, while kernels are loaded and unloaded, with the answers of one thread |
| NAIF's test suite | The calls of tspice, SPICELIB's test suite (400 families), recorded through [PEPPER](https://github.com/jacobwilliams/pepper) and replayed against Almanac.jl; a difference from SPICELIB is checked against CSPICE ([how to run it](@ref "Run NAIF's test suite"); not in the default run) |
| Allocation | A static proof (AllocCheck.jl) that 61 numeric routines cannot allocate: vector, matrix, rotation, coordinate and ellipsoid calls, conic elements, `unitim`, `deltet` and the table queries ([list](https://github.com/Orbitology/Almanac.jl/blob/main/test/test_alloccheck.jl)) |
| Package code | Aqua.jl: no method ambiguity, unbound type parameter, type piracy, stale dependency or dependency without compat bounds; ExplicitImports.jl: every imported name explicit and public, except four that Base keeps internal ([list](https://github.com/Orbitology/Almanac.jl/blob/main/test/test_static.jl)) |

| Measure | Value |
|---|---|
| Comparisons in the default run | about 40,000 |
| Comparisons in the exhaustive run (`ALMANAC_TEST_FULL=1`), with larger samples, longer search spans and larger files | about 100,000 |
| Line coverage of `src/` | 96.6 % (the coverage job of CI measures it on each push) |
| Julia versions tested | 1.10 and 1.13 |

## Run NAIF's test suite

```sh
julia --project=test test/tspice/run.jl
```

The first run builds PEPPER with gfortran 16 and fpm, which must be on `PATH`, and records the calls of the about 400 tspice families into `~/.cache/almanac-tspice` (about 2 GB, 15 to 40 minutes), with a manifest of the families and of the PEPPER commit. Later runs only replay them, in parallel processes. The recording keeps each distinct call of a routine without state once, and a query once between two changes of the loaded kernels or the settings; `TSPICE_ALL=1` keeps every call.

A run without a list of families is the full check. It stops if a family is missing from the recording, failed to record, or fails to replay. It fails if a mismatch of value, of error kind, or of the script itself is not in the baseline (`test/tspice/baseline_ids.tsv.gz`): each mismatch is identified by family, routine, kind and a hash of its inputs, so a new mismatch fails even when another one went away. It also fails if a routine has more calls skipped or compared with a tolerance than `test/tspice/coverage.tsv` records, or fewer calls compared. The mismatches are in `results/replay.tsv` and `results/ids.tsv`. Each routine of `test/tspice/baseline.tsv` has a class and a note: `documented`, a difference that a docstring or this manual states; `cspice`, Almanac follows CSPICE where SPICELIB differs; `gap`, a call that the replay cannot reproduce, such as a read of a file that the test wrote with routines outside the recording; `unexplained`. `--record` records again, `--update-baseline` makes the full run the new baseline and keeps the classes, and a list of families (`F_SPK08 F_CK06`) restricts the run to them, without the coverage check. `test/tspice/selftest.jl` checks the gate itself: on a copy of the recording it removes a mismatch from the baseline, deletes a log, damages a log and marks a family as failed, and requires the run to fail each time.

## Known limits of accuracy

Apart from the deliberate differences in [Differences from CSPICE](@ref), results differ from CSPICE by rounding. For states of the Moon and the planets from de440s.bsp, in inertial and body-fixed frames, the default run measures:

| Quantity | Aberration corrections | Largest relative difference from CSPICE |
|---|---|---|
| Position, velocity and light time | all nine (`NONE` to `XCN+S`) | 7.3e-16 |

Three limits hold in CSPICE as well, because they come from the size of a number or from the CSPICE method:

- The rotation angle of the Earth in `IAU_EARTH` grows to about 50,000 rad this century; its rounding alone moves the Sun by some meters in that frame.
- An aberration-corrected position has millimeter steps at the boundaries between SPK [records](@ref "Exact derivatives and the rates of CSPICE"), because target and observer are evaluated at different epochs and cross a boundary at different times.
- `spkezr` and `spkpos` differ by up to 2e-11 relative with the converged corrections (`CN`, `CN+S`), because they stop iterating by different rules.
