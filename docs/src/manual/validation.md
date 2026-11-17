# Validation

## How the results were checked

CSPICE itself is the reference. The tests call each routine and its CSPICE counterpart
(through SPICE.jl) on the same input and compare the two results, rather than
comparing with stored numbers. The inputs are NAIF kernels (planetary ephemerides,
lunar orientation, leap seconds, constants) and kernels that the tests write with the
CSPICE writers, so that every data type is exercised.

- **Files.** A kernel written by Almanac.jl is compared with the one CSPICE writes from
  the same input: SPK, CK and binary PCK files byte for byte, DSK and EK files word for
  word.
- **Derivatives.** The first derivative of a geometric state is compared with the
  velocity of CSPICE, and each higher derivative with the automatic differentiation of
  the one below it.
- **Errors.** The kind of each error is compared, and the text of the message where a
  program can depend on it.
- **Concurrency.** Answers from many threads, while kernels are loaded and unloaded,
  are compared with the answers of one thread.
- **NAIF's test suite.** The calls of tspice, the test suite of SPICELIB (400
  families), are recorded through [PEPPER](https://github.com/jacobwilliams/pepper)
  and replayed against Almanac.jl. A difference from SPICELIB is checked against
  CSPICE. This check is not part of the default run; see [Run NAIF's test
  suite](@ref).
- **Allocation.** AllocCheck.jl proves that 68 numeric routines cannot allocate:
  vector, matrix, rotation, coordinate and ellipsoid calls, conic elements, `unitim`,
  `deltet`, the table queries, and the seven calls of the orbit layer listed in
  [Validation of the orbit calls](@ref)
  ([list](https://github.com/Orbitology/Almanac.jl/blob/main/test/test_alloccheck.jl)).
- **Package code.** Aqua.jl finds no method ambiguity, unbound type parameter, type
  piracy, stale dependency or dependency without compat bounds. ExplicitImports.jl
  finds every imported name explicit and public, except four that Base keeps internal
  ([list](https://github.com/Orbitology/Almanac.jl/blob/main/test/test_static.jl)).

### Size of the test suite

- **Default run:** 58,254 tests passed on 2026-11-16 at commit 6fb27bd. The group
  `analysis_orbit` has grown since, from 16,470 to 18,633 tests, so the run now has
  about 60,400.
- **Exhaustive run** (`ALMANAC_TEST_FULL=1`, with larger samples, longer search spans
  and larger files): not measured since the orbit calls were added. The group
  `analysis_orbit` alone passed 371,971 tests at commit 1665511.
- **Line coverage of `src/`:** 96.6 % on 2026-11-16. The coverage job of CI measures
  it again on each push. The group `analysis_orbit` alone covers 98.8 % of `src/analysis/`.
- **Julia versions:** 1.10 and 1.13.

## Run NAIF's test suite

```sh
julia --project=test test/tspice/run.jl
```

The first run builds PEPPER, which needs gfortran 16 and fpm on `PATH`, and records the
calls of about 400 tspice families (about 2 GB, 15 to 40 minutes). Later runs only
replay the recording. A run fails if a mismatch appears that is not in the baseline of
known mismatches. The requirements, the options and the rules of the baseline are in
[`test/tspice/README.md`](https://github.com/Orbitology/Almanac.jl/blob/main/test/tspice/README.md).
