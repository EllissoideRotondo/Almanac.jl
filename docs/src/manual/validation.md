# Validation

CSPICE is the reference. The tests call each routine and its CSPICE counterpart
(through SPICE.jl) on the same input and compare the results, on NAIF kernels and on
kernels written by the CSPICE writers so that every data type is exercised. The
measured differences are in [Known limits of accuracy](@ref).

- **Files:** a kernel written by Almanac.jl equals the one CSPICE writes from the same
  input, byte for byte (SPK, CK, PCK) or word for word (DSK, EK).
- **Derivatives:** the first derivative of a state is compared with the velocity of
  CSPICE, and each higher derivative with automatic differentiation of the one below.
- **Errors:** the kind of each error matches CSPICE.
- **Threads:** answers from many threads, while kernels are loaded and unloaded, equal
  those of one thread.
- **Allocation:** AllocCheck.jl proves that the numeric hot paths cannot allocate.
- **Package hygiene:** Aqua.jl and ExplicitImports.jl.

The orbit calls have no CSPICE counterpart; see [Validation of the orbit calls](@ref).

## Run NAIF's test suite

NAIF's test suite of CSPICE N0067 (tspice_c, the SPICELIB test suite translated to C)
can be built with wrappers that record its calls, which are then replayed against
Almanac.jl. This check is not part of the default run:

```sh
julia --project=test test/tspice/run.jl
```

The first run downloads NAIF's archive, builds it with gcc and records the calls, about
2 GB; later runs only replay them. A run fails on a mismatch that is not in the baseline
of known mismatches. The `tspice` workflow of the repository runs this check, in about an
hour, on each push to `main` that changes the package. Details are in
[`test/tspice/README.md`](https://github.com/Orbitology/Almanac.jl/blob/main/test/tspice/README.md).
