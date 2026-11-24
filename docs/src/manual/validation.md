# Validation

CSPICE is the reference. The tests call each routine and its CSPICE counterpart
(through SPICE.jl) on the same input and compare the results, on NAIF kernels and on
kernels written by the CSPICE writers so that every data type is exercised. The
measured differences are in [Known limits of accuracy](@ref).

- **Files:** a kernel written by Almanac.jl holds the data, word for word, of the one
  CSPICE writes from the same input. The bytes differ only where CSPICE leaves them
  [unset](@ref "Fixes that change only the bytes of a file").
- **Derivatives:** the first derivative of a state is compared with the velocity of
  CSPICE, and each higher derivative with automatic differentiation of the one below.
- **Errors:** an input that CSPICE refuses gives the kind of the CSPICE error. Checks
  that CSPICE does not make give kinds of Almanac.jl, such as `:badarg`.
- **Threads:** answers from many threads, while kernels are loaded and unloaded, equal
  those of one thread.
- **Allocation:** AllocCheck.jl proves that the vector, matrix, coordinate, table and
  orbit calls cannot allocate. The tests measure the state and rotation calls of
  [Allocation](@ref) call by call.
- **Package hygiene:** Aqua.jl and ExplicitImports.jl.

The orbit calls have no CSPICE counterpart; see [Validation of the orbit calls](@ref).

## NAIF's test suite

The calls of NAIF's test suite of CSPICE N0067 (tspice_c, the SPICELIB test suite
translated to C) are recorded and replayed against Almanac.jl, and a replay fails on a
difference that is not in the baseline, the list of known differences from CSPICE.
