# Validation

Almanac.jl is a port of CSPICE N0067, and CSPICE is its reference: each routine follows
the arithmetic of its CSPICE counterpart operation by operation, including the order of
summation, and uses no fused multiply-add. The tests call each routine and its CSPICE
counterpart, through SPICE.jl, on the same input, with NAIF kernels and with kernels
written by the CSPICE writers so that every data type is exercised. For the states of
`de440s.bsp` with all nine aberration corrections, the states and light times agree bit
for bit in inertial frames and to 4e-15 relative in body-fixed frames, and a kernel
written by Almanac.jl holds the data, word for word, of the one that CSPICE writes from
the same input. NAIF's own test suite of CSPICE N0067 (tspice_c, the SPICELIB test
suite translated to C) is recorded and replayed against Almanac.jl as well, and a replay
fails on any difference that is not in its list of known differences from CSPICE.

## Known limits of accuracy

The elementary functions (`sin`, `cos`, `exp`, ...) are those of Julia, while CSPICE
calls the C library of the platform; the two differ in the last bit for about 3 % of
the arguments of `sin` and `cos`, so a result can differ from CSPICE in its last bit.
The difference grows where a computation carries it over a long span, as in a type 5
SPK segment propagated over six years (3e-8 km at 2e8 km from the Sun), and where the
geometry amplifies it, as at a limb point or a grazing intercept on a DSK shape (up to
2.5e-3 km in NAIF's test suite, where a change of the epoch by one unit in the last
place moves the limb point of CSPICE itself by 3e-3 km). With the `sin` and `cos` of
glibc in place of Julia's, these results equal those of CSPICE.

Three limits hold in CSPICE as well, because they come from the size of a number or from
the method of CSPICE:

- The rotation angle of the Earth in `IAU_EARTH` grows to about 50,000 rad this
  century. Its rounding alone moves the Sun by some meters in that frame.
- An aberration-corrected position has millimeter steps at the boundaries between SPK
  [records](@ref "Smoothness at record boundaries"), because target and observer are
  evaluated at different epochs and cross a boundary at different times.
- `spkezr` and `spkpos` differ by up to 2e-11 relative with the converged corrections
  (`CN`, `CN+S`), because they stop iterating by different rules.

The derivatives of states and rotations differ from the velocities and rates of CSPICE
where CSPICE approximates them; the tutorial section
[Differences from `spkezr` and `sxform`](@ref "Differences from spkezr and sxform") gives each case.

## Other checks

- **Errors:** an input that CSPICE refuses gives the kind of the CSPICE error. Checks
  that CSPICE does not make give kinds of Almanac.jl, such as `:badarg`.
- **Derivatives:** each order of a jet is compared with automatic differentiation of
  the order below, and every AD backend with ForwardDiff.
- **Threads:** answers from many threads, while kernels are loaded and unloaded, equal
  those of one thread.
- **Allocation:** AllocCheck.jl proves that the vector, matrix, coordinate and table
  calls cannot allocate for `Float64` inputs, and the tests measure the
  allocations of [`spkezr`](@ref), [`spkpos`](@ref), [`spkder`](@ref), [`pxform`](@ref),
  [`sxform`](@ref) and [`pxder`](@ref) call by call.
- **Package hygiene:** Aqua.jl and ExplicitImports.jl.
