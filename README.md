# Almanac.jl

[![CI](https://github.com/Orbitology/Almanac.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Orbitology/Almanac.jl/actions/workflows/CI.yml)
[![Manual](https://img.shields.io/badge/docs-dev-blue.svg)](https://orbitology.github.io/Almanac.jl/dev/)

Almanac.jl is a native Julia port of [CSPICE](https://naif.jpl.nasa.gov/naif/toolkit.html)
N0067, NASA's toolkit for spacecraft and planetary geometry. It reads SPICE kernels and
computes positions, orientations, times and events in pure Julia, with no C library behind
it. The calls keep the CSPICE names and, with few exceptions, the argument order, so code
written for CSPICE or [SPICE.jl](https://github.com/JuliaAstro/SPICE.jl) ports with few
changes.

## Installation

Requires Julia 1.10 or later.

```julia-repl
pkg> add https://github.com/Orbitology/Almanac.jl
```

## Example

With the NAIF generic kernels
[`naif0012.tls`](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls),
[`de440s.bsp`](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp) and
[`pck00011.tpc`](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc)
in the working directory:

```julia
using Almanac

furnsh("naif0012.tls")    # leap seconds
furnsh("de440s.bsp")      # planets, Moon and Sun
furnsh("pck00011.tpc")    # orientation and radii of the bodies

et = str2et("2025 MAR 01 12:00:00 UTC")              # seconds past J2000 TDB
state, lt = spkezr("MOON", et, "J2000", "LT+S", "EARTH")
state[1:3]                                           # position, km
```

[SPICEKernels.jl](https://github.com/JuliaAstro/EphemerisSources.jl/tree/main/lib/SPICEKernels)
downloads these files on first use; after `using SPICEKernels`, write
`furnsh(naif0012_tls())`, `furnsh(de440s())` and `furnsh(pck00011())`.

## Features

- CSPICE calls with Julia conventions: results are returned, not written to output
  arguments; indices start at 1; errors throw a `SpiceError`.
- Lock-free queries from any number of threads, also while another thread loads kernels.
- Independent kernel sets (`SpiceContext`), for example to compare two ephemerides.
- Time derivatives of any order (`spkder`, `pxder`) and ForwardDiff.jl dual-number epochs.
- Chebyshev tables of a trajectory or a rotation for inner loops and GPU kernels.
- Event searches, surface geometry (ellipsoids and DSK), and writers for SPK, CK, binary
  PCK, DSK and EK files.

Results match CSPICE to rounding, checked against CSPICE through SPICE.jl and against
NAIF's own test suite. The interface may still change.

## Documentation

The [manual](https://orbitology.github.io/Almanac.jl/dev/) has a tutorial, how-to guides,
the [differences from CSPICE](https://orbitology.github.io/Almanac.jl/dev/explanation/cspice/)
and the API reference. Each docstring links to the NAIF page of the routine it ports.

## Relationship to NAIF

Almanac.jl is an independent project, not developed, supported or endorsed by NAIF or JPL.
Report defects here, not to NAIF. The repository does not redistribute CSPICE or any kernel.

Publications that use SPICE should cite:

- C. H. Acton, "Ancillary data services of NASA's Navigation and Ancillary Information
  Facility", *Planetary and Space Science* 44 (1), 65–70, 1996.
- C. Acton, N. Bachman, B. Semenov, and E. Wright, "A look towards the future in the
  handling of space science mission geometry", *Planetary and Space Science* 150, 9–12,
  2018.

## Contributing

Questions and defect reports go to the [issues](https://github.com/Orbitology/Almanac.jl/issues).

```sh
julia --startup-file=no test/parallel.jl                     # tests
julia --project=docs -e 'using Pkg; Pkg.instantiate()'      # once
julia --project=docs docs/make.jl                            # manual, in docs/build/
```

## License

[MIT](LICENSE). The license covers only this repository and grants no rights over CSPICE.
