# Almanac.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://orbitology.github.io/Almanac.jl/dev/)
[![CI](https://github.com/Orbitology/Almanac.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Orbitology/Almanac.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/Orbitology/Almanac.jl/graph/badge.svg)](https://codecov.io/gh/Orbitology/Almanac.jl)

*NASA NAIF's SPICE toolkit in native Julia*

**Almanac.jl** is a pure-Julia port of [CSPICE](https://naif.jpl.nasa.gov/naif/toolkit.html)
(toolkit version N0067), the library of NASA's Navigation and Ancillary Information Facility
(NAIF) for spacecraft and planetary geometry. It reads SPICE kernels and computes positions,
velocities, orientations, time conversions, surface geometry and event times.

The CSPICE routines keep their names and input arguments, so NAIF's
[documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/index.html) applies to each
of them, and code written for CSPICE or [SPICE.jl](https://github.com/JuliaAstro/SPICE.jl)
carries over with few changes. Because it is plain Julia, Almanac.jl also offers lock-free
queries from many threads, automatic differentiation with respect to time, and ephemeris
tables that run on a GPU.

**Note:** This project is *not* affiliated with NASA, NAIF or JPL in any way.

**AI assistance:** the code for this project was 
written with the help of AI coding assistants (Claude Opus 5.5). 
The code was [validated](https://orbitology.github.io/Almanac.jl/dev/manual/validation/) against CSPICE using NAIF's test suite, [tspice](https://naif.jpl.nasa.gov/pub/naif/misc/tspice/N0067/). 

## Installation

```julia
import Pkg; Pkg.add("Almanac")
```

## Quickstart

```julia
using Almanac
using SPICEKernels

# Load the leap seconds kernel and a planetary ephemeris
furnsh(naif0012_tls())
furnsh(de440s())

# Convert the calendar date to ephemeris seconds past J2000
et = utc2et("2025-03-01T12:00:00")

# Position of the Moon as seen from the Earth, corrected for light time and stellar aberration
pos, lt = spkpos("MOON", et, "J2000", "LT+S", "EARTH")
```

## Features

- **The CSPICE routines.** States, frames, time systems, spacecraft clocks, the kernel pool
  and the mathematics routines, with every SPK and CK data type.
- **Geometry and events.** Sub-observer points, illumination angles, occultations and fields
  of view on ellipsoids and DSK shape models, and event searches with CSPICE's geometry finder
  (the `gf` routines).
- **Kernel writers.** SPK, CK, binary PCK, DSK and EK files, with the same contents as the
  files CSPICE writes.

Beyond CSPICE:

- **Thread safety.** Each query reads an immutable snapshot of the loaded kernels, so any
  number of threads can query at once, even while another thread loads kernels.
- **Independent kernel sets.** A `SpiceContext` holds its own kernels and settings, so a
  program can compare two ephemerides or keep two missions apart.
- **Derivatives.** `spkder` and `pxder` return time derivatives of any order, and states,
  rotations and surface points differentiate with respect to the epoch with ForwardDiff.jl,
  Enzyme.jl, Mooncake.jl, Zygote.jl, ReverseDiff.jl and Tracker.jl.
- **Fast tables.** `ephemeris_table` and `rotation_table` compress a trajectory or a rotation
  into Chebyshev coefficients, for inner loops and for code that runs on a GPU.

The test suite compares Almanac.jl with CSPICE through SPICE.jl and replays NAIF's own test
suite against Almanac.jl. Results agree to within rounding error, apart from a few
[documented differences](https://orbitology.github.io/Almanac.jl/dev/manual/validation/)
and the defects of CSPICE that Almanac.jl fixes.
The package is under active development, and its interface may still change.

## Citing

If you publish work that uses Almanac.jl, please cite the SPICE toolkit as NAIF requests:

- C. H. Acton, "Ancillary data services of NASA's Navigation and Ancillary Information
  Facility", *Planetary and Space Science* 44 (1), 65–70, 1996.
- C. Acton, N. Bachman, B. Semenov, and E. Wright, "A look towards the future in the handling
  of space science mission geometry", *Planetary and Space Science* 150, 9–12, 2018.
  
## Documentation

For information on using the package, see the
[documentation](https://orbitology.github.io/Almanac.jl/dev/). It opens with a
quickstart, followed by tutorials, pages on validation, and
the API reference, in which the docstring of each ported routine links to its NAIF page.

## Acknowledgements

Almanac.jl keeps the interface of [SPICE.jl](https://github.com/JuliaAstro/SPICE.jl), which
is also the test oracle. [PEPPER](https://github.com/jacobwilliams/pepper), a modern Fortran
reimagining of SPICELIB, inspired the port; the JSON kernels follow its
[jsonspice](https://github.com/jacobwilliams/jsonspice) format.

## License

Almanac.jl is released under the [MIT license](LICENSE), which covers this repository only.
The repository contains no CSPICE source and no kernels, and the license grants no rights
over CSPICE, which NAIF distributes under its own
[rules of use](https://naif.jpl.nasa.gov/naif/rules.html).
