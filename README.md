# Almanac.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://orbitology.github.io/Almanac.jl/dev/)
[![CI](https://github.com/Orbitology/Almanac.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Orbitology/Almanac.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/Orbitology/Almanac.jl/graph/badge.svg)](https://codecov.io/gh/Orbitology/Almanac.jl)

*NASA NAIF's SPICE toolkit, in native Julia*

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

**AI assistance:** Much of the code, the tests and the documentation of Almanac.jl was
written with AI coding assistants (Claude, by Anthropic), under the direction of the
maintainer. Its correctness rests on the comparisons with CSPICE and with NAIF's test suite
described [below](#features) and in the
[manual](https://orbitology.github.io/Almanac.jl/dev/manual/validation/). For results
that matter to your work, check them against CSPICE, and please report any difference in
the [issues](https://github.com/Orbitology/Almanac.jl/issues).

## Installation

Almanac.jl requires Julia 1.10 or later. It is not yet registered, so add it from its
repository:

```julia
using Pkg
Pkg.add(url="https://github.com/Orbitology/Almanac.jl")
```

## Quickstart

```julia
using Almanac
using SPICEKernels  # downloads the NAIF generic kernels and returns their paths

# Load the leap seconds kernel and a planetary ephemeris
furnsh(naif0012_tls())
furnsh(de440s())

# Convert the calendar date to ephemeris seconds past J2000
et = utc2et("2025-03-01T12:00:00")

# Position of the Moon as seen from the Earth, corrected for light time and stellar aberration
pos, lt = spkpos("MOON", et, "J2000", "LT+S", "EARTH")
```

`pos` is in kilometers and `lt`, the one-way light time, in seconds. As in SPICE.jl, outputs
are returned rather than written into arguments, and indices start at 1; errors throw a
`SpiceError`.

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
- **Independent kernel sets.** A `SpiceContext` holds its own kernels and settings, so you
  can compare two ephemerides or keep two missions apart.
- **Derivatives.** `spkder` and `pxder` return time derivatives of any order, and states,
  rotations and surface points differentiate with respect to the epoch with ForwardDiff.jl,
  Enzyme.jl, Mooncake.jl, Zygote.jl, ReverseDiff.jl and Tracker.jl.
- **Fast tables.** `ephemeris_table` and `rotation_table` compress a trajectory or a rotation
  into Chebyshev coefficients, for inner loops and for code that runs on a GPU.
- **Orbits.** Orbital elements, two-body propagation and relative motion in local frames,
  after [ANISE](https://github.com/nyx-space/anise)
  ([how-to](https://orbitology.github.io/Almanac.jl/dev/tutorials/orbits/)).

The test suite compares Almanac.jl with CSPICE through SPICE.jl and replays NAIF's own test
suite against Almanac.jl. Results agree to within rounding error, apart from a few
[documented differences](https://orbitology.github.io/Almanac.jl/dev/manual/numerics/).
The package is under active development, and its interface may still change.

## Tutorials and Documentation

For information on using the package, see the
[documentation](https://orbitology.github.io/Almanac.jl/dev/). It opens with a getting-started
page that goes from loading kernels to an event search, followed by tutorials, a manual
on the relation to CSPICE, and the API reference. The docstring of each ported routine links to its NAIF page. If you are porting a
program from CSPICE or SPICE.jl, start with the
[conventions](https://orbitology.github.io/Almanac.jl/dev/manual/conventions/).

## Citing

If you publish work that uses Almanac.jl, please cite the SPICE toolkit as NAIF requests:

- C. H. Acton, "Ancillary data services of NASA's Navigation and Ancillary Information
  Facility", *Planetary and Space Science* 44 (1), 65–70, 1996.
- C. Acton, N. Bachman, B. Semenov, and E. Wright, "A look towards the future in the handling
  of space science mission geometry", *Planetary and Space Science* 150, 9–12, 2018.

## Contributing

Bug reports and questions are welcome in the
[issues](https://github.com/Orbitology/Almanac.jl/issues). Please report bugs in Almanac.jl
here, not to NAIF. To run the tests and build the manual locally:

```sh
julia --startup-file=no test/parallel.jl                    # tests
julia --project=docs -e 'using Pkg; Pkg.instantiate()'      # once
julia --project=docs docs/make.jl                           # manual, in docs/build/
```

## License

Almanac.jl is released under the [MIT license](LICENSE), which covers this repository only.
The repository contains no CSPICE source and no kernels, and the license grants no rights
over CSPICE, which NAIF distributes under its own
[rules of use](https://naif.jpl.nasa.gov/naif/rules.html).
