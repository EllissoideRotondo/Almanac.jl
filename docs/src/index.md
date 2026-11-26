# Almanac.jl

Almanac.jl is a native Julia port of CSPICE N0067, the C toolkit of NASA's NAIF for
space geometry. It reads SPICE kernels and answers the questions CSPICE answers (where
a body is, how a frame is oriented, when an event occurs) in pure Julia.

## Installation

```julia-repl
pkg> add Almanac
```

## Quickstart

```julia
using Almanac
using SPICEKernels

furnsh(naif0012_tls())
furnsh(de440s())
furnsh(pck00011())

# Convert a UTC time to ephemeris time (TDB seconds past J2000), and back
et = str2et("2025-03-01T12:00:00")
et2utc(et, "ISOC", 3)

# State of the Moon relative to the Earth (km, km/s) and the one-way light time (s),
# corrected for light time and stellar aberration
state, lt = spkezr("MOON", et, "J2000", "LT+S", "EARTH")

# Rotate the position into the Earth-fixed frame
r = pxform("J2000", "IAU_EARTH", et)
mxv(r, state[1:3])
```

## Contents

### Getting started

```@contents
Pages = [
    "tutorials/contexts.md",
    "tutorials/derivatives.md",
    "tutorials/tables.md",
    "tutorials/providers.md",
    "tutorials/writing.md",
    "tutorials/orbits.md",
    "tutorials/events.md",
]
Depth = 1
```

### Validation


```@contents
Pages = ["manual/numerics.md", "manual/validation.md"]
Depth = 1
```

### API reference


```@contents
Pages = [
    "reference/kernels.md",
    "reference/time.md",
    "reference/states.md",
    "reference/frames.md",
    "reference/geometry.md",
    "reference/events.md",
    "reference/orbits.md",
    "reference/math/vectors.md",
    "reference/math/rotations.md",
    "reference/math/coordinates.md",
    "reference/math/ellipsoids.md",
    "reference/math/conics.md",
    "reference/files/daf.md",
    "reference/files/writers.md",
    "reference/files/das.md",
    "reference/files/ek.md",
]
Depth = 1
```

## Reporting issues

Report defects and ask questions in the
[issues](https://github.com/Orbitology/Almanac.jl/issues). Almanac.jl is not developed
or supported by NAIF, so do not send its defects to NAIF.
