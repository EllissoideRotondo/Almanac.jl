# Almanac.jl

Almanac.jl is a native Julia port of CSPICE N0067, the C toolkit of NASA's NAIF for
space geometry. It reads SPICE kernels and answers the questions CSPICE answers (where
a body is, how a frame is oriented, when an event occurs) in pure Julia.

The calls keep the names of CSPICE without the `_c` suffix, and its argument order with
[few exceptions](@ref "Names, arguments and results"); the docstring of each ported call
links to its NAIF page. Use the wrapper [SPICE.jl](https://github.com/JuliaAstro/SPICE.jl)
instead when results must be bit-identical to the CSPICE library on every platform.

## Installation

Almanac.jl needs Julia 1.10 or later. Install it from the Julia package manager:

```julia-repl
pkg> add Almanac
```

## Quickstart

[SPICEKernels.jl](https://github.com/JuliaAstro/EphemerisSources.jl/tree/main/lib/SPICEKernels)
downloads the NAIF generic kernels and returns their paths.

```julia
using Almanac
using SPICEKernels

furnsh(naif0012_tls())   # leap seconds
furnsh(de440s())         # planets, Moon and Sun, 1849 to 2150
furnsh(pck00011())       # orientation and radii of the bodies

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

Results are returned, not written into arguments. A call that cannot answer throws a
[`SpiceError`](@ref). Before porting a CSPICE or SPICE.jl program, read [Conventions](@ref).

## Contents

### Getting started

Each tutorial solves one task on the NAIF generic kernels, which SPICEKernels.jl
downloads.

```@contents
Pages = [
    "tutorials/events.md",
    "tutorials/contexts.md",
    "tutorials/derivatives.md",
    "tutorials/tables.md",
    "tutorials/providers.md",
    "tutorials/writing.md",
    "tutorials/orbits.md",
]
Depth = 1
```

### Porting from CSPICE

How Almanac.jl relates to CSPICE: its conventions, what it ports, what it adds, and
the defects of CSPICE that it copies or fixes.

```@contents
Pages = [
    "manual/conventions.md",
    "manual/coverage.md",
    "manual/additions.md",
    "manual/cspice_defects.md",
]
Depth = 1
```

### Accuracy and validation

How close the results are to those of CSPICE, and how they were tested.

```@contents
Pages = ["manual/numerics.md", "manual/validation.md", "manual/orbit_validation.md"]
Depth = 1
```

### API reference

The docstring of every exported call, grouped by subject.

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
