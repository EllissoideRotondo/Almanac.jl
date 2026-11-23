# Almanac.jl

Almanac.jl is a native Julia port of CSPICE N0067, the C toolkit of NASA's NAIF for
space geometry. It reads SPICE kernels and answers the questions CSPICE answers (where
a body is, how a frame is oriented, when an event occurs) in pure Julia.

The calls keep the names of CSPICE without the `_c` suffix, and its argument order with
[few exceptions](@ref "Names, arguments and results"); the docstring of each ported call
links to its NAIF page.

## Installation

Almanac.jl needs Julia 1.10 or later. Install it from the Julia package manager:

```julia-repl
pkg> add https://github.com/Orbitology/Almanac.jl
```

## Features

- Queries from any number of threads without locks: each query reads an immutable
  snapshot of the loaded kernels.
- Independent kernel sets, one per [`SpiceContext`](@ref).
- Event searches, surface geometry on ellipsoids and DSK shapes, and writers for SPK,
  CK, binary PCK, DSK and EK files.
- Beyond CSPICE: time derivatives of any order and automatic differentiation, tables
  for inner loops and GPU kernels, bodies and frames defined by Julia functions, and
  orbital elements after ANISE. [What Almanac.jl adds](@ref) lists them all.

Use the wrapper [SPICE.jl](https://github.com/JuliaAstro/SPICE.jl) instead when results
must be bit-identical to the CSPICE library on every platform.

## Contents

[Getting started](@ref) loads kernels, computes a state, a rotation and a sub-observer
point, and searches for events.

### Tutorials

Each tutorial solves one task, from several kernel sets and threads to orbits.

```@contents
Pages = [
    "tutorials/contexts.md",
    "tutorials/derivatives.md",
    "tutorials/tables.md",
    "tutorials/providers.md",
    "tutorials/events.md",
    "tutorials/writing.md",
    "tutorials/orbits.md",
]
Depth = 2
```

### Relation to CSPICE

These pages say where Almanac.jl follows CSPICE and where it departs from it; read
[Conventions](@ref) before porting a CSPICE or SPICE.jl program.

```@contents
Pages = [
    "manual/conventions.md",
    "manual/numerics.md",
    "manual/additions.md",
    "manual/cspice_defects.md",
]
Depth = 2
```

### Coverage and validation

What is ported, and how the results were checked against CSPICE and ANISE.

```@contents
Pages = [
    "manual/coverage.md",
    "manual/validation.md",
    "manual/orbit_validation.md",
]
Depth = 2
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

## The example kernels

The examples run on three small kernels that the manual writes itself: leap seconds,
constants of the Sun and the Earth with a spacecraft named `PROBE`, and two-body orbits
of the Earth and of `PROBE`, so the positions are not those of the real bodies.

## Reporting issues

Report defects and ask questions in the
[issues](https://github.com/Orbitology/Almanac.jl/issues). Almanac.jl is not developed
or supported by NAIF, so do not send its defects to NAIF.
