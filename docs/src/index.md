# Almanac.jl

Almanac.jl is a native Julia port of CSPICE N0067, the C toolkit of NASA's NAIF for
space geometry. It reads SPICE kernels and answers the questions CSPICE answers (where
a body is, how a frame is oriented, when an event occurs) in pure Julia.

The calls keep the names of CSPICE without the `_c` suffix, and its argument order with
[few exceptions](@ref "Names, arguments and results"). Each docstring links to the NAIF
page of the routine it ports.

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

## How this manual is organized

- [Getting started](@ref) loads kernels, computes a state, a rotation and a
  sub-observer point, and searches for events.
- The **Tutorials** each solve one task: [contexts and threads](@ref "Contexts and
  threads"), [time derivatives](@ref "Time derivatives"),
  [tables](@ref "Ephemeris and rotation tables"), [bodies and frames from your own
  functions](@ref "Bodies and frames from your own functions"), [event
  searches](@ref "Event searches"), [writing kernels](@ref "Writing kernels") and
  [orbits](@ref "Orbits and relative motion").
- The **Manual** explains how Almanac.jl relates to CSPICE, and how it was checked.
  Read [Conventions](@ref) before porting a CSPICE or SPICE.jl program.
  [Coverage](@ref), [Validation](@ref) and [Validation of the orbit calls](@ref) say
  what is available and how it was tested.
- The **API reference** documents every exported call, grouped by subject.

## The example kernels

The examples run on three small kernels that the manual writes itself: leap seconds,
constants of the Sun and the Earth with a spacecraft named `PROBE`, and two-body orbits
of the Earth and of `PROBE`. The positions are therefore not those of the real bodies.

## Contributing

Report defects and ask questions in the
[issues](https://github.com/Orbitology/Almanac.jl/issues). Almanac.jl is not developed
or supported by NAIF, so do not send its defects to NAIF.
