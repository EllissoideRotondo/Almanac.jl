# Almanac.jl

Almanac.jl is a native Julia port of CSPICE N0067. It reads SPICE kernels and answers the questions CSPICE answers (where a body is, how a frame is oriented, when an event occurs) in pure Julia. Calls keep the CSPICE names, without the `_c` suffix, and with [few exceptions](@ref "Call conventions") the argument order; each docstring links to the NAIF page of the routine it ports.

## Installation

```julia-repl
pkg> add https://github.com/Orbitology/Almanac.jl
```

Julia 1.10 or later.

## Features

- Lock-free queries from any number of threads; each query reads an immutable snapshot of the loaded kernels.
- Independent kernel sets, one per [`SpiceContext`](@ref).
- Time derivatives of any order ([`spkder`](@ref), [`pxder`](@ref)) and automatic differentiation with ForwardDiff.jl ([dual-number epochs](@ref "Dual-number epochs")), Enzyme.jl, Mooncake.jl, Zygote.jl, ReverseDiff.jl and Tracker.jl ([other AD packages](@ref "Other AD packages")).
- [Ephemeris and rotation tables](@ref) for inner loops and GPU kernels.
- Event searches, surface geometry on ellipsoids and DSK shapes, and writers for SPK, CK, binary PCK, DSK and EK files.
- Bodies and frames defined by Julia functions ([`add_ephemeris!`](@ref), [`add_frame!`](@ref)).
- Orbital elements, two-body propagation and relative motion in local frames ([`Orbit`](@ref)), after ANISE.

Use the wrapper [SPICE.jl](https://github.com/JuliaAstro/SPICE.jl) instead when results must be bit-identical to the CSPICE library on every platform.

## Where to start

- [First session](@ref): load kernels, compute a state, a rotation and a sub-observer point, search for events.
- [Differences from CSPICE](@ref): read before porting a CSPICE or SPICE.jl program.
- The how-to guides cover [threads and contexts](@ref "Contexts and threads"), [derivatives](@ref "Time derivatives"), [tables](@ref "Ephemeris and rotation tables"), [user-defined bodies](@ref "Bodies and frames from your own functions"), [event searches](@ref "Event searches") and [writing kernels](@ref "Writing kernels").
- [Orbits and relative motion](@ref): the calls beyond CSPICE.

The examples run on three small kernels that the manual writes itself: leap seconds, constants of the Sun and the Earth with a spacecraft named `PROBE`, and two-body orbits of the Earth and of `PROBE`. Positions are therefore not those of the real bodies.

## Contributing

Report defects and ask questions in the [issues](https://github.com/Orbitology/Almanac.jl/issues). Almanac.jl is not developed or supported by NAIF; do not send its defects to NAIF.
