# Almanac.jl

Almanac.jl is a pure-Julia port of CSPICE (toolkit version N0067), the library of NASA's
NAIF for spacecraft and planetary geometry, which reads SPICE kernels and computes
positions, velocities, orientations, time conversions, surface geometry and event times.

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

## Notable differences to CSPICE

- **No global state.** Loaded kernels, the kernel pool and settings such as
  [`timdef`](@ref) belong to a [`SpiceContext`](@ref), which a call takes as an
  optional first argument; without one, it uses [`default_context`](@ref). Queries
  take no lock, so any number of threads can query while another loads kernels.
- **Cells and windows.** A CSPICE cell is a `Vector` or a [`Window`](@ref), with no
  declared size. The calls that fill one, such as [`spkcov!`](@ref) and
  [`ckcov!`](@ref), take it as the first argument, where CSPICE takes it last.
- **Indices count from 1**, where CSPICE counts from 0: the interval of
  [`wnfetd`](@ref), the `which` of [`kdata`](@ref), the `start` of [`gnpool`](@ref)
  and [`gcpool`](@ref), and the rows and segments of an EK. An index copied unchanged
  from C code selects the element before the one intended.
- **No memory arguments.** String lengths, `room`, workspaces and counts such as
  `nintvls`, which serve only C memory, are absent: [`gdpool`](@ref) takes only the
  name of the variable.
- **Routines left out.** Error handling (`chkin`, `failed`, `reset`) is replaced by
  exceptions: a call throws a [`SpiceError`](@ref). The cell, set and string utilities
  are replaced by Base functions. Routines deprecated by NAIF are replaced by their
  successors ([`bodvrd`](@ref) for `bodvar`, [`subpnt`](@ref) for `subpt`). Loading by
  handle (`spklef`, `cklpf`, `pcklof`) is replaced by [`furnsh`](@ref).
