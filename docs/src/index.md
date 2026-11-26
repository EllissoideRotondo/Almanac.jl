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

## Reporting issues

Report defects and ask questions in the
[issues](https://github.com/Orbitology/Almanac.jl/issues). Almanac.jl is not developed
or supported by NAIF, so do not send its defects to NAIF.
