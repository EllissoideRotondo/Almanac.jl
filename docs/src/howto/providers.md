# Bodies and frames from your own functions

A trajectory computed by the program (a propagator, an optimizer, an interpolation) can be used directly, without writing an SPK. [`add_ephemeris!`](@ref) makes a function the source of a body's state, [`add_frame!`](@ref) of a frame's orientation. All calls that take a body or a frame then use them: any aberration correction, surface geometry, event searches, tables.

```@setup prov
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## A body

The function takes ET and returns position (km) and velocity (km/s) relative to the center, in the given frame. Here, a circular orbit of 8000 km around the Earth:

```@example prov
using Almanac

μ = 3.986004418e5                           # km³/s², the Earth
r = 8000.0                                  # km
n = sqrt(μ / r^3)                           # rad/s
circular(et) = [r * cos(n * et), r * sin(n * et), 0.0,
                -r * n * sin(n * et), r * n * cos(n * et), 0.0]

t0 = str2et("2025 MAR 01 00:00:00 UTC")
add_ephemeris!(circular, -1000, "EARTH", "J2000", t0, t0 + 86400.0; name = "my orbit")
state, light_time = spkezr(-1000, t0 + 600.0, "J2000", "LT+S", "PROBE")
round.(state; digits = 3)
```

The body acts as a file loaded at that moment: it takes priority over earlier files, those files still apply outside its interval, and [`unload`](@ref) with its name removes it.

```@example prov
kinfo("my orbit")                           # (type, meta-kernel, handle)
```

## A frame

The function returns the rotation from the new frame to its relative frame, the matrix of [`pxform`](@ref)`(new, relative, et)`.

```@example prov
spin(et) = (c = cos(n * et); s = sin(n * et); [c -s 0; s c 0; 0 0 1])
add_frame!(spin, "ORBIT_FIXED", 1_400_000, "EARTH", "J2000")

spkpos(-1000, t0 + 600.0, "ORBIT_FIXED", "NONE", "EARTH")[1]
```

In that frame the spacecraft stays on the X axis.

## Requirements on the function

The epoch passed in is not always a `Float64`: it is a `ForwardDiff.Dual` when the caller differentiates, and an internal Taylor number when a call needs derivatives of order 2 or more (an `LT+S` state of an observer, [`sxform`](@ref) of a frame). Plain arithmetic and `Base` functions (`sin`, `cos`, `sqrt`, ...) work with all of them. A function that branches on the epoch or calls `Float64`-only code works only in calls that need no derivative of order 2, such as geometric states and [`pxform`](@ref), and only with a `Float64` epoch.

The function is called from every thread that queries, so it must not mutate shared data. Beware of closures that assign a variable of the enclosing function: they share it. Declare such variables `local` in the closure, or give them names the enclosing function does not use.
