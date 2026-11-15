# Bodies and frames from your own functions

A trajectory that a program computes (the result of a propagator, an optimizer or an
interpolation) need not be written to an SPK file before the geometry calls can use it.
[`add_ephemeris!`](@ref) makes a function the source of the state of a body, and
[`add_frame!`](@ref) makes a function the source of the orientation of a frame. Every
call that takes a body or a frame then uses them: states with any aberration
correction, surface geometry, event searches, tables.

```@setup prov
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## A body

The function receives an epoch (ET, seconds past J2000 TDB) and returns the position
(km) and the velocity (km/s) relative to the center, in the frame of the call. Here a
second spacecraft moves in a circular orbit of 8000 km around the Earth.

```@example prov
using Almanac

μ = 3.986004418e5                           # km^3/s^2, the Earth
r = 8000.0                                  # km
n = sqrt(μ / r^3)                           # rad/s
circular(et) = [r * cos(n * et), r * sin(n * et), 0.0,
                -r * n * sin(n * et), r * n * cos(n * et), 0.0]

t0 = str2et("2025 MAR 01 00:00:00 UTC")
add_ephemeris!(circular, -1000, "EARTH", "J2000", t0, t0 + 86400.0; name = "my orbit")
state, light_time = spkezr(-1000, t0 + 600.0, "J2000", "LT+S", "PROBE")
round.(state; digits = 3)
```

The body behaves as a file loaded at that moment: it has priority over the files loaded
before it, the files apply outside its interval, and [`unload`](@ref) with its name
removes it.

```@example prov
kinfo("my orbit")                           # (type, meta-kernel, handle)
```

## A frame

The function returns the matrix that takes vectors from the new frame to the frame
named as relative, the direction of [`pxform`](@ref)`(new, relative, et)`.

```@example prov
spin(et) = (c = cos(n * et); s = sin(n * et); [c -s 0; s c 0; 0 0 1])
add_frame!(spin, "ORBIT_FIXED", 1_400_000, "EARTH", "J2000")

spkpos(-1000, t0 + 600.0, "ORBIT_FIXED", "NONE", "EARTH")[1]
```

In that frame the spacecraft stands still on the X axis.

## Derivatives and epochs that are not numbers

A call can pass the function an epoch that is not a `Float64`: a `ForwardDiff.Dual`
when the caller differentiates, and a Taylor number of this package for derivatives of
order 2 and more, which an `LT+S` state of an observer and [`sxform`](@ref) of a frame
need. A function that uses arithmetic and the functions of `Base` (`sin`, `cos`, `sqrt`,
...) works with all of them. A function that cannot (one that branches on the epoch, or
calls a library that takes only `Float64`) works only in the calls that need no
derivative of order 2, such as geometric states and positions and [`pxform`](@ref).

The function is called from any thread that makes a query, so it must not change data
that other threads read. A closure changes such data without showing it when it assigns
a variable that the enclosing function also assigns: the two share the variable. Declare
the variables of the closure `local`, or give them names of their own.
