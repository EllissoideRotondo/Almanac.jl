# Custom bodies and frames

A trajectory that the program computes (with a propagator, an optimizer, an interpolation) can be used directly, without writing an SPK. [`add_ephemeris!`](@ref) makes a function the source of the state of a body, and [`add_frame!`](@ref) of the orientation of a frame. Every call that takes a body or a frame then uses them: with any aberration correction, in surface geometry, in event searches and in tables.

## Define a body

The function takes an epoch (ET) and returns the position (km) and the velocity (km/s) relative to the center, in the given frame. Here, a circular orbit of radius 8000 km around the Earth, for a body with the ID -1000, which no NAIF kernel uses:

```@example prov
using Almanac, SPICEKernels

furnsh(naif0012_tls())
furnsh(de440s())
furnsh(pck00011())

μ = 3.986004418e5                           # km³/s², the Earth
r = 8000.0                                  # km
n = sqrt(μ / r^3)                           # rad/s
circular(et) = [r * cos(n * et), r * sin(n * et), 0.0,
                -r * n * sin(n * et), r * n * cos(n * et), 0.0]

t0 = str2et("2025 MAR 01 00:00:00 UTC")
add_ephemeris!(circular, -1000, "EARTH", "J2000", t0, t0 + 86400.0; name = "my orbit")
state, light_time = spkezr(-1000, t0 + 600.0, "J2000", "LT+S", "MOON")
round.(state; digits = 3)
```

The body acts as if it were a file loaded at the time of the call: queries that start after the call use it, and a search that runs on another thread does not. It takes priority over the files loaded earlier, which apply outside its interval, and [`unload`](@ref) with its name removes it.

```@example prov
kinfo("my orbit")                           # (file type, meta-kernel or "", handle 0)
```

The function is called at every query. When it is expensive, or when the states exist as samples, load them as a kernel in memory instead, as [Kernels in memory](@ref) shows.

## Define a frame

The function returns the rotation from the new frame to its relative frame, the matrix of [`pxform`](@ref)`(new, relative, et)`. The other arguments are the name and the ID code of the new frame, its center and its relative frame:

```@example prov
spin(et) = (c = cos(n * et); s = sin(n * et); [c -s 0; s c 0; 0 0 1])
add_frame!(spin, "ORBIT_FIXED", 1_400_000, "EARTH", "J2000")

spkpos(-1000, t0 + 600.0, "ORBIT_FIXED", "NONE", "EARTH")[1]
```

In that frame the spacecraft stays on the X axis.

## Epoch types and thread safety

This section is about the functions of [`add_ephemeris!`](@ref) and [`add_frame!`](@ref); the function of [`load_function!`](@ref) gets `Float64` epochs, during the call only. The epoch passed to the function is not always a `Float64`. It is a `ForwardDiff.Dual` when the caller differentiates. It is a Taylor number, a type internal to Almanac.jl that carries derivatives and defines no comparisons, when a call needs a derivative that the function does not return: order 2 or more of a body (an `LT+S` state of an observer, [`spkder`](@ref)), and order 1 or more of a frame ([`sxform`](@ref), a corrected state in the frame). Plain arithmetic and the functions of `Base` (`sin`, `cos`, `sqrt`, ...) work with all of them. A function that branches on the epoch, or calls code that accepts only `Float64`, works only with a `Float64` epoch, in the calls that need no such derivative, such as geometric states and [`pxform`](@ref).

The function is called from every thread that queries, so it must not change shared data. Beware of closures that assign to a variable of the enclosing function, because they share it. Declare such variables `local` in the closure, or give them names that the enclosing function does not use.
