# Getting started

This page loads kernels, converts a time, computes a state and a rotation, finds the point below a spacecraft and searches for perigees.

```@setup first
using Almanac
include(joinpath(@__DIR__, "..", "fixture.jl"))
```

## Load kernels

This page runs on [the example kernels](@ref "The example kernels") of the manual. `DocKernels.write`, defined in [`docs/fixture.jl`](https://github.com/Orbitology/Almanac.jl/blob/main/docs/fixture.jl) of the repository and not part of the package, writes them into a directory and returns their three paths: leap seconds, constants of the Sun and the Earth, and an SPK for the year 2025. The SPK holds the Sun at the solar system barycenter, the Earth on a two-body orbit around the Sun, and the spacecraft `PROBE` (NAIF ID -999) on a two-body orbit around the Earth. To run the page without that file, load the real kernels shown below and replace `PROBE` with `MOON`.

```@example first
using Almanac

kernels = DocKernels.write(mktempdir())
basename.(kernels)
```

[`furnsh`](@ref) loads a kernel of any kind, including a meta-kernel, and [`ktotal`](@ref) counts the loaded kernels:

```@example first
foreach(furnsh, kernels)
ktotal("ALL")
```

As in CSPICE, a relative path in the `KERNELS_TO_LOAD` of a meta-kernel is relative to the working directory, not to the meta-kernel. Use `PATH_SYMBOLS` and `PATH_VALUES` for other locations.

For the real bodies, [SPICEKernels.jl](https://github.com/JuliaAstro/EphemerisSources.jl/tree/main/lib/SPICEKernels) downloads the NAIF generic kernels and returns their paths:

```julia
using SPICEKernels

furnsh(naif0012_tls())   # leap seconds
furnsh(de440s())         # planets, Moon and Sun, 1849 to 2150
furnsh(pck00011())       # orientation and radii of the bodies
```

## Convert a time

The state, frame, geometry and event calls take epochs in seconds past J2000 TDB, called ET; the CK calls take encoded spacecraft clock. [`str2et`](@ref) parses a calendar string, and [`timout`](@ref) formats an epoch with a picture string.

```@example first
et = str2et("2025 MAR 01 12:00:00 UTC")
```

```@example first
timout(et, "YYYY Mon DD HR:MN:SC.### ::TDB")
```

The offset of 69.185 s is TDB − UTC at that date.

## Compute a state

[`spkezr`](@ref) returns the state of a target relative to an observer (km, km/s) and the one-way light time (s). The correction `"NONE"` gives the geometric state, and `"LT+S"` the apparent one.

```@example first
state, lt = spkezr("PROBE", et, "J2000", "NONE", "EARTH")
state
```

The state is an `SVector{6}`, and the call does not [allocate](@ref "Allocation"). [`spkpos`](@ref) returns only the position:

```@example first
sun, lt_sun = spkpos("SUN", et, "J2000", "LT+S", "PROBE")
lt_sun                        # s
```

## Rotate between frames

[`pxform`](@ref) returns the 3×3 rotation from one frame to another, and [`sxform`](@ref) the 6×6 state transformation. [`mxv`](@ref) applies a matrix to a vector with the arithmetic of CSPICE.

```@example first
r = pxform("J2000", "IAU_EARTH", et)
mxv(r, state[1:3])            # km, in the Earth-fixed frame
```

## Find the point below the spacecraft

[`subpnt`](@ref) returns the sub-observer point in the body-fixed frame, the epoch at the target, and the vector from the observer to the point. [`reclat`](@ref) converts the point to radius, longitude and latitude (rad).

```@example first
point, epoch, vector = subpnt("NEAR POINT/ELLIPSOID", "EARTH", et, "IAU_EARTH", "NONE", "PROBE")
radius, lon, lat = reclat(point)
(lon * dpr(), lat * dpr())    # degrees
```

## Search for an event

The `gf` calls search a [`Window`](@ref) of time for the epochs where a quantity meets a condition, and return a `Window`. The search below finds the local minima of the distance over one day, with the relation `"LOCMIN"`. After the relation come the reference value and the adjustment, both unused for `"LOCMIN"`, and the step of 600 s; [Choose the step](@ref) explains how to pick it.

```@example first
search = Window()
wninsd!(search, et, et + 86400.0)
perigees = gfdist("PROBE", "NONE", "EARTH", "LOCMIN", 0.0, 0.0, 600.0, search)

for i in 1:wncard(perigees)
    t, _ = wnfetd(perigees, i)
    distance = vnorm(spkpos("PROBE", t, "J2000", "NONE", "EARTH")[1])
    println(timout(t, "YYYY Mon DD HR:MN:SC ::TDB"), "   ", round(distance; digits=3), " km")
end
```

The example orbit has a perigee radius of 7078.137 km.

## Unload kernels

[`kclear`](@ref) removes all kernels. [`unload`](@ref) removes one file, or a meta-kernel and the files it loaded.

```@example first
kclear()
ktotal("ALL")
```

## Handle errors

A call that cannot answer throws a [`SpiceError`](@ref). Its `kind` is the CSPICE short message as a lower-case symbol (`SPICE(IDCODENOTFOUND)` becomes `:idcodenotfound`), and `msg` the text. With no kernels loaded after `kclear`, the state of `PROBE` is unknown:

```@example first
try
    spkezr("PROBE", et, "J2000", "NONE", "EARTH")
catch err
    err.kind
end
```

[Conventions](@ref "Errors and `found` flags") gives the other rules for errors, among them the subtype for missing SPK data.

## Next steps

- Before porting a CSPICE program, read [Conventions](@ref).
- The tutorials start with [Contexts and threads](@ref).
