# First session

This session loads kernels, converts a time, computes a state and a rotation, finds the point below a spacecraft and searches for perigees. Names are those of CSPICE without `_c`; a call that changes an argument ends in `!`.

```@setup first
using Almanac
include(joinpath(@__DIR__, "..", "fixture.jl"))
kernels = DocKernels.write(mktempdir())
```

## Load kernels

[`furnsh`](@ref) loads a kernel of any kind, including a meta-kernel; [`ktotal`](@ref) counts loaded kernels. Here `kernels` holds three paths: leap seconds, constants of the Sun and the Earth, and an SPK.

```@example first
using Almanac

foreach(furnsh, kernels)
ktotal("ALL")
```

As in CSPICE, relative paths in a meta-kernel's `KERNELS_TO_LOAD` are relative to the working directory, not to the meta-kernel; use `PATH_SYMBOLS` and `PATH_VALUES` for anything else.

For real bodies, [SPICEKernels.jl](https://github.com/JuliaAstro/EphemerisSources.jl/tree/main/lib/SPICEKernels) downloads the NAIF generic kernels and returns their paths:

```julia
using SPICEKernels

furnsh(naif0012_tls())   # leap seconds
furnsh(de440s())         # planets, Moon and Sun, 1849 to 2150
furnsh(pck00011())       # orientation and radii of the bodies
```

## Convert a time

State, frame, geometry and event calls take epochs in seconds past J2000 TDB (ET); the CK calls take encoded spacecraft clock. [`str2et`](@ref) parses a calendar string and [`timout`](@ref) formats an epoch with a picture string.

```@example first
et = str2et("2025 MAR 01 12:00:00 UTC")
```

```@example first
timout(et, "YYYY Mon DD HR:MN:SC.### ::TDB")
```

The 69.185 s offset is TDB − UTC at that date.

## State of a body

[`spkezr`](@ref) returns the state of a target relative to an observer (km, km/s) and the one-way light time (s). `"NONE"` gives the geometric state, `"LT+S"` the apparent one.

```@example first
state, lt = spkezr("PROBE", et, "J2000", "NONE", "EARTH")
state
```

The state is an `SVector{6}`. The call does not allocate when bodies and frame are integer IDs or upper-case names with single blanks (`"MARS BARYCENTER"`). [`spkpos`](@ref) returns the position only:

```@example first
sun, lt_sun = spkpos("SUN", et, "J2000", "LT+S", "PROBE")
lt_sun                        # s
```

## Rotation between frames

[`pxform`](@ref) returns the 3×3 rotation from one frame to another, [`sxform`](@ref) the 6×6 state transformation.

```@example first
r = pxform("J2000", "IAU_EARTH", et)
r * state[1:3]                # km, in the Earth-fixed frame
```

## The point below the spacecraft

[`subpnt`](@ref) returns the sub-observer point in the body-fixed frame, the epoch at the target and the observer-to-point vector. [`reclat`](@ref) converts it to radius, longitude and latitude (rad).

```@example first
point, epoch, vector = subpnt("NEAR POINT/ELLIPSOID", "EARTH", et, "IAU_EARTH", "NONE", "PROBE")
radius, lon, lat = reclat(point)
(lon * dpr(), lat * dpr())    # degrees
```

## Search for an event

The `gf` calls search a [`Window`](@ref) of time for the epochs where a quantity meets a condition, and return a `Window`. Below, the local minima of the distance over one day. The arguments after the relation are the reference value and the adjustment (both unused for `"LOCMIN"`) and the step, 600 s; [The step](@ref) says how to choose it.

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

## Unload

[`kclear`](@ref) removes all kernels; [`unload`](@ref) removes one file, or a meta-kernel and the files it loaded.

```@example first
kclear()
ktotal("ALL")
```

## Errors

A call that cannot answer throws a [`SpiceError`](@ref), usually its subtype [`SpiceKindError`](@ref). Its `kind` is the CSPICE short message as a lower-case symbol (`SPICE(IDCODENOTFOUND)` becomes `:idcodenotfound`), and `msg` the text.

```@example first
try
    spkezr("PROBE", et, "J2000", "NONE", "EARTH")
catch err
    err.kind
end
```

Missing SPK data throw the subtype [`SpkInsufficientDataError`](@ref), whose fields name the bodies, the epoch, and whether the missing body has data at other epochs.

Before porting a CSPICE program, read [Differences from CSPICE](@ref).
