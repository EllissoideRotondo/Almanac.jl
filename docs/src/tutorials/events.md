# Event searches

The geometry finder, the `gf` calls, answers the question "when is this quantity equal to, above or below a value, or at an extremum?" Each search takes a confinement window, which holds the times to search, and returns the window on which the condition holds.

## Load the kernels

The examples of this page need the leap seconds, the planetary ephemeris and the shapes and rotations of the bodies. [SPICEKernels.jl](https://github.com/JuliaAstro/EphemerisSources.jl/tree/main/lib/SPICEKernels) downloads the NAIF generic kernels on first use and returns their paths:

```@example gf
using Almanac, SPICEKernels

furnsh(naif0012_tls())
furnsh(de440s())
furnsh(pck00011())
```

## Define the search window

A [`Window`](@ref) is a sorted list of disjoint intervals. [`wninsd!`](@ref) adds an interval, [`wncard`](@ref) counts the intervals, and [`wnfetd`](@ref) returns one of them, counting from 1. The searches below cover March 2025:

```@example gf
start = str2et("2025 MAR 01 00:00:00 TDB")
search = Window()
wninsd!(search, start, start + 31 * 86400.0)
wncard(search)
```

## Choose the step

A search samples its quantity at the given step, then refines each change to the
tolerance `tol` (s), whose default comes from [`gfstol`](@ref) and is 1 µs unless
changed. A step that is too long misses events without warning. The step must be
shorter than:

- for a condition that holds or not (`"="`, `"<"`, `">"`, an occultation, a field of
  view): the shortest interval on which the condition holds, and the shortest gap
  between two such intervals;
- for an extremum (`"LOCMIN"`, `"LOCMAX"`, `"ABSMIN"`, `"ABSMAX"`): the shortest
  interval on which the quantity is monotone.

## Find perigees

[`gfdist`](@ref) searches on the distance between two bodies. The relation `"LOCMIN"` finds its local minima, here the perigees of the Moon. After the relation come the reference value and the adjustment, both unused for `"LOCMIN"`, and the step. The distance of the Moon is monotone for about two weeks at a time, so a step of one day is enough:

```@example gf
perigees = gfdist("MOON", "NONE", "EARTH", "LOCMIN", 0.0, 0.0, 86400.0, search)

for i in 1:wncard(perigees)
    t, _ = wnfetd(perigees, i)
    distance = vnorm(spkpos("MOON", t, "J2000", "NONE", "EARTH")[1])
    println(timout(t, "YYYY Mon DD HR:MN ::TDB"), "   ", round(distance; digits=1), " km")
end
```

## Find eclipses

[`gfoclt`](@ref) finds when one body hides another from an observer. The Moon is in the shadow of the Earth when the Earth occults the Sun as seen from the Moon. The arguments are the occultation type, the name, shape and body-fixed frame of the front body, the same for the back body, the correction, the observer, the step and the window. The type `"ANY"` finds when the center of the Moon is in the penumbra, and `"FULL"` when it is in the umbra:

```@example gf
for kind in ("ANY", "FULL")
    eclipse = gfoclt(kind, "EARTH", "ELLIPSOID", "IAU_EARTH",
        "SUN", "ELLIPSOID", "IAU_SUN", "LT", "MOON", 600.0, search)
    for i in 1:wncard(eclipse)
        t0, t1 = wnfetd(eclipse, i)
        println(rpad(kind, 6), timout(t0, "YYYY Mon DD HR:MN ::TDB"), "   ", round((t1 - t0) / 3600; digits=2), " h")
    end
end
```

This is the total lunar eclipse of 14 March 2025.

## Search on a coordinate

[`gfposc`](@ref) searches on one coordinate of a position, here a declination of the Moon above 20°. After the relation come the reference value, the adjustment (used only for absolute extrema) and the step. [`wnsumd`](@ref)`[1]` is the total duration.

```@example gf
north = gfposc("MOON", "J2000", "NONE", "EARTH", "RA/DEC", "DECLINATION",
    ">", 20 * rpd(), 0.0, 21600.0, search)
(wncard(north), wnsumd(north)[1] / 86400)    # intervals, days
```

## Search on a function of your own

[`gfuds`](@ref) searches on any scalar function of time. It also needs a function that tells whether the first one decreases at an epoch. [`uddc`](@ref) derives it from a numerical derivative over the given interval, here 60 s.

```@example gf
speed(t) = vnorm(spkezr("MOON", t, "J2000", "NONE", "EARTH")[1][4:6])
fast = gfuds(speed, t -> uddc(speed, t, 60.0), ">", 1.05, 0.0, 21600.0, search)
wncard(fast)                                # intervals with a speed above 1.05 km/s
```

## Search on several threads

Every search takes the keyword `threads`, `false` by default. With `threads=true` a search processes intervals on separate threads, if Julia started with `--threads` or `JULIA_NUM_THREADS`, and its result does not change.

- Every search runs the intervals of its confinement window in parallel.
- `"="`, `"<"`, `">"`, and `"ABSMIN"` or `"ABSMAX"` with a nonzero adjustment, also run in parallel on the monotone intervals that they find first, so they gain even on a window of one interval.
- Other extrema, occultations and field-of-view searches gain nothing on a window of one interval.

The functions given to `gfuds` or `gfudb` are then called concurrently, so they must be thread-safe.
