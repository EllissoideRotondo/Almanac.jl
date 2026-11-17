# Event searches

The geometry finder, the `gf` calls, answers the question "when is this quantity equal to, above or below a value, or at an extremum?" Each search takes a confinement window, which holds the times to search, and returns the window on which the condition holds.

```@setup gf
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Windows

A [`Window`](@ref) is a sorted list of disjoint intervals. [`wninsd!`](@ref) adds an interval, [`wncard`](@ref) counts the intervals, and [`wnfetd`](@ref) returns one of them, counting from 1.

```@example gf
using Almanac

start = str2et("2025 MAR 18 00:00:00 TDB")
search = Window()
wninsd!(search, start, start + 2 * 86400.0)
wncard(search)
```

## The step

A search samples its quantity at the given step, then refines each change to the
tolerance `tol` (s), whose default comes from [`gfstol`](@ref) and is 1 µs unless
changed. A step that is too long misses events without warning. The step must be
shorter than:

- for a condition that holds or not (`"="`, `"<"`, `">"`, an occultation, a field of
  view): the shortest interval on which the condition holds, and the shortest gap
  between two such intervals;
- for an extremum (`"LOCMIN"`, `"LOCMAX"`, `"ABSMIN"`, `"ABSMAX"`): the shortest
  interval on which the quantity is monotone.

## Eclipses

[`gfoclt`](@ref) finds when one body hides another from an observer. `PROBE` is in the shadow of the Earth when the Earth occults the Sun as seen from `PROBE`. The arguments are the occultation type (`"ANY"` for partial or full), the name, shape and body-fixed frame of the front body, the same for the back body, the correction, the observer, the step and the window. The first eclipse lasts less than two minutes, hence the step of 60 s; a step of 300 s can miss it.

```@example gf
eclipses = gfoclt("ANY", "EARTH", "ELLIPSOID", "IAU_EARTH",
    "SUN", "ELLIPSOID", "IAU_SUN", "LT", "PROBE", 60.0, search)

for i in 1:wncard(eclipses)
    t0, t1 = wnfetd(eclipses, i)
    println(timout(t0, "YYYY Mon DD HR:MN:SC ::TDB"), "   ", round(t1 - t0; digits=1), " s")
end
```

## A condition on a coordinate

[`gfposc`](@ref) searches on one coordinate of a position, here a latitude above 20°. After the relation come the reference value, the adjustment (used only for absolute extrema) and the step. [`wnsumd`](@ref)`[1]` is the total duration.

```@example gf
north = gfposc("PROBE", "IAU_EARTH", "NONE", "EARTH", "LATITUDINAL", "LATITUDE",
    ">", 20 * rpd(), 0.0, 600.0, search)
wnsumd(north)[1] / 3600                     # hours
```

## Search on a function of your own

[`gfuds`](@ref) searches on any scalar function of time. It also needs a function that tells whether the first one decreases at an epoch. [`uddc`](@ref) derives it from a numerical derivative over the given interval, here 10 s.

```@example gf
speed(t) = vnorm(spkezr("PROBE", t, "J2000", "NONE", "EARTH")[1][4:6])
fast = gfuds(speed, t -> uddc(speed, t, 10.0), ">", 8.0, 0.0, 600.0, search)
wncard(fast)                                # intervals with a speed above 8 km/s
```

## Search on several threads

Every search takes the keyword `threads=true`, which processes intervals on separate threads.

- Every search runs the intervals of its confinement window in parallel.
- `"="`, `"<"`, `">"`, and `"ABSMIN"` or `"ABSMAX"` with a nonzero adjustment, also run in parallel on the monotone intervals that they find first, so they gain even on a window of one interval.
- Other extrema, occultations and field-of-view searches gain nothing on a window of one interval.

The functions given to `gfuds` or `gfudb` are then called concurrently, so they must be thread-safe.
