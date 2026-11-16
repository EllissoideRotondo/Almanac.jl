# Event searches

The geometry finder (`gf` calls) answers "when is this quantity equal to, above or below a value, or at an extremum?". Each search takes a confinement window of times to search and returns the window where the condition holds.

```@setup gf
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Windows

A [`Window`](@ref) is a sorted list of disjoint intervals. [`wninsd!`](@ref) adds an interval, [`wncard`](@ref) counts them, [`wnfetd`](@ref) returns one (indices from 1).

```@example gf
using Almanac

start = str2et("2025 MAR 18 00:00:00 TDB")
search = Window()
wninsd!(search, start, start + 2 * 86400.0)
wncard(search)
```

## The step

A search samples its quantity at the given step, then refines each change to the tolerance `tol` (s); the default comes from [`gfstol`](@ref), 1 µs unless changed. A step that is too long misses events silently:

| Search | The step must be shorter than |
|---|---|
| A condition that holds or not (`"="`, `"<"`, `">"`, an occultation, a field of view) | the shortest interval on which the condition holds, and the shortest gap between two such intervals |
| An extremum (`"LOCMIN"`, `"LOCMAX"`, `"ABSMIN"`, `"ABSMAX"`) | the shortest interval on which the quantity is monotone |

## Eclipses

[`gfoclt`](@ref) finds when one body hides another from an observer. `PROBE` is in the Earth's shadow when the Earth occults the Sun as seen from `PROBE`. Arguments: the occultation type (`"ANY"`: partial or full), name, shape and body-fixed frame of the front body, then of the back body, correction, observer, step, window. The first eclipse lasts under two minutes, hence the 60 s step; a 300 s step can step over it.

```@example gf
eclipses = gfoclt("ANY", "EARTH", "ELLIPSOID", "IAU_EARTH",
    "SUN", "ELLIPSOID", "IAU_SUN", "LT", "PROBE", 60.0, search)

for i in 1:wncard(eclipses)
    t0, t1 = wnfetd(eclipses, i)
    println(timout(t0, "YYYY Mon DD HR:MN:SC ::TDB"), "   ", round(t1 - t0; digits=1), " s")
end
```

## A condition on a coordinate

[`gfposc`](@ref) searches on one coordinate of a position: here, latitude above 20°. After the relation: reference value, adjustment (absolute extrema only), step. [`wnsumd`](@ref)`[1]` is the total duration.

```@example gf
north = gfposc("PROBE", "IAU_EARTH", "NONE", "EARTH", "LATITUDINAL", "LATITUDE",
    ">", 20 * rpd(), 0.0, 600.0, search)
wnsumd(north)[1] / 3600                     # hours
```

## A function of your own

[`gfuds`](@ref) searches on any scalar function of time. It also needs a function that tells whether the first decreases at an epoch; [`uddc`](@ref) derives it from a numerical derivative over the given interval, here 10 s.

```@example gf
speed(t) = vnorm(spkezr("PROBE", t, "J2000", "NONE", "EARTH")[1][4:6])
fast = gfuds(speed, t -> uddc(speed, t, 10.0), ">", 8.0, 0.0, 600.0, search)
wncard(fast)                                # intervals with a speed above 8 km/s
```

## Several threads

Every search takes `threads=true`, which processes intervals on separate threads. The intervals are those of the confinement window and, for `"="`, `"<"`, `">"`, and `"ABSMIN"` or `"ABSMAX"` with a nonzero adjustment, the monotone intervals found first. Other extrema, occultations and field-of-view searches on a one-interval window gain nothing. Functions given to `gfuds` or `gfudb` are then called concurrently and must be thread-safe.
