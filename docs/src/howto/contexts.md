# Contexts and threads

Loaded kernels, the kernel pool and settings such as [`timdef`](@ref) belong to a [`SpiceContext`](@ref). Every call that reads them takes an optional context as first argument; without it, the call uses [`default_context`](@ref). Pure functions such as `vnorm` take no context. A program with one set of kernels never names a context.

```@setup ctx
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
kernels = DocKernels.write(mktempdir())
```

## Several kernel sets

Each context is independent: use one per ephemeris version or per mission.

```@example ctx
using Almanac

mission = SpiceContext()
foreach(k -> furnsh(mission, k), kernels)

(ktotal(mission, "ALL"), ktotal("ALL"))     # the default context is still empty
```

## Queries from many threads

A query reads an immutable snapshot of the context and takes no lock, so any number of threads can query one context at once:

```@example ctx
et = str2et(mission, "2025 MAR 01 12:00:00 UTC")
epochs = et .+ range(0, 86400; length=10_000)
distance = Vector{Float64}(undef, length(epochs))

Threads.@threads for i in eachindex(epochs)
    position, _ = spkpos(mission, "PROBE", epochs[i], "J2000", "LT+S", "EARTH")
    distance[i] = vnorm(position)
end

extrema(distance)                           # km
```

[`furnsh`](@ref), [`unload`](@ref) and the calls that change the pool ([`pdpool`](@ref), [`ldpool`](@ref), ...) may run while other threads query. Each builds a new snapshot and swaps it in atomically, so a query sees the set before or after the change, never a mixture; concurrent changes are serialized. Two consecutive queries can therefore see different sets.

## Resolve a query once

[`state_query`](@ref) and [`rotation_query`](@ref) resolve names, frame and correction once, on the current snapshot. Evaluating the query, `q(t)` for the state or `spkpos(q, t)` for the position, gives the result of [`spkezr`](@ref) or [`spkpos`](@ref) bit for bit (a rotation query that of [`pxform`](@ref)), skips the name lookup at each epoch, and keeps its snapshot even if another thread loads a kernel.

```@example ctx
q = state_query(mission, "PROBE", "J2000", "LT+S", "EARTH")

positions = [spkpos(q, t)[1] for t in epochs]
q(epochs[1]) == spkezr(mission, "PROBE", epochs[1], "J2000", "LT+S", "EARTH")
```

Make a new query to see kernels loaded later.

## What must not be shared

A file handle (from `dafopr`, `dasopr`, `spkopn`, ...) and a [`Window`](@ref) that is being changed must not be used from two threads at once. Queries ([`StateQuery`](@ref), [`RotationQuery`](@ref)) and tables ([`EphemerisTable`](@ref), [`RotationTable`](@ref)) are immutable and safe to share.
