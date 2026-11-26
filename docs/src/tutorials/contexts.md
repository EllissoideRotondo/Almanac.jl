# Contexts and threads

Every call that reads kernels takes an optional [`SpiceContext`](@ref) as its first argument, written `[ctx,]` in the reference, and uses [`default_context`](@ref) without one. `[ctx | ks,]` means that the call also accepts an [`Almanac.KernelSet`](@ref), the immutable snapshot of a context. A call that fills a container takes the context after it. Pure functions such as `vnorm` take no context. A program with one set of kernels never needs to name a context; this page is for programs with several sets, or with several threads.

## Create a context

Contexts are independent of each other. Use one per version of an ephemeris, or one per mission.

```@example ctx
using Almanac, SPICEKernels

kernels = [naif0012_tls(), de440s(), pck00011()]
mission = SpiceContext()
foreach(k -> furnsh(mission, k), kernels)

(ktotal(mission, "ALL"), ktotal("ALL"))     # the default context is still empty
```

## Query from many threads

A query reads an immutable snapshot of the context and takes no lock, so any number of threads can query one context at once. Julia uses several threads only when it starts with `--threads` or `JULIA_NUM_THREADS`:

```@example ctx
et = str2et(mission, "2025 MAR 01 12:00:00 UTC")
epochs = et .+ range(0, 86400; length=10_000)
distance = Vector{Float64}(undef, length(epochs))

Threads.@threads for i in eachindex(epochs)
    position, _ = spkpos(mission, "MOON", epochs[i], "J2000", "LT+S", "EARTH")
    distance[i] = vnorm(position)
end

extrema(distance)                           # km
```

[`furnsh`](@ref), [`unload`](@ref) and the calls that change the pool ([`pdpool`](@ref), [`ldpool`](@ref), ...) may run while other threads query. Each of them builds a new snapshot and swaps it in atomically, so a query sees the kernels before or after the change, never a mixture. Concurrent changes are serialized. A change can land between two consecutive queries, so they can see different kernels.

## Prepare a query for many epochs

[`state_query`](@ref) and [`rotation_query`](@ref) resolve the names, the frame and the correction once, on the current snapshot. Evaluating the query, as `q(t)` for the state or `spkpos(q, t)` for the position, gives the result of [`spkezr`](@ref) or [`spkpos`](@ref) bit for bit; a rotation query gives that of [`pxform`](@ref). The query skips the name lookup at each epoch, and keeps its snapshot even if another thread loads a kernel.

```@example ctx
q = state_query(mission, "MOON", "J2000", "LT+S", "EARTH")

positions = [spkpos(q, t)[1] for t in epochs]
q(epochs[1]) == spkezr(mission, "MOON", epochs[1], "J2000", "LT+S", "EARTH")
```

To see kernels loaded later, make a new query.

## Share objects between threads

A file handle (from `dafopr`, `dasopr`, `spkopn`, ...) and a [`Window`](@ref) that is being changed must not be used from two threads at once. Queries ([`StateQuery`](@ref), [`RotationQuery`](@ref)) and tables ([`EphemerisTable`](@ref), [`RotationTable`](@ref)) are immutable, and safe to share.
