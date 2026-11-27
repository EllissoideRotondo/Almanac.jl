# Benchmark against Ephemerides.jl: Almanac, Ephemerides.jl and CSPICE (SPICE.jl).
#
#     julia --startup-file=no --project=benchmark benchmark/ephemerides.jl
#
# Ephemerides.jl evaluates one link of a kernel: the state of a target relative to the
# center of its segment, in the frame of the segment, with no chain and no change of
# frame. Every case is such a link, so the three engines make the same computation:
# Almanac and CSPICE through `spkgps`, `spkgeo` and `spkder` in the frame of the segment,
# Ephemerides.jl through `ephem_vector3` to `ephem_vector12`. Almanac is also timed with
# a prepared query (`state_query`, `rotation_query`): `q(et)` returns the light time
# too, and `spkder(q, et, Val(1))` returns the position and the velocity only, as
# Ephemerides.jl does.
#
# Cases:
#
# - de440s.bsp (type 2): the Earth relative to the Earth-Moon barycenter at orders 0 to
#   3, and the ForwardDiff derivative of the position, at the 100 epochs of the ANISE
#   benchmark.
# - The other SPK types: the first J2000 segment of each type in the kernel of the test
#   group `spktypes` (`test/spktypefixtures.jl`), state only, at 100 epochs evenly
#   spaced over the segment.
# - Binary PCK type 2: the rotation from J2000 to MOON_PA_DE440 of
#   moon_pa_de440_200625.bpc. Almanac and CSPICE call `pxform`; Ephemerides.jl gives
#   the Euler angles, which `eul2m` of Almanac turns into the matrix.
# - de440s, the Moon relative to the Earth: a chain of two links, which Ephemerides.jl
#   does not follow; Almanac and CSPICE only.
# - The load of de440s.bsp.
#
# Before any timing, each segment of the type kernel is compared with CSPICE in the
# frame of the segment, and the result file lists the largest relative differences. A
# case where Ephemerides.jl throws or differs from CSPICE by more than `EPH_RTOL` is
# timed without Ephemerides.jl.
#
# The script runs on one thread: Ephemerides.jl 1.2.1 indexes its caches by
# `Threads.threadid()`.
#
# Ephemerides.jl is not in `benchmark/Project.toml`, which the AirspeedVelocity job of
# pull requests installs. Add it to the environment before a run:
#
#     julia --project=benchmark -e 'using Pkg; Pkg.add("Ephemerides")'

include("common.jl")
if Base.find_package("Ephemerides") === nothing
    printstyled("SKIP  Ephemerides.jl is not in the environment\n"; color = :yellow, bold = true)
    exit(0)
end
using Almanac
using Ephemerides: Ephemerides, EphemerisProvider, ephem_vector3, ephem_vector6,
    ephem_vector9, ephem_vector12, ephem_rotation3
import ForwardDiff
import SPICE
isdefined(@__MODULE__, :write_type_kernel) ||
    include(joinpath(@__DIR__, "..", "test", "spktypefixtures.jl"))

require_kernels(Kernels.DE440S, Kernels.MOON_BPC, Kernels.MOON_TF)

const EMB, MOON, J2000, MOON_PA = 3, 301, 1, 31008

"Largest relative difference from CSPICE for which a case is timed with Ephemerides.jl."
const EPH_RTOL = 1.0e-12

# One loop per call form. Every loop sums one component of the highest order, so no
# part of the call is removed.
function loop_alm_gps(ctx, epochs, targ, ref, obs)
    acc = 0.0
    for et in epochs
        acc += spkgps(ctx, targ, et, ref, obs)[1][1]
    end
    return acc
end
function loop_alm_geo(ctx, epochs, targ, ref, obs)
    acc = 0.0
    for et in epochs
        acc += spkgeo(ctx, targ, et, ref, obs)[1][4]
    end
    return acc
end
function loop_alm_der(ctx, epochs, targ, ref, obs, v::Val{N}) where {N}
    acc = 0.0
    for et in epochs
        acc += spkder(ctx, targ, et, ref, obs, v)[N + 1][1]
    end
    return acc
end
function loop_alm_ad(ctx, epochs, targ, ref, obs)
    acc = 0.0
    for et in epochs
        acc += ForwardDiff.derivative(t -> spkgps(ctx, targ, t, ref, obs)[1], et)[1]
    end
    return acc
end
function loop_alm_query(q, epochs)
    acc = 0.0
    for et in epochs
        acc += q(et)[1][4]
    end
    return acc
end
function loop_alm_query_der(q, epochs)
    acc = 0.0
    for et in epochs
        acc += spkder(q, et, Val(1))[2][1]
    end
    return acc
end
function loop_alm_query_pos(q, epochs)
    acc = 0.0
    for et in epochs
        acc += spkpos(q, et)[1][1]
    end
    return acc
end
function loop_alm_query_rot(q, epochs)
    acc = 0.0
    for et in epochs
        acc += q(et)[3, 1]
    end
    return acc
end
function loop_eph(f::F, eph, epochs, center, targ, k) where {F}
    acc = 0.0
    for et in epochs
        acc += f(eph, center, targ, et)[k]
    end
    return acc
end
function loop_eph_ad(eph, epochs, center, targ)
    acc = 0.0
    for et in epochs
        acc += ForwardDiff.derivative(t -> ephem_vector3(eph, center, targ, t), et)[1]
    end
    return acc
end
function loop_cspice_gps(epochs, targ, ref, obs)
    acc = 0.0
    for et in epochs
        acc += SPICE.spkgps(targ, et, ref, obs)[1][1]
    end
    return acc
end
function loop_cspice_geo(epochs, targ, ref, obs)
    acc = 0.0
    for et in epochs
        acc += SPICE.spkgeo(targ, et, ref, obs)[1][4]
    end
    return acc
end
function loop_alm_pxform(ctx, epochs, from, to)
    acc = 0.0
    for et in epochs
        acc += pxform(ctx, from, to, et)[3, 1]
    end
    return acc
end
function loop_eph_matrix(eph, epochs, from, to)
    acc = 0.0
    for et in epochs
        acc += pck_matrix(ephem_rotation3(eph, from, to, et))[3, 1]
    end
    return acc
end
function loop_cspice_pxform(epochs, from, to)
    acc = 0.0
    for et in epochs
        acc += SPICE.pxform(from, to, et)[3, 1]
    end
    return acc
end

"Rotation of the binary PCK type 2 angles `(phi, delta, w)`: `R3(w) R1(delta) R3(phi)`."
pck_matrix(a) = eul2m(a[3], a[2], a[1], 3, 1, 3)

load_almanac(path) = furnsh(SpiceContext(), path)
load_ephemerides(path) = EphemerisProvider(path)
load_cspice(path) = (SPICE.furnsh(path); SPICE.unload(path))

"100 epochs evenly spaced inside `[t0, t1]`, ends excluded."
inner_epochs(t0, t1) = [t0 + (t1 - t0) * (k + 0.5) / 100 for k in 0:99]

"Frame ID of the segment of `c`, from its descriptor."
function segment_frame(ctx, c::TypeCase)
    found = spksfs(ctx, c.body, (c.t0 + c.t1) / 2)
    @assert found !== nothing "no segment for $(c.label)"
    _, ic = dafus(found[2], 2, 6)
    return Int(ic[3])
end

"""
    max_reldiff(state, c, frame, epochs) -> Float64

Largest relative difference of `state(et)` (a 6-vector) from the CSPICE state of the
case `c` in `frame` over `epochs`: the larger of the position and velocity differences,
each over the size of its vector. `Inf` if `state` throws.
"""
function max_reldiff(state, c::TypeCase, frame::Integer, epochs)
    name = SPICE.frmnam(frame)
    worst = 0.0
    for et in epochs
        ref = SPICE.spkgeo(c.body, et, name, EARTH)[1]
        s = try
            collect(state(et))
        catch
            return Inf
        end
        dr = norm(s[1:3] .- ref[1:3]) / norm(ref[1:3])
        dv = norm(s[4:6] .- ref[4:6]) / norm(ref[4:6])
        worst = max(worst, dr, dv)
    end
    return worst
end

"First line of the error that `f()` throws, or `nothing`."
function first_error(f)
    try
        f()
    catch e
        return first(split(sprint(showerror, e), '\n'))
    end
    return nothing
end

"""
Speed-up of Almanac over Ephemerides.jl for each group that has both: the fastest
Ephemerides.jl row divided by the fastest Almanac row, minimum times.
"""
function ratio_table(rows::Vector{Row})
    io = IOBuffer()
    println(io, "| Case | Almanac min ns | Ephemerides.jl min ns | speed-up of Almanac |")
    println(io, "|---|---:|---:|---:|")
    for g in unique(r.group for r in rows)
        a = [r.min_ns for r in rows if r.group == g && engine(r) == "Almanac"]
        e = [r.min_ns for r in rows if r.group == g && engine(r) == "Ephemerides"]
        (isempty(a) || isempty(e)) && continue
        @printf(io, "| %s | %.1f | %.1f | %.2f |\n", g, minimum(a), minimum(e), minimum(e) / minimum(a))
    end
    return String(take!(io))
end

function main()
    @assert Threads.nthreads() == 1 "run on one thread: Ephemerides.jl caches by thread ID"
    typekernel = joinpath(mktempdir(), "spktypes.bsp")
    cases = write_type_kernel(typekernel, MersenneTwister(20261126))
    paths = (Kernels.DE440S, Kernels.MOON_TF, Kernels.MOON_BPC, typekernel)

    SPICE.kclear()
    foreach(SPICE.furnsh, paths)
    ctx = SpiceContext()
    foreach(p -> furnsh(ctx, p), paths)
    eph_de = EphemerisProvider(Kernels.DE440S)
    eph_moon = EphemerisProvider(Kernels.MOON_BPC)
    eph_types = EphemerisProvider(typekernel)

    # Agreement of every case of the type kernel, in the frame of its segment.
    agree = map(cases) do c
        frame = segment_frame(ctx, c)
        special = length(c.epochs) > 400 ? c.epochs[1:(length(c.epochs) ÷ 400 + 1):end] : c.epochs
        epochs = vcat(special, inner_epochs(c.t0, c.t1))
        alm = max_reldiff(et -> spkgeo(ctx, c.body, et, frame, EARTH)[1], c, frame, epochs)
        eph = max_reldiff(et -> ephem_vector6(eph_types, EARTH, c.body, et), c, frame, epochs)
        err = eph == Inf ? first_error(() -> ephem_vector6(eph_types, EARTH, c.body, epochs[1])) : nothing
        (; c, frame, alm, eph, err)
    end
    @assert all(a -> a.alm <= a.c.rtol, agree) "Almanac differs from CSPICE"

    epochs = anise_epochs()
    n = length(epochs)

    # de440s: the jets of the two engines, order by order, against each other, and the
    # orders 0 and 1 against CSPICE.
    jet_diff = zeros(4)
    for et in epochs
        a = spkder(ctx, EARTH, et, J2000, EMB, Val(3))
        e = ephem_vector12(eph_de, EMB, EARTH, et)
        for k in 0:3
            ek = SVector(e[3k + 1], e[3k + 2], e[3k + 3])
            jet_diff[k + 1] = max(jet_diff[k + 1], norm(ek - a[k + 1]) / norm(a[k + 1]))
        end
        c = SPICE.spkgeo(EARTH, et, "J2000", EMB)[1]
        @assert norm(a[1] - c[1:3]) / norm(c[1:3]) <= 1.0e-15 && norm(a[2] - c[4:6]) / norm(c[4:6]) <= 1.0e-15
    end
    rot_diff = zeros(2)
    for et in epochs
        c = SPICE.pxform("J2000", "MOON_PA_DE440", et)
        rot_diff[1] = max(rot_diff[1], maximum(abs.(pxform(ctx, J2000, MOON_PA, et) .- c)))
        rot_diff[2] = max(rot_diff[2], maximum(abs.(pck_matrix(ephem_rotation3(eph_moon, J2000, MOON_PA, et)) .- c)))
    end

    rows = Row[]
    g = "de440s, order 0"
    push!(rows, bench_row(g, "Almanac spkgps", loop_alm_gps, ctx, epochs, EARTH, J2000, EMB; n))
    q_earth = state_query(ctx, EARTH, J2000, "NONE", EMB)
    push!(rows, bench_row(g, "Almanac spkpos of a state_query", loop_alm_query_pos, q_earth, epochs; n))
    push!(rows, bench_row(g, "Ephemerides ephem_vector3", loop_eph, ephem_vector3, eph_de, epochs, EMB, EARTH, 1; n))
    push!(rows, bench_row(g, "CSPICE spkgps", loop_cspice_gps, epochs, EARTH, "J2000", EMB; n))
    g = "de440s, order 1"
    push!(rows, bench_row(g, "Almanac spkgeo", loop_alm_geo, ctx, epochs, EARTH, J2000, EMB; n))
    push!(rows, bench_row(g, "Almanac state_query", loop_alm_query, q_earth, epochs; n))
    push!(rows, bench_row(g, "Almanac spkder of a state_query, Val(1)", loop_alm_query_der, q_earth, epochs; n))
    push!(rows, bench_row(g, "Ephemerides ephem_vector6", loop_eph, ephem_vector6, eph_de, epochs, EMB, EARTH, 4; n))
    push!(rows, bench_row(g, "CSPICE spkgeo", loop_cspice_geo, epochs, EARTH, "J2000", EMB; n))
    g = "de440s, Moon relative to Earth (two links), order 1"
    push!(rows, bench_row(g, "Almanac spkgeo", loop_alm_geo, ctx, epochs, MOON, J2000, EARTH; n))
    push!(rows, bench_row(g, "Almanac state_query", loop_alm_query, state_query(ctx, MOON, J2000, "NONE", EARTH), epochs; n))
    push!(rows, bench_row(g, "CSPICE spkgeo", loop_cspice_geo, epochs, MOON, "J2000", EARTH; n))
    for (N, f) in ((2, ephem_vector9), (3, ephem_vector12))
        g = "de440s, order $N"
        push!(rows, bench_row(g, "Almanac spkder Val($N)", loop_alm_der, ctx, epochs, EARTH, J2000, EMB, Val(N); n))
        push!(rows, bench_row(g, "Ephemerides ephem_vector$(3N + 3)", loop_eph, f, eph_de, epochs, EMB, EARTH, 3N + 1; n))
    end
    g = "de440s, ForwardDiff of the position"
    push!(rows, bench_row(g, "Almanac spkgps, Dual epoch", loop_alm_ad, ctx, epochs, EARTH, J2000, EMB; n))
    push!(rows, bench_row(g, "Ephemerides ephem_vector3, Dual epoch", loop_eph_ad, eph_de, epochs, EMB, EARTH; n))

    timed = TypeCase[]
    for dtype in sort(unique(c.dtype for c in cases))
        i = findfirst(a -> a.c.dtype == dtype && a.frame == J2000, agree)
        i === nothing && continue
        a = agree[i]
        push!(timed, a.c)
        c = a.c
        ets = inner_epochs(c.t0, c.t1)
        g = "type $dtype, order 1 ($(c.label))"
        push!(rows, bench_row(g, "Almanac spkgeo", loop_alm_geo, ctx, ets, c.body, J2000, EARTH; n))
        q = state_query(ctx, c.body, J2000, "NONE", EARTH)
        push!(rows, bench_row(g, "Almanac state_query", loop_alm_query, q, ets; n))
        push!(rows, bench_row(g, "Almanac spkder of a state_query, Val(1)", loop_alm_query_der, q, ets; n))
        a.eph <= EPH_RTOL && push!(
            rows,
            bench_row(g, "Ephemerides ephem_vector6", loop_eph, ephem_vector6, eph_types, ets, EARTH, c.body, 4; n)
        )
        push!(rows, bench_row(g, "CSPICE spkgeo", loop_cspice_geo, ets, c.body, "J2000", EARTH; n))
    end

    g = "PCK type 2, J2000 to MOON_PA_DE440"
    push!(rows, bench_row(g, "Almanac pxform, names", loop_alm_pxform, ctx, epochs, "J2000", "MOON_PA_DE440"; n))
    push!(rows, bench_row(g, "Almanac pxform, frame IDs", loop_alm_pxform, ctx, epochs, J2000, MOON_PA; n))
    q_moon = rotation_query(ctx, J2000, MOON_PA)
    push!(rows, bench_row(g, "Almanac rotation_query", loop_alm_query_rot, q_moon, epochs; n))
    push!(rows, bench_row(g, "Ephemerides ephem_rotation3 and eul2m", loop_eph_matrix, eph_moon, epochs, J2000, MOON_PA; n))
    push!(rows, bench_row(g, "CSPICE pxform", loop_cspice_pxform, epochs, "J2000", "MOON_PA_DE440"; n))

    g = "load de440s.bsp"
    push!(rows, bench_row(g, "Almanac furnsh, new context", load_almanac, Kernels.DE440S; n = 1))
    push!(rows, bench_row(g, "Ephemerides EphemerisProvider", load_ephemerides, Kernels.DE440S; n = 1))
    SPICE.kclear()
    push!(rows, bench_row(g, "CSPICE furnsh and unload", load_cspice, Kernels.DE440S; n = 1))

    io = IOBuffer()
    println(io, "| Case | frame | Almanac | Ephemerides.jl | timed |")
    println(io, "|---|---|---:|---:|---|")
    for a in agree
        eph = a.err === nothing ? @sprintf("%.1e", a.eph) : "throws `$(a.err)`"
        how = !(a.c in timed) ? "" : a.eph <= EPH_RTOL ? "all three" : "without Ephemerides.jl"
        @printf(io, "| %s | %s | %.1e | %s | %s |\n", a.c.label, SPICE.frmnam(a.frame), a.alm, eph, how)
    end
    agree_text = String(take!(io))

    text = run_header("Benchmark against Ephemerides.jl") * """
        Ephemerides.jl $(pkgversion(Ephemerides)), Almanac $(pkgversion(Almanac)), SPICE.jl $(pkgversion(SPICE)).
        Kernels: `$(Kernels.DE440S)`, `$(Kernels.MOON_BPC)` with `$(Kernels.MOON_TF)`, and the
        SPK kernel of the test group `spktypes`, written for this run.

        ## Agreement

        de440s, the Earth relative to the Earth-Moon barycenter, largest relative difference
        of Ephemerides.jl from Almanac by order of the jet: $(join((@sprintf("%.1e", d) for d in jet_diff), ", ")).
        Almanac agrees with CSPICE to 1e-15 at orders 0 and 1.

        Rotation from J2000 to MOON_PA_DE440, largest difference of an element from CSPICE:
        Almanac $(@sprintf("%.1e", rot_diff[1])), Ephemerides.jl with `eul2m` $(@sprintf("%.1e", rot_diff[2])).

        SPK types: largest relative difference from CSPICE in the frame of the segment,
        at the knots and record bounds of the case and at 100 epochs inside it. A case is
        timed with Ephemerides.jl when its difference is at most $(EPH_RTOL).

        """ * agree_text * "\n## Time per query\n\n" * markdown_table(rows) *
        "\n## Almanac against Ephemerides.jl\n\n" * ratio_table(rows)
    println()
    println(text)
    write_results("ephemerides", text, rows)
    SPICE.kclear()
    return rows
end

main()
