# State benchmark: Almanac against CSPICE (SPICE.jl) against ANISE (anise_jll).
#
#     julia --startup-file=no --threads=8 --project=benchmark benchmark/states.jl
#
# ANISE is optional: `anise_jll` is a local build and is not in `benchmark/Project.toml`.
# Without it the ANISE rows are left out; `benchmark/anise.jl` says how to add it.
#
# The main case mirrors the benchmark of ANISE itself: 100 epochs evenly spaced from
# 1900-01-01 12:00 ET to 2099-01-01 12:00 ET (end excluded) on de440s.bsp, the Earth
# relative to the Moon, J2000, no aberration correction. Each engine loads de440s.bsp
# only. Before any timing the three engines are compared on the same queries, so the
# table compares the same computation.
#
# Further cases: the same epochs in random order (defeats a cache of the last record),
# 100 epochs 60 s apart (one record of the kernel, so CSPICE reads the file once and
# then answers from the record it keeps: the case least favourable to Almanac),
# the Moon relative to the Sun (longer chain: 301 -> 3 -> 0 <- 10), and the throughput
# of Almanac on all threads (through a context and on a `KernelSet` snapshot). CSPICE is
# not thread safe and is never called from a task.
# Set ALMANAC_BENCH_ANISE_THREADS=1 to add the threaded throughput of ANISE.

include("common.jl")
include("anise.jl")
using Almanac
import SPICE

require_kernels(Kernels.DE440S)

const EARTH, MOON, SUN, J2000 = 399, 301, 10, 1

# One loop per call form. Every loop sums one component, so the calls are not removed.
function loop_neo_names(ctx, epochs)
    acc = 0.0
    for et in epochs
        s, _ = spkezr(ctx, "EARTH", et, "J2000", "NONE", "MOON")
        acc += s[1]
    end
    return acc
end
function loop_neo_ids(ctx, epochs, targ, obs)
    acc = 0.0
    for et in epochs
        s, _ = spkezr(ctx, targ, et, J2000, "NONE", obs)
        acc += s[1]
    end
    return acc
end
function loop_neo_default(epochs)
    acc = 0.0
    for et in epochs
        s, _ = spkezr(EARTH, et, J2000, "NONE", MOON)
        acc += s[1]
    end
    return acc
end
function loop_neo_ks(ks, epochs)
    acc = 0.0
    for et in epochs
        s, _ = spkezr(ks, EARTH, et, J2000, "NONE", MOON)
        acc += s[1]
    end
    return acc
end
function loop_neo_der(ctx, epochs, v::Val{N}) where {N}
    acc = 0.0
    for et in epochs
        jet = spkder(ctx, EARTH, et, J2000, MOON, v)
        acc += jet[N + 1][1]
    end
    return acc
end
function loop_cspice_ezr(epochs, targ, obs)
    acc = 0.0
    for et in epochs
        s, _ = SPICE.spkezr(targ, et, "J2000", "NONE", obs)
        acc += s[1]
    end
    return acc
end
function loop_cspice_geo(epochs, targ, obs)
    acc = 0.0
    for et in epochs
        s, _ = SPICE.spkgeo(targ, et, "J2000", obs)
        acc += s[1]
    end
    return acc
end
function loop_anise(alm, epochs, targ, obs, out)
    acc = 0.0
    for et in epochs
        anise_spk_ezr(alm, targ, et, J2000, obs, out)
        acc += out[].position_km[1]
    end
    return acc
end

"""
    check_agreement(ctx, alm, epochs, targ, obs) -> NamedTuple

Largest absolute differences (km, km/s) of Almanac and of ANISE from CSPICE over
`epochs`. Throws if a difference is above 1e-4 km or 1e-9 km/s: the engines would then
not compute the same thing.
"""
function check_agreement(ctx, alm, epochs, targ, obs)
    dn = 0.0; da = 0.0; vn = 0.0; va = 0.0
    for et in epochs
        c, _ = SPICE.spkezr(string(targ), et, "J2000", "NONE", string(obs))
        n, _ = spkezr(ctx, targ, et, J2000, "NONE", obs)
        dn = max(dn, maximum(abs.(n[SOneTo(3)] .- c[1:3])))
        vn = max(vn, maximum(abs.(n[SVector(4, 5, 6)] .- c[4:6])))
        if alm !== nothing
            a = anise_state(alm, targ, et, J2000, obs)
            da = max(da, maximum(abs.(a[SOneTo(3)] .- c[1:3])))
            va = max(va, maximum(abs.(a[SVector(4, 5, 6)] .- c[4:6])))
        end
    end
    @assert dn <= 1.0e-4 && vn <= 1.0e-9 "Almanac differs from CSPICE: $dn km, $vn km/s"
    @assert da <= 1.0e-4 && va <= 1.0e-9 "ANISE differs from CSPICE: $da km, $va km/s"
    return (; neo_km = dn, neo_kms = vn, anise_km = da, anise_kms = va)
end

"Queries per second of `query(i)` for `total` queries on `nt` tasks (static chunks)."
function throughput(query, total::Int, nt::Int)
    @assert total > 0 && nt > 0
    chunk = total ÷ nt
    sums = zeros(nt)
    work(k) = begin
        acc = 0.0
        for i in 1:chunk
            acc += query(i + k)
        end
        sums[k] = acc
    end
    work(1)                                                # compile
    t = @elapsed begin
        if nt == 1
            work(1)
        else
            @sync for k in 1:nt
                Threads.@spawn work(k)
            end
        end
    end
    return chunk * nt / t
end

function main()
    epochs = anise_epochs()
    shuffled = shuffle(MersenneTwister(20261011), epochs)
    dense = [7.0e8 + 60.0 * k for k in 0:99]
    n = length(epochs)

    SPICE.kclear(); SPICE.furnsh(Kernels.DE440S)
    ctx = SpiceContext(); furnsh(ctx, Kernels.DE440S)
    kclear(); furnsh(Kernels.DE440S)
    ks = Almanac.kernelset(ctx)
    alm = try
        anise_almanac(Kernels.DE440S)
    catch e
        printstyled("ANISE is not available: ", sprint(showerror, e), "\n"; color = :yellow)
        nothing
    end
    out = Ref(AniseState())

    println("ANISE library: ", alm === nothing ? "none" : anise_version())
    agree = check_agreement(ctx, alm, epochs, EARTH, MOON)
    agree2 = check_agreement(ctx, alm, epochs, MOON, SUN)
    @printf(
        "agreement with CSPICE, Earth rel. Moon: Almanac %.2e km %.2e km/s; ANISE %.2e km %.2e km/s\n",
        agree.neo_km, agree.neo_kms, agree.anise_km, agree.anise_kms
    )
    @printf(
        "agreement with CSPICE, Moon rel. Sun:   Almanac %.2e km %.2e km/s; ANISE %.2e km %.2e km/s\n",
        agree2.neo_km, agree2.neo_kms, agree2.anise_km, agree2.anise_kms
    )

    rows = Row[]
    g = "Earth rel. Moon, ordered"
    push!(rows, bench_row(g, "Almanac spkezr, integer IDs, context", loop_neo_ids, ctx, epochs, EARTH, MOON; n))
    push!(rows, bench_row(g, "Almanac spkezr, names, context", loop_neo_names, ctx, epochs; n))
    push!(rows, bench_row(g, "Almanac spkezr, integer IDs, default context", loop_neo_default, epochs; n))
    push!(rows, bench_row(g, "Almanac spkezr, integer IDs, KernelSet", loop_neo_ks, ks, epochs; n))
    for N in 0:3
        push!(rows, bench_row(g, "Almanac spkder Val($N)", loop_neo_der, ctx, epochs, Val(N); n))
    end
    push!(rows, bench_row(g, "CSPICE spkezr (SPICE.jl, names)", loop_cspice_ezr, epochs, "EARTH", "MOON"; n))
    push!(rows, bench_row(g, "CSPICE spkgeo (SPICE.jl, integer IDs)", loop_cspice_geo, epochs, EARTH, MOON; n))
    alm === nothing ||
        push!(rows, bench_row(g, "ANISE spk_ezr (C interface)", loop_anise, alm, epochs, EARTH, MOON, out; n))

    g = "Earth rel. Moon, random order"
    push!(rows, bench_row(g, "Almanac spkezr, integer IDs, context", loop_neo_ids, ctx, shuffled, EARTH, MOON; n))
    push!(rows, bench_row(g, "CSPICE spkezr (SPICE.jl, names)", loop_cspice_ezr, shuffled, "EARTH", "MOON"; n))
    push!(rows, bench_row(g, "CSPICE spkgeo (SPICE.jl, integer IDs)", loop_cspice_geo, shuffled, EARTH, MOON; n))
    alm === nothing ||
        push!(rows, bench_row(g, "ANISE spk_ezr (C interface)", loop_anise, alm, shuffled, EARTH, MOON, out; n))

    g = "Earth rel. Moon, 60 s apart"
    push!(rows, bench_row(g, "Almanac spkezr, integer IDs, context", loop_neo_ids, ctx, dense, EARTH, MOON; n))
    push!(rows, bench_row(g, "CSPICE spkezr (SPICE.jl, names)", loop_cspice_ezr, dense, "EARTH", "MOON"; n))
    push!(rows, bench_row(g, "CSPICE spkgeo (SPICE.jl, integer IDs)", loop_cspice_geo, dense, EARTH, MOON; n))
    alm === nothing ||
        push!(rows, bench_row(g, "ANISE spk_ezr (C interface)", loop_anise, alm, dense, EARTH, MOON, out; n))

    g = "Moon rel. Sun, ordered"
    push!(rows, bench_row(g, "Almanac spkezr, integer IDs, context", loop_neo_ids, ctx, epochs, MOON, SUN; n))
    push!(rows, bench_row(g, "CSPICE spkezr (SPICE.jl, names)", loop_cspice_ezr, epochs, "MOON", "SUN"; n))
    push!(rows, bench_row(g, "CSPICE spkgeo (SPICE.jl, integer IDs)", loop_cspice_geo, epochs, MOON, SUN; n))
    alm === nothing ||
        push!(rows, bench_row(g, "ANISE spk_ezr (C interface)", loop_anise, alm, epochs, MOON, SUN, out; n))

    # Throughput on threads. Each task queries its own sequence of the 100 epochs.
    total = 4_000_000
    nt = Threads.nthreads()
    thr = IOBuffer()
    println(thr, "| Engine | threads | queries/s | speed-up |")
    println(thr, "|---|---:|---:|---:|")
    function throughput_rows(label, query)
        t1 = throughput(query, total, 1)
        @printf(thr, "| %s | 1 | %.3e | 1.00 |\n", label, t1)
        if nt > 1
            tn = throughput(query, total, nt)
            @printf(thr, "| %s | %d | %.3e | %.2f |\n", label, nt, tn, tn / t1)
        end
        return nothing
    end
    throughput_rows(
        "Almanac spkezr, integer IDs, context",
        i -> spkezr(ctx, EARTH, @inbounds(epochs[i % n + 1]), J2000, "NONE", MOON)[1][1]
    )
    # The same queries on the snapshot itself: no read of the context reference per query.
    throughput_rows(
        "Almanac spkezr, integer IDs, KernelSet",
        i -> spkezr(ks, EARTH, @inbounds(epochs[i % n + 1]), J2000, "NONE", MOON)[1][1]
    )
    if alm !== nothing && get(ENV, "ALMANAC_BENCH_ANISE_THREADS", "0") == "1"
        throughput_rows(
            "ANISE spk_ezr", function (i)
                o = Ref(AniseState())
                anise_spk_ezr(alm, EARTH, @inbounds(epochs[i % n + 1]), J2000, MOON, o)
                return o[].position_km[1]
            end
        )
    end
    nt == 1 && println("one thread only: run with --threads=8 for the throughput on threads")

    text = run_header("State benchmark: Almanac, CSPICE, ANISE") * """
        Kernel: `$(Kernels.DE440S)` only, in every engine. ANISE library $(alm === nothing ? "not available" : anise_version()).
        Epochs: 100, evenly spaced from 1900-01-01 12:00 ET to 2099-01-01 12:00 ET (end excluded).
        Frame J2000, no aberration correction.

        Largest difference from CSPICE on the timed queries:

        | Query | Almanac km | Almanac km/s | ANISE km | ANISE km/s |
        |---|---:|---:|---:|---:|
        | Earth rel. Moon | $(@sprintf("%.2e", agree.neo_km)) | $(@sprintf("%.2e", agree.neo_kms)) | $(@sprintf("%.2e", agree.anise_km)) | $(@sprintf("%.2e", agree.anise_kms)) |
        | Moon rel. Sun | $(@sprintf("%.2e", agree2.neo_km)) | $(@sprintf("%.2e", agree2.neo_kms)) | $(@sprintf("%.2e", agree2.anise_km)) | $(@sprintf("%.2e", agree2.anise_kms)) |

        ## Time per query

        """ * markdown_table(rows) * "\n## Throughput on threads ($(total) queries)\n\n" *
        String(take!(thr))
    println()
    println(text)
    write_results("states", text, rows)
    alm === nothing || anise_free(alm)
    SPICE.kclear(); kclear()
    return rows
end

main()
