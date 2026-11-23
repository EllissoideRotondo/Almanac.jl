# The Almanac.jl side of the speed comparison of PEPPER with ANISE.
#
# PEPPER (https://github.com/jacobwilliams/pepper) compares itself with ANISE in
# `test/core/bench_anise.f90` and `bench/anise/src/main.rs`. This script makes the same
# queries at the same epochs (100000 Park-Miller epochs over 20 years from J2000, seed
# 4242, and 100000 epochs 60 s apart) and prints the same lines,
# `key|n|microseconds per operation|checksum`, so that the three programs can be run on
# one machine and compared row by row. Each time is the fastest of `ALMANAC_BENCH_RUNS`
# runs (default 7), as in PEPPER's `bench/compare_anise.py`.
#
#     julia --project=benchmark benchmark/pepper.jl KERNEL_DIR
#
# KERNEL_DIR holds naif0012.tls, pck00011.tpc, de421.bsp, earth_latest_high_prec.bpc,
# moon_pa_de440_200625.bpc and moon_de440_250416.tf, the kernels of PEPPER's comparison.
#
# Keys without a suffix call the exported functions with names and strings, as a program
# written for CSPICE does. Keys ending in `_ids` call them on the `KernelSet` of the
# context with NAIF codes and a parsed `Aberration`, which is the nearest Almanac.jl has
# to PEPPER's prepared queries (`prepare_state`/`query_state`). Keys ending in `_table`
# evaluate an `EphemerisTable` or `RotationTable` built over the epochs; a table is an
# interpolation, so its checksum agrees with the others only to the error of the table.

using Almanac
using LinearAlgebra, Printf

const NQ = 100_000
const NSETS = 10
const YEAR = 365.25 * 86400.0
const STEP = 60.0
const RUNS = parse(Int, get(ENV, "ALMANAC_BENCH_RUNS", "7"))
const KERNELS = (
    "naif0012.tls", "pck00011.tpc", "de421.bsp", "earth_latest_high_prec.bpc",
    "moon_pa_de440_200625.bpc", "moon_de440_250416.tf",
)

"The epochs of PEPPER's `random_epochs`: span*(random(1000000)-1)/1e6, Park-Miller, seed 4242."
function random_epochs(n, span)
    seed = Int64(4242)
    return map(1:n) do _
        seed = (seed * 48271) % 2147483647
        span * ((seed % 1_000_000 + 1) - 1) / 1.0e6
    end
end

report(key, n, seconds, checksum) =
    @printf("%s|%d|%.5e|%.15e\n", key, n, 1.0e6 * seconds / n, checksum)

"Fastest of `RUNS` timings of `f(ets)`, which returns a checksum; the first call compiles."
function best(f, ets)
    sum = f(ets)
    t = Inf
    for _ in 1:RUNS
        t0 = time_ns()
        sum = f(ets)
        t = min(t, (time_ns() - t0) / 1.0e9)
    end
    return t, sum
end

function run(key, f, ets)
    t, s = best(f, ets)
    return report(key, length(ets), t, s)
end

# --- with names and strings, as from a CSPICE program ------------------------------------

function states_names(ctx, ets, sfx)
    run("spkezr_moon$sfx", e -> sum(et -> spkezr(ctx, "MOON", et, "J2000", "NONE", "EARTH")[1][1], e), ets)
    run("spkgeo_moon$sfx", e -> sum(et -> spkgeo(ctx, 301, et, "J2000", 399)[1][1], e), ets)
    run("spkezr_mars_lts$sfx", e -> sum(et -> spkezr(ctx, "MARS BARYCENTER", et, "ECLIPJ2000", "LT+S", "EARTH")[1][1], e), ets)
    run(
        "spkezr_mars_lt$sfx", e -> sum(e) do et
            s = spkezr(ctx, "MARS BARYCENTER", et, "ECLIPJ2000", "LT", "EARTH")[1]
            s[1] + s[4]
        end, ets
    )
    run("spkpos_mars_lts$sfx", e -> sum(et -> spkpos(ctx, "MARS BARYCENTER", et, "ECLIPJ2000", "LT+S", "EARTH")[1][1], e), ets)
    return run("spkpos_sun_cns$sfx", e -> sum(et -> spkpos(ctx, "SUN", et, "J2000", "CN+S", "EARTH")[1][1], e), ets)
end

function frames_names(ctx, ets, sfx)
    run("pxform_itrf93$sfx", e -> sum(et -> pxform(ctx, "J2000", "ITRF93", et)[1, 1], e), ets)
    run("sxform_itrf93$sfx", e -> sum(et -> sxform(ctx, "J2000", "ITRF93", et)[4, 1], e), ets)
    run("pxform_iau_mars$sfx", e -> sum(et -> pxform(ctx, "J2000", "IAU_MARS", et)[1, 1], e), ets)
    return run("pxform_moon_me$sfx", e -> sum(et -> pxform(ctx, "MOON_ME", "J2000", et)[1, 1], e), ets)
end

# --- on the snapshot, with codes and a parsed correction ------------------------------

function states_ids(ks, ets, sfx)
    none, lts, lt, cns = Aberration("NONE"), Aberration("LT+S"), Aberration("LT"), Aberration("CN+S")
    j2000, eclip = 1, Almanac.frame_id(ks, "ECLIPJ2000")
    run("spkezr_moon$sfx", e -> sum(et -> spkezr(ks, 301, et, j2000, none, 399)[1][1], e), ets)
    run("spkgeo_moon$sfx", e -> sum(et -> spkgeo(ks, 301, et, j2000, 399)[1][1], e), ets)
    run("spkezr_mars_lts$sfx", e -> sum(et -> spkezr(ks, 4, et, eclip, lts, 399)[1][1], e), ets)
    run(
        "spkezr_mars_lt$sfx", e -> sum(e) do et
            s = spkezr(ks, 4, et, eclip, lt, 399)[1]
            s[1] + s[4]
        end, ets
    )
    run("spkpos_mars_lts$sfx", e -> sum(et -> spkpos(ks, 4, et, eclip, lts, 399)[1][1], e), ets)
    return run("spkpos_sun_cns$sfx", e -> sum(et -> spkpos(ks, 10, et, j2000, cns, 399)[1][1], e), ets)
end

function frames_ids(ks, ets, sfx)
    itrf, mars, moon = (Almanac.frame_id(ks, n) for n in ("ITRF93", "IAU_MARS", "MOON_ME"))
    run("pxform_itrf93$sfx", e -> sum(et -> pxform(ks, 1, itrf, et)[1, 1], e), ets)
    run("sxform_itrf93$sfx", e -> sum(et -> sxform(ks, 1, itrf, et)[4, 1], e), ets)
    run("pxform_iau_mars$sfx", e -> sum(et -> pxform(ks, 1, mars, et)[1, 1], e), ets)
    return run("pxform_moon_me$sfx", e -> sum(et -> pxform(ks, moon, 1, et)[1, 1], e), ets)
end

# --- tables over the epochs -----------------------------------------------------------

function tables(ctx, ets, sfx)
    t0, t1 = extrema(ets)
    t1 += 1.0
    moon = ephemeris_table(ctx, "MOON", t0, t1, "J2000", "NONE", "EARTH")
    run("spkezr_moon$sfx", e -> sum(et -> table_state(moon, et)[1], e), ets)
    itrf = rotation_table(ctx, "J2000", "ITRF93", t0, t1)
    run("pxform_itrf93$sfx", e -> sum(et -> table_rotation(itrf, et)[1, 1], e), ets)
    mars = rotation_table(ctx, "J2000", "IAU_MARS", t0, t1)
    return run("pxform_iau_mars$sfx", e -> sum(et -> table_rotation(mars, et)[1, 1], e), ets)
end

function main(kdir)
    ets = random_epochs(NQ, 20 * YEAR)
    seq = [i * STEP for i in 0:(NQ - 1)]
    paths = [joinpath(kdir, k) for k in KERNELS]

    ctx = SpiceContext()
    load = function (_)
        kclear(ctx)
        foreach(p -> furnsh(ctx, p), paths)
        return 0.0
    end
    t, _ = best(load, 1:NSETS)          # one load per call; report per set
    report("load_set", 1, t, 0.0)

    states_names(ctx, ets, "")
    states_names(ctx, seq, "_seq")
    frames_names(ctx, ets, "")
    frames_names(ctx, seq, "_seq")

    ks = Almanac.kernelset(ctx)
    states_ids(ks, ets, "_ids")
    states_ids(ks, seq, "_ids_seq")
    frames_ids(ks, ets, "_ids")
    frames_ids(ks, seq, "_ids_seq")

    tables(ctx, ets, "_table")
    return tables(ctx, seq, "_table_seq")
end

isempty(ARGS) ? println("usage: julia --project=benchmark benchmark/pepper.jl KERNEL_DIR") : main(ARGS[1])
