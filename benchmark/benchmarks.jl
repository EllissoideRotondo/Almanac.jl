# The benchmark suite of AirspeedVelocity.jl: `SUITE` times the hot paths of Almanac, and
# `.github/workflows/benchmark_pr.yml` runs it on a pull request and on `main` and
# comments the ratios on the pull request. The script of the pull request runs on both
# revisions, so a benchmark must call only what `main` has too.
#
#     julia -e 'using Pkg; Pkg.add("AirspeedVelocity")'
#     ~/.julia/bin/benchpkg Almanac --rev=main,dirty --path=.        # from the root
#
# Each benchmark is a loop over 100 epochs, so one sample takes microseconds and needs
# no tuning; the times are per loop. The group "cspice" times CSPICE (SPICE.jl), which
# is the same on both revisions: its ratios show the noise of the runner. The comparison
# with CSPICE call by call is in the scripts beside this file (`run.jl`).
#
# The kernels are those of the tests, in ALMANAC_KERNEL_DIR, or else in a scratch
# directory of the depot; they are downloaded when missing.

haskey(ENV, "ALMANAC_KERNEL_DIR") ||
    (ENV["ALMANAC_KERNEL_DIR"] = joinpath(first(DEPOT_PATH), "scratchspaces", "almanac-benchmark-kernels"))
include(joinpath(@__DIR__, "common.jl"))         # `Kernels`, `anise_epochs`
using Almanac
import ForwardDiff, SPICE

let gone = Kernels.fetch_kernels(stderr) ∩ collect(Kernels.FULL_SET_FILES)
    isempty(gone) || error("kernels missing for the benchmarks: ", join(gone, ", "))
end

const SUITE = BenchmarkGroup()
const SECONDS = 1.0                 # per benchmark; two revisions of about 50 benchmarks

"Sum of one number of `f(x)` over `xs`, so the calls are not removed."
function sumloop(f::F, xs) where {F}
    acc = 0.0
    for x in xs
        acc += f(x)
    end
    return acc
end

"A benchmark of `sumloop(f, xs)`."
bench(f, xs) = @benchmarkable sumloop($f, $xs) seconds = SECONDS

const SET = Kernels.full_set()
const CTX = SpiceContext()
foreach(p -> furnsh(CTX, p), SET)
const KS = Almanac.kernelset(CTX)
SPICE.kclear(); foreach(SPICE.furnsh, SET)

# 100 epochs from 1900 to 2099, as in the benchmark of ANISE, and the same in random order.
const ETS = anise_epochs()
const SHUFFLED = shuffle(MersenneTwister(20261011), ETS)
const EARTH, MOON, SUN, J2000 = 399, 301, 10, 1
const GM_EARTH = 398600.435436
const DAY = 86400.0

let g = SUITE["spk"] = BenchmarkGroup()
    g["spkezr ids"] = bench(et -> spkezr(CTX, EARTH, et, J2000, "NONE", MOON)[1][1], ETS)
    g["spkezr ids, random order"] = bench(et -> spkezr(CTX, EARTH, et, J2000, "NONE", MOON)[1][1], SHUFFLED)
    g["spkezr names"] = bench(et -> spkezr(CTX, "EARTH", et, "J2000", "NONE", "MOON")[1][1], ETS)
    g["spkezr KernelSet"] = bench(et -> spkezr(KS, EARTH, et, J2000, "NONE", MOON)[1][1], ETS)
    g["spkezr Moon rel. Sun"] = bench(et -> spkezr(CTX, MOON, et, J2000, "NONE", SUN)[1][1], ETS)
    g["spkgeo"] = bench(et -> spkgeo(CTX, EARTH, et, J2000, MOON)[1][1], ETS)
    g["spkder Val(2)"] = bench(et -> spkder(CTX, EARTH, et, J2000, MOON, Val(2))[3][1], ETS)
end

let g = SUITE["aberration"] = BenchmarkGroup()
    for ab in ("LT", "LT+S", "CN+S")
        g["spkezr $ab"] = bench(et -> spkezr(CTX, EARTH, et, J2000, ab, MOON)[1][1], ETS)
    end
    g["spkpos CN+S"] = bench(et -> spkpos(CTX, EARTH, et, J2000, "CN+S", MOON)[1][1], ETS)
    g["spkder LT+S Val(2)"] = bench(et -> spkder(CTX, EARTH, et, J2000, "LT+S", MOON, Val(2))[3][1], ETS)
end

let g = SUITE["frames"] = BenchmarkGroup()
    for to in ("IAU_EARTH", "MOON_PA")
        id = namfrm(CTX, to)
        g["pxform $to ids"] = bench(et -> pxform(CTX, J2000, id, et)[1], ETS)
        g["sxform $to ids"] = bench(et -> sxform(CTX, J2000, id, et)[1], ETS)
        g["pxder $to Val(2)"] = bench(et -> pxder(CTX, J2000, id, et, Val(2))[3][1], ETS)
    end
    g["pxform IAU_EARTH names"] = bench(et -> pxform(CTX, "J2000", "IAU_EARTH", et)[1], ETS)
end

let g = SUITE["time"] = BenchmarkGroup()
    strings = [et2utc(CTX, et, "ISOC", 3) for et in ETS]
    pict = "YYYY-MON-DD HR:MN:SC.### ::TDB"
    g["str2et"] = bench(s -> str2et(CTX, s), strings)
    g["et2utc"] = bench(et -> Float64(length(et2utc(CTX, et, "C", 3))), ETS)
    g["timout"] = bench(et -> Float64(length(timout(CTX, et, pict))), ETS)
    g["unitim"] = bench(et -> unitim(CTX, et, "TDB", "JDTDB"), ETS)
end

let g = SUITE["pool"] = BenchmarkGroup()
    g["bodn2c"] = bench(_ -> Float64(bodn2c(CTX, "MOON")), ETS)
    g["namfrm"] = bench(_ -> Float64(namfrm(CTX, "IAU_EARTH")), ETS)
    g["bodvrd"] = bench(_ -> bodvrd(CTX, "EARTH", "RADII")[1], ETS)
end

let g = SUITE["conics"] = BenchmarkGroup()
    pairs = [(spkezr(CTX, MOON, et, J2000, "NONE", EARTH)[1], et) for et in ETS]
    epairs = [(oscelt(s, et, GM_EARTH), et) for (s, et) in pairs]
    g["oscelt"] = bench(((s, et),) -> oscelt(s, et, GM_EARTH)[2], pairs)
    g["conics"] = bench(((e, et),) -> conics(e, et)[1], epairs)
end

let g = SUITE["geometry"] = BenchmarkGroup()
    for method in ("NEAR POINT/ELLIPSOID", "INTERCEPT/ELLIPSOID")
        g["subpnt $method"] = bench(et -> subpnt(CTX, method, "EARTH", et, "IAU_EARTH", "LT+S", "MOON")[1][1], ETS)
    end
    rays = [(et, -spkpos(CTX, EARTH, et, J2000, "NONE", MOON)[1]) for et in ETS]
    g["sincpt"] = bench(((et, d),) -> sincpt(CTX, "ELLIPSOID", "EARTH", et, "IAU_EARTH", "CN+S", "MOON", "J2000", d)[1][1], rays)
    spoints = [(et, subpnt(CTX, "NEAR POINT/ELLIPSOID", "EARTH", et, "IAU_EARTH", "LT+S", "MOON")[1]) for et in ETS]
    g["ilumin"] = bench(((et, p),) -> ilumin(CTX, "ELLIPSOID", "EARTH", et, "IAU_EARTH", "LT+S", "MOON", p)[3], spoints)
    g["occult"] = bench(et -> Float64(occult(CTX, "MOON", "ELLIPSOID", "IAU_MOON", "SUN", "ELLIPSOID", "IAU_SUN", "CN", "EARTH", et)), ETS)
end

let g = SUITE["gf"] = BenchmarkGroup()
    t0 = ETS[60]; t1 = t0 + 365.25 * DAY
    g["gfdist LOCMAX one year"] = @benchmarkable(
        gfdist(CTX, "MOON", "NONE", "EARTH", "LOCMAX", 0.0, 0.0, 5 * DAY, Window([$t0, $t1])),
        seconds = SECONDS
    )
end

let g = SUITE["kernels"] = BenchmarkGroup()
    g["furnsh standard set"] = @benchmarkable(
        foreach(p -> furnsh(c, p), SET), setup = (c = SpiceContext()), evals = 1, seconds = SECONDS
    )
end

let g = SUITE["tables"] = BenchmarkGroup()
    t0 = ETS[60]; t1 = t0 + 30 * DAY
    span = collect(range(t0, t1; length = 100))
    moon = ephemeris_table(CTX, "MOON", t0, t1, "J2000", "NONE", "EARTH")
    earth = rotation_table(CTX, "J2000", "IAU_EARTH", t0, t1)
    g["ephemeris_table 30 days"] = @benchmarkable(
        ephemeris_table(CTX, "MOON", $t0, $t1, "J2000", "NONE", "EARTH"), seconds = SECONDS
    )
    g["table_state"] = bench(et -> table_state(moon, et)[1], span)
    g["table_rotation"] = bench(et -> table_rotation(earth, et)[1], span)
end

let g = SUITE["ad"] = BenchmarkGroup()
    g["ForwardDiff spkpos"] = bench(
        et -> ForwardDiff.derivative(t -> spkpos(CTX, EARTH, t, J2000, "NONE", MOON)[1], et)[1], ETS
    )
end

# The same on both revisions: the ratios of this group are the noise of the runner.
let g = SUITE["cspice"] = BenchmarkGroup()
    g["spkezr"] = bench(et -> SPICE.spkezr("EARTH", et, "J2000", "NONE", "MOON")[1][1], ETS)
    g["pxform"] = bench(et -> SPICE.pxform("J2000", "IAU_EARTH", et)[1], ETS)
end
