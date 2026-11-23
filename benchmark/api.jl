# Benchmark of the calls that the other scripts do not time: Almanac against CSPICE
# (SPICE.jl), one group per call, so that every group ends in its CSPICE row.
#
#     julia --startup-file=no --project=benchmark benchmark/api.jl
#
# Time conversion (`str2et`, `et2utc`, `timout`, `unitim`), names and the kernel pool
# (`bodn2c`, `namfrm`, `bodvrd`), conics (`oscelt`, `conics`), geometry on the IAU
# ellipsoids (`subpnt`, `sincpt`, `ilumin`, `occult`), one search of the geometry finder
# (`gfdist`) and the load of the standard kernel set. Both engines load the standard set
# of the tests (`Kernels.full_set()`), and before any timing the two are compared on the
# timed inputs. The SPICE.jl calls convert strings and arrays on every call, as a Julia
# program that uses CSPICE does; that cost is part of their rows.

include("common.jl")
using Almanac
import SPICE

require_kernels(Kernels.FULL_SET_FILES...)

const GM_EARTH = 398600.435436
const DAY = 86400.0

# Each loop sums one number of each result, so the calls are not removed.
function loop(f, xs)
    acc = 0.0
    for x in xs
        acc += f(x)
    end
    return acc
end

"Largest absolute difference of `a` and `b` over the inputs; throws above `tol`."
function agreement(label, a, b, xs, tol)
    d = maximum(x -> maximum(abs.(a(x) .- b(x))), xs)
    @assert d <= tol "Almanac differs from CSPICE in $label: $d"
    return @sprintf("| %s | %.2e |", label, d)
end

"The intervals of a window of either engine, as a flat vector."
window_values(w) = [v for i in 1:wncard(w) for v in wnfetd(w, i)]
cspice_window_values(w) = [v for i in 1:SPICE.wncard(w) for v in SPICE.wnfetd(w, i)]

function main()
    epochs = anise_epochs()
    n = length(epochs)
    set = Kernels.full_set()
    SPICE.kclear(); foreach(SPICE.furnsh, set)
    ctx = SpiceContext(); foreach(p -> furnsh(ctx, p), set)

    # Inputs that the timed calls share.
    strings = [SPICE.et2utc(et, "ISOC", 3) for et in epochs]
    moon = [spkezr(ctx, 301, et, 1, "NONE", 399)[1] for et in epochs]
    elts = [oscelt(s, et, GM_EARTH) for (s, et) in zip(moon, epochs)]
    pairs = collect(zip(moon, epochs))
    epairs = collect(zip(elts, epochs))
    # A ray from the Moon to the centre of the Earth, and the sub-Moon point.
    rays = [(et, -Vector(s[SOneTo(3)])) for (s, et) in pairs]
    spoints = [(et, subpnt(ctx, "NEAR POINT/ELLIPSOID", "EARTH", et, "IAU_EARTH", "LT+S", "MOON")[1]) for et in epochs]
    pict = "YYYY-MON-DD HR:MN:SC.### ::TDB"
    # The same inputs as `Vector`s for SPICE.jl, so that no conversion is timed.
    cpairs = [(Vector(s), et) for (s, et) in pairs]
    cepairs = [(Vector(e), et) for (e, et) in epairs]
    cspoints = [(et, Vector(p)) for (et, p) in spoints]

    agree = String[
        agreement("str2et", s -> str2et(ctx, s), SPICE.str2et, strings, 1.0e-6),
        agreement("unitim", et -> unitim(ctx, et, "TDB", "JDTDB"), et -> SPICE.unitim(et, "TDB", "JDTDB"), epochs, 1.0e-9),
        agreement("oscelt", ((s, et),) -> oscelt(s, et, GM_EARTH), ((s, et),) -> SPICE.oscelt(Vector(s), et, GM_EARTH), pairs, 1.0e-6),
        agreement("conics", ((e, et),) -> conics(e, et), ((e, et),) -> SPICE.conics(Vector(e), et), epairs, 1.0e-6),
        agreement(
            "subpnt NEAR POINT/ELLIPSOID",
            et -> subpnt(ctx, "NEAR POINT/ELLIPSOID", "EARTH", et, "IAU_EARTH", "LT+S", "MOON")[1],
            et -> SPICE.subpnt("NEAR POINT/ELLIPSOID", "EARTH", et, "IAU_EARTH", "LT+S", "MOON")[1], epochs, 1.0e-6
        ),
        agreement(
            "sincpt ELLIPSOID",
            ((et, d),) -> sincpt(ctx, "ELLIPSOID", "EARTH", et, "IAU_EARTH", "CN+S", "MOON", "J2000", d)[1],
            ((et, d),) -> SPICE.sincpt("ELLIPSOID", "EARTH", et, "IAU_EARTH", "CN+S", "MOON", "J2000", d)[1], rays, 1.0e-6
        ),
        agreement(
            "ilumin ELLIPSOID",
            ((et, p),) -> ilumin(ctx, "ELLIPSOID", "EARTH", et, "IAU_EARTH", "LT+S", "MOON", p)[3],
            ((et, p),) -> SPICE.ilumin("ELLIPSOID", "EARTH", et, "IAU_EARTH", "MOON", Vector(p); abcorr = "LT+S")[3], spoints, 1.0e-12
        ),
    ]
    @assert all(et -> et2utc(ctx, et, "C", 3) == SPICE.et2utc(et, "C", 3), epochs)
    @assert all(et -> timout(ctx, et, pict) == SPICE.timout(et, pict), epochs)
    @assert all(
        et -> occult(ctx, "MOON", "ELLIPSOID", "IAU_MOON", "SUN", "ELLIPSOID", "IAU_SUN", "CN", "EARTH", et) ==
            SPICE.occult("MOON", "ELLIPSOID", "IAU_MOON", "SUN", "ELLIPSOID", "IAU_SUN", "CN", "EARTH", et), epochs
    )

    # The search: the maxima of the Earth-Moon distance over one year, step 5 days.
    t0 = epochs[60]; t1 = t0 + 365.25 * DAY; step = 5 * DAY
    gf_almanac(_) = sum(window_values(gfdist(ctx, "MOON", "NONE", "EARTH", "LOCMAX", 0.0, 0.0, step, Window([t0, t1]))))
    function gf_cspice(_)
        cnfine = SPICE.SpiceDoubleCell(2)
        SPICE.wninsd!(cnfine, t0, t1)
        return sum(cspice_window_values(SPICE.gfdist("MOON", "NONE", "EARTH", "LOCMAX", 0.0, 0.0, step, 1000, cnfine)))
    end
    wa = window_values(gfdist(ctx, "MOON", "NONE", "EARTH", "LOCMAX", 0.0, 0.0, step, Window([t0, t1])))
    cn = SPICE.SpiceDoubleCell(2); SPICE.wninsd!(cn, t0, t1)
    wc = cspice_window_values(SPICE.gfdist("MOON", "NONE", "EARTH", "LOCMAX", 0.0, 0.0, step, 1000, cn))
    @assert length(wa) == length(wc) >= 24 "gfdist: $(length(wa)) and $(length(wc)) endpoints"
    push!(agree, @sprintf("| gfdist LOCMAX, one year (s) | %.2e |", maximum(abs.(wa .- wc))))
    @assert maximum(abs.(wa .- wc)) <= 1.0e-3

    rows = Row[]
    add(g, name, f, xs; n = length(xs)) = push!(rows, bench_row(g, name, loop, f, xs; n))

    add("str2et ISO string", "Almanac str2et", s -> str2et(ctx, s), strings)
    add("str2et ISO string", "CSPICE str2et (SPICE.jl)", SPICE.str2et, strings)
    add("et2utc C 3", "Almanac et2utc", et -> Float64(length(et2utc(ctx, et, "C", 3))), epochs)
    add("et2utc C 3", "CSPICE et2utc (SPICE.jl)", et -> Float64(length(SPICE.et2utc(et, "C", 3))), epochs)
    add("timout TDB picture", "Almanac timout", et -> Float64(length(timout(ctx, et, pict))), epochs)
    add("timout TDB picture", "CSPICE timout (SPICE.jl)", et -> Float64(length(SPICE.timout(et, pict))), epochs)
    add("unitim TDB to JDTDB", "Almanac unitim", et -> unitim(ctx, et, "TDB", "JDTDB"), epochs)
    add("unitim TDB to JDTDB", "CSPICE unitim (SPICE.jl)", et -> SPICE.unitim(et, "TDB", "JDTDB"), epochs)

    add("bodn2c MOON", "Almanac bodn2c", _ -> Float64(bodn2c(ctx, "MOON")), epochs)
    add("bodn2c MOON", "CSPICE bodn2c (SPICE.jl)", _ -> Float64(SPICE.bodn2c("MOON")), epochs)
    add("namfrm IAU_EARTH", "Almanac namfrm", _ -> Float64(namfrm(ctx, "IAU_EARTH")), epochs)
    add("namfrm IAU_EARTH", "CSPICE namfrm (SPICE.jl)", _ -> Float64(SPICE.namfrm("IAU_EARTH")), epochs)
    add("bodvrd EARTH RADII", "Almanac bodvrd", _ -> bodvrd(ctx, "EARTH", "RADII")[1], epochs)
    add("bodvrd EARTH RADII", "CSPICE bodvrd (SPICE.jl)", _ -> SPICE.bodvrd("EARTH", "RADII")[1], epochs)

    add("oscelt Moon about Earth", "Almanac oscelt", ((s, et),) -> oscelt(s, et, GM_EARTH)[2], pairs)
    add("oscelt Moon about Earth", "CSPICE oscelt (SPICE.jl)", ((s, et),) -> SPICE.oscelt(s, et, GM_EARTH)[2], cpairs)
    add("conics Moon about Earth", "Almanac conics", ((e, et),) -> conics(e, et)[1], epairs)
    add("conics Moon about Earth", "CSPICE conics (SPICE.jl)", ((e, et),) -> SPICE.conics(e, et)[1], cepairs)

    for method in ("NEAR POINT/ELLIPSOID", "INTERCEPT/ELLIPSOID")
        g = "subpnt $method, LT+S"
        add(g, "Almanac subpnt", et -> subpnt(ctx, method, "EARTH", et, "IAU_EARTH", "LT+S", "MOON")[1][1], epochs)
        add(g, "CSPICE subpnt (SPICE.jl)", et -> SPICE.subpnt(method, "EARTH", et, "IAU_EARTH", "LT+S", "MOON")[1][1], epochs)
    end
    g = "sincpt ELLIPSOID, CN+S"
    add(g, "Almanac sincpt", ((et, d),) -> sincpt(ctx, "ELLIPSOID", "EARTH", et, "IAU_EARTH", "CN+S", "MOON", "J2000", d)[1][1], rays)
    add(g, "CSPICE sincpt (SPICE.jl)", ((et, d),) -> SPICE.sincpt("ELLIPSOID", "EARTH", et, "IAU_EARTH", "CN+S", "MOON", "J2000", d)[1][1], rays)
    g = "ilumin ELLIPSOID, LT+S"
    add(g, "Almanac ilumin", ((et, p),) -> ilumin(ctx, "ELLIPSOID", "EARTH", et, "IAU_EARTH", "LT+S", "MOON", p)[3], spoints)
    add(g, "CSPICE ilumin (SPICE.jl)", ((et, p),) -> SPICE.ilumin("ELLIPSOID", "EARTH", et, "IAU_EARTH", "MOON", p; abcorr = "LT+S")[3], cspoints)
    g = "occult Moon and Sun from Earth, CN"
    add(g, "Almanac occult", et -> Float64(occult(ctx, "MOON", "ELLIPSOID", "IAU_MOON", "SUN", "ELLIPSOID", "IAU_SUN", "CN", "EARTH", et)), epochs)
    add(g, "CSPICE occult (SPICE.jl)", et -> Float64(SPICE.occult("MOON", "ELLIPSOID", "IAU_MOON", "SUN", "ELLIPSOID", "IAU_SUN", "CN", "EARTH", et)), epochs)

    g = "gfdist LOCMAX, one year, step 5 d"
    add(g, "Almanac gfdist", gf_almanac, 1:1)
    add(g, "CSPICE gfdist (SPICE.jl)", gf_cspice, 1:1)

    # Last, since the CSPICE row clears the pool of CSPICE.
    g = "load the standard set (6 files)"
    add(g, "Almanac furnsh, new context", _ -> (c = SpiceContext(); foreach(p -> furnsh(c, p), set); 0.0), 1:1)
    add(g, "CSPICE kclear + furnsh (SPICE.jl)", _ -> (SPICE.kclear(); foreach(SPICE.furnsh, set); 0.0), 1:1)

    text = run_header("Benchmark of further calls: Almanac, CSPICE") * """
        Kernels: the standard set of the tests (naif0012.tls, pck00011.tpc,
        moon_de440_200625.tf, the MOON_PA alias kernel, moon_pa_de440_200625.bpc,
        de440s.bsp), in both engines. Epochs: the 100 of the state benchmark; the
        search and the load are one query per call.

        Largest difference from CSPICE on the timed inputs (km, radians, days or
        seconds); `et2utc`, `timout` and `occult` agree exactly:

        | Call | max abs difference |
        |---|---:|
        $(join(agree, "\n"))

        ## Time per query

        """ * markdown_table(rows)
    println()
    println(text)
    write_results("api", text, rows)
    SPICE.kclear()
    return rows
end

main()
