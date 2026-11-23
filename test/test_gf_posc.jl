# Geometry finder, coordinates of a position vector: `gfposc` against `gfposc_c`, for
# every coordinate system and coordinate, every relation, several targets, frames,
# aberration corrections, steps and adjust values, and confinement windows of one and of
# several intervals. Longitude and right ascension are searched across the branch cut
# (at pi, and at 0 and 2 pi), with reference values on the cut and next to it. The
# number of intervals must be equal and each end point must agree within 2e-6 s (a root
# where the coordinate is almost stationary is accepted with a larger difference if the
# coordinate has the same value at the two epochs, to rounding); an inequality window is
# also checked against the value of the coordinate at random epochs, which does not
# depend on CSPICE.
#
# The default run compares one case in ten of the product of options and one case in
# thirteen at the branch cut (`gf_pick`), and does not make the table of search times;
# the exhaustive run (`ALMANAC_TEST_FULL=1`) compares every case and makes the table.

isdefined(@__MODULE__, :GFO) || include("gforacle.jl")

"""
    gf_inequality_holds(q, w, cn, rel, refval; n=150, margin=5.0) -> Bool

True if, at `n` random epochs of the window `cn` that are more than `margin` seconds
from every end point of `w` and of `cn`, the epoch is in `w` exactly when the value of
the quantity `q` is less than (`rel == "<"`) or greater than (`">"`) `refval`.
"""
function gf_inequality_holds(
        rng::AbstractRNG, q, w::Window, cn::Window, rel::AbstractString,
        refval::Real; n::Integer = 150, margin::Float64 = 5.0
    )
    ok = true
    ends = vcat(collect(w), collect(cn))
    for _ in 1:n
        k = rand(rng, 1:wncard(cn))
        t = cn[2k - 1] + (cn[2k] - cn[2k - 1]) * rand(rng)
        any(e -> abs(e - t) <= margin, ends) && continue
        v = Almanac._gf_value(q, t)
        holds = rel == "<" ? v < refval : v > refval
        ok &= (holds == wnelmd(w, t))
    end
    return ok
end

"Allocation of the second call of `Almanac._gf_coords(q, p)`; the first call compiles."
function gf_coords_alloc(q, p)
    Almanac._gf_coords(q, p)
    return @allocated Almanac._gf_coords(q, p)
end

@gf_light group("gf coordinates: gfposc, longitude and right ascension", Kernels.FULL_SET_FILES...) do
    oracle_load(Kernels.full_set())
    ctx = neo_context(Kernels.full_set())
    ks = Almanac.kernelset(ctx)
    rng = MersenneTwister(SEED + 6300)
    day = 86400.0

    # (target, observer, frame, corrections, steps, span, the frame center has radii)
    cases = (
        ("MOON", "EARTH", "J2000", GF_ABCORRS, (day, 2 * day), 90 * day, false),
        ("MOON", "EARTH", "IAU_EARTH", ("NONE", "LT+S", "CN", "XLT", "XCN+S"), (3 * 3600.0,), 4 * day, true),
        ("SUN", "EARTH", "J2000", ("NONE", "LT+S", "XCN"), (10 * day,), 800 * day, false),
        ("EARTH", "MOON", "MOON_ME", ("NONE", "LT", "CN+S"), (day,), 100 * day, true),
        ("SUN", "EARTH", "IAU_MARS", ("NONE",), (2 * 3600.0,), 3 * day, true),
    )

    function start_epoch()
        return -3.0e8 + 1.4e9 * rand(rng)
    end
    function confinement(span)
        t0 = start_epoch()
        return rand(rng) < 0.7 ? Window([t0, t0 + span]) :
            gf_random_window(rng, 3, t0, t0 + span)
    end

    @testset "the planetographic longitude of Mars is positive west" begin
        q = Almanac._gf_position_coordinate(
            ks, "SUN", "IAU_MARS", "NONE", "EARTH",
            "PLANETOGRAPHIC", "LONGITUDE"
        )
        @test q.sense == -1
        q = Almanac._gf_position_coordinate(
            ks, "MOON", "IAU_EARTH", "NONE", "EARTH",
            "PLANETOGRAPHIC", "LONGITUDE"
        )
        @test q.sense == 1
    end

    @testset "every system, coordinate and relation against CSPICE" begin
        nerr = 0
        i = 0
        fails = Failures()
        for (targ, obs, frame, abcorrs, steps, span, hasradii) in cases,
                (sys, crds) in GF_SYSTEMS, crd in crds, (rel, wadj) in GF_VARIANTS

            i += 1
            # (in every run: the error cases, and each relation on the longitude that is
            # positive west, which has a branch of its own for each relation)
            gf_pick(i, 10) || (!hasradii && sys in ("GEODETIC", "PLANETOGRAPHIC")) ||
                (frame == "IAU_MARS" && sys == "PLANETOGRAPHIC" && crd == "LONGITUDE") || continue
            abcorr = rand(rng, abcorrs)
            step = rand(rng, steps)
            cn = confinement(span)
            if !hasradii && sys in ("GEODETIC", "PLANETOGRAPHIC")
                # the center of the frame is not a body with radii
                rel == "<" || continue
                @test gf_kind(
                    () -> gfposc(
                        ctx, targ, frame, abcorr, obs, sys, crd, rel, 0.0, 0.0,
                        step, cn
                    )
                ) === :badframe
                @test GFO.kind(
                    () -> GFO.gfposc(
                        targ, frame, abcorr, obs, sys, crd, rel, 0.0, 0.0,
                        step, collect(cn)
                    )
                ) === :badframe
                nerr += 1
                continue
            end
            q = Almanac._gf_position_coordinate(ks, targ, frame, abcorr, obs, sys, crd)
            refval, adjust = gf_reference(rng, q, cn, rel, wadj)
            w = gfposc(ctx, targ, frame, abcorr, obs, sys, crd, rel, refval, adjust, step, cn)
            ref = GFO.gfposc(targ, frame, abcorr, obs, sys, crd, rel, refval, adjust, step, collect(cn))
            angular = (sys, crd) in GF_ANGULAR
            label = (angular ? "gfposc longitude/RA " : "gfposc other ") * gf_family(rel, adjust)
            scale = maximum(abs(Almanac._gf_value(q, Float64(t))) for t in range(cn[1], cn[end]; length = 20))
            ok = gf_compare(label, w, ref; value = t -> Almanac._gf_value(q, t), vtol = 64 * eps(scale))
            ok || @info "gfposc case" targ obs frame abcorr sys crd rel refval adjust step cn
            @test ok
            check!(fails, all(t -> wnelmd(cn, t), w), "gfposc inside the confinement window", targ, obs, frame, abcorr, sys, crd, rel, refval, adjust)
            if rel in ("<", ">")
                check!(fails, gf_inequality_holds(rng, q, w, cn, rel, refval), "gfposc inequality holds", targ, obs, frame, abcorr, sys, crd, rel, refval)
            end
        end
        @test isempty(fails)
        @test nerr == 2 * 6
    end

    @testset "longitude and right ascension at and next to the branch cut" begin
        refs = (
            0.0, Float64(pi), -Float64(pi), pi / 2, -pi / 2, 3pi / 2, 2pi, pi - 1.0e-3,
            -pi + 1.0e-3, 1.0e-3, 2pi - 1.0e-3, 4.0, -2.0, 7.0,
        )
        i = 0
        fails = Failures()
        for (targ, obs, frame, abcorrs, steps, span, hasradii) in cases, (sys, crd) in GF_ANGULAR
            (!hasradii && sys in ("GEODETIC", "PLANETOGRAPHIC")) && continue
            abcorr = first(abcorrs)
            step = first(steps)
            t0 = start_epoch()
            cn = Window([t0, t0 + span])
            q = Almanac._gf_position_coordinate(ks, targ, frame, abcorr, obs, sys, crd)
            for refval in refs, rel in ("<", "=", ">")
                gf_pick(i += 1, 13) || continue
                w = gfposc(ctx, targ, frame, abcorr, obs, sys, crd, rel, refval, 0.0, step, cn)
                ref = GFO.gfposc(targ, frame, abcorr, obs, sys, crd, rel, refval, 0.0, step, collect(cn))
                ok = gf_compare("gfposc branch cut " * gf_family(rel, 0.0), w, ref)
                ok || @info "gfposc branch-cut case" targ obs frame abcorr sys crd rel refval
                @test ok
                # A reference value outside the range of the coordinate is taken modulo
                # 2 pi by the proxy functions, so the direct check is for values in range.
                lo, hi = sys in ("LATITUDINAL", "SPHERICAL", "GEODETIC") ? (-pi, pi) : (0.0, 2pi)
                if rel != "=" && lo < refval < hi
                    check!(fails, gf_inequality_holds(rng, q, w, cn, rel, refval), "gfposc branch cut: inequality holds", targ, obs, frame, abcorr, sys, crd, rel, refval)
                end
            end
            for (rel, adjust) in (
                    ("LOCMIN", 0.0), ("LOCMAX", 0.0), ("ABSMIN", 0.0), ("ABSMAX", 0.0),
                    ("ABSMIN", 0.3), ("ABSMAX", 0.3), ("ABSMIN", 4.0), ("ABSMAX", 4.0),
                )
                gf_pick(i += 1, 13) || continue
                w = gfposc(ctx, targ, frame, abcorr, obs, sys, crd, rel, 0.0, adjust, step, cn)
                ref = GFO.gfposc(targ, frame, abcorr, obs, sys, crd, rel, 0.0, adjust, step, collect(cn))
                ok = gf_compare("gfposc branch cut " * gf_family(rel, adjust), w, ref)
                ok || @info "gfposc branch-cut case" targ obs frame abcorr sys crd rel adjust
                @test ok
            end
        end
        @test isempty(fails)
    end

    @testset "a right ascension that crosses 0 gives no false event" begin
        # The right ascension of the Moon goes through 2 pi -> 0 once in a month. The jump
        # is no root of "= pi" and no local extremum.
        t0 = SPICE.str2et("2024 MAR 01")
        cn = Window([t0, t0 + 120 * day])
        q = Almanac._gf_position_coordinate(
            ks, "MOON", "J2000", "NONE", "EARTH", "RA/DEC",
            "RIGHT ASCENSION"
        )
        ra(t) = Almanac._gf_value(q, t)
        roots = gfposc(
            ctx, "MOON", "J2000", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION", "=",
            pi, 0.0, day, cn
        )
        @test wncard(roots) in (4, 5)
        @test all(t -> abs(ra(t) - pi) < 1.0e-9, roots[1:2:end])
        for rel in ("LOCMIN", "LOCMAX")
            @test isempty(
                gfposc(
                    ctx, "MOON", "J2000", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION",
                    rel, 0.0, 0.0, day, cn
                )
            )
        end
        # A root next to the cut is found on both sides of it.
        for refval in (1.0e-4, 2pi - 1.0e-4)
            w = gfposc(
                ctx, "MOON", "J2000", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION", "=",
                refval, 0.0, day, cn
            )
            @test wncard(w) in (4, 5)
            @test all(t -> abs(ra(t) - refval) < 1.0e-9, w[1:2:end])
        end
        # The same events as a longitude: RA = 3 is longitude = 3, RA = 4 is longitude 4 - 2 pi.
        for refval in (3.0, 4.0)
            a = gfposc(
                ctx, "MOON", "J2000", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION", "=",
                refval, 0.0, day, cn
            )
            b = gfposc(
                ctx, "MOON", "J2000", "NONE", "EARTH", "LATITUDINAL", "LONGITUDE", "=",
                refval > pi ? refval - 2pi : refval, 0.0, day, cn
            )
            @test length(a) == length(b) && maximum(abs.(collect(a) .- collect(b))) < 1.0e-6
        end
        # '<' and '>' are complements in the confinement window, but for the end points.
        lt = gfposc(
            ctx, "MOON", "J2000", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION", "<", 2.0,
            0.0, day, cn
        )
        gt = gfposc(
            ctx, "MOON", "J2000", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION", ">", 2.0,
            0.0, day, cn
        )
        @test abs(wnsumd(lt)[1] + wnsumd(gt)[1] - 120 * day) < 1.0e-3
        @test wnsumd(wnintd(lt, gt))[1] < 1.0e-3
    end

    @testset "ABSMIN and ABSMAX of longitude with adjust and no extremum" begin
        # Regression for a defect of CSPICE N0067 (`zzgflong`): when no region holds an
        # extremum, CSPICE reads the first interval of the empty window to apply the
        # adjust value, so it is not called. At this epoch (a root of y found by bisection)
        # the Moon seen from the Earth is on the negative x axis of J2000 to 1e-14 rad, so
        # the epoch is in none of the regions where the extrema are searched (a sine
        # above sin(1e-11) or below its negative, or a positive cosine).
        t = 1.0125953371533352e8
        p = spkpos(ctx, "MOON", t, "J2000", "NONE", "EARTH")[1]
        @test p[1] < 0 && abs(p[2] / p[1]) < 1.0e-12
        for rel in ("ABSMIN", "ABSMAX")
            w = gfposc(ctx, "MOON", "J2000", "NONE", "EARTH", "LATITUDINAL", "LONGITUDE", rel, 0.0, 0.1, 3600.0, Window([t, t]))
            @test isempty(w)
        end
    end

    @testset "longitude and right ascension of reclat_, as zzgfcoq" begin
        # Regression: `zzgfcoq` and the cosine and sine proxies of `zzgflong` take the
        # longitude and the right ascension from the Fortran `reclat_` and `recrad_`,
        # the arc tangent of the components as they are. Almanac took them from
        # `reclat_c` and `recrad_c`, which divide the components by the largest one
        # first, and so gave another last bit for about one vector in sixteen. Julia's
        # `atan` can differ in the last bit from the `atan2` of the C library, so a
        # random vector is compared with the formula of `reclat_` bit for bit and with
        # the library within 4 ulps; the vectors with denormal x and y, whose longitude
        # is 0 in `reclat_c`, are compared with the library bit for bit. Own random
        # stream.
        r2 = MersenneTwister(SEED + 6301)
        qlon = Almanac._gf_position_coordinate(
            ks, "MOON", "J2000", "LT+S", "EARTH", "LATITUDINAL", "LONGITUDE"
        )
        qra = Almanac._gf_position_coordinate(
            ks, "MOON", "J2000", "LT+S", "EARTH", "RA/DEC", "RIGHT ASCENSION"
        )
        fails = Failures()
        for _ in 1:500
            p = SVector{3}(randn(r2, 3) .* 10.0 .^ rand(r2, 0:8))
            lon = atan(p[2], p[1])
            ra = lon < 0 ? lon + 2pi : lon
            a = Almanac._gf_coords(qlon, p)
            b = Almanac._gf_coords(qra, p)
            check!(fails, a[2] == lon, "longitude, the formula of reclat_", p)
            check!(fails, b[2] == ra, "right ascension, the formula of recrad_", p)
            check!(fails, isapprox(a[2], GFO.reclat_f(p)[2]; rtol = 4 * eps()), "longitude, reclat_", p)
            check!(fails, isapprox(b[2], GFO.recrad_f(p)[2]; rtol = 4 * eps()), "right ascension, recrad_", p)
        end
        for p in (SVector(5.0e-324, 5.0e-324, 1.0e10), SVector(0.0, -5.0e-324, 3.0e10))
            check!(fails, Almanac._gf_coords(qlon, p) == GFO.reclat_f(p), "denormal, reclat_", p)
            check!(fails, Almanac._gf_coords(qra, p) == GFO.recrad_f(p), "denormal, recrad_", p)
        end
        # the value of the quantity, from the position of the search
        for _ in 1:40
            t = -3.0e8 + 1.4e9 * rand(r2)
            p = spkpos(ctx, "MOON", t, "J2000", "LT+S", "EARTH")[1]
            lon = atan(p[2], p[1])
            check!(fails, Almanac._gf_value(qlon, t) == lon, "value, longitude", t)
            check!(fails, Almanac._gf_value(qra, t) == (lon < 0 ? lon + 2pi : lon), "value, right ascension", t)
        end
        @test isempty(fails)
    end

    @testset "planetographic coordinates without a lookup in the kernel pool" begin
        # Regression: each evaluation of a planetographic coordinate called `recpgr`,
        # which reads the sense of the longitude from the kernel pool (two lookups, each
        # with a new key string) although the sense is fixed for the search. The
        # coordinates are those of `recpgr`, bit for bit, and the conversion does not
        # allocate. Own random stream.
        r3 = MersenneTwister(SEED + 6302)
        for (targ, frame) in (("SUN", "IAU_MARS"), ("MOON", "IAU_EARTH"))
            q = Almanac._gf_position_coordinate(
                ks, targ, frame, "NONE", "EARTH", "PLANETOGRAPHIC", "LONGITUDE"
            )
            fails = Failures()
            for _ in 1:200
                p = SVector{3}(randn(r3, 3) .* 10.0 .^ rand(r3, 3:8))
                check!(fails, Almanac._gf_coords(q, p) == recpgr(ks, q.center, p, q.re, q.f), frame, p)
            end
            @test isempty(fails)
            @test gf_coords_alloc(q, SVector(1.0e4, -2.0e4, 3.0e3)) == 0
        end
    end

    @testset "errors, with the kinds of CSPICE" begin
        cn = Window([0.0, 30 * day])
        c = collect(cn)
        neo(a...) = gf_kind(() -> gfposc(ctx, a..., cn))
        ora(a...) = GFO.kind(() -> GFO.gfposc(a..., c))
        for (kind, a) in (
                (:notsupported, ("MOON", "J2000", "NONE", "EARTH", "POLAR", "RADIUS", "=", 1.0, 0.0, day)),
                (:notsupported, ("MOON", "J2000", "NONE", "EARTH", "RA/DEC", "LONGITUDE", "=", 1.0, 0.0, day)),
                (:notsupported, ("MOON", "J2000", "NONE", "EARTH", "RECTANGULAR", "W", "=", 1.0, 0.0, day)),
                (:idcodenotfound, ("NOSUCHBODY", "J2000", "NONE", "EARTH", "RECTANGULAR", "X", "=", 1.0, 0.0, day)),
                (:bodiesnotdistinct, ("EARTH", "J2000", "NONE", "399", "RECTANGULAR", "X", "=", 1.0, 0.0, day)),
                (:badframe, ("MOON", "J2000", "NONE", "EARTH", "GEODETIC", "ALTITUDE", "=", 1.0, 0.0, day)),
                (:noframe, ("MOON", "NOSUCHFRAME", "NONE", "EARTH", "GEODETIC", "ALTITUDE", "=", 1.0, 0.0, day)),
                (:unknownframe, ("MOON", "NOSUCHFRAME", "NONE", "EARTH", "RECTANGULAR", "X", "=", 1.0, 0.0, day)),
                (:invalidstep, ("MOON", "J2000", "NONE", "EARTH", "RECTANGULAR", "X", "=", 1.0, 0.0, 0.0)),
                (:notrecognized, ("MOON", "J2000", "NONE", "EARTH", "RECTANGULAR", "X", "!=", 1.0, 0.0, day)),
                (:notrecognized, ("MOON", "J2000", "NONE", "EARTH", "LATITUDINAL", "LONGITUDE", "!=", 1.0, 0.0, day)),
                (:valueoutofrange, ("MOON", "J2000", "NONE", "EARTH", "RECTANGULAR", "X", "ABSMAX", 1.0, -1.0, day)),
                (:valueoutofrange, ("MOON", "J2000", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION", "ABSMAX", 1.0, -1.0, day)),
                # `zzgfcslv` checks `adjust` before it sets up the quantity
                (:valueoutofrange, ("NOSUCHBODY", "J2000", "NONE", "EARTH", "RECTANGULAR", "X", "ABSMAX", 1.0, -1.0, day)),
                (:valueoutofrange, ("MOON", "J2000", "NONE", "EARTH", "POLAR", "RADIUS", "ABSMAX", 1.0, -1.0, day)),
                (:invalidvalue, ("MOON", "J2000", "NONE", "EARTH", "RECTANGULAR", "X", "<", 1.0, 1.0, day)),
                (:invalidvalue, ("MOON", "J2000", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION", "<", 1.0, 1.0, day)),
            )
            @test neo(a...) === kind
            @test ora(a...) === kind
        end
        # Regression: a NaN reference value is `:invalidvalue` for a longitude and a right
        # ascension too (`_gf_long`, `zzgflong`), as for the other coordinates
        # (`_gf_relx`), before the test for an empty window. CSPICE takes the NaN.
        for (sys, crd) in (("LATITUDINAL", "LONGITUDE"), ("RA/DEC", "RIGHT ASCENSION"), ("RECTANGULAR", "X")),
                rel in ("<", "="), w in (cn, Window())
            @test gf_kind(
                () -> gfposc(ctx, "MOON", "J2000", "NONE", "EARTH", sys, crd, rel, NaN, 0.0, day, w)
            ) === :invalidvalue
        end
        # names are read without regard to case and blanks
        a = gfposc(ctx, "MOON", "J2000", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION", "=", 1.0, 0.0, day, cn)
        @test gfposc(
            ctx, "moon", " j2000", "none", "earth", " ra / dec ", " right   ascension ", "=", 1.0,
            0.0, day, cn
        ) == a
        @test gfposc(
            ks, 301, 1, Aberration("NONE"), 399, "RA/DEC", "RIGHT ASCENSION", "=", 1.0, 0.0,
            day, cn
        ) == a
        @test isempty(
            gfposc(
                ctx, "MOON", "J2000", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION", "=",
                1.0, 0.0, day, Window()
            )
        )
    end

    FULL && @testset "search time against CSPICE (reported, not asserted)" begin
        t0 = SPICE.str2et("2020 JAN 01")
        cn = Window([t0, t0 + 3650 * day])
        c = collect(cn)
        cases2 = (
            (
                "gfposc Moon RA/DEC declination, 10 years, LOCMAX, LT+S",
                ("MOON", "J2000", "LT+S", "EARTH", "RA/DEC", "DECLINATION", "LOCMAX", 0.0, 0.0, day),
            ),
            (
                "gfposc Moon right ascension, 10 years, '=', LT+S",
                ("MOON", "J2000", "LT+S", "EARTH", "RA/DEC", "RIGHT ASCENSION", "=", 1.0, 0.0, day),
            ),
            (
                "gfposc Moon latitudinal longitude, 10 years, '<', NONE",
                ("MOON", "J2000", "NONE", "EARTH", "LATITUDINAL", "LONGITUDE", "<", 3.0, 0.0, day),
            ),
        )
        for (name, a) in cases2
            neo = () -> gfposc(ctx, a..., cn)
            ora = () -> GFO.gfposc(a..., c)
            @test gf_compare("timing", neo(), ora())
            tn, to = gf_time(neo), gf_time(ora)
            measured(
                @sprintf(
                    "  %-58s Almanac %8.2f ms  CSPICE %8.2f ms  ratio %.2f", name,
                    1.0e3 * tn, 1.0e3 * to, to / tn
                )
            )
        end
    end

    gf_report("gfposc")
end
