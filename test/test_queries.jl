# Group 53: prepared queries (`state_query`, `rotation_query`, src/queries/queries.jl).
#
# Invariant layer: a query gives what `spkezr`, `spkpos`, `pxform` and `sxform` give for
# the same arguments, to the last bit, for every aberration correction, an inertial and
# a rotating frame, names and codes. Then: the query keeps its kernels, dual-number
# epochs, threads, errors. Bodies and frames come from functions (`add_ephemeris!`) and
# from pool constants, so no kernel file of NAIF is needed.
#
# The chain plans (`src/queries/plans.jl`): chains of SPK segments and of frames whose
# links the query finds once give the same results and the same errors as the general
# path, for the shapes of chain that the plans distinguish (the observer or the target
# at the common body, chains on both sides, frames that differ from link to link, a
# list of segments with two centers, gaps in the coverage of a link, exact segment
# bounds; frame chains that reach the other frame, meet it at its first frame or at a
# later one, with and without a PCK frame). A last part repeats the comparison on
# de440s and the binary PCK of the Moon, also outside their coverage.

using ForwardDiff

@testset "prepared queries" begin
    earth(et) = [1.5e8 + 1.0e-3 * et, -2.0e7 + 29.0 * et, 3.0e6 - 0.5 * et, 1.0e-3, 29.0, -0.5]
    probe(et) = [
        7000.0 + 0.002 * et, -1200.0 + 7.4 * sin(1.0e-3 * et), 300.0, 0.002,
        7.4e-3 * cos(1.0e-3 * et), 0.0,
    ]
    t0, t1 = -86400.0, 86400.0
    ctx = SpiceContext()
    pdpool(ctx, "BODY399_POLE_RA", [10.0, 0.1, 0.0])
    pdpool(ctx, "BODY399_POLE_DEC", [80.0, -0.05, 0.0])
    pdpool(ctx, "BODY399_PM", [190.0, 360.9856, 0.0])
    add_ephemeris!(earth, ctx, 399, 0, "J2000", t0, t1)
    add_ephemeris!(probe, ctx, -1000, "EARTH", "J2000", t0, t1)
    boddef(ctx, "PROBE", -1000)
    ets = range(t0 + 1000, t1 - 1000; length = 5)

    @testset "the results of spkezr and spkpos, bit for bit" begin
        fails = Failures()
        for corr in ("NONE", "LT", "LT+S", "CN", "CN+S", "XLT", "XLT+S", "XCN", "XCN+S"),
                ref in ("J2000", "ECLIPJ2000", "IAU_EARTH"),
                (targ, obs) in (("PROBE", "EARTH"), (399, -1000), ("EARTH", "SSB"))
            q = state_query(ctx, targ, ref, corr, obs)
            for et in ets
                check!(fails, q(et) === spkezr(ctx, targ, et, ref, corr, obs), "state_query", corr, ref, targ, obs, et)
                check!(fails, spkpos(q, et) === spkpos(ctx, targ, et, ref, corr, obs), "spkpos of query", corr, ref, targ, obs, et)
                check!(fails, spkder(q, et, Val(1)) === spkder(ctx, targ, et, ref, corr, obs, Val(1)), "spkder of query", corr, ref, targ, obs, et)
            end
        end
        @test isempty(fails)
        q = state_query(ctx, "PROBE", "J2000", Aberration("LT+S"), "EARTH")
        @test spkezr.(q, ets) == [spkezr(ctx, "PROBE", et, "J2000", "LT+S", "EARTH") for et in ets]
    end

    @testset "the results of pxform and sxform, bit for bit" begin
        fails = Failures()
        for (from, to) in (("J2000", "IAU_EARTH"), ("IAU_EARTH", "ECLIPJ2000"), (1, 17))
            q = rotation_query(ctx, from, to)
            for et in ets
                check!(fails, q(et) === pxform(ctx, from, to, et), "rotation_query", from, to, et)
                check!(fails, sxform(q, et) === sxform(ctx, from, to, et), "sxform of query", from, to, et)
            end
        end
        @test isempty(fails)
    end

    "The value of `f()`, or the type, kind and message of the `SpiceError` it throws."
    outcome(f) = try
        f()
    catch e
        e isa SpiceError ? (typeof(e), e.kind, e.msg) : rethrow()
    end

    @testset "chain plans: the results and errors of the general path" begin
        moonlet(et) = [3.8e5, 1.0e3 * sin(1.0e-5 * et), -0.0, 0.0, 1.0e-2 * cos(1.0e-5 * et), 0.0]
        add_ephemeris!(moonlet, ctx, -2000, 399, "ECLIPJ2000", t0, t1)
        add_ephemeris!(et -> earth(et) .* 2, ctx, -3000, 0, "ECLIPJ2000", t0, t1)
        # Two centers: no plan.
        add_ephemeris!(probe, ctx, -4000, 399, "J2000", t0, 0.0)
        add_ephemeris!(et -> earth(et) .+ 5, ctx, -4000, 0, "J2000", 0.0, t1)
        # A gap in the coverage of one list.
        add_ephemeris!(probe, ctx, -5000, 399, "J2000", t0, -1000.0)
        add_ephemeris!(et -> probe(et) .- 3, ctx, -5000, 399, "ECLIPJ2000", 1000.0, t1)
        lmpool(
            ctx, [
                "FRAME_TKA = -100001", "FRAME_-100001_NAME = 'TKA'",
                "FRAME_-100001_CLASS = 4", "FRAME_-100001_CLASS_ID = -100001",
                "FRAME_-100001_CENTER = 399", "TKFRAME_-100001_RELATIVE = 'IAU_EARTH'",
                "TKFRAME_-100001_SPEC = 'ANGLES'", "TKFRAME_-100001_UNITS = 'DEGREES'",
                "TKFRAME_-100001_AXES = ( 3, 1, 3 )", "TKFRAME_-100001_ANGLES = ( 10.0, 20.0, 30.0 )",
                "FRAME_TKB = -100002", "FRAME_-100002_NAME = 'TKB'",
                "FRAME_-100002_CLASS = 4", "FRAME_-100002_CLASS_ID = -100002",
                "FRAME_-100002_CENTER = 399", "TKFRAME_-100002_RELATIVE = 'TKA'",
                "TKFRAME_-100002_SPEC = 'ANGLES'", "TKFRAME_-100002_UNITS = 'DEGREES'",
                "TKFRAME_-100002_AXES = ( 1, 2, 3 )", "TKFRAME_-100002_ANGLES = ( -5.0, 15.0, 70.0 )",
                "FRAME_TKC = -100003", "FRAME_-100003_NAME = 'TKC'",
                "FRAME_-100003_CLASS = 4", "FRAME_-100003_CLASS_ID = -100003",
                "FRAME_-100003_CENTER = 399", "TKFRAME_-100003_RELATIVE = 'ECLIPJ2000'",
                "TKFRAME_-100003_SPEC = 'ANGLES'", "TKFRAME_-100003_UNITS = 'DEGREES'",
                "TKFRAME_-100003_AXES = ( 3, 2, 1 )", "TKFRAME_-100003_ANGLES = ( 40.0, -25.0, 5.0 )",
            ]
        )
        epochs = [ets; -1000.0; 1000.0; 0.0; t0; t1; t1 + 10.0; -0.0]
        fails = Failures()
        for (targ, obs) in (
                    (-2000, -1000), (-2000, -3000), (-3000, -2000), (0, -2000), (-2000, 0),
                    (-4000, -1000), (-5000, -2000), (-2000, -5000),
                ), ref in ("J2000", "ECLIPJ2000", "IAU_EARTH", "TKB")
            q = state_query(ctx, targ, ref, "NONE", obs)
            for et in epochs
                check!(
                    fails, isequal(outcome(() -> q(et)), outcome(() -> spkezr(ctx, targ, et, ref, "NONE", obs))),
                    "state", targ, obs, ref, et
                )
                check!(
                    fails, isequal(outcome(() -> spkpos(q, et)), outcome(() -> spkpos(ctx, targ, et, ref, "NONE", obs))),
                    "position", targ, obs, ref, et
                )
                check!(
                    fails, isequal(outcome(() -> spkder(q, et, Val(2))), outcome(() -> spkder(ctx, targ, et, ref, obs, Val(2)))),
                    "jet", targ, obs, ref, et
                )
            end
        end
        # IAU_MARS has no constants in the pool: its link throws at each epoch.
        for (from, to) in (
                    ("TKB", "IAU_EARTH"), ("TKB", "TKA"), ("TKA", "TKB"), ("IAU_EARTH", "TKB"),
                    ("TKB", "ECLIPJ2000"), ("TKC", "J2000"), ("TKC", "TKB"), ("J2000", "TKC"),
                    ("J2000", "IAU_MARS"), ("TKB", "IAU_MARS"),
                ), et in epochs
            q = rotation_query(ctx, from, to)
            check!(
                fails, isequal(outcome(() -> q(et)), outcome(() -> pxform(ctx, from, to, et))),
                "pxform", from, to, et
            )
            check!(
                fails, isequal(outcome(() -> sxform(q, et)), outcome(() -> sxform(ctx, from, to, et))),
                "sxform", from, to, et
            )
        end
        @test isempty(fails)
    end

    @testset "the kernels of the query" begin
        q = state_query(ctx, "PROBE", "J2000", "NONE", "EARTH")
        before = q(0.0)
        local other = SpiceContext()
        add_ephemeris!(et -> probe(et) .+ 1, other, -1000, 399, "J2000", t0, t1)
        @test q(0.0) === before                       # another context: no effect
        add_ephemeris!(et -> probe(et) .+ 1, ctx, -1000, 399, "J2000", -10.0, 10.0; name = "later")
        @test q(0.0) === before                       # a later load: the query keeps its kernels
        @test state_query(ctx, "PROBE", "J2000", "NONE", "EARTH")(0.0)[1][1] == before[1][1] + 1
        unload(ctx, "later")
    end

    @testset "dual numbers and threads" begin
        q = state_query(ctx, "PROBE", "IAU_EARTH", "CN+S", "EARTH")
        g(et) = spkpos(q, et)[1]
        h(et) = spkpos(ctx, "PROBE", et, "IAU_EARTH", "CN+S", "EARTH")[1]
        @test ForwardDiff.derivative(g, 100.0) == ForwardDiff.derivative(h, 100.0)
        many = collect(range(t0 + 1000, t1 - 1000; length = 2000))
        serial = [q(et) for et in many]
        par = similar(serial)
        Threads.@threads for i in eachindex(many)
            par[i] = q(many[i])
        end
        @test par == serial
    end

    @testset "errors" begin
        kind(f) = try
            f()
            nothing
        catch e
            e isa SpiceError ? e.kind : rethrow()
        end
        @test kind(() -> state_query(ctx, "NO SUCH BODY", "J2000", "NONE", 399)) === :idcodenotfound
        @test kind(() -> state_query(ctx, -1000, "NO_FRAME", "NONE", 399)) === :unknownframe
        @test kind(() -> state_query(ctx, -1000, 1_234_567, "NONE", 399)) === :unknownframe
        @test kind(() -> rotation_query(ctx, "J2000", "NO_FRAME")) === :unknownframe
        @test state_query(ctx, -1000, "J2000", "NONE", 399) isa StateQuery
        q = state_query(ctx, -1000, "J2000", "NONE", 399)
        @test kind(() -> q(t1 + 10)) === kind(() -> spkezr(ctx, -1000, t1 + 10, "J2000", "NONE", 399))
    end
end

group("53. prepared queries on de440s and the binary PCK of the Moon", Kernels.FULL_SET_FILES...) do
    ctx = neo_context(Kernels.full_set())
    rng = MersenneTwister(SEED + 53)
    # Epochs inside and outside the coverage of de440s and of the binary PCK.
    epochs = [
        rand(rng, nsample(20, 200)) .* (DE440S_COVERAGE[2] - DE440S_COVERAGE[1]) .+ DE440S_COVERAGE[1];
        DE440S_COVERAGE[1] - 1.0e7; DE440S_COVERAGE[2] + 1.0e7; 0.0
    ]
    outcome(f) = try
        f()
    catch e
        e isa SpiceError ? (typeof(e), e.kind, e.msg) : rethrow()
    end
    fails = Failures()
    for (targ, obs) in ((399, 3), (301, 399), (399, 301), (4, 301), (0, 399), (10, 5)),
            ref in ("J2000", "ECLIPJ2000", "IAU_EARTH", "MOON_ME")
        q = state_query(ctx, targ, ref, "NONE", obs)
        for et in epochs
            check!(
                fails, isequal(outcome(() -> q(et)), outcome(() -> spkezr(ctx, targ, et, ref, "NONE", obs))),
                "state", targ, obs, ref, et
            )
            check!(
                fails, isequal(outcome(() -> spkpos(q, et)), outcome(() -> spkpos(ctx, targ, et, ref, "NONE", obs))),
                "position", targ, obs, ref, et
            )
            check!(
                fails, isequal(outcome(() -> spkder(q, et, Val(2))), outcome(() -> spkder(ctx, targ, et, ref, obs, Val(2)))),
                "jet", targ, obs, ref, et
            )
        end
    end
    for (from, to) in (("J2000", "MOON_PA"), ("MOON_ME", "IAU_EARTH"), ("MOON_PA", "MOON_ME")), et in epochs
        q = rotation_query(ctx, from, to)
        check!(fails, isequal(outcome(() -> q(et)), outcome(() -> pxform(ctx, from, to, et))), "pxform", from, to, et)
        check!(fails, isequal(outcome(() -> sxform(q, et)), outcome(() -> sxform(ctx, from, to, et))), "sxform", from, to, et)
    end
    @test isempty(fails)
end
