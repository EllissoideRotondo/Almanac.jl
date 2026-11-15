# Group 53: prepared queries (`state_query`, `rotation_query`, src/queries.jl).
#
# Invariant layer: a query gives what `spkezr`, `spkpos`, `pxform` and `sxform` give for
# the same arguments, to the last bit, for every aberration correction, an inertial and
# a rotating frame, names and codes. Then: the query keeps its kernels, dual-number
# epochs, threads, errors. Bodies and frames come from functions (`add_ephemeris!`) and
# from pool constants, so no kernel file of NAIF is needed.

using ForwardDiff

@testset "prepared queries" begin
    earth(et) = [1.5e8 + 1.0e-3 * et, -2.0e7 + 29.0 * et, 3.0e6 - 0.5 * et, 1.0e-3, 29.0, -0.5]
    probe(et) = [7000.0 + 0.002 * et, -1200.0 + 7.4 * sin(1e-3 * et), 300.0, 0.002,
        7.4e-3 * cos(1e-3 * et), 0.0]
    t0, t1 = -86400.0, 86400.0
    ctx = SpiceContext()
    pdpool(ctx, "BODY399_POLE_RA", [10.0, 0.1, 0.0])
    pdpool(ctx, "BODY399_POLE_DEC", [80.0, -0.05, 0.0])
    pdpool(ctx, "BODY399_PM", [190.0, 360.9856, 0.0])
    add_ephemeris!(earth, ctx, 399, 0, "J2000", t0, t1)
    add_ephemeris!(probe, ctx, -1000, "EARTH", "J2000", t0, t1)
    boddef(ctx, "PROBE", -1000)
    ets = range(t0 + 1000, t1 - 1000; length=5)

    @testset "the results of spkezr and spkpos, bit for bit" begin
        for corr in ("NONE", "LT", "LT+S", "CN", "CN+S", "XLT", "XLT+S", "XCN", "XCN+S"),
                ref in ("J2000", "ECLIPJ2000", "IAU_EARTH"),
                (targ, obs) in (("PROBE", "EARTH"), (399, -1000), ("EARTH", "SSB"))
            q = state_query(ctx, targ, ref, corr, obs)
            for et in ets
                @test q(et) === spkezr(ctx, targ, et, ref, corr, obs)
                @test spkpos(q, et) === spkpos(ctx, targ, et, ref, corr, obs)
            end
        end
        q = state_query(ctx, "PROBE", "J2000", Aberration("LT+S"), "EARTH")
        @test spkezr.(q, ets) == [spkezr(ctx, "PROBE", et, "J2000", "LT+S", "EARTH") for et in ets]
    end

    @testset "the results of pxform and sxform, bit for bit" begin
        for (from, to) in (("J2000", "IAU_EARTH"), ("IAU_EARTH", "ECLIPJ2000"), (1, 17))
            q = rotation_query(ctx, from, to)
            for et in ets
                @test q(et) === pxform(ctx, from, to, et)
                @test sxform(q, et) === sxform(ctx, from, to, et)
            end
        end
    end

    @testset "the kernels of the query" begin
        q = state_query(ctx, "PROBE", "J2000", "NONE", "EARTH")
        before = q(0.0)
        local other = SpiceContext()
        add_ephemeris!(et -> probe(et) .+ 1, other, -1000, 399, "J2000", t0, t1)
        @test q(0.0) === before                       # another context: no effect
        add_ephemeris!(et -> probe(et) .+ 1, ctx, -1000, 399, "J2000", -10.0, 10.0; name="later")
        @test q(0.0) === before                       # a later load: the query keeps its kernels
        @test state_query(ctx, "PROBE", "J2000", "NONE", "EARTH")(0.0)[1][1] == before[1][1] + 1
        unload(ctx, "later")
    end

    @testset "dual numbers and threads" begin
        q = state_query(ctx, "PROBE", "IAU_EARTH", "CN+S", "EARTH")
        g(et) = spkpos(q, et)[1]
        h(et) = spkpos(ctx, "PROBE", et, "IAU_EARTH", "CN+S", "EARTH")[1]
        @test ForwardDiff.derivative(g, 100.0) == ForwardDiff.derivative(h, 100.0)
        many = collect(range(t0 + 1000, t1 - 1000; length=2000))
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
