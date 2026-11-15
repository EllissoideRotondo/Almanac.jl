# Group 52: bodies and frames that come from functions (`add_ephemeris!`, `add_frame!`,
# src/providers.jl).
#
# Invariant layer: each provider describes linear motion, which an SPK of type 9 of
# degree 1 holds exactly, so a context with the providers and a context with the same
# motion written to an SPK file must give the same states, to rounding, for every
# correction, in an inertial and a rotating frame, and for the derivatives that the
# corrections need. Then the rules of a file: priority, coverage, unload, kclear, the
# listings; the errors; dual-number epochs; an event search; queries from threads.
# No kernel file of NAIF is needed.

using ForwardDiff, LinearAlgebra, StaticArrays

@testset "ephemeris providers" begin
    dir = mktempdir()
    # Earth relative to the barycenter and a probe relative to the Earth, as linear motion.
    earth(et) = [1.5e8 + 1.0e-3 * et, -2.0e7 + 29.0 * et, 3.0e6 - 0.5 * et, 1.0e-3, 29.0, -0.5]
    probe(et) = [7000.0 + 0.002 * et, -1200.0 + 7.4 * et, 300.0 + 0.1 * et, 0.002, 7.4, 0.1]
    t0, t1 = -86400.0, 86400.0
    function rotating!(ctx)          # IAU_EARTH from three pool variables
        pdpool(ctx, "BODY399_POLE_RA", [10.0, 0.1, 0.0])
        pdpool(ctx, "BODY399_POLE_DEC", [80.0, -0.05, 0.0])
        pdpool(ctx, "BODY399_PM", [190.0, 360.9856, 0.0])
        return ctx
    end

    spk = joinpath(dir, "linear.bsp")
    let h = spkopn(spk, "linear", 0), ep = [t0, t1]
        for (body, center, f) in ((399, 0, earth), (-1000, 399, probe))
            spkw09(h, body, center, "J2000", t0, t1, "linear", 1, reduce(hcat, f.(ep)), ep)
        end
        spkcls(h)
    end
    file = rotating!(SpiceContext())
    furnsh(file, spk)
    prov = rotating!(SpiceContext())
    add_ephemeris!(earth, prov, 399, 0, "J2000", t0, t1; name="earth")
    add_ephemeris!(probe, prov, -1000, "EARTH", "J2000", t0, t1)

    @testset "the same states as the SPK" begin
        # A corrected state is a difference of states relative to the barycenter, at
        # 1.6e8 km: the two contexts differ by rounding at that size, 2 ulps of 1.6e8 km
        # in position and 4e-12 km/s in velocity measured; the bounds are twice that.
        ap, av, al = 4 * eps(1.6e8), 1.0e-11, 4 * eps(1.6e8) / Almanac.CLIGHT
        ets = range(t0 + 1000, t1 - 1000; length=7)
        for corr in ("NONE", "LT", "LT+S", "CN", "CN+S", "XLT", "XLT+S", "XCN+S"),
                ref in ("J2000", "ECLIPJ2000", "IAU_EARTH"), (targ, obs) in ((-1000, 399), (399, -1000))
            for et in ets
                a, la = spkezr(file, targ, et, ref, corr, obs)
                b, lb = spkezr(prov, targ, et, ref, corr, obs)
                @test maximum(abs, b[1:3] - a[1:3]) <= ap && maximum(abs, b[4:6] - a[4:6]) <= av &&
                      abs(lb - la) <= al
            end
        end
        @test spkpos(prov, -1000, 10.0, "J2000", "NONE", 399)[1] ≈ probe(10.0)[1:3] rtol = 1e-15
        @test spkgeo(prov, -1000, 10.0, "J2000", 0)[1] ≈ probe(10.0) + earth(10.0) rtol = 1e-15
        # order 2 and 3 through Taylor numbers: the velocity is constant
        @test all(iszero, spkder(prov, -1000, 5.0, "J2000", 399, Val(3))[3])
        # the result types do not depend on the function (a regression: an untyped
        # provider made `spkgeo` of every SPK type type-unstable)
        @test (@inferred spkgeo(prov, -1000, 5.0, 1, 399)) isa Tuple{SVector{6,Float64},Float64}
        @test (@inferred spkder(prov, -1000, 5.0, 1, 399, Val(3))) isa NTuple{4,SVector{3,Float64}}
    end

    @testset "dual-number epochs" begin
        g(et) = spkpos(prov, -1000, et, "IAU_EARTH", "LT+S", 399)[1]
        h(et) = spkpos(file, -1000, et, "IAU_EARTH", "LT+S", 399)[1]
        @test ForwardDiff.derivative(g, 100.0) ≈ ForwardDiff.derivative(h, 100.0) rtol = 1e-12
    end

    @testset "the rules of a file" begin
        ctx = SpiceContext()
        furnsh(ctx, spk)
        @test spkezr(ctx, -1000, 0.0, "J2000", "NONE", 399)[1] == probe(0.0)
        shifted(et) = probe(et) .+ [1.0, 0, 0, 0, 0, 0]
        add_ephemeris!(shifted, ctx, -1000, 399, "J2000", -10.0, 10.0; name="shift")
        @test spkezr(ctx, -1000, 0.0, "J2000", "NONE", 399)[1][1] == probe(0.0)[1] + 1    # later wins
        @test spkezr(ctx, -1000, 20.0, "J2000", "NONE", 399)[1] == probe(20.0)            # outside: the file
        @test kinfo(ctx, "shift")[1] == "SPK" && ktotal(ctx, "SPK") == 2
        file2, descr, ident = spksfs(ctx, -1000, 0.0)
        @test (file2, ident) == ("shift", "shift")
        @test spkpvn(ctx, file2, descr, 0.0) == (1, shifted(0.0), 399)
        furnsh(ctx, spk)                                                                  # a later file wins
        @test spkezr(ctx, -1000, 0.0, "J2000", "NONE", 399)[1] == probe(0.0)
        unload(ctx, spk)
        @test spkezr(ctx, -1000, 0.0, "J2000", "NONE", 399)[1][1] == probe(0.0)[1] + 1
        unload(ctx, "shift")
        @test spkezr(ctx, -1000, 0.0, "J2000", "NONE", 399)[1] == probe(0.0)
        add_ephemeris!(shifted, ctx, -1000, 399, "J2000", -10.0, 10.0)
        @test kdata(ctx, 2)[1:2] == ("ephemeris of body -1000", "SPK")
        kclear(ctx)
        @test ktotal(ctx) == 0
        @test_throws SpiceError spkezr(ctx, -1000, 0.0, "J2000", "NONE", 399)
    end

    @testset "errors" begin
        ctx = SpiceContext()
        kind(f) = try
            f()
            nothing
        catch e
            e isa SpiceError ? e.kind : rethrow()
        end
        ok(et) = probe(et)
        @test kind(() -> add_ephemeris!(ok, ctx, "NO SUCH BODY", 399, "J2000", 0, 1)) === :idcodenotfound
        @test kind(() -> add_ephemeris!(ok, ctx, -1000, 399, "NO_FRAME", 0, 1)) === :unknownframe
        @test kind(() -> add_ephemeris!(ok, ctx, -1000, 399, "J2000", 1, 0)) === :baddescrtimes
        @test kind(() -> add_ephemeris!(ok, ctx, 399, "EARTH", "J2000", 0, 1)) === :barycenterephem
        @test ktotal(ctx) == 0
        add_ephemeris!(et -> [1.0, 2.0], ctx, -1000, 399, "J2000", 0, 1)
        @test kind(() -> spkezr(ctx, -1000, 0.5, "J2000", "NONE", 399)) === :badarraysize
        @test kind(() -> spkezr(ctx, -1000, 2.0, "J2000", "NONE", 399)) === :spkinsuffdata
    end

    @testset "an event search" begin
        # The probe's distance from the Earth falls below 300000 km once in [t0, t1].
        w = gfdist(prov, -1000, "NONE", 399, ">", 3.0e5, 0.0, 3600.0, Window([t0 + 1, t1 - 1]))
        crossing = w[2]
        @test norm(spkpos(file, -1000, crossing, "J2000", "NONE", 399)[1]) ≈ 3.0e5 rtol = 1e-9
    end

    @testset "queries from threads" begin
        ets = collect(range(t0 + 1000, t1 - 1000; length=4000))
        serial = [spkezr(prov, -1000, et, "IAU_EARTH", "CN+S", 399)[1] for et in ets]
        par = similar(serial)
        Threads.@threads for i in eachindex(ets)
            par[i] = spkezr(prov, -1000, ets[i], "IAU_EARTH", "CN+S", 399)[1]
        end
        @test par == serial
    end

    @testset "frames" begin
        ω = 2π / 86400                                        # one turn per day about Z
        # `local`: the testset assigns `s` below, and a closure that assigns an outer
        # local shares it between the threads that call it.
        spin(et) = (local c = cos(ω * et); local sn = sin(ω * et); [c -sn 0; sn c 0; 0 0 1])
        tilt = [1 0 0; 0 cosd(30) -sind(30); 0 sind(30) cosd(30)]
        ctx = SpiceContext()
        add_frame!(spin, ctx, "SPIN", 1_400_001, 399, "J2000")
        add_frame!(et -> tilt, ctx, "TILT", 1_400_002, 399, "SPIN"; path="tilt frame")
        lmpool(ctx, ["FRAME_TK_TILT = 1400003", "FRAME_1400003_NAME = 'TK_TILT'",
            "FRAME_1400003_CLASS = 4", "FRAME_1400003_CLASS_ID = 1400003",
            "FRAME_1400003_CENTER = 399", "TKFRAME_1400003_RELATIVE = 'SPIN'",
            "TKFRAME_1400003_SPEC = 'MATRIX'",
            "TKFRAME_1400003_MATRIX = ( $(join(vec(tilt), ' ')) )"])
        et = 5000.0
        @test pxform(ctx, "SPIN", "J2000", et) ≈ spin(et) rtol = 1e-15
        @test pxform(ctx, "J2000", "SPIN", et) ≈ transpose(spin(et)) rtol = 1e-15
        @test pxform(ctx, "TILT", "J2000", et) ≈ spin(et) * tilt rtol = 1e-15
        @test pxform(ctx, "TILT", "J2000", et) ≈ pxform(ctx, "TK_TILT", "J2000", et) rtol = 1e-15
        x = sxform(ctx, "TILT", "J2000", et)
        dspin = ω * [-sin(ω * et) -cos(ω * et) 0; cos(ω * et) -sin(ω * et) 0; 0 0 0]
        @test x[4:6, 1:3] ≈ dspin * tilt rtol = 1e-14
        @test x ≈ sxform(ctx, "TK_TILT", "J2000", et) rtol = 1e-14
        @test ForwardDiff.derivative(t -> pxform(ctx, "TILT", "J2000", t)[2, 1], et) ≈ x[5, 1] rtol = 1e-14
        @test frinfo(ctx, 1_400_001)[1:2] == (399, 7)
        @test namfrm(ctx, "SPIN") == 1_400_001 && frmnam(ctx, 1_400_002) == "TILT"
        # a state in the frame, a body defined in the frame
        add_ephemeris!(et -> [7000.0, 0, 0, 0, 0, 0], ctx, -1000, 399, "SPIN", t0, t1)
        r = spkpos(ctx, -1000, et, "J2000", "NONE", 399)[1]
        @test r ≈ spin(et) * [7000.0, 0, 0] rtol = 1e-15
        s = spkezr(ctx, -1000, et, "SPIN", "NONE", 399)[1]
        @test isapprox(s[1:3], [7000.0, 0, 0]; rtol=1e-15)
        @test maximum(abs, s[4:6]) < 1e-12
        # the load list
        @test kinfo(ctx, "tilt frame")[1] == "TEXT"
        unload(ctx, "tilt frame")
        @test_throws SpiceError pxform(ctx, "TILT", "J2000", et)
        @test pxform(ctx, "SPIN", "J2000", et) ≈ spin(et)
        # threads
        ets = collect(range(0.0, 86400.0; length=2000))
        serial = [sxform(ctx, "SPIN", "J2000", t) for t in ets]
        par = similar(serial)
        Threads.@threads for i in eachindex(ets)
            par[i] = sxform(ctx, "SPIN", "J2000", ets[i])
        end
        @test par == serial
        kclear(ctx)
        @test_throws SpiceError pxform(ctx, "SPIN", "J2000", et)
    end

    @testset "frame errors" begin
        ctx = SpiceContext()
        kind(f) = try
            f()
            nothing
        catch e
            e isa SpiceError ? e.kind : rethrow()
        end
        id3(et) = [1.0 0 0; 0 1 0; 0 0 1]
        @test kind(() -> add_frame!(id3, ctx, "IAU_EARTH", 1_400_001, 399, "J2000")) === :badframe
        @test kind(() -> add_frame!(id3, ctx, "X", 1, 399, "J2000")) === :badframe
        @test kind(() -> add_frame!(id3, ctx, "X", 1_400_001, 399, 1_400_001)) === :badframe
        @test kind(() -> add_frame!(id3, ctx, "X", 1_400_001, 399, "NO_FRAME")) === :unknownframe
        @test kind(() -> add_frame!(id3, ctx, "X", 1_400_001, "NO BODY", "J2000")) === :idcodenotfound
        @test kind(() -> add_frame!(id3, ctx, " ", 1_400_001, 399, "J2000")) === :emptystring
        @test ktotal(ctx) == 0
        add_frame!(et -> [1.0 0; 0 1], ctx, "BAD", 1_400_001, 399, "J2000")
        @test kind(() -> pxform(ctx, "BAD", "J2000", 0.0)) === :badarraysize
    end
end

