# Group 9: spacecraft clocks (SCLK type 1) against CSPICE.
#
# The SCLK kernel is written by `test/ckfixtures.jl`: clock -77 (two fields, two
# partitions, four coefficient records, parallel time TDB) and clock -78 (three fields
# with a first-field offset, three partitions, parallel time TDT).
#
# Error measures: seconds for ET, ticks for encoded SCLK; strings must be equal.

import ForwardDiff

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Almanac keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALMANAC_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

isdefined(@__MODULE__, :CKFix) || include("ckfixtures.jl")

# ET of the clocks is below 4e9 s (one unit in the last place 5e-7 s) at the end of the
# last partition and below 2e5 s (3e-11 s) in the part that the CK tests use. The two
# engines do the same operations in the same order, so the measured difference is 0;
# the thresholds leave room for a `sin` that differs in the last place (TDT clock).
const SCLK_TOL_ET = 1.0e-7          # seconds, over the whole clock range
const SCLK_TOL_TICKS = 1.0e-6       # ticks, continuous conversion over the whole clock range

"Largest encoded SCLK of the clock with partitions `pstart`, `pstop`."
sclk_mxtick(pstart, pstop) = sum(round.(pstop .- pstart))

"True if both engines throw, or both return and `same(a, b)` holds."
function agree(same, neo, oracle)
    a = try
        neo()
    catch e
        e isa Almanac.SpiceError || rethrow()
        nothing
    end
    b = try
        oracle()
    catch e
        e isa SPICE.SpiceError || rethrow()
        nothing
    end
    (a === nothing || b === nothing) && return a === nothing && b === nothing
    return same(a, b)
end

alloc_sct2e(ks, sc, x) = @allocated sct2e(ks, sc, x)
alloc_sce2c(ks, sc, x) = @allocated sce2c(ks, sc, x)
alloc_sce2t(ks, sc, x) = @allocated sce2t(ks, sc, x)

group("9. spacecraft clocks against CSPICE", Kernels.LSK) do
    rng = MersenneTwister(SEED + 9)
    set = [Kernels.LSK, CKFix.sclk_kernel()]
    oracle_load(set)
    ctx = neo_context(set)
    ks = Almanac.kernelset(ctx)
    clocks = (
        (CKFix.SC_TDB, CKFix.TDB_PSTART, CKFix.TDB_PSTOP, CKFix.TDB_COEF),
        (CKFix.SC_TDT, CKFix.TDT_PSTART, CKFix.TDT_PSTOP, CKFix.TDT_COEF),
    )

    @apart @testset "clock data" begin
        fails = Failures()
        for (sc, pstart, pstop, _) in clocks
            # CSPICE has no C wrapper of `sctype_`; the Fortran-style function is called.
            check!(fails, sctype(ctx, sc) == ccall((:sctype_, SPICE.libcspice), Cint, (Ref{Cint},), sc) == 1, "sctype", sc)
            ps, pe = scpart(ctx, sc)
            cs, ce = SPICE.scpart(sc)
            check!(fails, ps == cs == pstart, "partition starts", sc)
            check!(fails, pe == ce == pstop, "partition ends", sc)
        end
        @test isempty(fails)
        @test sort(collect(keys(ks.sclk))) == [CKFix.SC_TDT, CKFix.SC_TDB]
    end

    fails = Failures()
    for (sc, pstart, pstop, coef) in clocks
        mx = sclk_mxtick(pstart, pstop)
        # Random ticks over the whole range and over the part with the coefficient
        # records, and the special values: range ends, partition totals, record ticks.
        ticks = Float64[mx * rand(rng) for _ in 1:1500]
        append!(ticks, [coef[end, 1] * 1.2 * rand(rng) for _ in 1:1500])
        totals = cumsum(round.(pstop .- pstart))
        for s in (0.0, mx, totals..., coef[:, 1]...), d in (-1.0, -0.25, 0.0, 0.25, 1.0)
            0 <= s + d <= mx && push!(ticks, s + d)
        end
        det = Float64[]; dc = Float64[]; dt = Float64[]; rt = Float64[]
        strings_ok = true
        for t in ticks
            et = sct2e(ctx, sc, t)
            etc = SPICE.sct2e(sc, t)
            push!(det, abs(et - etc))
            push!(dc, abs(sce2c(ctx, sc, etc) - SPICE.sce2c(sc, etc)))
            push!(dt, abs(sce2t(ctx, sc, etc) - SPICE.sce2t(sc, etc)))
            push!(rt, abs(sce2c(ctx, sc, et) - t))
            s = scdecd(ctx, sc, t)
            strings_ok &= s == SPICE.scdecd(sc, t)
            strings_ok &= sce2s(ctx, sc, etc) == SPICE.sce2s(sc, etc)
            strings_ok &= scencd(ctx, sc, s) == SPICE.scencd(sc, s)
            strings_ok &= scs2e(ctx, sc, s) == SPICE.scs2e(sc, s)
            fields = last(split(s, '/'))
            strings_ok &= sctiks(ctx, sc, fields) == SPICE.sctiks(sc, fields)
            # `scfmt` of the tick count inside the partition
            k = sctiks(ctx, sc, fields)
            strings_ok &= scfmt(ctx, sc, k) == SPICE.scfmt(sc, k) == fields
        end
        a = report_quantiles("clock $sc: sct2e, seconds", det)
        b = report_quantiles("clock $sc: sce2c, ticks", dc)
        c = report_quantiles("clock $sc: sce2t, ticks", dt)
        r = report_quantiles("clock $sc: sce2c(sct2e(t)) - t, ticks", rt)
        within!(fails, a.max, SCLK_TOL_ET, "sct2e", sc)
        within!(fails, b.max, SCLK_TOL_TICKS, "sce2c", sc)
        check!(fails, c.max == 0.0, "sce2t", sc)
        # The round trip is exact to rounding, except in the 50 s jump of the parallel
        # time of clock -77 at its partition boundary, where two tick values share an ET.
        within!(fails, r.q50, SCLK_TOL_TICKS, "sce2c(sct2e(t)) - t, median", sc)
        check!(fails, strings_ok, "strings", sc)
    end
    @test isempty(fails)

    @apart @testset "clock strings as a user writes them" begin
        cases = (
            (
                CKFix.SC_TDB, [
                    "1/100.5", "100.5", "1/0000000100.005", "2/2000.000",
                    "2/ 2000 . 17", "  1/  77 12", "9999.255", "2/9999:255", "1/100", "100,3",
                    "2/2000-1", "1/.7", "1/7.", "2000.001", "1.5/100.5",
                ],
            ),
            (
                CKFix.SC_TDT, [
                    "1/1:0:0", "1:59:49", "1/000000854:20:01", "2/3:00:00",
                    "3/1:0:0", "2/5", "1/7 30", "1/7::3", "1/7:5:", "12:30",
                ],
            ),
        )
        fails = Failures()
        for (sc, list) in cases, s in list
            check!(fails, agree(==, () -> scencd(ctx, sc, s), () -> SPICE.scencd(sc, s)), "scencd", sc, repr(s))
            check!(fails, agree(==, () -> scs2e(ctx, sc, s), () -> SPICE.scs2e(sc, s)), "scs2e", sc, repr(s))
            f = String(last(split(s, '/')))
            check!(fails, agree(==, () -> sctiks(ctx, sc, f), () -> SPICE.sctiks(sc, f)), "sctiks", sc, repr(f))
        end
        @test isempty(fails)
        # At least one string of each clock is valid, so the comparison is not vacuous.
        @test scencd(ctx, CKFix.SC_TDB, "2/2000.001") == SPICE.scencd(CKFix.SC_TDB, "2/2000.001")
        @test scencd(ctx, CKFix.SC_TDT, "1/1:0:0") == 0.0
    end

    @apart @testset "errors" begin
        sc = CKFix.SC_TDB
        mx = sclk_mxtick(CKFix.TDB_PSTART, CKFix.TDB_PSTOP)
        @test throws_kind(() -> sct2e(ctx, sc, -1.0), :valueoutofrange)
        @test throws_kind(() -> sct2e(ctx, sc, mx + 1.0), :valueoutofrange)
        @test_throws SPICE.SpiceError SPICE.sct2e(sc, -1.0)
        @test_throws SPICE.SpiceError SPICE.sct2e(sc, mx + 1.0)
        @test throws_kind(() -> sce2c(ctx, sc, -10.0), :valueoutofrange)
        @test_throws SPICE.SpiceError SPICE.sce2c(sc, -10.0)
        @test throws_kind(() -> sce2c(ctx, sc, 1.0e12), :valueoutofrange)
        @test throws_kind(() -> sce2t(ctx, sc, 1.0e12), :valueoutofrange)
        # `sce2t` extends the first record to earlier epochs; the result is then negative.
        @test throws_kind(() -> sce2t(ctx, sc, -10.0), :valueoutofrange)
        @test throws_kind(() -> scdecd(ctx, sc, -1.0), :valueoutofrange)
        @test throws_kind(() -> scdecd(ctx, sc, mx + 1.0), :valueoutofrange)
        @test throws_kind(() -> scfmt(ctx, sc, -1.0), :valueoutofrange)
        @test throws_kind(() -> sctiks(ctx, sc, "  "), :invalidsclkstring)
        @test throws_kind(() -> sctiks(ctx, sc, "1.2.3"), :invalidsclkstring)
        @test throws_kind(() -> sctiks(ctx, sc, "12.x"), :invalidsclkstring)
        @test throws_kind(() -> sctiks(ctx, CKFix.SC_TDT, "0:0:0"), :invalidsclkstring)
        @test throws_kind(() -> scencd(ctx, sc, "/100.5"), :badpartnumber)
        @test throws_kind(() -> scencd(ctx, sc, "x/100.5"), :badpartnumber)
        @test throws_kind(() -> scencd(ctx, sc, "3/100.5"), :badpartnumber)
        @test throws_kind(() -> scencd(ctx, sc, "1/"), :invalidsclkstring)
        # no fields and no partition: the fields are checked first, as CSPICE (F_SCLK0)
        @test throws_kind(() -> scencd(ctx, sc, " /"), :invalidsclkstring)
        # Regression: a slash in the first column is BADPARTNUMBER before the fields are
        # read, and a blank partition number is read after the fields, as in scencd_.
        for str in ("/", "/  ", "  /x", "  /100.5")
            ccall((:scencd_c, SPICE.libcspice), Cvoid, (Cint, Cstring, Ref{Cdouble}), sc, str, Ref(0.0))
            @test throws_kind(() -> scencd(ctx, sc, str), cspice_kind())
        end
        # A character of several bytes in a field is an error of the string, not of the code.
        @test all(
            s -> try
                scencd(ctx, sc, s); false
            catch e
                e isa Almanac.SpiceError
            end,
            ("1/é:2", "1/100:4é", "é")
        )
        @test throws_kind(() -> scencd(ctx, sc, "1/10001.0"), :notinpart)
        @test throws_kind(() -> scencd(ctx, CKFix.SC_TDT, "2/1:0:0"), :notinpart)
        @test throws_kind(() -> scencd(ctx, CKFix.SC_TDT, "900000001:0:0"), :nopartition)
        @test throws_kind(() -> sct2e(ctx, -5, 0.0), :kernelvarnotfound)
        @test throws_kind(() -> sctype(ctx, -5), :kernelvarnotfound)
        @test throws_kind(() -> scpart(ctx, -5), :kernelvarnotfound)
        # A clock of another type than 1 has no conversion, as in CSPICE.
        other = neo_context([CKFix.sclk_kernel()])
        lmpool(other, ["SCLK_DATA_TYPE_90 = ( 2 )"])
        @test sctype(other, -90) == 2
        @test throws_kind(() -> sct2e(other, -90, 0.0), :notsupported)
        # A clock of type 1 whose partition start and end lists have different lengths
        # (`zzscup01` of CSPICE signals NUMPARTSUNEQUAL), and one that has its type
        # only. Both are set in the pool of the two engines and removed from that of
        # CSPICE afterwards.
        lines = [
            "SCLK_DATA_TYPE_91 = ( 1 )", "SCLK01_TIME_SYSTEM_91 = ( 1 )",
            "SCLK01_N_FIELDS_91 = ( 2 )", "SCLK01_MODULI_91 = ( 4294967296 256 )",
            "SCLK01_OFFSETS_91 = ( 0 0 )", "SCLK01_OUTPUT_DELIM_91 = ( 1 )",
            "SCLK_PARTITION_START_91 = ( 0.0 512000.0 )", "SCLK_PARTITION_END_91 = ( 2560000.0 )",
            "SCLK01_COEFFICIENTS_91 = ( 0.0 0.0 1.0 )", "SCLK_DATA_TYPE_92 = ( 1 )",
        ]
        lmpool(other, lines)
        SPICE.lmpool(lines)
        try
            for (id, kind) in ((-91, :numpartsunequal), (-92, :kernelvarnotfound))
                @test throws_kind(() -> sct2e(other, id, 0.0), kind)
                ccall((:sct2e_c, SPICE.libcspice), Cvoid, (Cint, Cdouble, Ref{Cdouble}), id, 0.0, Ref(0.0))
                @test cspice_kind() === kind
            end
        finally
            foreach(l -> SPICE.dvpool(String(strip(first(split(l, '='))))), lines)
        end
        # Regression: the type of a clock is the first value of SCLK_DATA_TYPE_<n>,
        # rounded, as `gipool_` of `scty01_` reads it; a clock was built only for one
        # value equal to 1.
        lines = [
            "SCLK_DATA_TYPE_93 = ( 1.2 7 )", "SCLK01_TIME_SYSTEM_93 = ( 1 )",
            "SCLK01_N_FIELDS_93 = ( 2 )", "SCLK01_MODULI_93 = ( 4294967296 256 )",
            "SCLK01_OFFSETS_93 = ( 0 0 )", "SCLK01_OUTPUT_DELIM_93 = ( 1 )",
            "SCLK_PARTITION_START_93 = ( 0.0 )", "SCLK_PARTITION_END_93 = ( 2560000.0 )",
            "SCLK01_COEFFICIENTS_93 = ( 0.0 0.0 1.0 )",
        ]
        lmpool(other, lines)
        SPICE.lmpool(lines)
        try
            @test sctype(other, -93) == 1
            @test sct2e(other, -93, 1000.0) == SPICE.sct2e(-93, 1000.0)
        finally
            foreach(l -> SPICE.dvpool(String(strip(first(split(l, '='))))), lines)
        end
        # The TDT clock needs the leap-second kernel; the TDB clock does not.
        bare = neo_context([CKFix.sclk_kernel()])
        @test throws_kind(() -> sct2e(bare, CKFix.SC_TDT, 1000.0), :missingtimeinfo)
        @test sct2e(bare, sc, 1000.0) == sct2e(ctx, sc, 1000.0)
    end

    @apart @testset "scfmt where ticks / modulus rounds up to a whole number" begin
        # A field of 1e9 (nanoseconds). Above about 9e18 ticks the quotient x / 1e9 can
        # round up to k, so that `DMOD` of f2c (x - y trunc(x / y)) is negative and its
        # magnitude is printed; the exact remainder counts y - 1024 ticks twice. Found
        # 2026-11-17: Almanac wrote "999999488" where CSPICE writes "000001024".
        lines = [
            "SCLK_KERNEL_ID = ( @2000-01-01 )", "SCLK_DATA_TYPE_91 = ( 1 )",
            "SCLK01_TIME_SYSTEM_91 = ( 1 )", "SCLK01_N_FIELDS_91 = ( 2 )",
            "SCLK01_MODULI_91 = ( 100000000000 1000000000 )", "SCLK01_OFFSETS_91 = ( 0 0 )",
            "SCLK01_OUTPUT_DELIM_91 = ( 2 )", "SCLK_PARTITION_START_91 = ( 0.0 )",
            "SCLK_PARTITION_END_91 = ( 1.0E20 )", "SCLK01_COEFFICIENTS_91 = ( 0.0 0.0 1.0 )",
        ]
        names = [String(strip(first(split(l, '=')))) for l in lines[2:end]]
        nano = neo_context(String[])
        lmpool(nano, lines)
        try
            SPICE.lmpool(lines)
            cases = Float64[]
            for k in 9_000_000_000:9_000_100_000
                x = prevfloat(k * 1.0e9)
                trunc(x / 1.0e9) * 1.0e9 > x && push!(cases, x)
                length(cases) == 20 && break
            end
            @test length(cases) == 20
            @test all(x -> scfmt(nano, -91, x) == SPICE.scfmt(-91, x), cases)
            @test scfmt(nano, -91, 8.723322218999999e18) == "08723322219:000001024"
        finally
            foreach(SPICE.dvpool, names)
        end
    end

    @apart @testset "derivatives pass through the conversions" begin
        fails = Failures()
        for (sc, _, _, coef) in clocks
            tikmsc = sc == CKFix.SC_TDB ? 256.0 : CKFix.TDT_TIKMSC
            for k in 1:(size(coef, 1) - 1)
                t = coef[k, 1] + (coef[k + 1, 1] - coef[k, 1]) * (0.2 + 0.6 * rand(rng))
                x = sct2e(ctx, sc, Almanac.taylor_var(t, Val(2)))
                e0, e1, e2 = Almanac.derivatives(x)
                check!(fails, e0 == sct2e(ctx, sc, t), "sct2e value", sc, t)
                h = 64.0
                fd = (sct2e(ctx, sc, t + h) - sct2e(ctx, sc, t - h)) / 2h
                check!(fails, isapprox(e1, fd; rtol = 1.0e-9), "sct2e derivative against differences", sc, t)
                # seconds per tick: the record rate (per count) over the ticks in a count
                check!(fails, isapprox(e1, coef[k, 3] / tikmsc; rtol = sc == CKFix.SC_TDB ? 1.0e-15 : 1.0e-7), "sct2e derivative = record rate", sc, t)
                sc == CKFix.SC_TDB && check!(fails, e2 == 0.0, "sct2e second derivative", sc, t)
                y = sce2c(ctx, sc, Almanac.taylor_var(e0, Val(2)))
                c0, c1, _ = Almanac.derivatives(y)
                check!(fails, isapprox(c0, t; atol = 1.0e-6), "sce2c(sct2e(t)) = t", sc, t)
                check!(fails, isapprox(c1 * e1, 1.0; rtol = 1.0e-12), "sce2c derivative = 1 / sct2e derivative", sc, t)
            end
        end
        @test isempty(fails)
        # A dual number passes as well, with the derivative of the Taylor number.
        sc, coef = clocks[1][1], clocks[1][4]
        t = (coef[1, 1] + coef[2, 1]) / 2
        d = ForwardDiff.derivative(s -> sct2e(ctx, sc, s), t)
        @test d ≈ Almanac.derivatives(sct2e(ctx, sc, Almanac.taylor_var(t, Val(1))))[2] rtol = 1.0e-15
    end

    @apart @testset "no allocation, inferred types" begin
        fails = Failures()
        for (sc, _, _, coef) in clocks
            t = coef[2, 1] * 0.7
            et = sct2e(ks, sc, t)
            tv = Almanac.taylor_var(t, Val(2))
            ev = Almanac.taylor_var(et, Val(2))
            check!(fails, CKFix.warm(alloc_sct2e, ks, sc, t) == 0, "sct2e", sc)
            check!(fails, CKFix.warm(alloc_sce2c, ks, sc, et) == 0, "sce2c", sc)
            check!(fails, CKFix.warm(alloc_sce2t, ks, sc, et) == 0, "sce2t", sc)
            check!(fails, CKFix.warm(alloc_sct2e, ks, sc, tv) == 0, "sct2e Taylor", sc)
            check!(fails, CKFix.warm(alloc_sce2c, ks, sc, ev) == 0, "sce2c Taylor", sc)
            @test (@inferred sct2e(ks, sc, t)) isa Float64
            @test (@inferred sce2c(ks, sc, et)) isa Float64
            @test (@inferred sce2t(ks, sc, et)) isa Float64
            @test (@inferred sct2e(ks, sc, tv)) isa Almanac.Taylor{2, Float64}
            @test (@inferred sce2c(ks, sc, ev)) isa Almanac.Taylor{2, Float64}
            @test (@inferred sct2e(ctx, sc, 1000)) isa Float64
        end
        @test isempty(fails)
    end

    @apart @testset "forms on the default context" begin
        kclear()
        foreach(furnsh, set)
        sc = CKFix.SC_TDT
        @test sct2e(sc, 5000.0) == sct2e(ctx, sc, 5000.0)
        @test sce2c(sc, 100.0) == sce2c(ctx, sc, 100.0)
        @test sce2t(sc, 100.0) == sce2t(ctx, sc, 100.0)
        @test sce2s(sc, 100.0) == sce2s(ctx, sc, 100.0)
        @test scs2e(sc, "1/7:30:00") == scs2e(ctx, sc, "1/7:30:00")
        @test scencd(sc, "1/7:30:00") == scencd(ctx, sc, "1/7:30:00")
        @test scdecd(sc, 5000.0) == scdecd(ctx, sc, 5000.0)
        @test sctiks(sc, "7:30:00") == sctiks(ctx, sc, "7:30:00")
        @test scfmt(sc, 5000.0) == scfmt(ctx, sc, 5000.0)
        @test scpart(sc) == scpart(ctx, sc)
        @test sctype(sc) == 1
        kclear()
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
