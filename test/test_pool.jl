# Group 5: kernel pool and names against CSPICE. The time conversions (`utc2et`, `et2utc`,
# `deltet`) are in the group `time`, which compares them with CSPICE around every leap
# second; this group keeps the error of a conversion without a leap-second kernel.
#
# Conventions for an absent item, as the prototype source documents them: `gdpool`,
# `gcpool`, `gipool` and `bodvrd` throw `SpiceError(:kernelvarnotfound)` (SPICE.jl returns
# `nothing`, or throws for `bodvrd`); `bodn2c`, `bods2c` and `bodc2n` return `nothing`;
# `namfrm` returns 0 and `frmnam` returns "", as CSPICE does.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Almanac keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALMANAC_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

"Pool variable names of CSPICE that match `template`."
function oracle_names(template::AbstractString)
    out = String[]
    start = 1
    while true
        part = SPICE.gnpool(template, start, 500)
        part === nothing && break
        append!(out, part)
        length(part) < 500 && break
        start += 500
    end
    return unique!(out)
end

group(
    "5. kernel pool and names against CSPICE", Kernels.TEXT_SET_FILES...,
    Kernels.MOON_TF_OLD
) do
    rng = MersenneTwister(SEED + 5)
    set = [Kernels.text_set()..., Kernels.MOON_TF_OLD]
    have_gm = isfile(Kernels.GM_TPC)
    have_gm ? push!(set, Kernels.GM_TPC) :
        skip_notice("5. GM values", "missing kernel: " * Kernels.GM_TPC)
    oracle_load(set)
    ctx = neo_context(set)

    @apart @testset "every pool variable of CSPICE" begin
        names = oracle_names("*")
        println("  ", length(names), " pool variables in CSPICE")
        @test length(names) > 500
        nnum = 0; nchr = 0; bad = String[]
        for name in names
            n, type = SPICE.dtpool(name)
            if type === :N
                ref = SPICE.gdpool(name; room = n)
                got = try
                    gdpool(ctx, name)
                catch e
                    e
                end
                (got isa Vector{Float64} && got == ref) || push!(bad, name)
                nnum += 1
            else
                ref = SPICE.gcpool(name; room = n, lenout = 1024)
                got = try
                    gcpool(ctx, name)
                catch e
                    e
                end
                (got isa Vector{String} && got == ref) || push!(bad, name)
                nchr += 1
            end
        end
        println("  numeric: ", nnum, ", character: ", nchr, ", different: ", length(bad))
        isempty(bad) || println("  first different variables: ", first(bad, 10))
        @test isempty(bad)
        @test nnum > 0 && nchr > 0
    end

    @apart @testset "gdpool, gcpool, gipool samples" begin
        @test gdpool(ctx, "BODY399_RADII") == [6378.1366, 6378.1366, 6356.7519]
        @test gdpool(ctx, "BODY399_RADII") isa Vector{Float64}
        fails = Failures()
        for name in (
                "BODY301_PM", "BODY499_POLE_RA", "BODY5_NUT_PREC_ANGLES",
                "DELTET/DELTA_AT", "DELTET/DELTA_T_A", "TKFRAME_31009_ANGLES",
                "BODY399_N_GEOMAG_CTR_DIPOLE_LAT",
            )
            check!(fails, gdpool(ctx, name) == SPICE.gdpool(name; room = 1000), "gdpool", name)
        end
        for name in (
                "FRAME_31008_NAME", "TKFRAME_31009_SPEC", "TKFRAME_31009_UNITS",
                "TKFRAME_31009_RELATIVE", "TKFRAME_31000_RELATIVE", "FRAME_31006_NAME",
            )
            check!(fails, gcpool(ctx, name) == SPICE.gcpool(name), "gcpool", name)
        end
        for name in (
                "FRAME_MOON_PA_DE440", "TKFRAME_31009_AXES", "FRAME_31008_CLASS",
                "BODY399_RADII",
            )
            check!(fails, gipool(ctx, name) == SPICE.gipool(name), "gipool", name)
        end
        @test isempty(fails)
        # `NINT` rounds a tie away from zero, inside the 32-bit range (ends included);
        # outside it, NaN included, CSPICE signals SPICE(INTOUTOFRANGE). Found 2026-11-17:
        # `round(Int, x)` gave 2 for 2.5, and 3000000000 for 3.0e9.
        ties = [2.5, -2.5, 0.5, -0.5, 2147483646.5, -2147483647.5, 2147483647.0, -2147483648.0]
        try
            pdpool(ctx, "GI_TIES", ties); SPICE.pdpool("GI_TIES", ties)
            @test gipool(ctx, "GI_TIES") == [3, -3, 1, -1, 2147483647, -2147483648, 2147483647, -2147483648]
            @test gipool(ctx, "GI_TIES") == SPICE.gipool("GI_TIES")
            for (k, x) in enumerate((2147483647.4, -2147483648.3, 3.0e9, NaN))
                pdpool(ctx, "GI_BIG$k", [1.0, x])
                @test throws_kind(() -> gipool(ctx, "GI_BIG$k"), :intoutofrange)
            end
        finally
            SPICE.dvpool("GI_TIES")
            foreach(n -> dvpool(ctx, n), ("GI_TIES", "GI_BIG1", "GI_BIG2", "GI_BIG3", "GI_BIG4"))
        end
        # The result is a copy: changing it does not change the pool.
        v = gdpool(ctx, "BODY399_RADII"); v[1] = -1.0
        @test gdpool(ctx, "BODY399_RADII")[1] == 6378.1366
        # Absent variable and wrong type.
        @test SPICE.gdpool("NO_SUCH_VARIABLE") === nothing
        @test throws_kind(() -> gdpool(ctx, "NO_SUCH_VARIABLE"), :kernelvarnotfound)
        @test throws_kind(() -> gcpool(ctx, "NO_SUCH_VARIABLE"), :kernelvarnotfound)
        @test throws_kind(() -> gipool(ctx, "NO_SUCH_VARIABLE"), :kernelvarnotfound)
        @test throws_kind(() -> gcpool(ctx, "BODY399_RADII"), :typemismatch)
        @test throws_kind(() -> gdpool(ctx, "TKFRAME_31009_SPEC"), :typemismatch)
    end

    @apart @testset "bodvrd" begin
        fails = Failures()
        for body in ("EARTH", "MOON", "SUN", "MARS", "JUPITER", "399", "earth")
            check!(fails, bodvrd(ctx, body, "RADII") == SPICE.bodvrd(body, "RADII"), "RADII", body)
        end
        @test bodvrd(ctx, 399, "RADII") == SPICE.bodvrd("399", "RADII")
        for (body, item) in (
                ("EARTH", "POLE_RA"), ("MOON", "NUT_PREC_PM"), ("MARS", "PM"),
                ("JUPITER", "NUT_PREC_RA"), ("SATURN", "POLE_DEC"),
            )
            check!(fails, bodvrd(ctx, body, item) == SPICE.bodvrd(body, item, 100), "bodvrd", body, item)
        end
        if have_gm
            for body in ("EARTH", "MOON", "SUN", "EARTH BARYCENTER", "JUPITER BARYCENTER", "3")
                check!(fails, bodvrd(ctx, body, "GM") == SPICE.bodvrd(body, "GM"), "GM", body)
            end
            @test bodvrd(ctx, "EARTH", "GM") isa Vector{Float64}
        end
        @test isempty(fails)
        @test_throws SPICE.SpiceError SPICE.bodvrd("EARTH", "NO_SUCH_ITEM")
        # A blank name is an error in CSPICE; Almanac reports "not found".
        @test bodn2c(ctx, "") === nothing
        @test namfrm(ctx, "") == 0
        @test throws_kind(() -> bodvrd(ctx, "EARTH", "NO_SUCH_ITEM"), :kernelvarnotfound)
        @test throws_kind(() -> bodvrd(ctx, "NO_SUCH_BODY", "RADII"), :notranslation)
    end

    @apart @testset "bodn2c, bods2c, bodc2n, bodc2s" begin
        names = [
            "EARTH", "earth", "  Earth  ", "MOON", "SUN", "SSB", "SOLAR SYSTEM BARYCENTER",
            "EMB", "EARTH BARYCENTER", "Earth-Moon   Barycenter", "MARS", "PHOBOS", "IO",
            "JUPITER BARYCENTER", "PLUTO", "CHARON", "CERES", "CASSINI", "VOYAGER 1", "MRO",
            "JUICE", "DSS-14", "TEMPEL 1", "NO_SUCH_BODY", "399", "-82",
        ]
        fails = Failures()
        for name in names
            check!(fails, bodn2c(ctx, name) == SPICE.bodn2c(name), "bodn2c", repr(name))
            check!(fails, bods2c(ctx, name) == SPICE.bods2c(name), "bods2c", repr(name))
        end
        codes = [
            0, 1, 3, 10, 199, 299, 301, 399, 401, 499, 501, 599, 699, 799, 899, 999, 901,
            2000001, -82, -31, -32, -74, 399014, 1000005, 123456789, -99999,
        ]
        for code in codes
            check!(fails, bodc2n(ctx, code) == SPICE.bodc2n(code), "bodc2n", code)
            check!(fails, bodc2s(ctx, code) == SPICE.bodc2s(code), "bodc2s", code)
        end
        @test isempty(fails)
        # scn2id and scid2n (`sctran` of SPICELIB, no C wrapper): the last "SCLK" among the
        # first 32 characters is removed.
        function c_scn2id(name)
            id, found = Ref{Cint}(0), Ref{Cint}(0)
            ccall((:scn2id_, SPICE.libcspice), Cint, (Ptr{UInt8}, Ref{Cint}, Ref{Cint}, Cint), name, id, found, ncodeunits(name))
            return found[] == 0 ? nothing : Int(id[])
        end
        function c_scid2n(id)
            buf, found = fill(UInt8(' '), 64), Ref{Cint}(0)
            ccall((:scid2n_, SPICE.libcspice), Cint, (Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Cint), id, buf, found, 64)
            return found[] == 0 ? nothing : String(rstrip(String(buf)))
        end
        clocks = [
            "CASSINI SCLK", "cassini sclk", "CASSINI", "  mgs   sclk ", "MARS GLOBAL SURVEYOR SCLK",
            "SCLK", "", "SCLK CASSINI", "CASSINISCLK", "CASSINI SCLK SCLK", "NO SUCH CLOCK SCLK",
            "VOYAGER 1" * " "^25 * "SCLK", "-82 SCLK", "MRO SCL",
        ]
        for name in clocks
            check!(fails, scn2id(ctx, name) == c_scn2id(name), "scn2id", repr(name))
        end
        for code in codes
            check!(fails, scid2n(ctx, code) == c_scid2n(code), "scid2n", code)
        end
        @test isempty(fails)
        @test bodc2s(ctx, typemin(Int32)) == SPICE.bodc2s(typemin(Int32)) == "-2147483648"
        # Round trip for a random sample of the built-in codes.
        fails = Failures()
        for code in rand(rng, 1:999, 200)
            name = SPICE.bodc2n(code)
            check!(fails, bodc2n(ctx, code) == name, "bodc2n", code)
            name === nothing || check!(fails, bodn2c(ctx, name) == SPICE.bodn2c(name), "bodn2c", code, repr(name))
        end
        @test isempty(fails)
    end

    @apart @testset "namfrm, frmnam" begin
        names = [
            "J2000", "j2000", "ECLIPJ2000", "B1950", "FK4", "GALACTIC", "MARSIAU",
            "ECLIPB1950", "DE-118", "IAU_EARTH", "IAU_MOON", "IAU_MARS", "IAU_JUPITER",
            "IAU_SUN", "IAU_PHOBOS", "ITRF93", "EARTH_FIXED", "MOON_PA", "MOON_ME",
            "MOON_PA_DE440", "MOON_ME_DE440_ME421", "MOON_PA_DE421", "MOON_ME_DE421",
            "NO_SUCH_FRAME",
        ]
        fails = Failures()
        for name in names
            check!(fails, namfrm(ctx, name) == SPICE.namfrm(name), "namfrm", name)
        end
        for id in [
                1:21..., 10001:10020..., 10081, 13000, 31000, 31001, 31006, 31007, 31008,
                31009, 0, -5, 123456,
            ]
            check!(fails, frmnam(ctx, id) == SPICE.frmnam(id), "frmnam", id)
        end
        @test isempty(fails)
    end

    @apart @testset "no leap-second kernel" begin
        empty = SpiceContext()
        @test_throws Almanac.SpiceError utc2et(empty, "2000-01-01T12:00:00")
    end
end

group("5b. body names that later pool entries take over, against CSPICE") do
    # Found by the replay of NAIF's tspice (F_BODCOD): a name given in the pool to a
    # second code takes the name away from the first code, as `zzbodker` of CSPICE does,
    # and a name of the pool takes it away from the code of `boddef` or of the built-in
    # table that had it (`zzbodc2n`).
    cases = [
        ["NAIF_BODY_CODE = ( 1011 1012 )", "NAIF_BODY_NAME = ( 'F_X' 'F_X' )"],
        ["NAIF_BODY_CODE = ( 1011 1012 1011 )", "NAIF_BODY_NAME = ( 'F_X' 'f_x' 'F_Y' )"],
        ["NAIF_BODY_CODE = ( 1011 1012 1013 )", "NAIF_BODY_NAME = ( 'A' 'B' 'A' )"],
        ["NAIF_BODY_CODE = ( 1012 )", "NAIF_BODY_NAME = ( 'EARTH' )"],
        # leading blanks of a name are kept (F_ZZBDTRN)
        ["NAIF_BODY_CODE = ( 1013 )", "NAIF_BODY_NAME = ( '   A' )"],
    ]
    fails = Failures()
    try
        for (k, cv) in enumerate(cases)
            SPICE.clpool()
            SPICE.lmpool(cv)
            ctx = SpiceContext()
            Almanac.lmpool(ctx, cv)
            for code in (399, 1011, 1012, 1013)
                expected = try
                    SPICE.bodc2n(code)
                catch
                    nothing
                end
                check!(fails, Almanac.bodc2n(ctx, code) == expected, "bodc2n", k, code)
            end
            for name in ("F_X", "F_Y", "A", "B", "EARTH")
                check!(fails, Almanac.bodn2c(ctx, name) == SPICE.bodn2c(name), "bodn2c", k, name)
            end
        end
        @test isempty(fails)
    finally
        SPICE.clpool()
    end
    # `boddef` cannot be undone in CSPICE, so these cases take the results that tspice
    # recorded from SPICELIB instead of the oracle.
    ctx = SpiceContext()
    Almanac.boddef(ctx, "F_M", 1013)
    Almanac.lmpool(ctx, ["NAIF_BODY_CODE = 1014", "NAIF_BODY_NAME = 'F_M'"])
    @test Almanac.bodc2n(ctx, 1013) === nothing && Almanac.bodn2c(ctx, "F_M") == 1014
    ctx = SpiceContext()
    Almanac.boddef(ctx, "EARTH", 1017)
    @test Almanac.bodc2n(ctx, 399) === nothing && Almanac.bodc2n(ctx, 1017) == "EARTH"
    Almanac.boddef(ctx, "EARTH", 399)
    @test Almanac.bodc2n(ctx, 399) == "EARTH"
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
