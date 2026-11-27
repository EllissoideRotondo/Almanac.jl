# Group 8: quality gates. Zero allocation and inferred return types for the integer-ID
# forms on an explicit context, after a warm-up call.
#
# Each measurement is inside a function, so that the arguments are not global variables
# (a call on untyped globals allocates on its own account).

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Almanac keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALMANAC_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

alloc_spkezr(ctx, et) = @allocated spkezr(ctx, 301, et, 1, "NONE", 399)
alloc_spkpos(ctx, et) = @allocated spkpos(ctx, 301, et, 1, "NONE", 399)
alloc_spkgeo(ctx, et) = @allocated spkgeo(ctx, 301, et, 1, 399)
alloc_spkgps(ctx, et) = @allocated spkgps(ctx, 301, et, 1, 399)
alloc_spkezr_chain(ctx, et) = @allocated spkezr(ctx, 301, et, 1, "NONE", 10)
alloc_spkezr_rot(ctx, et, ref) = @allocated spkezr(ctx, 301, et, ref, "NONE", 399)
alloc_spkder(ctx, et, ::Val{N}) where {N} = @allocated spkder(ctx, 301, et, 1, 399, Val(N))
alloc_spkder_rot(ctx, et, ref, ::Val{N}) where {N} = @allocated spkder(ctx, 301, et, ref, 399, Val(N))
alloc_pxform(ctx, from, to, et) = @allocated pxform(ctx, from, to, et)
alloc_sxform(ctx, from, to, et) = @allocated sxform(ctx, from, to, et)
alloc_pxder(ctx, from, to, et, ::Val{N}) where {N} = @allocated pxder(ctx, from, to, et, Val(N))
alloc_query(q, et) = @allocated q(et)
alloc_query_pos(q, et) = @allocated spkpos(q, et)
alloc_query_sx(q, et) = @allocated sxform(q, et)
alloc_names(ctx, et) = @allocated spkezr(ctx, "MOON", et, "J2000", "NONE", "EARTH")
alloc_pxform_names(ctx, et) = @allocated pxform(ctx, "J2000", "IAU_EARTH", et)
nalloc_mixed_names(ctx, et) = @allocations spkezr(ctx, "Moon", et, "J2000", "NONE", "Earth")
nalloc_normalize(s) = @allocations Almanac.normalize_name(s)
nalloc_str2et(ctx, s) = @allocations str2et(ctx, s)
nalloc_timout(ctx, et, pictur) = @allocations timout(ctx, et, pictur)

"Allocation of the second call of `f(args...)` (the first call compiles)."
function warm(f, args...)
    f(args...)
    return f(args...)
end

group("8. quality gates", Kernels.FULL_SET_FILES...) do
    ctx = neo_context(Kernels.full_set())
    et = 1.2345e8
    iau_earth = Int(namfrm(ctx, "IAU_EARTH"))
    moon_pa = Int(namfrm(ctx, "MOON_PA"))
    moon_me = Int(namfrm(ctx, "MOON_ME"))
    eclip = 17
    @test iau_earth == 10013 && moon_pa == 31000 && moon_me == 31001

    @apart @testset "zero allocation, integer-ID forms" begin
        @test warm(alloc_spkezr, ctx, et) == 0
        @test warm(alloc_spkpos, ctx, et) == 0
        @test warm(alloc_spkgeo, ctx, et) == 0
        @test warm(alloc_spkgps, ctx, et) == 0
        @test warm(alloc_spkezr_chain, ctx, et) == 0
        for ref in (eclip, iau_earth, moon_pa, moon_me)
            @test warm(alloc_spkezr_rot, ctx, et, ref) == 0
        end
        for N in 0:3
            @test warm(alloc_spkder, ctx, et, Val(N)) == 0
            @test warm(alloc_spkder_rot, ctx, et, iau_earth, Val(N)) == 0
            @test warm(alloc_spkder_rot, ctx, et, moon_pa, Val(N)) == 0
        end
        for (from, to) in (
                (1, eclip), (1, iau_earth), (iau_earth, 1), (1, moon_pa),
                (1, moon_me), (iau_earth, moon_pa),
            )
            @test warm(alloc_pxform, ctx, from, to, et) == 0
            @test warm(alloc_sxform, ctx, from, to, et) == 0
            for N in 0:3
                @test warm(alloc_pxder, ctx, from, to, et, Val(N)) == 0
            end
            q = rotation_query(ctx, from, to)
            @test warm(alloc_query, q, et) == 0
            @test warm(alloc_query_sx, q, et) == 0
        end
        for (targ, obs) in ((301, 399), (301, 10)), ref in (1, iau_earth, moon_me)
            q = state_query(ctx, targ, ref, "NONE", obs)
            @test warm(alloc_query, q, et) == 0
            @test warm(alloc_query_pos, q, et) == 0
        end
    end

    @apart @testset "allocation of the name forms (reported, not asserted)" begin
        measured(
            "  name form spkezr(ctx, \"MOON\", et, \"J2000\", \"NONE\", \"EARTH\"): " *
                "$(warm(alloc_names, ctx, et)) bytes"
        )
        measured(
            "  name form pxform(ctx, \"J2000\", \"IAU_EARTH\", et): " *
                "$(warm(alloc_pxform_names, ctx, et)) bytes"
        )
        @test true
    end

    @apart @testset "names in mixed case: the ljucrs rule, no allocation" begin
        nn = Almanac.normalize_name
        for (s, want) in (
                ("Earth", "EARTH"), (" solar  System   barycenter ", "SOLAR SYSTEM BARYCENTER"),
                ("j2000", "J2000"), ("  ", ""), ("", ""), ("MOON", "MOON"), ("Moon ", "MOON"),
            )
            @test nn(s) == want
            @test nn(nn(s)) == want
        end
        @test nn("Étoile  x") == "ÉTOILE X"
        @test bodn2c(ctx, " earth ") == bodn2c(ctx, "Earth") == 399
        @test warm(nalloc_normalize, "Earth") == 1
        # The lookup reads the name as given (it cost 14 allocations per name before).
        @test warm(nalloc_mixed_names, ctx, et) == 0
        # Names defined in the pool go through the same table.
        pctx = SpiceContext()
        lmpool(
            pctx, [
                "NAIF_BODY_NAME = ( ' My  Probe ', 'Étoile' )", "NAIF_BODY_CODE = ( -77, -78 )",
                "FRAME_MY_FRAME = -77001", "FRAME_-77001_NAME = 'My  Frame'",
                "FRAME_-77001_CLASS = 4", "FRAME_-77001_CLASS_ID = -77001",
                "FRAME_-77001_CENTER = 'my  probe'",
            ]
        )
        # A body name in a frame definition is resolved while the snapshot is built.
        @test frinfo(pctx, -77001) == (-77, 4, -77001)
        # The table of the snapshot holds the pairs of the pool, and no other.
        @test Dict(Almanac.kernelset(pctx).name2code) == Dict("MY PROBE" => -77, "ÉTOILE" => -78)
        @test get(Almanac.kernelset(pctx).name2code, :MY_PROBE, 0) == 0
        @test !haskey(Almanac.kernelset(pctx).name2code, :MY_PROBE)
        @test bodn2c(pctx, "my probe") == bodn2c(pctx, "  MY PROBE") == -77
        @test bodn2c(pctx, "MY  PROBEX") === nothing && bodn2c(pctx, "MY") === nothing
        @test bodn2c(pctx, "étoile") == bodn2c(pctx, " ÉTOILE") == -78   # not ASCII: `uppercase`
        @test namfrm(pctx, " my frame ") == -77001 && namfrm(pctx, "MYFRAME") == 0
        @test Almanac.name_lookup(Almanac.NameTable(), "EARTH") === nothing
    end

    @apart @testset "time strings: str2et allocates nothing, timout only its result" begin
        # The parser state (64 tokens, a picture of 320 bytes) and the line of `timout`
        # stay on the stack.
        for s in (
                "2020-01-01T00:00:00", "2020 JAN 01 12:00:00.000 TDB", "JD 2451545.0",
                "March 1, 2025 4:00 AM PST", "2025-03-01 12:00:00 UTC+5:30", "2025-060 // 12:00",
            )
            @test warm(nalloc_str2et, ctx, s) == 0
        end
        for pictur in (
                "YYYY-MON-DD HR:MN:SC.### ::TDB", "Mon DD, YYYY AP:MN ampm ::UTC-8 ::RND",
                "JULIAND.##### ::TDB", "Weekday ?ERA? YYYY-DOY HR:MN:SC.###### ::MCAL",
            )
            @test warm(nalloc_timout, ctx, et, pictur) == 1
        end
    end

    @apart @testset "inferred return types" begin
        @test (@inferred spkezr(ctx, 301, et, 1, "NONE", 399)) isa Tuple{SVector{6, Float64}, Float64}
        @test (@inferred spkpos(ctx, 301, et, 1, "NONE", 399)) isa Tuple{SVector{3, Float64}, Float64}
        @test (@inferred spkgeo(ctx, 301, et, 1, 399)) isa Tuple{SVector{6, Float64}, Float64}
        @test (@inferred spkezr(ctx, 301, et, iau_earth, "NONE", 399)) isa Tuple{SVector{6, Float64}, Float64}
        @test (@inferred spkezr(ctx, "MOON", et, "J2000", "NONE", "EARTH")) isa Tuple{SVector{6, Float64}, Float64}
        @test (@inferred spkder(ctx, 301, et, 1, 399, Val(0))) isa NTuple{1, SVector{3, Float64}}
        @test (@inferred spkder(ctx, 301, et, 1, 399, Val(1))) isa NTuple{2, SVector{3, Float64}}
        @test (@inferred spkder(ctx, 301, et, 1, 399, Val(2))) isa NTuple{3, SVector{3, Float64}}
        @test (@inferred spkder(ctx, 301, et, 1, 399, Val(3))) isa NTuple{4, SVector{3, Float64}}
        @test (@inferred spkder(ctx, 301, et, moon_pa, 399, Val(3))) isa NTuple{4, SVector{3, Float64}}
        @test (@inferred pxform(ctx, 1, iau_earth, et)) isa SMatrix{3, 3, Float64, 9}
        @test (@inferred pxform(ctx, 1, moon_me, et)) isa SMatrix{3, 3, Float64, 9}
        @test (@inferred pxform(ctx, "J2000", "IAU_EARTH", et)) isa SMatrix{3, 3, Float64, 9}
        @test (@inferred sxform(ctx, 1, iau_earth, et)) isa SMatrix{6, 6, Float64, 36}
        @test (@inferred pxder(ctx, 1, iau_earth, et, Val(0))) isa NTuple{1, SMatrix{3, 3, Float64, 9}}
        @test (@inferred pxder(ctx, 1, iau_earth, et, Val(1))) isa NTuple{2, SMatrix{3, 3, Float64, 9}}
        @test (@inferred pxder(ctx, 1, iau_earth, et, Val(2))) isa NTuple{3, SMatrix{3, 3, Float64, 9}}
        @test (@inferred pxder(ctx, 1, moon_pa, et, Val(3))) isa NTuple{4, SMatrix{3, 3, Float64, 9}}
        # An integer epoch and a Float32 epoch give Float64 results.
        @test (@inferred spkezr(ctx, 301, 0, 1, "NONE", 399)) isa Tuple{SVector{6, Float64}, Float64}
    end

    @apart @testset "package image and method table" begin
        # A method defined twice makes precompilation fail, and the package then loads
        # without its image (slow first calls) although every other test passes.
        project = dirname(Base.active_project())
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$project -e "using Almanac"`
        err = IOBuffer()
        ok = success(pipeline(cmd; stderr = err, stdout = devnull))
        text = String(take!(err))
        @test ok
        @test !occursin("overwritten", text)
        @test !occursin("precompilation", text)
        @test isempty(Test.detect_ambiguities(Almanac))
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
