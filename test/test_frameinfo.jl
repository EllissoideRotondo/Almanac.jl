# Group 19: frame information, the lists of frames, body constants, body-fixed rotations,
# binary PCK coverage, local solar time and L_s against CSPICE.
#
# Conventions for an absent item: `frinfo`, `cidfrm`, `cnmfrm` and `ccifrm` return
# `nothing` where CSPICE returns `found = false`; `tkfram` throws where `tkfram_c` signals
# an error. `tkfram_` has no wrapper in SPICE.jl and is called through `ccall`.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Almanac keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALMANAC_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

isdefined(@__MODULE__, :CKFix) || include("ckfixtures.jl")
isdefined(@__MODULE__, :DynFix) || include("dynfixtures.jl")
import ForwardDiff

"`tkfram_` of CSPICE: `(rot, relative frame ID)`, the error kind, or `nothing` if not found."
function oracle_tkfram(id::Integer)
    rot = zeros(9)
    frame = Ref{Cint}(0)
    found = Ref{Cint}(0)
    ccall(
        (:tkfram_, SPICE.libcspice), Cint, (Ref{Cint}, Ptr{Float64}, Ref{Cint}, Ref{Cint}),
        Cint(id), rot, frame, found
    )
    kind = cspice_kind()
    kind === :none || return kind
    found[] == 0 && return nothing
    return SMatrix{3, 3, Float64, 9}(rot), Int(frame[])
end

"`tkfram(ctx, id)`, or the kind of the `SpiceError` that it throws."
neo_tkfram(ctx, id) = try
    tkfram(ctx, id)
catch e
    e isa Almanac.SpiceError ? e.kind : rethrow()
end

group("19. frame information, body constants, local solar time", Kernels.FULL_SET_FILES...) do
    rng = MersenneTwister(SEED + 19)
    set = [Kernels.full_set()..., CKFix.frame_kernel(), DynFix.frame_kernel()]
    extra = joinpath(mktempdir(), "objects.tf")
    write(
        extra, """
        KPL/FK

        \\begindata
           OBJECT_MARS_FRAME      = 'IAU_MARS'
           OBJECT_599_FRAME       = 'J2000'
           OBJECT_EARTH_FRAME     = 13000
           OBJECT_-77_FRAME       = 'TESTSC_T3'
           OBJECT_301_FRAME       = 'MOON_ME'
           OBJECT_299_FRAME       = 31000
           NAIF_BODY_NAME        += ( 'TEST_PROBE' )
           NAIF_BODY_CODE        += ( -77 )
        \\begintext
        """
    )
    # The second pass adds the `OBJECT_*_FRAME` variables. Only `cidfrm` and `cnmfrm` read
    # them, so the other calls are compared in the first pass only.
    for with_objects in (false, true)
        with_objects && push!(set, extra)
        oracle_load(set)
        ctx = neo_context(set)
        label = with_objects ? " (with OBJECT_*_FRAME variables)" : ""

        frames = vcat(
            collect(-5:25), collect(10001:10150), collect(13000:13010),
            [1850000 + i for i in 0:25], [
                31000, 31001, 31002, 31003, 31004, 31005, 31006,
                31007, 31008, -77001, -77003, -77900, -77099, -99001, 1400499, 99999999,
            ]
        )
        with_objects || @testset "frinfo" begin
            nfound = 0
            fails = Failures()
            for id in frames
                want = SPICE.frinfo(id)
                got = frinfo(ctx, id)
                got == want || @info "frinfo" id got want
                check!(fails, got == want, "frinfo", id)
                want === nothing || (nfound += 1)
            end
            @test isempty(fails)
            println("  frinfo: ", length(frames), " frame IDs, ", nfound, " found")
            @test nfound > 150
            @test (@inferred Union{Nothing, NTuple{3, Int}} frinfo(ctx, 13000)) == SPICE.frinfo(13000)
        end

        @testset "cidfrm and cnmfrm$label" begin
            bodies = vcat(
                collect(-3:12), collect(199:100:999), collect(301:1:305),
                collect(401:402), collect(501:516), collect(601:618), [
                    10, 9511010, 2000001,
                    2000004, 2000433, 1000005, -77, -78, -31, 3, 4, 399, 301, 123456789,
                ]
            )
            nfound = 0
            fails = Failures()
            for b in bodies
                want = SPICE.cidfrm(b)
                got = cidfrm(ctx, b)
                got == want || @info "cidfrm" b got want
                check!(fails, got == want, "cidfrm", b)
                want === nothing || (nfound += 1)
            end
            println("  cidfrm: ", length(bodies), " bodies, ", nfound, " found")
            @test nfound > 50
            names = [
                "EARTH", "earth", "  Moon ", "MARS", "JUPITER", "SUN", "PHOBOS", "TEST_PROBE", "VENUS",
                "SOLAR SYSTEM BARYCENTER", "MARS BARYCENTER", "NOT A BODY", "TITAN", "CERES",
                "PLUTO", "IO", "Earth  Barycenter", "399", "-77", "499",
            ]
            for n in names
                want = SPICE.cnmfrm(n)
                got = cnmfrm(ctx, n)
                got == want || @info "cnmfrm" n got want
                check!(fails, got == want, "cnmfrm", n)
            end
            @test isempty(fails)
            # Regression: the frame name of `OBJECT_<body>_FRAME` is compared exactly with
            # the built-in names (`zzhscchk_`), so `iau_saturn` and ` IAU_SATURN` are not
            # the frame IAU_SATURN.
            if !with_objects
                try
                    for v in ("iau_saturn", " IAU_SATURN", "IAU_SATURN")
                        SPICE.pcpool("OBJECT_SATURN_FRAME", [v])
                        pcpool(ctx, "OBJECT_SATURN_FRAME", [v])
                        @test cidfrm(ctx, 699) == SPICE.cidfrm(699)
                        @test cnmfrm(ctx, "SATURN") == SPICE.cnmfrm("SATURN")
                    end
                    pcpool(ctx, "OBJECT_SATURN_FRAME", ["iau_saturn"])
                    @test cidfrm(ctx, 699) === nothing
                finally
                    SPICE.dvpool("OBJECT_SATURN_FRAME")
                    dvpool(ctx, "OBJECT_SATURN_FRAME")
                end
            end
        end

        with_objects || @testset "ccifrm" begin
            n = 0
            fails = Failures()
            for id in frames
                info = SPICE.frinfo(id)
                info === nothing && continue
                want = SPICE.ccifrm(info[2], info[3])
                got = ccifrm(ctx, info[2], info[3])
                got == want || @info "ccifrm" id info got want
                check!(fails, got == want, "ccifrm", id)
                n += 1
            end
            for (c, k) in (
                    (1, 0), (1, 22), (2, 398), (2, 31006), (3, -77001), (3, 5), (4, -77900),
                    (5, 1850001), (6, 1850021), (7, 1), (0, 0), (2, 10), (4, 31001),
                )
                want = SPICE.ccifrm(c, k)
                got = ccifrm(ctx, c, k)
                got == want || @info "ccifrm" c k got want
                check!(fails, got == want, "ccifrm by class", c, k)
            end
            @test isempty(fails)
            println("  ccifrm: ", n, " (class, class ID) pairs of loaded frames")
        end

        with_objects || @testset "tkfram" begin
            worst = 0.0
            nfound = 0
            for id in (-77900, 31000, 31001, 31002, 31003, 31006, 31007, 1850001, 10013, 5, -5)
                want = oracle_tkfram(id)
                got = neo_tkfram(ctx, id)
                want isa Tuple || (@test got === want; continue)
                nfound += 1
                @test got[2] == want[2]
                worst = max(worst, maximum(abs.(got[1] .- want[1])))
            end
            @test nfound >= 3
            @test worst <= 4 * eps()
            measured(
                @sprintf(
                    "  %-58s n=%-7d max=%.3e", "tkfram: element of the matrix",
                    nfound, worst
                )
            )
        end
    end

    @testset "tkfram and pxform: a TK frame definition that is not valid" begin
        # Regression: tkfram returned `nothing`, and pxform threw :framedatanotfound, for
        # an incomplete or inconsistent TK frame definition; tkfram_c and pxform_c signal
        # the error of tkfram_ for each. One case for each check of tkfram_.
        head = [
            "FRAME_TKA = -1000", "FRAME_-1000_NAME = 'TKA'", "FRAME_-1000_CLASS = 4",
            "FRAME_-1000_CLASS_ID = -1000", "FRAME_-1000_CENTER = 399",
        ]
        rel(name) = "TKFRAME_-1000_RELATIVE = '$name'"
        spec(s) = "TKFRAME_-1000_SPEC = '$s'"
        angles(units, axes) = [
            rel("J2000"); spec("ANGLES"); "TKFRAME_-1000_ANGLES = (1, 2, 3)";
            "TKFRAME_-1000_AXES = $axes"; "TKFRAME_-1000_UNITS = '$units'"
        ]
        unit = "TKFRAME_-1000_MATRIX = (1, 0, 0, 0, 1, 0, 0, 0, 1)"
        cases = (
            (-1000, head),                                                  # variablenotfound
            (-1000, [head; rel("NO_SUCH_FRAME"); spec("MATRIX"); unit]),    # badframespec
            (-1000, [head; rel("TKA"); spec("MATRIX"); unit]),              # badframespec2
            (-1000, [head; rel("J2000"); "TKFRAME_TKA_RELATIVE = 'J2000'"; spec("MATRIX"); unit]),
            (-1000, [head; rel("J2000"); spec("EULER")]),                   # unknownframespec
            (-1000, [head; rel("J2000"); spec(" MATRIX"); unit]),           # a leading blank counts
            (-1000, [head; rel("J2000"); spec("MATRIX"); "TKFRAME_-1000_MATRIX = (1, 0, 0)"]),
            (-1000, [head; angles("FURLONGS", "(3, 1, 3)")]),               # unitsnotrec
            (-1000, [head; angles("DEGREES", "(3, 4, 3)")]),                # badaxisnumbers
            (
                -2001, [
                    "FRAME_TKB = -1001", "FRAME_-1001_NAME = 'TKB'", "FRAME_-1001_CLASS = 4",
                    "FRAME_-1001_CLASS_ID = -2001", "FRAME_-1001_CENTER = 399",
                    "TKFRAME_-2001_RELATIVE = 'J2000'", "TKFRAME_-2001_SPEC = 'QUATERNION'",
                    "TKFRAME_-2001_Q = (1, 0, 0, 0)",
                ],
            ),                         # incompleteframe
        )
        r = zeros(9)
        try
            for (id, lines) in cases
                SPICE.kclear()
                SPICE.lmpool(lines)
                tctx = SpiceContext()
                lmpool(tctx, lines)
                want = oracle_tkfram(id)
                @test want isa Symbol && neo_tkfram(tctx, id) === want
                name = id == -1000 ? "TKA" : "TKB"
                ccall(
                    (:pxform_c, SPICE.libcspice), Cvoid, (
                        Cstring, Cstring, Cdouble,
                        Ptr{Float64},
                    ), name, "J2000", 0.0, r
                )
                @test throws_kind(() -> pxform(tctx, name, "J2000", 0.0), cspice_kind())
            end
            @test neo_tkfram(SpiceContext(), 0) === oracle_tkfram(0) === :zeroframeid
        finally
            oracle_load(set)
        end
    end

    set = [Kernels.full_set()..., CKFix.frame_kernel(), DynFix.frame_kernel()]
    oracle_load(set)
    ctx = neo_context(set)
    ks = Almanac.kernelset(ctx)

    # A definition without its center: KERNELVARNOTFOUND in CSPICE (tspice F_CCIFRM).
    @apart @testset "ccifrm: a frame without a center" begin
        lines = [
            "FRAME_-999_CLASS = 3", "FRAME_-999_CLASS_ID = -1",
            "FRAME_-999_NAME = 'BADFRAME'", "FRAME_BADFRAME = -999",
        ]
        c2 = SpiceContext()
        lmpool(c2, lines)
        @test throws_kind(() -> ccifrm(c2, 3, -1), :kernelvarnotfound)
        lmpool(c2, ["FRAME_-999_CENTER = 399"])
        @test ccifrm(c2, 3, -1) == (-999, "BADFRAME", 399)
        try
            SPICE.lmpool(lines)
            ccall(
                (:ccifrm_c, SPICE.libcspice), Cvoid, (
                    Cint, Cint, Cint, Ref{Cint}, Ptr{UInt8},
                    Ref{Cint}, Ref{Cint},
                ), 3, -1, 33, Ref{Cint}(0), zeros(UInt8, 33), Ref{Cint}(0),
                Ref{Cint}(0)
            )
            @test cspice_kind() === :kernelvarnotfound
        finally
            foreach(
                SPICE.dvpool, [
                    "FRAME_-999_CLASS", "FRAME_-999_CLASS_ID", "FRAME_-999_NAME",
                    "FRAME_BADFRAME",
                ]
            )
        end
    end

    @apart @testset "ccifrm: a frame class that is not a number" begin
        c2 = SpiceContext()
        lmpool(c2, ["FRAME_-999_CLASS = 'X'"])
        @test throws_kind(() -> ccifrm(c2, 3, -1), :invalidframedef)
        try
            SPICE.lmpool(["FRAME_-999_CLASS = 'X'"])
            @test_throws SPICE.SpiceError SPICE.ccifrm(3, -1)
        finally
            SPICE.dvpool("FRAME_-999_CLASS")
        end
    end

    # Regression for a defect of CSPICE N0067: `ccifrm_` finds the name variable by
    # replacing the first `_CLASS` of the key `FRAME_<name>_CLASS`, so for the frame
    # A_CLASS_B it reads FRAME_A_NAME_B_CLASS and signals INVALIDFRAMEDEF. Almanac
    # replaces the trailing `_CLASS`.
    @apart @testset "ccifrm: a frame name that contains _CLASS" begin
        lines = [
            "FRAME_A_CLASS_B = -1500001", "FRAME_-1500001_NAME = 'A_CLASS_B'",
            "FRAME_A_CLASS_B_NAME = 'A_CLASS_B'", "FRAME_A_CLASS_B_CLASS = 4",
            "FRAME_A_CLASS_B_CLASS_ID = -1500001", "FRAME_-1500001_CENTER = 399",
        ]
        c2 = SpiceContext()
        lmpool(c2, lines)
        @test ccifrm(c2, 4, -1500001) == (-1500001, "A_CLASS_B", 399)
        try
            SPICE.lmpool(lines)
            ccall(
                (:ccifrm_c, SPICE.libcspice), Cvoid, (
                    Cint, Cint, Cint, Ref{Cint}, Ptr{UInt8},
                    Ref{Cint}, Ref{Cint},
                ), 4, -1500001, 33, Ref{Cint}(0), zeros(UInt8, 33), Ref{Cint}(0),
                Ref{Cint}(0)
            )
            @test cspice_kind() === :invalidframedef
        finally
            foreach(l -> SPICE.dvpool(first(split(l, " = "))), lines)
        end
    end

    # Regression: an integer of the pool is read by `gipool_`, which rounds a tie away
    # from zero (`i_dnnt`); the frame items and the frame of a body were rounded to
    # even (class 2, center 398, frame 10014 here).
    @apart @testset "integer pool values with a fraction of one half" begin
        lines = [
            "FRAME_HALFTEST = -1977130.5", "FRAME_-1977131_NAME = 'HALFTEST'",
            "FRAME_-1977131_CLASS = 2.5", "FRAME_-1977131_CLASS_ID = -1977130.5",
            "FRAME_-1977131_CENTER = 398.5", "OBJECT_-1977132_FRAME = 10014.5",
        ]
        c2 = SpiceContext()
        lmpool(c2, lines)
        try
            SPICE.lmpool(lines)
            @test namfrm(c2, "HALFTEST") == SPICE.namfrm("HALFTEST") == -1977131
            @test frinfo(c2, -1977131) == SPICE.frinfo(-1977131) == (399, 3, -1977131)
            @test cidfrm(c2, -1977132) == SPICE.cidfrm(-1977132)
            @test cidfrm(c2, -1977132)[1] == 10015
        finally
            foreach(l -> SPICE.dvpool(first(split(l, " = "))), lines)
        end
    end

    @apart @testset "tipbod and tisbod" begin
        dp = Float64[]; ds = Float64[]; dd = Float64[]
        structure_ok = true
        # 31008 is the frame class of the lunar binary PCK; the others use text PCK constants.
        for body in (399, 301, 499, 599, 10, 199, 299, 699, 799, 899, 401, 502, 606, 31008),
                ref in ("J2000", "ECLIPJ2000", "B1950", "FK4", "GALACTIC", 1, 17)
            for et in random_epochs(rng, 6, -6.0e8, 6.0e8)
                refname = ref isa Integer ? frmnam(ctx, ref) : ref
                a = tipbod(ctx, ref, body, et)
                b = SPICE.tipbod(refname, body, et)
                push!(dp, maximum(abs.(a .- b)))
                s = tisbod(ctx, ref, body, et)
                t = SPICE.tisbod(refname, body, et)
                push!(ds, maximum(abs.(s[1:3, 1:3] .- t[1:3, 1:3])))
                scale = maximum(abs.(t[4:6, 1:3]))
                push!(dd, maximum(abs.(s[4:6, 1:3] .- t[4:6, 1:3])) / scale)
                structure_ok &= s[1:3, 4:6] == zeros(3, 3) && s[4:6, 4:6] == s[1:3, 1:3]
            end
        end
        sp = report_quantiles("tipbod: element of the matrix", dp)
        ss = report_quantiles("tisbod: element of the rotation block", ds)
        sd = report_quantiles("tisbod: derivative block (relative)", dd)
        @test sp.max <= 1.0e-13
        @test ss.max <= 1.0e-13
        @test sd.max <= 1.0e-12
        @test structure_ok          # [R 0; dR R]
        @test (@inferred tipbod(ks, 1, 399, 1.0e8)) isa SMatrix{3, 3, Float64, 9}
        @test (@inferred tisbod(ks, 1, 399, 1.0e8)) isa SMatrix{6, 6, Float64, 36}
        @test throws_kind(() -> tipbod(ctx, "IAU_EARTH", 399, 0.0), :irfnotrec)
        @test throws_kind(() -> tipbod(ctx, "J2000", 123456, 0.0), :framedatanotfound)
        # A rotation rate without a pole: CSPICE names the absent variable.
        c2 = SpiceContext()
        pdpool(c2, "BODY123456_PM", [10.0, 1.0])
        pdpool(c2, "BODY123456_POLE_DEC", [80.0])
        @test throws_kind(() -> tisbod(c2, "J2000", 123456, 0.0), :kernelvarnotfound)
    end

    # Regression for a defect of CSPICE N0067: for cubic angles
    # (`BODY<sys>_MAX_PHASE_DEGREE = 3`), `tisbod_` divides the rate of the cubic term by
    # the square of the seconds per century, so its derivative block misses that term.
    @apart @testset "tisbod: the derivative block of cubic nutation angles" begin
        # The Moon of pck00011.tpc as body 1301, which is its own system, with the 13
        # angles of the system made cubic: no quadratic term, 100 degrees per century
        # cubed.
        angles = [
            125.045, -1935.5364525, 250.089, -3871.072905, 260.008, 475263.3328725,
            176.625, 487269.629985, 357.529, 35999.0509575, 311.589, 964468.49931,
            134.963, 477198.869325, 276.617, 12006.300765, 34.226, 63863.5132425,
            15.134, -5806.6093575, 119.743, 131.84064, 239.961, 6003.1503825,
            25.053, 473327.79642,
        ]
        pool = [
            "BODY1301_POLE_RA" => [269.9949, 0.0031, 0.0],
            "BODY1301_POLE_DEC" => [66.5392, 0.013, 0.0],
            "BODY1301_PM" => [38.3213, 13.17635815, -1.4e-12],
            "BODY1301_NUT_PREC_RA" => [
                -3.8787, -0.1204, 0.07, -0.0172, 0.0, 0.0072, 0.0, 0.0, 0.0, -0.0052, 0.0, 0.0, 0.0043,
            ],
            "BODY1301_NUT_PREC_DEC" => [
                1.5419, 0.0239, -0.0278, 0.0068, 0.0, -0.0029, 0.0009, 0.0, 0.0, 0.0008, 0.0, 0.0, -0.0009,
            ],
            "BODY1301_NUT_PREC_PM" => [
                3.561, 0.1208, -0.0642, 0.0158, 0.0252, -0.0066, -0.0047, -0.0046, 0.0028, 0.0052,
                0.004, 0.0019, -0.0044,
            ],
            "BODY1301_NUT_PREC_ANGLES" => reduce(vcat, [[angles[2i - 1], angles[2i], 0.0, 100.0] for i in 1:13]),
            "BODY1301_MAX_PHASE_DEGREE" => [3.0],
        ]
        c2 = SpiceContext()
        foreach(kv -> pdpool(c2, kv[1], kv[2]), pool)
        # Five-point central difference of tipbod. At this epoch the prime meridian angle
        # is about 2700 rad, so each element of the matrix carries a rounding of about
        # 6e-13, and the difference with h = 1000 s about 1e-15; its truncation error is
        # below 1e-17. Measured: Almanac 6.8e-17, CSPICE 3.7e-12.
        et, h = 1.0e9, 1000.0
        r(t) = tipbod(c2, "J2000", 1301, t)
        fd = (8 * (r(et + h) - r(et - h)) - (r(et + 2h) - r(et - 2h))) / (12h)
        @test maximum(abs.(tisbod(c2, "J2000", 1301, et)[4:6, 1:3] - fd)) <= 1.0e-15
        try
            foreach(kv -> SPICE.pdpool(kv[1], kv[2]), pool)
            @test maximum(abs.(SPICE.tisbod("J2000", 1301, et)[4:6, 1:3] - fd)) >= 1.0e-12
        finally
            foreach(kv -> SPICE.dvpool(kv[1]), pool)
        end
    end

    # Regression: Almanac took the first of two forms of the epoch or frame variable, and
    # checked the degree of the phase angles only for a body with terms that use them;
    # `tisbod_` signals COMPETINGEPOCHSPEC, COMPETINGFRAMESPEC and DEGREEOUTOFRANGE.
    @apart @testset "tisbod: competing specifications and the phase degree, as CSPICE" begin
        base = [
            "BODY1399_POLE_RA" => [10.0, 0.1, 0.0], "BODY1399_POLE_DEC" => [80.0, 0.0, 0.0],
            "BODY1399_PM" => [5.0, 300.0, 0.0],
        ]
        cases = [
            :competingepochspec => ["BODY1399_CONSTANTS_JED_EPOCH" => [2451545.0], "BODY1399_CONSTS_JED_EPOCH" => [2451545.0]],
            :competingframespec => ["BODY1399_CONSTANTS_REF_FRAME" => [1.0], "BODY1399_CONSTS_REF_FRAME" => [1.0]],
            :degreeoutofrange => ["BODY1399_NUT_PREC_ANGLES" => [1.0, 2.0, 3.0, 4.0, 5.0], "BODY1399_MAX_PHASE_DEGREE" => [4.0]],
            :none => ["BODY1399_NUT_PREC_ANGLES" => [1.0, 2.0, 3.0], "BODY1399_MAX_PHASE_DEGREE" => [2.0]],
        ]
        for (want, extra) in cases
            pool = vcat(base, extra)
            c2 = SpiceContext()
            foreach(kv -> pdpool(c2, kv[1], kv[2]), pool)
            try
                foreach(kv -> SPICE.pdpool(kv[1], kv[2]), pool)
                ccall(
                    (:tisbod_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Cdouble, Ptr{Float64}),
                    "J2000", 1399, 0.0, zeros(36)
                )
                @test cspice_kind() === want
            finally
                foreach(kv -> SPICE.dvpool(kv[1]), pool)
            end
            for f in (tisbod, tipbod)
                @test throws_kind(() -> f(c2, "J2000", 1399, 0.0), want) || want === :none
            end
            want === :none && @test tisbod(c2, "J2000", 1399, 0.0) isa AbstractMatrix
        end
    end

    # CSPICE N0067 (`framex`) signals HASHISFULL past 5209 frames defined in kernels;
    # Almanac has no limit. CSPICE is not called.
    @apart @testset "more frames than CSPICE can hold" begin
        n = 5300
        lines = String[]
        for i in 1:n
            id = 1_400_000 + i
            append!(
                lines, [
                    "FRAME_TKF$i = $id", "FRAME_$(id)_NAME = 'TKF$i'", "FRAME_$(id)_CLASS = 4",
                    "FRAME_$(id)_CLASS_ID = $id", "FRAME_$(id)_CENTER = 399",
                ]
            )
        end
        c2 = SpiceContext()
        lmpool(c2, lines)
        for i in (1, 5210, n)
            @test namfrm(c2, "TKF$i") == 1_400_000 + i
            @test frmnam(c2, 1_400_000 + i) == "TKF$i"
        end
        @test frinfo(c2, 1_400_000 + n) == (399, 4, 1_400_000 + n)
    end

    @apart @testset "pckfrm and pckcov" begin
        want = sort!(Int[i for i in SPICE.pckfrm!(SPICE.SpiceIntCell(100), Kernels.MOON_BPC)])
        @test pckfrm(Kernels.MOON_BPC) == want
        @test pckfrm!([7, 3], Kernels.MOON_BPC) == sort!(unique!(vcat([7, 3], want)))
        # Regression: `spkobj!` did not sort its input, so that an input that was not
        # sorted came back with duplicates; `spkobj_c` inserts into a set.
        @test spkobj!([399, 3, 399], Kernels.DE440S) == sort!(unique!(vcat([399, 3], spkobj(Kernels.DE440S))))
        @test !isempty(want)
        for id in vcat(want, [399, 0])
            cell = SPICE.pckcov!(SPICE.SpiceDoubleCell(2000), Kernels.MOON_BPC, id)
            w = Float64[x for x in cell]
            got = pckcov(Kernels.MOON_BPC, id)
            @test got isa Window
            @test collect(got) == w
        end
        pre = Window([-1.0e12, -0.9e12])
        got = pckcov!(pre, Kernels.MOON_BPC, want[1])
        cell = SPICE.SpiceDoubleCell(2000)
        SPICE.wninsd!(cell, -1.0e12, -0.9e12)
        SPICE.pckcov!(cell, Kernels.MOON_BPC, want[1])
        @test collect(got) == Float64[x for x in cell]
        @test throws_kind(() -> pckfrm(Kernels.DE440S), :invalidfiletype)
        @test throws_kind(() -> pckcov(Kernels.LSK, 301), :invalidarchtype)
    end

    @apart @testset "bodfnd and bodvcd" begin
        n = 0
        fails = Failures()
        for body in (399, 301, 499, 10, 599, 1, 3, 9, -77, 2000004, 606),
                item in (
                    "RADII", "PM", "POLE_RA", "POLE_DEC", "GM", "NUT_PREC_RA", "LONG_AXIS",
                    "radii", "NOTHING", "NUT_PREC_ANGLES", "MAX_PHASE_DEGREE",
                )
            want = SPICE.bodfnd(body, item)
            check!(fails, bodfnd(ctx, body, item) == want, "bodfnd", body, item)
            want || continue
            check!(fails, bodvcd(ctx, body, item) == SPICE.bodvcd(body, item), "bodvcd", body, item)
            n += 1
        end
        @test isempty(fails)
        println("  bodvcd: ", n, " items equal to CSPICE")
        @test n > 30
        @test throws_kind(() -> bodvcd(ctx, 399, "NOTHING"), :kernelvarnotfound)
        @test (@inferred bodfnd(ctx, 399, "RADII")) === true
    end

    @apart @testset "et2lst" begin
        n = 0
        dsec = Int[]
        # Bodies with a segment in de440s.bsp; Venus rotates in the retrograde direction.
        for body in (399, 301, 299, 199, 10),
                typ in ("PLANETOCENTRIC", "PLANETOGRAPHIC", "planetographic", " Planetocentric ")
            for et in random_epochs(rng, 12, -6.0e8, 6.0e8)
                lon = (2 * rand(rng) - 1) * 2pi
                want = SPICE.et2lst(et, body, lon, strip(typ))
                got = et2lst(ctx, et, body, lon, typ)
                n += 1
                if got == want
                    push!(dsec, 0)
                    continue
                end
                # The second is truncated: one unit of rounding in the angle can move it.
                a = 3600 * got[1] + 60 * got[2] + got[3]
                b = 3600 * want[1] + 60 * want[2] + want[3]
                d = min(mod(a - b, 86400), mod(b - a, 86400))
                push!(dsec, d)
                d <= 1 || @info "et2lst" body typ et lon got want
            end
        end
        nexact = count(==(0), dsec)
        measured(
            @sprintf(
                "  %-58s n=%-7d equal=%d max=%d s", "et2lst: all five outputs", n,
                nexact, maximum(dsec)
            )
        )
        @test maximum(dsec) <= 1
        @test nexact >= n - 3
        @test et2lst(ctx, 0.0, 10, 1.0, "PLANETOCENTRIC") == SPICE.et2lst(0.0, 10, 1.0, "PLANETOCENTRIC")
        @test throws_kind(() -> et2lst(ctx, 0.0, 399, 0.0, "GEODETIC"), :unknownsystem)
        # `BODY<id>_PGR_POSITIVE_LON` sets the sense of a planetographic longitude.
        # Regression: `pgrrec_` compares the first four characters without blanks, in
        # upper case (`ljucrs_`), so "EASTERN" is east.
        try
            for sense in ("EAST", "WEST", "EASTERN", " w Est x")
                SPICE.pcpool("BODY299_PGR_POSITIVE_LON", [sense])
                pcpool(ctx, "BODY299_PGR_POSITIVE_LON", [sense])
                @test et2lst(ctx, 1.0e8, 299, 0.7, "PLANETOGRAPHIC") ==
                    SPICE.et2lst(1.0e8, 299, 0.7, "PLANETOGRAPHIC")
            end
            pcpool(ctx, "BODY299_PGR_POSITIVE_LON", ["NORTH"])
            @test throws_kind(() -> et2lst(ctx, 1.0e8, 299, 0.7, "PLANETOGRAPHIC"), :invalidoption)
        finally
            SPICE.dvpool("BODY299_PGR_POSITIVE_LON")
            dvpool(ctx, "BODY299_PGR_POSITIVE_LON")
        end
        @test throws_kind(() -> et2lst(ctx, 0.0, 1234567, 0.0, "PLANETOCENTRIC"), :cantfindframe)
        # Without `BODY<id>_PM` the sense of rotation is unknown. CSPICE N0067 signals the
        # short message cut to `SPICE(CANTGETROTATIONTYPE`, without its parenthesis.
        c2 = neo_context([Kernels.LSK, Kernels.DE440S])
        pcpool(c2, "OBJECT_301_FRAME", ["J2000"])
        @test throws_kind(() -> et2lst(c2, 1.0e8, 301, 0.3, "PLANETOCENTRIC"), :cantgetrotationtype)
        pm = SPICE.gdpool("BODY301_PM")
        try
            SPICE.pcpool("OBJECT_301_FRAME", ["J2000"])
            SPICE.dvpool("BODY301_PM")
            hr, mn, sc = Ref{Cint}(0), Ref{Cint}(0), Ref{Cint}(0)
            ccall(
                (:et2lst_c, SPICE.libcspice), Cvoid, (
                    Cdouble, Cint, Cdouble, Cstring, Cint, Cint,
                    Ref{Cint}, Ref{Cint}, Ref{Cint}, Ptr{UInt8}, Ptr{UInt8},
                ), 1.0e8, 301, 0.3, "PLANETOCENTRIC", 51, 51, hr, mn, sc, zeros(UInt8, 51),
                zeros(UInt8, 51)
            )
            buf = zeros(UInt8, 64)
            ccall((:getmsg_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Ptr{UInt8}), "SHORT", 64, buf)
            @test unsafe_string(pointer(buf)) == "SPICE(CANTGETROTATIONTYPE"
        finally
            ccall((:reset_c, SPICE.libcspice), Cvoid, ())
            SPICE.dvpool("OBJECT_301_FRAME")
            SPICE.pdpool("BODY301_PM", pm)
        end
        r = @inferred et2lst(ctx, 1.0e8, 399, 0.3, "PLANETOCENTRIC")
        @test r isa Tuple{Int, Int, Int, String, String}
    end

    @apart @testset "lspcn" begin
        d = Float64[]
        # Bodies with a segment in de440s.bsp and a pole in the PCK.
        for body in ("EARTH", 399, "MOON", 301, "VENUS", 199, "mercury"),
                ab in ("NONE", "LT", "LT+S", "CN", "CN+S")
            for et in random_epochs(rng, 10, -6.0e8, 6.0e8)
                want = SPICE.lspcn(string(body isa Integer ? bodc2n(ctx, body) : body), et, ab)
                got = lspcn(ctx, body, et, ab)
                x = abs(got - want)
                push!(d, min(x, 2pi - x))
            end
        end
        s = report_quantiles("lspcn: angle (rad)", d)
        @test s.max <= 2.0e-12
        @test 0 <= lspcn(ctx, "EARTH", 0.0, "NONE") < 2pi
        @test (@inferred lspcn(ks, 399, 1.0e8, "LT+S")) isa Float64
        # A plain epoch only, as et2lst: no derivatives, and no MethodError.
        @test throws_kind(() -> lspcn(ctx, 399, ForwardDiff.Dual(1.0e8, 1.0), "NONE"), :badarg)
    end

    @apart @testset "bltfrm and kplfrm" begin
        for cls in (-1, 1, 2, 3, 4, 5, 6)
            @test bltfrm(cls) == Int[i for i in SPICE.bltfrm(cls)]
            @test kplfrm(ctx, cls) == Int[i for i in SPICE.kplfrm(cls, 1000)]
        end
        @test length(kplfrm(ctx, -1)) > 20
        for cls in (0, 7, -2)
            @test throws_kind(() -> bltfrm(cls), :badframeclass)
            @test throws_kind(() -> kplfrm(ctx, cls), :badframeclass)
            @test_throws SPICE.SpiceError SPICE.kplfrm(cls, 10)
        end
        # The rules of `kplfrm_` on a pool of its own. CSPICE reads the variables, not the
        # frame definitions: KP_A has its class under its name; KP_B has no FRAME_KP_B and
        # KP_C no class, so both are left out; FRAME_KP_D is rounded to the ID of OTHER,
        # whose name does not match and which is listed; KP_E has a name variable in
        # lower case, which the template matches and the lookup of the ID does not.
        lines = [
            "FRAME_KP_A = 1700001", "FRAME_1700001_NAME = 'KP_A'", "FRAME_KP_A_CLASS = 3",
            "FRAME_1700002_NAME = 'KP_B'", "FRAME_1700002_CLASS = 4",
            "FRAME_KP_C = 1700003", "FRAME_1700003_NAME = 'KP_C'",
            "FRAME_1700005_NAME = 'KP_D'", "FRAME_KP_D = 1700004.4",
            "FRAME_1700004_NAME = 'OTHER'", "FRAME_1700004_CLASS = 4",
            "frame_1700006_name = 'KP_E'", "FRAME_KP_E = 1700006", "FRAME_1700006_CLASS = 4",
        ]
        try
            SPICE.kclear()
            SPICE.lmpool(lines)
            kc = SpiceContext()
            lmpool(kc, lines)
            for cls in (-1, 3, 4)
                @test kplfrm(kc, cls) == Int[i for i in SPICE.kplfrm(cls, 100)]
            end
            @test kplfrm(kc, -1) == [1700001, 1700004]
        finally
            SPICE.kclear()
        end
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
