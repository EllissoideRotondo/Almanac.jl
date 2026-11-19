# Group 2: rotations against CSPICE.
#
# Error measure: absolute difference per matrix element. Elements of a rotation matrix
# are at most 1 in magnitude, so the absolute difference is also the relative one. For
# the derivative block of `sxform` the difference is divided by `RATE_SCALE`, an upper
# bound of the angular rates involved, so both blocks share one threshold.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Almanac keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALMANAC_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

const ROT_FRAMES = [
    "J2000", "ECLIPJ2000", "B1950", "GALACTIC", "IAU_EARTH", "IAU_MOON",
    "IAU_MARS", "IAU_JUPITER", "MOON_PA", "MOON_ME",
]
const ROT_INERTIAL = ("J2000", "ECLIPJ2000", "B1950", "GALACTIC")

# Inertial to inertial: constant matrices, a few roundings of difference at most.
const ROT_TOL_INERTIAL = (q99 = 2.0e-15, max = 5.0e-15)
# With a body-fixed frame: see the note on TOL_BODYFIXED in test_states.jl. Measured
# maximum 3.3e-16.
const ROT_TOL_BODYFIXED = (q99 = 5.0e-16, max = 2.0e-15)
# Largest angular rate among the frames, rad/s (Jupiter 1.76e-4, pairs can add).
const RATE_SCALE = 4.0e-4
const ORTHO_TOL = 1.0e-14

group("2. rotations against CSPICE", Kernels.TEXT_SET_FILES..., Kernels.MOON_BPC) do
    rng = MersenneTwister(SEED + 2)
    set = [Kernels.text_set()..., Kernels.MOON_BPC]
    oracle_load(set)
    ctx = neo_context(set)
    epochs = random_epochs(rng, 40, DE440S_COVERAGE...)

    din = Float64[]; dbf = Float64[]          # pxform differences
    sin_ = Float64[]; sbf = Float64[]         # sxform derivative-block differences, scaled
    ortho = Float64[]
    structure_ok = true
    for from in ROT_FRAMES, to in ROT_FRAMES
        from == to && continue
        inertial = from in ROT_INERTIAL && to in ROT_INERTIAL
        for et in epochs
            R = pxform(ctx, from, to, et)
            X = sxform(ctx, from, to, et)
            Rc = SPICE.pxform(from, to, et)
            Xc = SPICE.sxform(from, to, et)
            append!(inertial ? din : dbf, abs.(R .- Rc))
            append!(inertial ? din : dbf, abs.(X[SOneTo(3), SOneTo(3)] .- Xc[1:3, 1:3]))
            append!(
                inertial ? sin_ : sbf,
                abs.(X[SVector(4, 5, 6), SOneTo(3)] .- Xc[4:6, 1:3]) ./ RATE_SCALE
            )
            push!(ortho, maximum(abs.(R * R' - I)))
            push!(ortho, abs(det(R) - 1))
            # sxform = [R 0; dR R], with R the pxform matrix to within rounding (the
            # two come from the order 1 and the order 0 evaluation; measured 1.1e-16).
            structure_ok &= maximum(abs.(X[SOneTo(3), SOneTo(3)] - R)) <= 4eps()
            structure_ok &= X[SVector(4, 5, 6), SVector(4, 5, 6)] == X[SOneTo(3), SOneTo(3)]
            structure_ok &= all(iszero, X[SOneTo(3), SVector(4, 5, 6)])
            inertial && (structure_ok &= all(iszero, X[SVector(4, 5, 6), SOneTo(3)]))
        end
    end
    a = report_quantiles("pxform, inertial pairs, absolute", din)
    b = report_quantiles("pxform, body-fixed pairs, absolute", dbf)
    c = report_quantiles("sxform dR block / $(RATE_SCALE), body-fixed pairs", sbf)
    o = report_quantiles("orthogonality |R R' - I|, |det R - 1|", ortho)
    @test a.q99 <= ROT_TOL_INERTIAL.q99
    @test a.max <= ROT_TOL_INERTIAL.max
    @test b.q99 <= ROT_TOL_BODYFIXED.q99
    @test b.max <= ROT_TOL_BODYFIXED.max
    @test c.q99 <= ROT_TOL_BODYFIXED.q99
    @test c.max <= ROT_TOL_BODYFIXED.max
    @test maximum(sin_) == 0.0
    @test o.max <= ORTHO_TOL
    @test structure_ok

    @apart @testset "forms and types" begin
        et = epochs[1]
        R = pxform(ctx, "J2000", "IAU_EARTH", et)
        @test R isa SMatrix{3, 3, Float64, 9}
        @test sxform(ctx, "J2000", "IAU_EARTH", et) isa SMatrix{6, 6, Float64, 36}
        @test pxform(ctx, 1, SPICE.namfrm("IAU_EARTH"), et) == R
        @test pxform(ctx, "j2000", " iau_earth ", et) == R
        @test pxform(ctx, "IAU_MARS", "IAU_MARS", et) == SMatrix{3, 3, Float64}(I)
        # The alias frames give the DE440 frames exactly as CSPICE resolves them.
        @test isapprox(pxform(ctx, "MOON_PA", "MOON_PA_DE440", et), I; atol = 1.0e-15)
        @test isapprox(pxform(ctx, "MOON_ME", "MOON_ME_DE440_ME421", et), I; atol = 1.0e-15)
        @test isapprox(
            pxform(ctx, "MOON_PA", "MOON_ME", et),
            SPICE.pxform("MOON_PA", "MOON_ME", et); atol = 1.0e-15
        )
    end

    # The same frames near J2000, where the angle W is small: this isolates the model
    # (pole, nutation-precession terms, Chebyshev angles) from the rounding of W.
    @apart @testset "near J2000" begin
        near = random_epochs(rng, 40, -86400.0 * 30, 86400.0 * 30)
        dn = Float64[]
        for from in ROT_FRAMES, to in ROT_FRAMES, et in near
            from == to && continue
            append!(dn, abs.(pxform(ctx, from, to, et) .- SPICE.pxform(from, to, et)))
        end
        n = report_quantiles("pxform, all pairs within 30 d of J2000", dn)
        @test n.q99 <= 5.0e-16
        @test n.max <= 2.0e-15
    end

    @apart @testset "built-in inertial pairs: one matrix on every path" begin
        # The pairs of the 21 built-in inertial frames skip the chain lookups; the value
        # is the same for a plain epoch, a series epoch and the jet, and the derivatives
        # are exact zeros.
        ks = Almanac.kernelset(ctx)
        et = 1.2345e8
        same = true; zeroder = true
        for i in Int32(1):Int32(21), j in Int32(1):Int32(21)
            R = Almanac.frame_rotation(ks, i, j, et)
            jet = Almanac.rotation_jet(ks, i, j, et, Val(2))
            tay = Almanac.frame_rotation(ks, i, j, Almanac.taylor_var(et, Val(2)))
            same &= jet[1] == R && map(x -> Almanac.primal(x), tay) == R
            zeroder &= all(m -> all(iszero, m), jet[2:3]) &&
                all(x -> all(iszero, Almanac.derivatives(x)[2:3]), tay)
        end
        @test same
        @test zeroder
    end

    @apart @testset "pxfrm2: rotation between two epochs" begin
        d = Float64[]
        for _ in 1:200
            from, to = rand(rng, ROT_FRAMES), rand(rng, ROT_FRAMES)
            e1 = rand(rng, epochs)
            e2 = e1 + 1.0e4 * randn(rng)
            append!(d, abs.(pxfrm2(ctx, from, to, e1, e2) .- SPICE.pxfrm2(from, to, e1, e2)))
        end
        r = report_quantiles("pxfrm2, all pairs at two epochs", d)
        @test r.q99 <= ROT_TOL_BODYFIXED.q99
        @test r.max <= ROT_TOL_BODYFIXED.max
        et = epochs[2]
        # with one epoch: the rotation of pxform, through J2000
        @test isapprox(
            pxfrm2(ctx, "IAU_EARTH", "IAU_MOON", et, et),
            pxform(ctx, "IAU_EARTH", "IAU_MOON", et); atol = 1.0e-15
        )
        @test pxfrm2(ctx, 1, SPICE.namfrm("IAU_MARS"), et, et + 60.0) ==
            pxfrm2(ctx, "J2000", "IAU_MARS", et - 1.0e5, et + 60.0)
        # measured: SPICE(UNKNOWNFRAME) for either name
        @test throws_kind(() -> pxfrm2(ctx, "NO_SUCH_FRAME", "J2000", et, et), :unknownframe)
        @test throws_kind(() -> pxfrm2(ctx, "J2000", "NO_SUCH_FRAME", et, et), :unknownframe)
    end

    # Regression: the rate of a binary PCK frame is that of `pckmat_` and `tisbod_`:
    # the angles and their rates from `chbint_` go to `eul2xf_`; it was the derivative
    # of a `Taylor` evaluation, which differed in the last bits of every derivative
    # block. `eul2xf_` uses `sin` and `cos`, which differ in the last bit between Julia
    # and the C library on some arguments; with the C library swapped in all agree to
    # the last bit. Measured with Julia's functions (500 epochs): 430 equal to the last
    # bit (old arithmetic 0). Own random stream.
    @apart @testset "binary PCK rates as tisbod_" begin
        r2 = MersenneTwister(SEED + 201)
        n = 200
        nexact = 0
        for et in random_epochs(r2, n, DE440S_COVERAGE...), from in ("J2000", "ECLIPJ2000")
            nexact += sxform(ctx, from, "MOON_PA", et) == SPICE.sxform(from, "MOON_PA", et)
        end
        @test nexact >= 0.6 * 2n
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
