# Group 9: the SPK data types other than 2 against CSPICE.
#
# Type 2 is the type of `de440s.bsp` and `de430.bsp` (groups 1 to 4). No kernel of the
# other types, type 3 included, is on this machine, so the group writes them: one SPK file
# with one segment per case, made by `write_type_kernel` of `spktypefixtures.jl` with the
# segment writers of CSPICE (`spkwNN_c`, the `spk14b_c`/`spk14a_c`/`spk14e_c` sequence,
# and the Fortran-interface routines `spkw01_`, `spkw19_`, `spkw21_` for the types that
# have no C wrapper). Each case has its own target body, relative to the Earth (399), so
# a query needs no other kernel.
#
# The segment writers are called directly with `ccall`: the wrappers of SPICE.jl for
# some of them take the state count from the wrong dimension.
#
# Checks of each case:
#   (a) the `spkgeo` state against CSPICE at random epochs, at the two ends of the
#       segment and at the special epochs of the case (knots, record bounds), where the
#       record or window selection of CSPICE is decided;
#   (b) no allocation and inferred result types of `spkgeo` and `spkder`;
#   (c) orders 2 and 3 of `spkder` against the nested `ForwardDiff.derivative` of the
#       velocity part of `spkgeo` (the generic path with a `Dual` epoch), and order 2
#       against a central finite difference of the velocity.
# A second file (`write_zero_kernel`, type 2 segments written by Almanac) holds a
# component -0.0: the states of chains of one and two links keep the sign of CSPICE.
#
# Tolerance of (a). The two engines make the same operations in the same order on the
# same numbers, so the states differ by a few units in the last place. The bound is
# `TYPE_RTOL` on |difference| / |vector|.
#
# Tolerance of (c). The Taylor path and the dual path differentiate the same function,
# so they differ in rounding only: `AD_VEL_RTOL` per order. The finite difference has
# the step `h` of the case and the bound `FD_VEL_RTOL`; it is the independent check that
# both paths carry the derivative through the Kepler solutions of types 5, 15 and 17.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Almanac keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALMANAC_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff

isdefined(@__MODULE__, :write_type_kernel) || include("spktypefixtures.jl")

const AD_VEL_RTOL = (1.0e-10, 1.0e-8)   # orders 2 and 3 of spkder
const FD_VEL_RTOL = 1.0e-5

# --- comparisons ----------------------------------------------------------------------

"Epochs of the state comparison of `c`: ends, special epochs, points just inside the ends, random points."
function case_epochs(rng, c::TypeCase, nrandom::Int)
    span = c.t1 - c.t0
    special = length(c.epochs) > 400 ? c.epochs[1:(length(c.epochs) ÷ 400 + 1):end] : c.epochs
    near = [c.t0 + span * 1.0e-9, c.t1 - span * 1.0e-9, c.t0 + span * 1.0e-3, c.t1 - span * 1.0e-3]
    return vcat(c.t0, c.t1, special, near, random_epochs(rng, nrandom, c.t0, c.t1))
end

"Largest relative difference of position and of velocity between Almanac and CSPICE for `c`."
function state_errors(ctx, c::TypeCase, epochs)
    worst_r = 0.0
    worst_v = 0.0
    for et in epochs
        ref, _ = SPICE.spkgeo(c.body, et, "J2000", EARTH)
        state, _ = spkgeo(ctx, c.body, et, "J2000", EARTH)
        worst_r = max(worst_r, norm(state[1:3] - ref[1:3]) / norm(ref[1:3]))
        worst_v = max(worst_v, norm(state[4:6] - ref[4:6]) / norm(ref[4:6]))
    end
    return worst_r, worst_v
end

"Allocation of the second call of `f(args...)` (the first call compiles)."
function warm_alloc(f, args...)
    f(args...)
    return f(args...)
end

"Relative difference `|a - b| / |b|`; for `b == 0` the absolute difference."
function type_reldiff(a, b)
    nb = norm(b)
    return nb == 0 ? norm(a - b) : norm(a - b) / nb
end

alloc_type_geo(ctx, body, et) = @allocated spkgeo(ctx, body, et, 1, EARTH)
alloc_type_der(ctx, body, et, ::Val{N}) where {N} = @allocated spkder(ctx, body, et, 1, EARTH, Val(N))

"Epochs of `c` for the derivative checks: random, and not within `100 h` of a special epoch."
function smooth_epochs(rng, c::TypeCase, n::Int)
    out = Float64[]
    while length(out) < n
        t = c.t0 + (c.t1 - c.t0) * (0.02 + 0.96 * rand(rng))
        all(s -> abs(s - t) > 100 * c.h, c.epochs) && push!(out, t)
    end
    return out
end

"Target of the segments with a component -0.0; `ZERO_BODY - 1` is relative to it."
const ZERO_BODY = -8000

"""
    write_zero_kernel(path) -> path

Write the SPK file `path` with two type 2 segments of degree 0 over [0, 200] s whose z
coefficient is -0.0: `ZERO_BODY` relative to the Earth and `ZERO_BODY - 1` relative to
`ZERO_BODY`. In the first half of a record the Chebyshev sum of CSPICE gives z = -0.0.
"""
function write_zero_kernel(path)
    h = spkopn(path, "signed zeros", 0)
    rec = [1000.0, 2000.0, -0.0]
    spkw02(h, ZERO_BODY, EARTH, "J2000", 0.0, 200.0, "zero", 100.0, [rec, rec], 0.0)
    spkw02(h, ZERO_BODY - 1, ZERO_BODY, "J2000", 0.0, 200.0, "zero 2", 100.0, [rec, rec], 0.0)
    spkcls(h)
    return path
end

group("9. SPK data types against CSPICE") do
    rng = MersenneTwister(SEED + 9)
    path = joinpath(mktempdir(), "almanac_types.bsp")
    cases = write_type_kernel(path, rng)
    zpath = write_zero_kernel(joinpath(mktempdir(), "almanac_zero.bsp"))
    oracle_load([path, zpath])
    ctx = neo_context([path, zpath])
    ks = Almanac.kernelset(ctx)
    println("  ", length(cases), " segments in ", path)
    @test sort(unique(c.dtype for c in cases)) == [1, 3, 5, 8, 9, 12, 13, 14, 15, 17, 18, 19, 20, 21]

    @apart @testset "states against CSPICE" begin
        by_type = Dict{Int, Tuple{Float64, Float64}}()
        for c in cases
            wr, wv = state_errors(ctx, c, case_epochs(rng, c, 300))
            measured(@sprintf("  %-44s position %.2e  velocity %.2e", c.label, wr, wv))
            @test wr <= c.rtol
            @test wv <= c.rtol
            old = get(by_type, c.dtype, (0.0, 0.0))
            by_type[c.dtype] = (max(old[1], wr), max(old[2], wv))
        end
        for dtype in sort(collect(keys(by_type)))
            measured(
                @sprintf(
                    "  type %-2d  max relative difference: position %.2e  velocity %.2e",
                    dtype, by_type[dtype]...
                )
            )
        end
    end

    @apart @testset "name forms and spkezr give the spkgeo state" begin
        c = cases[1]
        et = (c.t0 + c.t1) / 2
        byid = spkgeo(ctx, c.body, et, "J2000", EARTH)
        @test spkezr(ctx, string(c.body), et, "J2000", "NONE", "EARTH") == byid
        @test spkpos(ctx, c.body, et, "J2000", "NONE", EARTH)[1] == byid[1][SOneTo(3)]
    end

    # CSPICE (SPKGEO) subtracts the two states in their common segment frame and rotates
    # the difference; rotating each state to J2000 first changes the last bit, which is
    # 3e-8 km at the distance of Mars from the solar system barycenter.
    @apart @testset "a chain in a frame other than J2000: the CSPICE state to the bit" begin
        a, b = (only(c for c in cases if c.label == "type $t, degree 3, ECLIPJ2000") for t in (9, 13))
        lo, hi = max(a.t0, b.t0), min(a.t1, b.t1)
        fails = Failures()
        for et in random_epochs(MersenneTwister(SEED + 90), 20, lo, hi)   # own stream
            check!(fails, spkgeo(ctx, a.body, et, "J2000", b.body)[1] == SPICE.spkgeo(a.body, et, "J2000", b.body)[1], "spkgeo", et)
        end
        @test isempty(fails)
    end

    # SPKGEO takes the state of the first link of each chain as it is: a component -0.0
    # stays -0.0, where a sum that starts from zero gives +0.0.
    @apart @testset "signed zeros of the first link of a chain: those of CSPICE" begin
        @test signbit(SPICE.spkgeo(ZERO_BODY, 10.0, "J2000", EARTH)[1][3])
        fails = Failures()
        for (t, o) in (
                    (ZERO_BODY, EARTH), (ZERO_BODY - 1, EARTH), (EARTH, ZERO_BODY),
                    (ZERO_BODY - 1, ZERO_BODY),
                ), et in (10.0, 60.0, 150.0)
            ref = collect(SPICE.spkgeo(t, et, "J2000", o)[1])
            check!(fails, isequal(collect(spkgeo(ctx, t, et, "J2000", o)[1]), ref), "spkgeo", t, o, et)
            q = state_query(ctx, t, "J2000", "NONE", o)
            check!(fails, isequal(collect(q(et)[1]), ref), "state_query", t, o, et)
        end
        @test isempty(fails)
    end

    # `spke15` splits the mean anomaly with the rounded `d_mod` of f2c. The exact `rem`
    # gave a turn of the periapsis and the node that is off by ulps of the mean anomaly:
    # 3.1e-13 relative from CSPICE after 3e4 revolutions and 3.4e-12 after 3e5. With
    # `d_mod` the difference is the last bit of `cos` and `sin`, 7.6e-16 at most.
    @apart @testset "type 15 far from the periapsis epoch: the mean anomaly as d_mod" begin
        path15 = joinpath(mktempdir(), "type15_far.bsp")
        tp = [sind(30.0) * sind(50.0), -cosd(30.0) * sind(50.0), cosd(50.0)]
        pa = normalize(cross([0.0, 0.0, 1.0], tp))
        p, ecc, epoch, span = 8000.0, 0.05, 1234.5, 3.0e9      # span: about 95 years
        period = 2π * sqrt((p / (1 - ecc^2))^3 / GM_EARTH)
        bodies = (-8501, -8502, -8503)
        h = SPICE.spkopn(path15, "TYPE 15", 0)
        for (body, j2flg) in zip(bodies, (0, 1, 2))
            write_type15!(
                h, body, "J2000", epoch - span, epoch + span, epoch, tp, pa, p, ecc,
                j2flg, [0.0, 0.0, 1.0], 1.0826e-3, 6378.14
            )
        end
        SPICE.spkcls(h)
        ctx15 = neo_context([path15])
        SPICE.furnsh(path15)
        try
            worst = 0.0
            for body in bodies, revs in (3.0e4, -3.0e4, 3.0e5, -3.0e5), k in 1:5
                et = epoch + revs * period * (1 + k / 1000)
                ref = SPICE.spkgeo(body, et, "J2000", EARTH)[1]
                st = spkgeo(ctx15, body, et, "J2000", EARTH)[1]
                worst = max(
                    worst, type_reldiff(st[SOneTo(3)], ref[1:3]),
                    type_reldiff(st[SVector(4, 5, 6)], ref[4:6])
                )
            end
            measured(@sprintf("  type 15 after 3e4 and 3e5 revolutions: largest relative difference %.2e", worst))
            @test worst <= 1.0e-14
        finally
            SPICE.unload(path15)
        end
    end

    @apart @testset "zero allocation and inferred types" begin
        # `skip` does not evaluate the expression, so the allocations are not measured
        # when STACK_BUFFERS_ESCAPE is set.
        fails = Failures()
        for c in cases
            et = c.t0 + 0.37 * (c.t1 - c.t0)
            @testset "$(c.label)" begin
                STACK_BUFFERS_ESCAPE || check!(fails, warm_alloc(alloc_type_geo, ctx, c.body, et) == 0, "spkgeo allocates", c.label)
                for N in 0:3
                    STACK_BUFFERS_ESCAPE || check!(fails, warm_alloc(alloc_type_der, ctx, c.body, et, Val(N)) == 0, "spkder allocates", c.label, N)
                end
                @test (@inferred spkgeo(ctx, c.body, et, 1, EARTH)) isa Tuple{SVector{6, Float64}, Float64}
                @test (@inferred spkder(ctx, c.body, et, 1, EARTH, Val(3))) isa NTuple{4, SVector{3, Float64}}
            end
        end
        @test isempty(fails) skip = STACK_BUFFERS_ESCAPE
    end

    @apart @testset "spkder: lower orders, and higher orders against ForwardDiff and differences" begin
        fails = Failures()
        for c in cases
            vel(t) = spkgeo(ctx, c.body, t, 1, EARTH)[1][SVector(4, 5, 6)]
            acc(t) = ForwardDiff.derivative(vel, t)
            jerk(t) = ForwardDiff.derivative(acc, t)
            worst = zeros(3)            # order 2, order 3 against AD; order 2 against FD
            for et in smooth_epochs(rng, c, nsample(3, 8))
                jet = spkder(ctx, c.body, et, 1, EARTH, Val(3))
                state, _ = spkgeo(ctx, c.body, et, 1, EARTH)
                # Equal to the last bit in a normal run; with coverage instrumentation the
                # compiler arranges the sums differently, so one rounding is allowed.
                check!(fails, isapprox(jet[1], state[SOneTo(3)]; rtol = 4eps()), "spkder order 0", c.label, et)
                check!(fails, isapprox(jet[2], state[SVector(4, 5, 6)]; rtol = 4eps()), "spkder order 1", c.label, et)
                check!(fails, isapprox(spkder(ctx, c.body, et, 1, EARTH, Val(2))[3], jet[3]; rtol = 1.0e-12), "spkder Val(2) order 2", c.label, et)
                # The dual path gives the same plain values.
                dstate = spkgeo(ctx, c.body, ForwardDiff.Dual(et, 1.0), 1, EARTH)[1]
                check!(fails, isapprox(ForwardDiff.value.(dstate), state; rtol = 1.0e-14), "dual spkgeo value", c.label, et)
                fd = (vel(et + c.h) - vel(et - c.h)) / (2 * c.h)
                worst .= max.(
                    worst, (
                        type_reldiff(jet[3], acc(et)), type_reldiff(jet[4], jerk(et)),
                        type_reldiff(jet[3], fd),
                    )
                )
            end
            measured(@sprintf("  %-44s AD order 2 %.1e, 3 %.1e; FD order 2 %.1e", c.label, worst...))
            @test worst[1] <= AD_VEL_RTOL[1]
            @test worst[2] <= AD_VEL_RTOL[2]
            @test worst[3] <= FD_VEL_RTOL
        end
        @test isempty(fails)
    end

    @apart @testset "epoch number types: Taylor epoch, and the Kepler solutions at their special epochs" begin
        fails = Failures()
        for c in cases
            seg = ks.spk[Int32(c.body)][1]
            data = ks.daf[seg.file]
            et = c.t0 + 0.61 * (c.t1 - c.t0)
            jet = Almanac.spk_eval(data, seg, et, Val(3))
            # A Taylor epoch of order 1 with Val(2): series of series. The plain parts
            # are the jet of order 2 and the first-order parts the next orders.
            nested = Almanac.spk_eval(data, seg, Almanac.taylor_var(et, Val(1)), Val(2))
            check!(fails, all(map(x -> x.v, nested[k]) ≈ jet[k] for k in 1:3), "nested Taylor values", c.label)
            check!(fails, isapprox(map(x -> x.d[1], nested[2]), jet[3]; rtol = 1.0e-9), "nested Taylor order 2", c.label)
            check!(fails, isapprox(map(x -> x.d[1], nested[3]), jet[4]; rtol = 1.0e-7), "nested Taylor order 3", c.label)
        end
        # Type 15 at the periapsis epoch (true anomaly zero, time interval zero) and type
        # 5 at a state epoch: the derivatives are finite and agree with the differences.
        for c in cases
            (c.dtype == 15 || c.dtype == 17 || c.dtype == 5) || continue
            et = c.dtype == 5 ? c.epochs[2] : c.epochs[1]
            vel(t) = spkgeo(ctx, c.body, t, 1, EARTH)[1][SVector(4, 5, 6)]
            jet = spkder(ctx, c.body, et, 1, EARTH, Val(3))
            check!(fails, all(all(isfinite, j) for j in jet), "finite jet", c.label)
            if c.dtype != 5       # at a type 5 state epoch the pair of states changes
                fd = (vel(et + c.h) - vel(et - c.h)) / (2 * c.h)
                within!(fails, type_reldiff(jet[3], fd), FD_VEL_RTOL, "order 2 against FD", c.label)
            end
            within!(fails, type_reldiff(jet[3], ForwardDiff.derivative(vel, et)), AD_VEL_RTOL[1], "order 2 against AD", c.label)
        end
        @test isempty(fails)
    end

    @apart @testset "numeric helpers" begin
        tv = Almanac.taylor_var
        der = Almanac.derivatives
        # sinh and cosh: the derivatives alternate.
        x = 0.7
        @test collect(der(sinh(tv(x, Val(4))))) ≈ [sinh(x), cosh(x), sinh(x), cosh(x), sinh(x)]
        @test collect(der(cosh(tv(x, Val(4))))) ≈ [cosh(x), sinh(x), cosh(x), sinh(x), cosh(x)]
        @test der(sinh(2 * tv(x, Val(2))))[3] ≈ 4 * sinh(2x)
        # atan of a point that turns at unit rate while its radius changes: the angle is t.
        for t in (0.3, 1.5707963, 2.5, -2.9)
            tt = tv(t, Val(4))
            radius = 2.0 + 0.5 * tt * tt
            ang = atan(radius * sin(tt), radius * cos(tt))
            @test collect(der(ang)) ≈ [t, 1.0, 0.0, 0.0, 0.0] atol = 1.0e-12
        end
        # Orders that a number type carries, and the replacement of the plain value.
        D1 = ForwardDiff.Dual{Nothing, Float64, 1}
        @test Almanac._ad_order(Float64) == 0
        @test Almanac._ad_order(Almanac.Taylor{3, Float64}) == 3
        @test Almanac._ad_order(D1) == 1
        @test Almanac._ad_order(ForwardDiff.Dual{Nothing, D1, 1}) == 2
        @test Almanac._ad_order(Almanac.Taylor{2, D1}) == 3
        @test Almanac._ad_order(Rational{Int}) == 8
        @test Almanac._with_primal(3.0, 1.5) === 1.5
        y = Almanac._with_primal(tv(3.0, Val(2)), 1.5)
        @test y.c == (1.5, 1.0, 0.0)
        d = Almanac._with_primal(ForwardDiff.Dual(3.0, 2.0), 1.5)
        @test ForwardDiff.value(d) == 1.5 && ForwardDiff.partials(d, 1) == 2.0
        @test Almanac._with_primal(1 // 2, 0.25) == 0.25
        # Stumpff functions: the three branches meet the definitions.
        for z in (-9.0, -0.5, 0.0, 0.5, 9.0)
            c0, c1, c2, c3 = Almanac._stumpff(z)
            r = sqrt(abs(z))
            e0, e1 = z < 0 ? (cosh(r), sinh(r) / r) : z > 0 ? (cos(r), sin(r) / r) : (1.0, 1.0)
            @test c0 ≈ e0 && c1 ≈ e1
            @test 1 - z * c2 ≈ c0 && 1 - z * c3 ≈ c1
        end
        @test throws_kind(() -> Almanac._stumpff(-1.0e6), :badorbit)
        # Angle between vectors and rotation about an axis.
        ex, ey, ez = SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0)
        @test Almanac._m_vsep(ex, ey) == pi / 2
        @test Almanac._m_vsep(ex, -ex) ≈ pi
        @test Almanac._m_vsep(ex, SVector(1.0, 1.0, 0.0)) ≈ pi / 4
        @test Almanac._m_vsep(ex, zero(ex)) == 0.0
        @test Almanac._m_vrotv(ex, ez, pi / 2) ≈ ey
        @test Almanac._m_vrotv(ex, zero(ez), 1.0) == ex
        @test Almanac._m_vrotv(zero(ex), ez, 1.0) == zero(ex)
        @test Almanac._dsign(2.0, -0.0) == 2.0 && Almanac._dsign(-2.0, -1.0) == -2.0
        # Two-body propagation: zero time gives the state back; bad input is refused.
        r0, v0 = SVector(7000.0, 0.0, 0.0), SVector(0.0, 7.5, 1.0)
        @test Almanac._prop2b(GM_EARTH, r0, v0, 0.0) == (r0, v0)
        pd, vd = Almanac._prop2b(GM_EARTH, r0, v0, ForwardDiff.Dual(0.0, 1.0))
        @test ForwardDiff.value.(pd) ≈ r0 && ForwardDiff.partials.(pd, 1) ≈ v0
        @test ForwardDiff.partials.(vd, 1) ≈ orbit_acc(r0)
        for dt in (-5.0e4, 300.0, 9.0e5)
            ref = orbit_state(collect(vcat(r0, v0)), dt)
            p, v = Almanac._prop2b(GM_EARTH, r0, v0, dt)
            @test p ≈ ref[1:3] rtol = 1.0e-13
            @test v ≈ ref[4:6] rtol = 1.0e-13
        end
        # The kinds of `prop2b`, in its order.
        @test throws_kind(() -> Almanac._prop2b(-1.0, zero(r0), v0, 1.0), :nonpositivemass)
        @test throws_kind(() -> Almanac._prop2b(GM_EARTH, zero(r0), zero(v0), 1.0), :zeroposition)
        @test throws_kind(() -> Almanac._prop2b(GM_EARTH, r0, zero(v0), 1.0), :zerovelocity)
        @test throws_kind(() -> Almanac._prop2b(GM_EARTH, r0, 2.0 * r0, 1.0), :nonconicmotion)
        @test throws_kind(() -> Almanac._prop2b(GM_EARTH, r0, 3.0 * v0, 1.0e308), :dtoutofrange)
        @test Almanac._kepler_equinoctial(0.0, 0.3) == 0.0
        @test throws_kind(() -> Almanac._kepler_equinoctial(0.8, 0.8), :badorbit)
    end

    @apart @testset "errors" begin
        # Type 10 (two-line elements) is not available.
        seg10 = Almanac.Segment(
            Int32(-1), Int32(399), Int32(1), Int32(10), 0.0, 1.0, Int32(0),
            1, 20, 0.0, 0.0, Int32(0), Int32(0)
        )
        @test throws_kind(() -> Almanac.spk_eval(zeros(20), seg10, 0.5, Val(1)), :notsupported)
        # A segment over changed data takes the layout of those data, as a load does.
        relayout(bad, seg) = Almanac._segment_layout(bad, :SPK, seg)
        # A count at the end of a segment that does not agree with its length. (For
        # type 3 the count is checked when the file is opened, as for type 2.)
        for c in cases
            c.dtype in (3, 15, 17, 20) && continue
            seg = ks.spk[Int32(c.body)][1]
            bad = copy(ks.daf[seg.file])
            bad[seg.eaddr] += 1.0
            @test throws_kind(() -> Almanac.spk_eval(bad, relayout(bad, seg), (c.t0 + c.t1) / 2, Val(1)), :baddaf)
        end
        # One wrong value in the segment: `change(data, seg)` makes it, and the query
        # in the middle of the segment must throw an error of `kind`.
        function refuses(label, kind, change)
            c = cases[findfirst(x -> x.label == label, cases)]
            seg = ks.spk[Int32(c.body)][1]
            bad = copy(ks.daf[seg.file])
            change(bad, seg)
            return throws_kind(() -> Almanac.spk_eval(bad, relayout(bad, seg), (c.t0 + c.t1) / 2, Val(1)), kind)
        end
        e(seg) = seg.eaddr
        b(seg) = seg.baddr
        @test refuses("type 8, degree 7", :baddaf, (d, s) -> d[e(s) - 2] = 0.0)             # step
        @test refuses("type 8, degree 7", :baddaf, (d, s) -> d[e(s) - 1] = 40.0)            # window
        @test refuses("type 8, degree 7", :baddaf, (d, s) -> d[e(s)] = NaN)
        @test refuses("type 9, degree 9, J2000", :baddaf, (d, s) -> d[e(s) - 1] = 500.0)
        @test refuses(
            "type 9, degree 9, J2000", :baddaf,                                 # equal epochs
            (d, s) -> d[(b(s) + 6 * 230):(b(s) + 7 * 230 - 1)] .= 0.0
        )
        @test refuses(
            "type 13, degree 7, J2000", :baddaf,
            (d, s) -> d[(b(s) + 6 * 230):(b(s) + 7 * 230 - 1)] .= 0.0
        )
        # The values that the readers and evaluators of CSPICE check have their kinds.
        @test refuses("type 18, subtype 0, degree 7", :invalidvalue, (d, s) -> d[e(s) - 1] = 3.0)  # odd window
        @test refuses("type 18, subtype 0, degree 7", :invalidvalue, (d, s) -> d[e(s) - 1] = 10.0)
        @test refuses("type 18, subtype 0, degree 7", :notsupported, (d, s) -> d[e(s) - 2] = 2.0)  # subtype
        @test refuses("type 19, 130 intervals", :baddaf, (d, s) -> d[(e(s) - 2 - 131):(e(s) - 2)] .= 1.0)
        @test refuses("type 1, order 15", :baddaf, (d, s) -> d[b(s):(b(s) + 71 * 130 - 1)] .= 0.0)
        @test refuses("type 21, order 25", :difflinetoolarge, (d, s) -> d[e(s) - 1] = 26.0)
        @test refuses("type 14, degree 9", :notsupported, (d, s) -> d[e(s) - 12] = 2.0)       # RDRTYP
        @test refuses("type 14, degree 9", :baddaf, (d, s) -> d[e(s) - 2] = 7.0)             # PKTSZ
        @test refuses("type 20, degree 8", :baddaf, (d, s) -> d[e(s) - 5] = 0.0)             # TSCALE
        @test refuses("type 15, flag 0, ecc 0.2", :badlatusrectum, (d, s) -> d[b(s) + 7] = -1.0)   # P
        @test refuses("type 15, flag 0, ecc 0.2", :badvector, (d, s) -> d[(b(s) + 1):(b(s) + 3)] .= 0.0)
        @test refuses(
            "type 15, flag 0, ecc 0.2", :badinitstate,                          # PA along TP
            (d, s) -> d[(b(s) + 4):(b(s) + 6)] .= d[(b(s) + 1):(b(s) + 3)]
        )
        @test refuses("type 17, h 0.12, k -0.09", :badeccentricity, (d, s) -> d[b(s) + 2] = 0.95)  # H
        @test refuses("type 17, h 0.12, k -0.09", :badsemiaxis, (d, s) -> d[b(s) + 1] = 0.0)       # A
        @test refuses("type 5, ellipse", :nonpositivemass, (d, s) -> d[e(s) - 1] = 0.0)            # GM
        # Type 20: a record size of 3 leaves no coefficient (`spke20` INVALIDCOUNT). The
        # loader refuses such a directory, so the segment is built here.
        rec20 = [0.0, 0.0, 0.0, 1.0, 1.0, 2451545.0, 0.0, 1.0, 3.0, 1.0]
        seg20 = Almanac.Segment(
            Int32(-1), Int32(399), Int32(1), Int32(20), 0.0, 86400.0, Int32(0),
            1, 10, 0.0, 86400.0, Int32(3), Int32(1)
        )
        @test throws_kind(() -> Almanac.spk_eval(rec20, seg20, 100.0, Val(1)), :invalidcount)
    end

    # The layouts of types 14 and 19 are checked at load, but a bad segment fails at its
    # evaluation, as in CSPICE: the file loads, and a type 19 segment with one bad
    # mini-segment evaluates in the others.
    @apart @testset "bad meta data load and fail at evaluation" begin
        c14 = cases[findfirst(x -> x.label == "type 14, degree 9", cases)]
        c19 = cases[findfirst(x -> x.label == "type 19, 130 intervals", cases)]
        s14, s19 = ks.spk[Int32(c14.body)][1], ks.spk[Int32(c19.body)][1]
        @assert s14.file == s19.file
        bad = copy(ks.daf[s14.file])
        bad[s14.eaddr - 2] = 7.0                               # PKTSZ of type 14
        last = Int(bad[s19.eaddr - 2]) + s19.baddr - 2         # end of the last mini-segment
        bad[last - 2] = 5.0                                    # its subtype
        badpath = joinpath(mktempdir(), "bad_meta.bsp")
        write(badpath, bad)
        cbad = neo_context([badpath])
        @test throws_kind(() -> spkgeo(cbad, c14.body, (c14.t0 + c14.t1) / 2, 1, EARTH), :baddaf)
        @test spkgeo(cbad, c19.body, c19.t0, 1, EARTH) == spkgeo(ctx, c19.body, c19.t0, 1, EARTH)
        @test throws_kind(() -> spkgeo(cbad, c19.body, c19.t1, 1, EARTH), :notsupported)
    end

    # Regression: the evaluators threw `:baddaf` for the values that CSPICE refuses with
    # kinds of its own, and type 19 compared its selection flag with 1.0 where `spkr19`
    # rounds it. Each case writes a copy of the file with `change(words, seg)` applied to
    # the words of the file (word `a` is DAF address `a`) and queries both engines in the
    # middle of the segment.
    @apart @testset "segment values refused with the kinds of CSPICE" begin
        dir = mktempdir()
        function bad_file(label, tag, change)
            c = cases[findfirst(x -> x.label == label, cases)]
            bytes = read(path)
            change(reinterpret(Float64, bytes), ks.spk[Int32(c.body)][1])
            file = joinpath(dir, tag * ".bsp")
            write(file, bytes)
            return c, file
        end
        function both_kinds(label, tag, change)
            c, file = bad_file(label, tag, change)
            et = (c.t0 + c.t1) / 2
            kn = try
                spkgeo(neo_context([file]), c.body, et, "J2000", EARTH)
                :none
            catch err
                err isa Almanac.SpiceError || rethrow()
                err.kind
            end
            SPICE.furnsh(file)
            kc = try
                ccall(
                    (:spkgeo_c, LIBCSPICE), Cvoid,
                    (CInt, Cdouble, Cstring, CInt, Ptr{Cdouble}, Ref{Cdouble}),
                    c.body, et, "J2000", EARTH, zeros(6), Ref(0.0)
                )
                cspice_kind()
            finally
                SPICE.unload(file)
            end
            return kn, kc
        end
        e(seg) = seg.eaddr
        b(seg) = seg.baddr
        function zero_steps(d, s)                   # G(1) of every difference line
            m, nrec = Int(d[e(s) - 1]), Int(d[e(s)])
            for i in 0:(nrec - 1)
                d[b(s) + i * (4m + 11) + 1] = 0.0
            end
        end
        function mini_subtypes(d, s)                # the subtype of every mini-segment
            n = Int(d[e(s)])
            pbas = e(s) - 2 - (n + 1)
            for i in 1:n
                d[Int(d[pbas + i + 1]) + b(s) - 4] = 3.0
            end
        end
        for (label, tag, change, kind) in (
                ("type 21, order 25", "t21_dim", (d, s) -> d[e(s) - 1] = 26.0, :difflinetoolarge),
                ("type 21, order 25", "t21_step", zero_steps, :zerostep),
                ("type 18, subtype 0, degree 7", "t18_count", (d, s) -> d[e(s)] = 1.0, :toofewstates),
                ("type 18, subtype 0, degree 7", "t18_subtype", (d, s) -> d[e(s) - 2] = 2.0, :notsupported),
                ("type 18, subtype 0, degree 7", "t18_window", (d, s) -> d[e(s) - 1] = 10.0, :invalidvalue),
                ("type 19, 130 intervals", "t19_subtype", mini_subtypes, :notsupported),
                ("type 15, flag 0, ecc 0.2", "t15_p", (d, s) -> d[b(s) + 7] = -1.0, :badlatusrectum),
                ("type 17, h 0.12, k -0.09", "t17_a", (d, s) -> d[b(s) + 1] = 0.0, :badsemiaxis),
                ("type 5, ellipse", "t5_gm", (d, s) -> d[e(s) - 1] = 0.0, :nonpositivemass),
            )
            kn, kc = both_kinds(label, tag, change)
            @test (tag, kn) == (tag, kind)
            @test (tag, kc) == (tag, kind)
        end
        # Type 19: a flag of 1.2 rounds to 1 (the later interval at a bound), as in `spkr19`.
        c = cases[findfirst(x -> x.label == "type 19, 130 intervals", cases)]
        seg = ks.spk[Int32(c.body)][1]
        d0 = ks.daf[seg.file]
        n = Int(d0[e(seg)])
        ivbas = e(seg) - 2 - (n + 1) - n ÷ 100 - (n + 1)
        bound = d0[ivbas + n ÷ 2]                   # a bound between two intervals
        flag(x) = (d, s) -> d[e(s) - 1] = x
        _, f1 = bad_file(c.label, "t19_flag1", flag(1.0))
        _, f12 = bad_file(c.label, "t19_flag12", flag(1.2))
        s1 = spkgeo(neo_context([f1]), c.body, bound, "J2000", EARTH)[1]
        s12 = spkgeo(neo_context([f12]), c.body, bound, "J2000", EARTH)[1]
        @test s12 == s1
        SPICE.furnsh(f12)
        try
            ref = SPICE.spkgeo(c.body, bound, "J2000", EARTH)[1]
            @test type_reldiff(s12[SOneTo(3)], ref[1:3]) <= c.rtol
            @test type_reldiff(s12[SVector(4, 5, 6)], ref[4:6]) <= c.rtol
        finally
            SPICE.unload(f12)
        end
    end

    # Regression: the counts of a segment were rounded with ties away from zero, where
    # `i_dnnt` of f2c takes the floor of `x + 0.5`: 0.49999999999999994 + 0.5 rounds to
    # 1.0, so `spkr08_` reads the degree 1 (a window of two states), not 0.
    @apart @testset "segment counts as i_dnnt" begin
        c8 = cases[findfirst(x -> x.label == "type 8, degree 7", cases)]
        seg8 = ks.spk[Int32(c8.body)][1]
        whole, near = copy(ks.daf[seg8.file]), copy(ks.daf[seg8.file])
        whole[seg8.eaddr - 1] = 1.0
        near[seg8.eaddr - 1] = prevfloat(0.5)
        for et in (c8.t0 + (c8.t1 - c8.t0) / 3, c8.t1 - (c8.t1 - c8.t0) / 7)
            @test Almanac.spk_eval(near, seg8, et, Val(1)) == Almanac.spk_eval(whole, seg8, et, Val(1))
        end
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
