# Group 54: the analysis layer, part 1 (`Orbit`, its constructors, the Cartesian and
# kinematic elements, name dispatch; src/analysis/orbit.jl and elements.jl).
#
# The reference values are those of the GMAT case of the tests of ANISE
# (val_state_def_circ_inc) for an Earth GM of 398600.4415 km^3/s^2. The kernel part
# checks that `orbit` and `Orbit(ctx, ...)` read the same state, GM and radii as
# `spkezr` and `bodvcd`.
#
# The Keplerian elements are compared with the GMAT cases of ANISE (gmat_cases.jl),
# with Almanac's `oscltx` (validated against CSPICE in group math_conics) on random
# orbits, and on the degenerate orbits, where each expectation follows the ANISE source
# (src/astro/orbit.rs, b3aaadc).
#
# Run time of the group: about 18 s (17-18 s measured, one process; most of it is the
# compilation of the Keplerian code and of the CSPICE oracle wrapper).

import ForwardDiff
isdefined(@__MODULE__, :CS) || include("mathoracle.jl")
include("gmat_cases.jl")

# The tolerances against GMAT: the testset prints the largest difference of each case
# (`measured`); the largest is 4.0e-13 (circ_inc_kep, elliptical_kep). ANISE prints the
# values with 16 to 17 digits and compares with 1e-10 degrees; the difference is the
# rounding of values of 1e4 to 1e5 km (ulp 2e-12 to 1.5e-11 per km of magnitude 1e5, here
# scaled by max(|v|, 1)) and of angles near 2 pi, with the ANISE operations in the same
# order. The tolerance is 1e-12, 2.5 times the largest.
const RTOL_GMAT = 1e-12
const ATOL_GMAT = 1e-12

# A random orbit from the CSPICE conic elements (`keplerian_orbit` comes in a later task):
# periapsis 6500 to 50000 km, ecc in [0.005, 0.95] or [1.05, 5], inclination in
# [0.02, pi - 0.02], angles uniform. The bounds keep away from the circular and the
# equatorial orbit, where the angles are ill conditioned and the two codes may differ.
function random_orbit(rng, gm=398600.4415)
    ecc = rand(rng) < 0.5 ? 0.005 + 0.945 * rand(rng) : 1.05 + 3.95 * rand(rng)
    rp = 6500.0 + 43500.0 * rand(rng)
    inc = 0.02 + (π - 0.04) * rand(rng)
    m0 = ecc < 1 ? 2π * rand(rng) : 6 * rand(rng) - 3
    elts = [rp, ecc, inc, 2π * rand(rng), 2π * rand(rng), m0, 0.0, gm]
    return Orbit(conics(elts, 0.0), 0.0, 399, 1, gm)
end

angdiff(a, b) = abs(rem2pi(a - b, RoundNearest))

# Tolerances against `oscltx`: the largest differences measured in 200 random orbits
# (seed SEED + 31) are rp 9.6e-15, ecc 2.7e-15, inc 1.3e-15, raan 4.0e-13, aop 1.1e-14,
# ta 2.2e-12, ma 1.1e-11, sma and period 2.8e-15. The limit is about five times the
# measured value; raan, ta and ma are larger because the two codes compute them from the
# eccentricity and node vectors by different formulas (acos and atan2 near their
# unit-circle ends lose the digits of the argument).
const TOL_OSC = Dict(:rp => 5e-14, :ecc => 1.5e-14, :inc => 7e-15, :raan => 2e-12,
    :aop => 6e-14, :ta => 1e-11, :ma => 5e-11, :sma_period => 1.5e-14)
# Identities between elements of one orbit, 200 orbits: measured c3 1.6e-16, mean_motion
# 2.8e-16, ta_dot 2.4e-16, semi_minor 2.0e-14 (the square of a difference of ulps of
# sqrt(1-e^2)), ra 8.9e-16, dec 1.6e-15, vdec 4.4e-16, ha 8.9e-16: rounding of a few
# ulps, limit 1e-13. aol = aop + ta inherits the conditioning of aop and ta (2.2e-12
# measured), limit 1e-11.
const TOL_ID = 1e-13
const TOL_AOL = 1e-11
# B-plane of the GMAT case: the expected values are those of GMAT (as in ANISE's test), not
# of ANISE, and ANISE accepts 1 km and 0.01 deg, so the limits are ANISE's. Measured: 0.65
# km at most (br), 2e-5 rad for the angle; GMAT's ephemeris and axes (the state is in ICRS
# axes, read here as J2000) are not those of the test, so the formula cannot be tightened
# against these values.
const TOL_BPLANE = Dict(:bplane_br => 1.0, :bplane_bt => 1.0, :bplane_mag => 1.0,
    :bplane_angle => _d(0.01))
# Identities on random hyperbolas (measured: mag 2.2e-16, turn 4.4e-16, periapsis 1.2e-15;
# each is a few roundings of the same quantities), limit 1e-13 (TOL_ID); |b| against the
# semi-minor axis is the same number up to the rounding of `sqrt(1 - 1/e^2)` squared.
const TOL_BMAG = 1e-12
const TOL_VINF = 1e-13
const TOL_VINF_RP = 1e-13
# ANISE cases: measured turn angle 2.2e-16 rad, periapsis 5.5e-12 km (the radius divides
# by cos(rho) near 1, so the relative error of 1e-16 of ~1e4 km shows as 1e-12 km).
const TOL_VINF_KM = 1e-10
# b_plane_davis: ha 1.1e-11 rad measured; ANISE's 1e-9 deg.
const TOL_DAVIS_HA = _d(1e-9)

group("54. analysis: Orbit and Cartesian elements", Kernels.LSK, Kernels.DE440S,
        Kernels.GM_TPC, Kernels.PCK11) do
    @testset "Orbit: construction and Cartesian elements" begin
        s = SVector(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0)
        o = Orbit(s, 0.0, 399, 1, 398600.4415)
        @test o.state === s && o.center === Int32(399) && o.frame === Int32(1)
        @test (o.x, o.vz, o.rmag) == (s[1], s[6], norm(s[1:3]))
        @test o.position == s[1:3] && o.velocity == s[4:6]
        @test o.hvec == cross(s[1:3], s[4:6])
        @test o.hx ≈ 35065.806679607005 && o.hz ≈ 24796.2925419       # GMAT, val_state_def_circ_inc
        @test o.energy ≈ -25.842247282849137
        @test o.light_time == o.rmag / clight()
        @test element(o, :hz) === o.hz
        @test all(n -> n in propertynames(o), elements())
        # every element is a Float64, or the degenerate error of ANISE (:ha of an ellipse)
        @test all(n -> (y = outcome(() -> element(o, n)); y isa Float64 || y === :degeneratecase),
            elements())
        @test throws_kind(() -> element(o, :foo), :badarg)
        @test throws_kind(() -> o.foo, :badarg)
        @test throws_kind(() -> element(o, :position), :badarg)       # vectors are properties only
        @test all(isnan, o.radii)
    end

    @testset "Orbit: the number type follows the state and the epoch" begin
        s = SVector(7000.0, 100.0, -50.0, 0.1, 7.5, 1.0)
        o = Orbit(s, ForwardDiff.Dual(1.0, 1.0), 399, 1, 398600.4415)
        @test o isa Orbit{<:ForwardDiff.Dual} && o.rmag isa ForwardDiff.Dual
        o32 = Orbit(Float32.(s), 0.0f0, 399, 1, 398600.4415)
        @test o32 isa Orbit{Float64} && o32.state isa SVector{6,Float64}
        @test Orbit(SVector(7000, 0, 0, 0, 7, 0), 0, 399, 1, 1.0e5) isa Orbit{Float64}
        # An ellipse and a hyperbola, so that each element has a state where it is defined;
        # the degenerate error of ANISE (:ha of an ellipse, :ea of a hyperbola) is not a
        # type instability.
        sh = SVector(7000.0, 100.0, -50.0, 0.1, 11.5, 1.0)
        for orb in (Orbit(s, 0.0, 399, 1, 398600.4415), o,
                Orbit(sh, 0.0, 399, 1, 398600.4415), Orbit(sh, ForwardDiff.Dual(1.0, 1.0), 399, 1, 398600.4415)),
                n in elements()
            r = Ref(n)
            T = typeof(orb).parameters[1]
            y = try
                @inferred element(orb, r[])
            catch e
                e isa Almanac.SpiceError && e.kind === :degeneratecase || rethrow()
                nothing
            end
            y === nothing || @test y isa T
        end
        @test throws_kind(() -> Orbit(s[1:5], 0.0, 399, 1, 1.0), :badarg)
        @test throws_kind(() -> Orbit(s, 0.0, 399, 1, 1.0, [1.0, 2.0]), :badarg)
    end

    @testset "Orbit from kernels" begin
        ctx = neo_context([Kernels.LSK, Kernels.DE440S, Kernels.GM_TPC, Kernels.PCK11])
        et = 1.0e8
        o = orbit(ctx, "MOON", et, "J2000", "NONE", "EARTH")
        @test o.state == spkezr(ctx, "MOON", et, "J2000", "NONE", "EARTH")[1]
        @test o.gm == bodvcd(ctx, 399, "GM")[1] && o.radii == bodvcd(ctx, 399, "RADII")
        @test Orbit(ctx, o.state, et, "EARTH", "J2000") == o
        @test throws_kind(() -> Orbit(ctx, o.state, et, -999, "J2000"), :kernelvarnotfound)
    end

    @testset "Keplerian elements: GMAT cases ($(c.name))" for c in GMAT_CASES
        o = Orbit(c.state, 0.0, 399, 1, c.gm)
        worst = 0.0
        for (k, v) in c.expected
            x = element(o, k)
            @test isapprox(k in GMAT_ANGLES ? v + rem2pi(x - v, RoundNearest) : x, v;
                rtol=RTOL_GMAT, atol=ATOL_GMAT)
            worst = max(worst, k in GMAT_ANGLES ? angdiff(x, v) : abs(x - v) / max(abs(v), 1))
        end
        measured(@sprintf("  gmat %-24s max scaled difference %.3e", c.name, worst))
    end

    @testset "Keplerian elements against oscltx, and identities" begin
        rng = Xoshiro(SEED + 31)
        names = (:rp, :ecc, :inc, :raan, :aop, :ta, :ma, :sma_period,
            :aol, :c3, :mean_motion, :ta_dot, :semi_minor, :ra, :dec, :vdec, :ha)
        worst = Dict(k => 0.0 for k in names)
        note(k, d) = (worst[k] = max(worst[k], d); d)
        for _ in 1:nsample(200, 5000)
            o = random_orbit(rng)
            el = oscltx(o.state, o.et, o.gm)
            @test note(:rp, abs(o.periapsis - el[1]) / el[1]) < TOL_OSC[:rp]
            @test note(:ecc, abs(o.ecc - el[2])) < TOL_OSC[:ecc]
            @test note(:inc, angdiff(o.inc, el[3])) < TOL_OSC[:inc]
            @test note(:raan, angdiff(o.raan, el[4])) < TOL_OSC[:raan]
            @test note(:aop, angdiff(o.aop, el[5])) < TOL_OSC[:aop]
            @test note(:ta, angdiff(o.ta, el[9])) < TOL_OSC[:ta]
            if o.ecc < 1
                @test note(:ma, angdiff(o.ma, el[6])) < TOL_OSC[:ma]
                @test note(:sma_period, max(abs(o.sma - el[10]) / el[10],
                    abs(o.period - el[11]) / el[11])) < TOL_OSC[:sma_period]
                @test note(:mean_motion, abs(o.mean_motion * o.period - 2π) / 2π) < TOL_ID
                @test note(:ta_dot, abs(o.ta_dot - o.mean_motion) / o.mean_motion) < TOL_ID
            else
                @test note(:mean_motion, abs(o.mean_motion - sqrt(o.gm / abs(o.sma)^3)) / o.mean_motion) < TOL_ID
                @test o.period == 0.0 && o.ta_dot == Inf
                # the hyperbolic anomaly H of its defining relations with ecc and ta
                H = rem2pi(o.ha, RoundNearest)
                c, s = cos(o.ta), sin(o.ta)
                @test note(:ha, max(abs(cosh(H) - (o.ecc + c) / (1 + o.ecc * c)) / cosh(H),
                    abs(sinh(H) - sqrt(o.ecc^2 - 1) * s / (1 + o.ecc * c)) / max(abs(sinh(H)), 1))) < TOL_ID
            end
            @test note(:aol, angdiff(o.aol, o.aop + o.ta)) < TOL_AOL
            @test note(:c3, abs(o.c3 - 2 * o.energy) / abs(o.c3)) < TOL_ID
            @test note(:semi_minor, abs(o.semi_minor_axis^2 - abs(o.sma * o.semi_parameter)) /
                o.semi_minor_axis^2) < TOL_ID
            _, ra, dec = recrad(o.position)
            @test note(:ra, angdiff(o.right_ascension, ra)) < TOL_ID
            @test note(:dec, abs(o.declination - dec)) < TOL_ID
            @test note(:vdec, abs(o.velocity_declination - recrad(o.velocity)[3])) < TOL_ID
        end
        measured("  oscltx and identities, largest differences: " *
            join((@sprintf("%s %.1e", k, worst[k]) for k in names), ", "))
    end

    @testset "Keplerian elements: degenerate cases" begin
        gm = 398600.4415
        orb(v) = Orbit(SVector{6,Float64}(v), 0.0, 399, 1, gm)
        circ = orb((7000.0, 0, 0, 0, 0, sqrt(gm / 7000)))      # polar, ecc ~ 0 (rounding only)
        pro = orb((7000.0, 0, 0, 0, 8.0, 0))                   # equatorial, inc == 0
        retro = orb((7000.0, 0, 0, 0, -8.0, 0))                # equatorial, inc == pi
        para = orb((7000.0, 0, 0, 0, 0, sqrt(2gm / 7000)))     # |ecc - 1| < 1e-11
        @test circ.ecc < 1e-11
        @test throws_kind(() -> circ.ta, :degeneratecase)
        @test throws_kind(() -> circ.ea, :degeneratecase) && throws_kind(() -> circ.ma, :degeneratecase)
        @test throws_kind(() -> circ.tlong, :degeneratecase) && throws_kind(() -> circ.fpa, :degeneratecase)
        # ANISE aop_deg: acos(n.e / (|n| ecc)) with ecc ~ 1e-16 is rounding noise; it is
        # still a number in [0, 2 pi), never an error.
        @test 0 <= circ.aop < 2π
        for o in (pro, retro)
            # ANISE raan_deg: n = ẑ×ĥ is zero, 0/0 is NaN, and NaN becomes 0 (not an error).
            @test o.raan == 0.0
            @test throws_kind(() -> o.aol, :degeneratecase)       # |ẑ×ĥ| < 1e-6 deg
        end
        @test pro.inc == 0.0 && retro.inc == Float64(π)
        # ANISE aop_deg of an equatorial orbit: cos = 0/0, NaN, becomes 0 degrees.
        @test pro.aop == 0.0 && retro.aop == 0.0
        @test pro.ta == 0.0                                      # periapsis at the state
        # ANISE main: sma_km is -mu / (2 energy), with no check; a parabola has no
        # exception there (the brief expected one). Only `ma` is an error (ParabolicEccentricity).
        # This relies on the ecc of this state rounding to exactly 1.0 and its energy to
        # exactly 0 (sma = -Inf), which the next test line checks.
        @test para.ecc == 1.0 && para.energy == 0.0
        @test throws_kind(() -> para.ma, :degeneratecase)
        @test abs(para.ecc - 1) < 1e-11 && !isnan(para.sma) && abs(para.sma) > 1e9
        @test para.period == 0.0 && para.ta_dot == Inf           # ANISE: Duration::ZERO, 360/0
        @test throws_kind(() -> para.ha, :degeneratecase)        # NotHyperbolic for ecc <= 1
        # Zero radius or velocity: hvec and evec of ANISE are errors.
        zr = orb((0.0, 0, 0, 0, 8.0, 0))
        zv = orb((7000.0, 0, 0, 0, 0, 0))
        for n in (:inc, :raan, :aop, :aol, :ecc, :ta, :sma), o in (zr,)
            @test throws_kind(() -> element(o, n), :degeneratecase)
        end
        @test throws_kind(() -> zv.inc, :degeneratecase)         # hvec needs a velocity
        @test zv.evec isa SVector{3,Float64}                     # evec needs only a radius
        @test throws_kind(() -> zr.evec, :degeneratecase)
        # ANISE: energy needs a radius; hmag, hx, hy, hz go through hvec() (radius and velocity).
        @test throws_kind(() -> zr.energy, :degeneratecase)
        for n in (:hmag, :hx, :hy, :hz), o in (zr, zv)
            @test throws_kind(() -> element(o, n), :degeneratecase)
        end
        # A radial orbit with a velocity (r parallel to v): h = 0, and ANISE gives NaN
        # for inc and aol (0/0); Almanac throws.
        rect = orb((7000.0, 0, 0, 5.0, 0, 0))
        for n in (:inc, :aol)
            @test throws_kind(() -> element(rect, n), :degeneratecase)
        end
        @test rect.hmag == 0.0
        # zero velocity at ecc rounding above 1: ha = 0/0 (ANISE NaN)
        zvh = orb((7000.0, 0, 0, 0, 0, 1.0e-10))
        @test outcome(() -> zvh.ha) isa Symbol || !isnan(zvh.ha)
        # ANISE returns NaN for these (0/0 or Inf*0 at ecc = 1); Almanac throws instead.
        for n in (:semi_parameter, :semi_minor_axis, :periapsis)           # sma = -Inf, ecc = 1
            @test throws_kind(() -> element(para, n), :degeneratecase)
        end
        @test para.apoapsis == -Inf                              # ANISE: -Inf * 2, not NaN
        for n in (:ea, :fpa)                                      # ta = pi at ecc = 1: 0/0
            @test throws_kind(() -> element(zv, n), :degeneratecase)
        end
        for o in (circ, pro, retro, para, zr, zv, rect, zvh), n in elements()   # no NaN escapes
            y = outcome(() -> element(o, n))
            @test y isa Symbol || !isnan(y)
        end
        @test keplerian(pro) == (sma=pro.sma, ecc=pro.ecc, inc=pro.inc, raan=pro.raan,
            aop=pro.aop, ta=pro.ta)
        @test keys(keplerian(pro)) == (:sma, :ecc, :inc, :raan, :aop, :ta)
    end

    @testset "B-plane: GMAT case, geometry and degenerate cases" begin
        ctx = neo_context([Kernels.LSK, Kernels.DE440S, Kernels.GM_TPC, Kernels.PCK11])
        for c in BPLANE_CASES
            et = str2et(ctx, c.utc)
            earth = Orbit(c.state, et, 399, 1, c.gm)
            @test abs(earth.c3 - c.c3) < 1e-5                         # ANISE's check of the case
            moon = Orbit(ctx, c.state - spkezr(ctx, "MOON", et, "J2000", "NONE", "EARTH")[1],
                et, "MOON", "J2000")
            diffs = String[]
            for k in (:bplane_bt, :bplane_br, :bplane_mag, :bplane_angle)
                d = abs(element(moon, k) - c.expected[k])
                @test d < TOL_BPLANE[k]
                push!(diffs, @sprintf("%s %.3e", k, d))
            end
            measured(@sprintf("  bplane %-20s differences to the GMAT values: %s", c.name,
                join(diffs, ", ")))
            @test moon.bplane_ltof isa Float64                         # no value in ANISE's test
        end
        # Prograde and retrograde equatorial flybys (ANISE's documentation of `BPlane`).
        pro = Orbit(gmat_kep_state(-5000.0, 2.0, 0.0, 0.0, 0.0, _d(-45)), 0.0, 399, 1, GMAT_GM)
        retro = Orbit(gmat_kep_state(-5000.0, 2.0, π, 0.0, 0.0, _d(-45)), 0.0, 399, 1, GMAT_GM)
        @test pro.bplane_bt > 0 && abs(pro.bplane_br) < 1e-9
        @test retro.bplane_bt < 0 && abs(retro.bplane_br) < 1e-9
        @test abs(pro.bplane_mag - pro.semi_minor_axis) < 1e-9 * pro.semi_minor_axis
        @test pro.bplane_ltof > 0                                       # inbound: before the crossing
        # The incoming leg far away: the linearized time of flight tends to r / v_inf, with
        # an error of order (b/r)^2 and mu / (r v_inf^2), below 1e-6 at r ~ 1e8 b.
        far = Orbit(gmat_kep_state(-5000.0, 2.0, _d(28.5), _d(10), _d(20), -acos(-1 / 2) + 1e-8),
            0.0, 399, 1, GMAT_GM)
        @test far.rmag > 1e10
        @test abs(far.bplane_ltof * sqrt(far.c3) / far.rmag - 1) < 1e-6
        # elliptic orbit: ANISE NotHyperbolic
        ell = Orbit(gmat_kep_state(8000.0, 0.3, _d(30), 0.0, 0.0, 1.0), 0.0, 399, 1, GMAT_GM)
        for n in (:bplane_bt, :bplane_br, :bplane_mag, :bplane_angle, :bplane_ltof)
            @test throws_kind(() -> element(ell, n), :degeneratecase)
        end
        # incoming asymptote along the z axis: ANISE InfiniteValue (t = s x z = 0)
        ν = acos(1 / 2)
        e_hat = SVector(sin(ν), 0.0, -cos(ν))
        vp = sqrt(GMAT_GM * 3 / 10_000)
        polar = Orbit(vcat(10_000.0 * e_hat, vp * SVector(-cos(ν), 0.0, -sin(ν))), 0.0, 399, 1, GMAT_GM)
        @test polar.ecc ≈ 2
        for n in (:bplane_bt, :bplane_br, :bplane_mag, :bplane_angle, :bplane_ltof)
            @test throws_kind(() -> element(polar, n), :degeneratecase)
        end
        # c3 <= 0 with ecc > 1 cannot be built from a state; a zero speed has no hvec
        @test throws_kind(() -> Orbit(SVector(7000.0, 0, 0, 0, 0, 0), 0.0, 399, 1, GMAT_GM).bplane_bt,
            :degeneratecase)
    end

    @testset "B-plane and v-infinity identities on random hyperbolas" begin
        rng = Xoshiro(SEED + 32)
        names = (:mag, :angle, :turn, :periapsis)
        worst = Dict(k => 0.0 for k in names)
        note(k, d) = (worst[k] = max(worst[k], d); d)
        nhyp = ninbound = 0
        for _ in 1:nsample(200, 5000)
            o = random_orbit(rng)
            o.ecc > 1 || continue
            nhyp += 1
            @test note(:mag, abs(o.bplane_mag - hypot(o.bplane_bt, o.bplane_br)) / o.bplane_mag) < TOL_ID
            @test note(:angle, angdiff(o.bplane_angle, atan(o.bplane_br, o.bplane_bt))) < TOL_ID
            @test -π <= o.bplane_angle < π
            @test abs(o.bplane_mag - o.semi_minor_axis) / o.semi_minor_axis < TOL_BMAG
            @test note(:turn, abs(vinf_turn_angle(o, o.periapsis) - 2 * asin(1 / o.ecc))) < TOL_VINF
            rp = 7000.0 + 30000.0 * rand(rng)
            @test note(:periapsis, abs(vinf_periapsis(o, vinf_turn_angle(o, rp)) - rp) / rp) < TOL_VINF_RP
            if dot(o.position, o.velocity) < 0            # inbound leg: Review Focus 1
                ninbound += 1
                @test π < o.ta < 2π && o.ma < 0
            else
                @test 0 < o.ta < π && o.ma > 0
            end
        end
        @test ninbound > nsample(200, 5000) / 10 && nhyp > nsample(200, 5000) / 4
        measured("  b-plane and v-infinity identities, largest differences: " *
            join((@sprintf("%s %.1e", k, worst[k]) for k in names), ", "))
    end

    @testset "v-infinity: ANISE cases ($(c.name))" for c in VINF_CASES
        o = Orbit(c.state, 0.0, 399, 1, c.gm)
        if haskey(c, :ha)
            measured(@sprintf("  vinf %-24s ha difference %.3e rad", c.name, abs(o.ha - c.ha)))
            @test abs(o.ha - c.ha) < TOL_DAVIS_HA
            phi = vinf_turn_angle(o, c.rp)
            @test abs(vinf_periapsis(o, phi) - c.rp) < 1e-10          # ANISE's reciprocity check
            measured(@sprintf("  vinf %-24s reciprocity %.3e km", c.name, abs(vinf_periapsis(o, phi) - c.rp)))
        else
            expected = 2 * asin(1 / o.ecc)
            @test abs(o.ecc - c.ecc) < 1e-12 && abs(o.sma - c.sma) < 1e-7
            d1 = abs(vinf_turn_angle(o, o.periapsis) - expected)
            d2 = abs(vinf_periapsis(o, expected) - o.periapsis)
            measured(@sprintf("  vinf %-24s turn angle %.3e rad, periapsis %.3e km", c.name, d1, d2))
            @test d1 < TOL_VINF && d2 < TOL_VINF_KM
        end
    end

    @testset "v-infinity: degenerate cases" begin
        h = Orbit(SVector(7000.0, 0, 0, 0, 12.0, 0), 0.0, 399, 1, GMAT_GM)
        e = Orbit(SVector(7000.0, 0, 0, 0, 8.0, 0), 0.0, 399, 1, GMAT_GM)
        for f in (vinf_periapsis, vinf_turn_angle)
            @test throws_kind(() -> f(e, 1.0), :degeneratecase)       # ANISE NotHyperbolic
        end
        @test vinf_periapsis(h, Float64(π)) == 0.0                    # cos(0) = 1: (1/1 - 1) mu/c3
        @test vinf_periapsis(h, 0.0) > 1e15                           # 1/cos(pi/2), cos(pi/2) = 6e-17 in floating point
        @test throws_kind(() -> vinf_periapsis(h, NaN), :degeneratecase)
        # 1 + c3 rp / mu = 0 makes the argument of the arc cosine infinite (ANISE: NaN)
        @test throws_kind(() -> vinf_turn_angle(h, -GMAT_GM / h.c3), :degeneratecase)
        @test throws_kind(() -> vinf_turn_angle(h, NaN), :degeneratecase)
        @test vinf_turn_angle(h, 0.0) == Float64(π)                            # rho = 0: a head-on pass turns by pi
        @test vinf_turn_angle(h, 8000.0) isa Float64
    end
end
