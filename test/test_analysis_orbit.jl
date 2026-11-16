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
# Run time of the group: about 13 s (12.6 s measured, one process; most of it is the
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
end
