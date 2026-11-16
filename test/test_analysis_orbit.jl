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
# (src/astro/orbit.rs, b3aaadc). The altitudes and geodetic elements are compared with
# Vallado's examples as ANISE prints them and with `recgeo`; the equinoctial ones with their
# definitions; the Brouwer mean elements with the round trips of ANISE's own tests
# (orbit_mean_elements.rs). `keplerian_orbit` is compared with `conics` on random orbits,
# `with_element`, `with_apsides` with the loops of ANISE's verif_with_init, and `two_body` with
# `prop2b` and with ANISE's verif_orbit_at_epoch. The local frames are checked for orthonormality,
# against the axes as ANISE lists them, and the exact rate against ForwardDiff of `two_body`.
#
# Run time of the group: 16 to 30 s (one process; 23.8 to 29.6 s measured with the load of the
# machine; 29.6 s in the last run before the commit). About 12 s of it is the load of the
# oracle and the fixtures before the first testset.

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

# A random orbit from the CSPICE conic elements (not from `keplerian_orbit`, which would
# change the rounding of the sampled orbits and the earlier measured tolerances, and is itself
# tested against `conics` below):
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

# Conditioning of the shape: periapsis radius, eccentricity and semi-major axis come from
# cancellations near ecc = 1 (1 - e in a(1 - e) and in p / (1 - e); 1 - e^2 in the semi-minor
# axis), so the rounding is amplified by 1/|1 - e| (power 1) or 1/|1 - e^2| (power 2).
# A state rebuilt from its elements: the angles above, and sma and ecc through the semi-parameter
# p = a (1 - e^2), whose rounding is amplified by 1/|1 - e^2|. (The Kepler iteration is not in
# this path: `o` is the reference and the rebuilt states start from its true anomaly.)
cond_state(o) = cond_angle(o, :state) + eps() * 4 * o.ecc^2 / abs(1 - o.ecc^2)

# Error of a state that comes out of Kepler's equation: the Newton iteration stops at 1e-12 rad
# (`_MA_EPSILON` of ANISE, `conics` and `prop2b` of CSPICE stop at their own tolerances), and the
# state follows the anomaly with a factor that grows as 1/sqrt(|1 - e^2|) near the apsides.
cond_kepler(ecc) = 1e-12 * max(1, 1 / sqrt(abs(1 - ecc^2)))

# The semi-minor axis squared is compared with |a p|: b^2 = a^2 (1 - e^2) holds the rounding of
# e^2 (2 e^2 eps) divided by 1 - e^2, twice (the square, and the two sides).
cond_semiminor(o) = eps() * (1 + 4 * o.ecc^2 / abs(1 - o.ecc^2))

cond_shape(o, power=1) = eps() * (1 + 1 / abs(1 - o.ecc^power))
# Conditioning of the angles of an orbit (the model of the anise group, `cond_scale` of
# test/anise/runtests.jl). An element that is an arc cosine of a ratio of vectors has the rounding
# error eps of the ratio amplified by 1 / |sin angle|; the vectors carry a relative error of their
# own: eps / ecc for the eccentricity vector (ecc <= 1; a hyperbola has no small ecc in the
# samples) and eps / sin(inc) for the node vector. So the difference of two codes in
#   inc:  eps (1 + 1/sin inc)                      (acos of h_z / |h|)
#   raan: eps (1 + 1/(sin inc |sin raan|))         (acos of the node vector)
#   aop:  eps (1 + (1/ecc + 1/sin inc) / |sin aop|)
#   ta:   eps (1 + 1/(ecc |sin ta|))
#   ma:   the ta scale times max(1, (1 - e^2)^(3/2) / (1 + e cos ta)^2), the slope dM/dta
#   aol:  the larger of the aop and ta scales (aol = aop + ta)
#   state: the larger of the raan, aop and ta scales (a state rebuilt from its elements)
# is measured in units of that scale; the limits below are multiples of it (the `measured` lines
# print the largest multiple). A sample with an amplifier beyond 1e8 is bounded at 1e8.
_cs(x) = 1 / max(abs(sin(x)), 1e-8)
function cond_angle(o, which::Symbol)
    ae = 1 / min(o.ecc, 1.0)
    ai = 1 / max(sin(o.inc), 1e-8)
    which === :inc && return eps() * (1 + ai)
    which === :raan && return eps() * (1 + ai * _cs(o.raan))
    f_aop = 1 + (ae + ai) * _cs(o.aop)
    f_ta = 1 + ae * _cs(o.ta)
    f_raan = 1 + ai * _cs(o.raan)
    # ma = E - e sin E follows ta with dM/dnu = (1 - e^2)^(3/2) / (1 + e cos ta)^2, up to 12 at
    # e = 0.95 near apoapsis
    which === :ma && return eps() * f_ta * max(1, abs(1 - o.ecc^2)^1.5 / (1 + o.ecc * cos(o.ta))^2)
    return eps() * (which === :aop ? f_aop : which === :ta ? f_ta : which === :aol ? max(f_aop, f_ta) :
        max(f_aop, f_ta, f_raan))
end

# Tolerances against `oscltx`, in units of the conditioning models above (`cond_shape` for rp,
# ecc, sma_period; `cond_semiminor`; `cond_angle` for the angles). A limit in absolute terms
# failed at 5000 samples, because the worst of 200 samples has a condition number of 1e3 and
# that of 5000 one of 1e5. Largest multiples in the exhaustive tier (5000 orbits, seeds SEED
# and 1 to 5): rp 9.4, ecc 19, inc 1.9, raan 2.0, aop 3.3, ta 6.0, ma 5.9, sma_period 7.2, aol
# 3.0, semi_minor 8.1. Each limit is 4 to 5 times the largest (ecc 80, ma 25, semi_minor 35).
const TOL_OSC = Dict(:rp => 40, :ecc => 80, :inc => 8, :raan => 8,
    :aop => 14, :ta => 25, :ma => 25, :sma_period => 30)
# Identities between elements of one orbit, 200 orbits: measured c3 1.6e-16, mean_motion
# 2.8e-16, ta_dot 2.4e-16, ra 8.9e-16, dec 1.6e-15, vdec 4.4e-16, ha 8.9e-16 (also 5000 orbits:
# 4.4e-15 at most): rounding of a few ulps, limit 1e-13. The semi-minor axis is the square of
# a difference of ulps of sqrt(1-e^2), which `cond_semiminor` models (TOL_SEMIMINOR), and
# aol = aop + ta inherits the conditioning of aop and ta (TOL_AOL, units of cond_angle).
const TOL_ID = 1e-13
const TOL_AOL = 12
const TOL_SEMIMINOR = 35     # units of cond_semiminor
# B-plane of the GMAT case: the expected values are GMAT's (as in ANISE's test) and the
# limits are those of ANISE (1 km, 0.01 deg). Measured: 0.65 km at most (br), 2e-5 rad for
# the angle. The formula agrees with an independent B = b (S x h) form to 1e-10 km; a change
# of ephemeris (DE430 for DE440s) or of the epoch scale (A1 - TAI) moves the result by
# 0.013 km at most, and the rounding of the input values by 0.02 km at most. The rest of
# the difference lies in the reference values of GMAT, so the limits cannot be tightened.
const TOL_BPLANE = Dict(:bplane_br => 1.0, :bplane_bt => 1.0, :bplane_mag => 1.0,
    :bplane_angle => _d(0.01))
# Identities on random hyperbolas (measured: mag 2.2e-16, turn 4.4e-16, periapsis 1.2e-15;
# each is a few roundings of the same quantities), limit 1e-13 (TOL_ID); |b| against the
# semi-minor axis is the same number up to the rounding of `sqrt(1 - 1/e^2)` squared.
const TOL_BMAG = 1e-12      # measured 4.4e-16; |b| and the semi-minor axis differ by the rounding of sqrt(1 - 1/e^2)^2
const TOL_VINF = 1e-13
const TOL_VINF_RP = 1e-13
# ANISE cases: measured turn angle 2.2e-16 rad, periapsis 5.5e-12 km (the radius divides
# by cos(rho) near 1, so the relative error of 1e-16 of ~1e4 km shows as 1e-12 km).
const TOL_VINF_KM = 1e-10
# The sma printed by ANISE has 12 digits at 2e4 km: measured 4.0e-11 km (ulp 3.6e-12 km).
const TOL_VINF_SMA = 1e-9
# b_plane_davis: ha 1.1e-11 rad measured; ANISE's 1e-9 deg.
const TOL_DAVIS_HA = _d(1e-9)

# Geodetic elements against the values of ANISE (Vallado's examples, which ANISE prints with
# 16 digits): measured latitude 2.2e-16 rad, longitude 1.1e-16 rad, height 2.1e-12 km (a value of
# 5e3 km has an ulp of 9e-13), altitudes 4.2e-12 km; limits at about five times these.
const TOL_GEODETIC_ANG = 1e-15
const TOL_GEODETIC_KM = 2e-11
# Against `recgeo` (iteration of CSPICE), 900 points on three bodies at 1 to 4 equatorial radii:
# measured latitude 3.9e-16 rad, longitude 4.4e-16 rad, height 1.1e-15 equatorial radii, rounding
# of the two algorithms; limits at about five times these.
const TOL_RECGEO = (2e-15, 2e-15, 5e-15)
# The equinoctial identities are the same expressions as the elements (measured 0.0, the limit
# is TOL_ID). The conversion back to Keplerian elements (the sine and arc tangent of half the
# inclination and of p, q) has the rounding of the sum raan + aop + ma of three angles (up to
# 6 pi, ulp 3.6e-15 = 16 eps) and the acos of the inclination (1/sin inc): measured 1.9 at most
# in units of eps (16 + 1/sin inc) in 5000 orbits, six seeds; limit 8.
const TOL_EQ_BACK = 8
const _KEP6 = (:sma, :ecc, :inc, :raan, :aop, :ta)
# Building, editing and propagating (the `measured` lines of the testsets print the
# values; the limit is about five times the largest measured, 200 random orbits, seeds
# SEED + 51 to 53).
# keplerian_orbit against `conics`, relative, in units of `cond_kepler` (the Newton iteration of
# the Kepler equation stops at 1e-12 rad in ANISE, `conics` and `prop2b` stop at their own
# tolerance, and the state follows the anomaly with 1/sqrt(|1 - e^2|) at most): position 3.7 and
# velocity 1.9 at most (5000 orbits, six seeds; 1.2e-11 absolute at e = 0.949, where a fixed
# limit of 6e-12 failed); limit 25. Round trip of the state through the elements, the apsides
# or the true anomaly, in units of `cond_state`: 10.8 at most (limit 45); the Kepler iteration
# is not in that path (`o` is the reference and the rebuilt states start from its true
# anomaly). The elements read back: sma (relative) and ecc 1.0e-14 (limit 4e-14), the angles
# 9.2 in units of cond_angle (limit 40).
const TOL_BUILD = (25, 25, 45, 4e-14, 40)
# with_element: the asked sma (relative) or ecc (absolute: ecc = |e vector|, a few eps) in units
# of cond_shape of the new orbit: 25.5 at most (limit 100); the asked angle in units of cond_angle
# of the new orbit: 6.6 (limit 28); the other sma (relative) and ecc (absolute) 3.3e-14 (limit
# 1.4e-13); the other angles in units of the larger cond_angle of the old and the new orbit: 33
# (limit 140; a raan, aop or ta of a rebuilt state is an arc cosine). Six seeds, 5000 orbits.
const TOL_WITH = (100, 28, 1.4e-13, 140)
# verif_with_init of ANISE (loops of 900 or 360 values), largest 1.5e-11 (km or rad): the
# rounding of values of 8e3 to 1e4 km (ulp 1e-12 to 2e-12 each, summed over the conversion).
const TOL_INIT = 8e-11
# verif_orbit_at_epoch of ANISE: sma 3.6e-12 km (ANISE 1e-10), ecc 4.1e-16 (1e-7 deg), inc 0,
# raan 0, aop 3.3e-8 rad = 1.9e-6 deg (ANISE 3e-6 deg): the aop of an orbit of ecc 1e-8 is
# the angle of a vector of length 1e-8 whose rounding error is 1e-16.
const TOL_ATEPOCH = (2e-11, 1e-15, 1e-15, 1e-15, 5e-8)
# two_body against `prop2b`, relative, up to +-10 periods, in units of `cond_kepler`: position 8.1
# and velocity 5.5 at most (5000 orbits, six seeds; limit 25); sma 1.0e-14 (limit 3e-14, the
# mean motion divides by sma^1.5, a few ulps). The rounding of the mean longitude, a few ulps
# of up to 60 rad (1e-14), is 100 times smaller than the iteration tolerance.
const TOL_TWOBODY = (25, 25, 3e-14)
# A dual epoch: the derivative of the state is (v, -gm r/r^3), relative: 9e-15 at most over
# 1000 orbits (the derivative of the Kepler solver is exact; it was 3e-12 with the old lag).
const TOL_DUALDT = 4e-14
# The derivative of the true anomaly from the Kepler solver, against the closed forms, relative:
# d nu/dM 1.5e-14 (rounding), d nu/de (ellipses) 4.7e-12: the solver stops at 1e-12 rad, and the
# derivative is that of the implicit function at the converged E, not at the root. Limits 4 times
# those (seeds SEED and 1 to 4, 50000 samples).
const TOL_KEPLER_D = (6e-14, 2e-11)
# Local frames: the largest |R'R - I| and |det R - 1| are 8.9e-16 and 6.7e-16 in 5000 orbits (a
# few roundings of unit vectors, limit 3e-15). The rate of the frame of `two_body` by ForwardDiff
# against the exact rate at the state that `two_body` returns (the state of `o` differs from it
# by the rounding of the element conversions, up to 1e-12): 97 at most in units of eps |h|/r^2,
# the angular rate of the frame (limit 400). Before the Newton step of the Kepler solver
# (src/analysis/brouwer.jl) the difference was 3e-10 |h|/r^2, 1e6 eps.
const TOL_LF = (3e-15, 3e-15, 400)

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
        # (with radii: the altitudes, geodetic and Brouwer elements need them, see below)
        oe = Orbit(s, 0.0, 399, 1, 398600.4415, VALLADO_RADII)
        @test all(n -> (y = outcome(() -> element(oe, n)); y isa Float64 || y === :degeneratecase),
            elements())
        @test all(n -> outcome(() -> element(o, n)) in (:badarg, outcome(() -> element(oe, n))),
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
        gm = 398600.4415
        for orb in (Orbit(s, 0.0, 399, 1, gm, VALLADO_RADII),
                Orbit(s, ForwardDiff.Dual(1.0, 1.0), 399, 1, gm, VALLADO_RADII),
                Orbit(sh, 0.0, 399, 1, gm, VALLADO_RADII),
                Orbit(sh, ForwardDiff.Dual(1.0, 1.0), 399, 1, gm, VALLADO_RADII)),
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
            @test note(:rp, abs(o.periapsis - el[1]) / el[1] / cond_shape(o)) < TOL_OSC[:rp]
            @test note(:ecc, abs(o.ecc - el[2]) / cond_shape(o)) < TOL_OSC[:ecc]
            @test note(:inc, angdiff(o.inc, el[3]) / cond_angle(o, :inc)) < TOL_OSC[:inc]
            @test note(:raan, angdiff(o.raan, el[4]) / cond_angle(o, :raan)) < TOL_OSC[:raan]
            @test note(:aop, angdiff(o.aop, el[5]) / cond_angle(o, :aop)) < TOL_OSC[:aop]
            @test note(:ta, angdiff(o.ta, el[9]) / cond_angle(o, :ta)) < TOL_OSC[:ta]
            if o.ecc < 1
                @test note(:ma, angdiff(o.ma, el[6]) / cond_angle(o, :ma)) < TOL_OSC[:ma]
                @test note(:sma_period, max(abs(o.sma - el[10]) / el[10],
                    abs(o.period - el[11]) / el[11]) / cond_shape(o)) < TOL_OSC[:sma_period]
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
            @test note(:aol, angdiff(o.aol, o.aop + o.ta) / cond_angle(o, :aol)) < TOL_AOL
            @test note(:c3, abs(o.c3 - 2 * o.energy) / abs(o.c3)) < TOL_ID
            @test note(:semi_minor, abs(o.semi_minor_axis^2 - abs(o.sma * o.semi_parameter)) /
                o.semi_minor_axis^2 / cond_semiminor(o)) < TOL_SEMIMINOR
            _, ra, dec = recrad(o.position)
            @test note(:ra, angdiff(o.right_ascension, ra)) < TOL_ID
            @test note(:dec, abs(o.declination - dec)) < TOL_ID
            @test note(:vdec, abs(o.velocity_declination - recrad(o.velocity)[3])) < TOL_ID
        end
        measured("  oscltx and identities, largest differences: " *
            join((@sprintf("%s %.1e", k, worst[k]) for k in names), ", ") *
            " (rp, ecc, sma_period, semi_minor: units of cond_shape; inc, raan, aop, ta, ma, aol: cond_angle)")
    end

    @testset "Keplerian elements: degenerate cases" begin
        gm = 398600.4415
        orb(v) = Orbit(SVector{6,Float64}(v), 0.0, 399, 1, gm, VALLADO_RADII)
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
        measured(@sprintf("  bplane pro/retro: br %.3e, %.3e km; |b| - semi-minor %.3e km; retro angle %.17g",
            pro.bplane_br, retro.bplane_br, pro.bplane_mag - pro.semi_minor_axis, retro.bplane_angle))
        # Measured: br is exactly 0 in both (in the plane of the orbit), |b| equals the semi-minor
        # axis exactly; the limits of 1e-12 km and 1e-12 relative are roundings of values of 1e4 km.
        @test pro.bplane_bt > 0 && abs(pro.bplane_br) < 1e-12
        @test retro.bplane_bt < 0 && abs(retro.bplane_br) < 1e-12
        @test abs(pro.bplane_mag - pro.semi_minor_axis) < 1e-12 * pro.semi_minor_axis
        # br = +0 gives atan(0, bt<0) = π exactly, which [-π, π) maps to -π
        @test retro.bplane_angle === -Float64(π)
        @test pro.bplane_angle === 0.0
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
        names = (:mag, :angle, :turn, :periapsis, :bmag)
        worst = Dict(k => 0.0 for k in names)
        note(k, d) = (worst[k] = max(worst[k], d); d)
        nhyp = ninbound = 0
        for _ in 1:nsample(200, 5000)
            o = random_orbit(rng)
            o.ecc > 1 || continue
            nhyp += 1
            @test note(:mag, abs(o.bplane_mag - hypot(o.bplane_bt, o.bplane_br)) / o.bplane_mag) < TOL_ID
            @test note(:angle, angdiff(o.bplane_angle, atan(o.bplane_br, o.bplane_bt))) < TOL_ID
            @test -Float64(π) <= o.bplane_angle < Float64(π)
            @test note(:bmag, abs(o.bplane_mag - o.semi_minor_axis) / o.semi_minor_axis) < TOL_BMAG
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
            measured(@sprintf("  vinf %-24s ecc %.3e, sma %.3e km", c.name, abs(o.ecc - c.ecc), abs(o.sma - c.sma)))
            @test abs(o.ecc - c.ecc) < 1e-12 && abs(o.sma - c.sma) < TOL_VINF_SMA
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
    @testset "Geodetic and altitude elements: Vallado ($(c.name))" for c in VALLADO_CASES
        o = Orbit(vcat(c.pos, zeros(3)), 0.0, 399, 1, GMAT_GM, VALLADO_RADII)
        dl, dn, dh = abs(o.latitude - c.lat), angdiff(o.longitude, c.lon), abs(o.height - c.height)
        measured(@sprintf("  geodetic %-12s latitude %.3e rad, longitude %.3e rad, height %.3e km",
            c.name, dl, dn, dh))
        @test dl < TOL_GEODETIC_ANG && dn < TOL_GEODETIC_ANG && dh < TOL_GEODETIC_KM
        @test -π <= o.longitude < π
        if c.alts !== nothing
            d = (abs(o.sma_altitude - c.alts.sma), abs(o.apoapsis_altitude - c.alts.apoapsis),
                abs(o.periapsis_altitude - c.alts.periapsis))
            measured(@sprintf("  altitudes %-11s sma %.3e, apoapsis %.3e, periapsis %.3e km", c.name, d...))
            @test all(<(TOL_GEODETIC_KM), d)
            @test o.altitude == o.rmag - 6378.14
        end
    end

    @testset "Geodetic elements against recgeo" begin
        rng = Xoshiro(SEED + 41)
        bodies = ((6378.1366, 6378.1366, 6356.7519), (3396.19, 3396.19, 3376.2),
            (1000.0, 1000.0, 900.0))                       # the last is flattened far beyond any planet
        worst = zeros(3)
        for radii in bodies
            re = radii[1]
            f = (radii[1] - radii[3]) / radii[1]
            for _ in 1:nsample(300, 10000)
                u = randn(rng, 3)
                p = (1 + 3 * rand(rng)) * re * u / norm(u)     # 1 to 4 equatorial radii
                o = Orbit(vcat(p, 0.0, 0.0, 0.0), 0.0, 399, 1, GMAT_GM, collect(radii))
                lon, lat, alt = recgeo(p, re, f)
                d = (angdiff(o.latitude, lat), angdiff(o.longitude, lon), abs(o.height - alt) / re)
                worst .= max.(worst, d)
                @test d[1] < TOL_RECGEO[1] && d[2] < TOL_RECGEO[2] && d[3] < TOL_RECGEO[3]
            end
        end
        measured(@sprintf("  geodetic against recgeo, largest differences: latitude %.1e rad, longitude %.1e rad, height %.1e radii",
            worst...))
        # the pole (p = 0) gives no NaN
        pole = Orbit(SVector(0.0, 0.0, 6356.7519 + 100, 0, 0, 0), 0.0, 399, 1, GMAT_GM,
            [6378.1366, 6378.1366, 6356.7519])
        @test pole.latitude ≈ π / 2 && abs(pole.height - 100) < 1e-8 && pole.longitude == 0.0
    end

    @testset "Equinoctial elements" begin
        rng = Xoshiro(SEED + 42)
        worst = [0.0, 0.0]
        for _ in 1:nsample(200, 5000)
            o = random_orbit(rng)
            o.ecc < 1 || continue
            d = max(abs(o.equinoctial_h - o.ecc * sin(o.aop + o.raan)),
                abs(o.equinoctial_k - o.ecc * cos(o.aop + o.raan)),
                abs(o.equinoctial_p - sin(o.inc / 2) * sin(o.raan)),
                abs(o.equinoctial_q - sin(o.inc / 2) * cos(o.raan)),
                abs(o.equinoctial_lambda - (o.raan + o.aop + o.ma)))
            # ANISE equinoctial_to_keplerian returns the Keplerian elements
            sma, ecc, inc, raan, aop, ma = Almanac._equinoctial_to_kep(o.equinoctial_a, o.equinoctial_h,
                o.equinoctial_k, o.equinoctial_p, o.equinoctial_q, o.equinoctial_lambda)
            b = max(abs(ecc - o.ecc), abs(inc - o.inc), angdiff(raan, o.raan), angdiff(aop, o.aop),
                angdiff(ma, o.ma)) / (eps() * (16 + 1 / sin(o.inc)))
            worst .= max.(worst, (d, b))
            @test d < TOL_ID && o.equinoctial_a == o.sma && sma == o.sma && b < TOL_EQ_BACK
        end
        measured(@sprintf("  equinoctial identities %.1e, conversion back to Keplerian %.1e (units of eps (16 + 1/sin inc))", worst...))
        hyp = Orbit(gmat_kep_state(-20000.0, 1.5, _d(28.5), _d(10), _d(20), _d(30)), 0.0, 399, 1, GMAT_GM)
        @test hyp.equinoctial_h ≈ 1.5 * sin(hyp.aop + hyp.raan) && hyp.equinoctial_lambda ≈ hyp.raan + hyp.aop + hyp.ma
        # the errors are those of the Keplerian elements (a circular orbit has no ta, ma)
        circ = Orbit(SVector(7000.0, 0, 0, 0, 0, sqrt(GMAT_GM / 7000)), 0.0, 399, 1, GMAT_GM)
        @test throws_kind(() -> circ.equinoctial_lambda, :degeneratecase)
        @test isfinite(circ.equinoctial_h) && circ.equinoctial_a == circ.sma
    end

    @testset "Brouwer mean elements: ANISE round trips ($(c.name))" for c in BROUWER_CASES
        re, gm = VALLADO_RADII[1], GMAT_GM
        osc = Almanac._brouwer_mean_to_osc(c.mean..., re, Almanac.J2_EARTH)
        o = Orbit(Almanac._kep_ma_to_cart(osc..., gm), 0.0, 399, 1, gm, VALLADO_RADII)
        got = (o.brouwer_sma, o.brouwer_ecc, o.brouwer_inc, o.brouwer_raan, o.brouwer_aop, o.brouwer_ma)
        d = (abs(c.mean[1] - got[1]), abs(c.mean[2] - got[2]),
            (angdiff(c.mean[i], got[i]) for i in 3:6)...)
        measured(@sprintf("  brouwer %-14s sma %.2e km, ecc %.2e, inc %.2e, raan %.2e, aop %.2e, ma %.2e rad",
            c.name, d...))
        # ANISE: 4 km, 1e-2, 0.7 degrees; then the limits of this test
        @test d[1] < 4.0 && d[2] < 1e-2 && all(<(_d(0.7)), d[3:6])
        @test d[1] < c.tol[1] && d[2] < c.tol[2] && all(<(c.tol[3]), d[3:6])
    end

    @testset "Brouwer mean elements: validity" begin
        gm = GMAT_GM
        good = gmat_kep_state(7000.0, 0.01, _d(51.6), _d(30), _d(20), _d(10))
        o = Orbit(good, 0.0, 399, 1, gm, VALLADO_RADII)
        @test o.brouwer_sma isa Float64
        # no radii: the equatorial radius is unknown
        @test throws_kind(() -> Orbit(good, 0.0, 399, 1, gm).brouwer_sma, :badarg)
        # ANISE: Earth-centered only, periapsis above 3000 km, elliptic
        @test throws_kind(() -> Orbit(good, 0.0, 301, 1, gm, VALLADO_RADII).brouwer_sma, :degeneratecase)
        low = Orbit(gmat_kep_state(5000.0, 0.5, _d(30), 0.0, 0.0, 1.0), 0.0, 399, 1, gm, VALLADO_RADII)
        @test low.periapsis < 3000
        hyp = Orbit(gmat_kep_state(-20000.0, 1.5, _d(28.5), _d(10), _d(20), _d(30)), 0.0, 399, 1, gm, VALLADO_RADII)
        for n in (:brouwer_sma, :brouwer_ecc, :brouwer_inc, :brouwer_raan, :brouwer_aop, :brouwer_ma)
            @test throws_kind(() -> element(low, n), :degeneratecase)
            @test throws_kind(() -> element(hyp, n), :degeneratecase)
        end
        @test 0 <= o.brouwer_raan < 2π && 0 <= o.brouwer_aop < 2π && 0 <= o.brouwer_ma < 2π
        od = Orbit(good, ForwardDiff.Dual(0.0, 1.0), 399, 1, gm, VALLADO_RADII)
        @test od.brouwer_sma isa ForwardDiff.Dual && ForwardDiff.value(od.brouwer_sma) == o.brouwer_sma
    end

    @testset "Elements that need the radii of the center" begin
        ctx = neo_context([Kernels.LSK, Kernels.DE440S, Kernels.GM_TPC, Kernels.PCK11])
        o = Orbit(ctx, [7000.0, 0.0, 0.0, 0.0, 8.0, 0.0], 0.0, 3, "J2000")   # EMB: a GM, no radii
        @test all(isnan, o.radii) && o.sma isa Float64
        for n in (:altitude, :periapsis_altitude, :apoapsis_altitude, :sma_altitude, :latitude,
                :height)
            @test throws_kind(() -> element(o, n), :badarg)
        end
        # Brouwer: ANISE checks the center (Earth only) before it reads the shape
        @test throws_kind(() -> o.brouwer_sma, :degeneratecase)
        @test o.longitude == 0.0                    # ANISE longitude_deg needs no shape
        e = Orbit(ctx, [7000.0, 0.0, 0.0, 0.0, 8.0, 0.0], 0.0, 399, "J2000")
        @test e.radii == [6378.1366, 6378.1366, 6356.7519]
        @test e.altitude == e.rmag - 6378.1366 && e.sma_altitude == e.sma - 6378.1366
    end

    @testset "keplerian_orbit against conics and the round trip" begin
        rng = Xoshiro(SEED + 51)
        gm = 398600.4415
        worst = zeros(6)
        for _ in 1:nsample(200, 5000)
            ecc = rand(rng) < 0.5 ? 0.005 + 0.945 * rand(rng) : 1.05 + 3.95 * rand(rng)
            rp = 6500.0 + 43500.0 * rand(rng)
            inc = 0.02 + (π - 0.04) * rand(rng)
            raan, aop = 2π * rand(rng), 2π * rand(rng)
            m0 = ecc < 1 ? 2π * rand(rng) : 6 * rand(rng) - 3
            ref = conics([rp, ecc, inc, raan, aop, m0, 0.0, gm], 0.0)
            sma = rp / (1 - ecc)
            o = keplerian_orbit(0.0, 399, 1; sma, ecc, inc, raan, aop, ma=m0, gm)
            # the same orbit by the radii of the apsides and by the true anomaly, and the round trip
            ra = sma * (1 + ecc)
            o2 = ecc < 1 ? keplerian_orbit(0.0, 399, 1; periapsis=rp, apoapsis=ra, inc, raan, aop,
                ta=o.ta, gm) : o
            o3 = keplerian_orbit(0.0, 399, 1; sma, ecc, inc, raan, aop, ta=o.ta, gm)
            dr = norm(o.position - ref[1:3]) / norm(ref[1:3])
            dv = norm(o.velocity - ref[4:6]) / norm(ref[4:6])
            # round trip through the elements of the state
            k = keplerian(o)
            o4 = keplerian_orbit(0.0, 399, 1; k..., gm)
            ck = cond_kepler(ecc)
            dk = maximum(x -> norm(x.state - o.state) / norm(o.state), (o2, o3, o4)) / cond_state(o)
            worst[1:3] .= max.(worst[1:3], (dr / ck, dv / ck, dk))
            @test dr / ck < TOL_BUILD[1] && dv / ck < TOL_BUILD[2] && dk < TOL_BUILD[3]
            @test o.et === 0.0 && o.gm == gm && o.center === Int32(399) && o.frame === Int32(1)
            @test all(isnan, o.radii)
            kk = keplerian(o4)
            de = max(abs(kk.sma - k.sma) / abs(k.sma), abs(kk.ecc - k.ecc))
            da = max(angdiff(kk.inc, k.inc) / cond_angle(o, :inc),
                angdiff(kk.raan, k.raan) / cond_angle(o, :raan),
                angdiff(kk.aop, k.aop) / cond_angle(o, :aop), angdiff(kk.ta, k.ta) / cond_angle(o, :ta))
            worst[4] = max(worst[4], de)
            worst[5] = max(worst[5], da)
            @test de < TOL_BUILD[4] && da < TOL_BUILD[5]
        end
        measured(@sprintf("  keplerian_orbit against conics, relative, in units of cond_kepler: position %.2e,
    velocity %.2e; round trip of the state %.2e (units of cond_state), of sma and ecc %.2e,
    of the angles %.2e (units of cond_angle)", worst[1:5]...))
    end

    @testset "keplerian_orbit: signs, errors and number types" begin
        gm = 398600.4415
        kw = (inc=0.5, raan=1.0, aop=2.0, gm=gm)
        o = keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=0.1, ta=0.3, kw...)
        # ANISE changes the sign of a negative ecc and of an sma of the wrong sign
        @test keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=-0.1, ta=0.3, kw...).state == o.state
        @test keplerian_orbit(0.0, 399, 1; sma=-8000.0, ecc=0.1, ta=0.3, kw...).state == o.state
        h = keplerian_orbit(0.0, 399, 1; sma=-20000.0, ecc=1.5, ta=0.3, kw...)
        @test keplerian_orbit(0.0, 399, 1; sma=20000.0, ecc=1.5, ta=0.3, kw...).state == h.state
        # periapsis above apoapsis: the eccentricity is negative, and its sign changes
        @test keplerian_orbit(0.0, 399, 1; periapsis=9000.0, apoapsis=7000.0, ta=0.3, kw...).ecc ≈ 0.125 rtol=1e-12
        # argument errors
        for bad in (
                () -> keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=0.1, ta=0.3, ma=0.3, kw...),
                () -> keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=0.1, kw...),
                () -> keplerian_orbit(0.0, 399, 1; sma=8000.0, ta=0.3, kw...),
                () -> keplerian_orbit(0.0, 399, 1; ecc=0.1, ta=0.3, kw...),
                () -> keplerian_orbit(0.0, 399, 1; ta=0.3, kw...),
                () -> keplerian_orbit(0.0, 399, 1; periapsis=7000.0, ta=0.3, kw...),
                () -> keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=0.1, periapsis=7000.0,
                    apoapsis=9000.0, ta=0.3, kw...),
                () -> keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=0.1, apoapsis=9000.0, ta=0.3, kw...))
            @test throws_kind(bad, :badarg)
        end
        # the failures of ANISE: parabola, anomaly beyond the asymptotes, radii of apsis
        @test throws_kind(() -> keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=1.0, ta=0.3, kw...), :degeneratecase)
        @test throws_kind(() -> keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=1 + 5e-12, ta=0.3, kw...), :degeneratecase)
        @test keplerian_orbit(0.0, 399, 1; sma=-8000.0, ecc=1 + 2e-11, ta=0.3, kw...) isa Orbit
        @test throws_kind(() -> keplerian_orbit(0.0, 399, 1; sma=-8000.0, ecc=2.0, ta=2.5, kw...), :degeneratecase)
        @test keplerian_orbit(0.0, 399, 1; sma=-8000.0, ecc=2.0, ta=-2.0, kw...) isa Orbit    # inbound leg
        @test throws_kind(() -> keplerian_orbit(0.0, 399, 1; periapsis=0.0, apoapsis=9000.0, ta=0.3, kw...), :degeneratecase)
        @test throws_kind(() -> keplerian_orbit(0.0, 399, 1; periapsis=7000.0, apoapsis=-1.0, ta=0.3, kw...), :degeneratecase)
        @test throws_kind(() -> keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=1.0, ma=0.3, kw...), :degeneratecase)
        # no kernels with gm and integer names; the kernel forms read GM and the radii
        ctx0 = SpiceContext()
        @test keplerian_orbit(ctx0, 0.0, 399, 1; sma=8000.0, ecc=0.1, ta=0.3, kw...).state == o.state
        @test throws_kind(() -> keplerian_orbit(ctx0, 0.0, 399, 1; sma=8000.0, ecc=0.1, ta=0.3,
            inc=0.5, raan=1.0, aop=2.0), :kernelvarnotfound)
        ctx = neo_context([Kernels.LSK, Kernels.DE440S, Kernels.GM_TPC, Kernels.PCK11])
        ok = keplerian_orbit(ctx, 5.0, "EARTH", "J2000"; sma=8000.0, ecc=0.1, ta=0.3, inc=0.5, raan=1.0, aop=2.0)
        @test ok.gm == bodvcd(ctx, 399, "GM")[1] && ok.radii == bodvcd(ctx, 399, "RADII") && ok.et == 5.0
        @test ok == Orbit(ctx, keplerian_orbit(5.0, 399, 1; sma=8000.0, ecc=0.1, ta=0.3, inc=0.5, raan=1.0,
            aop=2.0, gm=ok.gm).state, 5.0, "EARTH", "J2000")
        ks = Almanac.kernelset(ctx)
        @test keplerian_orbit(ks, 5.0, 399, 1; sma=8000.0, ecc=0.1, ta=0.3, inc=0.5, raan=1.0, aop=2.0) == ok
        # dual numbers flow through the elements and the epoch
        d = ForwardDiff.derivative(x -> keplerian_orbit(0.0, 399, 1; sma=x, ecc=0.1, ta=0.3, kw...).sma, 8000.0)
        @test d ≈ 1 atol=1e-9
        dm = keplerian_orbit(ForwardDiff.Dual(0.0, 1.0), 399, 1; sma=8000.0, ecc=0.1, ma=0.3, kw...)
        @test dm isa Orbit{<:ForwardDiff.Dual} && ForwardDiff.partials(dm.et)[1] == 1.0
    end

    @testset "with_element and with_apsides" begin
        rng = Xoshiro(SEED + 52)
        worst = zeros(4)
        for _ in 1:nsample(100, 2000), name in (:sma, :ecc, :inc, :raan, :aop, :ta)
            o = random_orbit(rng)
            hyp = o.ecc > 1
            lim = hyp ? π - acos(1 / o.ecc) : π
            new = if name === :sma
                o.sma * (1 + 0.3 * rand(rng))
            elseif name === :ecc
                # a hyperbola keeps its true anomaly only inside the asymptotes of the new
                # eccentricity (1 + ecc cos(ta) > 0, as in `try_keplerian` of ANISE, which
                # refuses the others): the new value stays below -1/cos(ta), at 90% of the gap
                hyp ? o.ecc + 0.9 * (min(1.1 * o.ecc, cos(o.ta) < 0 ? -1 / cos(o.ta) : Inf) - o.ecc) * rand(rng) :
                0.005 + 0.945 * rand(rng)
            elseif name === :inc
                0.02 + (π - 0.04) * rand(rng)
            elseif name === :ta
                0.9 * lim * (2rand(rng) - 1)
            else
                2π * rand(rng)
            end
            p = with_element(o, name, new)
            before, after = keplerian(o), keplerian(p)
            # the asked element, and the other five unchanged (angles modulo a turn)
            got = getproperty(after, name)
            want = getproperty(before, name)
            # sma and ecc in relative and absolute difference, the angles in units of cond_angle
            # (of the new orbit for the asked element, the larger of the two for the others)
            d = name === :sma ? abs(got - new) / abs(new) / cond_shape(p) :
                name === :ecc ? abs(got - new) / cond_shape(p) : angdiff(got, new) / cond_angle(p, name)
            ia = name in (:sma, :ecc) ? 1 : 2
            worst[ia] = max(worst[ia], d)
            @test d < TOL_WITH[ia]
            for n in _KEP6
                n === name && continue
                x = n === :sma ? abs(after.sma - before.sma) / abs(before.sma) :
                    n === :ecc ? abs(after.ecc - before.ecc) :
                    angdiff(getproperty(after, n), getproperty(before, n)) /
                    max(cond_angle(o, n), cond_angle(p, n))
                j = n in (:sma, :ecc) ? 3 : 4
                worst[j] = max(worst[j], x)
                @test x < TOL_WITH[j]
            end
            @test isequal((p.et, p.center, p.frame, p.gm, p.radii), (o.et, o.center, o.frame, o.gm, o.radii))
        end
        measured(@sprintf("  with_element: asked sma (relative) or ecc (absolute) %.2e (units of cond_shape), angle %.2e;
    the others sma (relative) or ecc (absolute) %.2e, angles %.2e (units of cond_angle)", worst...))
        # ANISE verif_with_init (tests/astro/orbit.rs): the same loops, in rad
        kep = keplerian_orbit(0.0, 399, 1; sma=8191.93, ecc=0.0245, inc=_d(12.85), raan=_d(306.614),
            aop=_d(314.19), ta=_d(99.8877), gm=GMAT_GM)
        w = zeros(8)
        for incr in 100:999
            w[1] = max(w[1], abs(with_element(kep, :sma, kep.sma + incr).sma - (kep.sma + incr)))
        end
        for incr in 0:99
            e = kep.ecc + incr / 100
            w[2] = max(w[2], abs(with_element(kep, :ecc, e).ecc - e))
        end
        for incr in 0:359
            x = _d(incr)
            w[3] = max(w[3], angdiff(with_element(kep, :aop, mod(kep.aop + x, 2π)).aop, kep.aop + x))
            w[4] = max(w[4], angdiff(with_element(kep, :raan, mod(kep.raan + x, 2π)).raan, kep.raan + x))
            w[5] = max(w[5], angdiff(with_element(kep, :ta, mod(kep.ta + x, 2π)).ta, kep.ta + x))
            inc = abs(rem2pi(kep.inc + x, RoundNearest))
            w[6] = max(w[6], abs(with_element(kep, :inc, kep.inc + x).inc - inc))
        end
        for delta in 100:999
            ra, rp = kep.apoapsis + delta, kep.periapsis - delta
            q = with_apsides(kep, rp, ra)
            w[7] = max(w[7], abs(q.apoapsis - ra))
            w[8] = max(w[8], abs(q.periapsis - rp))
        end
        measured(@sprintf("  verif_with_init: sma %.1e km, ecc %.1e, aop %.1e, raan %.1e, ta %.1e, inc %.1e rad, ra %.1e, rp %.1e km", w...))
        # ANISE compares with 1e-10 (km or deg); the limit is that of the largest measured
        @test all(<(TOL_INIT), w)
        # the others stay when the apsides change
        q = with_apsides(kep, 7000.0, 9000.0)
        @test q.periapsis ≈ 7000.0 && q.apoapsis ≈ 9000.0 && angdiff(q.inc, kep.inc) < 1e-12 &&
              angdiff(q.raan, kep.raan) < 1e-12 && angdiff(q.aop, kep.aop) < 1e-11 && angdiff(q.ta, kep.ta) < 1e-11
        # errors: argument names, radii, and the circular orbit (no true anomaly)
        @test throws_kind(() -> with_element(kep, :foo, 1.0), :badarg)
        @test throws_kind(() -> with_element(kep, :ma, 1.0), :badarg)
        @test throws_kind(() -> with_apsides(kep, 0.0, 9000.0), :degeneratecase)
        @test throws_kind(() -> with_apsides(kep, 7000.0, 0.0), :degeneratecase)
        circ = Orbit(SVector(7000.0, 0, 0, 0, 0, sqrt(GMAT_GM / 7000)), 0.0, 399, 1, GMAT_GM)
        for n in (:sma, :ecc, :inc, :raan, :aop)
            @test throws_kind(() -> with_element(circ, n, 1.0), :degeneratecase)
        end
        @test throws_kind(() -> with_apsides(circ, 7000.0, 9000.0), :degeneratecase)
        @test with_element(circ, :ta, 1.0) isa Orbit      # ANISE set_ta_deg does not read ta
        @test with_element(kep, :sma, kep.sma) isa Orbit{Float64}
        dual = Orbit(kep.state, ForwardDiff.Dual(0.0, 1.0), 399, 1, GMAT_GM)
        @test with_element(dual, :aop, 1.0) isa Orbit{<:ForwardDiff.Dual}
    end

    @testset "derivative of the Kepler solver" begin
        # d nu / d M = (1 + e cos nu)^2 / |1 - e^2|^(3/2) (ellipse and hyperbola) and, on an ellipse,
        # d nu / d e = sin nu (2 + e cos nu) / (1 - e^2) at fixed M. Before the Newton step in the
        # number type (src/analysis/brouwer.jl), the dual part lagged the iteration by one step:
        # errors up to 5e-8. The samples include M within 1e-9 of pi, where the solver
        # returns E itself as the true anomaly (|E - pi| < 1e-8) and the derivative is
        # sqrt((1 - e)/(1 + e)) dE.
        rng = Xoshiro(SEED + 56)
        worst = zeros(3)
        for i in 1:nsample(2000, 50000)
            hyp = rand(rng) < 0.4
            ecc = hyp ? 1.05 + 3.95 * rand(rng) : 0.001 + 0.949 * rand(rng)
            m = hyp ? 6 * rand(rng) - 3 : (i % 5 == 0 ? π + 2e-9 * (rand(rng) - 0.5) : 6π * (2rand(rng) - 1))
            ta = Almanac._ma_to_ta(promote(m, ecc)...)
            @test ta isa Float64
            dta = ForwardDiff.derivative(x -> Almanac._ma_to_ta(promote(x, ecc)...), m)
            # in the branch |E - pi| < 1e-8 the solver returns E, which differs from nu by (1 - c) (E - pi)
            nu = !hyp && abs(ta - π) < 1e-8 ? π + sqrt((1 - ecc) / (1 + ecc)) * (ta - π) : ta
            want = (1 + ecc * cos(nu))^2 / abs(1 - ecc^2)^1.5
            e1 = abs(dta - want) / want
            worst[1] = max(worst[1], e1)
            @test e1 < TOL_KEPLER_D[1]
            if !hyp
                dte = ForwardDiff.derivative(x -> Almanac._ma_to_ta(promote(m, x)...), ecc)
                wante = sin(nu) * (2 + ecc * cos(nu)) / (1 - ecc^2)
                e2 = abs(dte - wante) / max(abs(wante), 1 / (1 - ecc^2))
                worst[2] = max(worst[2], e2)
                @test e2 < TOL_KEPLER_D[2]
            end
        end
        measured(@sprintf("  Kepler solver: d nu/dM relative error %.2e, d nu/de (ellipses) %.2e", worst[1:2]...))
    end

    @testset "two_body" begin
        rng = Xoshiro(SEED + 53)
        gm = 398600.4415
        worst = zeros(3)
        for _ in 1:nsample(200, 5000)
            o = random_orbit(rng)
            n = o.mean_motion
            dt = o.ecc < 1 ? 20 * (rand(rng) - 0.5) * 2π / n : (2rand(rng) - 1) * 3 / n
            ref = prop2b(gm, o.state, dt)
            t = two_body(o, o.et + dt)
            dr = norm(t.position - ref[1:3]) / norm(ref[1:3])
            dv = norm(t.velocity - ref[4:6]) / norm(ref[4:6])
            ds = abs(t.sma - o.sma) / abs(o.sma)
            # Kepler's equation is solved by a Newton iteration that stops at 1e-12 rad (`_MA_EPSILON`,
            # as in ANISE); the state follows the anomaly with the factor 1/sqrt(|1 - e^2|) at most
            # near the apsides of an eccentric orbit (rounding of the mean longitude, eps |n dt|
            # with |n dt| < 60, is 20 times smaller)
            cond = cond_kepler(o.ecc)
            worst .= max.(worst, (dr / cond, dv / cond, ds))
            @test dr / cond < TOL_TWOBODY[1] && dv / cond < TOL_TWOBODY[2] && ds < TOL_TWOBODY[3]
            @test t.et == o.et + dt && isequal((t.center, t.frame, t.gm, t.radii), (o.center, o.frame, o.gm, o.radii))
        end
        measured(@sprintf("  two_body against prop2b, relative: position %.2e, velocity %.2e (units of 1e-12 max(1, 1/sqrt(|1 - e^2|))); sma %.2e", worst...))
        # the same epoch: the orbit itself
        o = random_orbit(rng)
        @test two_body(o, o.et) === o
        # ANISE verif_orbit_at_epoch: half a period later only the anomaly has changed
        kep(args...) = keplerian_orbit(0.0, 399, 1; sma=8191.93, ecc=args[1], inc=_d(args[2]),
            raan=_d(306.614), aop=_d(314.19), ta=_d(args[3]), gm=GMAT_GM)
        rv, vv = 26559.381, 3.874002
        cases = (kep(1e-8, 12.85, 99.8877), kep(0.2, 12.85, -99.8877), kep(1e-8, 1e-5, -99.8877),
            (Orbit(SVector(a, b, 0, c, d, 0), 0.0, 399, 1, GMAT_GM) for (a, b, c, d) in
                ((rv, 0.0, 0.0, vv), (0.0, rv, -vv, 0.0), (-rv, 0.0, 0.0, -vv), (0.0, -rv, vv, 0.0)))...)
        w = zeros(5)
        for orb in cases
            f = two_body(orb, orb.et + 0.5 * orb.period)
            w .= max.(w, (abs(orb.sma - f.sma), abs(orb.ecc - f.ecc), angdiff(orb.inc, f.inc),
                angdiff(orb.raan, f.raan), angdiff(orb.aop, f.aop)))
        end
        measured(@sprintf("  verif_orbit_at_epoch: sma %.1e km, ecc %.1e, inc %.1e, raan %.1e, aop %.1e rad", w...))
        # ANISE: sma and raan 1e-10, ecc and inc 1e-7 deg, aop 3e-6 deg; the limits are of the measured
        @test w[1] < TOL_ATEPOCH[1] && w[2] < TOL_ATEPOCH[2] && w[3] < TOL_ATEPOCH[3] &&
              w[4] < TOL_ATEPOCH[4] && w[5] < TOL_ATEPOCH[5]
        # the true anomaly advanced as the mean anomaly: half a period of a circular-ish orbit
        c = cases[1]
        @test angdiff(two_body(c, c.et + c.period).ta, c.ta) < 1e-6
        # an exactly circular orbit has no mean anomaly: ANISE fails, so does this
        circ = Orbit(SVector(7000.0, 0, 0, 0, 0, sqrt(GMAT_GM / 7000)), 0.0, 399, 1, GMAT_GM)
        @test throws_kind(() -> two_body(circ, 100.0), :degeneratecase)
        @test two_body(circ, 0.0) === circ
        # a dual epoch: the derivative of the state is the velocity and the acceleration
        dd = zeros(2)
        for _ in 1:nsample(20, 1000)
            o = random_orbit(rng)
            t0 = 1000.0 * (2rand(rng) - 1)
            s = ForwardDiff.derivative(t -> Vector(two_body(o, t).state), t0)
            st = two_body(o, t0).state
            r = norm(st[1:3])
            d = (norm(s[1:3] - st[4:6]) / norm(st[4:6]), norm(s[4:6] + gm * st[1:3] / r^3) / (gm / r^2))
            dd .= max.(dd, d)
            @test all(<(TOL_DUALDT), d)
        end
        measured(@sprintf("  two_body with a dual epoch, relative: dr/dt - v %.2e, dv/dt + gm r/r^3 %.2e", dd...))
        o = random_orbit(rng)
        od = Orbit(o.state, ForwardDiff.Dual(0.0, 1.0), 399, 1, gm)
        @test two_body(od, 50.0) isa Orbit{<:ForwardDiff.Dual}
        @test two_body(o, ForwardDiff.Dual(o.et, 1.0)).state isa SVector{6,<:ForwardDiff.Dual}
    end

    @testset "local frames and relative states" begin
        rng = Xoshiro(SEED + 54)
        kinds = (:ric, :vnc, :rcn, :tnw)
        w = zeros(3)
        for _ in 1:nsample(6, 200)
            o = random_orbit(rng)
            r, v, h = o.position, o.velocity, o.hvec
            rh, vh, hh = r / norm(r), v / norm(v), h / norm(h)
            firsts = (:ric => rh, :vnc => vh, :rcn => rh, :tnw => vh)
            for (k, c1) in firsts
                R = local_frame(o, k)
                @test R isa SMatrix{3,3,Float64,9}
                w[1] = max(w[1], maximum(abs, R' * R - I))
                w[2] = max(w[2], abs(det(R) - 1))
                @test R[:, 1] ≈ c1 atol = 1e-15
                @test R[:, 3] ≈ (k in (:ric, :rcn, :tnw) ? hh : cross(vh, hh)) atol = 1e-15
                # the exact rate against the derivative of the frame of the propagated orbit
                S = local_frame_state(o, k)
                @test S isa SMatrix{6,6,Float64,36}
                @test S[1:3, 1:3] == R && S[4:6, 4:6] == R && all(iszero, S[1:3, 4:6])
                D = ForwardDiff.derivative(t -> local_frame(two_body(o, t), k), o.et)
                # two_body returns the state of `o` only up to the rounding of its element
                # conversions (cond_state, up to 1e-12), so the rate is compared at the state that
                # it returns: the derivative of the frame is the exact rate of that state
                tq = two_body(o, ForwardDiff.Dual(o.et, 1.0))
                q = Orbit(map(ForwardDiff.value, tq.state), o.et, 399, 1, o.gm)
                Sq = local_frame_state(q, k)
                w[3] = max(w[3], maximum(abs, Sq[4:6, 1:3] - D) / (norm(h) / norm(r)^2) / eps())
            end
        end
        measured(@sprintf("  local frames: |R'R - I| %.1e, |det R - 1| %.1e, rate against ForwardDiff %.1f (units of eps |h|/r^2, the angular rate of the frame)", w...))
        @test w[1] < TOL_LF[1] && w[2] < TOL_LF[2] && w[3] < TOL_LF[3]
        o = random_orbit(rng)
        @test local_frame(o, :rcn) == local_frame(o, :ric)
        # the axes as ANISE lists them
        r, v, h = o.position, o.velocity, o.hvec
        rh, vh, hh = normalize(r), normalize(v), normalize(h)
        @test local_frame(o, :ric) ≈ [rh cross(hh, rh) hh] atol = 1e-15
        @test local_frame(o, :vnc) ≈ [vh hh cross(vh, hh)] atol = 1e-15
        @test local_frame(o, :tnw) ≈ [vh cross(hh, vh) hh] atol = 1e-15
        # relative states
        for k in kinds
            @test relative_state(o, o, k) == zeros(6)
        end
        # (the position has the magnitude of 1e4 km, so the offset is exact to about 1e-11 km)
        radial = Orbit(o.state + SVector(o.position / norm(o.position)..., 0, 0, 0), o.et, 399, 1, o.gm)
        @test relative_state(radial, o, :ric)[1:3] ≈ SVector(1, 0, 0) atol = 1e-10
        # the inverse of the 6 by 6 matrix: it maps the relative state back to the difference
        # (ANISE ric_difference_unchecked), also for a difference that has a velocity part
        for _ in 1:nsample(5, 50), k in kinds
            c = random_orbit(rng)
            q = Orbit(c.state + SVector(1.5, -2.0, 0.7, 2e-3, -1e-3, 3e-3), c.et, 399, 1, c.gm)
            rs = relative_state(q, c, k)
            @test local_frame_state(c, k) * rs ≈ q.state - c.state rtol = 1e-12
        end
        # a pure position offset has a relative velocity: Rdot' dr, which is not zero
        S = local_frame_state(o, :ric)
        R, Rd = S[1:3, 1:3], S[4:6, 1:3]
        # (1 km along the radius: the frame turns at the angular rate |h|/r^2, so the relative
        # velocity is -|h|/r^2 along the in-track axis)
        dr = o.position / norm(o.position)
        pure = Orbit(o.state + SVector(dr..., 0, 0, 0), o.et, 399, 1, o.gm)
        rs = relative_state(pure, o, :ric)
        @test rs[4:6] ≈ Rd' * (pure.position - o.position) rtol = 1e-9
        @test rs[5] ≈ -norm(o.hvec) / norm(o.position)^2 rtol = 1e-9
        # dual numbers pass through
        od = Orbit(o.state, ForwardDiff.Dual(o.et, 1.0), 399, 1, o.gm)
        @test local_frame(od, :tnw) isa SMatrix{3,3,<:ForwardDiff.Dual}
        @test relative_state(pure, od, :vnc) isa SVector{6,<:ForwardDiff.Dual}
        # errors
        @test throws_kind(() -> local_frame(o, :xyz), :badarg)
        @test throws_kind(() -> local_frame_state(o, :xyz), :badarg)
        @test throws_kind(() -> relative_state(o, o, :xyz), :badarg)
        @test throws_kind(() -> relative_state(Orbit(o.state, o.et + 1, 399, 1, o.gm), o, :ric), :badarg)
        @test throws_kind(() -> relative_state(Orbit(o.state, o.et, 301, 1, o.gm), o, :ric), :badarg)
        @test throws_kind(() -> relative_state(Orbit(o.state, o.et, 399, 2, o.gm), o, :ric), :badarg)
        # ANISE fails through hvec/hmag for h = 0; no NaN leaves these calls
        for st in (SVector(7000.0, 0, 0, 7.0, 0, 0), SVector(7000.0, 0, 0, 0, 0, 0), SVector(0.0, 0, 0, 0, 7.0, 0))
            z = Orbit(st, 0.0, 399, 1, GMAT_GM)
            for k in kinds
                @test throws_kind(() -> local_frame(z, k), :degeneratecase)
                @test throws_kind(() -> local_frame_state(z, k), :degeneratecase)
                @test throws_kind(() -> relative_state(o, z, k), :degeneratecase)
            end
        end
    end
end
