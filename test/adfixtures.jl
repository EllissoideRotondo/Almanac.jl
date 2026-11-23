# Fixture of the groups `ad_enzyme`, `ad_mooncake`, `ad_zygote` and `ad_tracked`:
# automatic differentiation with the packages other than ForwardDiff (`src/ad.jl` and the
# extensions for ChainRulesCore, EnzymeCore, Mooncake, ReverseDiff and Tracker), driven
# through DifferentiationInterface. Each group loads one AD package, so that a process
# compiles only one of them (Enzyme and Mooncake take about 2 GB each).
#
# The reference is ForwardDiff. Group 3 ties the ForwardDiff derivatives to CSPICE, to
# `spkder` and `pxder` and to finite differences, so these groups check that every other
# backend gives the same derivative:
#
#   behaviour                                             layer      test set
#   first derivative of every derivative point, through   invariant  "first derivatives"
#     public calls, equals the ForwardDiff derivative
#     (and through the root solvers of nearpt and prop2b)
#   the value of a call is not changed by a backend       invariant  "first derivatives"
#   batched forward rule (two tangents at once)           invariant  "batched forward"
#   gradient over many epochs (reverse mode)              invariant  "gradient"
#   second derivative with ForwardDiff over a backend     invariant  "second derivative"
#     (not Mooncake or Tracker, which take no `Dual`)
#   a refused epoch type throws (Enzyme reverse with a    error      group file
#     `Dual`; ReverseDiff and Tracker in a plain-epoch call)
#   gradient of the Orbit calls with respect to a state  invariant  "Orbit gradients"
#     or elements (elements, local frame, relative state,
#     `two_body` and `keplerian_orbit` through the Kepler
#     solver, at E = pi too; the angles at pi through
#     `_ad_pin`) equals the ForwardDiff one
#   derivative of the local frame (exhaustive run: and    invariant  "first derivatives"
#     of `:sma`) of `orbit` with respect to the epoch
#     (through `spkezr`) equals the ForwardDiff one
#
# Tolerance. A backend and ForwardDiff take the same analytical jets, so they differ only
# in the order in which the chain rule adds and multiplies the parts. `AD_BACKEND_RTOL`
# (1e-13) is a margin of about 20 over the largest measured difference (5.6e-15,
# `subpnt`); `AD_BACKEND_RTOL2` (1e-11) is that of order 2 in group 3, where the rounding
# of the second derivative is amplified (measured 4.7e-13). The names differ from
# `AD_RTOL` of group 3, because all groups share `Main` in one process. `AD_ORBIT_RTOL`
# (5e-14) is that of the gradients of the Orbit calls with respect to a state or to
# elements, a margin of about 20 over the largest measured difference (2.1e-15, `:ta` with
# Tracker); there the backends differentiate the arithmetic of the Orbit layer themselves,
# and only the Kepler solver (`_ma_to_ta`, see `src/analysis/brouwer.jl`) and the angles
# (`_ad_pin`, see `src/ad.jl`) have rules.
#
# Size. The compilation of a backend for a call takes most of the time. The default run
# differentiates one function per backend that passes through every derivative point
# (`ad_combined`); the exhaustive run (`ALMANAC_TEST_FULL=1`) adds every call of the list
# in the manual one by one, and more epochs. The test code is compiled with full
# optimization, because Enzyme differentiates the closures of these files.

import ForwardDiff
import DifferentiationInterface as DI
using ADTypes: AutoForwardDiff

const AD_BACKEND_RTOL = 1.0e-13
const AD_BACKEND_RTOL2 = 1.0e-11
const AD_ORBIT_RTOL = 5.0e-14

"Relative difference `|a - b| / |b|`; for `b == 0` the absolute difference."
ad_reldiff(a, b) = (nb = norm(b); nb == 0 ? norm(a - b) : norm(a - b) / nb)

# Fixed weights that fold a vector or a matrix result into one number, so that every
# component reaches the derivative.
const W3 = SVector(0.3, -0.5, 0.8)
const W6 = SVector(0.3, -0.5, 0.8, 1.0e3, -2.0e3, 4.0e3)
const W33 = SMatrix{3, 3}(0.1, -0.2, 0.3, 0.4, -0.5, 0.6, -0.7, 0.8, 0.9)

"""
    ad_cases(ctx, tab, rtab, dvec) -> (default, full)

Scalar functions of the epoch, by name. The default cases reach each derivative point
of `src/ad.jl` once: `spk_jet` (orders 0 and 2), `abcorr_position`, `abcorr_state`,
`abcorr_jet`, `frame_rotation` (PCK and binary PCK frames), `frame_jet` (`sxform`), the
geometry of `subpnt`, the time scales, the tables and the local frame of `orbit` (the
state of `spkezr` through the Orbit layer); and the two root solvers that run on plain
numbers, that of `nearpt` (`subpnt` with `NEAR POINT`) and that of `prop2b`, whose
derivative a `Float64` takes through `_ad_pin`. The full cases add the other calls of
the list in the manual, `subslr` with `NEAR POINT`, `conics` and an element of `orbit`.
"""
function ad_cases(ctx, tab, rtab, dvec)
    default = [
        "spkpos NONE J2000" => t -> dot(W3, spkpos(ctx, "MOON", t, "J2000", "NONE", "EARTH")[1]),
        "spkpos LT+S IAU_EARTH" => t -> dot(W3, spkpos(ctx, "MOON", t, "IAU_EARTH", "LT+S", "EARTH")[1]),
        "spkezr CN+S state, lt" => t -> (
            s = spkezr(ctx, "MARS BARYCENTER", t, "ECLIPJ2000", "CN+S", "EARTH");
            dot(W6, s[1]) + 1.0e6 * s[2]
        ),
        "spkder LT order 2" => t -> dot(W3, spkder(ctx, "MOON", t, "J2000", "LT", "EARTH", Val(2))[3]),
        "spkder NONE order 2" => t -> dot(W3, spkder(ctx, "MOON", t, "J2000", "EARTH", Val(2))[3]),
        "pxform MOON_PA" => t -> sum(W33 .* pxform(ctx, "J2000", "MOON_PA", t)),
        "sxform IAU_EARTH" => t -> sum(sxform(ctx, "J2000", "IAU_EARTH", t)[SVector(4, 5, 6), SVector(1, 2, 3)] .* W33),
        "subpnt" => t -> dot(W3, subpnt(ctx, "INTERCEPT/ELLIPSOID", "EARTH", t, "IAU_EARTH", "LT+S", "MOON")[1]),
        "unitim" => t -> unitim(ctx, t, "ET", "TDT"),
        "table_position" => t -> dot(W3, table_position(tab, t)),
        "orbit local_frame" => t -> sum(W33 .* local_frame(orbit(ctx, "MOON", t, "J2000", "NONE", "EARTH"), :ric)),
        "subpnt NEAR POINT" => t -> dot(W3, subpnt(ctx, "NEAR POINT/ELLIPSOID", "EARTH", t, "IAU_EARTH", "LT+S", "MOON")[1]),
        # dt = t / 1e7 (a few revolutions over the span of the epochs); the factor 1e5
        # gives the derivative the size of the others in the combined sum
        "prop2b" => t -> 1.0e5 * dot(W6, prop2b(1.0, SVector(1.0, 0.0, 0.0, 0.0, 1.1, 0.1), 1.0e-7 * t)),
    ]
    full = [
        "spkgeo" => t -> dot(W6, spkgeo(ctx, 301, t, "J2000", 399)[1]),
        "spkezr XCN+S" => t -> dot(W6, spkezr(ctx, "MOON", t, "IAU_EARTH", "XCN+S", "EARTH")[1]),
        "pxfrm2" => t -> sum(W33 .* pxfrm2(ctx, "IAU_EARTH", "J2000", t, t - 2.0)),
        "pxder order 2" => t -> sum(W33 .* pxder(ctx, "J2000", "IAU_EARTH", t, Val(2))[3]),
        "sincpt" => t -> dot(W3, sincpt(ctx, "ELLIPSOID", "EARTH", t, "IAU_EARTH", "CN+S", "MOON", "J2000", dvec)[1]),
        "spkcpo" => t -> dot(
            W6, spkcpo(
                ctx, "MOON", t, "J2000", "OBSERVER", "LT+S",
                SVector(1000.0, 2000.0, 3000.0), "EARTH", "IAU_EARTH"
            )[1]
        ),
        "ilumin" => t -> sum(
            ilumin(
                ctx, "ELLIPSOID", "EARTH", t, "IAU_EARTH", "LT+S", "MOON",
                SVector(6378.0, 0.0, 0.0)
            )[3:5]
        ),
        "phaseq" => t -> phaseq(ctx, t, "MOON", "SUN", "EARTH", "LT+S"),
        "deltet" => t -> deltet(ctx, t, "ET"),
        "table_rotation" => t -> sum(W33 .* table_rotation(rtab, t)),
        "orbit sma" => t -> element(orbit(ctx, "MOON", t, "J2000", "NONE", "EARTH"), :sma),
        "subslr NEAR POINT" => t -> dot(W3, subslr(ctx, "NEAR POINT/ELLIPSOID", "MOON", t, "IAU_MOON", "LT+S", "EARTH")[1]),
        "conics" => t -> 1.0e5 * dot(W6, conics(SVector(1.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.0, 1.0), 1.0e-7 * t)),
    ]
    return default, full
end

"""
    ad_combined(cases) -> function

The sum of the case functions `cases`: one function of the epoch for a backend to
compile, which passes through every derivative point that the cases reach.
"""
function ad_combined(cases)
    fs = Tuple(map(last, cases))    # a tuple, so that the sum is type-stable
    return t -> sum(map(f -> f(t), fs))
end

"Context, epochs, tables and the ray of `sincpt` of the AD groups."
function ad_setup()
    rng = MersenneTwister(SEED + 61)
    ctx = neo_context([Kernels.full_set()..., Kernels.GM_TPC])    # GM_TPC for `orbit`
    epochs = record_safe_epochs(rng, nsample(1, 4); span = 3.0e8)
    t0, t1 = extrema(epochs)
    tab = ephemeris_table(ctx, "MOON", t0 - 1.0e5, t1 + 1.0e5, "J2000", "NONE", "EARTH")
    rtab = rotation_table(ctx, "J2000", "IAU_EARTH", t0 - 1.0e5, t1 + 1.0e5)
    # A ray from the Moon to the centre of the Earth, which meets the Earth at every epoch.
    dvec = spkpos(ctx, "EARTH", epochs[1], "J2000", "NONE", "MOON")[1]
    return ctx, epochs, ad_cases(ctx, tab, rtab, dvec)...
end

"""
    ad_check(bname, backend, ctx, epochs, default, full; gradient, batched, nested,
        combine=true, skip=())

The test sets of one backend: first derivatives and values (in the default run the
combined function, or with `combine=false` each default case on its own; in the
exhaustive run every case on its own; `skip` names default cases that the default run
leaves out), and, as the flags say, the
gradient over four epochs, a pushforward of two tangents, and the second derivative with
ForwardDiff over the backend.
"""
function ad_check(
        bname, backend, ctx, epochs, default, full; gradient::Bool,
        batched::Bool, nested::Bool, combine::Bool = true, skip = ()
    )
    pos(t) = spkpos(ctx, "MOON", t, "IAU_EARTH", "LT+S", "EARTH")[1]
    range_km(t) = norm(pos(t))
    @testset "$bname" begin
        mine = filter(c -> !(first(c) in skip), default)
        cases = FULL ? vcat(default, full) : combine ? ["combined" => ad_combined(mine)] : mine
        @testset "first derivatives" begin
            worst = 0.0
            for (cname, f) in cases, t in epochs
                ref = ForwardDiff.derivative(f, t)
                y, d = DI.value_and_derivative(f, backend, t)
                err = ad_reldiff(d, ref)
                worst = max(worst, err)
                @test err <= AD_BACKEND_RTOL || (@info "derivative" bname cname t d ref; false)
                # Within 2 eps, not bit for bit: with tracked numbers the closure's own
                # arithmetic (`dot`, `sum`) runs without fused multiply-adds.
                @test isapprox(Almanac.primal(y), f(t); rtol = 2eps())
            end
            measured(
                @sprintf(
                    "  %-18s %2d function(s), max relative difference to ForwardDiff %.2e",
                    bname, length(cases), worst
                )
            )
        end
        if batched
            @testset "batched forward" begin
                d = ForwardDiff.derivative(range_km, epochs[1])
                dd = DI.pushforward(range_km, backend, epochs[1], (1.0, 2.0))
                @test ad_reldiff(collect(dd), [d, 2d]) <= AD_BACKEND_RTOL
            end
        end
        if gradient
            @testset "gradient" begin
                ts = collect(epochs[1] .+ 86400.0 .* (0:3))
                total(x) = sum(range_km, x)
                @test ad_reldiff(DI.gradient(total, backend, ts), ForwardDiff.gradient(total, ts)) <= AD_BACKEND_RTOL
            end
        end
        if nested
            @testset "second derivative" begin
                t = epochs[1]
                ref = ForwardDiff.derivative(s -> ForwardDiff.derivative(range_km, s), t)
                d2 = DI.second_derivative(range_km, DI.SecondOrder(AutoForwardDiff(), backend), t)
                @test ad_reldiff(d2, ref) <= AD_BACKEND_RTOL2
            end
        end
    end
    return nothing
end

# --- the Orbit layer: gradients with respect to a state or to elements -------------------

const AD_GM = 398600.4415
const AD_RADII = SVector(6378.1366, 6378.1366, 6356.7519)
"The elements of the default element case; `:bplane_bt` is that of the hyperbola."
const AD_ELEMENTS = (:sma, :ecc, :inc, :raan, :aop, :ta, :ma, :latitude, :brouwer_sma)

ad_orbit(s) = Orbit(s, 0.0, 399, 1, AD_GM, AD_RADII)
ad_kep(x) = keplerian_orbit(
    0.0, 399, 1; sma = x[1], ecc = x[2], inc = x[3], raan = x[4], aop = x[5],
    ma = x[6], gm = AD_GM
).state
"The weighted state of `ad_kep`: one function, so that a backend compiles it once for the three points."
ad_kep_w(x) = dot(W6, ad_kep(x))

"""
    ad_orbit_cases() -> (default, full)

Scalar functions of a vector of six numbers, by name, each with the point at which it is
differentiated: an elliptic state (periapsis 7200 km, so that the Brouwer elements are
defined), a hyperbolic state, and Keplerian elements with the mean anomaly. The default
cases are those of `keplerian_orbit`, which reach every path of the Kepler solver
`_ma_to_ta`, the one derivative point of the Orbit layer: the ellipse, the hyperbola, and
the branch at E = pi (the mean anomaly pi, where ANISE returns E for the true anomaly and
the derivative is that of the true anomaly). The three are one function (`ad_kep_w`) at
three points, so that a backend compiles it once. The fourth case is the true longitude
(the sum of :raan, :aop and :ta, without the dispatch on the name of `element`) with the
three angles at pi, where each reaches `_ad_pin` at an arc cosine of about -1. The full
cases add `two_body`, the elements, the local frames and the relative state, which the
backends differentiate as plain arithmetic. They are left out of the default run for their compilation: `two_body` takes
Mooncake 8 s for each mode, and a backend other than Zygote compiles every branch of
`element` (it dispatches on the name at run time), 4 s (Tracker) to 44 s (Enzyme
reverse). Each function is divided by its value at the point.
"""
function ad_orbit_cases()
    se = Vector(
        keplerian_orbit(
            0.0, 399, 1; sma = 8000.0, ecc = 0.1, inc = 0.5, raan = 1.0, aop = 2.0,
            ta = 0.3, gm = AD_GM
        ).state
    )
    sh = Vector(
        keplerian_orbit(
            0.0, 399, 1; sma = -20000.0, ecc = 1.5, inc = 0.7, raan = 0.4,
            aop = 1.1, ta = 0.4, gm = AD_GM
        ).state
    )
    chief = ad_orbit(se)
    sd = se .+ [1.0, -2.0, 0.5, 1.0e-3, 2.0e-3, -1.0e-3]
    xpi = [8000.0, 0.6, 0.5, 1.0, 2.0, Float64(π)]
    xe, xh = [8000.0, 0.1, 0.5, 1.0, 2.0, 2.0], [-20000.0, 1.5, 0.7, 0.4, 1.1, 0.3]
    element_case(q, s0) = "element $q" => (s -> element(ad_orbit(s), q), s0)
    # raan, aop and ta at pi, where the arc cosines of ANISE have an infinite derivative.
    # The velocity is given in m/s: its partials are zero by symmetry here, and in km/s
    # their rounding (2e-16) would count against partials of position of 1e-4.
    sa = Vector(
        keplerian_orbit(
            0.0, 399, 1; sma = 8000.0, ecc = 0.1, inc = 0.5, raan = Float64(π),
            aop = Float64(π), ta = Float64(π), gm = AD_GM
        ).state .* SVector(1.0, 1.0, 1.0, 1.0e3, 1.0e3, 1.0e3)
    )
    ad_orbit_ms(s) = ad_orbit(vcat(s[1:3], s[4:6] ./ 1.0e3))
    default = [
        "keplerian_orbit ellipse" => (ad_kep_w, xe),
        "keplerian_orbit ma = pi" => (ad_kep_w, xpi),
        "keplerian_orbit hyperbola" => (ad_kep_w, xh),
        "tlong at the apoapsis" => (s -> Almanac._el_tlong(ad_orbit_ms(s)), sa),
    ]
    full = vcat(
        [
            "two_body ellipse" => (s -> dot(W6, two_body(ad_orbit(s), 3600.0).state), se),
            "two_body hyperbola" => (s -> dot(W6, two_body(ad_orbit(s), 3600.0).state), sh),
            "local_frame ric" => (s -> sum(W33 .* local_frame(ad_orbit(s), :ric)), se),
            "relative_state ric" => (s -> dot(W6, relative_state(ad_orbit(s), chief, :ric)), sd),
            "elements" => (
                s -> sum(
                    map(q -> element(ad_orbit(s), q), AD_ELEMENTS) ./
                        map(q -> abs(element(ad_orbit(se), q)), AD_ELEMENTS)
                ), se,
            ),
            element_case(:bplane_bt, sh),
        ], [element_case(q, se) for q in AD_ELEMENTS], [
            "local_frame vnc" => (s -> sum(W33 .* local_frame(ad_orbit(s), :vnc)), se),
        ]
    )
    # Each case divided by its value at its point.
    norm1(c) = (
        f = first(last(c)); x0 = last(last(c)); s0 = abs(f(x0));
        first(c) => (x -> f(x) / s0, x0)
    )
    return map(norm1, default), map(norm1, full)
end

"""
    ad_orbit_check(bname, backend; rules=false)

The gradients of the Orbit cases with `backend` against ForwardDiff, each case on its
own: Zygote compiles the sum of the cases as one function much more slowly than the cases
one by one (measured 70 s against 50 s, before the rules for `element` and `getproperty`
of the ChainRulesCore extension). With `rules = true` the default run adds an element
and a computed property, which reach those two rules (for Zygote; the exhaustive run has
them for every backend).
"""
function ad_orbit_check(bname, backend; rules::Bool = false)
    default, full = ad_orbit_cases()
    cases = FULL ? vcat(default, full) : default
    if rules && !FULL
        se = [7000.0, 500.0, 1200.0, -1.0, 7.0, 0.5]
        cases = vcat(
            cases, [
                "element sma" => (s -> element(ad_orbit(s), :sma) / 8000.0, se),
                "property ecc" => (s -> ad_orbit(s).ecc, se),
            ]
        )
    end
    @testset "$bname: Orbit gradients" begin
        worst = 0.0
        for (cname, (f, x)) in cases
            ref = ForwardDiff.gradient(f, x)
            y, g = DI.value_and_gradient(f, backend, x)
            err = ad_reldiff(g, ref)
            worst = max(worst, err)
            @test err <= AD_ORBIT_RTOL || (@info "gradient" bname cname g ref; false)
            @test isapprox(Almanac.primal(y), f(x); rtol = 2eps())
        end
        measured(
            @sprintf(
                "  %-18s %2d Orbit function(s), max relative difference to ForwardDiff %.2e",
                bname, length(cases), worst
            )
        )
    end
    return nothing
end
