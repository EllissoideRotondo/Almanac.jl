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
#
# Tolerance. A backend and ForwardDiff take the same analytical jets, so they differ only
# in the order in which the chain rule adds and multiplies the parts. `AD_BACKEND_RTOL`
# (1e-13) is a margin of about 20 over the largest measured difference (5.6e-15,
# `subpnt`); `AD_BACKEND_RTOL2` (1e-11) is that of order 2 in group 3, where the rounding
# of the second derivative is amplified (measured 4.7e-13). The names differ from
# `AD_RTOL` of group 3, because all groups share `Main` in one process.
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
geometry of `subpnt`, the time scales and the tables; and the two root solvers that run on plain
numbers, that of `nearpt` (`subpnt` with `NEAR POINT`) and that of `prop2b`, whose
derivative a `Float64` takes through `_ad_pin`. The full cases add the other calls of
the list in the manual, `subslr` with `NEAR POINT` and `conics`.
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
