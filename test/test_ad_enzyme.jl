# Group: automatic differentiation with Enzyme, forward and reverse mode. See
# `adfixtures.jl` for the behaviours and the tolerances.
#
# Two limits of Enzyme 0.13 on Julia 1.13 shape the checks. Its forward Jacobian of a
# function that allocates its result aborts the process (`x -> [sin(t) for t in x]` as
# well), so the batched forward rule is checked with a pushforward of two tangents. In
# reverse mode the cotangent of a `Dual` is not a `Dual`, so the rules refuse a `Dual`
# epoch, and ForwardDiff runs over forward mode only. Enzyme cannot prove the combined
# function of the fixture (a tuple of ten closures) read-only, so the backends declare it
# constant (`function_annotation=Enzyme.Const`); each closure on its own passes without
# it (measured), but compiling the ten one by one takes more than 4 GB. The default run
# also leaves out the three cases that Enzyme differentiates without a rule (`subpnt`,
# `unitim`, `table_position`: 9.5 min and 4.4 GB with them); each passes on its own, and
# the exhaustive run takes them. Run time of the group (default tier, one process, 2
# threads): 448 s and 3.6 GB, of which the Orbit cases about 50 s.
isdefined(@__MODULE__, :ad_check) || include("adfixtures.jl")

import Enzyme
using ADTypes: AutoEnzyme

group("ad enzyme", Kernels.FULL_SET_FILES..., Kernels.GM_TPC) do
    ctx, epochs, default, full = ad_setup()
    skip = ("subpnt", "unitim", "table_position")
    ad_check(
        "Enzyme forward", AutoEnzyme(; mode = Enzyme.Forward, function_annotation = Enzyme.Const), ctx, epochs, default, full;
        gradient = false, batched = true, nested = true, skip
    )
    rev = AutoEnzyme(; mode = Enzyme.Reverse, function_annotation = Enzyme.Const)
    ad_check(
        "Enzyme reverse", rev, ctx, epochs, default, full;
        gradient = true, batched = false, nested = false, skip
    )
    ad_orbit_check("Enzyme forward", AutoEnzyme(; mode = Enzyme.Forward, function_annotation = Enzyme.Const))
    ad_orbit_check("Enzyme reverse", rev)
    @testset "Enzyme reverse refuses a Dual epoch" begin
        f(t) = norm(spkpos(ctx, "MOON", t, "IAU_EARTH", "LT+S", "EARTH")[1])
        @test_throws ArgumentError DI.second_derivative(
            f, DI.SecondOrder(AutoForwardDiff(), rev),
            epochs[1]
        )
    end
    return nothing
end
