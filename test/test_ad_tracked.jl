# Group: automatic differentiation with the tracked numbers of ReverseDiff and Tracker.
# See `adfixtures.jl` for the behaviours and the tolerances. ReverseDiff runs without a
# compiled tape, because a compiled tape keeps the branches of its first epoch (segment,
# record). Run time of the group (default tier, one process, 2 threads): 40 s, of which
# the Orbit cases about 10 s.
isdefined(@__MODULE__, :ad_check) || include("adfixtures.jl")

import ReverseDiff, Tracker
using ADTypes: AutoReverseDiff, AutoTracker

group("ad tracked", Kernels.FULL_SET_FILES..., Kernels.GM_TPC) do
    ctx, epochs, default, full = ad_setup()
    # Tracker has no arithmetic between a `Dual` and a `TrackedReal`, so ForwardDiff runs
    # over ReverseDiff only.
    for (bname, backend, nested) in (("ReverseDiff", AutoReverseDiff(), true),
            ("Tracker", AutoTracker(), false))
        ad_check(bname, backend, ctx, epochs, default, full;
            gradient=true, batched=false, nested)
        ad_orbit_check(bname, backend)
        @testset "$bname: a plain-epoch call refuses a tracked epoch" begin
            lst(t) = et2lst(ctx, t, 399, 0.0, "PLANETOCENTRIC")[1]
            @test throws_kind(() -> DI.derivative(lst, backend, epochs[1]), :badarg)
        end
    end
    return nothing
end
