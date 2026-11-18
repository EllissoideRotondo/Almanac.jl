# Group: automatic differentiation with Zygote, which reads the ChainRulesCore rules.
# See `adfixtures.jl` for the behaviours and the tolerances. Run time of the group (default
# tier, one process, 2 threads): 72 s, of which the Orbit cases about 10 s.
isdefined(@__MODULE__, :ad_check) || include("adfixtures.jl")

import Zygote
using ADTypes: AutoZygote

group("ad zygote", Kernels.FULL_SET_FILES..., Kernels.GM_TPC) do
    ctx, epochs, default, full = ad_setup()
    ad_check(
        "Zygote", AutoZygote(), ctx, epochs, default, full;
        gradient = true, batched = false, nested = true
    )
    ad_orbit_check("Zygote", AutoZygote(); rules = true)
    return nothing
end
