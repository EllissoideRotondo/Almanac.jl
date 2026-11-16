# Group: automatic differentiation with Zygote, which reads the ChainRulesCore rules.
# See `adfixtures.jl` for the behaviours and the tolerances.
isdefined(@__MODULE__, :ad_check) || include("adfixtures.jl")

import Zygote
using ADTypes: AutoZygote

group("ad zygote", Kernels.FULL_SET_FILES...) do
    ctx, epochs, default, full = ad_setup()
    ad_check(
        "Zygote", AutoZygote(), ctx, epochs, default, full;
        gradient = true, batched = false, nested = true
    )
    return nothing
end
