# Group: automatic differentiation with Mooncake, forward and reverse mode. See
# `adfixtures.jl` for the behaviours and the tolerances. Mooncake takes no `Dual`
# argument, so ForwardDiff does not run over it.
isdefined(@__MODULE__, :ad_check) || include("adfixtures.jl")

import Mooncake
using ADTypes: AutoMooncake, AutoMooncakeForward

group("ad mooncake", Kernels.FULL_SET_FILES...) do
    ctx, epochs, default, full = ad_setup()
    ad_check("Mooncake forward", AutoMooncakeForward(), ctx, epochs, default, full;
        gradient=false, batched=true, nested=false)
    ad_check("Mooncake reverse", AutoMooncake(), ctx, epochs, default, full;
        gradient=true, batched=false, nested=false)
    return nothing
end
