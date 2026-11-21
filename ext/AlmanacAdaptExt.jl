module AlmanacAdaptExt

using Almanac: EphemerisTable, RotationTable
import Adapt

# Move the coefficient arrays of a table to another array type and keep the scalar
# fields. With these methods `adapt(CuArray, tab)` gives a table on the GPU, and CUDA.jl
# converts that table to an `isbits` value when it is an argument of a kernel or is
# captured by a broadcast. `cu(tab)` also converts the coefficients to `Float32`; the
# head of a mixed table must stay `Float64`, so for such a table `cu` throws and
# `adapt(CuArray, tab)` is the call to use.
Adapt.adapt_structure(to, tab::EphemerisTable) = EphemerisTable(
    Adapt.adapt(to, tab.coef),
    tab.t0, tab.dt, tab.nint, tab.ncoef, tab.target, tab.observer, tab.frame, tab.maxerr,
    Adapt.adapt(to, tab.head)
)
Adapt.adapt_structure(to, tab::RotationTable) = RotationTable(
    Adapt.adapt(to, tab.coef),
    tab.t0, tab.dt, tab.nint, tab.ncoef, tab.from, tab.to, tab.maxerr
)

end # module
