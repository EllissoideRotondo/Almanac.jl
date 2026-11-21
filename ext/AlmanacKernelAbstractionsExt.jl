module AlmanacKernelAbstractionsExt

# `table_positions!` and `table_states!` on GPU arrays: one KernelAbstractions kernel,
# one work item per epoch, on the backend of the arrays (CUDA, AMDGPU, Metal, oneAPI,
# OpenCL, JLArrays). The epochs were checked by the caller, so the kernel evaluates the
# table with `@inbounds`.

using Almanac: Almanac, EphemerisTable, SpiceError
using GPUArraysCore: AbstractGPUArray
import KernelAbstractions as KA
using KernelAbstractions: @kernel, @index

# The macros unqualified: KernelAbstractions 0.9 rewrites `@index` only under that name.
@kernel function _table_kernel!(out, tab, ets, v)
    j = @index(Global, Linear)
    @inbounds Almanac._table_row!(out, tab, ets, j, v)
end

_on_device(x::AbstractGPUArray, backend) = KA.get_backend(x) == backend
_on_device(::Nothing, backend) = true
_on_device(x, backend) = isbits(x)       # an `isbits` array (an SVector) goes into the kernel by value

function Almanac._table_batch!(out::AbstractGPUArray, tab::EphemerisTable, ets::AbstractGPUArray, v::Val)
    backend = KA.get_backend(out)
    (_on_device(ets, backend) && _on_device(tab.coef, backend) && _on_device(tab.head, backend)) || throw(
        SpiceError(
            :badarg,
            "table batch: out, ets and the table must be on one device; move the table with Adapt.adapt"
        )
    )
    _table_kernel!(backend)(out, tab, ets, v; ndrange = length(ets))
    return nothing
end

@noinline _throw_mixed() = throw(
    SpiceError(:badarg, "table batch: out and ets must both be GPU arrays of one device, or both CPU arrays")
)
Almanac._table_batch!(::AbstractGPUArray, ::EphemerisTable, ::AbstractVector, ::Val) = _throw_mixed()
Almanac._table_batch!(::AbstractMatrix, ::EphemerisTable, ::AbstractGPUArray, ::Val) = _throw_mixed()

end # module
