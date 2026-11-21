# Group: ephemeris and rotation tables on a GPU array type (`table_positions!`,
# `table_states!` through the KernelAbstractions extension, and `map` over a device
# array). JLArrays is the reference GPU array type of GPUArrays: its kernels go through
# GPUCompiler to SPIR-V and run on the CPU with POCL, so a table that a GPU cannot
# compile (a value that is not `isbits`, a dynamic call, an allocation) fails here.
#
# The oracle is the CPU evaluation of the same table. The device code is compiled by
# another compiler, which may contract a product and a sum into one fused operation, so
# the limits allow rounding: 4e-16 relative for Float64, 1e-6 for Float32, and 1e-12 for
# a mixed table, whose Float32 part is 1e-12 of the position. Measured on POCL: equal to
# the last bit. Time: 12 s with the packages precompiled, most of it compiling the
# kernels (local, one process).

import Adapt
import ForwardDiff
using JLArrays: JLArray
import KernelAbstractions as KA

group("table on a device", Kernels.FULL_SET_FILES...) do
    rng = MersenneTwister(SEED + 43)
    ctx = neo_context(Kernels.full_set())
    t0, t1 = 7.0e8 + 1234.5, 7.0e8 + 60 * 86400.0
    tab = ephemeris_table(ctx, 399, t0, t1, "J2000", "NONE", 10)
    lo, hi = table_span(tab)
    ets = vcat(lo, random_epochs(rng, 500, lo, hi), hi)
    dets = JLArray(ets)
    relmax(a, b) = maximum(norm(a[j, :] - b[j, :]) / norm(b[j, 1:3]) for j in axes(b, 1))

    @testset "batch kernel" begin
        @test Base.get_extension(Almanac, :AlmanacKernelAbstractionsExt) !== nothing
        for (name, t, limit) in [
                ("Float64", tab, 4.0e-16), ("Float32", Almanac.EphemerisTable{Float32}(tab), 1.0e-6),
                ("mixed", Almanac.EphemerisTable{Float32}(tab; mixed = true), 1.0e-12),
            ]
            dt = Adapt.adapt(JLArray, t)
            @test dt.coef isa JLArray && KA.get_backend(dt.coef) == KA.get_backend(dets)
            ref = table_states!(zeros(length(ets), 6), t, ets)
            out = table_states!(JLArray(zeros(length(ets), 6)), dt, dets)
            @test out isa JLArray
            e = relmax(Array(out), ref)
            measured(@sprintf("  table on device %-7s states: %.1e relative to the CPU", name, e))
            @test e <= limit
            pos = Array(table_positions!(JLArray(zeros(Float32, length(ets), 3)), dt, dets))
            @test maximum(abs.(pos - Float32.(ref[:, 1:3])) ./ norm(ref[1, 1:3])) <= 1.0e-6
        end
    end

    @testset "map over device epochs" begin
        dt = Adapt.adapt(JLArray, tab)
        ref = [table_state(tab, et) for et in ets]
        got = Array(map(et -> @inbounds(table_state(dt, et)), dets))
        @test maximum(j -> norm(got[j] - ref[j]) / norm(ref[j][1:3]), eachindex(ets)) <= 4.0e-16
        # Dual epochs on the device: the range rate.
        speed(t, et) = ForwardDiff.derivative(s -> norm(@inbounds(table_position(t, s))), et)
        rr = Array(map(et -> speed(dt, et), dets))
        @test maximum(j -> abs(rr[j] - speed(tab, ets[j])), eachindex(ets)) <= 1.0e-12
        rot = rotation_table(ctx, "J2000", "IAU_EARTH", t0, t0 + 2 * 86400.0)
        drot = Adapt.adapt(JLArray, rot)
        rets = random_epochs(rng, 100, table_span(rot)...)
        m = Array(map(et -> @inbounds(table_rotation(drot, et)), JLArray(rets)))
        @test maximum(j -> norm(m[j] - table_rotation(rot, rets[j])), eachindex(rets)) <= 4.0e-16
    end

    @testset "errors" begin
        dt = Adapt.adapt(JLArray, tab)
        out = JLArray(zeros(3, 6))
        # The epochs are checked on the device before the kernel runs.
        @test throws_kind(() -> table_states!(out, dt, JLArray([lo, hi + 1.0, lo])), :spkinsuffdata)
        @test throws_kind(() -> table_states!(out, dt, JLArray([lo, NaN, lo])), :spkinsuffdata)
        @test Array(out) == zeros(3, 6)
        @test throws_kind(() -> table_states!(JLArray(zeros(2, 6)), dt, JLArray([lo, lo, lo])), :badarg)
        # A table or an array that is not on the device.
        @test throws_kind(() -> table_states!(out, tab, JLArray([lo, lo, lo])), :badarg)
        @test throws_kind(() -> table_states!(out, Almanac.EphemerisTable{Float32}(tab; mixed = true), JLArray([lo, lo, lo])), :badarg)
        @test throws_kind(() -> table_states!(out, dt, [lo, lo, lo]), :badarg)
        @test throws_kind(() -> table_states!(zeros(3, 6), dt, JLArray([lo, lo, lo])), :badarg)
        # An `isbits` table goes into the kernel by value.
        small = Almanac.EphemerisTable(SVector{6, Float64}(1, 2, 3, 4, 5, 6), 0.0, 10.0, 1, 2, 1, 2, 1, 0.0)
        @test Array(table_positions!(JLArray(zeros(2, 3)), small, JLArray([5.0, 7.5]))) == [1.0 3.0 5.0; 2.0 5.0 8.0]
    end
end
