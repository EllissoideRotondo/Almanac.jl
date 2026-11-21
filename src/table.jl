# Ephemeris tables: the position of one body relative to another on a fixed interval of
# time, as one flat array of Chebyshev coefficients on intervals of equal length. A
# query needs no kernel set, no segment search and no chain of bodies, so it is suited
# to the inner loop of an optimizer and to a GPU kernel (the table is an `isbits` value
# when its array is one, see `ext/AlmanacAdaptExt.jl`).

"""
    EphemerisTable{T,A<:AbstractVector{T}}
    EphemerisTable(coef, t0, dt, nint, ncoef, target, observer, frame, maxerr)
    EphemerisTable{T}(tab; mixed=false) -> EphemerisTable{T,Vector{T}}

Position of one body relative to another on a fixed span of time, as piecewise
Chebyshev polynomials that need no kernel to evaluate. [`ephemeris_table`](@ref)
makes one; the constructor with nine arguments takes the fields of the table below.

A query does no segment search and follows no chain of bodies, so it suits the inner
loop of an optimizer. With Adapt.jl loaded, `adapt` moves the coefficients to another
array type; the table is then an `isbits` value on a GPU.

| Field              | Meaning                                                              |
|:-------------------|:---------------------------------------------------------------------|
| `target`, `observer`, `frame` | NAIF IDs of the bodies and ID of the frame                |
| `t0`, `dt`, `nint` | the table covers `t0 <= et <= t0 + nint * dt`, seconds past J2000 TDB |
| `ncoef`            | Chebyshev coefficients for each component and interval               |
| `coef`             | the coefficients, km                                                 |
| `maxerr`           | largest difference from `spkpos` that `ephemeris_table` measured, km   |
| `head`             | `nothing`, or the `Float64` part of a mixed table                    |

Interval `i` (from 0) covers `t0 + i dt` to `t0 + (i + 1) dt`. It holds `ncoef`
coefficients for each of x, y and z, in this order, from index `3 ncoef i + 1` of
`coef`: the layout of a record of SPK type 2 without its first two numbers.

`EphemerisTable{T}(tab)` returns a copy with coefficients of type `T` on the CPU. With
`Float32` the position keeps about seven digits (16 km at one astronomical unit);
`maxerr` is not measured again.

`EphemerisTable{T}(tab; mixed=true)` keeps the coefficients of degree 0 to 2 of each
series in `Float64` (the field `head`) and the others in `T`. The query returns
`Float64` values, computes the three largest terms in `Float64` and the rest in `T`.
On a GPU with slow `Float64` arithmetic, most of the work is then in `T`. Over 400 days
of `de440s.bsp` with `T = Float32`, the position of the Earth relative to the Sun stays
within 0.24 m of `Float64` (16 km for a `Float32` table), and the position of the Moon
relative to the Earth within 1.2 m. The higher coefficients are small in an inertial
frame; in a rotating frame such as `IAU_EARTH` they are not, and a mixed table gains
little (6 km for the Sun from the Earth). Move a mixed table to a GPU with
`adapt(CuArray, tab)`. `cu(tab)` converts the head to `Float32` and throws.

Throws `SpiceError(:badarg)` if the constructor gets `dt <= 0`, `nint < 1`,
`ncoef < 1` or `length(coef) != 3 * ncoef * nint`, or a head that is not `Float64`.

# Examples
```jldoctest
julia> coef = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];    # nint = 1 interval, ncoef = 2

julia> tab = EphemerisTable(coef, 0.0, 8.0, 1, 2, 301, 399, 1, 0.0);

julia> table_span(tab)
(0.0, 8.0)

julia> table_state(tab, 6.0)    # x = 1 + 2s, y = 3 + 4s, z = 5 + 6s at s = 0.5
6-element StaticArraysCore.SVector{6, Float64} with indices SOneTo(6):
 2.0
 5.0
 8.0
 0.5
 1.0
 1.5

julia> table_position(EphemerisTable{Float32}(tab; mixed = true), 6.0) == [2, 5, 8]
true
```

See also [`table_position`](@ref), [`table_state`](@ref), [`table_jet`](@ref),
[`table_positions!`](@ref), [`table_to_spk`](@ref), [`RotationTable`](@ref).
"""
struct EphemerisTable{T <: AbstractFloat, A <: AbstractVector{T}, H <: Union{Nothing, AbstractVector{Float64}}}
    coef::A
    t0::Float64
    dt::Float64
    nint::Int
    ncoef::Int
    target::Int32
    observer::Int32
    frame::Int32
    maxerr::Float64
    head::H
    rdt::Float64        # 1 / dt, so that a query multiplies where it would divide
    function EphemerisTable(
            coef::A, t0::Real, dt::Real, nint::Integer, ncoef::Integer,
            target::Integer, observer::Integer, frame::Integer, maxerr::Real,
            head::H = nothing
        ) where {T <: AbstractFloat, A <: AbstractVector{T}, H <: Union{Nothing, AbstractVector{Float64}}}
        (dt > 0 && nint >= 1 && ncoef >= 1 && length(coef) == 3 * ncoef * nint) ||
            throw(SpiceError(:badarg, "EphemerisTable: the sizes do not agree"))
        (head === nothing || length(head) == 3 * _TAB_HEAD * nint) ||
            throw(SpiceError(:badarg, "EphemerisTable: the sizes of the head do not agree"))
        return new{T, A, H}(coef, t0, dt, nint, ncoef, target, observer, frame, maxerr, head, 1 / dt)
    end
end

# A head of another element type: `cu(tab)` converts every Float64 array to Float32.
EphemerisTable(
    coef::AbstractVector{<:AbstractFloat}, t0::Real, dt::Real, nint::Integer, ncoef::Integer,
    target::Integer, observer::Integer, frame::Integer, maxerr::Real, head::AbstractVector
) = throw(
    SpiceError(
        :badarg,
        "EphemerisTable: the head of a mixed table must be Float64; move the table with Adapt.adapt, not with a conversion to Float32"
    )
)

const _TAB_HEAD = 3      # coefficients of each series in the head of a mixed table

"""
    _full_coef(tab) -> Vector{Float64}

The coefficients of `tab` in the layout of `coef`, on the CPU, as `Float64`. A mixed
table holds the coefficients of degree 0 to `_TAB_HEAD - 1` of each series in `head`
(in the layout of `coef` with `_TAB_HEAD` coefficients for each series) and zeros in
their places in `coef`.
"""
function _full_coef(tab::EphemerisTable)
    coef = Vector{Float64}(Array(tab.coef))
    tab.head === nothing && return coef
    head, n, h = Array(tab.head), tab.ncoef, _TAB_HEAD
    for i in 0:(tab.nint - 1), c in 0:2, k in 0:(min(n, h) - 1)   # invariant: the terms before (i, c, k) are complete
        coef[3n * i + n * c + k + 1] = head[3h * i + h * c + k + 1]
    end
    return coef
end

function EphemerisTable{T}(tab::EphemerisTable; mixed::Bool = false) where {T <: AbstractFloat}
    full = _full_coef(tab)
    mixed || return EphemerisTable(
        Vector{T}(full), tab.t0, tab.dt, tab.nint, tab.ncoef, tab.target,
        tab.observer, tab.frame, tab.maxerr
    )
    n, h = tab.ncoef, _TAB_HEAD
    head = zeros(3 * h * tab.nint)
    for i in 0:(tab.nint - 1), c in 0:2, k in 0:(min(n, h) - 1)   # invariant: the terms before (i, c, k) are split
        a = 3n * i + n * c + k + 1
        head[3h * i + h * c + k + 1] = full[a]
        full[a] = 0.0
    end
    return EphemerisTable(
        Vector{T}(full), tab.t0, tab.dt, tab.nint, tab.ncoef, tab.target,
        tab.observer, tab.frame, tab.maxerr, head
    )
end

"""
    table_span(tab) -> (first, last)

Return the first and the last epoch of the table `tab`, in seconds past J2000 TDB.

`tab` is an [`EphemerisTable`](@ref) or a [`RotationTable`](@ref). A query is valid
for `first <= et <= last`. A table that [`ephemeris_table`](@ref) made from the records
of the kernels covers whole records, so its span can be wider than the requested one.
"""
@inline table_span(tab::EphemerisTable) = (tab.t0, tab.t0 + tab.nint * tab.dt)

Base.show(io::IO, tab::EphemerisTable{T}) where {T} = print(
    io, "EphemerisTable{", T, "}(", tab.head === nothing ? "" : "Float64 head, ",
    tab.target, " relative to ", tab.observer, ", frame ", tab.frame, ", ET ", tab.t0,
    " to ", table_span(tab)[2], ", ", tab.nint, " intervals of ", tab.dt, " s, ",
    tab.ncoef, " coefficients, largest error ", tab.maxerr, " km)"
)

@noinline _throw_table_range() =
    throw(SpiceError(:spkinsuffdata, "the epoch is outside the ephemeris table"))

@inline _tab_arg(::Type{T}, s::AbstractFloat) where {T} = T(s)
@inline _tab_arg(::Type, s::Real) = s

"""
    table_jet(tab, et, Val(N)) -> NTuple{N+1,SVector{3}}

Return the position from the ephemeris table `tab` at `et`, in km, and its time
derivatives to order `N`.

Element `k + 1` is the derivative of order `k`, in km/s^k; every order is a
derivative of the position polynomial. `et` is in seconds past J2000 TDB and is a
plain number or a dual number. The interval and the normalized time are computed in
`Float64` also when the coefficients are `Float32`; the polynomial is then summed in
`Float32`, except the head of a mixed table (see [`EphemerisTable`](@ref)), and
a mixed table returns `Float64`. The call reads no kernel and does not allocate.

Throws `SpiceError(:spkinsuffdata)` if `et` is outside [`table_span`](@ref)`(tab)`;
this is the kind of [`spkezr`](@ref) for an epoch outside the loaded data. The check is
a bounds check. Under `@inbounds`, as in a GPU kernel, it is removed, and an
epoch outside the table then gives the extrapolation of the nearest interval.

See also [`table_position`](@ref), [`table_state`](@ref).
"""
Base.@propagate_inbounds function table_jet(tab::EphemerisTable, et::Real, v::Val)
    i, s, rr = _tab_locate(tab, et)
    return _tab_eval(tab.head, tab, i, s, rr, v)
end

@inline _tab_eval(::Nothing, tab::EphemerisTable{T}, i, s, rr, v::Val) where {T} =
    _tab_scale(_cheb3(tab.coef, 3 * tab.ncoef * i + 1, tab.ncoef, _tab_arg(T, s), one(T), v), _tab_arg(T, rr))

# A mixed table: the head (degrees 0 to 2) is summed in Float64 at `s`, the rest in `T`
# at `s` rounded to `T`. Measured on de440s over 400 days with a Float32 tail: 0.24 m
# for the Earth relative to the Sun, against 16 km for a Float32 table and 14 m with a
# head of degrees 0 and 1 only, where the Float32 rounding of the Clenshaw sum of the
# degree-2 term (4.5e4 km) dominated. A first-order correction for the rounding of `s`
# gained 30 % for twice the Float32 work, and was left out.
@inline function _tab_eval(head::AbstractVector{Float64}, tab::EphemerisTable{T}, i, s, rr, v::Val{N}) where {T, N}
    tail = _tab_scale(_cheb3(tab.coef, 3 * tab.ncoef * i + 1, tab.ncoef, _tab_arg(T, s), one(T), v), _tab_arg(T, rr))
    top = _tab_scale(_cheb3(head, 3 * _TAB_HEAD * i + 1, _TAB_HEAD, s, 1.0, v), rr)
    R = promote_type(eltype(top[1]), eltype(tail[1]))
    return ntuple(m -> Vec3{R}(top[m] + tail[m]), Val(N + 1))
end

"""
    _tab_locate(tab, et) -> (i, s, 1 / radius)

The interval `i` (from 0) of `et`, the normalized time `s` in it, and the reciprocal
of the half length of an interval. The interval is chosen on `primal(et)`; the bounds
check throws `:spkinsuffdata` outside the table. The epoch arithmetic is `Float64`
(or the type of `et`) and multiplies by `tab.rdt` in place of a division, which on
a GPU costs a sequence of `Float64` operations.
"""
Base.@propagate_inbounds function _tab_locate(tab, et::Real)
    t = primal(et)
    @boundscheck (tab.t0 <= t <= tab.t0 + tab.nint * tab.dt) || _throw_table_range()
    i = clamp(unsafe_trunc(Int, (t - tab.t0) * tab.rdt), 0, tab.nint - 1)
    rr = 2 * tab.rdt
    s = (et - (tab.t0 + (2i + 1) * (tab.dt / 2))) * rr
    return i, s, rr
end

"Order `m` of the jet (derivatives with respect to `s`) times `rr^m`: derivatives with respect to time."
@inline function _tab_scale(jet::NTuple{M}, rr) where {M}
    return ntuple(Val(M)) do m
        v = jet[m]
        for _ in 2:m      # invariant: `v` has been multiplied by `rr` once per pass
            v = v * rr
        end
        v
    end
end

"""
    table_position(tab, et) -> SVector{3}

Return the position from the ephemeris table `tab` at `et`, in km.

`et` is in seconds past J2000 TDB. The bodies and the frame are `tab.target`,
`tab.observer` and `tab.frame`. The rules for the epoch, for `@inbounds` and for the
number types are those of [`table_jet`](@ref).

Throws `SpiceError(:spkinsuffdata)` if `et` is outside the table.
"""
Base.@propagate_inbounds table_position(tab::EphemerisTable, et::Real) = table_jet(tab, et, Val(0))[1]

"""
    table_state(tab, et) -> SVector{6}

Return the position (km) and the velocity (km/s) from the ephemeris table `tab` at
`et`.

The velocity is the derivative of the position polynomial. For SPK data whose
velocity is the derivative of the position, such as type 2, it is therefore order 1 of
[`spkder`](@ref) to the accuracy of the table; with an `LT` correction, that differs
from the velocity of [`spkezr`](@ref). The rules for the epoch, for `@inbounds` and for
the number types are those of [`table_jet`](@ref).

Throws `SpiceError(:spkinsuffdata)` if `et` is outside the table.
"""
Base.@propagate_inbounds function table_state(tab::EphemerisTable, et::Real)
    r, v = table_jet(tab, et, Val(1))
    return vcat(r, v)
end

"""
    table_positions!(out, tab, ets) -> out

Write the position (km) from the ephemeris table `tab` at each epoch of `ets` into
`out`, an `n` by 3 matrix for `n = length(ets)`. Row `j` holds the position at
`ets[j]`, one component per column. This call has no CSPICE counterpart.

The values are those of [`table_position`](@ref), converted to the element type of
`out`. If `out`, `ets` and the coefficients of `tab` are GPU
arrays of one device (`tab` moved with `Adapt.adapt`), the call runs one
KernelAbstractions kernel on that device and returns without waiting for it. It
needs `KernelAbstractions` and `GPUArraysCore` in the environment, which the GPU
packages load. On the CPU it is a loop that does not allocate. A column per
component lets neighbouring threads of a GPU write neighbouring addresses.

Throws:
- `SpiceError(:badarg)` if the size of `out` does not agree with `ets`, or if a GPU
  array is mixed with an array or a table that is not on the same device;
- `SpiceError(:spkinsuffdata)` if an epoch is outside [`table_span`](@ref)`(tab)` or
  is `NaN`. The epochs are checked before anything is written, also on a GPU;
- `ArgumentError` if the indices of `out` or `ets` do not start at 1.

# Examples
```jldoctest
julia> tab = EphemerisTable([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 0.0, 10.0, 1, 2, 301, 399, 1, 0.0);

julia> table_positions!(zeros(2, 3), tab, [5.0, 7.5])
2×3 Matrix{Float64}:
 1.0  3.0  5.0
 2.0  5.0  8.0
```

See also [`table_states!`](@ref).
"""
function table_positions!(out::AbstractMatrix, tab::EphemerisTable, ets::AbstractVector)
    _check_batch(out, tab, ets, 3)
    isempty(ets) || _table_batch!(out, tab, ets, Val(0))
    return out
end

"""
    table_states!(out, tab, ets) -> out

Write the position and the velocity (km, km/s) from the ephemeris table `tab` at each
epoch of `ets` into `out`, an `n` by 6 matrix for `n = length(ets)`: the values of
[`table_state`](@ref), one row per epoch. Devices, errors and the order of the checks
are those of [`table_positions!`](@ref).
"""
function table_states!(out::AbstractMatrix, tab::EphemerisTable, ets::AbstractVector)
    _check_batch(out, tab, ets, 6)
    isempty(ets) || _table_batch!(out, tab, ets, Val(1))
    return out
end

function _check_batch(out, tab, ets, width)
    Base.require_one_based_indexing(out, ets)
    size(out) == (length(ets), width) || throw(
        SpiceError(:badarg, "table batch: out must have length(ets) rows and $width columns, not size $(size(out))")
    )
    isempty(ets) && return nothing
    lo, hi = extrema(primal, ets)
    (tab.t0 <= lo && hi <= tab.t0 + tab.nint * tab.dt) || _throw_table_range()
    return nothing
end

# Row `j` of `out` from the jet of order N (position, or position and velocity).
Base.@propagate_inbounds function _table_row!(out, tab, ets, j, v::Val{N}) where {N}
    jet = table_jet(tab, ets[j], v)
    for o in 0:N, c in 1:3      # invariant: the columns before (o, c) of row j are set
        out[j, 3o + c] = jet[o + 1][c]
    end
    return nothing
end

# CPU loop; the KernelAbstractions extension adds the methods for GPU arrays.
function _table_batch!(out::AbstractMatrix, tab, ets::AbstractVector, v::Val)
    @inbounds for j in eachindex(ets)   # invariant: the rows before j are set
        _table_row!(out, tab, ets, j, v)
    end
    return nothing
end

# ---------------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------------

"Segments from `body` to the root of its chain at `t`, or `nothing` if the chain is too long."
function _tab_up(ks::KernelSet, body::Int32, t::Float64)
    list = Segment[]
    while true     # invariant: `list` is the chain from the first body to `body`
        segs = get(ks.spk, body, nothing)
        segs === nothing && return list
        i = find_segment(segs, t)
        i == 0 && return list
        length(list) == MAXCHAIN && return nothing
        push!(list, segs[i])
        body = segs[i].center
    end
    return
end

"""
    _tab_chain(ks, targ, obs, t) -> Union{Nothing,Vector{Segment}}

The segments that `spk_jet` adds and subtracts for `targ` relative to `obs` at `t`,
or `nothing` if the two bodies are not connected at `t`.
"""
function _tab_chain(ks::KernelSet, targ::Int32, obs::Int32, t::Float64)
    a, b = _tab_up(ks, targ, t), _tab_up(ks, obs, t)
    (a === nothing || b === nothing) && return nothing
    root(list, body) = isempty(list) ? body : last(list).center
    root(a, targ) == root(b, obs) || return nothing
    while !isempty(a) && !isempty(b) && last(a) == last(b)   # invariant: the tails removed are equal
        pop!(a)
        pop!(b)
    end
    return append!(a, b)
end

"""
    _tab_native_grid(ks, targ, ref, obs, t0, t1) -> Union{Nothing,Tuple{Float64,Float64,Int,Int}}

`(start, dt, nint, ncoef)` of the grid on which the position of `targ` relative to
`obs` in frame `ref` is one polynomial of degree `ncoef - 1` for each interval, or
`nothing` if there is no such grid. A grid exists when `ref` is inertial, every
segment of the chain is SPK type 2 in an inertial frame, the record boundaries of
all segments fall on the boundaries of the shortest record, and the same segments
serve the whole of `start` to `start + nint * dt`, which contains `t0` to `t1`.
"""
function _tab_native_grid(ks::KernelSet, targ::Int32, ref::Int32, obs::Int32, t0::Float64, t1::Float64)
    _is_inertial(ks, ref) || return nothing
    chain = _tab_chain(ks, targ, obs, (t0 + t1) / 2)
    (chain === nothing || isempty(chain)) && return nothing
    all(s -> s.dtype == Int32(2) && _is_inertial(ks, s.frame), chain) || return nothing
    short = chain[argmin(map(s -> s.intlen, chain))]
    dt, origin = short.intlen, short.init
    whole(x) = abs(x - round(x)) <= 1.0e-9
    all(s -> whole((s.init - origin) / dt) && whole(s.intlen / dt), chain) || return nothing
    i0 = floor(Int, (t0 - origin) / dt)
    i1 = max(i0 + 1, ceil(Int, (t1 - origin) / dt))
    start, nint = origin + i0 * dt, i1 - i0
    slack = 1.0e-9 * dt
    all(s -> s.t0 <= start + slack && s.t1 >= start + nint * dt - slack, chain) || return nothing
    for i in 0:(nint - 1)      # invariant: the chain of every interval before `i` is `chain`
        _tab_chain(ks, targ, obs, start + (i + 0.5) * dt) == chain || return nothing
    end
    ncoef = maximum(s -> (Int(s.rsize) - 2) ÷ 3, chain)
    return start, dt, nint, ncoef
end

"""
    _tab_fit(f, start, dt, nint, ncoef) -> (coef::Vector{Float64}, maxerr, scale, velerr)

Chebyshev coefficients of `f(et)::SVector{3}` on `nint` intervals of length `dt` from
`start`: interpolation at the `ncoef` zeros of the Chebyshev polynomial of degree
`ncoef`. An epoch near 1e9 s has a resolution of 1e-7 s, so a zero is not an exact
epoch; the function solves the interpolation conditions at the normalized times of
the epochs that `f` received, which keeps a polynomial `f` exact to rounding.
`maxerr` is the largest norm of the difference between `f` and the interpolant at the
`ncoef - 1` points between the zeros (the points where the error of the interpolation
is largest); `scale` is the largest norm of `f` at the zeros.
Requires `dt > 0`, `nint >= 1`, `ncoef >= 2`.
"""
_tab_fit(f::F, start::Float64, dt::Float64, nint::Int, ncoef::Int) where {F} =
    _tab_fit(f, start, dt, nint, ncoef, Val(1))

"The value at `s` of the `C` vector series of one record of `ncoef` coefficients from index `a`, as one vector."
@inline function _tab_value(coef, a::Int, n::Int, s, radius, ::Val{C}) where {C}
    cols = ntuple(c -> _cheb3(coef, a + 3 * n * (c - 1), n, s, radius, Val(0))[1], Val(C))
    return vcat(cols...)
end

"""
    _tab_fit(f, start, dt, nint, ncoef, Val(C); split=false) -> (coef, maxerr, scale, velerr)

The fit above for a function `f(et)` of `3C` components (for `C = 3`, the columns of a
matrix one after the other): a record holds `ncoef` coefficients for each component,
in the order of the components, and the two norms are norms of all `3C` components.
With `split` (`C = 2`, a position and a velocity), `maxerr` and `scale` are norms of the
first three components and `velerr` is the error of the last three; else `velerr` is
`NaN`.
"""
function _tab_fit(
        f::F, start::Float64, dt::Float64, nint::Int, ncoef::Int, vc::Val{C};
        split::Bool = false
    ) where {F, C}
    @assert !split || C == 2
    @assert dt > 0 && nint >= 1 && ncoef >= 2
    n = ncoef
    x = [cospi((j - 0.5) / n) for j in 1:n]
    xc = [cospi(j / n) for j in 1:(n - 1)]
    w = 3 * C * n                       # numbers in one record
    coef = zeros(w * nint)
    vals = Matrix{Float64}(undef, n, 3 * C)
    basis = Matrix{Float64}(undef, n, n)
    radius = dt / 2
    scale = 0.0
    for i in 0:(nint - 1)      # invariant: the coefficients of the intervals before `i` are set
        mid = start + (2i + 1) * radius
        for j in 1:n
            et = mid + radius * x[j]
            s = (et - mid) / radius
            r = f(et)
            vals[j, :] .= r
            scale = max(scale, split ? norm(view(r, 1:3)) : norm(r))
            basis[j, 1] = 1.0
            basis[j, 2] = s
            for k in 3:n       # invariant: basis[j, m] is the polynomial of degree m - 1 at s
                basis[j, k] = 2 * s * basis[j, k - 1] - basis[j, k - 2]
            end
        end
        coef[(w * i + 1):(w * (i + 1))] .= vec(basis \ vals)
    end
    maxerr = 0.0
    velerr = split ? 0.0 : NaN
    for i in 0:(nint - 1), sc in xc
        mid = start + (2i + 1) * radius
        et = mid + radius * sc
        d = _tab_value(coef, w * i + 1, n, (et - mid) / radius, radius, vc) - f(et)
        if split
            maxerr = max(maxerr, norm(view(d, 1:3)))
            velerr = max(velerr, norm(view(d, 4:6)))
        else
            maxerr = max(maxerr, norm(d))
        end
    end
    @assert length(coef) == w * nint
    return coef, maxerr, scale, velerr
end

"""
    _tab_refine(f, a, b, tol, ncoef, maxintervals, Val(C), who, unit; split=false)
        -> (coef, dt, nint, maxerr, velerr)

Fit `f` on `a` to `b` with `_tab_fit` and double the number of intervals until the
measured error is at most `tol`. Without `tol` (`nothing`) the goal is `1e-12` times
the largest norm of `f`, and the loop stops also when the error is small and one more
doubling did not divide it by four: the source then limits the accuracy. `split` is that
of `_tab_fit`: the goal is then on the position, and `velerr` is the velocity error of
the accepted fit. Throws `SpiceError(:noconvergence)` if `maxintervals` intervals do not
reach the goal.
"""
function _tab_refine(
        f::F, a::Float64, b::Float64, tol::Union{Nothing, Real}, ncoef::Int,
        maxintervals::Integer, vc::Val{C}, who::AbstractString, unit::AbstractString;
        split::Bool = false
    ) where {F, C}
    nint, previous = 1, Inf
    while true     # invariant: no number of intervals below `nint` was accepted
        dt = (b - a) / nint
        coef, maxerr, scale, velerr = _tab_fit(f, a, dt, nint, ncoef, vc; split = split)
        goal = tol === nothing ? 1.0e-12 * scale : Float64(tol)
        floor_reached = tol === nothing && maxerr <= 1.0e-6 * scale && maxerr > previous / 4
        (maxerr <= goal || floor_reached) && return coef, dt, nint, maxerr, velerr
        2 * nint > maxintervals && throw(
            SpiceError(
                :noconvergence,
                "$who: $nint intervals give an error of $maxerr$unit, above the tolerance of $goal$unit"
            )
        )
        nint, previous = 2 * nint, maxerr
    end
    return
end

"""
    ephemeris_table([ctx,] targ, t0, t1, ref, abcorr, obs;
                    tol=nothing, ncoef=16, maxintervals=65536) -> EphemerisTable

Return a table of the position of `targ` relative to `obs` in frame `ref` for
`t0 <= et <= t1`, which [`table_position`](@ref), [`table_state`](@ref) and
[`table_jet`](@ref) evaluate without the kernels. This call has no CSPICE counterpart.

`t0` and `t1` are in seconds past J2000 TDB. `targ`, `ref`, `abcorr` and `obs` are as
for [`spkpos`](@ref). The table is a snapshot: a later change of the kernels does not
change it. The result is an `EphemerisTable{Float64,Vector{Float64}}`, and its field
`maxerr` is the measured error in km.

Two methods give the coefficients:

- **Records of the kernels.** If `abcorr` is `"NONE"`, `ref` is inertial, and every
  segment between the two bodies is SPK type 2 in an inertial frame with record
  boundaries on one common grid, the table takes the record length and the degree of
  the kernels. It then agrees with `spkpos` to rounding, and `tol`, `ncoef` and
  `maxintervals` are not used. The table covers whole records of the shortest record
  length, at least `t0` to `t1`: it can start up to one record before `t0` and end up
  to one record after `t1`, and a query in that margin does not throw.
- **Interpolation.** In all other cases the call interpolates on `t0` to `t1` with
  `ncoef` coefficients for each interval, and it doubles the number of intervals
  until the measured error is at most `tol` (km). With `tol`, the result has
  `maxerr <= tol` or the call throws. Without `tol` the goal is `1e-12` times the
  largest distance in the table, and the call also stops when the error is below
  `1e-6` times that distance and one more doubling does not divide it by four.
  `spkpos` then limits the accuracy: rounding of a rotation angle of many turns, or the
  small steps of an aberration-corrected position at the record boundaries.

Throws `SpiceError(:badarg)` if the epochs are not finite with `t0 < t1`, if
`ncoef < 2`, `tol <= 0` or `maxintervals < 1`. Throws `SpiceError(:noconvergence)` if
`maxintervals` intervals do not reach the goal, and the errors of [`spkpos`](@ref) if
the kernels do not cover `t0` to `t1`.

# Examples
```julia
furnsh("naif0012.tls"); furnsh("de440s.bsp")
t0, t1 = str2et("2025-01-01"), str2et("2026-01-01")
tab = ephemeris_table("MARS BARYCENTER", t0, t1, "ECLIPJ2000", "NONE", "SUN")
tab.maxerr                              # km
r = table_position(tab, t0 + 86400.0)   # km, no kernel access
x = table_state(tab, t0 + 86400.0)      # km, km/s
tab32 = EphemerisTable{Float32}(tab)    # smaller, about seven digits
```

See also [`EphemerisTable`](@ref), [`rotation_table`](@ref), [`table_to_spk`](@ref).
"""
function ephemeris_table(
        ks::KernelSet, targ::NameOrID, t0::Real, t1::Real, ref::NameOrID,
        abcorr::Union{AbstractString, Aberration}, obs::NameOrID;
        tol::Union{Nothing, Real} = nothing, ncoef::Integer = 16, maxintervals::Integer = 65536
    )
    a, b = Float64(t0), Float64(t1)
    (isfinite(a) && isfinite(b) && a < b) ||
        throw(SpiceError(:badarg, "ephemeris_table: the epochs must be finite, with t0 < t1"))
    ncoef >= 2 || throw(SpiceError(:badarg, "ephemeris_table: ncoef must be at least 2"))
    (tol === nothing || tol > 0) || throw(SpiceError(:badarg, "ephemeris_table: tol must be positive"))
    maxintervals >= 1 || throw(SpiceError(:badarg, "ephemeris_table: maxintervals must be at least 1"))
    t, f, o = body_code(ks, targ), frame_id(ks, ref), body_code(ks, obs)
    corr = Aberration(abcorr)
    source = corr.geometric ? (et -> spk_jet(ks, t, et, f, o, Val(0))[1]) :
        (et -> abcorr_position(ks, t, et, f, corr, o)[1])
    source(a), source(b)       # data at both ends, or the error of the SPK layer
    grid = corr.geometric ? _tab_native_grid(ks, t, f, o, a, b) : nothing
    if grid !== nothing
        start, dt, nint, n = grid
        coef, maxerr, _ = _tab_fit(source, start, dt, nint, n)
        tab = EphemerisTable(coef, start, dt, nint, n, t, o, f, maxerr)
    else
        coef, dt, nint, maxerr = _tab_refine(
            source, a, b, tol, Int(ncoef), maxintervals, Val(1),
            "ephemeris_table", " km"
        )
        tab = EphemerisTable(coef, a, dt, nint, ncoef, t, o, f, maxerr)
    end
    lo, hi = table_span(tab)
    @assert lo <= a && b <= hi + 1.0e-9 * tab.dt "the table does not cover the request"
    return tab
end
ephemeris_table(
    ctx::SpiceContext, targ::NameOrID, t0::Real, t1::Real, ref::NameOrID,
    abcorr::Union{AbstractString, Aberration}, obs::NameOrID; kw...
) =
    ephemeris_table(kernelset(ctx), targ, t0, t1, ref, abcorr, obs; kw...)
ephemeris_table(
    targ::NameOrID, t0::Real, t1::Real, ref::NameOrID,
    abcorr::Union{AbstractString, Aberration}, obs::NameOrID; kw...
) =
    ephemeris_table(default_context(), targ, t0, t1, ref, abcorr, obs; kw...)

# ---------------------------------------------------------------------------------
# Export of a table as an SPK segment
# ---------------------------------------------------------------------------------

"""
    table_to_spk([ctx,] handle, tab, segid; first=table_span(tab)[1],
                 last=table_span(tab)[2]) -> Nothing

Add the ephemeris table `tab` to the SPK file `handle` as one segment of type 2.

The records of the table are the records of the segment. [`spkpos`](@ref) on the new
file therefore gives the position of the table to rounding, and the velocity is its
derivative. `handle` comes from [`spkopn`](@ref), [`spkopa`](@ref) or [`spkops`](@ref); `segid` is the
segment identifier; `first` and `last` are the bounds of the segment in seconds past
J2000 TDB. The segment is for `tab.target` relative to `tab.observer` in the frame
`tab.frame`, which must be a frame known in `ctx` (default: the default context; see
[`frinfo`](@ref)). The coefficients are written as `Float64`; a table on a GPU is copied
to the CPU.

A segment has no aberration correction. Write only a table that was made with
`abcorr = "NONE"`; the call cannot check this.

Throws `SpiceError(:invalidrefframe)` if `tab.frame` is not a frame of `ctx`, as
[`spkw02_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw02_c.html) for an unknown frame, and the other errors of [`spkw02`](@ref). Among its
conditions: `tab.ncoef <= 28`, and `first` and `last` inside the table.

# Examples
```julia
tab = ephemeris_table("MOON", t0, t1, "J2000", "NONE", "EARTH")
h = spkopn("moon.bsp", "moon.bsp", 0)
table_to_spk(h, tab, "Moon from table")
spkcls(h)
```
"""
function table_to_spk(
        ctx::SegContext, h::DAFHandle, tab::EphemerisTable,
        segid::AbstractString; first::Real = table_span(tab)[1], last::Real = table_span(tab)[2]
    )
    frinfo(_seg_ks(ctx), tab.frame) === nothing && throw(
        SpiceError(
            :invalidrefframe,
            "table_to_spk: the reference frame $(tab.frame) of the table is not known"
        )
    )
    cdata = reshape(_full_coef(tab), 3 * tab.ncoef, tab.nint)
    spkw02(
        ctx, h, tab.target, tab.observer, Int(tab.frame), first, last, segid, tab.dt,
        cdata, tab.t0
    )
    return nothing
end
table_to_spk(
    h::DAFHandle, tab::EphemerisTable, segid::AbstractString;
    first::Real = table_span(tab)[1], last::Real = table_span(tab)[2]
) =
    table_to_spk(default_context(), h, tab, segid; first = first, last = last)

# ---------------------------------------------------------------------------------
# Rotation tables
# ---------------------------------------------------------------------------------

"""
    RotationTable{T,A<:AbstractVector{T}}
    RotationTable{T}(tab) -> RotationTable{T,Vector{T}}

Rotation matrix from one frame to another on a fixed span of time, as piecewise
Chebyshev polynomials that need no kernel to evaluate: the counterpart of
[`EphemerisTable`](@ref) for [`pxform`](@ref). [`rotation_table`](@ref) makes one.

| Field              | Meaning                                                              |
|:-------------------|:---------------------------------------------------------------------|
| `from`, `to`       | frame IDs; the matrix rotates a vector from `from` to `to`           |
| `t0`, `dt`, `nint` | the table covers `t0 <= et <= t0 + nint * dt`, seconds past J2000 TDB |
| `ncoef`            | Chebyshev coefficients for each matrix element and interval          |
| `coef`             | the coefficients                                                     |
| `maxerr`           | largest Frobenius norm of the difference from `pxform` that `rotation_table` measured   |

Interval `i` (from 0) holds `ncoef` coefficients for each of the nine elements of the
matrix, column after column, from index `9 ncoef i + 1` of `coef`. A matrix of the
table is orthogonal to the accuracy `maxerr`; it is not made orthogonal again. A query
outside the table throws `SpiceError(:spkinsuffdata)`, as for an `EphemerisTable`.

`RotationTable{T}(tab)` returns a copy with coefficients of type `T` on the CPU;
`maxerr` is not measured again.

Throws `SpiceError(:badarg)` if the constructor gets `dt <= 0`, `nint < 1`,
`ncoef < 1` or `length(coef) != 9 * ncoef * nint`.

See also [`table_rotation`](@ref), [`table_rotation_jet`](@ref), [`table_span`](@ref).
"""
struct RotationTable{T <: AbstractFloat, A <: AbstractVector{T}}
    coef::A
    t0::Float64
    dt::Float64
    nint::Int
    ncoef::Int
    from::Int32
    to::Int32
    maxerr::Float64
    rdt::Float64        # 1 / dt
    function RotationTable(
            coef::A, t0::Real, dt::Real, nint::Integer, ncoef::Integer,
            from::Integer, to::Integer, maxerr::Real
        ) where {T <: AbstractFloat, A <: AbstractVector{T}}
        (dt > 0 && nint >= 1 && ncoef >= 1 && length(coef) == 9 * ncoef * nint) ||
            throw(SpiceError(:badarg, "RotationTable: the sizes do not agree"))
        return new{T, A}(coef, t0, dt, nint, ncoef, from, to, maxerr, 1 / dt)
    end
end

RotationTable{T}(tab::RotationTable) where {T <: AbstractFloat} =
    RotationTable(
    Vector{T}(tab.coef), tab.t0, tab.dt, tab.nint, tab.ncoef, tab.from, tab.to,
    tab.maxerr
)

@inline table_span(tab::RotationTable) = (tab.t0, tab.t0 + tab.nint * tab.dt)

Base.show(io::IO, tab::RotationTable{T}) where {T} = print(
    io, "RotationTable{", T, "}(frame ",
    tab.from, " to frame ", tab.to, ", ET ", tab.t0, " to ", table_span(tab)[2], ", ", tab.nint,
    " intervals of ", tab.dt, " s, ", tab.ncoef, " coefficients, largest error ", tab.maxerr, ")"
)

"""
    table_rotation_jet(tab, et, Val(N)) -> NTuple{N+1,SMatrix{3,3}}

Return the rotation matrix from the rotation table `tab` at `et` and its time
derivatives to order `N`.

Element `k + 1` is the derivative of order `k`, in 1/s^k; every order is a derivative
of the polynomials of the table. Order 1 thus follows order 1 of [`pxder`](@ref), which
lists the frame classes for which the derivative block of [`sxform`](@ref) differs.
`et` is in seconds past J2000 TDB. The rules for the
epoch, for `@inbounds` and for the number types are those of [`table_jet`](@ref). The
call reads no kernel and does not allocate.

Throws `SpiceError(:spkinsuffdata)` if `et` is outside [`table_span`](@ref)`(tab)`.
"""
Base.@propagate_inbounds function table_rotation_jet(tab::RotationTable{T}, et::Real, v::Val{N}) where {T, N}
    i, s0, rr0 = _tab_locate(tab, et)
    s, rr = _tab_arg(T, s0), _tab_arg(T, rr0)
    a = 9 * tab.ncoef * i + 1
    c1 = _tab_scale(_cheb3(tab.coef, a, tab.ncoef, s, one(T), v), rr)
    c2 = _tab_scale(_cheb3(tab.coef, a + 3 * tab.ncoef, tab.ncoef, s, one(T), v), rr)
    c3 = _tab_scale(_cheb3(tab.coef, a + 6 * tab.ncoef, tab.ncoef, s, one(T), v), rr)
    return ntuple(k -> hcat(c1[k], c2[k], c3[k]), Val(N + 1))
end

"""
    table_rotation(tab, et) -> SMatrix{3,3}

Return the rotation matrix from frame `tab.from` to frame `tab.to` at `et`, from the
rotation table `tab`.

`et` is in seconds past J2000 TDB. The rules for the epoch, for `@inbounds` and for
the number types are those of [`table_jet`](@ref).

Throws `SpiceError(:spkinsuffdata)` if `et` is outside the table.
"""
Base.@propagate_inbounds table_rotation(tab::RotationTable, et::Real) =
    table_rotation_jet(tab, et, Val(0))[1]

"""
    rotation_table([ctx,] from, to, t0, t1; tol=nothing, ncoef=16, maxintervals=65536)
        -> RotationTable

Return a table of the rotation from frame `from` to frame `to` for `t0 <= et <= t1`,
which [`table_rotation`](@ref) and [`table_rotation_jet`](@ref) evaluate without the
kernels. This call has no CSPICE counterpart.

`t0` and `t1` are in seconds past J2000 TDB; the frames are as for [`pxform`](@ref).
The result is a `RotationTable{Float64,Vector{Float64}}` with
`table_span(result) == (t0, t1)` to rounding, and its field `maxerr` is the measured
error.

The call interpolates the nine elements of the matrix with `ncoef` coefficients for
each interval, and it doubles the number of intervals until the measured error
(Frobenius norm, no unit) is at most `tol`. With `tol`, the result has
`maxerr <= tol` or the call throws. Without `tol` the goal is `1.7e-12`, and the call
also stops when the error is below `1.7e-6` and one more doubling does not divide it
by four. `pxform` then limits the accuracy: the rotation angle of many turns has a
rounding error above the goal. Read `maxerr` for the error that was reached.

Throws `SpiceError(:badarg)` if the epochs are not finite with `t0 < t1`, if
`ncoef < 2`, `tol <= 0` or `maxintervals < 1`. Throws `SpiceError(:noconvergence)` if
`maxintervals` intervals do not reach the goal, and the errors of [`pxform`](@ref) if
the kernels do not cover `t0` to `t1`.

# Examples
```julia
furnsh("naif0012.tls"); furnsh("pck00011.tpc")
t0, t1 = str2et("2025-01-01"), str2et("2025-02-01")
tab = rotation_table("J2000", "IAU_MARS", t0, t1)
R = table_rotation(tab, t0 + 3600.0)
```

See also [`RotationTable`](@ref), [`ephemeris_table`](@ref).
"""
function rotation_table(
        ks::KernelSet, from::NameOrID, to::NameOrID, t0::Real, t1::Real;
        tol::Union{Nothing, Real} = nothing, ncoef::Integer = 16, maxintervals::Integer = 65536
    )
    a, b = Float64(t0), Float64(t1)
    (isfinite(a) && isfinite(b) && a < b) ||
        throw(SpiceError(:badarg, "rotation_table: the epochs must be finite, with t0 < t1"))
    ncoef >= 2 || throw(SpiceError(:badarg, "rotation_table: ncoef must be at least 2"))
    (tol === nothing || tol > 0) || throw(SpiceError(:badarg, "rotation_table: tol must be positive"))
    maxintervals >= 1 || throw(SpiceError(:badarg, "rotation_table: maxintervals must be at least 1"))
    f, t = frame_id(ks, from), frame_id(ks, to)
    source = et -> SVector{9, Float64}(frame_rotation(ks, f, t, et))
    source(a), source(b)       # data at both ends, or the error of the frame layer
    coef, dt, nint, maxerr = _tab_refine(
        source, a, b, tol, Int(ncoef), maxintervals, Val(3),
        "rotation_table", ""
    )
    tab = RotationTable(coef, a, dt, nint, ncoef, f, t, maxerr)
    @assert table_span(tab)[1] == a && abs(table_span(tab)[2] - b) <= 1.0e-9 * (b - a)
    return tab
end
rotation_table(ctx::SpiceContext, from::NameOrID, to::NameOrID, t0::Real, t1::Real; kw...) =
    rotation_table(kernelset(ctx), from, to, t0, t1; kw...)
rotation_table(from::NameOrID, to::NameOrID, t0::Real, t1::Real; kw...) =
    rotation_table(default_context(), from, to, t0, t1; kw...)
