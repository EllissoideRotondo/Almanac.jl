# Ephemeris interface: geometric states from the SPK segments of a `KernelSet`, with
# time derivatives of any order. The aberration corrections are in `aberration.jl`.

"Speed of light in km/s, the value of `clight_c`."
const CLIGHT = 299792.458

"Number of bodies of a target chain that `spk_jet` holds; a longer chain goes to `spk_jet_long`."
const MAXCHAIN = 10

"""
    rotate_jet(Rjet::NTuple{M}, rjet::NTuple{M}) -> NTuple{M}

Return the jet of the product `R(t) * r(t)` from the jets of `R` and `r` (ordinary
derivatives, orders 0 to M-1), by the Leibniz rule.
"""
@inline function rotate_jet(Rjet::NTuple{M, Any}, rjet::NTuple{M, Any}) where {M}
    return ntuple(Val(M)) do i
        k = i - 1
        acc = _m_mxv(Rjet[1], rjet[i])
        for j in 1:k
            acc += binomial(k, j) * _m_mxv(Rjet[j + 1], rjet[i - j])
        end
        acc
    end
end

# Order 1: each velocity component sums its six products from left to right, those of
# the derivative block first, as `mxvg_` applies the state transformation of `frmchg_`
# in `spkgeo_` and `spkez_`; the sum of the general method above rounds differently.
@inline function rotate_jet(Rjet::NTuple{2, Any}, rjet::NTuple{2, Any})
    R, dR = Rjet
    r, v = rjet
    w = SVector(
        ntuple(Val(3)) do i
            dR[i, 1] * r[1] + dR[i, 2] * r[2] + dR[i, 3] * r[3] +
                R[i, 1] * v[1] + R[i, 2] * v[2] + R[i, 3] * v[3]
        end
    )
    return (_m_mxv(R, r), w)
end

"""
    frame_jet(ks, from::Int32, to::Int32, et, ::Val{N}) -> NTuple{N+1,Mat3}

Return the rotation matrix from frame `from` to frame `to` at `et` and its time
derivatives to order N, for any epoch number type: `rotation_jet` for a plain epoch,
`frame_rotation` of a `Taylor` epoch otherwise.
"""
@inline frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Float64, v::Val) =
    rotation_jet(ks, from, to, et, v)
# The method for `Float64` and `Val{0}` resolves the ambiguity of the two above and below.
@inline frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Float64, v::Val{0}) =
    rotation_jet(ks, from, to, et, v)
@inline frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Real, ::Val{0}) =
    (frame_rotation(ks, from, to, et),)
@inline function frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Real, v::Val{N}) where {N}
    d = map(derivatives, frame_rotation(ks, from, to, taylor_var(et, v)))
    return ntuple(k -> map(x -> x[k], d), Val(N + 1))
end

"Jet of order N of `seg` at `et` in the segment frame `seg.frame`."
@inline function segment_jet_native(ks::KernelSet, seg::Segment, et::T, ::Val{N}) where {T <: Real, N}
    seg.dtype == PROVIDER_DTYPE &&
        return provider_jet(ks.loaded[seg.baddr].content::EphemerisProvider, et, Val(N))
    seg.dtype == Int32(10) && return _t10_jet(_lsk(ks), ks.daf[seg.file], seg, et, Val(N))
    return spk_eval(ks.daf[seg.file], seg, et, Val(N))
end

"`jet` from frame `from` to frame `to` at `et`; `jet` itself if the frames are equal."
@inline _hop(ks::KernelSet, from::Int32, to::Int32, et, jet::NTuple{M, Any}) where {M} =
    from == to ? jet : rotate_jet(frame_jet(ks, from, to, et, Val(M - 1)), jet)

"""
    _chain_add(ks, et, acc, accframe, seg, ::Val{N}) -> (acc, frame)

Add the jet of `seg` to the sum `acc` in frame `accframe` as `SPKGEO` of CSPICE does:
in the segment frame, after `acc` is rotated to it when the frames differ. With this
order a chain of segments in a frame other than J2000 gives the CSPICE state to the
last bit.
"""
@inline function _chain_add(ks::KernelSet, et, acc, accframe::Int32, seg::Segment, v::Val)
    jet = segment_jet_native(ks, seg, et, v)
    return map(+, _hop(ks, accframe, seg.frame, et, acc), jet), seg.frame
end

"""
    _chain_combine(ks, et, tjet, tframe, ojet, oframe, ref) -> jet

`tjet` (frame `tframe`) minus `ojet` (frame `oframe`) in frame `ref`, with the frame
changes of `SPKGEO` of CSPICE.
"""
@inline function _chain_combine(
        ks::KernelSet, et, tjet, tframe::Int32, ojet, oframe::Int32,
        ref::Int32
    )
    if tframe == oframe
        jet = map(-, tjet, ojet)
    elseif tframe == ref
        jet = map(-, tjet, _hop(ks, oframe, ref, et, ojet))
        oframe = ref
    else
        jet = map(-, _hop(ks, tframe, oframe, et, tjet), ojet)
    end
    return _hop(ks, oframe, ref, et, jet)
end

@inline function chain_index(chain::NTuple{MAXCHAIN, Int32}, n::Int, body::Int32)
    for k in 1:n
        chain[k] == body && return k
    end
    return 0
end

"Text that lists the SPK coverage intervals of `body` in `ks` (ET seconds)."
function coverage_text(ks::KernelSet, body::Integer)
    segs = get(ks.spk, Int32(body), nothing)
    segs === nothing && return "no SPK segment is loaded for this body"
    parts = ["[$(s.t0), $(s.t1)] (center $(s.center))" for s in Iterators.take(segs, 8)]
    more = length(segs) > 8 ? " and $(length(segs) - 8) more" : ""
    return "available coverage: " * join(parts, ", ") * more
end

"""
    throw_no_connection(ks, targ, obs, t, tend, oend) -> Union{}

Throw the [`SpkInsufficientDataError`](@ref) for two chains that do not meet at `t`;
`tend` and `oend` are the last bodies of the target chain and of the observer chain.
CSPICE signals NOLOADEDFILES with no SPK file and SPKINSUFFDATA for every other case.
"""
@noinline function throw_no_connection(ks::KernelSet, targ, obs, t::Float64, tend, oend)
    isempty(ks.spk) && throw(
        SpkInsufficientDataError(
            :noloadedfiles, "no SPK file is " *
                "loaded (needed for target $targ relative to observer $obs at ET $t s); load " *
                "one with furnsh", targ, obs, t, nothing, false
        )
    )
    for body in (oend, tend)
        if body != 0 && haskey(ks.spk, Int32(body))
            throw(
                SpkInsufficientDataError(
                    :spkinsuffdata, "no SPK segment covers body " *
                        "$body at ET $t s (needed for target $targ relative to observer $obs); " *
                        coverage_text(ks, body), targ, obs, t, body, true
                )
            )
        end
    end
    # No end of a chain is in the files at another epoch: the end that is not the solar
    # system barycenter is the body without data.
    missing_body = tend != 0 ? tend : oend != 0 ? oend : nothing
    throw(
        SpkInsufficientDataError(
            :spkinsuffdata, "no chain of SPK segments connects " *
                "target $targ and observer $obs at ET $t s: the target chain ends at body $tend " *
                "and the observer chain ends at body $oend", targ, obs, t, missing_body, false
        )
    )
end

# CSPICE does not return on a cycle; `spk_jet` stops after as many steps as there are
# bodies with SPK data, which no chain without a cycle exceeds.
@noinline throw_chain_cycle(targ, obs, t::Float64) = throw(
    SpkInsufficientDataError(
        :spkinsuffdata, "the SPK chain between target $targ and observer $obs at ET $t s " *
            "has a cycle", targ, obs, t, nothing, false
    )
)

"""
    spk_jet(ks, targ::Int32, et, ref::Int32, obs::Int32, ::Val{N}) -> NTuple{N+1,Vec3}

Return the geometric position of `targ` relative to `obs` in frame `ref` at `et` and
its time derivatives to order N (no aberration correction). The function does not
allocate, except for a target chain of more than MAXCHAIN bodies. Throws
[`SpkInsufficientDataError`](@ref).
"""
function spk_jet(
        ks::KernelSet, targ::Int32, et::T, ref::Int32, obs::Int32,
        ::Val{N}
    ) where {T <: Real, N}
    t = primal(et)
    # Pass 1: bodies and segment indices of the target chain up to the observer or the
    # root, without evaluation.
    tbody = ntuple(_ -> Int32(0), Val(MAXCHAIN))
    tseg = ntuple(_ -> 0, Val(MAXCHAIN))
    tbody = Base.setindex(tbody, targ, 1)
    nt = 1
    body = targ
    while body != 0 && body != obs
        # invariant: tbody[1:nt] are the bodies reached, tseg[1:nt-1] their segments
        nt == MAXCHAIN && return spk_jet_long(ks, targ, et, ref, obs, Val(N))
        segs = get(ks.spk, body, nothing)
        segs === nothing && break
        i = find_segment(segs, t)
        i == 0 && break
        tseg = Base.setindex(tseg, i, nt)
        body = segs[i].center
        nt += 1
        tbody = Base.setindex(tbody, body, nt)
    end
    # Pass 2: walk the observer chain until it meets the target chain.
    ojet = ntuple(_ -> zero(Vec3{T}), Val(N + 1))
    oframe = Int32(1)
    body = obs
    depth = 0
    common = chain_index(tbody, nt, body)
    while common == 0
        # `spkgeo` ends the observer chain at the barycenter (0), even with a segment
        # whose target is 0.
        body == 0 && throw_no_connection(ks, targ, obs, t, tbody[nt], body)
        segs = get(ks.spk, body, nothing)
        segs === nothing && throw_no_connection(ks, targ, obs, t, tbody[nt], body)
        i = find_segment(segs, t)
        i == 0 && throw_no_connection(ks, targ, obs, t, tbody[nt], body)
        depth += 1
        depth > length(ks.spk) && throw_chain_cycle(targ, obs, t)
        seg = segs[i]
        depth == 1 && (oframe = seg.frame)
        ojet, oframe = _chain_add(ks, et, ojet, oframe, seg, Val(N))
        body = seg.center
        common = chain_index(tbody, nt, body)
    end
    # Pass 3: target chain up to the common body, summed as the observer chain.
    tjet = ntuple(_ -> zero(Vec3{T}), Val(N + 1))
    tframe = oframe
    for k in 1:(common - 1)
        seg = ks.spk[tbody[k]][tseg[k]]
        k == 1 && (tframe = seg.frame)
        tjet, tframe = _chain_add(ks, et, tjet, tframe, seg, Val(N))
    end
    depth == 0 && (oframe = tframe)   # no observer segment: CSPICE takes the target frame
    return _chain_combine(ks, et, tjet, tframe, ojet, oframe, ref)
end

"""
    spk_jet_long(ks, targ::Int32, et, ref::Int32, obs::Int32, ::Val{N}) -> NTuple{N+1,Vec3}

`spk_jet` for a target chain of more than MAXCHAIN bodies; it allocates. As in
`SPKGEO`, the length of a chain has no limit.
"""
@noinline function spk_jet_long(
        ks::KernelSet, targ::Int32, et::T, ref::Int32, obs::Int32,
        ::Val{N}
    ) where {T <: Real, N}
    t = primal(et)
    tbody = Int32[targ]
    tseg = Segment[]
    body = targ
    while body != 0 && body != obs
        # invariant: tseg[k] is the segment from tbody[k] to tbody[k+1]
        segs = get(ks.spk, body, nothing)
        segs === nothing && break
        i = find_segment(segs, t)
        i == 0 && break
        length(tseg) < length(ks.spk) || throw_chain_cycle(targ, obs, t)
        push!(tseg, segs[i])
        body = segs[i].center
        push!(tbody, body)
    end
    ojet = ntuple(_ -> zero(Vec3{T}), Val(N + 1))
    oframe = Int32(1)
    body = obs
    depth = 0
    common = findfirst(==(body), tbody)
    while common === nothing
        body == 0 && throw_no_connection(ks, targ, obs, t, last(tbody), body)
        segs = get(ks.spk, body, nothing)
        segs === nothing && throw_no_connection(ks, targ, obs, t, last(tbody), body)
        i = find_segment(segs, t)
        i == 0 && throw_no_connection(ks, targ, obs, t, last(tbody), body)
        depth += 1
        depth > length(ks.spk) && throw_chain_cycle(targ, obs, t)
        depth == 1 && (oframe = segs[i].frame)
        ojet, oframe = _chain_add(ks, et, ojet, oframe, segs[i], Val(N))
        body = segs[i].center
        common = findfirst(==(body), tbody)
    end
    tjet = ntuple(_ -> zero(Vec3{T}), Val(N + 1))
    tframe = oframe
    for k in 1:(common - 1)
        k == 1 && (tframe = tseg[k].frame)
        tjet, tframe = _chain_add(ks, et, tjet, tframe, tseg[k], Val(N))
    end
    depth == 0 && (oframe = tframe)   # no observer segment: CSPICE takes the target frame
    return _chain_combine(ks, et, tjet, tframe, ojet, oframe, ref)
end

@inline function check_epoch(et::Real)
    isfinite(primal(et)) || throw(SpiceError(:badarg, "the epoch is not finite"))
    return float(et)
end

"""
    _plain_epoch(et, who) -> et

Return `et` unchanged if it is a plain real number. A call whose result has no
derivatives (`limbpt`, `occult`, `et2lst`, ...) takes such an epoch only; throws
`SpiceError(:badarg)` for a number that carries derivatives (a `ForwardDiff.Dual`, a
`Taylor`), which would otherwise fail later with a `MethodError`.
"""
function _plain_epoch(et::Real, who::AbstractString)
    (!(et isa Taylor) && float(et) isa AbstractFloat) || throw(
        SpiceError(
            :badarg,
            "$who: the epoch must be a plain real number such as a Float64, not a " *
                "$(nameof(typeof(et))); this call does not give derivatives"
        )
    )
    return et
end

"A body or a frame given by name or by integer ID."
const NameOrID = Union{AbstractString, Integer}

function spkder(
        ks::KernelSet, targ::NameOrID, et::Real, ref::NameOrID, obs::NameOrID,
        ::Val{N}
    ) where {N}
    N isa Int && N >= 0 || throw(SpiceError(:badarg, "spkder: the order must be an integer >= 0"))
    return spk_jet(
        ks, body_code(ks, targ), check_epoch(et), frame_id(ks, ref),
        body_code(ks, obs), Val(N)
    )
end
spkder(ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID, obs::NameOrID, v::Val) =
    spkder(kernelset(ctx), targ, et, ref, obs, v)
spkder(targ::NameOrID, et::Real, ref::NameOrID, obs::NameOrID, v::Val) =
    spkder(default_context(), targ, et, ref, obs, v)

"""
    spkgeo([ctx,] targ::Integer, et, ref, obs::Integer) -> (state, lt)

Return the geometric state of body `targ` relative to body `obs` in frame `ref` at
`et`, and the one-way light time, as
`spkgeo_c`.

`state` is an `SVector{6}`: position in km, then velocity in km/s. `lt` is the
distance divided by the speed of light, in seconds. `et` is in seconds past J2000
TDB. `targ` and `obs` are NAIF ID codes; `ref` is a frame name or a frame ID.

Throws `SpiceError(:unknownframe)` if `ref` is not a known frame name,
`SpiceError(:badarg)` if `et` is not finite, and [`SpkInsufficientDataError`](@ref) if
no SPK file is loaded (`:noloadedfiles`) or the loaded SPK files do not connect `targ`
and `obs` at `et` (`:spkinsuffdata`). As in CSPICE, `targ == obs` gives a zero state
and a zero light time before `ref` and `et` are checked; so does [`spkezr`](@ref) with
the correction `"NONE"`.

See also [`spkezr`](@ref), [`spkssb`](@ref).

# References

- [NAIF documentation of `spkgeo_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkgeo_c.html)
"""
function spkgeo(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID, obs::Integer)
    # As in CSPICE, a body relative to itself is zero before `ref` and `et` are checked.
    targ == obs && return zero(SVector{6, typeof(float(et))}), zero(float(et))
    r, v = spk_jet(ks, Int32(targ), check_epoch(et), frame_id(ks, ref), Int32(obs), Val(1))
    return vcat(r, v), _m_norm(r) / CLIGHT
end
spkgeo(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID, obs::Integer) =
    spkgeo(kernelset(ctx), targ, et, ref, obs)
spkgeo(targ::Integer, et::Real, ref::NameOrID, obs::Integer) =
    spkgeo(default_context(), targ, et, ref, obs)

"""
    spkgps([ctx,] targ::Integer, et, ref, obs::Integer) -> (pos, lt)

Return the geometric position of body `targ` relative to body `obs` in frame `ref` at
`et`, and the one-way light time, as
`spkgps_c`.

`pos` is an `SVector{3}` in km and `lt` is in seconds. The arguments and the errors are
those of [`spkgeo`](@ref). As in CSPICE, `targ == obs` gives a zero position and a zero
light time before `ref` and `et` are checked.

See also [`spkezp`](@ref), [`spkpos`](@ref).

# References

- [NAIF documentation of `spkgps_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkgps_c.html)
"""
function spkgps(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID, obs::Integer)
    targ == obs && return zero(Vec3{typeof(float(et))}), zero(float(et))
    (r,) = spk_jet(ks, Int32(targ), check_epoch(et), frame_id(ks, ref), Int32(obs), Val(0))
    return r, _m_norm(r) / CLIGHT
end
spkgps(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID, obs::Integer) =
    spkgps(kernelset(ctx), targ, et, ref, obs)
spkgps(targ::Integer, et::Real, ref::NameOrID, obs::Integer) =
    spkgps(default_context(), targ, et, ref, obs)

"""
    spkezr([ctx,] targ, et, ref, abcorr, obs) -> (state, lt)

Return the state of `targ` relative to `obs` in frame `ref` at `et`, and the one-way
light time, as
`spkezr_c`.

`state` is an `SVector{6}`: position in km, then velocity in km/s. `lt` is in seconds.

- `targ`, `obs`: body name, NAIF ID code, or a string that holds a code (`"-77"`).
- `et`: epoch of the observer in seconds past J2000 TDB (see [`str2et`](@ref)).
- `ref`: frame name or frame ID. For a frame that is not inertial and a correction
  other than `"NONE"`, the frame is evaluated at `et` corrected for the light time
  between the observer and the frame center, as in CSPICE.
- `abcorr`: `"NONE"` (geometric), `"LT"`, `"LT+S"`, `"CN"`, `"CN+S"` for reception, or
  the same with the prefix `X` for transmission. Letter case and blanks are ignored.
  An [`Aberration`](@ref) value is also accepted.

The call does not allocate, with integer IDs or with names in any case, unless it throws or
a chain has more than 10 SPK segments. With ForwardDiff loaded, `et`
can be a dual number (see [Differentiate with ForwardDiff](@ref)).

The SPK segments connect a body to its center, and that center to the next one. The
call needs such a chain from `targ` and from `obs` to a common body.

Throws `SpiceError`:
- `:idcodenotfound` if a body name is not known, `:unknownframe` if a frame name is
  not known;
- `:invalidoption` if `abcorr` is not valid, `:badarg` if `et` is not finite;
- `:noloadedfiles` if no SPK file is loaded;
- `:spkinsuffdata` if the loaded SPK files hold no chain between `targ` and `obs` at
  the epoch needed: a body is in no file, or the epoch is outside its segments;
- the kinds of [`pxform`](@ref) if `ref` is not inertial and its orientation data are
  not loaded.

`:noloadedfiles` and `:spkinsuffdata` are a [`SpkInsufficientDataError`](@ref), which
holds the body and the epoch.

# Examples
```julia
furnsh("naif0012.tls")
furnsh("de440s.bsp")
et = str2et("2025-03-01T12:00:00")
state, lt = spkezr("MARS BARYCENTER", et, "ECLIPJ2000", "LT+S", "EARTH")
r, v = state[1:3], state[4:6]                  # km, km/s
state, lt = spkezr(4, et, 17, Aberration("LT+S"), 399)   # the same with IDs
```

See also [`spkpos`](@ref), [`spkgeo`](@ref), [`spkder`](@ref), [`furnsh`](@ref).

# References

- [NAIF documentation of `spkezr_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html)
"""
spkezr(
    ks::KernelSet, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString,
    obs::NameOrID
) = spkezr(ks, targ, et, ref, Aberration(abcorr), obs)
spkezr(
    ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString,
    obs::NameOrID
) = spkezr(kernelset(ctx), targ, et, ref, Aberration(abcorr), obs)
spkezr(targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString, obs::NameOrID) =
    spkezr(default_context(), targ, et, ref, abcorr, obs)

"""
    spkpos([ctx,] targ, et, ref, abcorr, obs) -> (pos, lt)

Return the position of `targ` relative to `obs` in frame `ref` at `et`, and the
one-way light time, as
`spkpos_c`.

`pos` is an `SVector{3}` in km and `lt` is in seconds. The arguments and the errors
are those of [`spkezr`](@ref), except that an `abcorr` that is not valid throws
`SpiceError(:spkinvalidoption)`, as in CSPICE.

As in CSPICE, the `CN` corrections make three light-time iterations here;
[`spkezr`](@ref) makes up to five, with the convergence test of [`spkltc`](@ref). The
position can therefore differ from the position part of `spkezr` by up to 2e-11
relative (see [Known limits of accuracy](@ref)).

# Examples
```julia
et = str2et("2025-03-01T12:00:00")
pos, lt = spkpos("SUN", et, "J2000", "LT+S", "EARTH")    # km, s
```

See also [`spkezr`](@ref), [`spkder`](@ref).

# References

- [NAIF documentation of `spkpos_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkpos_c.html)
"""
spkpos(
    ks::KernelSet, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString,
    obs::NameOrID
) = spkpos(ks, targ, et, ref, _spk_abcorr(abcorr), obs)
spkpos(
    ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString,
    obs::NameOrID
) = spkpos(kernelset(ctx), targ, et, ref, _spk_abcorr(abcorr), obs)
spkpos(targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString, obs::NameOrID) =
    spkpos(default_context(), targ, et, ref, abcorr, obs)
