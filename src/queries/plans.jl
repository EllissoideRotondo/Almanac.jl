# Chain plans of the prepared queries. `state_query` and `rotation_query` find once the
# SPK segment lists and the frame links that a chain passes through, when the chain does
# not depend on the epoch; the evaluation then skips the search of the chain and the
# lookups by body and frame ID. A plan makes the operations of `spk_jet` and
# `frame_rotation` in the same order, so the results are the same to the last bit, and
# it hands an epoch at which a link has no data to the general path, which throws or
# finds what those calls find.

# --- SPK chains ---------------------------------------------------------------------------

"""
    _SPKPlan

The SPK chain of the geometric state of `targ` relative to `obs` in frame `ref`.
`links[1:nt]` are the segment lists of the target chain up to the first body that the
chain of the observer meets, and the other elements those of the chain of the observer.
`frame` is the frame of every segment of the lists, or 0 if they have more than one.
"""
struct _SPKPlan
    ks::KernelSet
    targ::Int32
    ref::Int32
    obs::Int32
    links::Vector{SegmentList{Segment}}
    nt::Int
    frame::Int32
end

"The center of every segment of `segs`, or `nothing` if the segments have more than one."
function _one_center(segs::SegmentList{Segment})
    isempty(segs) && return nothing
    c = segs[1].center
    return all(s -> s.center == c, segs) ? c : nothing
end

"""
    _spk_plan(ks, targ::Int32, ref::Int32, obs::Int32, corr) -> Union{Nothing,_SPKPlan}

Plan of the geometric chain from `targ` to `obs`, or `nothing`. It exists when `corr` is geometric,
`targ != obs`, every body on the two chains has segments of a single center, the target
chain is shorter than `MAXCHAIN` bodies and the plan has at most `MAXCHAIN` links. The
two chains are then the same at every epoch, so the segments that `spk_jet` evaluates
at an epoch are those of the plan's lists that cover it.
"""
function _spk_plan(ks::KernelSet, targ::Int32, ref::Int32, obs::Int32, corr::Aberration)
    (corr.geometric && targ != obs) || return nothing
    bodies = Int32[targ]
    links = SegmentList{Segment}[]
    body = targ
    # The walk of pass 1 of `spk_jet`, without epoch.
    while body != 0 && body != obs
        # invariant: links[k] holds the segments from bodies[k] to bodies[k+1]
        length(bodies) == MAXCHAIN && return nothing     # `spk_jet_long`
        segs = get(ks.spk, body, nothing)
        segs === nothing && break
        c = _one_center(segs)
        (c === nothing || c in bodies) && return nothing
        push!(links, segs)
        body = c
        push!(bodies, body)
    end
    # The walk of pass 2: a chain of the observer that does not meet the target chain
    # throws at every epoch, from `spk_jet`.
    olinks = SegmentList{Segment}[]
    body = obs
    common = findfirst(==(body), bodies)
    while common === nothing
        # invariant: olinks lead from obs to body, and no body before body is in `bodies`
        body == 0 && return nothing
        segs = get(ks.spk, body, nothing)
        segs === nothing && return nothing
        c = _one_center(segs)
        (c === nothing || length(olinks) == MAXCHAIN) && return nothing
        push!(olinks, segs)
        body = c
        common = findfirst(==(body), bodies)
    end
    chain = vcat(links[1:(common - 1)], olinks)
    length(chain) <= MAXCHAIN || return nothing
    f = chain[1][1].frame
    frame = all(l -> all(s -> s.frame == f, l), chain) ? f : Int32(0)
    return _SPKPlan(ks, targ, ref, obs, chain, common - 1, frame)
end

"""
    _spk_plan_jet(p::_SPKPlan, et, ::Val{N}) -> NTuple{N+1,Vec3}

`spk_jet(p.ks, p.targ, et, p.ref, p.obs, Val(N))` from the plan, for a `Float64` epoch.
The segments of all links are found before any evaluation; at an epoch where a list has
no segment, `spk_jet` takes the call and throws or finds what it finds. The chain of the
observer and the target chain are then summed as in `spk_jet`. An epoch of another
number type goes to `spk_jet`.
"""
function _spk_plan_jet(p::_SPKPlan, et::Float64, v::Val{N}) where {N}
    if length(p.links) == 1 && p.frame == p.ref
        # One link in the frame of the query: the target relative to the center of its
        # segment less nothing, or nothing less the observer relative to the target.
        segs = @inbounds p.links[1]
        i = find_segment(segs, et)
        i == 0 && return spk_jet(p.ks, p.targ, et, p.ref, p.obs, v)
        jet = segment_jet_native(p.ks, @inbounds(segs[i]), et, v)
        return p.nt == 1 ? jet : map(-, ntuple(_ -> zero(Vec3{Float64}), Val(N + 1)), jet)
    end
    return _plan_links(p, et, v)
end
_spk_plan_jet(p::_SPKPlan, et::Real, v::Val) = spk_jet(p.ks, p.targ, et, p.ref, p.obs, v)

# A plan of more than one link, or of one link in another frame. Not inlined, so that
# the path of one link stays small.
@noinline function _plan_links(p::_SPKPlan, et::Float64, v::Val{N}) where {N}
    idx = ntuple(_ -> 0, Val(MAXCHAIN))
    for k in eachindex(p.links)
        # invariant: idx[1:k-1] are the segments of links 1 to k-1 that cover et
        i = find_segment(@inbounds(p.links[k]), et)
        i == 0 && return spk_jet(p.ks, p.targ, et, p.ref, p.obs, v)
        idx = Base.setindex(idx, i, k)
    end
    return p.frame != 0 ? _plan_one_frame(p, idx, et, v) : _plan_frames(p, idx, et, v)
end

@inline _plan_segment(p::_SPKPlan, idx, k::Int) = @inbounds p.links[k][idx[k]]

# All segments in the frame `p.frame`: every frame change of `spk_jet` but the last is
# from a frame to itself, and the sums are those of `_chain_first` and `_chain_add`. The
# links of the observer are evaluated first.
@inline function _plan_one_frame(p::_SPKPlan, idx, et::Float64, ::Val{N}) where {N}
    ks = p.ks
    nt = p.nt
    zero_jet = ntuple(_ -> zero(Vec3{Float64}), Val(N + 1))
    ojet = zero_jet
    for k in (nt + 1):length(p.links)
        # invariant: ojet sums the links of the observer before k
        jet = segment_jet_native(ks, _plan_segment(p, idx, k), et, Val(N))
        ojet = k == nt + 1 ? jet : map(+, ojet, jet)
    end
    tjet = zero_jet
    for k in 1:nt
        jet = segment_jet_native(ks, _plan_segment(p, idx, k), et, Val(N))
        tjet = k == 1 ? jet : map(+, tjet, jet)
    end
    jet = map(-, tjet, ojet)
    return p.frame == p.ref ? jet : _hop(ks, p.frame, p.ref, et, jet)
end

# Segments in more than one frame: passes 2 and 3 of `spk_jet` and `_chain_combine`.
@noinline function _plan_frames(p::_SPKPlan, idx, et::Float64, ::Val{N}) where {N}
    ks = p.ks
    nt = p.nt
    ojet = ntuple(_ -> zero(Vec3{Float64}), Val(N + 1))
    oframe = Int32(1)
    for k in (nt + 1):length(p.links)
        seg = _plan_segment(p, idx, k)
        ojet, oframe = k == nt + 1 ? _chain_first(ks, et, ojet, seg, Val(N)) :
            _chain_add(ks, et, ojet, oframe, seg, Val(N))
    end
    tjet = ntuple(_ -> zero(Vec3{Float64}), Val(N + 1))
    tframe = oframe
    for k in 1:nt
        seg = _plan_segment(p, idx, k)
        tjet, tframe = k == 1 ? _chain_first(ks, et, tjet, seg, Val(N)) :
            _chain_add(ks, et, tjet, tframe, seg, Val(N))
    end
    nt == length(p.links) && (oframe = tframe)
    return _chain_combine(ks, et, tjet, tframe, ojet, oframe, p.ref)
end

# --- frame chains -------------------------------------------------------------------------

"""
    _FrameLink

One link of a frame chain of a plan: frame `id` of class 1 or 4, whose matrix to its
parent is the constant `rot`, or of class 2, whose matrix comes from the binary PCK
segments `pck` of its class ID (`_NO_SEGMENTS` for none), or else from the text PCK
constants.
"""
struct _FrameLink
    id::Int32
    class::Int32
    classid::Int32
    pck::SegmentList{Segment}
    rot::Mat3{Float64}
end

"""
    _FramePlan

The frame chains of the rotation from `from` to `to`: `up` from `from` to J2000 or to
`to`, and `down` from `to` to the frame number `meet` of the chain of `from` (1 for
`from`; `meet == 0` when the chain of `from` reaches `to`). `state` is true when
`sxform` takes the state transformation of CSPICE (`_ck_frame_jet`).
"""
struct _FramePlan
    ks::KernelSet
    from::Int32
    to::Int32
    up::Vector{_FrameLink}
    down::Vector{_FrameLink}
    meet::Int
    state::Bool
end

"""
    _frame_link(ks, id) -> Union{Nothing,Tuple{Int32,_FrameLink}}

The parent of frame `id` and its link, for a built-in inertial frame, a PCK frame or a
fixed-offset frame with data; `nothing` for any other frame, whose parent can depend on
the epoch, and for a frame that `_frame_step` refuses.
"""
function _frame_link(ks::KernelSet, id::Int32)
    fd = _frame_def_or_nothing(ks, id)
    fd === nothing && return nothing
    if fd.class == 1 && 1 <= id <= N_INERTIAL_FRAMES
        rot = Mat3{Float64}(transpose(inertial_rotation(id)))
        return _J2000_FRAME, _FrameLink(id, 1, fd.classid, _NO_SEGMENTS, rot)
    elseif fd.class == 2
        pck = get(ks.pck, fd.classid, _NO_SEGMENTS)
        return _J2000_FRAME, _FrameLink(id, 2, fd.classid, pck, one(Mat3{Float64}))
    elseif fd.class == 4
        tk = get(ks.tk, id, nothing)
        tk === nothing && return nothing
        return tk.relative, _FrameLink(id, 4, fd.classid, _NO_SEGMENTS, tk.rot)
    end
    return nothing
end

"""
    _frame_plan(ks, from::Int32, to::Int32) -> Union{Nothing,_FramePlan}

Plan of the chains between the frames `from` and `to`, which `rotation_query` has
checked, or `nothing`. It exists when every frame on them is a built-in inertial frame, a PCK frame or
a fixed-offset frame, the chain of `from` has fewer than `MAX_FRAME_CHAIN` frames, and
the two frames are not equal or both built-in inertial frames (whose rotation is cheap
without a plan). The parents of those frames do not depend on the epoch.
"""
function _frame_plan(ks::KernelSet, from::Int32, to::Int32)
    (from == to || _both_inertial(from, to)) && return nothing
    nodes = Int32[from]
    up = _FrameLink[]
    cur = from
    while cur != _J2000_FRAME && cur != to
        # invariant: up[k] is the link from nodes[k] to nodes[k+1]
        length(nodes) == MAX_FRAME_CHAIN && return nothing   # `_chain_rotation` folds links
        l = _frame_link(ks, cur)
        l === nothing && return nothing
        cur, link = l
        cur in nodes && return nothing
        push!(up, link)
        push!(nodes, cur)
    end
    down = _FrameLink[]
    meet = 0
    if cur != to
        cur = to
        k = findfirst(==(cur), nodes)
        while k === nothing
            # invariant: down leads from to to cur, and no frame before cur is in `nodes`
            length(down) == _FRAME_CHAIN_LIMIT && return nothing
            l = _frame_link(ks, cur)
            l === nothing && return nothing
            cur, link = l
            push!(down, link)
            k = findfirst(==(cur), nodes)
        end
        meet = k
    end
    state = try
        !_exact_rates(ks) && (_ck_in_chain(ks, from) || _ck_in_chain(ks, to))
    catch e
        e isa SpiceError || rethrow()
        return nothing
    end
    return _FramePlan(ks, from, to, up, down, meet, state)
end

"Matrix of the link `l` in the number type of the mode, as `_chain_step` gives it."
@inline function _link_rotation(mode, ks::KernelSet, l::_FrameLink, et::T) where {T <: Real}
    l.class == 2 && return _link_body(mode, ks, l, et)
    return Mat3{_ck_numtype(mode, T)}(l.rot)
end
@inline _link_body(::_CKExact, ks::KernelSet, l::_FrameLink, et::T) where {T <: Real} =
    Mat3{T}(transpose(_body_rotation(ks, l.pck, l.classid, et)))
@inline _link_body(::_CKState, ks::KernelSet, l::_FrameLink, et::Float64) =
    Mat3{Taylor{1, Float64}}(transpose(_body_state(ks, l.pck, l.classid, et)))

"""
    _plan_chain(mode, p::_FramePlan, et) -> Mat3

`_chain_rotation(mode, p.ks, p.from, p.to, et)` from the plan: the products of the
links in the same order. As there, the chain of `from` is evaluated to its end, also
past the frame where the chain of `to` meets it.
"""
function _plan_chain(mode, p::_FramePlan, et::T) where {T <: Real}
    S = _ck_numtype(mode, T)
    ks = p.ks
    rot = one(Mat3{S})
    up = rot
    for k in eachindex(p.up)
        # invariant: rot is the matrix from `from` to the parent of link k - 1
        step = _link_rotation(mode, ks, @inbounds(p.up[k]), et)
        rot = k == 1 ? step : _ck_mxm(mode, step, rot)
        k + 1 == p.meet && (up = rot)
    end
    p.meet == 0 && return rot
    down = one(Mat3{S})
    for k in eachindex(p.down)
        step = _link_rotation(mode, ks, @inbounds(p.down[k]), et)
        down = k == 1 ? step : _ck_mxm(mode, step, down)
    end
    p.meet == 1 && return Mat3{S}(transpose(down))
    return Mat3{S}(_ck_mtxm(mode, down, up))
end

"""
    _plan_rotation(p::_FramePlan, et) -> Mat3

`frame_rotation(p.ks, p.from, p.to, et)` from the plan, for a `Float64` epoch; another
number type goes to `frame_rotation`.
"""
_plan_rotation(p::_FramePlan, et::Float64) = _plan_chain(_CKExact(), p, et)
_plan_rotation(p::_FramePlan, et::Real) = frame_rotation(p.ks, p.from, p.to, et)

"""
    _plan_frame_jet(p::_FramePlan, et, ::Val{N}) -> NTuple{N+1,Mat3}

`frame_jet(p.ks, p.from, p.to, et, Val(N))` from the plan, for order 1 and a `Float64`
epoch: the state transformation of `_ck_frame_jet` or the jet of `rotation_jet`, as the
plan's `state` says. Other orders and number types go to `frame_jet`.
"""
function _plan_frame_jet(p::_FramePlan, et::Float64, v::Val{1})
    if p.state
        x = _plan_chain(_CKState(), p, et)
        return map(e -> e.v, x), map(e -> e.d[1], x)
    end
    d = map(derivatives, _plan_chain(_CKExact(), p, taylor_var(et, v)))
    return ntuple(k -> map(x -> x[k], d), Val(2))
end
_plan_frame_jet(p::_FramePlan, et::Real, v::Val) = frame_jet(p.ks, p.from, p.to, et, v)
