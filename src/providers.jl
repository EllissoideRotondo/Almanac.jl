# Bodies whose states come from functions of the program (ephemeris providers). A
# provider enters the load list as a file with one SPK segment of the data type
# `PROVIDER_DTYPE`, so every call that finds SPK segments (states, corrections, the
# event searches, tables) uses it, with the priority of a file loaded at that moment.

"""
    provider_jet(p::EphemerisProvider, et, ::Val{N}) -> NTuple{N+1,Vec3}

Contract: the jet of the position of `p.body` relative to `p.center` in `p.frame`. Order
0 and 1 come from the state that `p.f(et)` returns; order `k >= 2` is derivative `k-1`
of its velocity, from `p.f` called on a Taylor number (as for SPK types that store the
velocity, `spk_eval`). Throws `SpiceError(:badarraysize)` if `p.f` does not return six
numbers.
"""
@noinline function provider_jet(p::EphemerisProvider, et::T, ::Val{N}) where {T <: Real, N}
    # `p.f` is not inferred: every value is converted to `T` with an assertion, so that
    # the jet has a concrete type and the callers of `segment_jet` stay type-stable.
    # `@noinline` keeps this branch out of the inlined `segment_jet` of the file types.
    if N <= 1
        s = p.f(et)
        _provider_check(p, s, et)
        r = Vec3{T}(_tc(T, s[1]), _tc(T, s[2]), _tc(T, s[3]))
        v = Vec3{T}(_tc(T, s[4]), _tc(T, s[5]), _tc(T, s[6]))
        return ntuple(k -> k == 1 ? r : v, Val(N + 1))::NTuple{N + 1, Vec3{T}}
    end
    s = p.f(taylor_var(et, Val(N - 1)))
    _provider_check(p, s, et)
    # a name of its own: a variable that both branches assign is boxed in the closures
    r0 = Vec3{T}(_tc(T, _tvalue(s[1])), _tc(T, _tvalue(s[2])), _tc(T, _tvalue(s[3])))
    dx = _tderivs(T, s[4], Val(N - 1))
    dy = _tderivs(T, s[5], Val(N - 1))
    dz = _tderivs(T, s[6], Val(N - 1))
    return ntuple(k -> k == 1 ? r0 : Vec3{T}(dx[k - 1], dy[k - 1], dz[k - 1]), Val(N + 1))::NTuple{N + 1, Vec3{T}}
end

@inline _tc(::Type{T}, x) where {T} = convert(T, x)::T
_tvalue(x::Taylor) = x.v
_tvalue(x) = x
_tderivs(::Type{T}, x::Taylor, ::Val{M}) where {T, M} = map(y -> _tc(T, y), derivatives(x))::NTuple{M + 1, T}
_tderivs(::Type{T}, x, ::Val{M}) where {T, M} = ntuple(k -> k == 1 ? _tc(T, x) : zero(T), Val(M + 1))

@noinline function _provider_check(p::EphemerisProvider, s, et)
    (applicable(length, s) && length(s) == 6) || throw(
        SpiceError(
            :badarraysize,
            "the provider '$(p.name)' of body $(p.body) returned $(repr(s)) at ET " *
                "$(primal(et)) s; it must return the six numbers of a state"
        )
    )
    return nothing
end

"""
    _kernel_bodies(who, bodies) -> Vector{Tuple{String,Int32}}

Contract: the pairs `name => code` of the iterable `bodies` as a vector of name and code,
with the blanks around each name removed. Throws `SpiceError(:badarg)` for a blank name,
as [`boddef`](@ref) does, and for a code that is not a 32-bit integer.
"""
function _kernel_bodies(who::AbstractString, bodies)
    out = Tuple{String, Int32}[]
    for p in bodies
        nm = String(strip(first(p)))
        isempty(nm) && throw(SpiceError(:badarg, "$who: a body name of `bodies` is blank"))
        code = last(p)
        (code isa Integer && typemin(Int32) <= code <= typemax(Int32)) || throw(
            SpiceError(:badarg, "$who: the code $(repr(code)) of $(repr(nm)) is not a 32-bit integer")
        )
        push!(out, (nm, Int32(code)))
    end
    return out
end

"""
    _resolve_body(ks, x, bodies) -> Int32

Contract: the code of the body `x`: a name of the pairs `bodies` (the later pair first),
else the code that [`body_code`](@ref) gives in `ks`.
"""
function _resolve_body(ks::KernelSet, x::NameOrID, bodies::Vector{Tuple{String, Int32}})
    if x isa AbstractString
        n = normalize_name(x)
        for i in length(bodies):-1:1
            normalize_name(bodies[i][1]) == n && return bodies[i][2]
        end
    end
    return body_code(ks, x)
end

"""
    add_ephemeris!(f, [ctx,] body, center, frame, t0, t1; name="", bodies=()) -> Nothing

Make `body` a body whose state relative to `center` in `frame` is `f(et)` for
`t0 <= et <= t1` (ET, seconds past J2000 TDB). `f` returns six numbers: the position
(km) and the velocity (km/s).

The body is then available to every call that reads SPK data: [`spkezr`](@ref) with
any aberration correction, [`spkpos`](@ref), the geometry calls, the event searches,
[`ephemeris_table`](@ref). It enters the load list as a file named `name` with one SPK
segment, so the rules of a file hold: it has priority over the files loaded before it
and not over those loaded after it, outside `[t0, t1]` the files that cover the body
apply, [`unload`](@ref) with `name` removes it and [`kclear`](@ref) removes it with
the rest. [`kdata`](@ref) lists it with the type `"SPK"`. An empty `name` stands for
`"ephemeris of body <code>"`, with the NAIF code of `body`.

`bodies`, pairs such as `["MYSAT" => -1001]`, are body names that hold while the
provider is loaded, as for [`furnsh`](@ref) of a handle; `body` and `center` may be
one of them.

`f` is called from any thread that makes a query, so it must not change shared data. It
receives the epoch type of the query: `Float64`, a `ForwardDiff.Dual` for a dual-number
epoch, and for a derivative of order 2 or more (the acceleration that an `LT+S` state
of an observer needs, [`spkder`](@ref)) a Taylor number of this package, on which `f`
may use arithmetic and the functions of `Base` but no comparison. `f` may query `ctx`
for other bodies, but not for `body` in the interval it covers.

`body` and `center` are NAIF codes or names, and `frame` a frame name or code, all
known to `ctx` at the call. This call has no CSPICE counterpart.

Throws `SpiceError`:

- `:idcodenotfound` for an unknown body or center name;
- `:unknownframe` for an unknown frame;
- `:baddescrtimes` if `t1 < t0`;
- `:barycenterephem` if `body` and `center` are the same body.

A query throws `SpiceError(:badarraysize)` if `f` returns something other than six
numbers.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> add_ephemeris!(ctx, -1000, 399, "J2000", 0.0, 3600.0) do et
           [7000.0 + 0.001et, 0.0, 0.0, 0.001, 0.0, 0.0]
       end

julia> spkezr(ctx, -1000, 60.0, "J2000", "NONE", 399)[1] ≈ [7000.06, 0, 0, 0.001, 0, 0]
true
```

See also [`unload`](@ref), [`ephemeris_table`](@ref).
"""
function add_ephemeris!(
        f, ctx::SpiceContext, body::NameOrID, center::NameOrID,
        frame::NameOrID, t0::Real, t1::Real; name::AbstractString = "", bodies = ()
    )
    t0 <= t1 || throw(
        SpiceError(
            :baddescrtimes,
            "add_ephemeris!: the interval [$t0, $t1] ends before it starts"
        )
    )
    kb = _kernel_bodies("add_ephemeris!", bodies)
    lock(ctx.lock) do
        ks = kernelset(ctx)
        b, c = _resolve_body(ks, body, kb), _resolve_body(ks, center, kb)
        fr = frame_id(ks, frame)
        b == c && throw(
            SpiceError(
                :barycenterephem,
                "add_ephemeris!: body $b cannot be its own center"
            )
        )
        n = isempty(name) ? "ephemeris of body $b" : String(name)
        p = EphemerisProvider(n, b, c, fr, Float64(t0), Float64(t1), f, :SPK, kb)
        list = copy(ks.loaded)
        push!(list, Loaded(n, "", p))
        _pool_x_swap_loaded!(ctx, list, Loaded[])
    end
    return nothing
end
add_ephemeris!(
    f, body::NameOrID, center::NameOrID, frame::NameOrID, t0::Real, t1::Real;
    kw...
) = add_ephemeris!(f, default_context(), body, center, frame, t0, t1; kw...)

# ---------------------------------------------------------------------------------------
# Frames
# ---------------------------------------------------------------------------------------

"""
    frame_provider(ks, fd::FrameDef) -> Union{Nothing,FrameProvider}

Contract: the provider of the class 7 frame `fd`, or `nothing` if `fd` comes from the
kernel pool (a frame kernel may give any class; CSPICE has no class 7).
"""
@inline function frame_provider(ks::KernelSet, fd::FrameDef)
    i = fd.classid
    1 <= i <= length(ks.loaded) || return nothing
    c = ks.loaded[i].content
    return (c isa FrameProvider && c.id == fd.id) ? c : nothing
end

"""
    provider_rotation(p::FrameProvider, et::T) -> Mat3{T}

Contract: the matrix from frame `p.id` to frame `p.relative` that `p.f(et)` returns.
Throws `SpiceError(:badarraysize)` if it does not have 3 x 3 elements.
"""
@inline function provider_rotation(p::FrameProvider, et::T) where {T <: Real}
    m = p.f(et)
    (applicable(size, m) && size(m) == (3, 3)) || throw(
        SpiceError(
            :badarraysize,
            "the provider of frame '$(p.name)' returned $(repr(m)) at ET $(primal(et)) s; " *
                "it must return a 3 x 3 rotation matrix"
        )
    )
    return Mat3{T}(ntuple(k -> _tc(T, m[k]), Val(9)))::Mat3{T}
end

"""
    add_frame!(f, [ctx,] name, id, center, relative; path="") -> Nothing

Define the frame `name` with the frame ID `id`, centered on the body `center`, whose
orientation relative to the frame `relative` is `f(et)`: `f` returns the 3 x 3 matrix
that takes the coordinates of a vector in frame `name` to its coordinates in frame
`relative` at `et` (ET, seconds past J2000 TDB), the direction of the matrix of
[`pxform`](@ref)`(name, relative, et)`.

The frame is then available to every call that takes a frame: [`pxform`](@ref),
[`sxform`](@ref), states in the frame, [`add_ephemeris!`](@ref) and other frames
relative to it, [`frinfo`](@ref) (class 7, a class that CSPICE does not have). It
enters the load list as a frame kernel named `path` (`"frame <name>"` if `path` is
empty): [`unload`](@ref) with `path` removes it and [`kclear`](@ref) removes it with the rest. A frame of a text kernel with
the same name or ID loses to it, whatever the load order.

`f` is called from any thread that makes a query, so it must not change shared data. It
receives the epoch type of the query: `Float64`, a `ForwardDiff.Dual`, and for the
derivative of the matrix that [`sxform`](@ref) and corrected states in the frame need, a
Taylor number of this package, on which `f` may use arithmetic and the functions of
`Base` (`sin`, `cos`, ...) but no comparison. The matrix is used as returned; it is not
made orthogonal. `f` may query `ctx` for other frames and bodies, but not for this
frame.

`center` is a NAIF code or name and `relative` a frame name or ID, known to `ctx` at
the call. This call has no CSPICE counterpart.

Throws `SpiceError`:

- `:emptystring` if `name` is empty or blank;
- `:idcodenotfound` for an unknown center name;
- `:unknownframe` for an unknown relative frame;
- `:badframe` if `name` or `id` is a frame built into the package, or `id` is the ID
  of `relative`.

A query throws `SpiceError(:badarraysize)` if `f` returns something other than a 3 x 3
matrix.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> add_frame!(ctx, "SPIN", 1_400_001, 399, "J2000") do et
           s, c = sincos(1e-3 * et)
           [c -s 0; s c 0; 0 0 1]
       end

julia> pxform(ctx, "SPIN", "J2000", 0.0) == [1 0 0; 0 1 0; 0 0 1]
true

julia> sxform(ctx, "SPIN", "J2000", 0.0)[5, 1]    # derivative of the (2, 1) element, 1/s
0.001
```

See also [`add_ephemeris!`](@ref), [`unload`](@ref).
"""
function add_frame!(
        f, ctx::SpiceContext, name::AbstractString, id::Integer,
        center::NameOrID, relative::NameOrID; path::AbstractString = ""
    )
    nn = normalize_name(name)
    isempty(nn) && throw(SpiceError(:emptystring, "add_frame!: the frame name is empty"))
    (haskey(BUILTIN_FRAMENAME2ID, nn) || haskey(BUILTIN_FRAMES, Int32(id))) &&
        throw(SpiceError(:badframe, "add_frame!: '$name' or ID $id is a built-in frame"))
    lock(ctx.lock) do
        ks = kernelset(ctx)
        c, r = body_code(ks, center), frame_id(ks, relative)
        r == Int32(id) && throw(
            SpiceError(
                :badframe,
                "add_frame!: frame $id cannot be relative to itself"
            )
        )
        frame_def(ks, r)                    # the relative frame must be defined
        p = isempty(path) ? "frame $name" : String(path)
        list = copy(ks.loaded)
        push!(list, Loaded(p, "", FrameProvider(String(name), p, Int32(id), c, r, f, :TEXT)))
        _pool_x_swap_loaded!(ctx, list, Loaded[])
    end
    return nothing
end
add_frame!(f, name::AbstractString, id::Integer, center::NameOrID, relative::NameOrID; kw...) =
    add_frame!(f, default_context(), name, id, center, relative; kw...)
