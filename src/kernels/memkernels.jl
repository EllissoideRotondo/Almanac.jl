# Kernels built in memory: `furnsh` of a DAF handle of `dafops`, `spkops`, `ckops` or
# `pckops`, and the shortcuts `load_samples!` and `load_function!`.
#
# A handle in memory holds the bytes that the writers would put in a file. A load copies
# them into a new vector of words and parses it as `open_daf` parses a file, so the
# loaded kernel is an ordinary `DAF` in the load list, and writing to the handle later
# never changes it (the snapshot rule of `KernelSet`).

"""
    _daf_of_handle(h, name, who) -> DAF

Return the `DAF`, named `name`, of the bytes that the open handle `h` holds now, with the
checks of `open_daf`. An array that is being added is not in the result. Throws
`SpiceError(:dafnosuchhandle)` if `h` is closed.
"""
function _daf_of_handle(h::DAFHandle, name::AbstractString, who::AbstractString)
    _daf_check_open(h, who)
    bytes = if _daf_inmemory(h)
        _daf_bytes(h)
    else
        seekstart(h.io)
        read(h.io)
    end
    native() = copy(reinterpret(Float64, view(bytes, 1:(8 * (length(bytes) ÷ 8)))))
    return _daf_parse(IOBuffer(bytes), bytes, name, native)
end

"""
    _daf_from_handle(h, name, who) -> DAF

Return `_daf_of_handle` for the load of an in-memory handle. Throws `SpiceError` with kind
`:dafnosuchhandle` if `h` is closed, `:dafinvalidaccess` if it has a file, and
`:dafnewconflict` if an array is being added.
"""
function _daf_from_handle(h::DAFHandle, name::AbstractString, who::AbstractString)
    _daf_check_open(h, who)
    _daf_inmemory(h) || throw(
        SpiceError(
            :dafinvalidaccess,
            "$who: the DAF $(h.path) is a file; close it and load the file by its path"
        )
    )
    h.adding && throw(
        SpiceError(
            :dafnewconflict,
            "$who: an array is being added to the DAF in memory; end it with dafena first"
        )
    )
    return _daf_of_handle(h, name, who)
end

"""
    furnsh([ctx,] handle::DAFHandle; name, bodies=()) -> Nothing

Load the segments that are in the in-memory DAF `handle` (from [`spkops`](@ref),
[`ckops`](@ref), [`pckops`](@ref) or [`dafops`](@ref)) into `ctx` as an SPK, CK or
binary PCK kernel, with no file. This method has no CSPICE counterpart.

The kernel enters the load list under `name` with the rules of a file, so
[`kdata`](@ref), [`unload`](@ref) and [`kclear`](@ref) take `name`. The call copies the
segments: writing to `handle` afterwards does not change the loaded kernel, and `handle`
stays open.

`bodies`, pairs such as `["MYSAT" => -1001]`, are body names that hold while the kernel
is loaded. They rank below the names of the kernel pool (`NAIF_BODY_NAME`) and above
those of [`boddef`](@ref) and the built-in names; among kernels, the later one wins.

Throws `SpiceError`:

- `:emptystring` or `:blankfilename` if `name` is empty or blank, `:badarg` if a name of
  `bodies` is blank or its code is not a 32-bit integer;
- `:dafnosuchhandle` if `handle` is closed;
- `:dafinvalidaccess` if `handle` has a file (load the file by its path);
- `:dafnewconflict` if an array is being added (`dafbna` without `dafena`);
- `:baddaf` if `handle` is not an SPK, CK or binary PCK (see [`dafops`](@ref));
- `:nosegmentsfound` if `handle` holds no segment.

A load that throws leaves `ctx` unchanged.

# Examples
```jldoctest
julia> h = spkops();

julia> spkw09(h, -1000, 399, "J2000", 0.0, 10.0, "line", 1,
              [[7000.0, 0, 0, 0, 1.0, 0], [7000.0, 10.0, 0, 0, 1.0, 0]], [0.0, 10.0])

julia> ctx = SpiceContext();

julia> furnsh(ctx, h; name = "line");

julia> spkpos(ctx, -1000, 5.0, "J2000", "NONE", 399)[1] == [7000.0, 5.0, 0.0]
true

julia> spkcls(h)
```

See also [`save_kernel`](@ref), [`load_samples!`](@ref), [`load_function!`](@ref).
"""
function furnsh(ctx::SpiceContext, h::DAFHandle; name::AbstractString, bodies = ())
    _chk_filename("furnsh", name)
    kb = _kernel_bodies("furnsh", bodies)
    d = _daf_from_handle(h, name, "furnsh")
    daf = DAF(d.path, d.kind, d.data, d.segments, d.cksegments, d.names, kb)
    (isempty(daf.segments) && isempty(daf.cksegments)) && throw(
        SpiceError(
            :nosegmentsfound,
            "furnsh: the DAF in memory $(repr(name)) holds no segment"
        )
    )
    lock(ctx.lock) do
        list = copy(kernelset(ctx).loaded)
        push!(list, Loaded(String(name), "", daf))
        _pool_x_swap_loaded!(ctx, list, Loaded[])
    end
    return nothing
end
furnsh(h::DAFHandle; kw...) = furnsh(default_context(), h; kw...)

# ---------------------------------------------------------------------------------------
# Coverage of a handle, and saving
# ---------------------------------------------------------------------------------------

"Return the DAF of the handle `h` for the coverage call `who`, which needs the kind `kind`."
function _cov_daf(h::DAFHandle, who::AbstractString, kind::Symbol)
    name = _daf_inmemory(h) ? "(DAF in memory)" : h.path
    daf = _daf_of_handle(h, name, who)
    daf.kind === kind || throw(
        SpiceError(
            :invalidfiletype,
            "$who: the DAF $name is a binary $(daf.kind) kernel; the file must be a binary $kind file"
        )
    )
    return daf
end

"`spkobj!` and `pckfrm!` of a handle: the targets of the segments, for a DAF of `kind`."
_cov_ids!(ids::AbstractVector{<:Integer}, h::DAFHandle, who::AbstractString, kind::Symbol) =
    _add_ids!(ids, Int[s.target for s in _cov_daf(h, who, kind).segments])

"`spkcov!` and `pckcov!` of a handle: the spans of the segments of `idcode`, in file order."
function _cov_window!(cover::Window, h::DAFHandle, idcode::Integer, who::AbstractString, kind::Symbol)
    for s in _cov_daf(h, who, kind).segments
        s.target == idcode && wninsd!(cover, s.t0, s.t1)
    end
    return cover
end

spkobj!(ids::AbstractVector{<:Integer}, h::DAFHandle) = _cov_ids!(ids, h, "spkobj", :SPK)
spkobj(h::DAFHandle) = spkobj!(Int[], h)
spkcov!(cover::Window, h::DAFHandle, idcode::Integer) = _cov_window!(cover, h, idcode, "spkcov", :SPK)
spkcov(h::DAFHandle, idcode::Integer) = spkcov!(Window(), h, idcode)
pckfrm!(ids::AbstractVector{<:Integer}, h::DAFHandle) = _cov_ids!(ids, h, "pckfrm", :PCK)
pckfrm(h::DAFHandle) = pckfrm!(Int[], h)
pckcov!(cover::Window, h::DAFHandle, idcode::Integer) = _cov_window!(cover, h, idcode, "pckcov", :PCK)
pckcov(h::DAFHandle, idcode::Integer) = pckcov!(Window(), h, idcode)

# `ckobj`, `ckobj!`, `ckcov` and `ckcov!` of `ck_api.jl` reach a handle through this method.
_ck_open(f, h::DAFHandle, who::AbstractString) = f(_cov_daf(h, who, :CK))

"""
    save_kernel(handle, path) -> Nothing

Write the in-memory DAF `handle` (from [`spkops`](@ref), [`ckops`](@ref),
[`pckops`](@ref) or [`dafops`](@ref)) to the new file `path`: the bytes that the same
write calls put in a file made by [`spkopn`](@ref). An array that is being added is not
in the file. An existing file is never replaced, and the handle stays open. This call has
no CSPICE counterpart.

Throws `SpiceError`:

- `:emptystring` or `:blankfilename` if `path` is empty or blank;
- `:dafnosuchhandle` if `handle` is closed, `:dafinvalidaccess` if it has a file;
- `:fileopenfail` if the file exists or cannot be created, `:dafwritefail` if the write
  fails (the file is then removed).

See also [`furnsh`](@ref).
"""
function save_kernel(h::DAFHandle, path::AbstractString)
    _chk_filename("save_kernel", path)
    _daf_check_open(h, "save_kernel")
    _daf_inmemory(h) || throw(
        SpiceError(
            :dafinvalidaccess,
            "save_kernel: the DAF $(h.path) is a file, not a DAF in memory"
        )
    )
    bytes = _daf_bytes(h)
    fd = try
        Base.Filesystem.open(
            path, Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL |
                Base.Filesystem.JL_O_WRONLY, 0o644
        )
    catch e
        _io_failure(e) || rethrow()
        throw(
            SpiceError(
                :fileopenfail,
                "save_kernel: attempt to create the new file $path failed: $(sprint(showerror, e))"
            )
        )
    end
    try
        write(fd, bytes)
        close(fd)
    catch e
        close(fd)
        rm(path; force = true)           # no partial file stays behind
        throw(
            SpiceError(
                :dafwritefail,
                "save_kernel: the write to $path failed: $(sprint(showerror, e))"
            )
        )
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Shortcuts
# ---------------------------------------------------------------------------------------

"""
    _load_segment!(write!, ctx, name, segid, bodies) -> Nothing

Write one segment with `write!(h, segid)` into a new in-memory SPK and load it into `ctx`
under `name` with the names `bodies`. A call that throws leaves `ctx` unchanged.
"""
function _load_segment!(write!, ctx::SpiceContext, name::AbstractString, segid::AbstractString, bodies)
    h = spkops(name)
    try
        write!(h, first(segid, 40))
        furnsh(ctx, h; name = name, bodies = bodies)
    finally
        dafcls(h)
    end
    return nothing
end

"""
    load_samples!([ctx,] name, body, center, frame, ets, states;
                  method=:hermite, degree=7, segid=name, bodies=()) -> Nothing

Load the states `states` (km, km/s) of `body` relative to `center` in `frame` at the epochs
`ets` (seconds past J2000 TDB, increasing) into `ctx` as a kernel named `name`, with no
file; queries interpolate between them. This call has no CSPICE counterpart.

The kernel holds one SPK segment from `ets[1]` to `ets[end]`: type 13 (Hermite
interpolation, which uses the velocities) for `method = :hermite`, type 9 (Lagrange
interpolation) for `:lagrange`, with polynomials of degree `degree`. Type 13 needs an odd
`degree` in `1:27` and at least `(degree + 1) ÷ 2` states; type 9 a `degree` in `1:27` and
more than `degree` states. `states` is a 6 by n matrix, a vector of 6-vectors, or a flat
vector, as for [`spkw13`](@ref). The kernel is the one that [`spkw13`](@ref) or
[`spkw09`](@ref) writes to a file, so queries give the same results as on that file.
`segid` is cut to its first 40 characters. `bodies` are body names of the kernel, as for
[`furnsh`](@ref) of a handle; `body` and `center` are NAIF codes or names, and may be
names of `bodies`. A trajectory with jumps (burns) needs one call per arc.

Throws `SpiceError`:

- `:emptystring` or `:blankfilename` if `name` is empty or blank;
- `:invalidargument` if `method` is not `:hermite` or `:lagrange`;
- `:toofewstates` if `ets` is empty;
- `:idcodenotfound` for an unknown body name;
- the errors of [`spkw13`](@ref) or [`spkw09`](@ref) (degree, order of the epochs, ...).

A call that throws leaves `ctx` unchanged.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> ets = [0.0, 60.0, 120.0];

julia> states = [[7000.0 + 0.5t, 0, 0, 0.5, 0, 0] for t in ets];

julia> load_samples!(ctx, "probe", "PROBE", 399, "J2000", ets, states; degree = 3,
                     bodies = ["PROBE" => -1000])

julia> spkpos(ctx, "PROBE", 30.0, "J2000", "NONE", "EARTH")[1] ≈ [7015.0, 0, 0]
true
```

See also [`load_function!`](@ref), [`furnsh`](@ref), [`spkops`](@ref).
"""
function load_samples!(
        ctx::SpiceContext, name::AbstractString, body::NameOrID, center::NameOrID,
        frame::NameOrID, ets::AbstractVector{<:Real}, states; method::Symbol = :hermite,
        degree::Integer = 7, segid::AbstractString = name, bodies = ()
    )
    who = "load_samples!"
    _chk_filename(who, name)
    method in (:hermite, :lagrange) || throw(
        SpiceError(
            :invalidargument,
            "$who: the method is $(repr(method)); it must be :hermite or :lagrange"
        )
    )
    isempty(ets) && throw(SpiceError(:toofewstates, "$who: there are no epochs"))
    kb = _kernel_bodies(who, bodies)
    ks = kernelset(ctx)
    b, c = _resolve_body(ks, body, kb), _resolve_body(ks, center, kb)
    writer = method === :hermite ? spkw13 : spkw09
    _load_segment!(ctx, name, segid, kb) do h, id
        writer(ctx, h, b, c, frame, ets[1], ets[end], id, degree, states, ets)
    end
    return nothing
end
load_samples!(name::AbstractString, body::NameOrID, center::NameOrID, frame::NameOrID, ets::AbstractVector{<:Real}, states; kw...) =
    load_samples!(default_context(), name, body, center, frame, ets, states; kw...)

"""
    load_function!([ctx,] name, f, body, center, frame, t0, t1;
                   tol=nothing, ncoef=16, maxintervals=65536, segid=name, bodies=())
        -> (poserr, velerr)

Fit the trajectory `f(et)` of `body` relative to `center` in `frame` on `t0 <= et <= t1`
(seconds past J2000 TDB) with Chebyshev polynomials, and load the fit into `ctx` as a kernel
named `name`, with no file. `f` is called with `Float64` epochs, on the calling task, during
the call only, so a query costs what a query of an SPK file costs. This call has no CSPICE
counterpart.

`f` returns six numbers, the position (km) and the velocity (km/s), or three, the position
only. With six the kernel has an SPK segment of type 3 (a series for the position and one
for the velocity); with three, of type 2 (the velocity is the derivative of the position).
The fit takes `ncoef` coefficients on each of `n` intervals of equal length and doubles `n`
until the position error is at most `tol` (km). Without `tol` the goal is `1e-12` times the
largest distance, with the early stop of [`ephemeris_table`](@ref). The result is the
measured largest position error (km) and velocity error (km/s, `NaN` for three numbers).
`body`, `center`,
`segid` and `bodies` are as for [`load_samples!`](@ref).

An `ODESolution` goes in through a closure; see
[Kernels in memory](@ref).

Throws `SpiceError`:

- `:emptystring` or `:blankfilename` if `name` is empty or blank;
- `:badarg` unless `t0 < t1` are finite, `2 <= ncoef <= 28`, `tol > 0` and
  `maxintervals >= 1`;
- `:badarraysize` if `f` returns something other than three or six numbers;
- `:noconvergence` if `maxintervals` intervals do not reach the tolerance;
- `:idcodenotfound` for an unknown body name, and the errors of `f`.

A call that throws leaves `ctx` unchanged.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> gm, pv0 = 398600.4418, [7000.0, 0, 0, 0, 7.5, 1.0];

julia> poserr, velerr = load_function!(ctx, "orbit", et -> prop2b(gm, pv0, et), -1000, 399,
                                        "J2000", 0.0, 86400.0; tol = 1e-6);

julia> poserr <= 1e-6
true

julia> vnorm(spkpos(ctx, -1000, 5000.0, "J2000", "NONE", 399)[1] - prop2b(gm, pv0, 5000.0)[1:3]) <= 1e-6
true
```

See also [`load_samples!`](@ref), [`ephemeris_table`](@ref), [`add_ephemeris!`](@ref).
"""
function load_function!(
        ctx::SpiceContext, name::AbstractString, f, body::NameOrID, center::NameOrID,
        frame::NameOrID, t0::Real, t1::Real; tol::Union{Nothing, Real} = nothing,
        ncoef::Integer = 16, maxintervals::Integer = _TAB_MAXINTERVALS, segid::AbstractString = name,
        bodies = ()
    )
    who = "load_function!"
    _chk_filename(who, name)
    a, b = _tab_checks(who, t0, t1, ncoef, tol, maxintervals, 28)
    kb = _kernel_bodies(who, bodies)
    ks = kernelset(ctx)
    bc, cc = _resolve_body(ks, body, kb), _resolve_body(ks, center, kb)
    n = length(f(a))
    function checked(et, ::Val{N}) where {N}
        r = f(et)
        length(r) == N || throw(
            SpiceError(
                :badarraysize,
                "$who: the function returned $(length(r)) numbers at ET $et s; it must return $n"
            )
        )
        return SVector{N, Float64}(ntuple(k -> Float64(r[k]), Val(N)))
    end
    if n == 6
        coef, dt, nint, poserr, velerr = _tab_refine(
            et -> checked(et, Val(6)), a, b, tol, Int(ncoef), maxintervals, Val(2), who,
            " km"; split = true
        )
    elseif n == 3
        coef, dt, nint, poserr, velerr = _tab_refine(
            et -> checked(et, Val(3)), a, b, tol, Int(ncoef), maxintervals, Val(1), who, " km"
        )
    else
        throw(
            SpiceError(
                :badarraysize,
                "$who: the function returned $n numbers at ET $a s; it must return 3 or 6"
            )
        )
    end
    cdata = reshape(coef, n * Int(ncoef), nint)
    writer = n == 6 ? spkw03 : spkw02
    _load_segment!(ctx, name, segid, kb) do h, id
        writer(ctx, h, bc, cc, frame, a, b, id, dt, cdata, a)
    end
    return poserr, velerr
end
load_function!(name::AbstractString, f, body::NameOrID, center::NameOrID, frame::NameOrID, t0::Real, t1::Real; kw...) =
    load_function!(default_context(), name, f, body, center, frame, t0, t1; kw...)
