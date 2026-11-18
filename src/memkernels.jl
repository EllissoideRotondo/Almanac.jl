# Kernels built in memory: `furnsh` of a DAF handle of `dafops`, `spkops`, `ckops` or
# `pckops`, and the shortcuts `load_samples!` and `load_function!`.
#
# A handle in memory holds the bytes that the writers would put in a file. A load copies
# them into a new vector of words and parses it as `open_daf` parses a file, so the
# loaded kernel is an ordinary `DAF` in the load list, and writing to the handle later
# never changes it (the snapshot rule of `KernelSet`).

"""
    _daf_of_handle(h, name, who) -> DAF

Contract: the `DAF`, named `name`, of the bytes that the open handle `h` holds now (in
memory or in its file), with the checks of `open_daf`. An array that is being added has
no summary yet and is not in the result. Throws `SpiceError(:dafnosuchhandle)` if `h` is
closed.
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

Contract: `_daf_of_handle` for the load of an in-memory handle. Throws `SpiceError` with
kind `:dafnosuchhandle` if `h` is closed, `:dafinvalidaccess` if it has a file, and
`:dafnewconflict` if an array is being added (`dafbna` without `dafena`).
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

The kernel enters the load list under `name`, with the rules of a file: it has priority
over the kernels loaded before it, [`kdata`](@ref) lists `name` with its type,
[`unload`](@ref) with `name` and [`kclear`](@ref) remove it, and loading the same name
twice makes two entries. The call copies the segments: writing to `handle` afterwards
does not change the loaded kernel, and `handle` stays open, so it can take more segments
and be loaded again.

`bodies`, pairs such as `["MYSAT" => -1001]`, are body names that hold while the kernel
is loaded and go with it. A name that the kernel pool defines (`NAIF_BODY_NAME`) keeps the
definition of the pool; among kernels, the later one wins; a kernel name has priority over
one of [`boddef`](@ref) and over the built-in names. After [`unload`](@ref) the next
definition in this order holds again.

Throws `SpiceError`:

- `:emptystring` or `:blankfilename` if `name` is empty or blank, `:badarg` if a name of
  `bodies` is blank;
- `:dafnosuchhandle` if `handle` is closed;
- `:dafinvalidaccess` if `handle` has a file (load the file by its path);
- `:dafnewconflict` if an array is being added (`dafbna` without `dafena`);
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
    daf = DAF(d.path, d.kind, d.data, d.segments, d.cksegments, kb)
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

"The DAF of the handle `h` for the coverage call `who`, which needs the kind `kind`."
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

function spkobj!(ids::AbstractVector{<:Integer}, h::DAFHandle)
    daf = _cov_daf(h, "spkobj", :SPK)
    return _add_ids!(ids, Int[s.target for s in daf.segments])
end
spkobj(h::DAFHandle) = spkobj!(Int[], h)

function spkcov!(cover::Window, h::DAFHandle, idcode::Integer)
    for s in _cov_daf(h, "spkcov", :SPK).segments   # in file order, as the summaries
        s.target == idcode && wninsd!(cover, s.t0, s.t1)
    end
    return cover
end
spkcov(h::DAFHandle, idcode::Integer) = spkcov!(Window(), h, idcode)

function pckfrm!(ids::AbstractVector{<:Integer}, h::DAFHandle)
    daf = _cov_daf(h, "pckfrm", :PCK)
    return _add_ids!(ids, Int[s.target for s in daf.segments])
end
pckfrm(h::DAFHandle) = pckfrm!(Int[], h)

function pckcov!(cover::Window, h::DAFHandle, idcode::Integer)
    for s in _cov_daf(h, "pckcov", :PCK).segments
        s.target == idcode && wninsd!(cover, s.t0, s.t1)
    end
    return cover
end
pckcov(h::DAFHandle, idcode::Integer) = pckcov!(Window(), h, idcode)

# `ckobj`, `ckobj!`, `ckcov` and `ckcov!` of `ck_api.jl` reach a handle through this method.
_ck_open(f, h::DAFHandle, who::AbstractString) = f(_cov_daf(h, who, :CK))

"""
    save_kernel(handle, path) -> Nothing

Write the in-memory DAF `handle` (from [`spkops`](@ref), [`ckops`](@ref),
[`pckops`](@ref) or [`dafops`](@ref)) to the new file `path`: the bytes that the same
write calls would have put in a file made by [`spkopn`](@ref) (or `ckopn`, `pckopn`,
`dafonw`). The file is a valid kernel; an array that is being added (`dafbna` without
`dafena`) is not in it. An existing file is never replaced. The handle stays open. This
call has no CSPICE counterpart.

Throws `SpiceError`:

- `:emptystring` or `:blankfilename` if `path` is empty or blank;
- `:dafnosuchhandle` if `handle` is closed, `:dafinvalidaccess` if it has a file;
- `:fileopenfailed` if the file exists or cannot be created.

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
    try
        fd = Base.Filesystem.open(
            path, Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL |
                Base.Filesystem.JL_O_WRONLY, 0o644
        )
        try
            write(fd, bytes)
        finally
            close(fd)
        end
    catch e
        throw(
            SpiceError(
                :fileopenfailed,
                "save_kernel: attempt to create the new file $path failed: $(sprint(showerror, e))"
            )
        )
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Shortcuts
# ---------------------------------------------------------------------------------------

"""
    _load_segment!(write!, ctx, name, segid, bodies, who) -> Nothing

Contract: write one segment with `write!(h, segid)` into a new in-memory SPK and load it
into `ctx` under `name` with the names `bodies`. The handle is closed on return. A call
that throws leaves `ctx` unchanged.
"""
function _load_segment!(write!, ctx::SpiceContext, name::AbstractString, segid::AbstractString, bodies, who::AbstractString)
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
file. This call has no CSPICE counterpart.

The kernel holds one SPK segment from `ets[1]` to `ets[end]`: type 13 (Hermite
interpolation, which uses the velocities) for `method = :hermite`, type 9 (Lagrange
interpolation) for `:lagrange`, with polynomials of degree `degree`. `states` is a 6 by n
matrix, a vector of 6-vectors, or a flat vector, as for [`spkw13`](@ref). The kernel is the
one that [`spkw13`](@ref) or [`spkw09`](@ref) would write to a file, so queries give the
same results as on that file. `segid` is cut to its first 40 characters. `body`, `center`
and the names `bodies` are as for [`furnsh`](@ref) of a handle, and the rules of the load
list are those of a file. A trajectory with jumps (burns) needs one call per arc, or
several segments through [`spkops`](@ref).

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
    _load_segment!(ctx, name, segid, kb, who) do h, id
        writer(ctx, h, b, c, frame, ets[1], ets[end], id, degree, states, ets)
    end
    return nothing
end
load_samples!(name::AbstractString, body::NameOrID, center::NameOrID, frame::NameOrID, ets::AbstractVector{<:Real}, states; kw...) =
    load_samples!(default_context(), name, body, center, frame, ets, states; kw...)
