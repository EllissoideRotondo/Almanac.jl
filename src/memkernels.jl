# Kernels built in memory: `furnsh` of a DAF handle of `dafops`, `spkops`, `ckops` or
# `pckops`, and the shortcuts `load_samples!` and `load_function!`.
#
# A handle in memory holds the bytes that the writers would put in a file. A load copies
# them into a new vector of words and parses it as `open_daf` parses a file, so the
# loaded kernel is an ordinary `DAF` in the load list, and writing to the handle later
# never changes it (the snapshot rule of `KernelSet`).

"""
    _daf_from_handle(h, name, who) -> DAF

Contract: the `DAF` of the current bytes of the open in-memory handle `h`, named `name`,
with the checks of `open_daf`. Throws `SpiceError` with kind `:dafnosuchhandle` if `h` is
closed, `:dafinvalidaccess` if it has a file, and `:dafnewconflict` if an array is being
added (`dafbna` without `dafena`).
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
    bytes = _daf_bytes(h)
    native() = copy(reinterpret(Float64, view(bytes, 1:(8 * (length(bytes) ÷ 8)))))
    return _daf_parse(IOBuffer(bytes), bytes, name, native)
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
