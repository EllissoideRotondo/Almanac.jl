# Kernel-pool, name, frame-name and time interface on a context.

"Return the value of pool variable `name`; throws `SpiceError(:kernelvarnotfound)` if it is absent."
function pool_value(ks::KernelSet, name::AbstractString)
    v = get(ks.pool, strip(name), nothing)          # a `SubString` key: no copy of the name
    v === nothing && throw(SpiceError(:kernelvarnotfound, "kernel-pool variable $name is not loaded"))
    return v
end

"""
    gdpool([ctx,] name) -> Vector{Float64}

Return the numeric values of the kernel-pool variable `name`, as
[`gdpool_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gdpool_c.html).

The result is a new vector with every value of the variable. Different from CSPICE,
the call has no `start` and `room` arguments and no `found` flag: an absent variable
is an error. Use [`expool`](@ref) or [`dtpool`](@ref) to test first. Blanks around
`name` are ignored; letter case is significant.

Throws `SpiceError(:kernelvarnotfound)` if the variable is not in the pool and
`SpiceError(:typemismatch)` if it holds character values.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pdpool(ctx, "BODY399_RADII", [6378.1366, 6378.1366, 6356.7519]);

julia> gdpool(ctx, "BODY399_RADII")
3-element Vector{Float64}:
 6378.1366
 6378.1366
 6356.7519
```

See also [`gipool`](@ref), [`gcpool`](@ref), [`bodvrd`](@ref), [`pdpool`](@ref).
"""
function gdpool(ks::KernelSet, name::AbstractString)
    v = pool_value(ks, name)
    v isa Vector{Float64} || throw(SpiceError(:typemismatch, "kernel-pool variable $name is not numeric"))
    return copy(v)
end

"""
    gipool([ctx,] name) -> Vector{Int}

Return the values of the numeric kernel-pool variable `name` rounded to the nearest
integers, a tie away from zero, as
[`gipool_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gipool_c.html).

The result holds every value; there is no `found` flag (see [`gdpool`](@ref)).

Throws `SpiceError(:kernelvarnotfound)` if the variable is not in the pool,
`SpiceError(:typemismatch)` if it holds character values, and
`SpiceError(:intoutofrange)` if a value is outside the 32-bit integers or NaN.
"""
function gipool(ks::KernelSet, name::AbstractString)
    v = gdpool(ks, name)
    for (i, x) in enumerate(v)
        # `gipool` of CSPICE: DBLE(INTMIN()) <= x <= DBLE(INTMAX()), then NINT
        _gipool_value(x) === nothing && _gipool_range_error(name, i, x)
    end
    return Int[_nint(x) for x in v]
end

"""
    gcpool([ctx,] name) -> Vector{String}

Return the character values of the kernel-pool variable `name`, as
[`gcpool_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gcpool_c.html).

The result is a new vector with every value; there is no `found` flag (see
[`gdpool`](@ref)).

Throws `SpiceError(:kernelvarnotfound)` if the variable is not in the pool and
`SpiceError(:typemismatch)` if it holds numeric values.
"""
function gcpool(ks::KernelSet, name::AbstractString)
    v = pool_value(ks, name)
    v isa Vector{String} || throw(SpiceError(:typemismatch, "kernel-pool variable $name is not character"))
    return copy(v)
end

"""
    _pool_key(prefix, mid, s1, s2="", upper=false) -> String

The pool name `prefix * mid * s1 * s2`, built in one buffer: `mid` is an integer (its
decimal digits) or a string (in upper case if `upper`).
"""
function _pool_key(
        prefix::String, mid::Union{Integer, AbstractString}, s1::AbstractString,
        s2::AbstractString = "", upper::Bool = false
    )
    if mid isa Integer
        n = abs(Int(mid))
        nmid = ndigits(n) + (mid < 0)
    else
        nmid = ncodeunits(mid)
    end
    np, n1, n2 = ncodeunits(prefix), ncodeunits(s1), ncodeunits(s2)
    buf = Base.StringVector(np + nmid + n1 + n2)
    copyto!(buf, 1, codeunits(prefix), 1, np)
    if mid isa Integer
        mid < 0 && (buf[np + 1] = UInt8('-'))
        # Loop invariant: the `k` last digits of `n` are written, the last at `np + nmid`.
        for k in 0:(ndigits(n) - 1)
            n, d = divrem(n, 10)
            buf[np + nmid - k] = UInt8('0') + d
        end
    else
        for k in 1:nmid
            c = codeunit(mid, k)
            buf[np + k] = upper ? _upper(c) : c
        end
    end
    copyto!(buf, np + nmid + 1, codeunits(s1), 1, n1)
    copyto!(buf, np + nmid + n1 + 1, codeunits(s2), 1, n2)
    return String(buf)
end

"The pool name `BODY<code>_<item>`, with the blanks around `item` removed."
_body_key(code::Integer, item::AbstractString) = _pool_key("BODY", code, "_", strip(item))

# `bodvrd_c` signals SPICE(NOTRANSLATION) for a name without a code, where the state
# calls signal SPICE(IDCODENOTFOUND).
function _bodvrd_code(ks::KernelSet, body)
    try
        return body_code(ks, body)
    catch e
        (e isa SpiceError && e.kind === :idcodenotfound) || rethrow()
        throw(SpiceError(:notranslation, "bodvrd: the body name '$body' has no NAIF ID code"))
    end
end

"""
    bodvrd([ctx,] body, item) -> Vector{Float64}

Return the values of the kernel-pool variable `BODY<id>_<item>` for `body`, as
[`bodvrd_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html).

`body` is a body name, a string that holds an integer (`"399"`) or an integer NAIF ID
code. `item` is case sensitive, as in CSPICE;
blanks around it are ignored. The values have the units of the kernel: km for
`"RADII"`, km^3/s^2 for `"GM"` in the NAIF kernels. There is no `maxn` argument.

Throws `SpiceError(:notranslation)` if the body name is not known,
`SpiceError(:kernelvarnotfound)` if the variable is not in the pool and `SpiceError(:typemismatch)`
if it holds character values.

# Examples
```julia
furnsh("pck00011.tpc")
a, b, c = bodvrd("EARTH", "RADII")    # km
```

See also [`bodvcd`](@ref), [`gdpool`](@ref).
"""
bodvrd(ks::KernelSet, body::NameOrID, item::AbstractString) =
    gdpool(ks, _body_key(_bodvrd_code(ks, body), item))

"""
    bods2c([ctx,] name) -> Union{Int,Nothing}

Return the NAIF ID code of the body `name`, or `nothing` if the name is not known, as
[`bods2c_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bods2c_c.html).

`name` is a body name or a string that holds an integer, such as `"399"`. Letter case
and blanks around the name are ignored. Names from the kernel pool and from
[`boddef`](@ref) have priority over the built-in names. `nothing` takes the place of
the `found` flag of CSPICE.

# Examples
```jldoctest
julia> bods2c("Earth"), bods2c("-77"), bods2c("NO SUCH BODY")
(399, -77, nothing)
```

See also [`bodn2c`](@ref), [`bodc2n`](@ref).
"""
function bods2c(ks::KernelSet, name::AbstractString)
    try
        return Int(body_code(ks, name))
    catch e
        e isa SpiceError && e.kind === :idcodenotfound && return nothing
        rethrow()
    end
end

"""
    bodn2c([ctx,] name) -> Union{Int,Nothing}

Return the NAIF ID code of the body name `name`, or `nothing` if the name is not
known, as
[`bodn2c_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodn2c_c.html).

Different from [`bods2c`](@ref), a string that holds an integer is not a name.
`nothing` takes the place of the `found` flag of CSPICE.

# Examples
```jldoctest
julia> bodn2c("MOON"), bodn2c("301")
(301, nothing)
```
"""
function bodn2c(ks::KernelSet, name::AbstractString)
    code = bods2c(ks, name)
    code === nothing && return nothing
    # `bods2c` also accepts integer strings (`_parse_int_name`, as `beint`); accept one
    # only if it is a defined name.
    if _parse_int_name(name) !== nothing
        haskey(ks.name2code, normalize_name(name)) || return nothing
    end
    return code
end

"""
    bodc2n([ctx,] code) -> Union{String,Nothing}

Return the name of the body with the NAIF ID `code`, or `nothing` if the code has no
name, as
[`bodc2n_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodc2n_c.html).

A name from the kernel pool or from [`boddef`](@ref) has priority over the built-in
name.

# Examples
```jldoctest
julia> bodc2n(399), bodc2n(123456789)
("EARTH", nothing)
```

See also [`bodn2c`](@ref).
"""
bodc2n(ks::KernelSet, code::Integer) = body_name(ks, code)

"""
    bodc2s([ctx,] code) -> String

Return the name of the body with the NAIF ID `code`, or the decimal text of `code` if
the code has no name, as
[`bodc2s_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodc2s_c.html).

The name is the one of [`bodc2n`](@ref).

# Examples
```jldoctest
julia> bodc2s(399), bodc2s(123456789)
("EARTH", "123456789")
```

See also [`bods2c`](@ref).
"""
function bodc2s(ks::KernelSet, code::Integer)
    name = body_name(ks, code)
    return name === nothing ? string(code) : name
end

"""
    scn2id([ctx,] clknam) -> Union{Int,Nothing}

Return the NAIF ID code of the spacecraft clock `clknam`, or `nothing` if the name is
not known, as
[`scn2id`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/scn2id.html)
of SPICELIB. CSPICE has no C wrapper for it.

The name is usually the body name of the spacecraft with the word `SCLK` after it, as
in `"MGS SCLK"`. Only the first 32 characters count, as in SPICELIB; letter case is
ignored, the last `SCLK` among them is blanked (it need not be a separate word), and the
rest is looked up with [`bodn2c`](@ref), so a bare body name gives the same code.

# Examples
```jldoctest
julia> scn2id("CASSINI SCLK"), scn2id("cassini"), scn2id("NO SUCH CLOCK")
(-82, -82, nothing)
```

See also [`scid2n`](@ref).
"""
function scn2id(ks::KernelSet, clknam::AbstractString)
    # `ucase` into a 32-character buffer, then the last "SCLK" becomes blanks.
    s = codeunits(String(clknam))
    b = UInt8[_upper(c) for c in view(s, 1:min(length(s), 32))]
    for i in (length(b) - 3):-1:1
        if b[i] == UInt8('S') && b[i + 1] == UInt8('C') && b[i + 2] == UInt8('L') &&
                b[i + 3] == UInt8('K')
            b[i:(i + 3)] .= UInt8(' ')
            break
        end
    end
    name = String(b)
    isempty(strip(name)) && return nothing
    return bodn2c(ks, name)
end

"""
    scid2n([ctx,] clkid) -> Union{String,Nothing}

Return the name of the spacecraft clock with the NAIF ID `clkid`, or `nothing` if the
code has no body name, as
[`scid2n`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/scid2n.html)
of SPICELIB. CSPICE has no C wrapper for it.

The name is the body name of [`bodc2n`](@ref) with `" SCLK"` after it.

# Examples
```jldoctest
julia> scid2n(-82), scid2n(123456789)
("CASSINI SCLK", nothing)
```

See also [`scn2id`](@ref).
"""
function scid2n(ks::KernelSet, clkid::Integer)
    name = body_name(ks, clkid)
    return name === nothing ? nothing : string(rstrip(name), " SCLK")
end

"""
    namfrm([ctx,] name) -> Int

Return the ID of the frame `name`, or 0 if the frame is not known, as
[`namfrm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/namfrm_c.html).

Built-in frames and frames defined in the kernel pool are known.

# Examples
```jldoctest
julia> namfrm("J2000"), namfrm("eclipj2000"), namfrm("NO_SUCH_FRAME")
(1, 17, 0)
```

See also [`frmnam`](@ref), [`frinfo`](@ref).
"""
function namfrm(ks::KernelSet, name::AbstractString)
    try
        return Int(frame_id(ks, name))
    catch e
        e isa SpiceError && e.kind === :unknownframe && return 0
        rethrow()
    end
end

"""
    frmnam([ctx,] id) -> String

Return the name of the frame with ID `id`, or `""` if the frame is not known, as
[`frmnam_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/frmnam_c.html).

# Examples
```jldoctest
julia> frmnam(17), frmnam(-123456)
("ECLIPJ2000", "")
```

See also [`namfrm`](@ref).
"""
function frmnam(ks::KernelSet, id::Integer)
    name = frame_name(ks, id)
    return name === nothing ? "" : name
end

# Context and default-context forms of the functions above.
for f in (:gdpool, :gipool, :gcpool, :bods2c, :bodn2c, :scn2id, :namfrm)
    @eval begin
        $f(ctx::SpiceContext, name::AbstractString) = $f(kernelset(ctx), name)
        $f(name::AbstractString) = $f(default_context(), name)
    end
end
for f in (:bodc2n, :bodc2s, :scid2n, :frmnam)
    @eval begin
        $f(ctx::SpiceContext, code::Integer) = $f(kernelset(ctx), code)
        $f(code::Integer) = $f(default_context(), code)
    end
end
bodvrd(ctx::SpiceContext, body::NameOrID, item::AbstractString) =
    bodvrd(kernelset(ctx), body, item)
bodvrd(body::NameOrID, item::AbstractString) = bodvrd(default_context(), body, item)
