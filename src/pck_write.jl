# Binary PCK writer: `pckopn`, `pckw02` and `pckcls`. A binary PCK is a DAF of the type
# PCK whose summaries hold 2 doubles (the coverage) and 5 integers (frame class ID of the
# body-fixed frame, base frame, data type, first and last address). `pckw02` makes the
# input checks of SPICELIB `pckw02` in the same order and stores the records of
# `spkw02` (`_seg_cheb_records`), so the file is equal byte for byte to the one that
# CSPICE writes from the same input.

"""
    _pck_irfnum(frame) -> Int

ID of the built-in inertial frame `frame` (a name or an ID), or 0. A name is compared
as `irfnum_` compares it: letter case and blanks do not count.
"""
function _pck_irfnum(frame::SegFrame)
    frame isa Integer && return 1 <= frame <= N_INERTIAL_FRAMES ? Int(frame) : 0
    key = uppercase(replace(frame, ' ' => ""))
    for (i, def) in enumerate(INERTIAL_FRAME_DEFS)   # invariant: no frame before i has the name
        def[1] == key && return i
    end
    return 0
end

"""
    pckopn(path, ifname, ncomch=0) -> DAFHandle

Create the binary PCK file `path` and open it for writing, as
[`pckopn_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pckopn_c.html).

`ifname` is the internal file name and `ncomch` is the number of comment characters to
make room for. The file is on disk when the call returns, byte for byte the file that
`pckopn_c` creates. An existing file is never replaced. Add segments with
[`pckw02`](@ref) and close the handle with [`pckcls`](@ref). Load the file with
[`furnsh`](@ref) only after the close.

Throws `SpiceError(:emptystring)` if `path` or `ifname` is empty and `:blankfilename` if
`path` is blank, as CSPICE, and `:fileopenfailed` if the file exists or cannot be created.

# Examples
Write a segment for the frame IAU_MOON (frame class ID 301) whose angles are constant,
then load the file:
```jldoctest
julia> path = joinpath(mktempdir(), "example.bpc");

julia> h = pckopn(path, "example PCK");

julia> cdata = [[0.5, 0.25, 1.0]];    # phi, delta, w (rad), degree 0

julia> pckw02(h, 301, "J2000", 0.0, 100.0, "moon", 100.0, cdata, 0.0)

julia> pckcls(h)

julia> ctx = SpiceContext();

julia> furnsh(ctx, path)

julia> pxform(ctx, "J2000", "IAU_MOON", 50.0) ≈ eul2m(1.0, 0.25, 0.5, 3, 1, 3)
true
```

See also [`dafac`](@ref) to add comments.
"""
function pckopn(path::AbstractString, ifname::AbstractString, ncomch::Integer = 0)
    _seg_chkfstr("pckopn", "path", path)
    _seg_chkfstr("pckopn", "ifname", ifname)
    ncomr = ncomch > 0 ? (Int(ncomch) - 1) ÷ _DAF_CRECLEN + 1 : 0
    return dafonw(path, "PCK", 2, 5, ifname, ncomr)
end

"""
    pckops(ifname="", ncomch=0) -> DAFHandle

Open a new PCK file that lives in memory, for writing: [`pckopn`](@ref) without the
path. No file is made. The writers (such as [`pckw02`](@ref)) take the handle as they take
one of a file and write the same bytes. [`furnsh`](@ref) loads its segments as a kernel,
and [`save_kernel`](@ref) writes them to a file; [`pckcls`](@ref) discards them. This call
has no CSPICE counterpart; see [`dafops`](@ref).

Throws the errors of [`dafops`](@ref).
"""
function pckops(ifname::AbstractString = "", ncomch::Integer = 0)
    ncomr = ncomch > 0 ? (Int(ncomch) - 1) ÷ _DAF_CRECLEN + 1 : 0
    return dafops("PCK", 2, 5, ifname, ncomr)
end

"""
    pckw02(handle, clssid, frame, first, last, segid, intlen, cdata, btime) -> Nothing

Write a binary PCK segment of type 2 (Chebyshev polynomials for three Euler angles, on
intervals of equal length) to the file of `handle`, as
[`pckw02_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pckw02_c.html).

- `handle`: a [`DAFHandle`](@ref) from [`pckopn`](@ref).
- `clssid`: the frame class ID of the body-fixed frame, the `classid` of
  [`frinfo`](@ref), for example 301 for IAU_MOON.
- `frame`: the base frame, a built-in inertial frame given by name or by frame ID (1 to
  21, which [`bltfrm`](@ref) lists). Letter case and blanks in the name do not count.
- `first`, `last`: the coverage of the segment, seconds past J2000 TDB.
- `segid`: the segment identifier, at most 40 printable characters.
- `intlen`, `btime`: the length of every interval (s) and the start of the first interval
  (seconds past J2000 TDB).
- `cdata`: `n` records, one for each interval. All records have the same length, a
  multiple of 3: the coefficients (rad) of the angle `phi`, then of `delta`, then of
  `w`, a third of the record for each, so the polynomial degree is a third of the
  length minus 1. The matrix from the base frame to the body-fixed frame is
  `R3(w) * R1(delta) * R3(phi)`, which is `eul2m(w, delta, phi, 3, 1, 3)`.

`cdata` is a vector of vectors or a matrix with one record in each column, as in
[`spkw02`](@ref). The record count `n` and the degree, which CSPICE takes as arguments,
follow from its shape.

The segment is on disk when the call returns. A call that throws writes nothing.

Throws `SpiceError`:
- `:dafnosuchhandle` if the handle is closed, `:dafinvalidaccess` if it is open for
  reading;
- `:emptystring` if `frame` or `segid` is empty;
- `:numcoeffsnotpos` if `cdata` has no record (an empty vector of vectors or a matrix
  without columns);
- `:badarg` if the records differ in length, or if their length is not a positive
  multiple of 3;
- `:intlennotpos` unless `intlen > 0`;
- `:invalidrefframe` if `frame` is not a built-in inertial frame;
- `:baddescrtimes` if `first > last`, or unless [`first`, `last`] is inside [`btime`,
  `btime + n * intlen`] (no tolerance);
- `:idstringtoolong` or `:nonprintablechars` if `segid` is not valid.

See [`pckopn`](@ref) for an example.
"""
function pckw02(
        h::DAFHandle, clssid::Integer, frame::SegFrame, first::Real, last::Real,
        segid::AbstractString, intlen::Real, cdata, btime::Real
    )
    who = "pckw02"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    cd, n, ninrec = _seg_flat(who, "cdata", cdata, 0)
    first, last, intlen, btime = Float64(first), Float64(last), Float64(intlen), Float64(btime)
    n <= 0 && throw(
        SpiceError(
            :numcoeffsnotpos,
            "$who: the number of sets of Euler angle coefficients is not positive: $n"
        )
    )
    ninrec % 3 == 0 && ninrec > 0 || throw(
        SpiceError(
            :badarg, "$who: every record of " *
                "`cdata` must hold 3 * (degree + 1) coefficients, it holds $ninrec"
        )
    )
    intlen <= 0.0 && throw(SpiceError(:intlennotpos, "$who: the interval length is not positive: $intlen"))
    refcod = _pck_irfnum(frame)
    refcod == 0 && throw(
        SpiceError(
            :invalidrefframe,
            "$who: the reference frame $frame is not supported"
        )
    )
    first > last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is greater than the segment end time $last"
        )
    )
    first < btime && throw(
        SpiceError(
            :baddescrtimes, "$who: the segment descriptor " *
                "start time $first is less than the beginning time of the segment data $btime"
        )
    )
    ltime = btime + n * intlen
    last > ltime && throw(
        SpiceError(
            :baddescrtimes, "$who: the segment descriptor end " *
                "time $last is greater than the end time of the segment data $ltime"
        )
    )
    _seg_chckid(who, segid, "PCK")
    _daf_add_array!(
        h, Float64[first, last], Int[clssid, refcod, 2, 0, 0], segid,
        _seg_cheb_records(cd, n, ninrec, intlen, btime)
    )
    return nothing
end

"""
    pckw20([ctx,] handle, clssid, frame, first, last, segid, intlen, cdata, ascale,
           tscale, initjd, initfr) -> Nothing

Write a binary PCK segment of type 20 (Chebyshev polynomials for the rates of three
Euler angles, on intervals of equal length) to the file of `handle`, as
[`pckw20`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/pckw20.html) of SPICELIB.

Each record of `cdata` holds `polydg + 1` coefficients of the rate of `phi`, then the
value of `phi` at the midpoint of the interval, then the same for `delta` and for `w`
(the angles of [`pckw02`](@ref)). The angles are in units of `ascale` radians and the
rates in units of `ascale` radians per `tscale` s. `intlen` is the length of every
interval in days, and `initjd + initfr` is the TDB Julian date of the start of the first
interval, as in [`spkw20`](@ref). Different from [`pckw02`](@ref), whose base frame must
be a built-in inertial frame, `frame` is any frame that `ctx` knows, by name or ID, as
`pckw20_` of CSPICE looks it up with `namfrm`. See [`pckw02`](@ref) for the other
arguments.

The segment is on disk when the call returns. A call that throws writes nothing.

Throws `SpiceError`:
- `:dafnosuchhandle` if the handle is closed, `:dafinvalidaccess` if it is open for
  reading;
- `:emptystring` if `frame` or `segid` is empty;
- `:invalidcount` if `cdata` has no record;
- `:badarg` if the records differ in length, or if their length is not 3 * (degree + 2);
- `:invaliddegree` if the degree is greater than 50;
- `:intlennotpos` unless `intlen > 0`;
- `:invalidrefframe` if `frame` is a name that `ctx` does not know;
- `:baddescrtimes` unless `first < last`;
- `:nonpositivescale` unless `ascale > 0` and `tscale > 0`;
- `:coveragegap` unless [`first`, `last`] is inside the intervals of the data, with the
  tolerance of [`spkw20`](@ref);
- `:idstringtoolong` or `:nonprintablechars` if `segid` is not valid.
"""
function pckw20(
        ctx::SegContext, h::DAFHandle, clssid::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, intlen::Real, cdata, ascale::Real, tscale::Real,
        initjd::Real, initfr::Real
    )
    who = "pckw20"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    cd, n, ninrec = _seg_flat(who, "cdata", cdata, 0)
    first, last, intlen = Float64(first), Float64(last), Float64(intlen)
    ascale, tscale, initjd, initfr = Float64(ascale), Float64(tscale), Float64(initjd), Float64(initfr)
    n <= 0 && throw(
        SpiceError(
            :invalidcount,
            "$who: the number of sets of coordinate coefficients is not positive: $n"
        )
    )
    (ninrec % 3 == 0 && ninrec >= 6) || throw(
        SpiceError(
            :badarg, "$who: every record of " *
                "`cdata` must hold 3 * (degree + 2) values, it holds $ninrec"
        )
    )
    polydg = ninrec ÷ 3 - 2
    polydg > 50 && throw(
        SpiceError(
            :invaliddegree, "$who: the interpolating polynomials " *
                "have degree $polydg; the valid degree range is [0, 50]"
        )
    )
    intlen <= 0.0 && throw(SpiceError(:intlennotpos, "$who: the interval length is not positive: $intlen"))
    refcod = _seg_refcode(who, ctx, frame)
    first >= last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is not less than the segment end time $last"
        )
    )
    ascale <= 0.0 && throw(SpiceError(:nonpositivescale, "$who: the angle scale is not positive: $ascale"))
    tscale <= 0.0 && throw(SpiceError(:nonpositivescale, "$who: the time scale is not positive: $tscale"))
    btime = 86400.0 * (initjd - 2451545.0 + initfr)
    ltime = btime + n * intlen * 86400.0
    tol = max(abs(btime), abs(ltime)) * 1.0e-13
    _seg_coverage_gap(who, first, last, btime, ltime, tol)
    _seg_chckid(who, segid, "PCK")
    data = append!(cd, (ascale, tscale, initjd, initfr, intlen, Float64(ninrec), Float64(n)))
    _daf_add_array!(h, Float64[first, last], Int[clssid, refcod, 20, 0, 0], segid, data)
    return nothing
end
pckw20(
    h::DAFHandle, clssid::Integer, frame::SegFrame, first::Real, last::Real,
    segid::AbstractString, intlen::Real, cdata, ascale::Real, tscale::Real, initjd::Real,
    initfr::Real
) = pckw20(
    default_context(), h, clssid, frame, first, last, segid, intlen,
    cdata, ascale, tscale, initjd, initfr
)

"""
    pckcls(handle) -> Nothing

Close the binary PCK file of `handle`, as
[`pckcls_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pckcls_c.html). The call does nothing if the handle
is closed already. The close writes no data: every write call puts its bytes on disk
before it returns. Load the file with [`furnsh`](@ref) after the close.

Throws `SpiceError(:nosegmentsfound)` if the file is open for writing and holds no
segment. The handle then stays open; close it with [`dafcls`](@ref).
"""
pckcls(h::DAFHandle) = _seg_close(h, "PCK")
