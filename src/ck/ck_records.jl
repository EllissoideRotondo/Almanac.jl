# Records of CK segments of types 1 to 6 through a DAF handle: `cknrNN`, `ckgrNN`,
# `cknm06` and `ckmp06`, as the SPICELIB routines of the same names. The descriptor is
# the packed summary of `dafgs` (2 doubles, 6 integers); CSPICE does not check that it
# belongs to the file of the handle, and neither does this file.

"Integer components of the CK descriptor `descr`; throws `SpiceError(:ckwrongdatatype)` unless the data type is `dtype`."
function _ck_descr_ic(who::AbstractString, descr::AbstractVector{<:Real}, dtype::Int)
    _, ic = dafus(descr, 2, 6)
    ic[3] == dtype || throw(
        SpiceError(
            :ckwrongdatatype, "$who: the data type of the " *
                "segment should be $dtype; the descriptor shows type $(ic[3])"
        )
    )
    return ic
end

function _ck_check_recno(who::AbstractString, recno::Integer, nrec::Int)
    1 <= recno <= nrec || throw(
        SpiceError(
            :cknonexistrec, "$who: the requested record " *
                "number ($recno) does not exist; there are $nrec records in the segment"
        )
    )
    return nothing
end

"""
    cknr02(handle, descr) -> Int

Return the number of pointing records in the CK type 2 segment with the descriptor
`descr`, as
[`cknr02_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/cknr02_c.html).

`descr` is the packed descriptor of the segment (5 doubles), as [`dafgs`](@ref) returns
it during a search on `handle`. The count follows from the addresses in the descriptor,
so the file is not read, and a closed handle gives the count too.

Throws `SpiceError(:ckwrongdatatype)` if the segment is not of type 2, and
`:dafinvalidparams` if `descr` has fewer than 5 doubles.

# Examples
```julia
h = dafopr("mission.bc")
dafbfs(h)
while daffna(h)
    descr = dafgs(h)
    _, ic = dafus(descr, 2, 6)
    ic[3] == 2 || continue
    for i in 1:cknr02(h, descr)
        start, stop, rate = ckgr02(h, descr, i)[1:3]
    end
end
dafcls(h)
```

See also [`ckgr02`](@ref), [`cknr03`](@ref).
"""
function cknr02(h::DAFHandle, descr::AbstractVector{<:Real})
    ic = _ck_descr_ic("cknr02", descr, 2)
    return _nint(((ic[6] - ic[5] + 1) * 100.0 + 1.0) / 1001.0)
end

"""
    ckgr02(handle, descr, recno) -> Vector{Float64}

Return the pointing record `recno` of the CK type 2 segment with the descriptor `descr`
in the file of `handle`, as
[`ckgr02_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckgr02_c.html).
Records count from 1, as in CSPICE.

The record has 10 doubles: the start and the stop of its interval (ticks), the clock rate
(seconds per tick), the quaternion of the C-matrix at the start (SPICE form, scalar
first) and the constant angular velocity (rad/s, components in the reference frame of the
segment). `descr` is as in [`cknr02`](@ref).

Throws `SpiceError`:
- `:ckwrongdatatype` if the segment is not of type 2;
- `:cknonexistrec` unless `1 <= recno <= cknr02(handle, descr)`;
- `:dafnosuchhandle` if the handle is closed;
- `:dafinvalidparams` if `descr` has fewer than 5 doubles.
"""
function ckgr02(h::DAFHandle, descr::AbstractVector{<:Real}, recno::Integer)
    ic = _ck_descr_ic("ckgr02", descr, 2)
    nrec = cknr02(h, descr)
    _ck_check_recno("ckgr02", recno, nrec)
    beg = ic[5]
    addr = beg + 8 * (recno - 1)
    prec = dafgda(h, addr, addr + 7)             # quaternion, angular velocity, rate
    start = dafgda(h, beg + 8nrec + recno - 1, beg + 8nrec + recno - 1)[1]
    stop = dafgda(h, beg + 9nrec + recno - 1, beg + 9nrec + recno - 1)[1]
    record = Float64[
        start, stop, prec[8], prec[1], prec[2], prec[3], prec[4], prec[5],
        prec[6], prec[7],
    ]
    @assert length(record) == 10
    return record
end

"""
    cknr03(handle, descr) -> Int

Return the number of pointing instances in the CK type 3 segment with the descriptor
`descr` in the file of `handle`, as
[`cknr03_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/cknr03_c.html).
`descr` is as in [`cknr02`](@ref).

Throws `SpiceError(:ckwrongdatatype)` if the segment is not of type 3,
`:dafnosuchhandle` if the handle is closed, and `:dafinvalidparams` if `descr` has fewer
than 5 doubles.
"""
function cknr03(h::DAFHandle, descr::AbstractVector{<:Real})
    ic = _ck_descr_ic("cknr03", descr, 3)
    return _nint(dafgda(h, ic[6], ic[6])[1])
end

"""
    ckgr03(handle, descr, recno) -> Vector{Float64}

Return the pointing instance `recno` of the CK type 3 segment with the descriptor `descr`
in the file of `handle`, as
[`ckgr03_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckgr03_c.html).
Records count from 1, as in CSPICE.

The record is the time tag (ticks) and the quaternion of the C-matrix (SPICE form, scalar
first): 5 doubles for a segment without angular velocity, and 8 doubles, the last three
the angular velocity (rad/s, components in the reference frame of the segment), for a
segment with it. The fourth integer of the descriptor (see [`dafus`](@ref)) is 1 for a
segment with angular velocity. `descr` is as in [`cknr02`](@ref).

Throws `SpiceError`:
- `:ckwrongdatatype` if the segment is not of type 3;
- `:cknonexistrec` unless `1 <= recno <= cknr03(handle, descr)`;
- `:dafnosuchhandle` if the handle is closed;
- `:dafinvalidparams` if `descr` has fewer than 5 doubles.
"""
function ckgr03(h::DAFHandle, descr::AbstractVector{<:Real}, recno::Integer)
    ic = _ck_descr_ic("ckgr03", descr, 3)
    psiz = ic[4] == 1 ? 7 : 4
    nrec = _nint(dafgda(h, ic[6], ic[6])[1])
    _ck_check_recno("ckgr03", recno, nrec)
    addr = ic[5] + psiz * (recno - 1)
    tick = ic[5] + psiz * nrec + recno - 1
    record = [dafgda(h, tick, tick); dafgda(h, addr, addr + psiz - 1)]
    @assert length(record) == psiz + 1
    return record
end

# --- types 1, 4, 5 and 6 -------------------------------------------------------------------

"""
    cknr01(handle, descr) -> Int

Return the number of pointing instances in the CK type 1 segment with the descriptor
`descr` in the file of `handle`, as
[`cknr01`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/cknr01.html)
of SPICELIB. CSPICE has no C wrapper for it. `descr` is the packed descriptor of the
segment (5 doubles), as [`dafgs`](@ref) returns it during a search of the file.

The record readers of the types 1 and 4 to 6 throw `:ckwrongdatatype` for a segment of
another type, `:dafnosuchhandle` for a closed handle, `:dafinvalidparams` for a short
descriptor and, the `ckgr` readers, `:cknonexistrec` for a record number out of range.
"""
function cknr01(h::DAFHandle, descr::AbstractVector{<:Real})
    ic = _ck_descr_ic("cknr01", descr, 1)
    return _cint(dafgda(h, ic[6], ic[6])[1])   # `(integer) n`: truncation
end

"""
    ckgr01(handle, descr, recno) -> Vector{Float64}

Return the pointing instance `recno` of the CK type 1 segment with the descriptor `descr`
in the file of `handle`, as
[`ckgr01`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/ckgr01.html)
of SPICELIB. CSPICE has no C wrapper for it. Records count from 1.

The record has the layout of [`ckgr03`](@ref): the time tag (ticks), the quaternion
(SPICE form, scalar first) and, for a segment with angular velocity, the angular velocity
(rad/s), 5 or 8 doubles. The errors are those of [`cknr01`](@ref).
"""
function ckgr01(h::DAFHandle, descr::AbstractVector{<:Real}, recno::Integer)
    ic = _ck_descr_ic("ckgr01", descr, 1)
    nrec = _cint(dafgda(h, ic[6], ic[6])[1])
    psiz = ic[4] == 1 ? 7 : 4
    _ck_check_recno("ckgr01", recno, nrec)
    addr = ic[5] + psiz * (recno - 1)
    tick = ic[5] + psiz * nrec + recno - 1
    return [dafgda(h, tick, tick); dafgda(h, addr, addr + psiz - 1)]
end

# The meta data of the generic segment from `baddr` to `eaddr` in the file of `h`, as
# `sgmeta`: 17 items, the base addresses (items 1, 3, 6, 8, 11, 13) made absolute.
function _sgmeta(h::DAFHandle, baddr::Int, eaddr::Int)
    metasz = _nint(dafgda(h, eaddr, eaddr)[1])
    metasz < 15 && throw(
        SpiceError(
            :invalidmetadata, "sgmeta: this segment reports that it has $metasz meta " *
                "data items; every generic segment must have at least 15"
        )
    )
    m, _ = _sg_meta_base(metasz, eaddr)
    metasz = clamp(metasz, 16, 17)
    x = dafgda(h, m + 1, m + metasz)
    meta = Int[_nint(v) for v in x]
    length(meta) == 16 && push!(meta, 0)
    meta[17] = metasz
    for i in metasz:16      # as `sgmeta`, which also clears item 16 of a 16-item segment
        meta[i] = 0
    end
    for i in (1, 3, 6, 8, 11, 13)
        meta[i] += baddr - 1
    end
    return meta
end

"""
    cknr04(handle, descr) -> Int

Return the number of pointing records in the CK type 4 segment with the descriptor
`descr` in the file of `handle`, as
[`cknr04`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/cknr04.html)
of SPICELIB. CSPICE has no C wrapper for it. The errors are those of [`cknr01`](@ref),
and `:invalidmetadata` if the segment has fewer than 15 meta data items.
"""
function cknr04(h::DAFHandle, descr::AbstractVector{<:Real})
    ic = _ck_descr_ic("cknr04", descr, 4)
    return _sgmeta(h, ic[5], ic[6])[12]
end

"""
    ckgr04(handle, descr, recno) -> Vector{Float64}

Return the pointing record `recno` of the CK type 4 segment with the descriptor `descr`
in the file of `handle`, as
[`ckgr04`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/ckgr04.html)
of SPICELIB. CSPICE has no C wrapper for it. Records count from 1.

The record has the layout that [`ckw04a`](@ref) takes: the midpoint and the radius of
the interval (ticks), the seven coefficient counts (quaternion components q0 to q3, then
the angular velocity components), and the Chebyshev coefficients, `9 + sum(counts)`
doubles. The errors are those of [`cknr04`](@ref).
"""
function ckgr04(h::DAFHandle, descr::AbstractVector{<:Real}, recno::Integer)
    ic = _ck_descr_ic("ckgr04", descr, 4)
    meta = _sgmeta(h, ic[5], ic[6])
    nrec = meta[12]
    _ck_check_recno("ckgr04", recno, nrec)
    pdrtyp, pktoff, pktsiz, pktbas = meta[10], meta[16], meta[15], meta[11]
    if pdrtyp == 0                              # fixed-size packets, as `sgfpkt`
        b = pktbas + (recno - 1) * (pktsiz + pktoff) + 1 + pktoff
        e = b + pktsiz - 1
    else
        d = dafgda(h, meta[8] + recno, meta[8] + recno + 1)
        b = pktbas + _cint(d[1])
        e = b + _cint(d[2]) - _cint(d[1]) - pktoff - 1
    end
    pkt = dafgda(h, b, e)
    # The seven counts are packed in one double, base 128 (`zzck4d2i`).
    packed = pkt[3]
    x = 128.0^6
    counts = zeros(Int, 7)
    for k in 7:-1:1
        counts[k] = _cint(packed / x)
        packed -= counts[k] * x
        x /= 128.0
    end
    numall = sum(counts)
    return [pkt[1]; pkt[2]; Float64.(counts); pkt[4:(3 + numall)]]
end

"""
    cknr05(handle, descr) -> Int

Return the number of pointing instances in the CK type 5 segment with the descriptor
`descr` in the file of `handle`, as
[`cknr05`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/cknr05.html)
of SPICELIB. CSPICE has no C wrapper for it. The errors are those of [`cknr01`](@ref).
"""
function cknr05(h::DAFHandle, descr::AbstractVector{<:Real})
    ic = _ck_descr_ic("cknr05", descr, 5)
    return _nint(dafgda(h, ic[6], ic[6])[1])
end

"""
    ckgr05(handle, descr, recno) -> Vector{Float64}

Return the pointing instance `recno` of the CK type 5 segment with the descriptor `descr`
in the file of `handle`, as
[`ckgr05`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/ckgr05.html)
of SPICELIB. CSPICE has no C wrapper for it. Records count from 1.

The record is the time tag (ticks), the subtype of the segment and the packet of the
subtype (8, 4, 14 or 7 doubles for the subtypes 0 to 3; see [`ckw05`](@ref)). The errors
are those of [`cknr01`](@ref), and `:notsupported` for a subtype other than 0 to 3.
"""
function ckgr05(h::DAFHandle, descr::AbstractVector{<:Real}, recno::Integer)
    ic = _ck_descr_ic("ckgr05", descr, 5)
    beg, en = ic[5], ic[6]
    sub = dafgda(h, en - 3, en - 3)[1]
    subtyp = _nint(sub)
    0 <= subtyp <= 3 || throw(
        SpiceError(:notsupported, "ckgr05: unexpected CK type 5 subtype $subtyp found in type 5 segment")
    )
    packsz = _CK5_PACKSZ[subtyp + 1]
    nrec = _nint(dafgda(h, en, en)[1])
    _ck_check_recno("ckgr05", recno, nrec)
    addr = beg + packsz * (recno - 1)
    tick = beg + packsz * nrec + recno - 1
    return [dafgda(h, tick, tick); sub; dafgda(h, addr, addr + packsz - 1)]
end

"""
    cknm06(handle, descr) -> Int

Return the number of mini-segments in the CK type 6 segment with the descriptor `descr`
in the file of `handle`, as
[`cknm06`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/cknm06.html)
of SPICELIB. CSPICE has no C wrapper for it. The errors are those of [`cknr01`](@ref).

See also [`ckmp06`](@ref), [`ckgr06`](@ref).
"""
function cknm06(h::DAFHandle, descr::AbstractVector{<:Real})
    ic = _ck_descr_ic("cknm06", descr, 6)
    return _nint(dafgda(h, ic[6], ic[6])[1])
end

function _ck06_check_msno(who, msno, nintvl)
    1 <= msno <= nintvl || throw(
        SpiceError(:indexoutofrange, "$who: mini-segment index must be in range 1:$nintvl but was $msno")
    )
    return nothing
end

"""
    ckmp06(handle, descr, msno) -> (rate, subtyp, winsiz, nrec, ivlbds, lstepc)

Return the parameters of the mini-segment `msno` of the CK type 6 segment with the
descriptor `descr` in the file of `handle`, as
[`ckmp06`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/ckmp06.html)
of SPICELIB. CSPICE has no C wrapper for it. Mini-segments count from 1.

`rate` is the clock rate (seconds per tick), `subtyp` the subtype, `winsiz` the
interpolation window size, `nrec` the number of records (these three are `Int`), `ivlbds`
a vector of the two bounds of the interval of the mini-segment (ticks) and `lstepc` the
last epoch of its records (ticks).

Throws `SpiceError(:indexoutofrange)` unless `1 <= msno <= cknm06(handle, descr)`; the
other errors are those of [`cknr01`](@ref).
"""
function ckmp06(h::DAFHandle, descr::AbstractVector{<:Real}, msno::Integer)
    ic = _ck_descr_ic("ckmp06", descr, 6)
    baddr, eaddr = ic[5], ic[6]
    nintvl = _nint(dafgda(h, eaddr, eaddr)[1])
    _ck06_check_msno("ckmp06", msno, nintvl)
    ptrbas = eaddr - 2 - (nintvl + 1)
    ivlbas = ptrbas - nintvl ÷ 100 - (nintvl + 1)
    ivlbds = dafgda(h, ivlbas + msno, ivlbas + msno + 1)
    minie = baddr - 1 + _nint(dafgda(h, ptrbas + msno + 1, ptrbas + msno + 1)[1]) - 1
    buf = dafgda(h, minie - 3, minie)
    nrec = _nint(buf[4])
    epaddr = minie - 4 - (nrec - 1) ÷ 100
    return buf[1], _nint(buf[2]), _nint(buf[3]), nrec, ivlbds, dafgda(h, epaddr, epaddr)[1]
end

"""
    ckgr06(handle, descr, msno, recno) -> Vector{Float64}

Return the record `recno` of the mini-segment `msno` of the CK type 6 segment with the
descriptor `descr` in the file of `handle`, as
[`ckgr06`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/ckgr06.html)
of SPICELIB. CSPICE has no C wrapper for it. Mini-segments and records count from 1.

The record is the time tag (ticks), the subtype, the clock rate (seconds per tick) and the
packet of the subtype (8, 4, 14 or 7 doubles for the subtypes 0 to 3; see
[`ckw06`](@ref)).

Throws `SpiceError`:
- `:indexoutofrange` unless `1 <= msno <= cknm06(handle, descr)`;
- `:notsupported` for a subtype other than 0 to 3;
- `:cknonexistrec` unless `1 <= recno <= nrec` of [`ckmp06`](@ref);
- the errors of [`cknr01`](@ref).
"""
function ckgr06(h::DAFHandle, descr::AbstractVector{<:Real}, msno::Integer, recno::Integer)
    ic = _ck_descr_ic("ckgr06", descr, 6)
    baddr, eaddr = ic[5], ic[6]
    nintvl = _nint(dafgda(h, eaddr, eaddr)[1])
    _ck06_check_msno("ckgr06", msno, nintvl)
    ptrbas = eaddr - 2 - (nintvl + 1)
    ptr = dafgda(h, ptrbas + msno, ptrbas + msno + 1)
    minib = baddr - 1 + _nint(ptr[1])
    minie = baddr - 1 + _nint(ptr[2]) - 1
    buf = dafgda(h, minie - 3, minie)
    rate, subtyp, nrec = buf[1], _nint(buf[2]), _nint(buf[4])
    0 <= subtyp <= 3 || throw(
        SpiceError(:notsupported, "ckgr06: unexpected CK type 6 subtype $subtyp found in mini-segment $msno")
    )
    pktsiz = _CK5_PACKSZ[subtyp + 1]
    1 <= recno <= nrec || throw(
        SpiceError(:cknonexistrec, "ckgr06: record index must be in range 1:$nrec but was $recno")
    )
    epaddr = minie - 4 - (nrec - 1) ÷ 100 - nrec + recno
    b = minib + (recno - 1) * pktsiz
    return [dafgda(h, epaddr, epaddr); Float64(subtyp); rate; dafgda(h, b, b + pktsiz - 1)]
end
