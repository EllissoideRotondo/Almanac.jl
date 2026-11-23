# Running an encoded EK query on a `KernelSet`: the search of `eksrch` (`ekqmgr.c`).
#
# The WHERE clause of an encoded query is a disjunction of conjunctions of constraints.
# For one conjunction the search does this, as CSPICE:
#
#   1. For each table of the FROM clause and each loaded segment of the table (load
#      order), find the rows that satisfy the constraints on that table alone. If a
#      constraint compares an indexed column with a value, the index of one such column
#      (the "key") gives a range of candidate rows, which are then visited in index
#      order; otherwise all rows are visited in row order (`zzekkey`).
#   2. Join the tables from left to right (`zzekjoin`): the rows of tables 1 to t-1 with
#      the rows of table t, under the constraints that compare columns of two of these
#      tables. An equality constraint gives a sort-merge join, an inequality a join over
#      two sorted lists, else a nested loop.
#
# The results of the conjunctions are put together; a row combination that an earlier
# conjunction gave is removed (`zzekweed`). Then the ORDER BY clause sorts the result
# with the Shell sort of `zzekjsrt`, so equal keys come out in the order of CSPICE.
#
# CSPICE keeps these row sets in a global scratch file. Here they are local arrays, and
# the search is a pure function of the snapshot and the query: it runs on many threads.
#
# Four defects of CSPICE N0067 change results; the search has a switch `compat` that
# reproduces them, and the tests use it to show that every difference from CSPICE comes
# from these defects. With `compat = false` (the public calls) the results are correct.
#
#   (a) `zzekweed` compares each candidate row with the rows of every segment vector of
#       an earlier conjunction (it reads the candidate's segment vector twice), and the
#       comparison includes the position of the segment vector in its row set. So a row
#       that two conjunctions give can stay twice in the result (the positions differ),
#       and two different rows with the same record pointer in two segments can be taken
#       for one (the positions agree).
#   (b) `zzekillt`/`zzekille` round a double precision value to the nearest integer when
#       the key column has integer type, so `ID < 17.4` on an indexed column loses the
#       rows with `ID = 17`; for a double precision key column and a value that is not of
#       integer or double precision type they use an integer that is not set.
#   (c) When `zzekweed` has deleted all rows of one segment vector of a conjunction,
#       `zzekjsqz` removes that segment vector but keeps its row list, so each later
#       segment vector of the conjunction gets the rows of the one before it: the rows
#       of the next segment are lost, and the rows after it are read from a wrong
#       segment.
#   (d) `zzekqcon` converts an integer value of a query with `i_dnnt`, which gives the
#       smallest 32-bit integer for a value outside the 32-bit range, and `zzekscmp`
#       compares with that integer: `DCOL > 10000000000` holds for every row that is
#       not null. The public search compares such a value as a double.
#
# A fifth difference is not a defect of the search. The scanner of CSPICE reads a number
# of a query with `nparsd`, whose result can be one unit in the last place away from the
# correctly rounded double ("2.756" gives the double after 2.756), so that a row with
# the value 2.756 is not equal to the number 2.756 of a query. The public calls compare
# with the correctly rounded double (`dexact` of a constraint); with `compat` the search
# compares with the double of `nparsd` (`dval`).

# Operator codes (`ekopcd.inc`) and the codes of a comparison result.
const _EK_EQ = 1
const _EK_GE = 2
const _EK_GT = 3
const _EK_LE = 4
const _EK_LT = 5
const _EK_NE = 6
const _EK_LIKE = 7
const _EK_UNLIKE = 8
const _EK_ISNULL = 9
const _EK_NOTNUL = 10
const _EK_MAXSTR = 1024           # characters of a string that take part in a comparison

"One segment as the search uses it: its file, the segment, and its columns in table order."
struct _EKQSeg
    file::EKFile
    seg::_EKSegment
    cols::Vector{_EKColumn}
end

# ---------------------------------------------------------------------------------------
# Comparisons
# ---------------------------------------------------------------------------------------

"""
    _ek_bcmp(a, b, na, nb) -> Int

Order of the first `na` bytes of `a` and the first `nb` bytes of `b` as Fortran compares
strings (`l_lt`, `l_gt`): byte by byte, the shorter one continued with blanks. Returns
-1, 0 or 1.
"""
function _ek_bcmp(a::AbstractVector{UInt8}, b::AbstractVector{UInt8}, na::Int, nb::Int)
    for i in 1:max(na, nb)               # invariant: the characters before i are equal
        x = i <= na ? a[i] : UInt8(' ')
        y = i <= nb ? b[i] : UInt8(' ')
        x == y || return x < y ? -1 : 1
    end
    return 0
end

"Relation code (`_EK_LT`, `_EK_EQ`, `_EK_GT`) of two strings; the first 1024 characters count."
function _ek_fcmp(a::AbstractString, b::AbstractString)
    ca = codeunits(a)
    cb = codeunits(b)
    r = _ek_bcmp(ca, cb, min(length(ca), _EK_MAXSTR), min(length(cb), _EK_MAXSTR))
    return r < 0 ? _EK_LT : r > 0 ? _EK_GT : _EK_EQ
end

@inline _ek_ncmp(x, y) = x < y ? _EK_LT : x > y ? _EK_GT : _EK_EQ

@inline function _ek_relop(op::Int, rel::Int)
    op == _EK_EQ && return rel == _EK_EQ
    op == _EK_LT && return rel == _EK_LT
    op == _EK_LE && return rel != _EK_GT
    op == _EK_GT && return rel == _EK_GT
    op == _EK_GE && return rel != _EK_LT
    op == _EK_NE && return rel != _EK_EQ
    throw(SpiceError(:unnaturalrelation, "The relational operator $op was not recognized."))
end

@noinline _ek_noelt(col::_EKColumn, recptr::Integer) = throw(
    SpiceError(
        :invalidindex,
        "COLUMN = $(col.name); ROW = $recptr; ELTIDX = 1. Column entry element was not found."
    )
)

"The first `_EK_MAXSTR` characters of the element `r` of the store `f` as a `String`."
_ek_ctrunc(f, r::_EKCRef) = _ek_cstring(f.chars, _ek_cfirst(r, _EK_MAXSTR))

"""
    _ek_vrel(f, col, recptr, dtype, cval, dval, ival, keynull, writer) -> (rel, enull, cref)

Relation code (`_EK_LT`, `_EK_EQ`, `_EK_GT`) of the first element of the entry of
column `col` in the record `recptr` of the store `f` to a value of type code `dtype`
(the bytes `cval`, `dval` or `ival`), as `zzekscmp` finds it; `enull` tells that the
entry is null, and `cref` is the place of the element of a character column. A null
entry is less than every value. With `keynull` the value is null, which is equal to a
null entry and less than the others (the index searches of the writers). Of a string,
the first 1024 characters count.

Throws `SpiceError(:bug)` for types that cannot be compared, `:invalidindex` for an
entry without elements, and the errors of the readers. A column type other than the
four is read as integer, or with `writer` throws `:invaliddatatype`.
"""
function _ek_vrel(
        f, col::_EKColumn, recptr::Integer, dtype::Int, cval::AbstractVector{UInt8},
        dval::Float64, ival::Int, keynull::Bool, writer::Bool
    )
    coltyp = Int(col.dtype)
    cref = _EKCRef(0)
    eltd = 0.0
    elti = 0
    if coltyp == _EK_CHR
        cref, enull, found = _ek_cref(f, col, recptr, 1)
    elseif coltyp == _EK_DP || coltyp == _EK_TIME
        eltd, enull, found = _ek_elt_d(f, col, recptr, 1)
    elseif coltyp == _EK_INT || !writer
        elti, enull, found = _ek_elt_i(f, col, recptr, 1)
    else
        throw(SpiceError(:invaliddatatype, "Data type code $coltyp not recognized."))
    end
    found || _ek_noelt(col, recptr)
    enull && return (keynull ? _EK_EQ : _EK_LT), true, cref
    keynull && return _EK_GT, false, cref
    badtype() = throw(SpiceError(:bug, "Column type is $coltyp; value type is $dtype."))
    rel = if coltyp == _EK_CHR
        dtype == _EK_CHR || badtype()
        r = _ek_ccmp(
            f.chars, _ek_cfirst(cref, _EK_MAXSTR), cval,
            _EKCRef(min(length(cval), _EK_MAXSTR))
        )
        r < 0 ? _EK_LT : r > 0 ? _EK_GT : _EK_EQ
    elseif coltyp == _EK_TIME
        (dtype == _EK_TIME || dtype == _EK_DP) || badtype()
        _ek_ncmp(eltd, dval)
    elseif coltyp == _EK_DP
        dtype == _EK_CHR && badtype()
        _ek_ncmp(eltd, dtype == _EK_INT ? Float64(ival) : dval)
    else
        (dtype == _EK_INT || dtype == _EK_DP) || badtype()
        _ek_ncmp(Float64(elti), dtype == _EK_INT ? Float64(ival) : dval)
    end
    return rel, false, cref
end

"True if the entry of `col` in the record `recptr` is null; `:invalidindex` if it has no elements."
function _ek_isnull(f, col::_EKColumn, recptr::Integer)
    coltyp = Int(col.dtype)
    _, enull, found = coltyp == _EK_CHR ? _ek_cref(f, col, recptr, 1) :
        (coltyp == _EK_DP || coltyp == _EK_TIME) ? _ek_elt_d(f, col, recptr, 1) :
        _ek_elt_i(f, col, recptr, 1)
    found || _ek_noelt(col, recptr)
    return enull
end

"""
    _ek_scmp(op, f, col, recptr, dtype, cval, dval, ival) -> Bool

Port of `zzekscmp`: true if the first element of the entry of column `col` in the record
`recptr` of the store `f` has the relation `op` to a value of type code `dtype`
(`cval`, `dval` or `ival`), with the rules of [`_ek_vrel`](@ref). LIKE and UNLIKE are
false for a null entry. Throws as `_ek_vrel`, and `:unnaturalrelation` for an operator
that the types do not have.
"""
function _ek_scmp(
        op::Int, f, col::_EKColumn, recptr::Integer, dtype::Int, cval::String,
        dval::Float64, ival::Int
    )
    op == _EK_ISNULL && return _ek_isnull(f, col, recptr)
    op == _EK_NOTNUL && return !_ek_isnull(f, col, recptr)
    rel, enull, cref = _ek_vrel(f, col, recptr, dtype, codeunits(cval), dval, ival, false, false)
    op <= _EK_NE && return _ek_relop(op, rel)
    if (op == _EK_LIKE || op == _EK_UNLIKE) && dtype == _EK_CHR
        enull && return false
        m = _matchi(_ek_ctrunc(f, cref), cval)
        return op == _EK_LIKE ? m : !m
    end
    throw(
        SpiceError(
            :unnaturalrelation, "The relational operator $op was not recognized " *
                "or was not applicable for data type $dtype."
        )
    )
end

"""
    _ek_ecmp(f1, col1, row1, f2, col2, row2) -> Int

Port of `zzekecmp`: relation code of the first elements of two column entries. Two null
entries are equal, and a null entry is less than an entry that is not null. An integer
can be compared with a double; time only with time and character only with character
(`SpiceError(:bug)` otherwise).
"""
function _ek_ecmp(f1, col1::_EKColumn, row1::Integer, f2, col2::_EKColumn, row2::Integer)
    lt = Int(col1.dtype)
    rt = Int(col2.dtype)
    bad() = throw(SpiceError(:bug, "LHS data type is $lt; RHSTYP is $rt."))
    if lt == _EK_CHR
        rt == _EK_CHR || bad()
        a, n1, fnd1 = _ek_cref(f1, col1, row1, 1)
        b, n2, fnd2 = _ek_cref(f2, col2, row2, 1)
        fnd1 || _ek_noelt(col1, row1)
        fnd2 || _ek_noelt(col2, row2)
        (n1 || n2) && return n1 && n2 ? _EK_EQ : n1 ? _EK_LT : _EK_GT
        r = _ek_ccmp(
            f1.chars, _ek_cfirst(a, _EK_MAXSTR), f2.chars, _ek_cfirst(b, _EK_MAXSTR)
        )
        return r < 0 ? _EK_LT : r > 0 ? _EK_GT : _EK_EQ
    end
    if lt == _EK_TIME
        rt == _EK_TIME || bad()
    else
        (rt == _EK_INT || rt == _EK_DP) || bad()
    end
    x, n1 = _ek_number(f1, col1, row1)
    y, n2 = _ek_number(f2, col2, row2)
    (n1 || n2) && return n1 && n2 ? _EK_EQ : n1 ? _EK_LT : _EK_GT
    (lt == _EK_INT && rt == _EK_INT) && return _ek_ncmp(Int(x), Int(y))
    return _ek_ncmp(x, y)
end

@inline function _ek_number(f, col::_EKColumn, row::Integer)
    if col.dtype == _EK_INT
        v, isnull, found = _ek_elt_i(f, col, row, 1)
        found || _ek_noelt(col, row)
        return Float64(v), isnull
    end
    v, isnull, found = _ek_elt_d(f, col, row, 1)
    found || _ek_noelt(col, row)
    return v, isnull
end

"""
    _ek_colmatch(op, f1, col1, row1, f2, col2, row2) -> Bool

One constraint between two columns, as `zzekvmch`: the relation of `_ek_ecmp`, or for
LIKE and UNLIKE the match of the left string against the right string as template (a
null entry counts as a blank string).
"""
function _ek_colmatch(
        op::Int, f1, col1::_EKColumn, row1::Integer, f2, col2::_EKColumn,
        row2::Integer
    )
    op <= _EK_NE && return _ek_relop(op, _ek_ecmp(f1, col1, row1, f2, col2, row2))
    ((op == _EK_LIKE || op == _EK_UNLIKE) && col1.dtype == _EK_CHR) || throw(
        SpiceError(
            :unnaturalrelation, "The relational operator $op was not recognized."
        )
    )
    a, n1, fnd1 = _ek_cref(f1, col1, row1, 1)
    b, n2, fnd2 = _ek_cref(f2, col2, row2, 1)
    fnd1 || _ek_noelt(col1, row1)
    fnd2 || _ek_noelt(col2, row2)
    m = _matchi(n1 ? "" : _ek_ctrunc(f1, a), n2 ? "" : _ek_ctrunc(f2, b))
    return op == _EK_LIKE ? m : !m
end

# ---------------------------------------------------------------------------------------
# Indexes
# ---------------------------------------------------------------------------------------

"""
    _ek_ixlk(f, col, key, ix = _EK_NOIX) -> Int

Record pointer at position `key` of the index of column `col`, as `zzekixlk`. An index
of type 1 is a tree; an index of type 2 (fast writers) is an order vector in adjacent
integer pages. `ix` is the index decoded at load ([`_ek_ixvals`](@ref)); if it is
empty the file is read. Throws `SpiceError(:invalidtype)` for another index type.
"""
function _ek_ixlk(f, col::_EKColumn, key::Int, ix::Vector{Int32} = _EK_NOIX)
    isempty(ix) || return Int(ix[key])
    col.ixtype == 1 && return _ek_trdp(f.ints, col.ixptr, key)
    col.ixtype == 2 || throw(
        SpiceError(
            :invalidtype,
            "The index type $(col.ixtype) is not supported."
        )
    )
    return _ek_ordvec_value(f.ints, col.ixptr, key)
end

"""
    _ek_ixvals(qs, col) -> Vector{Int32}

The record pointers of the index of the column `col` of the segment of `qs` in index
order, as `_ek_segment` decoded them at load; empty if they were not decoded.
"""
function _ek_ixvals(qs::_EKQSeg, col::_EKColumn)
    cols = qs.seg.columns
    j = Int(col.ord)
    (1 <= j <= length(cols) && cols[j] == col) || (j = something(findfirst(==(col), cols), 0))
    return j == 0 ? _EK_NOIX : qs.seg.ixvals[j]
end

"""
    _ek_ixlast(f, col, nrows, strict, dtype, cval, dval, ival, ix = _EK_NOIX) -> Int

Position of the last entry in the index of `col` that is less than (`strict`) or not
greater than the value, 0 if there is none: the binary search of `zzeklltc`/`zzekllec`
and their integer and double forms. `ix` as in [`_ek_ixlk`](@ref). Precondition: `col`
has an index over `nrows` rows.
"""
function _ek_ixlast(
        f, col::_EKColumn, nrows::Int, strict::Bool, dtype::Int, cval::String,
        dval::Float64, ival::Int, ix::Vector{Int32} = _EK_NOIX
    )
    nrows == 0 && return 0
    notbelow = strict ? _EK_GE : _EK_GT       # the first entry is too large
    below = strict ? _EK_LT : _EK_LE
    cmp(op, pos) = _ek_scmp(op, f, col, _ek_ixlk(f, col, pos, ix), dtype, cval, dval, ival)
    cmp(notbelow, 1) && return 0
    cmp(below, nrows) && return nrows
    b = 1
    e = nrows
    while e > b + 1                      # invariant: entry b is below the value, entry e is not
        m = (b + e) ÷ 2
        if cmp(below, m)
            b = m
        else
            e = m
        end
    end
    return b
end

"""
    _ek_ixbound(f, col, nrows, strict, con, compat, ix) -> Int

`zzekillt` (`strict`) or `zzekille` for the value of the constraint `con` and the key
column `col`, whose decoded index is `ix`. With `compat` the value is converted as CSPICE converts it (see defect
(b) in the file header); else it is compared as it is.
"""
function _ek_ixbound(
        f, col::_EKColumn, nrows::Int, strict::Bool, con, compat::Bool, ix::Vector{Int32}
    )
    coltyp = Int(col.dtype)
    dtype = Int(con.dtype)
    coltyp == _EK_CHR && return _ek_ixlast(f, col, nrows, strict, _EK_CHR, con.cval, 0.0, 0, ix)
    if !compat
        vt = dtype == _EK_INT ? _EK_INT : coltyp == _EK_TIME ? _EK_TIME : _EK_DP
        return _ek_ixlast(f, col, nrows, strict, vt, "", con.dval, con.ival, ix)
    end
    if coltyp == _EK_DP
        dnum = dtype == _EK_DP ? con.dval : Float64(con.ival)
        return _ek_ixlast(f, col, nrows, strict, _EK_DP, "", dnum, 0, ix)
    elseif coltyp == _EK_TIME
        return _ek_ixlast(f, col, nrows, strict, _EK_DP, "", con.dval, 0, ix)
    end
    inum = dtype == _EK_DP ? round(Int, con.dval, RoundNearestTiesAway) : con.ival
    return _ek_ixlast(f, col, nrows, strict, _EK_INT, "", 0.0, inum, ix)
end

"""
    _ek_key(qs, cons, active, compat) -> (found, keycol, begidx, endidx)

Port of `zzekkey`: choose the key column of the segment `qs` for the value constraints
`cons` that are `active`. Among the indexed columns that have an active constraint (in
column order), the key is the one whose constraints with EQ, LT, LE, GT, GE leave the
smallest range `begidx:endidx` of index positions; a range that is empty ends the
search. Postcondition: if `found`, the EQ, LT, LE, GT and GE constraints on the key
column are not active any more, because the range satisfies them.
"""
function _ek_key(qs::_EKQSeg, cons::Vector, active::Vector{Bool}, compat::Bool)
    nrows = Int(qs.seg.nrows)
    f = qs.file
    idxset = Int[]
    for (i, c) in enumerate(cons)
        (active[i] && qs.cols[c.lhs.colidx].ixtype != -1) && push!(idxset, c.lhs.colidx)
    end
    isempty(idxset) && return false, qs.cols[1], 1, nrows
    idxset = sort!(unique!(idxset))
    begidx = 1
    endidx = nrows
    best = idxset[1]
    elim = false
    col = best
    for outer col in idxset
        kc = qs.cols[col]
        ix = _ek_ixvals(qs, kc)
        minptr = 1
        maxptr = nrows
        for (i, c) in enumerate(cons)    # invariant: minptr:maxptr satisfies constraints before i
            (c.lhs.colidx == col && active[i]) || continue
            op = c.opcode
            if op == _EK_LT
                maxptr = min(_ek_ixbound(f, kc, nrows, true, c, compat, ix), maxptr)
                elim = maxptr == 0
            elseif op == _EK_LE
                maxptr = min(_ek_ixbound(f, kc, nrows, false, c, compat, ix), maxptr)
                elim = maxptr == 0
            elseif op == _EK_EQ
                lastlt = _ek_ixbound(f, kc, nrows, true, c, compat, ix)
                lastle = _ek_ixbound(f, kc, nrows, false, c, compat, ix)
                if lastlt < lastle
                    minptr = max(lastlt + 1, minptr)
                    maxptr = min(lastle, maxptr)
                else
                    elim = true
                end
            elseif op == _EK_GT
                lastle = _ek_ixbound(f, kc, nrows, false, c, compat, ix)
                minptr = max(lastle + 1, minptr)
                elim = lastle == nrows
            elseif op == _EK_GE
                lastlt = _ek_ixbound(f, kc, nrows, true, c, compat, ix)
                minptr = max(lastlt + 1, minptr)
                elim = lastlt == nrows
            end
            minptr > maxptr && (elim = true)
            elim && break
        end
        elim && break
        if maxptr - minptr + 1 < endidx - begidx + 1
            best = col
            begidx = minptr
            endidx = maxptr
        end
    end
    key = best
    if elim
        key = col
        begidx = 1
        endidx = 0
    end
    for (i, c) in enumerate(cons)
        if active[i] && c.lhs.colidx == key && c.opcode in (_EK_LT, _EK_LE, _EK_EQ, _EK_GE, _EK_GT)
            active[i] = false
        end
    end
    return true, qs.cols[key], begidx, endidx
end

# ---------------------------------------------------------------------------------------
# Row sets
# ---------------------------------------------------------------------------------------

"""
    _EKRowSet

A "join row set" of CSPICE: row combinations of `nt` tables. `segvecs[k]` names one
segment for each table (indices into the segment list of the search), and `rows[k]`
holds the row vectors that belong to it: `nt` record pointers for each row vector, one
after the other. The order of the segment vectors and of the row vectors is the order
of the result.
"""
struct _EKRowSet
    nt::Int
    segvecs::Vector{Vector{Int}}
    rows::Vector{Vector{Int32}}
end

_ek_nrows(rs::_EKRowSet) = sum(r -> length(r) ÷ rs.nt, rs.rows; init = 0)

"Remove the segment vectors without rows from `rs`, as `zzekjsqz`."
function _ek_squeeze!(rs::_EKRowSet)
    keep = findall(!isempty, rs.rows)
    length(keep) == length(rs.rows) && return rs
    kept_sv = rs.segvecs[keep]
    kept_rows = rs.rows[keep]
    resize!(rs.segvecs, length(keep)) .= kept_sv
    resize!(rs.rows, length(keep)) .= kept_rows
    return rs
end

"""
    _ek_shellsort!(ord, jle) -> ord

The Shell sort of `zzekjsrt` on the order vector `ord`: `jle(a, b)` tells that item `a`
can stay before item `b`. The sequence of comparisons and swaps is that of CSPICE, so
items with equal keys end in the same order.
"""
function _ek_shellsort!(ord::Vector{Int}, jle)
    n = length(ord)
    gap = n ÷ 2
    while gap > 0
        for i in (gap + 1):n
            j = i - gap
            while j > 0
                jg = j + gap
                if jle(ord[j], ord[jg])
                    j = 0
                else
                    ord[j], ord[jg] = ord[jg], ord[j]
                end
                j -= gap
            end
        end
        gap ÷= 2
    end
    return ord
end

"""
    _EKKeys

The first elements of one column for a list of rows, read once so that a sort does not
read the file for each comparison: `null[i]`, and the value in `ints` (integer column),
`dps` (double or time) or `strs` (character, at most 1024 characters). `dtype` is the
data type code of the column.
"""
struct _EKKeys
    dtype::Int
    null::Vector{Bool}
    ints::Vector{Int}
    dps::Vector{Float64}
    strs::Vector{String}
end

_EKKeys(dtype::Int) = _EKKeys(dtype, Bool[], Int[], Float64[], String[])

function _ek_pushkey!(keys::_EKKeys, f, col::_EKColumn, recptr::Integer)
    col.dtype == keys.dtype || throw(
        SpiceError(
            :bug, "Data type mismatch for the column " *
                "$(col.name): type $(col.dtype) in one segment, $(keys.dtype) in another"
        )
    )
    if keys.dtype == _EK_CHR
        r, isnull, found = _ek_cref(f, col, recptr, 1)
        found || _ek_noelt(col, recptr)
        push!(keys.strs, _ek_ctrunc(f, r))
    elseif keys.dtype == _EK_INT
        v, isnull, found = _ek_elt_i(f, col, recptr, 1)
        found || _ek_noelt(col, recptr)
        push!(keys.ints, v)
    else
        v, isnull, found = _ek_elt_d(f, col, recptr, 1)
        found || _ek_noelt(col, recptr)
        push!(keys.dps, v)
    end
    push!(keys.null, isnull)
    return keys
end

"Relation code of item `i` of `a` and item `j` of `b`, with the rules of `_ek_ecmp`."
function _ek_keycmp(a::_EKKeys, i::Int, b::_EKKeys, j::Int)
    n1 = a.null[i]
    n2 = b.null[j]
    (n1 || n2) && return n1 && n2 ? _EK_EQ : n1 ? _EK_LT : _EK_GT
    if a.dtype == _EK_CHR
        b.dtype == _EK_CHR || throw(SpiceError(:bug, "LHS data type is 1; RHSTYP is $(b.dtype)."))
        return _ek_fcmp(a.strs[i], b.strs[j])
    end
    ok = a.dtype == _EK_TIME ? b.dtype == _EK_TIME : (b.dtype == _EK_INT || b.dtype == _EK_DP)
    ok || throw(SpiceError(:bug, "LHS data type is $(a.dtype); RHSTYP is $(b.dtype)."))
    (a.dtype == _EK_INT && b.dtype == _EK_INT) && return _ek_ncmp(a.ints[i], b.ints[j])
    x = a.dtype == _EK_INT ? Float64(a.ints[i]) : a.dps[i]
    y = b.dtype == _EK_INT ? Float64(b.ints[j]) : b.dps[j]
    return _ek_ncmp(x, y)
end

# ---------------------------------------------------------------------------------------
# One table
# ---------------------------------------------------------------------------------------

"""
    _ek_table_rows(segs, segids, t, cons, compat) -> _EKRowSet

The rows of table number `t` of the FROM clause that satisfy the constraints of the
conjunction `cons` on that table alone (step 1 of the file header). `segids` are the
segments of the table in load order, as indices into `segs`.
"""
function _ek_table_rows(
        segs::Vector{_EKQSeg}, segids::Vector{Int}, t::Int, cons::Vector,
        compat::Bool
    )
    rs = _EKRowSet(1, Vector{Int}[], Vector{Int32}[])
    n = length(cons)
    activc = Vector{Bool}(undef, n)
    activv = Vector{Bool}(undef, n)
    for sid in segids
        qs = segs[sid]
        f = qs.file
        for (i, c) in enumerate(cons)
            activc[i] = c.cnstyp == 1 && c.lhs.tabidx == t && c.rhs.tabidx == t
            activv[i] = c.cnstyp == 2 && c.lhs.tabidx == t
        end
        keyfnd, keycol, begidx, endidx = _ek_key(qs, cons, activv, compat)
        ix = keyfnd ? _ek_ixvals(qs, keycol) : _EK_NOIX
        rows = Int32[]
        for r in begidx:endidx           # invariant: `rows` holds the matches before position r
            recptr = keyfnd ? _ek_ixlk(f, keycol, r, ix) : Int(qs.seg.recptrs[r])
            ok = true
            for (i, c) in enumerate(cons)
                activv[i] || continue
                ok = _ek_scmp(
                    c.opcode, f, qs.cols[c.lhs.colidx], recptr, c.dtype, c.cval,
                    c.dval, c.ival
                )
                ok || break
            end
            ok || continue
            for (i, c) in enumerate(cons)
                activc[i] || continue
                ok = _ek_colmatch(
                    c.opcode, f, qs.cols[c.lhs.colidx], recptr, f,
                    qs.cols[c.rhs.colidx], recptr
                )
                ok || break
            end
            ok && push!(rows, Int32(recptr))
        end
        push!(rs.segvecs, [sid])
        push!(rs.rows, rows)
    end
    return _ek_squeeze!(rs)
end

# ---------------------------------------------------------------------------------------
# Join
# ---------------------------------------------------------------------------------------

"Keys of the column `colidx` of table `tab` for the `nr` row vectors `rows` (`nt` pointers each) of one segment vector."
function _ek_join_keys(
        segs::Vector{_EKQSeg}, segvec::Vector{Int}, tab::Int, colidx::Int,
        rows::Vector{Int32}, nt::Int, off::Int
    )
    qs = segs[segvec[tab]]
    col = qs.cols[colidx]
    keys = _EKKeys(Int(col.dtype))
    for i in 1:(length(rows) ÷ nt)
        _ek_pushkey!(keys, qs.file, col, rows[(i - 1) * nt + tab - off])
    end
    return keys
end

"""
    _ek_join(segs, rs1, rs2, cons, active) -> _EKRowSet

Port of `zzekjoin` with `zzekjprp`/`zzekjnxt`: the row combinations of `rs1` (tables 1
to `nt1`) and `rs2` (the next table) that satisfy the `active` constraints of `cons`,
which compare a column of one row set with a column of the other. The order of the
result is that of CSPICE: segment vectors of `rs1` outside, of `rs2` inside; in a pair
of segment vectors, the order that the join method gives (see the file header).
"""
function _ek_join(
        segs::Vector{_EKQSeg}, rs1::_EKRowSet, rs2::_EKRowSet, cons::Vector,
        active::Vector{Bool}
    )
    nt1 = rs1.nt
    nt2 = rs2.nt
    nt3 = nt1 + nt2
    out = _EKRowSet(nt3, Vector{Int}[], Vector{Int32}[])
    # The constraint that selects the join method: the first equality, else the first
    # inequality with LT, LE, GE or GT.
    equal = findfirst(i -> active[i] && cons[i].opcode == _EK_EQ, eachindex(cons))
    ordered = equal !== nothing ? nothing :
        findfirst(
            i -> active[i] && cons[i].opcode in (_EK_LT, _EK_LE, _EK_GE, _EK_GT),
            eachindex(cons)
        )
    cnstr = equal !== nothing ? equal : ordered
    case = equal !== nothing ? 1 : ordered !== nothing ? 2 : 3
    rowvec = Vector{Int32}(undef, nt3)
    for (s1, sv1) in enumerate(rs1.segvecs), (s2, sv2) in enumerate(rs2.segvecs)
        segvec = vcat(sv1, sv2)
        rows1 = rs1.rows[s1]
        rows2 = rs2.rows[s2]
        nr1 = length(rows1) ÷ nt1
        nr2 = length(rows2) ÷ nt2
        res = Int32[]
        # All active constraints except `skip` hold for the combination in `rowvec`.
        function others(skip::Int)
            for (j, c) in enumerate(cons)
                (active[j] && j != skip) || continue
                l = c.lhs.tabidx
                r = c.rhs.tabidx
                ql = segs[segvec[l]]
                qr = segs[segvec[r]]
                _ek_colmatch(
                    c.opcode, ql.file, ql.cols[c.lhs.colidx], rowvec[l], qr.file,
                    qr.cols[c.rhs.colidx], rowvec[r]
                ) || return false
            end
            return true
        end
        function candidate(i1::Int, i2::Int, skip::Int)
            copyto!(rowvec, 1, rows1, (i1 - 1) * nt1 + 1, nt1)
            copyto!(rowvec, nt1 + 1, rows2, (i2 - 1) * nt2 + 1, nt2)
            others(skip) && append!(res, rowvec)
            return nothing
        end
        if nr1 == 0 || nr2 == 0
            # nothing to join
        elseif case == 3
            for i1 in 1:nr1, i2 in 1:nr2
                candidate(i1, i2, 0)
            end
        else
            c = cons[cnstr]
            # "Left" is the row set that holds the table of the left side of the constraint.
            lfirst = c.lhs.tabidx <= nt1
            lkeys = lfirst ? _ek_join_keys(segs, segvec, c.lhs.tabidx, c.lhs.colidx, rows1, nt1, 0) :
                _ek_join_keys(segs, segvec, c.lhs.tabidx, c.lhs.colidx, rows2, nt2, nt1)
            rfirst = c.rhs.tabidx <= nt1
            rkeys = rfirst ? _ek_join_keys(segs, segvec, c.rhs.tabidx, c.rhs.colidx, rows1, nt1, 0) :
                _ek_join_keys(segs, segvec, c.rhs.tabidx, c.rhs.colidx, rows2, nt2, nt1)
            nl = length(lkeys.null)
            nr = length(rkeys.null)
            lov = _ek_shellsort!(collect(1:nl), (a, b) -> _ek_keycmp(lkeys, a, lkeys, b) != _EK_GT)
            rov = _ek_shellsort!(collect(1:nr), (a, b) -> _ek_keycmp(rkeys, a, rkeys, b) != _EK_GT)
            emit(li, ri) = lfirst ? candidate(li, ri, cnstr) : candidate(ri, li, cnstr)
            lptr = 1
            lcur = 1
            rptr = 1
            done = false
            if case == 1
                while !done
                    li = lov[lcur]
                    ri = rov[rptr]
                    rel = _ek_keycmp(lkeys, li, rkeys, ri)
                    if rel == _EK_LT
                        if lptr < nl
                            lptr += 1
                            lcur = lptr
                        else
                            done = true
                        end
                    elseif rel == _EK_EQ
                        emit(li, ri)
                        if lcur < nl
                            lcur += 1
                        elseif lcur == nl && rptr < nr
                            rptr += 1
                            lcur = lptr
                        else
                            done = true
                        end
                    else
                        if rptr < nr
                            rptr += 1
                            lcur = lptr
                        else
                            done = true
                        end
                    end
                end
            else
                lsmall = c.opcode == _EK_LT || c.opcode == _EK_LE
                while !done
                    li = lov[lptr]
                    ri = rov[rptr]
                    if _ek_relop(c.opcode, _ek_keycmp(lkeys, li, rkeys, ri))
                        emit(li, ri)
                        if lsmall
                            if rptr < nr
                                rptr += 1
                            elseif lptr < nl
                                lptr += 1
                                rptr = 1
                            else
                                done = true
                            end
                        else
                            if lptr < nl
                                lptr += 1
                            elseif rptr < nr
                                rptr += 1
                                lptr = 1
                            else
                                done = true
                            end
                        end
                    elseif lsmall
                        rptr < nr ? (rptr += 1) : (done = true)
                    else
                        lptr < nl ? (lptr += 1) : (done = true)
                    end
                end
            end
        end
        push!(out.segvecs, segvec)
        push!(out.rows, res)
    end
    return _ek_squeeze!(out)
end

# ---------------------------------------------------------------------------------------
# Union of the conjunctions
# ---------------------------------------------------------------------------------------

"""
    _ek_squeeze_compat!(rs, counts) -> rs

`zzekjsqz` as CSPICE does it after `zzekweed` deleted row vectors: `counts[k]` is the
number of row vectors of segment vector `k` before the deletion, and `rs.rows[k]` holds
the ones that are left. CSPICE removes a segment vector that has lost all its row
vectors, but keeps its (now empty) list of row vectors and removes the last list in its
place; so each later segment vector gets the row vectors of the one before it (defect
(c) of the file header).
"""
function _ek_squeeze_compat!(rs::_EKRowSet, counts::Vector{Int})
    keep = findall(k -> counts[k] > 0 && !isempty(rs.rows[k]), eachindex(counts))
    length(keep) == length(counts) && return rs
    lists = rs.rows[findall(>(0), counts)]
    kept_sv = rs.segvecs[keep]
    resize!(rs.segvecs, length(keep)) .= kept_sv
    resize!(rs.rows, length(keep)) .= lists[1:length(keep)]
    return rs
end

"""
    _ek_weed!(sets, compat) -> sets

Port of `zzekweed`: remove from each row set of `sets` the row vectors that an earlier
row set holds, then the segment vectors and the row sets that are empty. Two row
vectors are the same if their record pointers agree and their segment vectors name the
same segments; with `compat` the second condition is that the segment vectors have the
same position in their row sets, which is what CSPICE tests (defect (a) of the file
header), and the empty segment vectors are removed as CSPICE does it (defect (c)).
"""
function _ek_weed!(sets::Vector{_EKRowSet}, compat::Bool)
    isempty(sets) && return sets
    nt = sets[1].nt
    seen = Set{Tuple{Vector{Int}, Vector{Int32}}}()
    for (c, rs) in enumerate(sets)
        # invariant: `seen` holds the row vectors of the row sets before c
        newkeys = Tuple{Vector{Int}, Vector{Int32}}[]
        counts = Int[length(r) ÷ nt for r in rs.rows]
        for (pos, sv) in enumerate(rs.segvecs)
            rows = rs.rows[pos]
            kept = Int32[]
            tag = compat ? [pos] : sv
            for i in 1:counts[pos]
                rv = rows[((i - 1) * nt + 1):(i * nt)]
                if c > 1 && (tag, rv) in seen
                    continue
                end
                append!(kept, rv)
                push!(newkeys, (tag, rv))
            end
            rs.rows[pos] = kept
        end
        union!(seen, newkeys)
        compat ? _ek_squeeze_compat!(rs, counts) : _ek_squeeze!(rs)
    end
    filter!(rs -> _ek_nrows(rs) > 0, sets)
    return sets
end

# ---------------------------------------------------------------------------------------
# The search
# ---------------------------------------------------------------------------------------

"""
    EKResult

The result of one [`ekfind`](@ref): the rows that satisfy the query, in the order of
the ORDER BY clause, and the columns of the SELECT clause. It takes the place of the
global result of the last query in CSPICE.

Read it with [`ekgc`](@ref), [`ekgd`](@ref), [`ekgi`](@ref), [`eknelt`](@ref) and
[`ekcolumn`](@ref), with rows and SELECT columns that count from 1.

- `nmrows`: the number of matching rows; `length(result)` gives the same.
- `seltype`: the data type of each SELECT column, `:CHR`, `:DP`, `:INT` or `:TIME`.
- `selname`: the name of each SELECT column, as `TABLE.COLUMN` in upper case.

The other fields are internal. A result stays valid after the EK files are unloaded,
and many tasks can read one result at the same time. It is never modified.

# Examples
```julia
result, failed, errmsg = ekfind("SELECT NAME, RADIUS FROM MOONS ORDER BY NAME")
failed && throw(ArgumentError(errmsg))
for row in 1:length(result)
    name, _, _ = ekgc(result, 1, row)
    radius, isnull, _ = ekgd(result, 2, row)
    println(name, " ", isnull ? "unknown" : radius)
end
```
"""
struct EKResult
    nmrows::Int
    ntab::Int
    segs::Vector{_EKQSeg}
    rowseg::Matrix{Int32}
    rowptr::Matrix{Int32}
    seltab::Vector{Int}
    selcol::Vector{Int}
    seltype::Vector{Symbol}
    selname::Vector{String}
end

Base.length(r::EKResult) = r.nmrows
Base.show(io::IO, r::EKResult) = print(
    io, "EKResult(", r.nmrows, " rows, columns ",
    join(r.selname, ", "), ")"
)

_ek_empty_result() = EKResult(
    0, 0, _EKQSeg[], Matrix{Int32}(undef, 0, 0),
    Matrix{Int32}(undef, 0, 0), Int[], Int[], Symbol[], String[]
)

"""
    _ek_exact_constraint(c) -> _EKQConstraint

`c` as the public search compares it: with the correctly rounded value `c.dexact` in
place of `c.dval`, and an integer value outside the 32-bit range as a double. CSPICE
compares such a value through `i_dnnt`, which gives the smallest 32-bit integer (defect
(d) of the file header).
"""
function _ek_exact_constraint(c::_EKQConstraint)
    big = c.dtype == _EK_INT && !(-2147483648.0 <= c.dexact <= 2147483647.0)
    (c.dval === c.dexact && !big) && return c
    return _EKQConstraint(
        c.cnstyp, c.lhs, c.opcode, c.rhs, big ? _EK_DP : c.dtype, c.cval, c.dexact,
        c.ival, c.dexact
    )
end

"""
    _ek_search(ks, q; compat=false) -> EKResult

Port of `eksrch`: run the encoded query `q` (the result of `_ek_encode_query` for the
tables of `ks`) on the snapshot `ks`. Pure function of its arguments. `compat = true`
reproduces four defects of CSPICE (see the file header); it exists for the tests.
"""
function _ek_search(ks::KernelSet, q; compat::Bool = false)
    ntab = length(q.tables)
    segs = _EKQSeg[]
    tabsegs = Vector{Vector{Int}}(undef, ntab)
    first_of = Dict{Int, Vector{Int}}()
    for t in 1:ntab
        tabsegs[t] = get!(first_of, q.tabcat[t]) do
            table = ks.ektables[q.tabcat[t]]
            ids = Int[]
            for (k, (fno, segno)) in enumerate(table.segs)
                file = ks.ekfiles[fno]
                push!(segs, _EKQSeg(file, file.segments[segno], table.cols[k]))
                push!(ids, length(segs))
            end
            ids
        end
    end
    sets = _EKRowSet[]
    nconj = length(q.conjsizes)
    cjend = 0
    for conj in 1:max(1, nconj)
        cjsize = nconj == 0 ? 0 : q.conjsizes[conj]
        cons = q.constraints[(cjend + 1):(cjend + cjsize)]
        # The public search compares with the correctly rounded number of the query.
        compat || (cons = map(_ek_exact_constraint, cons))
        cjend += cjsize
        res = _ek_table_rows(segs, tabsegs[1], 1, cons, compat)
        for t in 2:ntab
            # invariant: `res` holds the combinations of tables 1 to t-1
            rst = _ek_table_rows(segs, tabsegs[t], t, cons, compat)
            active = Bool[
                c.cnstyp == 1 && 1 <= c.lhs.tabidx <= t && 1 <= c.rhs.tabidx <= t &&
                    c.lhs.tabidx != c.rhs.tabidx && (c.lhs.tabidx == t || c.rhs.tabidx == t)
                    for c in cons
            ]
            res = _ek_join(segs, res, rst, cons, active)
        end
        push!(sets, res)
        # CSPICE removes the duplicates after each conjunction; the result differs from
        # one pass at the end only through its defects.
        compat && _ek_weed!(sets, true)
    end
    compat || _ek_weed!(sets, false)
    n = sum(_ek_nrows, sets; init = 0)
    rowseg = Matrix{Int32}(undef, ntab, n)
    rowptr = Matrix{Int32}(undef, ntab, n)
    k = 0
    for rs in sets, (pos, sv) in enumerate(rs.segvecs)
        rows = rs.rows[pos]
        for i in 1:(length(rows) ÷ ntab)
            k += 1
            for t in 1:ntab
                rowseg[t, k] = sv[t]
                rowptr[t, k] = rows[(i - 1) * ntab + t]
            end
        end
    end
    if !isempty(q.orderby) && n > 0
        ord = _ek_order(segs, rowseg, rowptr, q.orderby)
        rowseg = rowseg[:, ord]
        rowptr = rowptr[:, ord]
    end
    seltab = Int[c.tabidx for c in q.select]
    selcol = Int[c.colidx for c in q.select]
    seltype = Symbol[ks.ektables[q.tabcat[c.tabidx]].attrs[c.colidx].dtype for c in q.select]
    selname = String[string(q.tables[c.tabidx], ".", c.name) for c in q.select]
    return EKResult(n, ntab, segs, rowseg, rowptr, seltab, selcol, seltype, selname)
end

"""
    _ek_order(segs, rowseg, rowptr, orderby) -> Vector{Int}

The order vector of `zzekjsrt` for the rows of a search: position `i` of the sorted
result is row `ord[i]`. `orderby` holds the column and the sense code (0 ascending,
1 descending) of each key. A null entry is less than every value. The Shell sort is
that of CSPICE, so rows with equal keys end in the same order.
"""
function _ek_order(
        segs::Vector{_EKQSeg}, rowseg::Matrix{Int32}, rowptr::Matrix{Int32},
        orderby::Vector
    )
    n = size(rowseg, 2)
    keys = _EKKeys[]
    senses = Int[]
    for (c, sense) in orderby
        t = c.tabidx
        k = _EKKeys(Int(segs[rowseg[t, 1]].cols[c.colidx].dtype))
        for i in 1:n
            qs = segs[rowseg[t, i]]
            _ek_pushkey!(k, qs.file, qs.cols[c.colidx], rowptr[t, i])
        end
        push!(keys, k)
        push!(senses, sense)
    end
    function jle(a::Int, b::Int)
        for (k, sense) in zip(keys, senses)   # invariant: the keys before this one are equal
            rel = _ek_keycmp(k, a, k, b)
            rel == _EK_EQ && continue
            sense == 1 && (rel = rel == _EK_LT ? _EK_GT : _EK_LT)
            return rel != _EK_GT
        end
        return true
    end
    return _ek_shellsort!(collect(1:n), jle)
end
