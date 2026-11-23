# Text kernels: architecture sniffing, the kernel-pool grammar of `pool.c`/`zzrvar.c`,
# meta-kernel expansion as `keeper.c` does it, and the number and date parsers
# (`nparsd.c`, `tparse.c`) that fix the exact value of every numeric pool variable.

const _TK_LINLEN = 132      # LINLEN of pool.c: longer lines are cut
const _TK_MAXNAME = 32      # MAXLEN: longest variable name
const _TK_MAXSTR = 80       # MAXCHR: longest string value

# ---------------------------------------------------------------------------------
# Architecture
# ---------------------------------------------------------------------------------

"""
    kernel_arch(path) -> Symbol

Return `:SPK`, `:PCK` or `:CK` for a binary DAF of that type (the old `NAIF/DAF` too, as
`open_daf` reads it), `:DSK` or `:EK` for a DAS file of that type, and `:TEXT` for any
other file. Throws `SpiceError(:filenotfound)` for a missing file and
`SpiceError(:notsupported)` for another binary kernel type.
"""
function kernel_arch(path::AbstractString)::Symbol
    _isfile(path) || throw(SpiceError(:filenotfound, "kernel file '$path' does not exist"))
    head = open(io -> read(io, 16), path)
    n = length(head)
    idword = String(head[1:min(n, 8)])
    if startswith(idword, "DAF/")
        typ = strip(idword[5:end])
        typ == "SPK" && return :SPK
        typ == "PCK" && return :PCK
        typ == "CK" && return :CK
        throw(SpiceError(:notsupported, "DAF kernel type '$typ' of '$path' is not supported"))
    elseif idword == "NAIF/DAF"
        n >= 16 || throw(SpiceError(:badfile, "DAF file '$path' is truncated"))
        nd = reinterpret(Int32, head[9:12])[1]
        ni = reinterpret(Int32, head[13:16])[1]
        for (d, i) in ((nd, ni), (bswap(nd), bswap(ni)))
            if d == 2 && i == 6
                daf = open_daf(path)
                _unmap!(daf.data)            # only the kind is needed; see `_unmap!`
                return daf.kind
            end
            d == 2 && i == 5 && return :PCK
        end
        throw(
            SpiceError(
                :notsupported,
                "old-style DAF '$path' has a summary format that is neither SPK nor PCK"
            )
        )
    elseif startswith(idword, "DAS/")
        strip(idword[5:end]) == "DSK" && return :DSK
        strip(idword[5:end]) == "EK" && return :EK
        throw(SpiceError(:notsupported, "DAS kernel '$path' is not supported"))
    end
    return :TEXT
end

# ---------------------------------------------------------------------------------
# Numbers: port of NPARSD
# ---------------------------------------------------------------------------------

const _NPARSD_LOOKUP = (1.0, 10.0, 100.0, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8, 1.0e9, 1.0e10)

"""
    _nparsd(s) -> Float64

Return the double that `nparsd` computes for the number `s` (digits, one decimal point,
`E` or `D` exponent, signs, `PI`; blanks and commas ignored), which is not always the
correctly rounded one. Throws `SpiceError(:badnumber)` for a string that is not a number.
"""
function _nparsd(s::Union{AbstractString, AbstractVector{UInt8}})::Float64
    bytes = s isa AbstractString ? codeunits(s) : s
    bad(i) = throw(
        SpiceError(
            :badnumber,
            "'$(String(copy(bytes)))' is not a number (character $i)"
        )
    )
    intbnd = 1.0e15
    maxexp = 308.0
    minexp = -309.0
    dpmax = floatmax(Float64)
    smlbnd = dpmax / 1.0e10
    pntok = signok = expok = bpiok = doint = roundd = roundi = true
    dodec = doexp = mantsa = sigchr = epiok = zeroi = false
    intval = decval = expval = ecount = 0.0
    divisr = factor = 1.0
    signdx = 1
    sgn1 = sgn2 = 1.0
    l = length(bytes)
    while l > 0 && bytes[l] == UInt8(' ')
        l -= 1
    end
    l == 0 && throw(SpiceError(:badnumber, "a blank string is not a number"))
    @inbounds for i in 1:l
        c = bytes[i]
        if _isdigit(c)
            v = Float64(c - UInt8('0'))
            bpiok = epiok = signok = false
            sigchr = true
            if doint
                mantsa = true
                if intval < intbnd
                    intval = intval * 10.0 + v
                else
                    ecount += 1.0
                    factor /= 10.0
                    if roundi
                        roundi = false
                        v > 5.0 && (intval += 1.0)
                    end
                end
            elseif dodec
                mantsa = true
                if zeroi
                    if decval < intbnd
                        decval = decval * 10.0 + v
                        ecount -= 1.0
                    elseif roundd
                        roundd = false
                        v >= 5.0 && (decval += 1.0)
                    end
                elseif divisr < intbnd
                    decval = decval * 10.0 + v
                    divisr *= 10.0
                end
            elseif doexp
                if expval + ecount > maxexp
                    throw(
                        SpiceError(
                            :badnumber,
                            "'$(String(copy(bytes)))' is too large for a double"
                        )
                    )
                elseif expval + ecount < minexp
                else
                    expval = expval * 10.0 + sgn2 * v
                end
            else
                bad(i)
            end
        elseif c == UInt8('.')
            pntok || bad(i)
            bpiok = epiok = pntok = signok = false
            dodec = true
            doint = doexp = false
            zeroi = intval == 0.0
        elseif c == UInt8('E') || c == UInt8('D') || c == UInt8('e') || c == UInt8('d')
            sigchr = true
            expok || bad(i)
            bpiok = epiok = expok = pntok = dodec = doint = false
            doexp = signok = true
            signdx = 2
        elseif c == UInt8('+') || c == UInt8('-')
            signok || bad(i)
            v = c == UInt8('-') ? -1.0 : 1.0
            signdx == 1 ? (sgn1 = v) : (sgn2 = v)
            signok = false
        elseif c == UInt8('p') || c == UInt8('P')
            sigchr = true
            bpiok || bad(i)
            doint = dodec = doexp = expok = pntok = bpiok = signok = false
            epiok = true
        elseif c == UInt8('i') || c == UInt8('I')
            epiok || bad(i)
            doint = dodec = doexp = expok = pntok = bpiok = signok = epiok = false
            mantsa = true
            intval = Float64(pi)
        elseif c == UInt8(' ') || c == UInt8(',')
        else
            bad(i)
        end
    end
    epiok && bad(l)
    value = if mantsa
        intval + decval / divisr * factor
    elseif sigchr
        1.0
    else
        bad(l)
    end
    expval += ecount
    if expval < 0.0
        while expval < -10.0
            value /= 1.0e10
            expval += 10.0
        end
        value /= _NPARSD_LOOKUP[Int(-trunc(expval)) + 1]
    elseif expval > 0.0
        while expval > 10.0
            value >= smlbnd && throw(
                SpiceError(
                    :badnumber,
                    "'$(String(copy(bytes)))' is too large for a double"
                )
            )
            value *= 1.0e10
            expval -= 10.0
        end
        f = _NPARSD_LOOKUP[Int(round(expval)) + 1]
        value < dpmax / f || throw(
            SpiceError(
                :badnumber,
                "'$(String(copy(bytes)))' is too large for a double"
            )
        )
        value *= f
    end
    return sgn1 * value
end

# ---------------------------------------------------------------------------------
# Integers: NPARSI, and the rounding of GIPOOL
# ---------------------------------------------------------------------------------

const _NPARSI_RANGE = "NPARSI: Value entered is beyond the bounds of representable integers."

"""
    _nparsi(s) -> Int

Return the integer of `s`, as `nparsi`: the number of `_nparsd` cut toward zero. Throws
`SpiceError(:badnumber)` for a string that is not a number and for a value outside the
32-bit integers.
"""
function _nparsi(s::Union{AbstractString, AbstractVector{UInt8}})::Int
    t = trunc(_nparsd(s))
    -2147483648.0 <= t <= 2147483647.0 || throw(SpiceError(:badnumber, _NPARSI_RANGE))
    return Int(t)
end

"The number of `s` (`_nparsd`), or `nothing` where `nparsd` signals an error."
function _try_nparsd(s::Union{AbstractString, AbstractVector{UInt8}})::Union{Float64, Nothing}
    return try
        _nparsd(s)
    catch err
        err isa SpiceError || rethrow()
        nothing
    end
end

"The integer of `s` (`_nparsi`), or `nothing` where `nparsi` signals an error."
function _try_nparsi(s::Union{AbstractString, AbstractVector{UInt8}})::Union{Int, Nothing}
    return try
        _nparsi(s)
    catch err
        err isa SpiceError || rethrow()
        nothing
    end
end

"""
Return `x` rounded to the nearest integer, a tie away from zero, when it is in the range
that `gipool_` accepts; `nothing` otherwise (NaN included), where `gipool_` signals
SPICE(INTOUTOFRANGE).
"""
@inline _gipool_value(x::Float64)::Union{Int32, Nothing} =
    -2147483648.0 <= x <= 2147483647.0 ? Int32(_nint(x)) : nothing

@noinline _gipool_range_error(name::AbstractString, i::Integer, x::Float64) = throw(
    SpiceError(
        :intoutofrange, "the value associated with index $i of the kernel " *
            "variable $name is outside the range of integers; the value stored was $x"
    )
)

"""
    _gipool_first(v, name; room = 1, single = false, strict = true) -> Union{Int32,Nothing}

Return the first value of the pool variable `name` whose values are `v`, as
`gipool_(name, 1, room, n, ivals, found)` gives it, or `nothing` if `v` is not numeric or
is empty. Each of the first `room` values must be in the range of `gipool_`, else
`SpiceError(:intoutofrange)`, or `nothing` with `strict = false`. With `single`, the
result is `nothing` unless `v` has one value.
"""
function _gipool_first(
        v, name::AbstractString; room::Int = 1, single::Bool = false, strict::Bool = true
    )::Union{Int32, Nothing}
    (v isa Vector{Float64} && !isempty(v)) || return nothing
    for i in 1:min(room, length(v))
        if _gipool_value(@inbounds v[i]) === nothing
            strict || return nothing
            _gipool_range_error(name, i, v[i])
        end
    end
    single && length(v) != 1 && return nothing
    return _gipool_value(@inbounds v[1])
end

"Return the numeric values of the pool variable `name`, or `nothing` if it is missing or holds strings."
function _pool_numbers(pool::AbstractDict, name::AbstractString)
    v = get(pool, name, nothing)
    return v isa Vector{Float64} ? v : nothing
end

# ---------------------------------------------------------------------------------
# Dates: the `@` forms, as TPARSE converts them
# ---------------------------------------------------------------------------------

const _MONTH_ABBR = (
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT",
    "NOV", "DEC",
)

"""
    _tparse(s, yearbound=1969, check=false) -> Float64

Return the seconds past J2000 that `tparse` gives for the kernel date `s` (without the
`@`), as `zzrvar_` reads it. `yearbound` and `check` are the settings of `tsetyr` and
`tparch` of the context at the load. Throws `SpiceError(:badtimespec)` for a string that
`tparse` rejects.
"""
function _tparse(s::AbstractString, yearbound::Int = 1969, check::Bool = false)::Float64
    @assert yearbound > 0
    d = TimeDefaults()
    try
        return tparse(s, TimeDefaults(d.calendar, d.system, d.zone, yearbound, check))
    catch e
        e isa SpiceError && e.kind in (:unparsedtime, :badtimestring) || rethrow()
        throw(SpiceError(:badtimespec, e.msg))
    end
end

# ---------------------------------------------------------------------------------
# Grammar: port of the scanner of ZZRVAR
# ---------------------------------------------------------------------------------

# Token types of ZZRVAR.
const _TOK_QUOTED = 1
const _TOK_VALUE = 2
const _TOK_LPAREN = 3
const _TOK_RPAREN = 4
const _TOK_EQUAL = 5
const _TOK_PLUSEQ = 6

"""
    _tokenize(line) -> (begs, ends, types)

Return the token boundaries (1-based, inclusive) and types of the data line `line`, as
the scanner loop of `zzrvar` finds them.
"""
function _tokenize(line::Vector{UInt8})
    begs = Int[]
    ends = Int[]
    types = Int[]
    even = true
    intokn = false
    inquot = false
    insepf = true
    @inbounds for i in eachindex(line)
        c = line[i]
        if c == UInt8(' ') || c == UInt8(',') || c == UInt8('(') || c == UInt8(')') ||
                c == UInt8('=') || c == 0x09
            if intokn || (inquot && even)
                inquot = false
                intokn = false
                insepf = true
            end
            if insepf
                if c == UInt8('=')
                    push!(begs, i)
                    push!(ends, i)
                    push!(types, _TOK_EQUAL)
                    if i > 1 && line[i - 1] == UInt8('+')
                        n = length(begs)
                        types[n] = _TOK_PLUSEQ
                        begs[n] = i - 1
                        if n >= 2 && begs[n - 1] == ends[n - 1]
                            pop!(begs)
                            pop!(ends)
                            pop!(types)
                            begs[n - 1] = i - 1
                            ends[n - 1] = i
                            types[n - 1] = _TOK_PLUSEQ
                        elseif n >= 2
                            ends[n - 1] -= 1
                        end
                    end
                elseif c == UInt8(')')
                    push!(begs, i)
                    push!(ends, i)
                    push!(types, _TOK_RPAREN)
                elseif c == UInt8('(')
                    push!(begs, i)
                    push!(ends, i)
                    push!(types, _TOK_LPAREN)
                end
            end
        elseif c == UInt8('\'')
            even = !even
            if !inquot
                insepf = false
                intokn = false
                inquot = true
                push!(begs, i)
                push!(ends, i)
                push!(types, _TOK_QUOTED)
            end
            ends[end] = i
        else
            if insepf || (inquot && even)
                inquot = false
                insepf = false
                intokn = true
                push!(begs, i)
                push!(ends, i)
                push!(types, _TOK_VALUE)
            end
            ends[end] = i
        end
    end
    return begs, ends, types
end

"Return the value of the quoted token `line[b:e]` as `zzrvar` stores it: quotes removed, `''` reduced to `'`, at most 80 characters, trailing blanks removed."
function _string_value(line::Vector{UInt8}, b::Int, e::Int)::String
    line[e] == UInt8('\'') || (e += 1)
    out = UInt8[]
    j = b + 1
    while j < e
        line[j] == UInt8('\'') && (j += 1)
        j <= length(line) || break
        length(out) < _TK_MAXSTR && push!(out, line[j])
        j += 1
    end
    while !isempty(out) && out[end] == UInt8(' ')
        pop!(out)
    end
    return String(out)
end

"""
    _clean_line(raw, tabs=true) -> Vector{UInt8}

Return the physical line `raw` as `rdkdat` hands it on: cut at 132 characters, tabs
replaced by blanks if `tabs`, blanks removed at both ends. The lines of `lmpool` keep
their tabs (`tabs = false`), as `zzrvbf` reads them.
"""
function _clean_line(raw::AbstractString, tabs::Bool = true)::Vector{UInt8}
    cu = codeunits(raw)
    n = length(cu)
    while n > 0 && (cu[n] == UInt8('\r') || cu[n] == UInt8('\n'))
        n -= 1
    end
    n = min(n, _TK_LINLEN)
    blank(c) = c == UInt8(' ') || (tabs && c == 0x09)
    hi = n
    while hi > 0 && blank(cu[hi])
        hi -= 1
    end
    lo = 1
    while lo <= hi && blank(cu[lo])
        lo += 1
    end
    # One copy of the bytes `lo:hi`, with the tabs replaced.
    b = Vector{UInt8}(undef, hi - lo + 1)
    for i in lo:hi            # invariant: b[1:i-lo] holds the cleaned bytes lo:i-1
        c = cu[i]
        b[i - lo + 1] = (tabs && c == 0x09) ? UInt8(' ') : c
    end
    return b
end

const _BEGINDATA = Vector{UInt8}(codeunits("\\begindata"))
const _BEGINTEXT = Vector{UInt8}(codeunits("\\begintext"))

# The kinds are the short messages of `zzrvar`: BADVARASSIGN, NONPRINTINGCHAR, BADVARNAME,
# TYPEMISMATCH, NUMBEREXPECTED, BADTIMESPEC.
@noinline _tk_fail(path, lineno::Int, msg, kind::Symbol = :badvarassign) =
    throw(SpiceError(kind, "$msg (line $lineno of '$path')"))

function _tk_number(path, lineno, bytes)
    try
        return _nparsd(bytes)
    catch e
        (e isa SpiceError && e.kind === :badnumber) || rethrow()
        _tk_fail(path, lineno, e.msg, :numberexpected)
    end
end

"Add the assignment of `name` with the values of `kind` (1 strings, 2 numbers) to `out`; nothing for kind 0."
function _tk_finish!(out::Vector{Assignment}, kind::Int, name, append::Bool, strs, nums)
    @assert 0 <= kind <= 2
    kind == 1 && push!(out, Assignment(name, append, strs))
    kind == 2 && push!(out, Assignment(name, append, nums))
    return nothing
end

"""
    parse_text_lines(path, lines, yearbound=1969, check=false; buffer=false) -> TextKernel

Return the assignments of the physical lines `lines` of a text kernel, in file order, with
the values that CSPICE puts in the kernel pool. `path` names the result and the messages;
`yearbound` and `check` are the settings of `tsetyr` and `tparch` for `@` dates. With
`buffer`, the lines are those of `lmpool`, read as `zzrvbf` reads them: every line is
data and a tab stays a tab. Throws `SpiceError` with the kind of `zzrvar` for a malformed
assignment.
"""
function parse_text_lines(
        path::AbstractString, lines, yearbound::Int = 1969,
        check::Bool = false; buffer::Bool = false
    )::TextKernel
    out = Assignment[]
    indata = buffer
    invec = false
    name = ""
    append = false
    nums = Float64[]
    strs = String[]
    kind = 0                       # 0 unknown, 1 string, 2 numeric
    # The kinds of the variables that the earlier assignments of these lines set: `zzrvar_`
    # stores each assignment before it reads the next one, so a `+=` takes the kind of the
    # variable and a value of the other kind is TYPEMISMATCH on its own line.
    kinds = Dict{String, Int}()
    lineno = 0
    for raw in lines
        lineno += 1
        line = _clean_line(raw, !buffer)
        if !buffer && line == _BEGINTEXT
            indata = false
            continue
        elseif !buffer && line == _BEGINDATA
            indata = true
            continue
        end
        (indata && !isempty(line)) || continue
        begs, ends, types = _tokenize(line)
        count = length(begs)
        nxt = 1
        if !invec
            count >= 3 || _tk_fail(
                path, lineno, "a kernel variable assignment must have the form " *
                    "'<name> [+]= <values>'"
            )
            # `lastpc_` of the name: the last character must be printable (33 to 126).
            0x21 <= line[ends[1]] <= 0x7e || _tk_fail(
                path, lineno, "the variable name ends in a nonprinting character " *
                    "(ASCII code $(Int(line[ends[1]])))", :nonprintingchar
            )
            ends[1] - begs[1] + 1 <= _TK_MAXNAME ||
                _tk_fail(
                path, lineno, "a variable name is longer than $_TK_MAXNAME characters",
                :badvarname
            )
            name = String(line[begs[1]:ends[1]])
            types[2] == _TOK_EQUAL || types[2] == _TOK_PLUSEQ ||
                _tk_fail(path, lineno, "the assignment operator of '$name' is neither '=' nor '+='")
            append = types[2] == _TOK_PLUSEQ
            nums = Float64[]
            strs = String[]
            kind = append ? get(kinds, name, 0) : 0
            if types[3] == _TOK_LPAREN
                nxt = 4
            else
                nxt = 3
                push!(begs, 0)
                push!(ends, 0)
                push!(types, _TOK_RPAREN)
                count += 1
            end
        end
        # Loop invariant: tokens 1..nxt-1 of this line are consumed and every value
        # among them is in `nums` or `strs`.
        while nxt <= count && types[nxt] != _TOK_RPAREN
            t, b, e = types[nxt], begs[nxt], ends[nxt]
            if kind == 0
                kind = t == _TOK_QUOTED ? 1 : t == _TOK_VALUE ? 2 :
                    _tk_fail(
                        path, lineno, "the first item after the assignment operator of '$name' " *
                        "is not a value"
                    )
            end
            if kind == 1
                t == _TOK_QUOTED || _tk_fail(
                    path, lineno, "'$name' is a string variable but a value " *
                        "is not a string", :typemismatch
                )
                b + 1 >= e && _tk_fail(
                    path, lineno, "there is a quoted string with no characters",
                    :typemismatch
                )
                push!(strs, _string_value(line, b, e))
            else
                t == _TOK_VALUE || _tk_fail(
                    path, lineno, "'$name' is a numeric variable but a value " *
                        "is not a number or a time", :typemismatch
                )
                if line[b] == UInt8('@')
                    e == b && _tk_fail(path, lineno, "'@' is not followed by a time", :badtimespec)
                    push!(nums, _tparse(String(line[(b + 1):e]), yearbound, check))
                else
                    push!(nums, _tk_number(path, lineno, view(line, b:e)))
                end
            end
            nxt += 1
        end
        # `( )`: a closing parenthesis before the first value
        nxt <= count && isempty(strs) && isempty(nums) &&
            _tk_fail(path, lineno, "the vector assigned to '$name' has no values")
        invec = nxt > count
        if !invec
            _tk_finish!(out, kind, name, append, strs, nums)
            kind == 0 || (kinds[name] = kind)
        end
    end
    invec && _tk_finish!(out, kind, name, append, strs, nums)
    return TextKernel(String(path), out)
end

"""
    parse_text_kernel(path, yearbound=1969, check=false) -> TextKernel

Return the assignments of the text kernel `path` in file order. String elements keep
their continuation markers (see `pool_string_joined`). A file whose name ends in `.json`
is read by `parse_json_kernel`.
"""
function parse_text_kernel(
        path::AbstractString, yearbound::Int = 1969,
        check::Bool = false
    )::TextKernel
    _isfile(path) || throw(SpiceError(:filenotfound, "kernel file '$path' does not exist"))
    is_json_kernel(path) && return parse_json_kernel(path, yearbound, check)
    text = read(path, String)
    return parse_text_lines(path, _physical_lines(text), yearbound, check)
end

"""
    _physical_lines(text::String) -> Vector{SubString{String}}

Return the lines of `text` split at CR-LF, LF or a lone CR, as
`eachsplit(text, r"\\r\\n|\\n|\\r")` does, without a regex search per line.
"""
function _physical_lines(text::String)
    cu = codeunits(text)
    n = length(cu)
    out = SubString{String}[]
    start = 1
    i = 1
    @inbounds while i <= n
        c = cu[i]
        if c == UInt8('\n') || c == UInt8('\r')
            push!(out, SubString(text, start, prevind(text, i)))
            i += (c == UInt8('\r') && i < n && cu[i + 1] == UInt8('\n')) ? 2 : 1
            start = i
        else
            i += 1
        end
    end
    push!(out, SubString(text, start, lastindex(text)))
    return out
end

# ---------------------------------------------------------------------------------
# Pool updates and meta-kernels
# ---------------------------------------------------------------------------------

"""
    apply_assignments!(pool, tk) -> pool

Apply the assignments of `tk` to `pool` in order: `=` replaces a variable and `+=`
appends to it. A stored vector is never mutated, because older snapshots share it.
Throws `SpiceError(:typemismatch)` when `+=` mixes numeric and string values.
"""
function apply_assignments!(pool::Dict{String, PoolValue}, tk::TextKernel)
    for a in tk.assignments
        old = a.append ? get(pool, a.name, nothing) : nothing
        if old === nothing
            pool[a.name] = a.value
        elseif old isa Vector{Float64} && a.value isa Vector{Float64}
            pool[a.name] = vcat(old, a.value)
        elseif old isa Vector{String} && a.value isa Vector{String}
            pool[a.name] = vcat(old, a.value)
        else
            throw(
                SpiceError(
                    :typemismatch, "'+=' on kernel variable '$(a.name)' in " *
                        "'$(tk.path)' mixes numeric and string values"
                )
            )
        end
    end
    return pool
end

"""
    pool_strings_joined(values, contin) -> Vector{String}

Return the logical strings of the raw elements `values` of a string variable, as `stpool`
joins them: an element that ends in the marker `contin` is joined with the next one.
"""
function pool_strings_joined(values::Vector{String}, contin::AbstractString)::Vector{String}
    out = String[]
    buf = IOBuffer()
    open_string = false
    for v in values
        part = rstrip(v, ' ')
        if endswith(part, contin)
            print(buf, chop(part; tail = length(contin)))
            open_string = true
        else
            print(buf, part)
            push!(out, String(take!(buf)))
            open_string = false
        end
    end
    open_string && push!(out, String(take!(buf)))
    return out
end

"""
    pool_string_joined(values, nth, contin) -> Union{String,Nothing}

Return the `nth` logical string of `pool_strings_joined`, or `nothing` if there
is none.
"""
function pool_string_joined(values::Vector{String}, nth::Integer, contin::AbstractString)
    all = pool_strings_joined(values, contin)
    return 1 <= nth <= length(all) ? all[nth] : nothing
end

"""
    metakernel_files(tk) -> Vector{String}

Return the files of `KERNELS_TO_LOAD` of `tk` in load order, with every `\$SYMBOL` of
`PATH_SYMBOLS` replaced by its `PATH_VALUES` entry; empty if `tk` is not a meta-kernel.
A relative path is relative to the current directory, as in CSPICE.
"""
function metakernel_files(tk::TextKernel)::Vector{String}
    pool = apply_assignments!(Dict{String, PoolValue}(), tk)
    ktl = get(pool, "KERNELS_TO_LOAD", nothing)
    ktl isa Vector{String} || return String[]
    files = pool_strings_joined(ktl, "+")
    syms = get(pool, "PATH_SYMBOLS", nothing)
    syms isa Vector{String} || return files
    vals = get(pool, "PATH_VALUES", nothing)
    vals isa Vector{String} || throw(
        SpiceError(
            :badkernel,
            "meta-kernel '$(tk.path)' has PATH_SYMBOLS but no PATH_VALUES"
        )
    )
    values = pool_strings_joined(vals, "+")
    length(values) == length(syms) || throw(
        SpiceError(
            :badkernel,
            "meta-kernel '$(tk.path)' has $(length(syms)) path symbols and " *
                "$(length(values)) path values"
        )
    )
    return [_substitute_symbols(f, syms, values) for f in files]
end

"Return `file` with each `\$SYMBOL` replaced by its value; the longest matching symbol wins, and inserted text is not scanned again."
function _substitute_symbols(file::String, syms::Vector{String}, values::Vector{String})
    out = IOBuffer()
    i = firstindex(file)
    while i <= lastindex(file)
        c = file[i]
        if c == '$'
            rest = SubString(file, nextind(file, i))
            use = 0
            for (k, s) in enumerate(syms)
                !isempty(s) && startswith(rest, s) &&
                    (use == 0 || length(s) > length(syms[use])) && (use = k)
            end
            if use > 0
                print(out, values[use])
                i = nextind(file, i, length(syms[use]) + 1)
                continue
            end
        end
        print(out, c)
        i = nextind(file, i)
    end
    return String(take!(out))
end
