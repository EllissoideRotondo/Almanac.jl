# JSON text kernels: a text kernel written as one JSON object, in the format of
# jsonspice (https://github.com/jacobwilliams/jsonspice), which PEPPER also reads.
# The parser turns the file into the same `TextKernel` that `parse_text_lines` makes,
# so that `furnsh`, `ldpool`, `unload`, meta-kernels and the pool watchers treat it
# as any text kernel.

"""
    is_json_kernel(path) -> Bool

Contract: true if the name of `path` ends in `.json`, in any case. The content is not
inspected.
"""
is_json_kernel(path::AbstractString) = endswith(lowercase(path), ".json")

"Parser state: the text, the byte position, and the path for the messages."
mutable struct _JSONReader
    text::Vector{UInt8}
    pos::Int
    path::String
end

function _json_line(r::_JSONReader, at::Int=r.pos)
    return 1 + count(==(UInt8('\n')), view(r.text, 1:min(at, length(r.text)) - 1))
end

@noinline _json_fail(r::_JSONReader, msg, at::Int=r.pos) =
    throw(SpiceError(:badkernel, "JSON kernel: $msg (line $(_json_line(r, at)) of '$(r.path)')"))

_json_eof(r::_JSONReader) = r.pos > length(r.text)

"Skip blanks and `//` comments."
function _json_skip!(r::_JSONReader)
    t = r.text
    while r.pos <= length(t)
        c = t[r.pos]
        if c == UInt8(' ') || c == UInt8('\t') || c == UInt8('\n') || c == UInt8('\r')
            r.pos += 1
        elseif c == UInt8('/') && r.pos < length(t) && t[r.pos+1] == UInt8('/')
            while r.pos <= length(t) && t[r.pos] != UInt8('\n')
                r.pos += 1
            end
        elseif c == UInt8('/')
            _json_fail(r, "'/' that does not start a '//' comment (block comments are not allowed)")
        else
            return nothing
        end
    end
    return nothing
end

function _json_expect!(r::_JSONReader, c::UInt8, what)
    _json_skip!(r)
    (!_json_eof(r) && r.text[r.pos] == c) || _json_fail(r, "expected $what")
    r.pos += 1
    return nothing
end

function _json_hex4(r::_JSONReader)
    r.pos + 3 <= length(r.text) || _json_fail(r, "a '\\u' escape needs four hexadecimal digits")
    s = String(r.text[r.pos:r.pos+3])
    v = tryparse(UInt16, s; base=16)
    v === nothing && _json_fail(r, "'\\u$s' is not a valid escape")
    r.pos += 4
    return v
end

"Read a JSON string at `r.pos` (on the opening quote) and return it."
function _json_string!(r::_JSONReader)::String
    t = r.text
    start = r.pos
    r.pos += 1
    io = IOBuffer()
    while true
        _json_eof(r) && _json_fail(r, "a string is not closed", start)
        c = t[r.pos]
        if c == UInt8('"')
            r.pos += 1
            break
        elseif c == UInt8('\\')
            r.pos += 1
            _json_eof(r) && _json_fail(r, "a string ends in '\\'", start)
            e = t[r.pos]
            r.pos += 1
            if e == UInt8('"') || e == UInt8('\\') || e == UInt8('/')
                write(io, e)
            elseif e == UInt8('b')
                write(io, '\b')
            elseif e == UInt8('f')
                write(io, '\f')
            elseif e == UInt8('n')
                write(io, '\n')
            elseif e == UInt8('r')
                write(io, '\r')
            elseif e == UInt8('t')
                write(io, '\t')
            elseif e == UInt8('u')
                u = _json_hex4(r)
                if 0xd800 <= u <= 0xdbff       # a surrogate pair
                    (r.pos + 1 <= length(t) && t[r.pos] == UInt8('\\') && t[r.pos+1] == UInt8('u')) ||
                        _json_fail(r, "a high surrogate '\\u$(string(u; base=16))' is not followed by a low one")
                    r.pos += 2
                    lo = _json_hex4(r)
                    0xdc00 <= lo <= 0xdfff || _json_fail(r, "'\\u$(string(lo; base=16))' is not a low surrogate")
                    write(io, Char(0x10000 + ((UInt32(u) - 0xd800) << 10) + (UInt32(lo) - 0xdc00)))
                else
                    0xdc00 <= u <= 0xdfff && _json_fail(r, "a low surrogate without a high one")
                    write(io, Char(u))
                end
            else
                _json_fail(r, "'\\$(Char(e))' is not a JSON escape")
            end
        elseif c < 0x20
            _json_fail(r, "a control character in a string (write it as an escape)")
        else
            write(io, c)
            r.pos += 1
        end
    end
    s = String(take!(io))
    isvalid(s) || _json_fail(r, "a string is not valid UTF-8", start)
    return s
end

"Read a JSON number at `r.pos` (RFC 8259 grammar) and return it, correctly rounded."
function _json_number!(r::_JSONReader)::Float64
    t = r.text
    start = r.pos
    isdig(i) = i <= length(t) && UInt8('0') <= t[i] <= UInt8('9')
    i = r.pos
    i <= length(t) && t[i] == UInt8('-') && (i += 1)
    if i <= length(t) && t[i] == UInt8('0')
        i += 1
    elseif isdig(i)
        while isdig(i); i += 1; end
    else
        _json_fail(r, "a value is not a number, a string, an array, true or false", start)
    end
    if i <= length(t) && t[i] == UInt8('.')
        i += 1
        isdig(i) || _json_fail(r, "a number has no digit after its '.'", start)
        while isdig(i); i += 1; end
    end
    if i <= length(t) && (t[i] == UInt8('e') || t[i] == UInt8('E'))
        i += 1
        i <= length(t) && (t[i] == UInt8('+') || t[i] == UInt8('-')) && (i += 1)
        isdig(i) || _json_fail(r, "a number has no digit in its exponent", start)
        while isdig(i); i += 1; end
    end
    r.pos = i
    return parse(Float64, String(t[start:i-1]))
end

"""
One scalar: `(kind, value)` with kind 1 for a string, 2 for a number (`true` and `false`
are 1 and 0), and 3 for a time (a string that starts with `@`, without the `@`).
"""
function _json_scalar!(r::_JSONReader)
    _json_skip!(r)
    _json_eof(r) && _json_fail(r, "the file ends where a value is expected")
    c = r.text[r.pos]
    if c == UInt8('"')
        s = _json_string!(r)
        return startswith(s, '@') ? (3, s[2:end]) : (1, s)
    elseif c == UInt8('t') || c == UInt8('f') || c == UInt8('n')
        for (word, v) in (("true", 1.0), ("false", 0.0), ("null", nothing))
            w = codeunits(word)
            if r.pos + length(w) - 1 <= length(r.text) && view(r.text, r.pos:r.pos+length(w)-1) == w
                v === nothing && _json_fail(r, "null is not a kernel value")
                r.pos += length(w)
                return (2, v)
            end
        end
        _json_fail(r, "a value is not a number, a string, an array, true or false")
    elseif c == UInt8('{')
        _json_fail(r, "a nested object is not a kernel value")
    elseif c == UInt8('[')
        _json_fail(r, "an array inside an array is not a kernel value")
    end
    return (2, _json_number!(r))
end

"Check `name` as `parse_text_lines` checks a variable name."
function _json_check_name(r::_JSONReader, name::String, at::Int)
    isempty(name) && _json_fail(r, "a variable name is empty", at)
    ncodeunits(name) <= _TK_MAXNAME ||
        _json_fail(r, "the variable name '$name' is longer than $_TK_MAXNAME characters", at)
    all(c -> isascii(c) && isprint(c) && c != ' ', name) ||
        _json_fail(r, "the variable name '$name' has a blank or a character that is not printable", at)
    return nothing
end

"Check and convert a string value as a text kernel stores it."
function _json_string_value(r::_JSONReader, s::String, at::Int)
    v = rstrip(s, ' ')
    isempty(v) && _json_fail(r, "a string value is empty", at)
    ncodeunits(v) <= _TK_MAXSTR ||
        _json_fail(r, "a string value is longer than $_TK_MAXSTR characters, the limit of the kernel pool", at)
    return String(v)
end

"""
    parse_json_kernel(path, yearbound=1969, check=false) -> TextKernel

Contract: `path` names a readable JSON kernel; returns its assignments in file order,
with the values that the equivalent text kernel would put in the kernel pool.
`yearbound` and `check` are those of `parse_text_kernel`, for the `@` times.

The file is one JSON object (RFC 8259, plus `//` line comments and trailing commas;
a UTF-8 byte order mark is skipped). Each member is one assignment:

- a key `"NAME"` is `NAME = ...` and `"+NAME"` is `NAME += ...`;
- a number, `true` (1) or `false` (0) is a numeric value, read correctly rounded;
- a string is a string value, and a string that starts with `@` is a time, as `@` in a
  text kernel;
- an array of such values is a vector; its values are all numbers and times, or all
  strings.

Throws `SpiceError(:badkernel)` for a syntax error, `null`, a nested object or array,
an empty or mixed array, a variable name that a text kernel does not accept, or a
string value that is empty or longer than 80 characters. Nothing of the file is
returned on an error, so a load that fails changes nothing.
"""
function parse_json_kernel(path::AbstractString, yearbound::Int=1969,
        check::Bool=false)::TextKernel
    isfile(path) || throw(SpiceError(:filenotfound, "kernel file '$path' does not exist"))
    r = _JSONReader(read(path), 1, String(path))
    length(r.text) >= 3 && r.text[1:3] == [0xef, 0xbb, 0xbf] && (r.pos = 4)
    out = Assignment[]
    _json_expect!(r, UInt8('{'), "'{' at the start of a JSON kernel")
    isfirst = true
    while true
        _json_skip!(r)
        _json_eof(r) && _json_fail(r, "the object is not closed with '}'")
        if r.text[r.pos] == UInt8('}')
            r.pos += 1
            break
        end
        if !isfirst
            _json_expect!(r, UInt8(','), "',' or '}' after a member")
            _json_skip!(r)
            if !_json_eof(r) && r.text[r.pos] == UInt8('}')      # a trailing comma
                r.pos += 1
                break
            end
        end
        isfirst = false
        _json_skip!(r)
        keyat = r.pos
        (!_json_eof(r) && r.text[r.pos] == UInt8('"')) || _json_fail(r, "expected a quoted variable name")
        key = _json_string!(r)
        append = startswith(key, '+')
        name = append ? key[2:end] : key
        _json_check_name(r, name, keyat)
        _json_expect!(r, UInt8(':'), "':' after the name '$key'")
        _json_skip!(r)
        valat = r.pos
        items = Tuple{Int,Any}[]
        if !_json_eof(r) && r.text[r.pos] == UInt8('[')
            r.pos += 1
            while true
                _json_skip!(r)
                _json_eof(r) && _json_fail(r, "the array of '$key' is not closed", valat)
                if r.text[r.pos] == UInt8(']')
                    r.pos += 1
                    break
                end
                if !isempty(items)
                    _json_expect!(r, UInt8(','), "',' or ']' in the array of '$key'")
                    _json_skip!(r)
                    if !_json_eof(r) && r.text[r.pos] == UInt8(']')
                        r.pos += 1
                        break
                    end
                end
                push!(items, _json_scalar!(r))
            end
            isempty(items) && _json_fail(r, "the array of '$key' is empty", valat)
        else
            push!(items, _json_scalar!(r))
        end
        strings = count(it -> it[1] == 1, items)
        0 < strings < length(items) &&
            _json_fail(r, "the array of '$key' mixes strings with numbers or times", valat)
        if strings > 0
            push!(out, Assignment(name, append,
                String[_json_string_value(r, it[2], valat) for it in items]))
        else
            nums = Float64[]
            for (k, v) in items
                if k == 3
                    isempty(strip(v)) && _json_fail(r, "'@' is not followed by a time", valat)
                    push!(nums, _tparse(String(v), yearbound, check))
                else
                    push!(nums, v)
                end
            end
            push!(out, Assignment(name, append, nums))
        end
    end
    _json_skip!(r)
    _json_eof(r) || _json_fail(r, "text after the closing '}'")
    return TextKernel(String(path), out)
end
