# ASCII string helpers of SPICELIB: upper and lower case of a byte (`ucase`, `lcase`),
# comparison without regard to case and blanks (`eqstr`) and the template match of
# `matchi`. They read the code units of a string and do not allocate, except
# `_squeeze`, which returns a new string.

@inline _upper(c::UInt8) = (UInt8('a') <= c <= UInt8('z')) ? c - 0x20 : c
@inline _lower(c::UInt8) = (UInt8('A') <= c <= UInt8('Z')) ? c + 0x20 : c
@inline _isdigit(c::UInt8) = UInt8('0') <= c <= UInt8('9')

# `cmprss_` with zero blanks, then `ucase_`.
_squeeze(s::AbstractString) = uppercase(replace(s, ' ' => ""))

"""
    _cut_bytes(s, n) -> SubString

Return `s` cut to its first `n` bytes, as a Fortran string of length `n` holds it; a
character of several bytes that does not fit whole is left out.
"""
_cut_bytes(s::AbstractString, n::Int) =
    ncodeunits(s) <= n ? SubString(s) : SubString(s, 1, prevind(s, thisind(s, n + 1)))

"""
    _repmc(msg, value) -> String

Return `msg` with its first `#` replaced by `value`, as `repmc_` and `repmi_`: a string
loses its leading and trailing blanks, and a blank string becomes one blank.
"""
function _repmc(msg::String, value::Union{AbstractString, Integer})::String
    r = findfirst(==(UInt8('#')), codeunits(msg))
    r === nothing && return msg
    v = value isa Integer ? string(value) : String(strip(value, ' '))
    isempty(v) && (v = " ")
    c = codeunits(msg)
    return string(String(c[1:(r - 1)]), v, String(c[(r + 1):end]))
end

"""
    _ljucase(s) -> String

Return `s` without leading and trailing blanks and with ASCII letters in upper case, as
`ljust_` and `ucase_` give it before a comparison with a constant.
"""
function _ljucase(s::AbstractString)
    t = strip(s, ' ')
    any(c -> 'a' <= c <= 'z', t) || return String(t)
    return map(c -> 'a' <= c <= 'z' ? c - 32 : c, t)
end

"""
    _ljucrs4(s) -> NTuple{4,UInt8}

Return the first four non-blank characters of `s` in upper case, padded with blanks, as
`ljucrs_(0, s, out)` writes them into an `out` of four characters.
"""
function _ljucrs4(s::AbstractString)
    b = (UInt8(' '), UInt8(' '), UInt8(' '), UInt8(' '))
    j = 0
    for c in codeunits(s)                # invariant: b[1:j] holds the first j non-blanks
        c == UInt8(' ') && continue
        j += 1
        b = Base.setindex(b, _upper(c), j)
        j == 4 && break
    end
    return b
end

"Return true if `a` and `b` are equal without regard to ASCII case and blanks, as `eqstr_`."
function _eqstr(a::AbstractString, b::AbstractString)
    ca = codeunits(a)
    cb = codeunits(b)
    i = 1
    j = 1
    while true                           # invariant: the non-blank characters before i and j agree
        while i <= length(ca) && ca[i] == UInt8(' ')
            i += 1
        end
        while j <= length(cb) && cb[j] == UInt8(' ')
            j += 1
        end
        (i > length(ca) || j > length(cb)) && break
        _upper(ca[i]) == _upper(cb[j]) || return false
        i += 1
        j += 1
    end
    return i > length(ca) && j > length(cb)
end

"""
    _eqci(a, b, trail=false) -> Bool

Return true if `a` equals the upper-case ASCII string `b` without regard to ASCII case
and, with `trail`, to trailing blanks of `a`.
"""
function _eqci(a::AbstractString, b::String, trail::Bool = false)::Bool
    ca = codeunits(a)
    cb = codeunits(b)
    n = length(cb)
    (trail ? length(ca) >= n : length(ca) == n) || return false
    @inbounds for i in 1:n
        _upper(ca[i]) == cb[i] || return false
    end
    @inbounds for i in (n + 1):length(ca)
        ca[i] == UInt8(' ') || return false
    end
    return true
end

"""
    _matchi(string, templ) -> Bool

Return true if `string` matches the template `templ` of `matchi_` without regard to ASCII
case and to leading and trailing blanks. `*` stands for any number of characters and `%`
for exactly one.
"""
function _matchi(string::AbstractString, templ::AbstractString)
    s = codeunits(string)
    t = codeunits(templ)
    wstr = UInt8('*')
    wchr = UInt8('%')
    blank = UInt8(' ')
    sfirst = something(findfirst(!=(blank), s), 0)
    slast = something(findlast(!=(blank), s), 0)
    tfirst = something(findfirst(!=(blank), t), 0)
    tlast = something(findlast(!=(blank), t), 0)
    scur = max(1, sfirst)
    tcur = tfirst
    tlast == 0 && return slast == 0
    samei(i, j) = 1 <= i <= length(s) && 1 <= j <= length(t) && _upper(s[i]) == _upper(t[j])
    iswstr(j) = 1 <= j <= length(t) && t[j] == wstr
    iswchr(j) = 1 <= j <= length(t) && t[j] == wchr
    # The part of the template before the first wild string must match the start.
    while tcur <= tlast && scur <= slast && !iswstr(tcur)
        (_upper(t[tcur]) != _upper(s[scur]) && t[tcur] != wchr) && return false
        tcur += 1
        scur += 1
    end
    (tcur > tlast && scur > slast) && return true
    (tcur > tlast && scur <= slast) && return false
    if tcur <= tlast && scur > slast
        return all(j -> t[j] == wstr, tcur:tlast)
    end
    while tcur <= tlast && iswstr(tcur)
        tcur += 1
    end
    tcur > tlast && return true
    left = tcur
    while tcur <= tlast && !iswstr(tcur)
        tcur += 1
    end
    right = tcur - 1
    while true
        # invariant: the template before `left` matches the string before `scur`, and
        # `left:right` (with `left <= tlast`) is the next part of the template without a
        # wild string; each pass returns or moves to the next part
        slast - scur < right - left && return false
        if right == tlast
            i = slast
            j = tlast
            while j >= left
                (iswchr(j) || samei(i, j)) || return false
                j -= 1
                i -= 1
            end
            return true
        end
        nosubm = true
        while nosubm
            i = scur
            j = left
            while j <= right && (samei(i, j) || iswchr(j))
                i += 1
                j += 1
            end
            if j > right
                scur = i
                nosubm = false
            else
                scur += 1
                slast - scur < right - left && return false
            end
        end
        while tcur <= tlast && iswstr(tcur)
            tcur += 1
        end
        tcur > tlast && break
        left = tcur
        while tcur <= tlast && !iswstr(tcur)
            tcur += 1
        end
        right = tcur - 1
    end
    return true
end
