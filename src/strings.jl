# ASCII string helpers of SPICELIB: upper and lower case of a byte (`ucase`, `lcase`),
# comparison without regard to case and blanks (`eqstr`) and the template match of
# `matchi`. They read the code units of a string and do not allocate, except
# `_squeeze`, which returns a new string.

"Contract: byte `c` in upper case if it is an ASCII lower-case letter, else `c`."
@inline _upper(c::UInt8) = (UInt8('a') <= c <= UInt8('z')) ? c - 0x20 : c
"Contract: byte `c` in lower case if it is an ASCII upper-case letter, else `c`."
@inline _lower(c::UInt8) = (UInt8('A') <= c <= UInt8('Z')) ? c + 0x20 : c

"`s` without blanks and in upper case, as `cmprss_` with zero blanks and `ucase_`."
_squeeze(s::AbstractString) = uppercase(replace(s, ' ' => ""))

"""
    _ljucrs4(s) -> NTuple{4,UInt8}

Contract: the first four characters of `s` after its blanks are removed, ASCII letters in
upper case and blanks after the end, as `ljucrs_(0, s, out)` writes them into an `out` of
four characters. Does not allocate.
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

"Contract: true if `a` and `b` are equal without regard to ASCII case and blanks, as `eqstr_`. Does not allocate."
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
        (i > length(ca) || j > length(cb)) && return i > length(ca) && j > length(cb)
        _upper(ca[i]) == _upper(cb[j]) || return false
        i += 1
        j += 1
    end
end

"""
    _eqci(a, b, trail=false) -> Bool

Contract: `b` is upper-case ASCII; true if `a` equals `b` when ASCII case is ignored
and, with `trail`, trailing blanks of `a` are ignored. Does not allocate.
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

Port of `matchi_` with the wild string `*` and the wild character `%`: true if `string`
matches the template `templ` without regard to ASCII case; leading and trailing blanks
of both are ignored; `*` stands for any number of characters and `%` for exactly one.
Does not allocate.
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
        tcur > tlast && return true
        left = tcur
        while tcur <= tlast && !iswstr(tcur)
            tcur += 1
        end
        right = tcur - 1
    end
end
