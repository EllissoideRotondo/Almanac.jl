# The byte writers of the number and time formats, and `_ByteBuf`, the buffer in which
# `et2utc` and `etcal` build their result without touching the heap. A writer takes a
# buffer and returns it: a `Vector{UInt8}` (the buffer of `timout`) is changed in place,
# a `_ByteBuf` is immutable and comes back as a new value. `_putn`, `_puts` and
# `_vsbstr` take either one; the others take a `_ByteBuf` only. `_inline` turns a
# `_ByteBuf` into the `String63` (InlineStrings.jl) that the calls return.

# At most 63 bytes, the capacity of a `String63`. Byte `i` is byte `(i - 1) % 8` of word
# `(i - 1) ÷ 8 + 1`; the bytes past `len` are undefined, and `_inline` never reads them.
struct _ByteBuf <: AbstractString      # only so that the constructor of `String63` reads it
    words::NTuple{8, UInt64}
    len::Int
    function _ByteBuf(words::NTuple{8, UInt64}, len::Int)
        0 <= len <= 63 || _bytebuf_overflow(len)
        return new(words, len)
    end
end

@noinline _bytebuf_overflow(len) = throw(
    SpiceError(:buffertoosmall, "the result has $len bytes; at most 63 fit in a String63")
)

_ByteBuf() = _ByteBuf(ntuple(_ -> zero(UInt64), Val(8)), 0)

# The interface of an `AbstractString` that the constructor of `String63` uses. It is
# faster than `InlineStrings.addcodeunit` byte by byte, which shifts all 64 bytes.
Base.ncodeunits(b::_ByteBuf) = b.len
Base.codeunit(::_ByteBuf) = UInt8
Base.codeunit(b::_ByteBuf, i::Int) = _bget(b, i)
Base.isvalid(b::_ByteBuf, i::Int) = 1 <= i <= b.len
Base.iterate(b::_ByteBuf, i::Int = 1) = i > b.len ? nothing : (Char(_bget(b, i)), i + 1)

"The bytes of `b` as a `String63`; they fit, as `b` holds at most 63."
_inline(b::_ByteBuf) = String63(b)

# The buffer operations. A writer that grows a buffer with `_bresize` writes every new byte.
_blen(b::Vector{UInt8}) = length(b)
_blen(b::_ByteBuf) = b.len
@inline _bget(b::_ByteBuf, i::Int) = (b.words[((i - 1) >> 3) + 1] >> (((i - 1) & 7) << 3)) % UInt8
@inline _bset(b::Vector{UInt8}, i::Int, c::UInt8) = (b[i] = c; b)
@inline function _bset(b::_ByteBuf, i::Int, c::UInt8)
    k = ((i - 1) >> 3) + 1
    sh = ((i - 1) & 7) << 3
    w = b.words
    words = ntuple(j -> j == k ? (w[j] & ~(UInt64(0xff) << sh)) | (UInt64(c) << sh) : w[j], Val(8))
    return _ByteBuf(words, b.len)
end
@inline _bresize(b::Vector{UInt8}, n::Int) = resize!(b, n)
@inline _bresize(b::_ByteBuf, n::Int) = _ByteBuf(b.words, n)
@inline _bpush(b::Vector{UInt8}, c::UInt8) = push!(b, c)
@inline _bpush(b::_ByteBuf, c::UInt8) = _bset(_ByteBuf(b.words, b.len + 1), b.len + 1, c)

"Copy the `n` bytes at `src` to `dst` in the same buffer; the ranges may overlap."
function _bmove(b::_ByteBuf, dst::Int, src::Int, n::Int)
    if dst < src
        for k in 0:(n - 1)
            b = _bset(b, dst + k, _bget(b, src + k))
        end
    elseif dst > src
        for k in (n - 1):-1:0
            b = _bset(b, dst + k, _bget(b, src + k))
        end
    end
    return b
end

"Delete the bytes `from:to` of `b` (nothing if `to < from`)."
function _bdelete(b::_ByteBuf, from::Int, to::Int)
    to < from && return b
    n = _blen(b)
    b = _bmove(b, from, to + 1, n - to)
    return _bresize(b, n - (to - from + 1))
end

"Insert the byte `c` before position `i` of `b`."
function _binsert(b::_ByteBuf, i::Int, c::UInt8)
    n = _blen(b)
    b = _bresize(b, n + 1)
    b = _bmove(b, i + 1, i, n - i + 1)
    return _bset(b, i, c)
end

"The first position from `i` on that holds no blank; `_blen(b) + 1` if there is none."
function _skip_blanks(b::_ByteBuf, i::Int)
    while i <= _blen(b) && _bget(b, i) == UInt8(' ')
        i += 1
    end
    return i
end

"Append the byte `c` `n` times to `b` (nothing for `n <= 0`)."
@inline function _putn(b, c::UInt8, n::Int)
    for _ in 1:n
        b = _bpush(b, c)
    end
    return b
end

"Append the bytes of `s`."
@inline _puts(b::Vector{UInt8}, s::AbstractString) = append!(b, codeunits(s))
@inline function _puts(b::_ByteBuf, s::AbstractString)
    for c in codeunits(s)
        b = _bpush(b, c)
    end
    return b
end

"Append the decimal digits of `n >= 0`, zero padded on the left to `width` (`lpad(string(n), width, '0')`)."
function _putint(b::_ByteBuf, n::Integer, width::Int = 1)
    nd = ndigits(n)
    b = _putn(b, UInt8('0'), width - nd)
    base = _blen(b)
    b = _bresize(b, base + nd)
    # Loop invariant: the `k` last digits of `n` are written.
    for k in 0:(nd - 1)
        n, d = divrem(n, 10)
        b = _bset(b, base + nd - k, UInt8('0') + UInt8(d))
    end
    return b
end
