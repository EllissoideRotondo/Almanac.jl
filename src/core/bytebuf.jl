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

# The part of the `AbstractString` interface that the constructor of `String63` reads,
# in InlineStrings 1.x and 2.x: the length and the bytes. It is faster than
# `InlineStrings.addcodeunit` byte by byte, which shifts all 64 bytes.
Base.ncodeunits(b::_ByteBuf) = b.len
Base.codeunit(::_ByteBuf) = UInt8
Base.codeunit(b::_ByteBuf, i::Int) = _bget(b, i)

"The bytes of `b` as a `String63`; they fit, as `b` holds at most 63."
_inline(b::_ByteBuf) = String63(b)

# The buffer operations. A writer that grows a buffer with `_bresize` writes every new byte.
_blen(b::Vector{UInt8}) = length(b)
_blen(b::_ByteBuf) = b.len
@inline _bget(b::Vector{UInt8}, i::Int) = b[i]
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

# `_StackVec`, a vector of fixed length for the state of a parser. It is an `MVector`
# whose index check throws without a reference to the vector: when every function
# that takes it is inlined, the compiler keeps it on the stack, which a `BoundsError`
# that holds the vector would prevent. The check is made with `@inbounds` too.
struct _StackVec{N, T}
    data::MVector{N, T}
end

"A `_StackVec` of `N` copies of `x`."
@inline function _StackVec{N, T}(x::T) where {N, T}
    v = _StackVec{N, T}(MVector{N, T}(undef))
    # Loop invariant: elements `1:i-1` are `x`.
    for i in 1:N
        @inbounds v.data[i] = x
    end
    return v
end

"A `_StackVec` whose elements are undefined."
@inline _StackVec{N, T}(::UndefInitializer) where {N, T} = _StackVec{N, T}(MVector{N, T}(undef))

Base.length(::_StackVec{N}) where {N} = N
@inline function Base.getindex(v::_StackVec{N}, i::Int) where {N}
    1 <= i <= N || _stackvec_index_error(i, N)
    return @inbounds v.data[i]
end
@inline function Base.setindex!(v::_StackVec{N}, x, i::Int) where {N}
    1 <= i <= N || _stackvec_index_error(i, N)
    @inbounds v.data[i] = x
    return v
end
@noinline _stackvec_index_error(i::Int, N::Int) = throw(BoundsError(1:N, i))

"The first `n` bytes of `v` as a `String`, the one allocation."
@inline function _string_of(v::_StackVec{N, UInt8}, n::Int) where {N}
    0 <= n <= N || _stackvec_index_error(n, N)
    m = v.data
    return GC.@preserve m unsafe_string(pointer(m), n)
end

"""
    _StackBuf{N}

A byte buffer of at most `N` bytes with the interface of the writers (`_blen`,
`_bget`, `_bset`, `_bresize`, `_bpush`, `_puts`), changed in place, and with `length`
and indexing, which read its first `_blen` bytes. It stays on the stack where every
function that takes it is inlined.
"""
struct _StackBuf{N}
    data::_StackVec{N, UInt8}
    len::_StackVec{1, Int}
end

"An empty `_StackBuf` of `N` bytes. Its bytes are undefined: the writers write a byte before they read it."
@inline _StackBuf{N}() where {N} = _StackBuf{N}(_StackVec{N, UInt8}(undef), _StackVec{1, Int}(0))

@inline _blen(b::_StackBuf) = b.len[1]
@inline _bget(b::_StackBuf, i::Int) = b.data[i]
@inline _bset(b::_StackBuf, i::Int, c::UInt8) = (b.data[i] = c; b)
@inline function _bresize(b::_StackBuf{N}, n::Int) where {N}
    0 <= n <= N || _stackvec_index_error(n, N)
    b.len[1] = n
    return b
end
@inline function _bpush(b::_StackBuf, c::UInt8)
    n = _blen(b) + 1
    _bresize(b, n)
    b.data[n] = c
    return b
end
"Append the bytes of `s`, a string or anything with `length` and indexing."
@inline function _puts(b::_StackBuf, s)
    u = s isa AbstractString ? codeunits(s) : s
    for k in 1:length(u)
        _bpush(b, u[k])
    end
    return b
end
Base.length(b::_StackBuf) = _blen(b)
@inline function Base.getindex(b::_StackBuf, i::Int)
    1 <= i <= _blen(b) || _stackvec_index_error(i, _blen(b))
    return b.data[i]
end

"The bytes `from:from+len-1` of the buffer `buf` (anything with indexing), as a vector with `length` and indexing that holds no copy."
struct _Span{B}
    buf::B
    from::Int
    len::Int
end
Base.length(s::_Span) = s.len
@inline function Base.getindex(s::_Span, k::Int)
    1 <= k <= s.len || _stackvec_index_error(k, s.len)
    return s.buf[s.from + k - 1]
end
