# Layouts of the SPK and binary PCK segments whose evaluator reads the places of its
# parts from meta data at the end of the segment: SPK type 10 and type 14 (generic
# segments, `sgmeta`), SPK type 19 (mini-segments) and binary PCK type 3 (the layout of
# SPK type 14).
#
# CSPICE checks these meta data at each evaluation (`spkr14`, `spkr10`, `spkr19`, with
# `sgmeta`), so a file with a bad segment loads and only the evaluation of that segment
# fails. `_daf_parse` calls `_segment_layout` once per segment: when the checks pass, it
# stores the layout in `Segment.layout` and the evaluator reads the layout from there;
# when they fail, the layout stays zero and the evaluator makes the checks again at each
# call, so that it throws the same error at the same time as before. A `Segment` that
# another constructor builds (a descriptor of `spkpvn`, a test) has a zero layout too.

@noinline _badsegment(seg::Segment, msg::String) = throw(
    SpiceError(
        :baddaf,
        "SPK type $(seg.dtype) segment at $(seg.baddr):$(seg.eaddr): $msg"
    )
)

"""
    _segerror(kind, seg, msg) -> Union{}

Throws `SpiceError(kind)` for a value of the segment `seg` that the evaluator of CSPICE
refuses with this kind.
"""
@noinline _segerror(kind::Symbol, seg::Segment, msg::String) = throw(
    SpiceError(kind, "SPK type $(seg.dtype) segment at $(seg.baddr):$(seg.eaddr): $msg")
)

"""
    _seg_count(seg, x::Float64, lo::Int, hi::Int, what) -> Int

Nearest whole number to `x` (as `i_dnnt` of f2c); throws `SpiceError(:baddaf)` that
names `what` unless `x` is finite and the result is in `lo:hi`.
"""
@inline function _seg_count(seg::Segment, x::Float64, lo::Int, hi::Int, what::String)
    (isfinite(x) && lo - 0.5 < x < hi + 0.5) || _badsegment(seg, _count_msg(what, x, lo, hi))
    n = _nint(x)
    (lo <= n <= hi) || _badsegment(seg, _count_msg(what, x, lo, hi))
    return n
end
@noinline _count_msg(what::String, x::Float64, lo::Int, hi::Int) = "$what = $x is not in $lo:$hi"

"Largest whole number that a segment count may be (the 32-bit address range of a DAF)."
const _MAXCOUNT = Int(typemax(Int32))

"""
    _sg_meta_base(nmeta, e) -> (m, haspktoff)

Addresses of the meta data of a generic segment that ends at the address `e` and
reports `nmeta >= 15` items, as `sgmeta_` reads them: item k is at `m + k`. A segment
of 15 items was written as 16 doubles. Only a segment of 17 or more items has a packet
offset (item 16): `sgmeta_` clears item 16 of a segment of 15 or 16 items. Items after
the 17th are not used.
"""
@inline _sg_meta_base(nmeta::Int, e::Int) = e - (nmeta == 15 ? 16 : nmeta), nmeta >= 17

"True if `seg.layout` holds the layout that `_segment_layout` found valid at load."
@inline _has_layout(seg::Segment) = @inbounds(seg.layout[1]) != Int32(0)

# --- SPK type 14 and binary PCK type 3 ------------------------------------------------------

"""
    _spk14_check(data, seg) -> NamedTuple

The parts of the SPK type 14 or binary PCK type 3 segment `seg` (the same layout) from
its meta data, as offsets from the address before the segment: `conbas`, `refbas`,
`pktbas`, the counts `nref` and `npkt`, the packet size `pktsz`, the packet offset
`pktoff`, `stride = pktsz + pktoff`, and the coefficient count `ncoef` of a series.
Throws `SpiceError(:notsupported)` for a layout that `spk14b` does not write and
`SpiceError(:baddaf)` for one that does not fit.
"""
@inline function _spk14_check(data::Vector{Float64}, seg::Segment)
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    nmeta = _seg_count(seg, @inbounds(data[e]), 15, len, "meta data count")
    mb, hasoff = _sg_meta_base(nmeta, e)        # meta item i is data[mb+i]
    mb >= b - 1 || _badsegment(seg, "the meta data do not fit the segment")
    @inbounds begin
        # A base of `len` would put the first constant after the segment.
        conbas = _seg_count(seg, data[mb + 1], 0, len - 1, "constants base")
        ncon = _seg_count(seg, data[mb + 2], 1, len, "constants count")
        rdrtyp = _seg_count(seg, data[mb + 5], 0, 4, "reference directory type")
        refbas = _seg_count(seg, data[mb + 6], 0, len, "reference base")
        nref = _seg_count(seg, data[mb + 7], 1, len, "reference count")
        pdrtyp = _seg_count(seg, data[mb + 10], 0, 1, "packet directory type")
        pktbas = _seg_count(seg, data[mb + 11], 0, len, "packet base")
        npkt = _seg_count(seg, data[mb + 12], 1, len, "packet count")
        pktsz = _seg_count(seg, data[mb + 15], 8, len, "packet size")
        pktoff = hasoff ? _seg_count(seg, data[mb + 16], 0, len, "packet offset") : 0
    end
    (rdrtyp == 3 && pdrtyp == 0) || throw(
        SpiceError(
            :notsupported,
            "$(seg.dtype == 14 ? "SPK type 14" : "binary PCK type 3") segment with " *
                "reference directory type $rdrtyp and packet " *
                "directory type $pdrtyp is not supported"
        )
    )
    stride = pktsz + pktoff
    (nref == npkt && refbas + nref <= len && Int128(pktbas) + Int128(npkt) * stride <= len) ||
        _badsegment(seg, "the meta data do not agree with the length")
    ncoef = _seg_count(seg, @inbounds(data[b - 1 + conbas + 1]), 1, len, "coefficient count")
    6 * ncoef + 2 == pktsz || _badsegment(seg, "the packet size is not 2 + 6*NCOEF")
    return (; conbas, refbas, nref, pktbas, npkt, pktsz, pktoff, stride, ncoef)
end

"`Segment.layout` of a type 14 or PCK type 3 segment with the parts `lay` of `_spk14_check`."
@inline _spk14_encode(lay) = (
    Int32(lay.npkt), Int32(lay.refbas), Int32(lay.pktbas), Int32(lay.stride),
    Int32(lay.pktoff), Int32(lay.ncoef),
)

"""
    _spk14_layout(data, seg) -> NamedTuple

The parts of the SPK type 14 or binary PCK type 3 segment `seg` that the evaluator
needs: `refbas`, `nref`, `pktbas`, `npkt`, `pktoff`, `stride` and `ncoef`, with the
meaning of [`_spk14_check`](@ref). They come from `seg.layout`, or from the checks of
`_spk14_check` (and their errors) when the segment has no layout.
"""
@inline function _spk14_layout(data::Vector{Float64}, seg::Segment)
    l = _has_layout(seg) ? seg.layout : _spk14_encode(_spk14_check(data, seg))
    @inbounds n, refbas, pktbas, stride, pktoff, ncoef = Int.(l)
    return (; refbas, nref = n, pktbas, npkt = n, pktoff, stride, ncoef)
end

# --- SPK type 10 ------------------------------------------------------------------------------

@noinline _t10_error(seg::Segment, msg::String) = throw(
    SpiceError(
        :baddaf,
        "SPK type 10 segment at $(seg.baddr):$(seg.eaddr): $msg"
    )
)

"""
    _t10_check(data, seg) -> NamedTuple

Addresses of the parts of the generic segment `seg`, from its meta data (`sgmeta`):
`con` is the address before the first constant, `pkt` the address before the first
packet record, `ref` the address before the first reference epoch, `n` the count of
packets, `pksz` the packet size (10 or 14) and `stride` the size of a packet record
(the packet and the `pktoff` doubles that precede it).

The last double of the segment is the count of meta items. A segment with 15 or 16
items has no packet offset (item 16), as `sgmeta` reads the old files; with 17 or more
items the first 17 are used. Throws `SpiceError(:baddaf)` unless the segment has
fixed-size packets of 10 or 14 doubles, one reference epoch for each packet with the
reference index type 4 (closest epoch), at least 8 constants, and all its parts
inside the segment.
"""
function _t10_check(data::Vector{Float64}, seg::Segment)
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    (b >= 1 && e <= length(data) && len >= 16) || _t10_error(seg, "the segment is too short")
    nmeta = _seg_count(seg, @inbounds(data[e]), 15, _MAXCOUNT, "meta item count")
    m, hasoff = _sg_meta_base(nmeta, e)      # meta item k is data[m + k]
    span = e - m                             # doubles from the first meta item to the end
    span <= len || _t10_error(seg, "the meta data is longer than the segment")
    item(k, lo, what) = _seg_count(seg, @inbounds(data[m + k]), lo, _MAXCOUNT, what)
    conbas = item(1, 0, "constants base")
    ncon = item(2, 0, "constants count")
    rdrtyp = item(5, 0, "reference index type")
    refbas = item(6, 0, "reference base")
    nref = item(7, 0, "reference count")
    pdrtyp = item(10, 0, "packet directory type")
    pktbas = item(11, 0, "packet base")
    npkt = item(12, 1, "packet count")
    pksz = item(15, 1, "packet size")
    pktoff = hasoff ? item(16, 0, "packet offset") : 0
    ncon >= 8 || _t10_error(seg, "the segment has $ncon constants, 8 are needed")
    pdrtyp == 0 || _t10_error(seg, "the packets do not have a fixed size")
    (pksz == 10 || pksz == 14) || _t10_error(seg, "the packet size $pksz is not 10 or 14")
    rdrtyp == 4 || _t10_error(seg, "the reference index type $rdrtyp is not 4")
    nref == npkt || _t10_error(seg, "the counts of packets and reference epochs differ")
    stride = pksz + pktoff
    data_len = len - span
    (
        conbas + 8 <= data_len && pktoff <= data_len && npkt <= (data_len - pktbas) ÷ stride &&
            pktbas <= data_len && refbas <= data_len && nref <= data_len - refbas
    ) ||
        _t10_error(seg, "a part of the segment is outside the segment")
    return (;
        con = b - 1 + conbas, pkt = b - 1 + pktbas + pktoff, ref = b - 1 + refbas, n = npkt,
        pksz, stride,
    )
end

"`Segment.layout` of a type 10 segment `seg` with the parts `lay` of `_t10_check`: the addresses relative to the one before the segment."
@inline function _t10_encode(seg::Segment, lay)
    a = seg.baddr - 1
    return (
        Int32(lay.n), Int32(lay.con - a), Int32(lay.pkt - a), Int32(lay.ref - a),
        Int32(lay.pksz), Int32(lay.stride),
    )
end

"""
    _t10_layout(data, seg) -> NamedTuple

The parts `con`, `pkt`, `ref`, `n`, `pksz` and `stride` of the type 10 segment `seg`,
with the meaning of [`_t10_check`](@ref). They come from `seg.layout`, or from the
checks of `_t10_check` (and their errors) when the segment has no layout.
"""
@inline function _t10_layout(data::Vector{Float64}, seg::Segment)
    l = _has_layout(seg) ? seg.layout : _t10_encode(seg, _t10_check(data, seg))
    a = seg.baddr - 1
    @inbounds n, con, pkt, ref, pksz, stride = Int.(l)
    return (; con = a + con, pkt = a + pkt, ref = a + ref, n, pksz, stride)
end

# --- SPK type 19 ------------------------------------------------------------------------------

"Packet size of a subtype of types 18 and 19: 12 for Hermite with two derivative sets, else 6."
@inline _packet_size(subtyp::Int) = subtyp == 0 ? 12 : 6

"""
    _spk19_head(data, seg) -> (later, n, pbas, ivbas)

The selection flag `later` (an epoch on a bound takes the later interval), the
interval count `n` of the type 19 segment `seg`, and the addresses before its
mini-segment pointers (`pbas`) and before its interval bounds (`ivbas`). Throws
`SpiceError(:baddaf)` unless the flag is a finite number and the count agrees with the
length.
"""
@inline function _spk19_head(data::Vector{Float64}, seg::Segment)
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    len >= 6 || _badsegment(seg, "the segment is too short")
    # `spkr19` rounds the flag (ISEL = i_dnnt); the writer stores -1 for "earlier".
    later = _seg_count(seg, @inbounds(data[e - 1]), -_MAXCOUNT, _MAXCOUNT, "selection flag") == 1
    n = _seg_count(seg, @inbounds(data[e]), 1, _MAXCOUNT, "interval count")
    pbas = e - 2 - (n + 1)                      # pointer i is data[pbas+i]
    ivbas = pbas - n ÷ 100 - (n + 1)            # interval bound i is data[ivbas+i]
    # The last pointer is the place after the last mini-segment: the first bound.
    (ivbas >= b && @inbounds(data[pbas + n + 1]) == ivbas - b + 2) ||
        _badsegment(seg, "the interval count does not agree with the length")
    return later, n, pbas, ivbas
end

"""
    _spk19_mini(data, seg, pbas, ivbas, mini) -> (minib, minie, subtyp, wnd, npkt)

First and last address, subtype, window size and packet count of mini-segment `mini`
of the type 19 segment `seg` with the addresses `pbas` and `ivbas` of `_spk19_head`.
Throws `SpiceError(:baddaf)` unless the pointers and the packet count agree with the
length, and, in the order of `spkr19`, `:notsupported` for a subtype not in `0:2` and
`:invalidvalue` for a window size that is odd or not in `2:14` (subtype 0) or `2:28`.
"""
@inline function _spk19_mini(data::Vector{Float64}, seg::Segment, pbas::Int, ivbas::Int, mini::Int)
    b = seg.baddr
    len = seg.eaddr - b + 1
    @inbounds p1, p2 = data[pbas + mini], data[pbas + mini + 1]
    minib = _seg_count(seg, p1, 1, len, "mini-segment pointer") + b - 1
    minie = _seg_count(seg, p2, 1, len, "mini-segment pointer") + b - 2
    (minie - minib >= 3 && minie <= ivbas) ||
        _badsegment(seg, "the pointers of mini-segment $mini are not valid")
    # The checks and kinds of `spkr19`, in its order: the subtype, then the window.
    subtyp = _seg_count(seg, @inbounds(data[minie - 2]), -_MAXCOUNT, _MAXCOUNT, "subtype")
    wnd = _seg_count(seg, @inbounds(data[minie - 1]), -_MAXCOUNT, _MAXCOUNT, "window size")
    0 <= subtyp <= 2 || _segerror(:notsupported, seg, "the subtype $subtyp is not supported")
    maxwnd = subtyp == 0 ? 14 : 28
    (2 <= wnd <= maxwnd && iseven(wnd)) || _segerror(
        :invalidvalue, seg,
        "the window size $wnd is not even and in 2:$maxwnd for subtype $subtyp"
    )
    npkt = _seg_count(seg, @inbounds(data[minie]), 1, _MAXCOUNT, "packet count")
    npkt * (_packet_size(subtyp) + 1) + (npkt - 1) ÷ 100 + 3 == minie - minib + 1 ||
        _badsegment(seg, "the packet count of mini-segment $mini does not agree with its length")
    return minib, minie, subtyp, wnd, npkt
end

"""
    _spk19_check(data, seg) -> NTuple{6,Int32}

`Segment.layout` of the type 19 segment `seg`: the interval count, the addresses of
`_spk19_head` relative to the address before the segment, and the selection flag (1
for the later interval, else 0). Throws the `SpiceError` of `_spk19_head` or
`_spk19_mini` unless the head and every mini-segment pass their checks.
"""
function _spk19_check(data::Vector{Float64}, seg::Segment)
    later, n, pbas, ivbas = _spk19_head(data, seg)
    for mini in 1:n                             # invariant: mini-segments 1:mini-1 are valid
        _spk19_mini(data, seg, pbas, ivbas, mini)
    end
    a = seg.baddr - 1
    return (Int32(n), Int32(pbas - a), Int32(ivbas - a), Int32(later), Int32(0), Int32(0))
end

# --- at load ----------------------------------------------------------------------------------

"""
    _segment_layout(data, kind, seg) -> Segment

`seg`, a segment of a DAF of `kind` (`:SPK` or `:PCK`) over `data`, with its `layout`:
the layout of an SPK type 10, 14 or 19 or a binary PCK type 3 segment whose meta data
pass the checks of its evaluator, else zero. Never throws a `SpiceError`: a segment that
fails the checks throws at evaluation, as in CSPICE.
"""
function _segment_layout(data::Vector{Float64}, kind::Symbol, seg::Segment)
    @assert kind === :SPK || kind === :PCK
    zero = ntuple(_ -> Int32(0), Val(6))
    d = seg.dtype
    layout = if seg.eaddr < seg.baddr          # no data: the evaluation throws
        zero
    else
        try
            if kind === :SPK && d == 14 || kind === :PCK && d == 3
                _spk14_encode(_spk14_check(data, seg))
            elseif kind === :SPK && d == 10
                _t10_encode(seg, _t10_check(data, seg))
            elseif kind === :SPK && d == 19
                _spk19_check(data, seg)
            else
                zero
            end
        catch err
            err isa SpiceError || rethrow()
            zero
        end
    end
    @assert layout === zero || layout[1] >= 1
    return Segment(
        seg.target, seg.center, seg.frame, seg.dtype, seg.t0, seg.t1, seg.file,
        seg.baddr, seg.eaddr, seg.init, seg.intlen, seg.rsize, seg.nrec, layout
    )
end
