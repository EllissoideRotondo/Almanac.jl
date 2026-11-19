# Replay of the calls that NAIF's test suite (tspice) made to SPICELIB, against
# Almanac.jl, and comparison of the results.
#
#     julia --project=test test/tspice/replay.jl LOGDIR [FAMILY ...]
#
# LOGDIR holds one log per tspice family (`F_XXX.log.gz` or `F_XXX.log`, written by
# record.sh through the instrumented SPICELIB interface of PEPPER, see instrument.py),
# with the copies of the working directory that the log names. Each family is replayed in a new
# `SpiceContext`, call by call, in the order of the log. A call is:
#
# - `match`: Almanac gives the outputs of SPICELIB within `RTOL`, or throws the
#   `SpiceError` whose kind is the short message that SPICELIB signaled;
# - `cspice`: Almanac differs from SPICELIB, and so does CSPICE, which Almanac agrees with
#   (TSPICE_ORACLE=1);
# - `near`: the outputs differ by more than `RTOL` and less than `NEAR` relative (to the
#   largest magnitude in each group of three numbers);
# - `spread`: the outputs differ by more, SPICELIB and CSPICE differ from each other, and
#   Almanac is no farther from CSPICE than twice SPICELIB's distance (TSPICE_ORACLE=1);
# - `value`: the outputs differ by more;
# - `error`: SPICELIB signaled an error and Almanac did not, or the reverse, or the kinds
#   differ;
# - `adapter`: the call could not be made (a method error of the mapping below), which
#   is a fault of this script, not of Almanac;
# - `skipped`: made by another recorded routine, or made while SPICELIB was in an error
#   state (it then does nothing), or with an array too long for the log, or after a
#   load or pool edit with such an array, or a routine with no mapping.
#
# The output is a table per routine and the first mismatches of each kind, and
# `LOGDIR/replay.tsv` with one line per mismatch.

using Almanac, Printf, LinearAlgebra, StaticArrays

# With TSPICE_ORACLE=1 each mismatch is also computed by CSPICE (through SPICE.jl, with
# the loads of the family mirrored into it): a call where Almanac agrees with CSPICE and
# both differ from SPICELIB, the Fortran library that tspice ran, is counted as
# `cspice` instead (CSPICE is the reference of Almanac, and some of its routines are
# not translations of SPICELIB's: vperp_c, the time defaults).
const ORACLE = get(ENV, "TSPICE_ORACLE", "0") == "1"
ORACLE && import SPICE

const RTOL = parse(Float64, get(ENV, "TSPICE_RTOL", "1e-12"))
const ATOL = parse(Float64, get(ENV, "TSPICE_ATOL", "1e-14"))
const NEAR = parse(Float64, get(ENV, "TSPICE_NEAR", "1e-9"))

"Marks an argument whose values the log does not hold."
const TRUNCATED = :truncated

struct Arg
    name::String
    kind::Char          # 'd', 'i', 'l', 'c'
    values::Vector{Any}
end

struct Call
    name::String
    line::Int
    depth::Int
    failed_before::Bool
    inputs::Vector{Arg}
    outputs::Vector{Arg}
    failed::Bool
    short::String
    snapshot::String
    "A writer that the log does not record (a `W` line): only its handles are logged."
    unrecorded::Bool
end

"The log at `path`, decompressed if it ends in `.gz`."
open_log(f, path) = endswith(path, ".gz") ? open(f, `gzip -dc $path`) : open(f, path)

"""
Call `f(call)` for each call of the log at `path`, in order, reading the log as a stream.
A call of a routine for which `want(name)` is false is passed with no arguments, which
saves the parsing of its values.
"""
function foreach_call(f, path; want = _ -> true)
    return open_log(path) do io
        pending = nothing                   # a line read ahead
        lineno = 0
        nextline() = (pending === nothing ? (lineno += 1; readline(io)) : (l = pending; pending = nothing; l))
        while !eof(io) || pending !== nothing
            l = nextline()
            startswith(l, "@ ") || continue
            parts = split(l)
            length(parts) == 4 || continue
            _, name, depth, fb = parts
            start = lineno
            ins, outs, snap, unrec, bin = Arg[], Arg[], "", false, ""
            keep = want(lowercase(name))
            failed, short = false, ""
            while !eof(io)
                l = nextline()
                if startswith(l, "! ")
                    p = split(l, ' '; limit = 3)
                    failed = p[2] == "T"
                    short = length(p) == 3 ? strip(p[3]) : ""
                    break
                elseif startswith(l, "S ")
                    # the snapshot directory as recorded, or beside this log if the
                    # recording was moved
                    snap = l[3:end]
                    m = match(r"([^/]+\.log\.files/[^/]+)$", snap)
                    (!isdir(snap) && m !== nothing) && (snap = joinpath(dirname(path), m[1]))
                elseif startswith(l, "B ")
                    # the values of the next argument, too long for the log, in a binary file
                    bin = l[3:end]
                    m = match(r"([^/]+\.log\.files/bin/[^/]+)$", bin)
                    (!isfile(bin) && m !== nothing) && (bin = joinpath(dirname(path), m[1]))
                elseif l == "W"
                    unrec = keep = true     # its handles follow
                elseif !keep
                    continue
                elseif startswith(l, "<") || startswith(l, ">")
                    p = split(l)
                    length(p) == 3 || continue
                    aname, k, n = p
                    n = parse(Int, n)
                    if n < 0                # more values than the log keeps
                        vals = isempty(bin) || !isfile(bin) ? Any[TRUNCATED] : read_bin(bin, k[1])
                        bin = ""
                        push!((l[1] == '<' ? ins : outs), Arg(aname[2:end], k[1], vals))
                        continue
                    end
                    vals = Vector{Any}(undef, n)
                    for j in 1:n
                        v = nextline()
                        vals[j] = k == "d" ? parse(Float64, strip(v)) :
                            k == "i" ? parse(Int, strip(v)) :
                            k == "l" ? strip(v) == "T" :
                            # `%`: a string with characters outside printable ASCII, in hex
                            v[1] == '%' ? String(hex2bytes(v[2:end])) : String(v[2:end])
                    end
                    push!((l[1] == '<' ? ins : outs), Arg(aname[2:end], k[1], vals))
                elseif startswith(l, "@ ")
                    pending = l             # a nested call before this one ended
                    break
                end
            end
            f(Call(name, start, parse(Int, depth), fb == "T", ins, outs, failed, short, snap, unrec))
        end
    end
end

"""
The values of a binary file of the recorder (`B` line): doubles, 4-byte integers, or
strings of the length that the name gives (`.c80.gz`), without their trailing blanks as
the log has them.
"""
function read_bin(path, kind)
    bytes = read(`gzip -dc $path`)
    kind == 'd' && return Any[reinterpret(Float64, bytes)...]
    kind == 'i' && return Any[Int.(reinterpret(Int32, bytes))...]
    m = match(r"\.c(\d+)\.gz$", path)
    m === nothing && return Any[TRUNCATED]
    l = parse(Int, m[1])
    l == 0 && return Any["" for _ in 1:0]
    return Any[rstrip(String(bytes[(k + 1):(k + l)]), ' ') for k in 0:l:(length(bytes) - l)]
end

# --- argument conversion -----------------------------------------------------------------

"A logged value as Almanac takes it: a scalar, a static vector or a 3x3 or 6x6 matrix."
function jl(a::Arg)
    v = a.values
    a.kind == 'c' && return length(v) == 1 ? v[1] : String.(v)
    n = length(v)
    T = a.kind == 'd' ? Float64 : a.kind == 'i' ? Int : Bool
    n == 1 && return T(v[1])
    n == 9 && return SMatrix{3, 3, T}(v...)
    n == 36 && return SMatrix{6, 6, T}(v...)
    # past the length of the log (values from a binary file): a plain vector
    n > 2000 && return T.(v)
    return SVector{n, T}(v...)
end

"The values of an Almanac result, flattened as SPICELIB's output arguments are."
flat(x::Nothing) = nothing
flat(x::AbstractString) = Any[String(x)]
flat(x::Symbol) = Any[String(x)]
flat(x::Number) = Any[x]
flat(x::AbstractArray{<:Number}) = Any[vec(collect(x))...]
flat(x::AbstractArray{<:AbstractString}) = Any[String.(x)...]
flat(x::AbstractArray) = Any[reduce(vcat, flat.(x))...]
flat(x::SpiceEllipse) = Any[x.center..., x.semimajor..., x.semiminor...]
flat(x::SpicePlane) = Any[x.normal..., x.constant]
flat(x::Tuple) = Any[reduce(vcat, [flat(v) for v in x]; init = Any[])...]
# the descriptors in the order of the arrays of SPICELIB (DLADSZ = 8, DSKDSZ = 24)
flat(x::Almanac.DLADescr) = Any[getfield(x, f) for f in fieldnames(typeof(x))]
flat(x::Almanac.DSKDescr) = Any[
    x.surfce, x.center, x.dclass, x.dtype, x.frmcde, x.corsys,
    x.corpar..., x.co1min, x.co1max, x.co2min, x.co2max, x.co3min, x.co3max, x.start, x.stop,
]
# the descriptors of CSPICE have the same fields
ORACLE && @eval begin
    flat(x::SPICE.DLADescr) = Any[getfield(x, f) for f in fieldnames(typeof(x))]
    flat(x::SPICE.Plane) = Any[x.normal..., x.constant]
    flat(x::SPICE.Ellipse) = Any[x.center..., x.semi_major..., x.semi_minor...]
    flat(x::SPICE.DSKDescr) = Any[
        x.surfce, x.center, x.dclass, x.dtype, x.frmcde, x.corsys,
        x.corpar..., x.co1min, x.co1max, x.co2min, x.co2max, x.co3min, x.co3max, x.start, x.stop,
    ]
end

close_enough(a::Real, b::Real) = a == b || abs(a - b) <= ATOL + RTOL * max(abs(a), abs(b))
close_enough(a::AbstractString, b::AbstractString) = rstrip(a) == rstrip(b)
close_enough(a::Bool, b::Bool) = a == b
close_enough(a, b) = false

"Largest relative difference of two value lists (for the report)."
function maxdiff(exp, got)
    m = 0.0
    for (a, b) in zip(exp, got)
        (a isa Real && b isa Real) || continue
        m = max(m, abs(a - b) / max(abs(a), abs(b), 1.0e-300))
    end
    return m
end

# --- the mapping of each routine ---------------------------------------------------------
#
# A mapping takes the context and the input arguments by name, and returns the list of
# expected/actual output pairs to compare: `[(name, expected values, actual values)]`.
# The default mapping calls the Almanac function of the same name with the inputs in the
# order of SPICELIB, with the context first when the function has a method for that, and
# pairs its results with the outputs of SPICELIB in order. A `found` output of SPICELIB is
# `nothing` in Almanac.


# Input arguments that SPICELIB needs and Almanac does not (sizes of output arrays).
const DROP = Dict(
    "gdpool" => ["start", "room"], "gipool" => ["start", "room"],
    "gcpool" => ["start", "room"], "bodvrd" => ["maxn"], "bodvcd" => ["maxn"],
    "pdpool" => ["n"], "pipool" => ["n"], "pcpool" => ["n"], "lmpool" => ["n"], "ekbseg" => ["ncols"]
)

for g in (
        "vnormg", "vminug", "vsubg", "vaddg", "vsclg", "vequg", "vzerog", "vdotg", "vlcomg", "dasadd", "dasadi",
        "vhatg", "vdistg", "vsepg", "vrelg", "vprojg", "vlcomg", "unormg",
    )
    DROP[g] = ["ndim"]
end

for (r, names) in (
        "qderiv" => ["ndim"], "srfnrm" => ["npts"], "latsrf" => ["npts"],
        "hrmint" => ["n", "work"], "lgrint" => ["n", "work"], "lgrind" => ["n", "work"],
        "dskrb2" => ["nv", "np"], "hrmesp" => ["n", "work"], "lgresp" => ["n", "work"],
        "dafada" => ["n"], "vlcomg" => ["n"], "ekacei" => ["nvals"], "ekaced" => ["nvals"],
        "ekacec" => ["nvals"], "ekucei" => ["nvals"], "ekuced" => ["nvals"], "ekucec" => ["nvals"],
        "dasadd" => ["n"], "dasadi" => ["n"], "dskxsi" => ["nsurf", "maxd", "maxi"],
        "dskxv" => ["nsurf", "nrays"], "chbigr" => ["degp"], "chbint" => ["degp"],
        "chbval" => ["degp"], "chbder" => ["degp"], "dafps" => ["nd", "ni"],
        "ckw01" => ["avflag", "nrec"], "ckw02" => ["nrec"], "ckw03" => ["avflag", "nrec", "nints"],
        "pltvol" => ["nv", "np"], "pltar" => ["nv", "np"], "pltexp" => [], "dskw02" => ["nv", "np"],
    )
    DROP[r] = names
end

"Routines whose arrays of 9 or 36 numbers are vectors, not matrices, in Almanac."
const VECTOR_ROUTINES = Set(
    [
        "vnormg", "vminug", "vsubg", "vaddg", "vsclg", "vequg",
        "vzerog", "vdotg", "vhatg", "vdistg", "vsepg", "vrelg", "vprojg", "vlcomg", "unormg",
        "qderiv", "eqncpv", "hrmint", "lgrint", "lgrind", "hrmesp", "lgresp", "dafada", "dasadd", "dasadi", "chbigr", "chbint",
        "chbval", "chbder", "dafps", "dskxsi",
    ]
)

"Routines whose array of points is a list of points in Almanac: name => (argument, size)."
const POINT_LISTS = Dict("srfnrm" => ("srfpts", 3), "latsrf" => ("lonlat", 2))

"A matrix argument of a routine of general dimensions: column-major, `nr` by `nc`."
gmat(a::Arg, nr, nc) = reshape(Float64.(a.values), nr, nc)
gmatn(a::Arg, nr, nc) = reshape(Float64.(a.values)[1:(nr * nc)], nr, nc)

# Output arguments that Almanac does not return (counts of returned arrays).
const NOTOUT = Dict(
    "gdpool" => ["n"], "gipool" => ["n"], "gcpool" => ["n"],
    "dskp02" => ["n"], "dskv02" => ["n"], "dski02" => ["n"], "dskd02" => ["n"],
    "ekrcec" => ["segno", "nvals"], "ekrced" => ["segno", "nvals"], "ekrcei" => ["segno", "nvals"],
    "bodvrd" => ["dim"], "bodvcd" => ["dim"], "dtpool" => ["found"],
    "lgrint" => ["work"], "lgrind" => ["work"], "hrmint" => ["work"], "lgresp" => ["work"],
    "hrmesp" => ["work"], "getfvn" => ["n"], "getfov" => ["n"]
)

"An input by its name: the planes and ellipses of SPICELIB are arrays of 4 and 9 doubles."
function jl_named(a::Arg)
    # the DLA descriptor of SPICELIB, 8 integers, is a DLADescr in Almanac
    if a.kind == 'i' && length(a.values) == 8 && a.name in ("dladsc", "nxtdsc", "prvdsc")
        return Almanac.DLADescr(Int.(a.values)...)
    end
    if a.kind == 'd' && length(a.values) == 4 && a.name in ("plane", "projpl", "invpl")
        v = Float64.(a.values)
        return SpicePlane(SVector{3}(v[1], v[2], v[3]), v[4])
    elseif a.kind == 'd' && length(a.values) == 9 && a.name in ("ellips", "ellipse", "limb")
        v = Float64.(a.values)
        return SpiceEllipse(SVector{3}(v[1:3]...), SVector{3}(v[4:6]...), SVector{3}(v[7:9]...))
    end
    return jl(a)
end
"Segment writers: the count `n` is the length of the arrays in Almanac."
const WRITERS = r"^(spkw|ckw|pckw)\d\d$"
"Array arguments of the segment writers, vectors even when they hold one number."
const WRITER_ARRAYS = Set(
    [
        "states", "epochs", "sclkdp", "quats", "avvs", "starts", "start", "stop",
        "cdata", "packts", "records", "rates", "coeffs", "dlines", "elems", "ivlbds",
    ]
)

"Routines whose vertices and plates are 3 x n matrices in Almanac."
const PLATE_ROUTINES = Set(["pltar", "pltvol", "dskw02", "dskmi2", "dskrb2"])

function jl_named(a::Arg, routine::AbstractString)
    if routine in PLATE_ROUTINES && a.name in ("vrtces", "plates") && a.kind != 'c'
        v = a.kind == 'd' ? Float64.(a.values) : Int.(a.values)
        return reshape(v, 3, length(v) ÷ 3)
    end
    routine == "dskw02" && a.name == "corpar" && length(a.values) > 10 &&
        return Float64.(a.values[1:10])
    if occursin(WRITERS, routine) && a.name in WRITER_ARRAYS && a.kind == 'd'
        return Float64.(a.values)
    end
    # the values of an EK entry are a vector, also of one value
    if occursin(r"^ek[au]ce[cdi]$", routine) && a.name in ("ivals", "dvals", "cvals")
        return a.kind == 'c' ? String.(a.values) : a.kind == 'i' ? Int.(a.values) : Float64.(a.values)
    end
    # the log holds the declared size of CONSTS(8) plus one
    routine == "spkw10" && a.name == "consts" && length(a.values) > 8 &&
        return Float64.(a.values[1:8])
    pl = get(POINT_LISTS, routine, nothing)
    if pl !== nothing && a.name == pl[1]
        k = pl[2]
        v = Float64.(a.values)
        return [SVector{k}(v[i:(i + k - 1)]) for i in 1:k:(length(v) - k + 1)]
    end
    routine in VECTOR_ROUTINES && a.kind != 'c' &&
        (
        length(a.values) > 1 || a.name in (
            "xvals", "yvals", "f0", "f2", "v1", "v2", "data",
            "cp", "srflst", "dc", "ic", "vrtces", "plates",
        )
    ) &&
        return length(a.values) > 2000 ? (a.kind == 'd' ? Float64.(a.values) : Int.(a.values)) :
        a.kind == 'd' ? SVector{length(a.values), Float64}(a.values...) :
        SVector{length(a.values), Int}(a.values...)
    return jl_named(a)
end

# --- files opened by handle ------------------------------------------------------------
#
# The integer handles of SPICELIB in the log are mapped to the handles that Almanac
# returned for the same open. A file is opened in the working directory of the replay
# (a new directory for each family), not in the snapshot of the call: a snapshot is a
# record of SPICELIB's directory and stays unchanged. A file that the open reads is
# first copied from the snapshot.

const OPENERS = r"^(daf|das|spk|ck|pck|ek|dsk|dla)(opn|opr|opw|onw|opa)$"
const CLOSERS = r"^(daf|das|spk|ck|pck|ek|dsk)cls$"
"Opens that create a file; SPICELIB's succeeded, so a file of that name did not exist."
const CREATORS = r"(opn|onw)$"
const HANDLE_ARGS = Set(["handle", "inhan", "han", "handl", "newh"])

"A handle of a file that SPICELIB opened in a routine that is not recorded (a tspice
utility such as TSTSPK): the call is skipped; the file is read from the snapshots."
struct UnmappedHandle <: Exception end

"Paths of the files opened by handle, and the files whose writes were not all replayed."
const HPATH = Dict{Int, String}()
const TAINTED = Set{String}()

"Almanac's handle of each handle of SPICELIB."
const HMAP = Dict{Int, Any}()

handle_value(a::Arg) = a.kind == 'i' && length(a.values) == 1 && a.name in HANDLE_ARGS ?
    get(HMAP, a.values[1], nothing) : nothing

"""
A reader handle of Almanac for a file that SPICELIB loaded (`cklpf`, `spklef`): the DAF
calls take a handle where the queries take the path. Opened once, in the directory where
the file was loaded.
"""
const READERS = Dict{String, Any}()
"The directory in which each file loaded by handle was loaded, and that of the last furnsh."
const LOADDIR = Dict{String, String}()
const LAST_LOAD = Ref("")
function daf_reader(path::AbstractString, c::Call)
    return get!(READERS, path) do
        dir = get(LOADDIR, path, isempty(c.snapshot) ? pwd() : c.snapshot)
        p = isabspath(path) ? path : joinpath(dir, path)
        first(Almanac.getfat(p)) == "DAS" ? Almanac.dasopr(p) : Almanac.dafopr(p)
    end
end

function open_by_handle(ctx, c::Call)
    name = lowercase(c.name)
    path = String(c.inputs[1].values[1])
    if occursin(CREATORS, name)
        c.failed || (rm(path; force = true); delete!(TAINTED, abspath(path)))
    else
        snap = isempty(c.snapshot) ? "" : joinpath(c.snapshot, path)
        if (!isfile(path) || abspath(path) in TAINTED) && isfile(snap)
            mkpath(dirname(abspath(path)))
            cp(snap, path; force = true)
            delete!(TAINTED, abspath(path))
        end
        abspath(path) in TAINTED && throw(UnmappedHandle())
    end
    f = getfield(Almanac, Symbol(name))
    # SPICELIB takes an empty internal file name and CSPICE does not (`emptystring`);
    # the name is not under test, and without the file the family could not go on.
    args = Any[a.name == "ifname" && isempty(strip(jl(a))) ? "replay" : jl(a) for a in c.inputs]
    any(a -> a.name == "ifname" && isempty(strip(jl(a))), c.inputs) ? push!(IFNAME_REPLAY, abspath(path)) :
        occursin(CREATORS, name) && delete!(IFNAME_REPLAY, abspath(path))
    h = f(args...)
    k = findfirst(o -> o.name in HANDLE_ARGS, c.outputs)
    k === nothing || (HMAP[c.outputs[k].values[1]] = h; HPATH[c.outputs[k].values[1]] = abspath(path))
    (k !== nothing && occursin(r"(opn|onw|opw|opa)$", name)) && push!(WRITE_OPEN, c.outputs[k].values[1])
    return []
end

"What the close of a file that SPICELIB wrote did with it, for files.tsv."
const FILE_CHECK = Ref("")
"Why the replay gave up each handle (drop_handle), for files.tsv."
const DROPPED = Dict{Int, String}()

"The files that Almanac created with the internal name `replay` for SPICELIB's blank one."
const IFNAME_REPLAY = Set{String}()

"The EK fast loads that Almanac's ekifld began, by SPICELIB's handle and segment."
const EKFAST = Dict{Tuple{Int, Int}, Any}()

"""
SPKW19 and CKW06: the counts of the mini-segments are vectors in Almanac (also of one
value), and the count `nname` is their length.
"""
function minisegment_writer(ctx, c, f, nname, dbl)
    a = Dict(x.name => x for x in c.inputs)
    any(x -> any(v -> v isa Symbol, x.values), values(a)) && throw(UnmappedHandle())
    h = handle_value(a["handle"]); h === nothing && throw(UnmappedHandle())
    n = jl(a[nname])
    n < 0 && throw(UnmappedHandle())
    args = Any[]
    for x in c.inputs
        x.name in ("handle", nname) && continue
        push!(
            args, x.name in ("npkts", "subtps", "degres") ? Int.(x.values)[1:n] :
                x.name in dbl ? Float64.(x.values) : jl(x)
        )
    end
    applicable(f, ctx, h, args...) ? f(ctx, h, args...) : f(h, args...)
    return []
end

"EKACLC, EKACLD, EKACLI: the column of the fast load that EKIFLD began."
function ek_column(c, f, T)
    a = Dict(x.name => x for x in c.inputs)
    any(x -> any(v -> v isa Symbol, x.values), values(a)) && throw(UnmappedHandle())
    h = handle_value(a["handle"]); h === nothing && throw(UnmappedHandle())
    segno = jl(a["segno"])
    rc = get(EKFAST, (a["handle"].values[1], segno), nothing)
    rc === nothing && throw(UnmappedHandle())
    vk = only(k for k in keys(a) if k in ("cvals", "dvals", "ivals"))
    vals = T === String ? String.(a[vk].values) : T.(a[vk].values)
    f(h, segno, jl(a["column"]), vals, Int.(a["entszs"].values), Bool.(a["nlflgs"].values), rc)
    return []
end

"DAFWS and DAFRS: the summary, logged unpacked, of the current array of the search."
function daf_summary_writer(c, f)
    a = Dict(x.name => x for x in c.inputs)
    h = CURRENT_DAF[]
    # a search that the replay does not have: Almanac's file would miss the write
    h === nothing && (c.failed || foreach(drop_handle, collect(WRITE_OPEN)); throw(UnmappedHandle()))
    f(h, Almanac.dafps(Float64.(a["dc"].values), Int.(a["ic"].values)))
    return []
end

"The type 14 segment writers of Almanac, by SPICELIB's handle."
const SPK14W = Dict{Int, Any}()

"SPICELIB's handles of the files open for write."
const WRITE_OPEN = Set{Int}()

"""
DASLLC closes the DAS file on a logical unit, after the test wrote it with Fortran I/O,
which is not in the log: each file open for write is closed in Almanac, and the later
opens read SPICELIB's file from the snapshots.
"""
function dasllc_mapping(ctx, c)
    foreach(drop_handle, collect(WRITE_OPEN))
    return []
end

"""
Give up SPICELIB's handle `hid`: Almanac's file misses data that SPICELIB's has. The
file is closed in Almanac, so that a later open does not find it open, and later opens
read SPICELIB's file from the snapshots.
"""
function drop_handle(hid, why = "a write that the replay could not make")
    haskey(HMAP, hid) && (DROPPED[hid] = why)
    haskey(HPATH, hid) && push!(TAINTED, HPATH[hid])
    h = get(HMAP, hid, nothing)
    h isa Almanac.DASHandle && try
        Almanac.dascls(h)
    catch end
    h isa Almanac.DAFHandle && try
        Almanac.dafcls(h)
    catch end
    CURRENT_DAF[] === h && (CURRENT_DAF[] = nothing)
    delete!(HMAP, hid)
    return delete!(WRITE_OPEN, hid)
end

# SPICELIB keeps a current DAF for the array search and the array being added; Almanac
# takes the handle of each call.
const DAF_SELECT = Set(["dafbfs", "dafbbs", "dafbna", "dafcs", "dafcad"])
const DAF_CURRENT = Set(["daffna", "daffpa", "dafgs", "dafgn", "dafena", "dafada", "dafws", "dafrs", "dafrn"])
const CURRENT_DAF = Ref{Any}(nothing)

"""
The characters that DASADC and DASUDC of SPICELIB take from the strings `data`, the
positions `bpos:epos` of each string in turn, `n` of them: Almanac takes them as one
string.
"""
function das_chars(data, n, bpos, epos)
    io = IOBuffer()
    for s in data
        n <= 0 && break
        sub = rpad(s, epos)[bpos:epos]
        k = min(n, length(sub))
        print(io, sub[1:k])
        n -= k
    end
    return String(take!(io))
end

"The writers of records of one size, with the array of records and the size of a record."
const RECORD_WRITERS = Dict(
    "spkw02" => ("cdata", a -> 3 * (jl(a["polydg"]) + 1)),
    "spkw03" => ("cdata", a -> 6 * (jl(a["polydg"]) + 1)),
    "spkw20" => ("cdata", a -> 3 * (jl(a["polydg"]) + 2)),
    "pckw02" => ("cdata", a -> 3 * (jl(a["polydg"]) + 1)),
    "pckw20" => ("cdata", a -> 3 * (jl(a["polydg"]) + 2)),
    "spkw21" => ("dlines", a -> jl(a["dlsize"]))
)

const EK_QUERY = Set(["ekfind", "eknelt", "ekgc", "ekgd", "ekgi"])
const CURRENT_EK = Ref{Any}(nothing)   # result of the last ekfind

function default_mapping(ctx, c::Call)
    name = lowercase(c.name)
    if name == "dasrdc"
        a = Dict(x.name => x for x in c.inputs)
        h = handle_value(a["handle"]); h === nothing && throw(UnmappedHandle())
        f, l = jl(a["first"]), jl(a["last"])
        got = dasrdc(h, f, l)
        (c.failed || isempty(c.outputs) || any(v -> v isa Symbol, c.outputs[1].values)) && return []
        exp = das_chars(c.outputs[1].values, l - f + 1, jl(a["bpos"]), jl(a["epos"]))
        return [("data", Any[exp], Any[got])]
    end
    if name in ("dasadc", "dasudc")
        a = Dict(x.name => x for x in c.inputs)
        any(x -> x.kind != 'c' && any(v -> v isa Symbol, x.values), values(a)) && throw(UnmappedHandle())
        h = handle_value(a["handle"]); h === nothing && throw(UnmappedHandle())
        data = a["data"].values
        b, e = jl(a["bpos"]), jl(a["epos"])
        if name == "dasadc"
            dasadc(h, das_chars(data, jl(a["n"]), b, e))
        else
            f, l = jl(a["first"]), jl(a["last"])
            dasudc(h, f, l, das_chars(data, l - f + 1, b, e))
        end
        return []
    end
    if haskey(RECORD_WRITERS, name)
        # SPICELIB: N records of a size that the degree (or DLSIZE) gives, in one array;
        # Almanac: a vector of the records, whose size gives the degree
        data, sizeof = RECORD_WRITERS[name]
        a = Dict(x.name => x for x in c.inputs)
        any(x -> x.kind != 'c' && any(v -> v isa Symbol, x.values), values(a)) && throw(UnmappedHandle())
        h = handle_value(a["handle"]); h === nothing && throw(UnmappedHandle())
        n, m = jl(a["n"]), sizeof(a)
        v = Float64.(a[data].values)
        # a size or count that a vector of records cannot carry, or a log without all
        # values (then Almanac's file misses the segment)
        (m <= 0 || n < 0) && throw(UnmappedHandle())
        # degree -1 of the type 20 writers (records of 3), or no record of SPKW21: the vector of records gives Almanac
        # the degree or the size of a difference line, which it then cannot carry
        (name in ("pckw20", "spkw20") && jl(a["polydg"]) < 0) && throw(UnmappedHandle())
        (name == "spkw21" && n == 0) && throw(UnmappedHandle())
        length(v) < n * m && (c.failed || drop_handle(a["handle"].values[1]); throw(UnmappedHandle()))
        recs = [v[((k - 1) * m + 1):(k * m)] for k in 1:n]
        drop = ("handle", "n", "polydg", "dlsize")
        args = Any[x.name == data ? recs : jl_named(x, name) for x in c.inputs if !(x.name in drop)]
        f = getfield(Almanac, Symbol(name))
        applicable(f, ctx, h, args...) ? f(ctx, h, args...) : f(h, args...)
        return []
    end
    if name == "dskxsi"
        # Almanac returns the point, the segment and the plate; SPICELIB the point, the
        # handle, the descriptors and the plate as IC(1): the point and the plate are
        # compared, and whether there is an intercept
        a = Dict(x.name => x for x in c.inputs)
        any(x -> x.kind != 'c' && any(v -> v isa Symbol, x.values), values(a)) && throw(UnmappedHandle())
        # a negative surface count or output arrays too small for one value: Almanac
        # takes a vector and returns a tuple, so the call has no counterpart
        (jl(a["nsurf"]) < 0 || jl(a["maxd"]) < 1 || jl(a["maxi"]) < 1) && throw(UnmappedHandle())
        srf = Int.(a["srflst"].values)[1:max(0, jl(a["nsurf"]))]
        r = dskxsi(
            ctx, jl(a["pri"]), jl(a["target"]), srf, jl(a["et"]), jl(a["fixref"]),
            jl(a["vertex"]), jl(a["raydir"])
        )
        out = Dict(o.name => o for o in c.outputs)
        fexp = out["found"].values[1]
        r === nothing && return [("found", Any[fexp], Any[false])]
        fexp || return [("found", Any[false], Any[true])]
        return [("xpt", out["xpt"].values, flat(r[1])), ("ic", out["ic"].values[1:1], Any[r[3]])]
    end
    if name in EK_QUERY
        a = Dict(x.name => x for x in c.inputs)
        out = Dict(o.name => o for o in c.outputs)
        if name == "ekfind"
            res, failed, msg = ekfind(ctx, jl(a["query"]))
            CURRENT_EK[] = failed ? nothing : res
            pairs = Any[("error", out["error"].values, Any[failed])]
            failed || push!(pairs, ("nmrows", out["nmrows"].values, Any[res.nmrows]))
            return pairs
        end
        # the fetch routines read the result of the last EKFIND; Almanac takes it as an argument
        CURRENT_EK[] === nothing && throw(UnmappedHandle())
        res, s, r = CURRENT_EK[], jl(a["selidx"]), jl(a["row"])
        name == "eknelt" && return [("nelt", out["nelt"].values, Any[eknelt(res, s, r)])]
        v, isnull, found = getfield(Almanac, Symbol(name))(res, s, r, jl(a["elment"]))
        pairs = Any[("found", out["found"].values, Any[found])]
        found || return pairs
        push!(pairs, ("null", out["null"].values, Any[isnull]))
        isnull && return pairs
        dk = only(k for k in keys(out) if k in ("cdata", "ddata", "idata"))
        exp = dk == "cdata" ? Any[rstrip(String(out[dk].values[1]))] : out[dk].values
        push!(pairs, (dk, exp, Any[v]))
        return pairs
    end
    if name == "ekcii"
        # ATTDSC: class, type code, length, size, index flag, null flag; the length of a
        # character column differs from CSPICE as EKColumnAttr documents and is not compared
        a = Dict(x.name => x for x in c.inputs)
        out = Dict(o.name => o for o in c.outputs)
        col, at = ekcii(ctx, jl(a["table"]), jl(a["cindex"]))
        flag(b) = b ? 1 : -1
        d = out["attdsc"].values
        got = Any[at.cclass, findfirst(==(at.dtype), (:CHR, :DP, :INT, :TIME)), at.size, flag(at.indexd), flag(at.nullok)]
        return [
            ("column", Any[rstrip(String(out["column"].values[1]))], Any[col]),
            ("attdsc", Any[d[1], d[2], d[4], d[5], d[6]], got),
        ]
    end
    occursin(OPENERS, name) && return open_by_handle(ctx, c)
    if name in DAF_SELECT
        h = handle_value(c.inputs[1])
        h isa AbstractString && (h = daf_reader(h, c))
        # a file that Almanac does not have: the search that follows is skipped
        CURRENT_DAF[] = h
        h === nothing && throw(UnmappedHandle())
        name in ("dafcs", "dafcad") && return []
    end
    name == "dafgh" && return map_handles([("handle", c.outputs[1].values, Any[CURRENT_DAF[]])])
    any(a -> a.name in HANDLE_ARGS && a.kind == 'i' && handle_value(a) === nothing, c.inputs) &&
        throw(UnmappedHandle())
    f = getfield(Almanac, Symbol(name))
    drop = occursin(WRITERS, name) ? vcat(get(DROP, name, String[]), ["n"]) : get(DROP, name, String[])
    args = Any[
        (h = handle_value(a); h === nothing ? jl_named(a, name) : h)
            for a in c.inputs if !(a.name in drop)
    ]
    # the angular velocities of a CK writer without them are `nothing` in Almanac
    if name in ("ckw01", "ckw03")
        av = findfirst(a -> a.name == "avflag", c.inputs)
        k = findfirst(a -> a.name == "avvs", [a for a in c.inputs if !(a.name in drop)])
        (av !== nothing && k !== nothing && !c.inputs[av].values[1]) && (args[k] = nothing)
    end
    # a loaded file is known by its path; a call that takes no path takes a reader handle
    if !applicable(f, ctx, args...) && !applicable(f, args...)
        args = Any[x isa AbstractString && x in values(HMAP) ? daf_reader(x, c) : x for x in args]
    end
    if name in DAF_CURRENT
        if CURRENT_DAF[] === nothing
            # the array was begun by a routine outside the log: a file open for write
            # (negative handle) misses data that SPICELIB's has
            name in ("dafada", "dafena") &&
                foreach(h -> h < 0 && push!(TAINTED, HPATH[h]), keys(HPATH))
            throw(UnmappedHandle())
        end
        pushfirst!(args, CURRENT_DAF[])
    end
    if occursin(CLOSERS, name)
        hid = c.inputs[1].values[1]
        try
            f(args...)
        catch e
            # SPICELIB closed a file with segments that a routine outside the log wrote
            # (the CK and PCK writers of some families are not recorded): the later
            # opens read SPICELIB's file from the snapshots.
            (e isa SpiceError && e.kind === :nosegmentsfound && !c.failed) || rethrow()
            try
                Almanac.dafcls(args[1])
            catch end    # the close left it open
            haskey(HPATH, hid) && push!(TAINTED, HPATH[hid])
            delete!(HMAP, hid)
            FILE_CHECK[] = "not compared: Almanac's file has no segment (written outside the log)"
            throw(UnmappedHandle())
        end
        filter!(kv -> kv.second !== args[1], READERS)
        args[1] isa AbstractString || delete!(HMAP, hid)
        written = hid in WRITE_OPEN
        delete!(WRITE_OPEN, hid)
        # a file that SPICELIB wrote: the recorder copied it after the close (logs from
        # before that have no copy, and then the close compares nothing)
        (written && !isempty(c.snapshot) && !c.failed && haskey(HPATH, hid)) || return []
        path = HPATH[hid]
        theirs = joinpath(c.snapshot, basename(path))
        # a write that the replay could not make, or no copy (past the cap of the
        # recorder): the files are not compared
        FILE_CHECK[] = path in TAINTED ? "not compared: a write outside the replay" :
            !isfile(theirs) ? "not compared: no copy of SPICELIB's file (past the cap of the recorder, or deleted before the close)" : "compared"
        FILE_CHECK[] == "compared" || throw(UnmappedHandle())
        return compare_files(path, theirs; ifname = path in IFNAME_REPLAY ? "replay" : nothing)
    end
    r = try
        applicable(f, ctx, args...) ? f(ctx, args...) : f(args...)
    catch e
        # An EK that SPICELIB filled through routines outside the log (the segments of
        # the tspice EK utilities): a segment or record that Almanac's file lacks.
        if e isa SpiceError && e.kind === :invalidindex && !c.failed && startswith(name, "ek")
            for a in c.inputs
                (a.name in HANDLE_ARGS && haskey(HPATH, a.values[1])) || continue
                drop_handle(a.values[1])
            end
            throw(UnmappedHandle())
        end
        rethrow()
    end
    notout = get(NOTOUT, name, String[])
    outs = [o for o in c.outputs if !(o.name in notout)]
    # the arguments that PEPPER declares in-out and that the EK routines leave unchanged
    occursin(r"^ek(ace|uce|rce)[cdi]$|^ekdelr$", name) &&
        filter!(o -> !any(a -> a.name == o.name && a.values == o.values, c.inputs), outs)
    # a routine without outputs (a setter that Almanac makes return its value), and the
    # SET action of TIMDEF, whose VALUE SPICELIB leaves as given
    (isempty(outs) || (name == "timdef" && uppercase(strip(jl(c.inputs[1]))) == "SET")) &&
        return []
    found = findfirst(o -> o.name == "found", outs)
    if found !== nothing
        fexp = outs[found].values[1]
        deleteat!(outs, found)
        r isa Bool && return [("found", Any[fexp], Any[r])]
        r === nothing && return [("found", Any[fexp], Any[false])]
        fexp || return [("found", Any[false], Any[true])]
    end
    vals = r isa Tuple ? collect(r) : Any[r]
    # a null EK entry: SPICELIB leaves the values undefined, Almanac returns none; only
    # the null flag is compared
    if occursin(r"^ekrce[cdi]$", name) && length(vals) == length(outs)
        k = findfirst(o -> o.name == "isnull", outs)
        if k !== nothing && outs[k].values == Any[true]
            keep = [i for i in eachindex(outs) if !(outs[i].name in ("cvals", "dvals", "ivals"))]
            outs, vals = outs[keep], vals[keep]
        end
    end
    WRITE_HANDLE[] = Int[a.values[1] for a in c.inputs if a.name in HANDLE_ARGS && a.kind == 'i']
    pairs = counted(name, map_handles(output_pairs(name, outs, vals)))
    # kinfo and kdata give the handle of a loaded file with its name: from here the handle
    # stands for the path, as for a file that a loader loaded
    if name in ("kinfo", "kdata") && !c.failed
        fi = findfirst(a -> a.name == "file", c.inputs)
        fo = findfirst(o -> o.name == "file", c.outputs)
        ho = findfirst(o -> o.name == "handle", c.outputs)
        path = fi !== nothing ? jl(c.inputs[fi]) : fo !== nothing ? jl(c.outputs[fo]) : nothing
        if path isa AbstractString && ho !== nothing
            HMAP[c.outputs[ho].values[1]] = String(path)
            # a file that a meta-kernel listed: read where the last furnsh read
            haskey(LOADDIR, path) || isempty(LAST_LOAD[]) || (LOADDIR[String(path)] = LAST_LOAD[])
        end
    end
    return pairs
end

# --- the files that the writers wrote -------------------------------------------------
#
# At the close of a file that SPICELIB wrote, the file that Almanac wrote is compared with
# the copy of SPICELIB's that the recorder made after the close, through Almanac's
# readers: the content, not the bytes, since Almanac may lay out the records otherwise
# (a DAS file is segregated at the close, the reserved records of a DAF are not kept).
# The result is that of the close: a match, or a value mismatch whose name says what
# differs.

"""
Add the pair `name` of SPICELIB's values `e` and Almanac's `g` to `pairs`: the lengths
if they agree, else the first element that differs, with its index in the name.
"""
function file_pair!(pairs, name, e, g)
    e, g = Any[e...], Any[g...]
    if length(e) != length(g)
        push!(pairs, ("$name: length", Any[length(e)], Any[length(g)]))
        return pairs
    end
    k = findfirst(i -> !close_enough(e[i], g[i]), eachindex(e))
    push!(pairs, k === nothing ? (name, Any[length(e)], Any[length(g)]) : ("$name[$k of $(length(e))]", Any[e[k]], Any[g[k]]))
    return pairs
end

"The arrays of a DAF in the order of the file: (dc, ic, name, data)."
function daf_arrays(h)
    nd, ni = Almanac.dafhsf(h)
    out = []
    Almanac.dafbfs(h)
    while Almanac.daffna(h)
        dc, ic = Almanac.dafus(Almanac.dafgs(h), nd, ni)
        data = ic[end] >= ic[end - 1] ? Almanac.dafgda(h, ic[end - 1], ic[end]) : Float64[]
        push!(out, (dc, ic, Almanac.dafgn(h), data))
    end
    return out
end

function compare_daf!(pairs, mine, theirs)
    hm, ht = Almanac.dafopr(mine), Almanac.dafopr(theirs)
    try
        # the summary format and the internal name; the record numbers and the free
        # address depend on the layout
        file_pair!(pairs, "nd,ni,ifname", Any[Almanac.dafrfr(ht)[1:3]...], Any[Almanac.dafrfr(hm)[1:3]...])
        file_pair!(pairs, "comments", Almanac.dafec(ht), Almanac.dafec(hm))
        at, am = daf_arrays(ht), daf_arrays(hm)
        file_pair!(pairs, "number of arrays", Any[length(at)], Any[length(am)])
        for (k, ((dt, it, nt, xt), (dm, im, nm, xm))) in enumerate(zip(at, am))
            file_pair!(pairs, "array $k name", Any[nt], Any[nm])
            file_pair!(pairs, "array $k dc", dt, dm)
            # the integers but the addresses, which depend on the layout
            file_pair!(pairs, "array $k ic", it[1:(end - 2)], im[1:(end - 2)])
            file_pair!(pairs, "array $k data", xt, xm)
        end
    finally
        Almanac.dafcls(hm); Almanac.dafcls(ht)
    end
    return pairs
end

function compare_das!(pairs, mine, theirs)
    hm, ht = Almanac.dasopr(mine), Almanac.dasopr(theirs)
    try
        file_pair!(pairs, "idword,ifname,nresvr,nresvc,ncomr,ncomc", Any[Almanac.dasrfr(ht)...], Any[Almanac.dasrfr(hm)...])
        file_pair!(pairs, "comments", Almanac.dasec(ht), Almanac.dasec(hm))
        lt, lm = Almanac.dashfs(ht)[6], Almanac.dashfs(hm)[6]
        file_pair!(pairs, "last addresses (c, d, i)", lt, lm)
        for (k, kind) in enumerate(("character", "double", "integer"))
            n = min(lt[k], lm[k])
            n > 0 || continue
            rd = (Almanac.dasrdc, Almanac.dasrdd, Almanac.dasrdi)[k]
            vt, vm = rd(ht, 1, n), rd(hm, 1, n)
            k == 1 && ((vt, vm) = (codeunits(vt), codeunits(vm)))
            file_pair!(pairs, "$kind data", vt, vm)
        end
    finally
        Almanac.dascls(hm); Almanac.dascls(ht)
    end
    return pairs
end

ek_attrs(x) = Any[x.cclass, String(x.dtype), x.strlen, x.size, x.indexd, x.nullok]

function compare_ek!(pairs, mine, theirs)
    hm, ht = Almanac.ekopr(mine), Almanac.ekopr(theirs)
    try
        rt, rm_ = Almanac.dasrfr(ht), Almanac.dasrfr(hm)
        file_pair!(pairs, "ifname,ncomc", Any[rt[2], rt[6]], Any[rm_[2], rm_[6]])
        file_pair!(pairs, "comments", Almanac.dasec(ht), Almanac.dasec(hm))
        nt, nm = Almanac.eknseg(ht), Almanac.eknseg(hm)
        file_pair!(pairs, "number of segments", Any[nt], Any[nm])
        for s in 1:min(nt, nm)
            st, sm = Almanac.ekssum(ht, s), Almanac.ekssum(hm, s)
            file_pair!(pairs, "segment $s table,nrows,ncols", Any[st.tabnam, st.nrows, st.ncols], Any[sm.tabnam, sm.nrows, sm.ncols])
            file_pair!(pairs, "segment $s columns", st.cnames, sm.cnames)
            for (k, (at, am)) in enumerate(zip(st.cdescrs, sm.cdescrs))
                file_pair!(pairs, "segment $s column $(st.cnames[k]) attributes", ek_attrs(at), ek_attrs(am))
            end
            (st.nrows == sm.nrows && st.cnames == sm.cnames) || continue
            for (k, col) in enumerate(st.cnames)
                dt = st.cdescrs[k].dtype
                rd = dt === :CHR ? Almanac.ekrcec : dt === :INT ? Almanac.ekrcei : Almanac.ekrced
                et, gt = Any[], Any[]
                for r in 1:st.nrows
                    (vt, nlt), (vm, nlm) = rd(ht, s, r, col), rd(hm, s, r, col)
                    push!(et, nlt, length(vt), vt...)
                    push!(gt, nlm, length(vm), vm...)
                end
                file_pair!(pairs, "segment $s column $col (null flag, size, values of each row)", et, gt)
            end
        end
    finally
        Almanac.dascls(hm); Almanac.dascls(ht)
    end
    return pairs
end

"The content of the file `mine` that Almanac wrote and of SPICELIB's copy `theirs`, as pairs."
function compare_files(mine, theirs; ifname = nothing)
    pairs = Tuple{String, Vector{Any}, Any}[]
    ft, fm = Almanac.getfat(theirs), Almanac.getfat(mine)
    file_pair!(pairs, "architecture,type", Any[ft...], Any[fm...])
    ft == fm || return pairs
    ft[1] == "DAF" && compare_daf!(pairs, mine, theirs)
    ft == ("DAS", "EK") && compare_ek!(pairs, mine, theirs)
    ft[1] == "DAS" && ft[2] != "EK" && compare_das!(pairs, mine, theirs)
    # the internal name `replay` that the open gave for SPICELIB's blank one (open_by_handle)
    ifname === nothing || replace!(
        p -> startswith(p[1], "nd,ni,ifname") || startswith(p[1], "idword,ifname") ||
            startswith(p[1], "ifname,") ? (p[1], p[2], Any[x == ifname ? "" : x for x in p[3]]) : p, pairs
    )
    return pairs
end

"""
A writer that the log does not record (a `W` mark): Almanac's copy of the file misses
its write. The files open for write that the call names (SPICELIB's handles of those are
negative) are given up; if it names none, it writes the current DAF (DAFWS, the SPKSnn
subset writers, whose handle is of the file they read), or else any file open for write.
"""
function unrecorded_writer(c::Call)
    hs = [a.values[1] for a in c.inputs if a.name in HANDLE_ARGS && a.kind == 'i' && length(a.values) == 1]
    w = filter(<(0), hs)
    if isempty(w)
        cur = current_daf_ids()
        w = isempty(cur) ? collect(WRITE_OPEN) : cur
    end
    foreach(h -> drop_handle(h, "a writer outside the log ($(lowercase(c.name)))"), w)
    return
end

"SPICELIB's handles of the current DAF of Almanac."
current_daf_ids() = CURRENT_DAF[] === nothing ? Int[] : [k for (k, v) in HMAP if v === CURRENT_DAF[]]

# the handles of the current call, for `counted`
const WRITE_HANDLE = Ref(Int[])

"Outputs that hold a value only if a count output is at least a number: SPICELIB leaves them undefined otherwise."
const COUNTED = Dict(
    "inelpl" => ("nxpts", ["xpt1" => 1, "xpt2" => 2]),
    "inrypl" => ("nxpts", ["xpt" => 1])
)

function counted(name, pairs)
    # for a file open for write Almanac's dashfs describes the segregated file that dascls
    # writes (documented); SPICELIB's describes the records as they are now
    name == "dashfs" && return filter(
        p -> !(p[1] in ("free", "lastrc", "lastwd")) ||
            !any(<(0), WRITE_HANDLE[]), pairs
    )
    rule = get(COUNTED, name, nothing)
    rule === nothing && return pairs
    k = findfirst(p -> p[1] == rule[1], pairs)
    k === nothing && return pairs
    n = pairs[k][2][1]
    return filter(p -> !any(r -> r.first == p[1] && n < r.second, rule[2]), pairs)
end

"""
Pairs with a handle output compared through `HMAP`: equal when Almanac returned the
handle (or path) that is mapped to SPICELIB's handle; an unmapped handle is mapped.
"""
function map_handles(pairs)
    return map(pairs) do (n, e, g)
        (n in HANDLE_ARGS && length(e) == 1 && g !== nothing && length(g) == 1) || return (n, e, g)
        h = get!(HMAP, e[1], g[1])
        return (n, e, h == g[1] ? e : g)
    end
end

# --- the CSPICE oracle ----------------------------------------------------------------

"Arguments as SPICE.jl takes them: plain arrays."
plain(x::StaticArray) = Array(x)
# an argument of mixed numbers from the log; SPICE.jl takes only concrete arrays
plain(x::AbstractArray{Any}) = all(v -> v isa Integer, x) ? Int.(x) :
    all(v -> v isa Real, x) ? Float64.(x) : all(v -> v isa AbstractString, x) ? String.(x) : x
plain(x) = x

"An input as SPICE.jl takes it: the planes and ellipses of SPICELIB (4 and 9 doubles) as its types."
function oracle_arg(a::Arg)
    if a.kind == 'd' && length(a.values) == 4 && a.name in ("plane", "projpl", "invpl")
        v = Float64.(a.values)
        return SPICE.Plane((v[1], v[2], v[3]), v[4])
    elseif a.kind == 'd' && length(a.values) == 9 && a.name in ("ellips", "ellipse", "limb")
        v = Float64.(a.values)
        return SPICE.Ellipse((v[1], v[2], v[3]), (v[4], v[5], v[6]), (v[7], v[8], v[9]))
    end
    return plain(jl(a))
end

"The routines that change the kernels or the settings, mirrored into CSPICE."
const STATEFUL = Set(
    [
        "furnsh", "unload", "kclear", "ldpool", "clpool", "pdpool", "pcpool",
        "pipool", "lmpool", "dvpool", "timdef", "spklef", "cklpf", "pcklof",
    ]
)

# SPICE.jl puts only the long message of CSPICE into its error; the short message, which
# is the kind, is added in front here.
ORACLE && Core.eval(
    SPICE, quote
        function handleerror()
            ccall((:failed_c, libcspice), Bool, ()) || return nothing
            buf = Array{UInt8}(undef, 1841)
            ccall((:getmsg_c, libcspice), Cvoid, (Cstring, Cint, Ref{UInt8}), "SHORT", 1841, buf)
            short = String(buf[1:(findfirst(iszero, buf) - 1)])
            ccall((:getmsg_c, libcspice), Cvoid, (Cstring, Cint, Ref{UInt8}), "LONG", 1841, buf)
            long = String(buf[1:(findfirst(iszero, buf) - 1)])
            ccall((:reset_c, libcspice), Cvoid, ())
            throw(SpiceError(short * " " * long))
        end
    end
)

"The kind of an error of SPICE.jl: the short message of CSPICE, which `handleerror` puts first."
function cspice_kind(e)
    m = match(r"SPICE\((\w+)\)", sprint(showerror, e))
    return m === nothing ? :unknown : Symbol(lowercase(m.captures[1]))
end

"Call the SPICE.jl function of the routine of `c` with its inputs; an exception is returned."
# Routines that open or write files are not called in CSPICE: an inconsistent file
# makes it stop the process.
const NO_ORACLE = r"^(ek|das|daf|dla|spkw|ckw|pckw|spko|cko|pcko|dsko|dskw|dskcls|spkcls|ckcls|pckcls)"

# CSPICE routines that SPICE.jl does not wrap (tangpt), or whose wrapper returns arrays of
# the full size `maxn` (limbpt, termpt) or takes `abcorr` as a keyword (ilumin,
# illumg), called as SPICE.jl calls the others.
const SpD = Float64
function cspice_tangpt(method, target, et, fixref, abcorr, corloc, obsrvr, dref, dvec)
    tanpt, srfpt, srfvec = zeros(3), zeros(3), zeros(3)
    alt, range, trgepc = Ref(0.0), Ref(0.0), Ref(0.0)
    ccall(
        (:tangpt_c, SPICE.libcspice), Cvoid,
        (
            Cstring, Cstring, SpD, Cstring, Cstring, Cstring, Cstring, Cstring, Ref{SpD},
            Ref{SpD}, Ref{SpD}, Ref{SpD}, Ref{SpD}, Ref{SpD}, Ref{SpD},
        ),
        method, target, et, fixref, abcorr, corloc, obsrvr, dref, Vector{SpD}(dvec),
        tanpt, alt, range, srfpt, trgepc, srfvec
    )
    SPICE.handleerror()
    return tanpt, alt[], range[], srfpt, trgepc[], srfvec
end
function cspice_termpt(
        method, ilusrc, target, et, fixref, abcorr, corloc, obsrvr, refvec,
        rolstp, ncuts, schstp, soltol, maxn
    )
    npts = zeros(Int32, ncuts)
    points, epochs, trmvcs = zeros(3, maxn), zeros(maxn), zeros(3, maxn)
    ccall(
        (:termpt_c, SPICE.libcspice), Cvoid,
        (
            Cstring, Cstring, Cstring, SpD, Cstring, Cstring, Cstring, Cstring, Ref{SpD}, SpD,
            Int32, SpD, SpD, Int32, Ref{Int32}, Ref{SpD}, Ref{SpD}, Ref{SpD},
        ),
        method, ilusrc, target, et, fixref, abcorr, corloc, obsrvr, Vector{SpD}(refvec),
        rolstp, ncuts, schstp, soltol, maxn, npts, points, epochs, trmvcs
    )
    SPICE.handleerror()
    n = Int(sum(npts))
    return Int.(npts), [points[:, i] for i in 1:n], epochs[1:n], [trmvcs[:, i] for i in 1:n]
end
function cspice_limbpt(
        method, target, et, fixref, abcorr, corloc, obsrvr, refvec,
        rolstp, ncuts, schstp, soltol, maxn
    )
    npts, points, epochs, tangts = SPICE.limbpt(
        method, target, et, fixref, abcorr, corloc, obsrvr, refvec,
        rolstp, ncuts, schstp, soltol, maxn
    )
    n = Int(sum(npts))
    return Int.(npts), points[1:n], epochs[1:n], tangts[1:n]
end
cspice_ilumin(method, target, et, fixref, abcorr, obsrvr, spoint) =
    SPICE.ilumin(method, target, et, fixref, obsrvr, spoint; abcorr)
cspice_illumg(method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint) =
    SPICE.illumg(method, target, ilusrc, et, fixref, obsrvr, spoint; abcorr)
const CSPICE_EXTRA = Dict{String, Function}(
    "tangpt" => cspice_tangpt, "termpt" => cspice_termpt,
    "limbpt" => cspice_limbpt, "ilumin" => cspice_ilumin, "illumg" => cspice_illumg,
)

# evsgp4 is not wrapped by SPICE.jl; its azlcpo wrapper takes method, azccw and elplsz as
# keywords, so the positional call of the oracle failed.
function cspice_evsgp4(et, geophs, elems)
    state = zeros(6)
    ccall(
        (:evsgp4_c, SPICE.libcspice), Cvoid, (SpD, Ref{SpD}, Ref{SpD}, Ref{SpD}),
        et, Vector{SpD}(geophs), Vector{SpD}(elems), state
    )
    SPICE.handleerror()
    return state
end
cspice_azlcpo(method, target, et, abcorr, azccw, elplsz, obspos, obsctr, obsref) =
    SPICE.azlcpo(target, et, abcorr, obspos, obsctr, obsref; method, azccw, elplsz)
# SPICE.jl's timdef has a buffer of 10 characters: a longer value (a zone) throws
function cspice_timdef(action, item, value)
    buf = zeros(UInt8, 256)
    v = Vector{UInt8}(value); copyto!(buf, v)
    ccall((:timdef_c, SPICE.libcspice), Cvoid, (Cstring, Cstring, Cint, Ptr{UInt8}), action, item, 256, buf)
    SPICE.handleerror()
    return String(buf[1:(findfirst(iszero, buf) - 1)])
end
CSPICE_EXTRA["timdef"] = cspice_timdef
CSPICE_EXTRA["evsgp4"] = cspice_evsgp4
CSPICE_EXTRA["azlcpo"] = cspice_azlcpo

# Routines whose SPICE.jl wrapper does not take the arguments as the log has them: the
# general matrix routines take Julia matrices of their own shape (the log has the column
# major array of the declared size and the dimensions), spkaps takes an acceleration of
# 6 numbers (SPICELIB's has 3), latsrf, srfnrm and dskrb2 take lists of points (the log
# has the flat array), and vequ takes an output vector. The matrix routines are called
# at the level of the f2c translation of SPICELIB (`mxvg_` ...), which takes the column
# major arrays as SPICELIB does; the others through the C interface.
const SpI = Cint
"An array of the log as a vector of doubles; an array of one value can come as a number."
_dv(x) = Float64.(vec(vcat(x)))
_cm(m, nr, nc) = reshape(_dv(m)[1:(nr * nc)], nr, nc)
function cspice_mxvg(m1, v2, nr1, nc1r2)
    vout = zeros(nr1)
    ccall(
        (:mxvg_, SPICE.libcspice), Cint, (Ref{SpD}, Ref{SpD}, Ref{SpI}, Ref{SpI}, Ref{SpD}),
        _cm(m1, nr1, nc1r2), _dv(v2)[1:nc1r2], nr1, nc1r2, vout
    )
    return vout
end
function cspice_mtxvg(m1, v2, nc1, nr1r2)
    vout = zeros(nc1)
    ccall(
        (:mtxvg_, SPICE.libcspice), Cint, (Ref{SpD}, Ref{SpD}, Ref{SpI}, Ref{SpI}, Ref{SpD}),
        _cm(m1, nr1r2, nc1), _dv(v2)[1:nr1r2], nc1, nr1r2, vout
    )
    return vout
end
function cspice_mxmg(m1, m2, nr1, nc1r2, nc2)
    mout = zeros(nr1, nc2)
    ccall(
        (:mxmg_, SPICE.libcspice), Cint, (Ref{SpD}, Ref{SpD}, Ref{SpI}, Ref{SpI}, Ref{SpI}, Ref{SpD}),
        _cm(m1, nr1, nc1r2), _cm(m2, nc1r2, nc2), nr1, nc1r2, nc2, mout
    )
    return mout
end
function cspice_mtxmg(m1, m2, nc1, nr1r2, nc2)
    mout = zeros(nc1, nc2)
    ccall(
        (:mtxmg_, SPICE.libcspice), Cint, (Ref{SpD}, Ref{SpD}, Ref{SpI}, Ref{SpI}, Ref{SpI}, Ref{SpD}),
        _cm(m1, nr1r2, nc1), _cm(m2, nr1r2, nc2), nc1, nr1r2, nc2, mout
    )
    return mout
end
function cspice_mxmtg(m1, m2, nr1, nc1c2, nr2)
    mout = zeros(nr1, nr2)
    ccall(
        (:mxmtg_, SPICE.libcspice), Cint, (Ref{SpD}, Ref{SpD}, Ref{SpI}, Ref{SpI}, Ref{SpI}, Ref{SpD}),
        _cm(m1, nr1, nc1c2), _cm(m2, nr2, nc1c2), nr1, nc1c2, nr2, mout
    )
    return mout
end
function cspice_xposeg(matrix, nrow, ncol)
    xposem = zeros(ncol, nrow)
    ccall(
        (:xposeg_, SPICE.libcspice), Cint, (Ref{SpD}, Ref{SpI}, Ref{SpI}, Ref{SpD}),
        _cm(matrix, nrow, ncol), nrow, ncol, xposem
    )
    return xposem
end
# the C interface takes the matrix in row major order: its transpose in Julia's order
function cspice_vtmvg(v1, matrix, v2, nrow, ncol)
    r = ccall(
        (:vtmvg_c, SPICE.libcspice), SpD, (Ref{SpD}, Ref{SpD}, Ref{SpD}, SpI, SpI),
        _dv(v1)[1:nrow], permutedims(_cm(matrix, nrow, ncol)),
        _dv(v2)[1:ncol], nrow, ncol
    )
    SPICE.handleerror()
    return r
end
function cspice_vequ(vin)
    vout = zeros(3)
    ccall((:vequ_c, SPICE.libcspice), Cvoid, (Ref{SpD}, Ref{SpD}), _dv(vin), vout)
    return vout
end
function cspice_spkaps(targ, et, ref, abcorr, stobs, accobs)
    starg = zeros(6)
    lt, dlt = Ref(0.0), Ref(0.0)
    ccall(
        (:spkaps_c, SPICE.libcspice), Cvoid,
        (SpI, SpD, Cstring, Cstring, Ref{SpD}, Ref{SpD}, Ref{SpD}, Ref{SpD}, Ref{SpD}),
        targ, et, ref, abcorr, _dv(stobs), _dv(accobs), starg, lt, dlt
    )
    SPICE.handleerror()
    return starg, lt[], dlt[]
end
# the points are rows of a C array, columns of a Julia matrix
function cspice_latsrf(method, target, et, fixref, lonlat)
    ll = reshape(_dv(lonlat), 2, :)
    npts = size(ll, 2)
    srfpts = zeros(3, npts)
    ccall(
        (:latsrf_c, SPICE.libcspice), Cvoid,
        (Cstring, Cstring, SpD, Cstring, SpI, Ref{SpD}, Ref{SpD}),
        method, target, et, fixref, npts, ll, srfpts
    )
    SPICE.handleerror()
    return srfpts
end
function cspice_srfnrm(method, target, et, fixref, srfpts)
    p = reshape(_dv(srfpts), 3, :)
    npts = size(p, 2)
    normls = zeros(3, npts)
    ccall(
        (:srfnrm_c, SPICE.libcspice), Cvoid,
        (Cstring, Cstring, SpD, Cstring, SpI, Ref{SpD}, Ref{SpD}),
        method, target, et, fixref, npts, p, normls
    )
    SPICE.handleerror()
    return normls
end
function cspice_dskrb2(vrtces, plates, corsys, corpar)
    v = reshape(_dv(vrtces), 3, :)
    p = reshape(SpI.(vec(vcat(plates))), 3, :)
    mncor3, mxcor3 = Ref(0.0), Ref(0.0)
    ccall(
        (:dskrb2_c, SPICE.libcspice), Cvoid,
        (SpI, Ref{SpD}, SpI, Ref{SpI}, SpI, Ref{SpD}, Ref{SpD}, Ref{SpD}),
        size(v, 2), v, size(p, 2), p, corsys, _dv(corpar), mncor3, mxcor3
    )
    SPICE.handleerror()
    return mncor3[], mxcor3[]
end
function cspice_vequg(vin)
    v = _dv(vin)
    vout = zeros(length(v))
    ccall((:vequg_c, SPICE.libcspice), Cvoid, (Ref{SpD}, SpI, Ref{SpD}), v, length(v), vout)
    return vout
end
# SPICE.jl's polyds takes no degree (the log has the declared size of the coefficients)
# a negative degree has no coefficients in the log; CSPICE reads none of them then
function cspice_polyds(coeffs, deg, nderiv, t)
    p = zeros(max(nderiv, 0) + 1)
    c = vcat(_dv(coeffs), zeros(max(deg + 1, 1)))[1:max(deg + 1, 1)]
    ccall(
        (:polyds_c, SPICE.libcspice), Cvoid, (Ref{SpD}, SpI, SpI, SpD, Ref{SpD}),
        c, deg, nderiv, t, p
    )
    SPICE.handleerror()
    return p
end
# the counts nv and np are dropped: the arrays give them
_plates(vrtces, plates) = reshape(_dv(vrtces), 3, :), reshape(SpI.(vec(vcat(plates))), 3, :)
function cspice_pltar(vrtces, plates)
    v, p = _plates(vrtces, plates)
    r = ccall(
        (:pltar_c, SPICE.libcspice), SpD, (SpI, Ref{SpD}, SpI, Ref{SpI}),
        size(v, 2), v, size(p, 2), p
    )
    SPICE.handleerror()
    return r
end
function cspice_pltvol(vrtces, plates)
    v, p = _plates(vrtces, plates)
    r = ccall(
        (:pltvol_c, SPICE.libcspice), SpD, (SpI, Ref{SpD}, SpI, Ref{SpI}),
        size(v, 2), v, size(p, 2), p
    )
    SPICE.handleerror()
    return r
end
CSPICE_EXTRA["pltar"] = cspice_pltar
CSPICE_EXTRA["pltvol"] = cspice_pltvol
# swpool has a count of names; the watchers it sets decide what cvpool returns, so it
# is mirrored into CSPICE as the loads are
cspice_swpool(agent, nnames, names) =
    SPICE.swpool(agent, String.(vcat(names))[1:max(nnames, 0)])
CSPICE_EXTRA["swpool"] = cspice_swpool
push!(STATEFUL, "swpool")
# the azimuth and elevation conversions of SPICE.jl take azccw and elplsz as keywords
cspice_azlrec(range, az, el, azccw, elplsz) = SPICE.azlrec(range, az, el; azccw, elplsz)
cspice_drdazl(range, az, el, azccw, elplsz) = SPICE.drdazl(range, az, el; azccw, elplsz)
cspice_recazl(rectan, azccw, elplsz) = SPICE.recazl(rectan; azccw, elplsz)
cspice_dazldr(x, y, z, azccw, elplsz) = SPICE.dazldr(x, y, z; azccw, elplsz)
for (name, f) in (
        "azlrec" => cspice_azlrec, "drdazl" => cspice_drdazl, "recazl" => cspice_recazl,
        "dazldr" => cspice_dazldr, "vequg" => cspice_vequg, "polyds" => cspice_polyds,
        "mxvg" => cspice_mxvg, "mtxvg" => cspice_mtxvg, "mxmg" => cspice_mxmg,
        "mtxmg" => cspice_mtxmg, "mxmtg" => cspice_mxmtg, "xposeg" => cspice_xposeg,
        "vtmvg" => cspice_vtmvg, "vequ" => cspice_vequ, "spkaps" => cspice_spkaps,
        "latsrf" => cspice_latsrf, "srfnrm" => cspice_srfnrm, "dskrb2" => cspice_dskrb2,
    )
    CSPICE_EXTRA[name] = f
end
# body names that a family defines with boddef (F_SINCPT maps SUN to 399) reach CSPICE
push!(STATEFUL, "boddef")

"The answer of `cspice_call` for a call that CSPICE is not asked (`nothing` is the result of a wrapper that found nothing)."
struct NoOracle end

function cspice_call(c::Call)
    name = lowercase(c.name)
    occursin(NO_ORACLE, name) && return NoOracle()
    any(a -> handle_value(a) !== nothing, c.inputs) && return NoOracle()
    name in ("spklef", "cklpf", "pcklof") && (name = "furnsh")
    f = haskey(CSPICE_EXTRA, name) ? CSPICE_EXTRA[name] :
        isdefined(SPICE, Symbol(name)) ? getfield(SPICE, Symbol(name)) : nothing
    f === nothing && return NoOracle()
    drop = get(DROP, lowercase(c.name), String[])
    args = Any[oracle_arg(a) for a in c.inputs if !(a.name in drop)]
    v = isempty(c.inputs) ? Any[] : c.inputs[end].values
    name == "pdpool" && (args = Any[args[1], Float64.(v)])
    name == "pipool" && (args = Any[args[1], Int.(v)])
    name == "pcpool" && (args = Any[args[1], String.(v)])
    name == "lmpool" && (args = Any[String.(c.inputs[1].values)])
    try
        r = in_snapshot(() -> f(args...), c)
        # some wrappers of SPICE.jl (ltime, polyds, chbder, szpool, ...) do not check for
        # an error: it would return garbage, and leave CSPICE in the error state, in
        # which the next call, a load mirrored into CSPICE among them, does nothing
        SPICE.handleerror()
        return r
    catch e
        # an exception that is not an error of CSPICE (a MethodError of a call that does
        # not fit SPICE.jl's signature) means the oracle did not answer: it is counted,
        # so that a mapping fault does not pass for a disagreement with CSPICE
        e isa SPICE.SpiceError || note_oracle_failure(name, e)
        return e
    end
end

"Calls for which the oracle failed, per routine: count and first message (oracle_failures.tsv)."
const ORACLE_FAILURES = Dict{String, Tuple{Int, String}}()

function note_oracle_failure(name, e)
    n, msg = get(ORACLE_FAILURES, name, (0, ""))
    ORACLE_FAILURES[name] = (n + 1, n == 0 ? replace(first(sprint(showerror, e), 200), r"\s+" => " ") : msg)
    return
end

function cspice_reset()
    ccall((:reset_c, SPICE.libcspice), Cvoid, ())
    SPICE.kclear()
    SPICE.timdef("SET", "CALENDAR", "GREGORIAN")
    return SPICE.timdef("SET", "SYSTEM", "UTC")
end

"""
True if CSPICE gives for `c` what Almanac gave and SPICELIB did not: the same error kind,
or outputs that agree within `NEAR` with `almanac_result`, the outputs of Almanac paired
with SPICELIB's, where CSPICE also differs from SPICELIB in each output that Almanac
differs in. Without that last condition a fault of Almanac below `NEAR` would pass as
`cspice` wherever CSPICE and SPICELIB agree. For a call that SPICELIB failed and CSPICE
did not, SPICELIB's outputs are undefined: Almanac's outputs are compared with CSPICE's.
"""
function agrees_with_cspice(c::Call, almanac_result)
    r = cspice_call(c)
    r isa NoOracle && return false
    if almanac_result isa SpiceError
        return r isa Exception && cspice_kind(r) === almanac_result.kind
    end
    r isa Exception && return false
    # SPICE.jl returns `nothing` where CSPICE found nothing; Almanac's result is then
    # the output `found` alone, and it must be false
    r === nothing && return length(almanac_result) == 1 &&
        almanac_result[1][1] == "found" && almanac_result[1][3] == Any[false]
    name = lowercase(c.name)
    notout = get(NOTOUT, name, String[])
    outs = [o for o in c.outputs if !(o.name in notout && o.name != "found")]
    filter!(o -> o.name != "found", outs)
    vals = r isa Tuple ? collect(r) : Any[r]
    oracle = output_pairs(name, outs, vals)
    length(oracle) == length(almanac_result) || return false
    for ((_, e, a), (_, _, o)) in zip(almanac_result, oracle)
        if a === nothing || o === nothing
            a === o || return false
            continue
        end
        length(a) == length(o) || return false
        near_output(a, o) || return false
        # SPICELIB signalled an error and CSPICE did not: SPICELIB's outputs are the
        # values the caller left there, so Almanac's are compared with CSPICE's alone
        c.failed && continue
        length(e) == length(a) && all(close_enough.(e, a)) && continue
        length(e) == length(o) && all(close_enough.(e, o)) && return false
    end
    return true
end

"""
True if the outputs `a` and `b` agree within `NEAR`, relative to the largest magnitude
(at least 1) in each group of three numbers: a 3-vector is compared as a whole, so a
small component of a large vector does not need a tighter agreement than the vector
has, and the position and the velocity of a state are compared apart.
"""
function near_output(a, b)
    length(a) == length(b) || return false
    k = length(a) % 3 == 0 ? 3 : max(length(a), 1)
    for i in 1:k:length(a)
        r = i:min(i + k - 1, length(a))
        scale = 1.0
        for x in Iterators.flatten((a[r], b[r]))
            x isa Real && isfinite(x) && (scale = max(scale, abs(Float64(x))))
        end
        all(
            t -> close_enough(t[1], t[2]) ||
                (t[1] isa Real && t[2] isa Real && abs(t[1] - t[2]) <= NEAR * scale), zip(a[r], b[r])
        ) ||
            return false
    end
    return true
end

"""
True if SPICELIB and CSPICE disagree with each other, and Almanac is no farther from
CSPICE than twice the distance of SPICELIB, in every output: the result is within the
scatter of the references. An ill-conditioned geometry gives this: a ray that grazes a
plate from 2e8 km moves its intercept by 1e-5 km for one unit in the last place of the
rotation to the body frame, whose sine and cosine differ in the last bit between Julia
and the C library.
"""
function within_spread(c::Call, pairs)
    r = cspice_call(c)
    (r isa NoOracle || r === nothing || r isa Exception) && return false
    name = lowercase(c.name)
    notout = get(NOTOUT, name, String[])
    outs = filter(o -> !(o.name in notout) && o.name != "found", c.outputs)
    oracle = output_pairs(name, outs, r isa Tuple ? collect(r) : Any[r])
    length(oracle) == length(pairs) || return false
    for ((_, e, g), (_, _, o)) in zip(pairs, oracle)
        (g === nothing || o === nothing || length(g) != length(o) || length(e) != length(o)) && return false
        all(close_enough.(e, g)) && continue
        # an output where CSPICE agrees with SPICELIB has no scatter to be within
        all(close_enough.(e, o)) && return false
        # integers (handles, counts, indices) have no rounding scatter: they must agree
        if all(x -> x isa AbstractFloat, Iterators.flatten((e, g, o)))
            dac = maximum(abs.(Float64.(g) .- Float64.(o)); init = 0.0)
            dsc = maximum(abs.(Float64.(e) .- Float64.(o)); init = 0.0)
            dac <= 2dsc || near_output(g, o) || return false
        else
            near_output(g, o) || return false
        end
    end
    return true
end

function output_pairs(name, outs, vals)
    pairs = Tuple{String, Vector{Any}, Any}[]
    if length(vals) != length(outs)
        # a single array output that SPICELIB splits into several, or the reverse
        flatvals = reduce(vcat, [flat(v) for v in vals if v !== nothing]; init = Any[])
        exp = reduce(vcat, [o.values for o in outs]; init = Any[])
        push!(pairs, (join([o.name for o in outs], "+"), exp, flatvals))
        return pairs
    end
    for (o, v) in zip(outs, vals)
        push!(pairs, (o.name, o.values, flat(v)))
    end
    return pairs
end

"SPICELIB's handles of the files loaded by handle in this family, and their paths."
const HANDLES = Dict{Int, String}()

function load_by_handle(ctx, c)
    path = jl(c.inputs[1])
    in_snapshot(c) do
        furnsh(ctx, path)
    end
    # Almanac keeps the path as given, and `spksfs` returns it as the handle.
    c.failed || (
        HANDLES[c.outputs[1].values[1]] = path; HMAP[c.outputs[1].values[1]] = path;
        LOADDIR[path] = isempty(c.snapshot) ? pwd() : c.snapshot
    )
    return []
end

function unload_by_handle(ctx, c)
    path = get(HANDLES, jl(c.inputs[1]), nothing)
    path === nothing && return []      # not loaded by handle (F_SGMETA unloads its write handles)
    unload(ctx, path)
    delete!(HMAP, jl(c.inputs[1]))
    return []
end

# Routines whose arguments do not follow the default.
const MAPPINGS = Dict{String, Function}(
    "spklef" => load_by_handle, "cklpf" => load_by_handle, "pcklof" => load_by_handle,
    "eklef" => load_by_handle, "ekuef" => unload_by_handle, "dasllc" => dasllc_mapping,
    "dasops" => (ctx, c) -> (h = dasops(); c.failed || (HMAP[c.outputs[1].values[1]] = h); []),
    # SPICELIB keeps the type 14 segment being written for its handle; Almanac returns
    # a writer from spk14b and takes it in spk14a and spk14e
    "spk14b" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        h = handle_value(a["handle"]); h === nothing && throw(UnmappedHandle())
        SPK14W[a["handle"].values[1]] = spk14b(
            ctx, h, jl(a["segid"]), jl(a["body"]), jl(a["center"]),
            jl(a["frame"]), jl(a["first"]), jl(a["last"]), jl(a["chbdeg"])
        )
        return []
    end,
    "spk14a" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        w = get(SPK14W, a["handle"].values[1], nothing); w === nothing && throw(UnmappedHandle())
        spk14a(w, jl(a["ncsets"]), Float64.(a["coeffs"].values), Float64.(a["epochs"].values))
        return []
    end,
    "spk14e" => function (ctx, c)
        # a file given up (SPK14A is not recorded) is closed in Almanac
        haskey(HMAP, c.inputs[1].values[1]) || throw(UnmappedHandle())
        w = get(SPK14W, c.inputs[1].values[1], nothing); w === nothing && throw(UnmappedHandle())
        spk14e(w)
        return []
    end,
    "xposeg" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        nr, nc = jl(a["nrow"]), jl(a["ncol"])
        r = xposeg(reshape(Float64.(a["matrix"].values)[1:(nr * nc)], nr, nc))
        return [("xposem", c.outputs[1].values, flat(r))]
    end,
    "swpool" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        names = String.(a["names"].values)[1:max(0, jl(a["nnames"]))]
        swpool(ctx, jl(a["agent"]), names)
        return []
    end,
    "badkpv" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        r = badkpv(
            ctx, jl(a["caller"]), jl(a["name"]), jl(a["comp"]), jl(a["size"]),
            jl(a["divby"]), only(jl(a["type"]))
        )
        return output_pairs("badkpv", c.outputs, Any[r])
    end,
    "ekops" => (ctx, c) -> (h = ekops(); c.failed || (HMAP[c.outputs[1].values[1]] = h); []),
    "ekpsel" => function (ctx, c)
        # SPICELIB: the item types and classes as words; byte positions from 1 in both
        r = ekpsel(ctx, jl(c.inputs[1]))
        out = Dict(o.name => o for o in c.outputs)
        pairs = Any[("error", out["error"].values, Any[r.failed])]
        r.failed && return pairs
        n = r.n
        push!(pairs, ("n", out["n"].values, Any[n]))
        for (k, v) in (("xbegs", r.xbegs), ("xends", r.xends))
            push!(pairs, (k, out[k].values[1:min(n, end)], Any[v...]))
        end
        for (k, v) in (("xtypes", r.xtypes), ("xclass", r.xclass), ("tabs", r.tabs), ("cols", r.cols))
            push!(pairs, (k, Any[rstrip(String(x)) for x in out[k].values[1:min(n, end)]], Any[string(x) for x in v]))
        end
        return pairs
    end,
    "spkw19" => (ctx, c) -> minisegment_writer(ctx, c, spkw19, "nintvl", ("packts", "epochs", "ivlbds")),
    "ckw06" => (ctx, c) -> minisegment_writer(ctx, c, ckw06, "nmini", ("packts", "rates", "sclkdp", "ivlbds")),
    # the EK fast load: Almanac's EKIFLD returns the state of the load, which the column
    # writers and EKFFLD take where SPICELIB takes the record pointers
    "ekifld" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        h = handle_value(a["handle"]); h === nothing && throw(UnmappedHandle())
        n = max(0, jl(a["ncols"]))
        segno, rc = ekifld(
            h, jl(a["tabnam"]), jl(a["nrows"]), String.(a["cnames"].values)[1:n],
            String.(a["decls"].values)[1:n]
        )
        EKFAST[(a["handle"].values[1], segno)] = rc
        out = Dict(o.name => o for o in c.outputs)
        return [("segno", out["segno"].values, Any[segno])]
    end,
    "ekaclc" => (ctx, c) -> ek_column(c, ekaclc, String),
    "ekacld" => (ctx, c) -> ek_column(c, ekacld, Float64),
    "ekacli" => (ctx, c) -> ek_column(c, ekacli, Int),
    "ekffld" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        h = handle_value(a["handle"]); h === nothing && throw(UnmappedHandle())
        rc = get(EKFAST, (a["handle"].values[1], jl(a["segno"])), nothing)
        rc === nothing && throw(UnmappedHandle())
        ekffld(h, jl(a["segno"]), rc)
        return []
    end,
    # DAFWS and DAFRS write the summary of the current array of the search
    "dafws" => (ctx, c) -> daf_summary_writer(c, dafws),
    "dafrs" => (ctx, c) -> daf_summary_writer(c, dafrs),
    # the workspace is not logged; SPICELIB's integer component has the declared size
    # SPXISZ, of which the index uses the length that Almanac returns
    "dskmi2" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        any(x -> any(v -> v isa Symbol, x.values), values(a)) && throw(UnmappedHandle())
        nv, np = jl(a["nv"]), jl(a["np"])
        (nv < 0 || np < 0) && throw(UnmappedHandle())
        v = reshape(Float64.(a["vrtces"].values)[1:(3nv)], 3, nv)
        p = reshape(Int.(a["plates"].values)[1:(3np)], 3, np)
        d, i = dskmi2(
            v, p, jl(a["finscl"]), jl(a["corscl"]), jl(a["worksz"]), jl(a["voxpsz"]),
            jl(a["voxlsz"]), jl(a["makvtl"]), jl(a["spxisz"])
        )
        out = Dict(o.name => o for o in c.outputs)
        ei = out["spaixi"].values
        return [("spaixd", out["spaixd"].values, flat(d)), ("spaixi", ei[1:min(end, length(i))], flat(i))]
    end,
    "dafcs" => default_mapping, "dafgh" => default_mapping,
    # DAFCAD selects the DAF to which DAFADA and DAFENA add (Almanac takes the handle)
    "dafcad" => default_mapping,
    # the log holds the summary unpacked (the integers packed into doubles would not
    # survive the decimal log)
    "dafbna" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        any(x -> any(v -> v isa Symbol, x.values), values(a)) && throw(UnmappedHandle())
        h = handle_value(a["handle"])
        CURRENT_DAF[] = h
        h === nothing && throw(UnmappedHandle())
        dafbna(h, Almanac.dafps(Float64.(a["dc"].values), Int.(a["ic"].values)), jl(a["name"]))
        return []
    end,
    "spkuef" => unload_by_handle, "ckupf" => unload_by_handle, "pckuof" => unload_by_handle,
    "tparse" => function (ctx, c)
        errmsg = c.outputs[2].values[1]
        r = try
            tparse(ctx, jl(c.inputs[1]))
        catch e
            e isa SpiceError || rethrow()
            e
        end
        isempty(errmsg) || return [("errmsg", Any["error: " * errmsg], Any[r isa SpiceError ? "error: " * errmsg : string(r)])]
        r isa SpiceError && throw(r)
        return [("sp2000", c.outputs[1].values, Any[r])]
    end,
    "gdpool" => (ctx, c) -> pool_get(ctx, c, gdpool),
    "gipool" => (ctx, c) -> pool_get(ctx, c, gipool),
    "gcpool" => (ctx, c) -> pool_get(ctx, c, gcpool),
    "dtpool" => function (ctx, c)
        found, n, typ = c.outputs[1].values[1], c.outputs[2].values[1], c.outputs[3].values[1]
        r = dtpool(ctx, jl(c.inputs[1]))
        r[2] === :X && return [("found", Any[found], Any[false])]   # Almanac: (0, :X) when absent
        found || return [("found", Any[false], Any[true])]
        return [("n,type", Any[n, typ], Any[r[1], string(r[2])])]
    end,
    "furnsh" => function (ctx, c)
        path = jl(c.inputs[1])
        in_snapshot(() -> furnsh(ctx, path), c)
        # a later call that takes the handle of the file (kinfo) reads it where it was loaded
        isempty(c.snapshot) || (LOADDIR[path] = c.snapshot; LAST_LOAD[] = c.snapshot)
        return []
    end,
    "ldpool" => (ctx, c) -> (
        in_snapshot(c) do;
            ldpool(ctx, jl(c.inputs[1]))
        end; []
    ),
    "lmpool" => (ctx, c) -> (lmpool(ctx, String.(c.inputs[1].values)); []),
    "pdpool" => (ctx, c) -> (pdpool(ctx, jl(c.inputs[1]), Float64.(c.inputs[3].values)); []),
    "pipool" => (ctx, c) -> (pipool(ctx, jl(c.inputs[1]), Int.(c.inputs[3].values)); []),
    "pcpool" => (ctx, c) -> (pcpool(ctx, jl(c.inputs[1]), String.(c.inputs[3].values)); []),
    "mxvg" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        nr, nc = a["nr1"].values[1], a["nc1r2"].values[1]
        r = mxvg(gmat(a["m1"], nr, nc), Float64.(a["v2"].values))
        return output_pairs("mxvg", c.outputs, Any[r])
    end,
    "mtxvg" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        nc, nr = a["nc1"].values[1], a["nr1r2"].values[1]
        r = mtxvg(gmat(a["m1"], nr, nc), Float64.(a["v2"].values))
        return output_pairs("mtxvg", c.outputs, Any[r])
    end,
    # the logged arrays have the declared size; the matrices are their first values
    "polyds" => function (ctx, c)
        # COEFFS(0:DEG): also a vector when it holds one number
        a = Dict(x.name => x for x in c.inputs)
        r = polyds(Float64.(a["coeffs"].values)[1:(jl(a["deg"]) + 1)], jl(a["nderiv"]), jl(a["t"]))
        return output_pairs("polyds", c.outputs, Any[r])
    end,
    "mtxmg" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        nc1, nr, nc2 = jl(a["nc1"]), jl(a["nr1r2"]), jl(a["nc2"])
        r = mtxmg(gmatn(a["m1"], nr, nc1), gmatn(a["m2"], nr, nc2))
        return output_pairs("mtxmg", c.outputs, Any[r])
    end,
    "mxmtg" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        nr1, nc, nr2 = jl(a["nr1"]), jl(a["nc1c2"]), jl(a["nr2"])
        r = mxmtg(gmatn(a["m1"], nr1, nc), gmatn(a["m2"], nr2, nc))
        return output_pairs("mxmtg", c.outputs, Any[r])
    end,
    "vtmvg" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        nr, nc = jl(a["nrow"]), jl(a["ncol"])
        r = vtmvg(Float64.(a["v1"].values)[1:nr], gmatn(a["matrix"], nr, nc), Float64.(a["v2"].values)[1:nc])
        return output_pairs("vtmvg", c.outputs, Any[r])
    end,
    "dasac" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        h = handle_value(a["handle"]); h === nothing && throw(UnmappedHandle())
        dasac(h, String.(a["buffer"].values)[1:max(0, jl(a["n"]))])
        return []
    end,
    "mxmg" => function (ctx, c)
        a = Dict(x.name => x for x in c.inputs)
        nr1, nc1r2, nc2 = a["nr1"].values[1], a["nc1r2"].values[1], a["nc2"].values[1]
        r = mxmg(gmat(a["m1"], nr1, nc1r2), gmat(a["m2"], nc1r2, nc2))
        return output_pairs("mxmg", c.outputs, Any[r])
    end,
    "gnpool" => function (ctx, c)
        # Almanac returns the names in ascending order, CSPICE in the order of its hash
        # table: compare the sets when SPICELIB returned all of them.
        a = Dict(x.name => jl(x) for x in c.inputs)
        out = Dict(o.name => o.values for o in c.outputs)
        # the checks of the arguments (room), then all the names
        gnpool(ctx, a["name"], 1, a["room"])
        r = gnpool(ctx, a["name"])
        found, n = out["found"][1], out["n"][1]
        isempty(r) && return [("found", Any[found], Any[false])]
        found || return [("found", Any[false], Any[true])]
        # past the start: the count of the names that are left, up to the room
        (a["start"] == 1 && n < a["room"]) ||
            return [("n", Any[n], Any[clamp(length(r) - max(a["start"], 1) + 1, 0, a["room"])])]
        return [("cvals", sort(out["cvals"][1:n]), Any[sort(r)...])]
    end,
    "ktotal" => (ctx, c) -> [("count", c.outputs[1].values, Any[ktotal(ctx, jl(c.inputs[1]))])],
    "et2lst" => function (ctx, c)
        a = Dict(x.name => jl(x) for x in c.inputs)
        r = et2lst(ctx, a["et"], a["body"], a["lon"], a["type"])
        return output_pairs("et2lst", c.outputs, collect(r))
    end,
)

function pool_get(ctx, c, f)
    a = Dict(x.name => jl(x) for x in c.inputs)
    out = Dict(o.name => o.values for o in c.outputs)
    r = f(ctx, a["name"])
    found = out["found"][1]
    r === nothing && return [("found", Any[found], Any[false])]
    found || return [("found", Any[false], Any[true])]
    n = out["n"][1]
    got = Any[r[a["start"]:min(end, a["start"] + n - 1)]...]
    key = c.name == "GDPOOL" ? "values" : c.name == "GIPOOL" ? "ivals" : "cvals"
    return [(key, out[key][1:n], got)]
end

"Run `f` in the copy of the working directory that the call names (or where we are)."
in_snapshot(f, c::Call) = isempty(c.snapshot) ? f() : cd(f, c.snapshot)

# --- replay ------------------------------------------------------------------------------

# TSPICE_AUDIT=1 writes audit.tsv: per routine, the matches, those that compared no
# value, those that compared fewer values than SPICELIB's outputs hold (`found` and the
# outputs of NOTOUT aside), and the outputs that those left out. A routine with many of
# either is a routine whose matches say little.
const AUDIT = get(ENV, "TSPICE_AUDIT", "0") == "1"
const AUDITS = Dict{String, Vector{Int}}()

const AUDIT_MISSING = Dict{String, Dict{String, Int}}()

function audit!(name, c, result)
    ncmp = sum((length(e) for (_, e, _) in result); init = 0)
    notout = get(NOTOUT, name, String[])
    outs = filter(o -> o.name != "found" && !(o.name in notout), c.outputs)
    nlog = sum((length(o.values) for o in outs); init = 0)
    a = get!(() -> [0, 0, 0], AUDITS, name)
    a[1] += 1
    ncmp == 0 && (a[2] += 1)
    if 0 < ncmp < nlog
        a[3] += 1
        seen = Set(Iterators.flatten(split(n, '+') for (n, _, _) in result))
        m = get!(() -> Dict{String, Int}(), AUDIT_MISSING, name)
        for o in outs
            o.name in seen || (m[o.name] = get(m, o.name, 0) + 1)
        end
    end
    return
end

# TSPICE_ORACLE_ALL=1 (with TSPICE_ORACLE=1) calls CSPICE for every call that returns
# values, not only for the mismatches, and writes cspice_agreement.tsv: per routine, how
# many calls of Almanac and of SPICELIB agree with CSPICE exactly, within `RTOL`, within
# `NEAR`, or not. It measures how close each is to CSPICE; it does not change the counts.
const ORACLE_ALL = ORACLE && get(ENV, "TSPICE_ORACLE_ALL", "0") == "1"
const AGREEMENT = Dict{String, Vector{Int}}()     # n, Almanac exact/rtol/near/far, SPICELIB the same

function _agreement_level(x, o)
    (x === nothing || o === nothing) && return x === o ? 1 : 4
    length(x) == length(o) || return 4
    all(t -> t[1] == t[2] || (t[1] isa AbstractString && t[2] isa AbstractString && rstrip(t[1]) == rstrip(t[2])), zip(x, o)) && return 1
    all(close_enough.(x, o)) && return 2
    near_output(x, o) && return 3
    return 4
end

function compare_all!(name, c, result)
    r = cspice_call(c)
    (r isa NoOracle || r === nothing || r isa Exception) && return
    notout = get(NOTOUT, name, String[])
    outs = filter(o -> o.name != "found" && !(o.name in notout), c.outputs)
    oracle = output_pairs(name, outs, r isa Tuple ? collect(r) : Any[r])
    length(oracle) == length(result) || return
    la = lb = 1
    for ((_, e, g), (_, _, o)) in zip(result, oracle)
        la = max(la, _agreement_level(g, o))
        lb = max(lb, _agreement_level(e, o))
    end
    v = get!(() -> zeros(Int, 9), AGREEMENT, name)
    v[1] += 1
    v[1 + la] += 1
    v[5 + lb] += 1
    return
end

# The reasons of the skipped calls, per routine, for skips.tsv.
const SKIP_REASONS = ("nested", "error_state", "truncated", "after_truncated_load", "no_mapping", "unmapped_handle")
const SKIPS = Dict{String, Dict{String, Int}}()
note_skip(name, why) = (d = get!(() -> Dict{String, Int}(), SKIPS, name); d[why] = get(d, why, 0) + 1; nothing)

# For selftest.jl: TSPICE_FAULT=value scales each floating-point output of Almanac by
# 1 + 1e-10, below `NEAR`; TSPICE_FAULT=kind changes the kind of each SpiceError.
const FAULT = get(ENV, "TSPICE_FAULT", "")

function inject_fault(result)
    if FAULT == "value" && result isa AbstractVector
        scale(g) = g === nothing ? g : Any[x isa AbstractFloat ? x * (1 + 1.0e-10) : x for x in g]
        return [(n, e, scale(g)) for (n, e, g) in result]
    elseif FAULT == "kind" && result isa SpiceError
        return SpiceError(result.kind === :invalidoption ? :badframe : :invalidoption, "injected")
    end
    return result
end

mutable struct Tally
    match::Int
    cspice::Int
    near::Int
    spread::Int
    value::Int
    error::Int
    adapter::Int
    skipped::Int
end
Tally() = Tally(0, 0, 0, 0, 0, 0, 0, 0)

function replay_family(path, tallies, mismatches)
    family = replace(basename(path), r"\.log(\.gz)?$" => "")
    ctx = SpiceContext()
    empty!(HANDLES)
    empty!(HMAP)
    CURRENT_DAF[] = nothing; CURRENT_EK[] = nothing
    empty!(HPATH); empty!(TAINTED); empty!(WRITE_OPEN); empty!(SPK14W); empty!(IFNAME_REPLAY); empty!(DROPPED); empty!(EKFAST)
    foreach(
        h -> try
            Almanac.dafcls(h)
        catch end, values(READERS)
    ); empty!(READERS)
    empty!(LOADDIR); LAST_LOAD[] = ""
    cd(mktempdir())
    ORACLE && cspice_reset()
    # A load or pool edit whose values the log does not hold leaves Almanac without
    # the kernels of SPICELIB: the rest of the family is skipped.
    desync = false
    return foreach_call(path; want = n -> haskey(MAPPINGS, n) || isdefined(Almanac, Symbol(n))) do c
        name = lowercase(c.name)
        t = get!(Tally, tallies, name)
        truncated = any(a -> a.values == Any[TRUNCATED], c.inputs) ||
            any(a -> a.values == Any[TRUNCATED], c.outputs)
        truncated && name in STATEFUL && (desync = true)
        # a write that is not replayed leaves the file of Almanac incomplete: the later
        # calls on its handle are skipped
        for a in c.inputs
            (truncated && a.name in HANDLE_ARGS && a.kind == 'i') || continue
            drop_handle(a.values[1], "a call with more values than the log keeps ($name)")
        end
        # the same for a write to the current DAF (DAFADA of more values than the log keeps)
        (truncated && name in DAF_CURRENT && c.depth == 1 && !c.failed_before) &&
            foreach(h -> drop_handle(h, "a call with more values than the log keeps ($name)"), current_daf_ids())
        if c.unrecorded
            c.depth == 1 && !c.failed_before && unrecorded_writer(c)
            t.skipped += 1
            return
        end
        if desync || c.depth > 1 || c.failed_before || truncated ||
                !(haskey(MAPPINGS, name) || isdefined(Almanac, Symbol(name)))
            t.skipped += 1
            note_skip(
                name, c.depth > 1 ? "nested" : c.failed_before ? "error_state" :
                    truncated ? "truncated" : desync ? "after_truncated_load" : "no_mapping"
            )
            return
        end
        mapping = get(MAPPINGS, name, default_mapping)
        FILE_CHECK[] = ""
        m0 = t.match
        result = try
            mapping(ctx, c)
        catch e
            e
        end
        FAULT == "" || (result = inject_fault(result))
        ORACLE && name in STATEFUL && cspice_call(c)
        oracle_ok() = ORACLE && !(name in STATEFUL) && !(result isa Exception && !(result isa SpiceError)) &&
            agrees_with_cspice(c, result)
        expkind = c.failed ? Symbol(lowercase(replace(c.short, r"^SPICE\(|\)$" => ""))) : nothing
        if result isa UnmappedHandle
            t.skipped += 1
            note_skip(name, "unmapped_handle")
        elseif result isa SpiceError
            if expkind === result.kind
                t.match += 1
            elseif oracle_ok()
                t.cspice += 1
                record!(mismatches, (family, c, "cspice", "SPICELIB: $(c.failed ? c.short : "no error"); Almanac and CSPICE: $(result.kind)"))
            else
                t.error += 1
                record!(
                    mismatches, (
                        family, c, "error", "SPICELIB: $(c.failed ? c.short : "no error"); " *
                            "Almanac: $(result.kind): $(first(sprint(showerror, result), 160))",
                    )
                )
            end
        elseif result isa Exception
            t.adapter += 1
            record!(mismatches, (family, c, "adapter", first(sprint(showerror, result), 200)))
        elseif c.failed && oracle_ok()
            t.cspice += 1
            record!(mismatches, (family, c, "cspice", "SPICELIB: $(c.short); Almanac and CSPICE: no error"))
        elseif c.failed
            t.error += 1
            record!(mismatches, (family, c, "error", "SPICELIB: $(c.short); Almanac: no error"))
        else
            bad = [
                (n, e, g) for (n, e, g) in result if
                    g === nothing || length(e) != length(g) || !all(close_enough.(e, g))
            ]
            ORACLE_ALL && !(name in STATEFUL) && compare_all!(name, c, result)
            if isempty(bad)
                t.match += 1
                AUDIT && audit!(name, c, result)
            elseif oracle_ok()
                t.cspice += 1
                n, e, g = bad[1]
                record!(mismatches, (family, c, "cspice", "$n: SPICELIB $(e) Almanac and CSPICE $(first(string(g), 300))"))
            else
                n, e, g = bad[1]
                near = all(((n, e, g),) -> g !== nothing && length(e) == length(g) && near_output(e, g), bad)
                kind = near ? "near" : ORACLE && within_spread(c, result) ? "spread" : "value"
                kind == "near" ? (t.near += 1) : kind == "spread" ? (t.spread += 1) : (t.value += 1)
                d = g === nothing ? "nothing" : string(g)
                record!(
                    mismatches, (
                        family, c, kind, "$n: SPICELIB $(e) Almanac $(first(d, 300))" *
                            (g === nothing ? "" : @sprintf(" (max rel diff %.2e)", maxdiff(e, g))),
                    )
                )
            end
        end
        # files.tsv: each file that SPICELIB wrote and closed (the recorder copied it),
        # and whether the close compared it
        if occursin(CLOSERS, name) && !isempty(c.snapshot) && mismatches.files !== nothing
            why = result isa UnmappedHandle && isempty(FILE_CHECK[]) ?
                "not compared: " * get(DROPPED, c.inputs[1].values[1], "the handle is not replayed (an open outside the log)") :
                FILE_CHECK[] == "compared" ? "compared: " * (t.match > m0 ? "match" : "mismatch") :
                isempty(FILE_CHECK[]) ? "not compared: " * (result isa Exception ? "the close failed" : "SPICELIB's close failed") :
                FILE_CHECK[]
            println(mismatches.files, family, '\t', c.line, '\t', name, '\t', why)
        end
    end
end

"""
The mismatches go to `replay.tsv` as they are found, the first `MAXROWS` of each routine
and kind; the others are only counted. Every mismatch (each kind but `match` and
`skipped`) also goes to `ids.tsv` with its identity: family, routine, kind and a hash of
all its inputs, which does not depend on the line of the log. run.jl compares these
identities with the baseline, so that a new mismatch fails even when another one went
away.
"""
mutable struct Mismatches
    io::IO
    rows::Dict{Tuple{String, String}, Int}
    ids::Union{IO, Nothing}
    files::Union{IO, Nothing}
end
Mismatches(io::IO, rows) = Mismatches(io, rows, nothing, nothing)
Mismatches(io::IO, rows, ids) = Mismatches(io, rows, ids, nothing)
const MAXROWS = 200
const GATED_KINDS = ("cspice", "near", "spread", "value", "error", "adapter")

"FNV-1a hash of a string, 64 bits: stable across Julia versions, unlike `hash`."
function fnv64(s::AbstractString)
    h = 0xcbf29ce484222325
    for b in codeunits(s)
        h = (h ⊻ b) * 0x00000100000001b3
    end
    return h
end

"The identity of a mismatch of the call `c`: a hash of its name and all its inputs."
call_id(c::Call) = string(
    fnv64(c.name * "(" * join(["$(a.name)=$(repr(a.values))" for a in c.inputs], ",") * ")");
    base = 16, pad = 16
)

function record!(m::Mismatches, (fam, c, kind, msg))
    if m.ids !== nothing && kind in GATED_KINDS
        println(m.ids, fam, '\t', lowercase(c.name), '\t', kind, '\t', call_id(c))
    end
    key = (c.name, kind)
    n = get(m.rows, key, 0) + 1
    m.rows[key] = n
    n <= MAXROWS || return
    return println(m.io, fam, '\t', c.line, '\t', kind, '\t', describe(c), '\t', replace(msg, '\n' => ' '))
end

function describe(c::Call)
    ins = join(["$(a.name)=$(length(a.values) == 1 ? repr(a.values[1]) : repr(a.values))" for a in c.inputs], ", ")
    return "$(c.name)($(first(ins, 240)))"
end

"The logs of LOGDIR (`F_XXX.log` or `F_XXX.log.gz`), of the families named, or all."
function family_logs(logdir, families)
    logs = sort(filter(f -> occursin(r"^F_.*\.log(\.gz)?$", basename(f)), readdir(logdir; join = true)))
    wanted = Set("F_" * uppercase(replace(fam, r"^F_"i => "")) for fam in families)
    isempty(families) || filter!(f -> uppercase(replace(basename(f), r"\.log(\.gz)?$" => "")) in wanted, logs)
    return logs
end

const TALLY_FIELDS = fieldnames(Tally)

"`tallies` as a table, one routine per line: the name and the counts of `Tally`."
function write_tallies(path, tallies)
    return open(path, "w") do io
        println(io, join(["routine"; String.(collect(TALLY_FIELDS))], '\t'))
        for name in sort(collect(keys(tallies)))
            t = tallies[name]
            println(io, join([name; [string(getfield(t, f)) for f in TALLY_FIELDS]], '\t'))
        end
    end
end

# The results go to TSPICE_OUT (default LOGDIR): replay.tsv, the mismatches, and
# tally.tsv, the counts per routine (read by run.jl).
function main(args)
    logdir = abspath(args[1])
    outdir = abspath(get(ENV, "TSPICE_OUT", logdir))
    mkpath(outdir)
    # The writers that the families call create files in the working directory.
    cd(mktempdir())
    logs = family_logs(logdir, args[2:end])
    tallies = Dict{String, Tally}()
    mismatches = Mismatches(
        open(joinpath(outdir, "replay.tsv"), "w"),
        Dict{Tuple{String, String}, Int}(), open(joinpath(outdir, "ids.tsv"), "w"),
        open(joinpath(outdir, "files.tsv"), "w")
    )
    # families.tsv: the status of each family, `done` or `failed` with the error. A family
    # that fails is not a mismatch that the baseline can allow: the exit status is 1.
    failed = 0
    status = open(joinpath(outdir, "families.tsv"), "w")
    for path in logs
        fam = replace(basename(path), r"\.log(\.gz)?$" => "")
        println(stderr, basename(path))
        flush(stderr)
        try
            replay_family(path, tallies, mismatches)
            println(status, fam, "\tdone")
        catch e
            failed += 1
            msg = replace(first(sprint(showerror, e), 300), r"\s+" => " ")
            println("family $fam FAILED: ", msg)
            println(status, fam, "\tfailed\t", msg)
        end
        flush(status)
    end
    close(status)
    close(mismatches.ids)
    close(mismatches.files)
    println(rpad("routine", 10), lpad("match", 10), lpad("cspice", 8), lpad("near", 8), lpad("spread", 8), lpad("value", 8), lpad("error", 8), lpad("adapter", 10), lpad("skipped", 10))
    tot = Tally()
    for name in sort(collect(keys(tallies)))
        t = tallies[name]
        println(rpad(name, 10), lpad(t.match, 10), lpad(t.cspice, 8), lpad(t.near, 8), lpad(t.spread, 8), lpad(t.value, 8), lpad(t.error, 8), lpad(t.adapter, 10), lpad(t.skipped, 10))
        for f in fieldnames(Tally)
            setfield!(tot, f, getfield(tot, f) + getfield(t, f))
        end
    end
    println(rpad("all", 10), lpad(tot.match, 10), lpad(tot.cspice, 8), lpad(tot.near, 8), lpad(tot.spread, 8), lpad(tot.value, 8), lpad(tot.error, 8), lpad(tot.adapter, 10), lpad(tot.skipped, 10))
    close(mismatches.io)
    write_tallies(joinpath(outdir, "tally.tsv"), tallies)
    # skips.tsv: per routine, the skipped calls by reason (a call made by another recorded
    # routine is `nested`, made while SPICELIB was in an error state `error_state`, ...)
    open(joinpath(outdir, "skips.tsv"), "w") do io
        println(io, join(["routine"; collect(SKIP_REASONS)], '\t'))
        for (r, d) in sort(collect(SKIPS); by = first)
            println(io, join([r; [string(get(d, k, 0)) for k in SKIP_REASONS]], '\t'))
        end
    end
    open(joinpath(outdir, "oracle_failures.tsv"), "w") do io
        for (r, (n, msg)) in sort(collect(ORACLE_FAILURES); by = first)
            println(io, r, '\t', n, '\t', msg)
        end
    end
    isempty(ORACLE_FAILURES) ||
        println("the oracle failed for ", sum(first, values(ORACLE_FAILURES)), " calls: see oracle_failures.tsv")
    ORACLE_ALL && open(joinpath(outdir, "cspice_agreement.tsv"), "w") do io
        println(io, "routine\tcalls\talmanac_exact\talmanac_rtol\talmanac_near\talmanac_far\tspicelib_exact\tspicelib_rtol\tspicelib_near\tspicelib_far")
        for (r, v) in sort(collect(AGREEMENT); by = first)
            println(io, r, '\t', join(v, '\t'))
        end
    end
    AUDIT && open(joinpath(outdir, "audit.tsv"), "w") do io
        println(io, "routine\tmatch\tnone_compared\tfewer_compared\tnot_compared")
        for (r, v) in sort(collect(AUDITS); by = first)
            m = get(AUDIT_MISSING, r, Dict{String, Int}())
            println(io, r, '\t', join(v, '\t'), '\t', join(["$k=$n" for (k, n) in sort(collect(m))], ' '))
        end
    end
    return failed
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS) == 0 ? 0 : 1)
end
