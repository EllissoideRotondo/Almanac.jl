# Replay of the calls that NAIF's test suite (tspice) made to SPICELIB, against
# Almanac.jl, and comparison of the results.
#
#     julia --project=test test/tspice/replay.jl LOGDIR [FAMILY ...]
#
# LOGDIR holds one log per tspice family (`F_XXX.log`, written by record.sh through the
# instrumented SPICELIB interface of PEPPER, see instrument.py), with the copies of the
# working directory that the log names. Each family is replayed in a new
# `SpiceContext`, call by call, in the order of the log. A call is:
#
# - `match`: Almanac gives the outputs of SPICELIB within `RTOL`, or throws the
#   `SpiceError` whose kind is the short message that SPICELIB signaled;
# - `cspice`: Almanac differs from SPICELIB and agrees with CSPICE (TSPICE_ORACLE=1);
# - `near`: the outputs differ by more than `RTOL` and less than `NEAR` relative;
# - `value`: the outputs differ by more;
# - `error`: SPICELIB signaled an error and Almanac did not, or the reverse, or the kinds
#   differ;
# - `adapter`: the call could not be made (a method error of the mapping below), which
#   is a fault of this script, not of Almanac;
# - `skipped`: made by another recorded routine, or made while SPICELIB was in an error
#   state (it then does nothing), or with an array too long for the log, or a routine
#   with no mapping.
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
end

"""
Call `f(call)` for each call of the log at `path`, in order, reading the log as a stream.
A call of a routine for which `want(name)` is false is passed with no arguments, which
saves the parsing of its values.
"""
function foreach_call(f, path; want=_ -> true)
    open(path) do io
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
            ins, outs, snap = Arg[], Arg[], ""
            keep = want(lowercase(name))
            failed, short = false, ""
            while !eof(io)
                l = nextline()
                if startswith(l, "! ")
                    p = split(l, ' '; limit=3)
                    failed = p[2] == "T"
                    short = length(p) == 3 ? strip(p[3]) : ""
                    break
                elseif startswith(l, "S ")
                    snap = l[3:end]
                elseif !keep
                    continue
                elseif startswith(l, "<") || startswith(l, ">")
                    p = split(l)
                    length(p) == 3 || continue
                    aname, k, n = p
                    n = parse(Int, n)
                    if n < 0                # more values than the log keeps
                        push!((l[1] == '<' ? ins : outs), Arg(aname[2:end], k[1], Any[TRUNCATED]))
                        continue
                    end
                    vals = Vector{Any}(undef, n)
                    for j in 1:n
                        v = nextline()
                        vals[j] = k == "d" ? parse(Float64, strip(v)) :
                                  k == "i" ? parse(Int, strip(v)) :
                                  k == "l" ? strip(v) == "T" : String(v[2:end])
                    end
                    push!((l[1] == '<' ? ins : outs), Arg(aname[2:end], k[1], vals))
                elseif startswith(l, "@ ")
                    pending = l             # a nested call before this one ended
                    break
                end
            end
            f(Call(name, start, parse(Int, depth), fb == "T", ins, outs, failed, short, snap))
        end
    end
end

# --- argument conversion -----------------------------------------------------------------

"A logged value as Almanac takes it: a scalar, a static vector or a 3x3 or 6x6 matrix."
function jl(a::Arg)
    v = a.values
    a.kind == 'c' && return length(v) == 1 ? v[1] : String.(v)
    n = length(v)
    T = a.kind == 'd' ? Float64 : a.kind == 'i' ? Int : Bool
    n == 1 && return T(v[1])
    n == 9 && return SMatrix{3,3,T}(v...)
    n == 36 && return SMatrix{6,6,T}(v...)
    return SVector{n,T}(v...)
end

"The values of an Almanac result, flattened as SPICELIB's output arguments are."
flat(x::Nothing) = nothing
flat(x::AbstractString) = Any[String(x)]
flat(x::Symbol) = Any[String(x)]
flat(x::Number) = Any[x]
flat(x::AbstractArray{<:Number}) = Any[vec(collect(x))...]
flat(x::AbstractArray{<:AbstractString}) = Any[String.(x)...]
flat(x::AbstractArray) = Any[reduce(vcat, flat.(x))...]
flat(x::SpiceEllipse) = Any[x.center..., x.semi_major..., x.semi_minor...]

close_enough(a::Real, b::Real) = a == b || abs(a - b) <= ATOL + RTOL * max(abs(a), abs(b))
close_enough(a::AbstractString, b::AbstractString) = rstrip(a) == rstrip(b)
close_enough(a::Bool, b::Bool) = a == b
close_enough(a, b) = false

"Largest relative difference of two value lists (for the report)."
function maxdiff(exp, got)
    m = 0.0
    for (a, b) in zip(exp, got)
        (a isa Real && b isa Real) || continue
        m = max(m, abs(a - b) / max(abs(a), abs(b), 1e-300))
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
const DROP = Dict("gdpool" => ["start", "room"], "gipool" => ["start", "room"],
    "gcpool" => ["start", "room"], "bodvrd" => ["maxn"], "bodvcd" => ["maxn"],
    "pdpool" => ["n"], "pipool" => ["n"], "pcpool" => ["n"], "lmpool" => ["n"])

# Output arguments that Almanac does not return (counts of returned arrays).
const NOTOUT = Dict("gdpool" => ["n"], "gipool" => ["n"], "gcpool" => ["n"],
    "bodvrd" => ["dim"], "bodvcd" => ["dim"], "dtpool" => ["found"])

function default_mapping(ctx, c::Call)
    name = lowercase(c.name)
    f = getfield(Almanac, Symbol(name))
    drop = get(DROP, name, String[])
    args = Any[jl(a) for a in c.inputs if !(a.name in drop)]
    r = applicable(f, ctx, args...) ? f(ctx, args...) : f(args...)
    notout = get(NOTOUT, name, String[])
    outs = [o for o in c.outputs if !(o.name in notout)]
    found = findfirst(o -> o.name == "found", outs)
    if found !== nothing
        fexp = outs[found].values[1]
        deleteat!(outs, found)
        r === nothing && return [("found", Any[fexp], Any[false])]
        fexp || return [("found", Any[false], Any[true])]
    end
    vals = r isa Tuple ? collect(r) : Any[r]
    return output_pairs(name, outs, vals)
end

# --- the CSPICE oracle ----------------------------------------------------------------

"Arguments as SPICE.jl takes them: plain arrays."
plain(x::StaticArray) = Array(x)
plain(x) = x

"The routines that change the kernels or the settings, mirrored into CSPICE."
const STATEFUL = Set(["furnsh", "unload", "kclear", "ldpool", "clpool", "pdpool", "pcpool",
    "pipool", "lmpool", "dvpool", "timdef", "spklef", "cklpf", "pcklof"])

"The kind of an error of SPICE.jl (its message starts with the short message of CSPICE)."
function cspice_kind(e)
    m = match(r"SPICE\((\w+)\)", sprint(showerror, e))
    return m === nothing ? :unknown : Symbol(lowercase(m.captures[1]))
end

"Call the SPICE.jl function of the routine of `c` with its inputs; an exception is returned."
function cspice_call(c::Call)
    name = lowercase(c.name)
    name in ("spklef", "cklpf", "pcklof") && (name = "furnsh")
    isdefined(SPICE, Symbol(name)) || return nothing
    f = getfield(SPICE, Symbol(name))
    drop = get(DROP, lowercase(c.name), String[])
    args = Any[plain(jl(a)) for a in c.inputs if !(a.name in drop)]
    name in ("pdpool", "pipool", "pcpool") && (args = Any[args[1], plain(jl(c.inputs[end]))])
    name == "lmpool" && (args = Any[String.(c.inputs[1].values)])
    try
        return in_snapshot(() -> f(args...), c)
    catch e
        return e
    end
end

function cspice_reset()
    SPICE.kclear()
    SPICE.timdef("SET", "CALENDAR", "GREGORIAN")
    SPICE.timdef("SET", "SYSTEM", "UTC")
end

"""
True if CSPICE gives for `c` what Almanac gave: the same error kind, or outputs that
agree within `NEAR` with `got`, the outputs of Almanac paired with SPICELIB's.
"""
function agrees_with_cspice(c::Call, almanac_result)
    r = cspice_call(c)
    r === nothing && return false
    if almanac_result isa SpiceError
        return r isa Exception && cspice_kind(r) === almanac_result.kind
    end
    r isa Exception && return false
    name = lowercase(c.name)
    notout = get(NOTOUT, name, String[])
    outs = [o for o in c.outputs if !(o.name in notout && o.name != "found")]
    filter!(o -> o.name != "found", outs)
    vals = r isa Tuple ? collect(r) : Any[r]
    oracle = output_pairs(name, outs, vals)
    length(oracle) == length(almanac_result) || return false
    for ((_, _, a), (_, _, o)) in zip(almanac_result, oracle)
        (a === nothing || o === nothing) && return a === o
        length(a) == length(o) || return false
        all(t -> close_enough(t[1], t[2]) || (t[1] isa Real && t[2] isa Real &&
            abs(t[1] - t[2]) <= NEAR * max(abs(t[1]), abs(t[2]), 1.0)), zip(a, o)) || return false
    end
    return true
end

function output_pairs(name, outs, vals)
    pairs = Tuple{String,Vector{Any},Any}[]
    if length(vals) != length(outs)
        # a single array output that SPICELIB splits into several, or the reverse
        flatvals = reduce(vcat, [flat(v) for v in vals if v !== nothing]; init=Any[])
        exp = reduce(vcat, [o.values for o in outs]; init=Any[])
        push!(pairs, (join([o.name for o in outs], "+"), exp, flatvals))
        return pairs
    end
    for (o, v) in zip(outs, vals)
        push!(pairs, (o.name, o.values, flat(v)))
    end
    return pairs
end

"SPICELIB's handles of the files loaded by handle in this family, and their paths."
const HANDLES = Dict{Int,String}()

function load_by_handle(ctx, c)
    path = jl(c.inputs[1])
    in_snapshot(c) do
        furnsh(ctx, path)
    end
    c.failed || (HANDLES[c.outputs[1].values[1]] = path)
    return []
end

function unload_by_handle(ctx, c)
    path = get(HANDLES, jl(c.inputs[1]), nothing)
    path === nothing || unload(ctx, path)
    return []
end

# Routines whose arguments do not follow the default.
const MAPPINGS = Dict{String,Function}(
    "spklef" => load_by_handle, "cklpf" => load_by_handle, "pcklof" => load_by_handle,
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
        r === nothing && return [("found", Any[found], Any[false])]
        found || return [("found", Any[false], Any[true])]
        return [("n,type", Any[n, typ], Any[r[1], string(r[2])])]
    end,
    "furnsh" => (ctx, c) -> (in_snapshot(c) do; furnsh(ctx, jl(c.inputs[1])); end; []),
    "ldpool" => (ctx, c) -> (in_snapshot(c) do; ldpool(ctx, jl(c.inputs[1])); end; []),
    "lmpool" => (ctx, c) -> (lmpool(ctx, String.(c.inputs[1].values)); []),
    "pdpool" => (ctx, c) -> (pdpool(ctx, jl(c.inputs[1]), Float64.(c.inputs[3].values)); []),
    "pipool" => (ctx, c) -> (pipool(ctx, jl(c.inputs[1]), Int.(c.inputs[3].values)); []),
    "pcpool" => (ctx, c) -> (pcpool(ctx, jl(c.inputs[1]), String.(c.inputs[3].values)); []),
    "gnpool" => function (ctx, c)
        # Almanac returns the names in ascending order, CSPICE in the order of its hash
        # table: compare the sets when SPICELIB returned all of them.
        a = Dict(x.name => jl(x) for x in c.inputs)
        out = Dict(o.name => o.values for o in c.outputs)
        r = gnpool(ctx, a["name"])
        found, n = out["found"][1], out["n"][1]
        isempty(r) && return [("found", Any[found], Any[false])]
        found || return [("found", Any[false], Any[true])]
        (a["start"] == 1 && n < a["room"]) || return [("n", Any[min(n, 1)], Any[1])]
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

mutable struct Tally
    match::Int
    cspice::Int
    near::Int
    value::Int
    error::Int
    adapter::Int
    skipped::Int
end
Tally() = Tally(0, 0, 0, 0, 0, 0, 0)

function replay_family(path, tallies, mismatches)
    family = replace(basename(path), ".log" => "")
    ctx = SpiceContext()
    empty!(HANDLES)
    ORACLE && cspice_reset()
    foreach_call(path; want=n -> haskey(MAPPINGS, n) || isdefined(Almanac, Symbol(n))) do c
        name = lowercase(c.name)
        t = get!(Tally, tallies, name)
        truncated = any(a -> a.values == Any[TRUNCATED], c.inputs) ||
                    any(a -> a.values == Any[TRUNCATED], c.outputs)
        if c.depth > 1 || c.failed_before || truncated ||
                !(haskey(MAPPINGS, name) || isdefined(Almanac, Symbol(name)))
            t.skipped += 1
            return
        end
        mapping = get(MAPPINGS, name, default_mapping)
        result = try
            mapping(ctx, c)
        catch e
            e
        end
        ORACLE && name in STATEFUL && cspice_call(c)
        oracle_ok() = ORACLE && !(name in STATEFUL) && !(result isa Exception && !(result isa SpiceError)) &&
                      agrees_with_cspice(c, result)
        expkind = c.failed ? Symbol(lowercase(replace(c.short, r"^SPICE\(|\)$" => ""))) : nothing
        if result isa SpiceError
            if expkind === result.kind
                t.match += 1
            elseif oracle_ok()
                t.cspice += 1
            else
                t.error += 1
                record!(mismatches, (family, c, "error", "SPICELIB: $(c.failed ? c.short : "no error"); " *
                    "Almanac: $(result.kind): $(first(sprint(showerror, result), 160))"))
            end
        elseif result isa Exception
            t.adapter += 1
            record!(mismatches, (family, c, "adapter", first(sprint(showerror, result), 200)))
        elseif c.failed && oracle_ok()
            t.cspice += 1
        elseif c.failed
            t.error += 1
            record!(mismatches, (family, c, "error", "SPICELIB: $(c.short); Almanac: no error"))
        else
            bad = [(n, e, g) for (n, e, g) in result if
                   g === nothing || length(e) != length(g) || !all(close_enough.(e, g))]
            if isempty(bad)
                t.match += 1
            elseif oracle_ok()
                t.cspice += 1
            else
                n, e, g = bad[1]
                near = all(((n, e, g),) -> g !== nothing && length(e) == length(g) &&
                    all(t -> close_enough(t[1], t[2]) || (t[1] isa Real && t[2] isa Real &&
                        abs(t[1] - t[2]) <= NEAR * max(abs(t[1]), abs(t[2]), 1.0)), zip(e, g)), bad)
                near ? (t.near += 1) : (t.value += 1)
                d = g === nothing ? "nothing" : string(g)
                record!(mismatches, (family, c, near ? "near" : "value", "$n: SPICELIB $(e) Almanac $(first(d, 300))" *
                    (g === nothing ? "" : @sprintf(" (max rel diff %.2e)", maxdiff(e, g)))))
            end
        end
    end
end

"""
The mismatches go to `replay.tsv` as they are found, the first `MAXROWS` of each routine
and kind; the others are only counted.
"""
mutable struct Mismatches
    io::IO
    rows::Dict{Tuple{String,String},Int}
end
const MAXROWS = 200

function record!(m::Mismatches, (fam, c, kind, msg))
    key = (c.name, kind)
    n = get(m.rows, key, 0) + 1
    m.rows[key] = n
    n <= MAXROWS || return
    println(m.io, fam, '\t', c.line, '\t', kind, '\t', describe(c), '\t', replace(msg, '\n' => ' '))
end

function describe(c::Call)
    ins = join(["$(a.name)=$(length(a.values) == 1 ? repr(a.values[1]) : repr(a.values))" for a in c.inputs], ", ")
    return "$(c.name)($(first(ins, 240)))"
end

function main(args)
    logdir = args[1]
    families = args[2:end]
    logs = sort(filter(f -> endswith(f, ".log"), readdir(logdir; join=true)))
    wanted = Set("F_" * uppercase(replace(fam, r"^F_"i => "")) * ".LOG" for fam in families)
    isempty(families) || filter!(f -> uppercase(basename(f)) in wanted, logs)
    tallies = Dict{String,Tally}()
    mismatches = Mismatches(open(joinpath(logdir, "replay.tsv"), "w"), Dict{Tuple{String,String},Int}())
    for path in logs
        println(stderr, basename(path))
        flush(stderr)
        try
            replay_family(path, tallies, mismatches)
        catch e
            println("family $(basename(path)): ", first(sprint(showerror, e), 300))
        end
    end
    println(rpad("routine", 10), lpad("match", 8), lpad("cspice", 8), lpad("near", 7), lpad("value", 7), lpad("error", 7), lpad("adapter", 9), lpad("skipped", 9))
    tot = Tally()
    for name in sort(collect(keys(tallies)))
        t = tallies[name]
        println(rpad(name, 10), lpad(t.match, 8), lpad(t.cspice, 8), lpad(t.near, 7), lpad(t.value, 7), lpad(t.error, 7), lpad(t.adapter, 9), lpad(t.skipped, 9))
        for f in fieldnames(Tally)
            setfield!(tot, f, getfield(tot, f) + getfield(t, f))
        end
    end
    println(rpad("all", 10), lpad(tot.match, 8), lpad(tot.cspice, 8), lpad(tot.near, 7), lpad(tot.value, 7), lpad(tot.error, 7), lpad(tot.adapter, 9), lpad(tot.skipped, 9))
    close(mismatches.io)
    return tot
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
