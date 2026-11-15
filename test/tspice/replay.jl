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
# - `value`: the outputs differ;
# - `error`: SPICELIB signaled an error and Almanac did not, or the reverse, or the kinds
#   differ;
# - `adapter`: the call could not be made (a method error of the mapping below), which
#   is a fault of this script, not of Almanac;
# - `skipped`: made by another recorded routine, or made while SPICELIB was in an error
#   state (it then does nothing), or a routine with no mapping.
#
# The output is a table per routine and the first mismatches of each kind, and
# `LOGDIR/replay.tsv` with one line per mismatch.

using Almanac, Printf, LinearAlgebra, StaticArrays

const RTOL = parse(Float64, get(ENV, "TSPICE_RTOL", "1e-12"))
const ATOL = parse(Float64, get(ENV, "TSPICE_ATOL", "1e-14"))

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

function parse_log(path)
    calls = Call[]
    lines = readlines(path)
    i = 1
    while i <= length(lines)
        l = lines[i]
        if !startswith(l, "@ ")
            i += 1
            continue
        end
        _, name, depth, fb = split(l)
        start = i
        ins, outs, snap = Arg[], Arg[], ""
        failed, short = false, ""
        i += 1
        while i <= length(lines)
            l = lines[i]
            if startswith(l, "! ")
                parts = split(l, ' '; limit=3)
                failed = parts[2] == "T"
                short = length(parts) == 3 ? strip(parts[3]) : ""
                i += 1
                break
            elseif startswith(l, "S ")
                snap = l[3:end]
                i += 1
            elseif startswith(l, "<") || startswith(l, ">")
                aname, k, n = split(l)
                n = parse(Int, n)
                vals = Any[]
                for j in 1:n
                    v = lines[i+j]
                    push!(vals, k == "d" ? parse(Float64, replace(strip(v), "D" => "E")) :
                                k == "i" ? parse(Int, strip(v)) :
                                k == "l" ? strip(v) == "T" : String(v[2:end]))
                end
                push!((l[1] == '<' ? ins : outs), Arg(aname[2:end], k[1], vals))
                i += n + 1
            elseif startswith(l, "@ ")
                break                       # a nested call before this one ended
            else
                i += 1
            end
        end
        push!(calls, Call(name, start, parse(Int, depth), fb == "T", ins, outs, failed, short, snap))
    end
    return calls
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
# order of SPICELIB, with the context first when the function has a method for it, and
# pairs its results with the outputs of SPICELIB in order. A `found` output of SPICELIB is
# `nothing` in Almanac.

const PURE = Set(["twovec", "latrec", "reclat", "georec", "recgeo", "sphrec", "recsph",
    "cylrec", "reccyl", "radrec", "recrad", "eul2m", "m2eul", "q2m", "m2q", "axisar",
    "raxisa", "xf2eul", "eul2xf", "rav2xf", "xf2rav", "invert", "vrotv", "vsep", "vnorm",
    "vhat", "ucrss", "vcrss", "nearpt", "surfpt", "edlimb", "inelpl", "nplnpt", "npedln",
    "oscelt", "oscltx", "conics", "prop2b", "stelab", "stlabx", "twopi", "halfpi", "rpd",
    "dpr", "spd", "jyear", "j2000", "b1950"])

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
    r = name in PURE ? f(args...) : f(ctx, args...)
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

# Routines whose arguments do not follow the default.
const MAPPINGS = Dict{String,Function}(
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
    value::Int
    error::Int
    adapter::Int
    skipped::Int
end
Tally() = Tally(0, 0, 0, 0, 0)

function replay_family(path, tallies, mismatches)
    family = replace(basename(path), ".log" => "")
    ctx = SpiceContext()
    for c in parse_log(path)
        name = lowercase(c.name)
        t = get!(Tally, tallies, name)
        if c.depth > 1 || c.failed_before || !(haskey(MAPPINGS, name) || isdefined(Almanac, Symbol(name)))
            t.skipped += 1
            continue
        end
        mapping = get(MAPPINGS, name, default_mapping)
        result = try
            mapping(ctx, c)
        catch e
            e
        end
        expkind = c.failed ? Symbol(lowercase(replace(c.short, r"^SPICE\(|\)$" => ""))) : nothing
        if result isa SpiceError
            if expkind === result.kind
                t.match += 1
            else
                t.error += 1
                push!(mismatches, (family, c, "error", "SPICELIB: $(c.failed ? c.short : "no error"); " *
                    "Almanac: $(result.kind): $(first(sprint(showerror, result), 160))"))
            end
        elseif result isa Exception
            t.adapter += 1
            push!(mismatches, (family, c, "adapter", first(sprint(showerror, result), 200)))
        elseif c.failed
            t.error += 1
            push!(mismatches, (family, c, "error", "SPICELIB: $(c.short); Almanac: no error"))
        else
            bad = [(n, e, g) for (n, e, g) in result if
                   g === nothing || length(e) != length(g) || !all(close_enough.(e, g))]
            if isempty(bad)
                t.match += 1
            else
                t.value += 1
                n, e, g = bad[1]
                d = g === nothing ? "nothing" : string(g)
                push!(mismatches, (family, c, "value", "$n: SPICELIB $(e) Almanac $(first(d, 300))" *
                    (g === nothing ? "" : @sprintf(" (max rel diff %.2e)", maxdiff(e, g)))))
            end
        end
    end
end

function describe(c::Call)
    ins = join(["$(a.name)=$(length(a.values) == 1 ? repr(a.values[1]) : repr(a.values))" for a in c.inputs], ", ")
    return "$(c.name)($(first(ins, 240)))"
end

function main(args)
    logdir = args[1]
    families = args[2:end]
    logs = sort(filter(f -> endswith(f, ".log"), readdir(logdir; join=true)))
    isempty(families) || filter!(f -> any(fam -> occursin(Regex("^F_" * uppercase(replace(fam, r"^F_"i => "")) * r"\.log$"),
        uppercase(basename(f))), families), logs)
    tallies = Dict{String,Tally}()
    mismatches = Any[]
    for path in logs
        try
            replay_family(path, tallies, mismatches)
        catch e
            println("family $(basename(path)): ", first(sprint(showerror, e), 300))
        end
    end
    println(rpad("routine", 10), lpad("match", 8), lpad("value", 7), lpad("error", 7), lpad("adapter", 9), lpad("skipped", 9))
    tot = Tally()
    for name in sort(collect(keys(tallies)))
        t = tallies[name]
        println(rpad(name, 10), lpad(t.match, 8), lpad(t.value, 7), lpad(t.error, 7), lpad(t.adapter, 9), lpad(t.skipped, 9))
        for f in fieldnames(Tally)
            setfield!(tot, f, getfield(tot, f) + getfield(t, f))
        end
    end
    println(rpad("all", 10), lpad(tot.match, 8), lpad(tot.value, 7), lpad(tot.error, 7), lpad(tot.adapter, 9), lpad(tot.skipped, 9))
    open(joinpath(logdir, "replay.tsv"), "w") do io
        for (fam, c, kind, msg) in mismatches
            println(io, fam, '\t', c.line, '\t', kind, '\t', describe(c), '\t', replace(msg, '\n' => ' '))
        end
    end
    return tot
end

abspath(PROGRAM_FILE) == @__FILE__ && main(ARGS)
