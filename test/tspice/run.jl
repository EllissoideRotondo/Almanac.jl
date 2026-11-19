# The tspice test: record the calls of NAIF's test suite to SPICELIB (record.sh) if no
# recording is there yet, replay them against Almanac.jl in parallel processes
# (replay.jl), and compare the mismatches with the baseline of known ones.
#
#     julia --project=test test/tspice/run.jl [options] [FAMILY ...]
#
#     --dir DIR          the recording (DIR/logs) and the results (DIR/results);
#                        default ~/.cache/almanac-tspice
#     --record           record again, even if DIR/logs has a recording
#     --jobs N           parallel processes for recording and replay (default: the CPU
#                        threads, at most one per 3 GB of free memory)
#     --no-oracle        do not check the mismatches with CSPICE (then the baseline,
#                        made with CSPICE, is not applied; failures of the harness count)
#     --update-baseline  make this run the baseline: baseline.tsv (counts, class and note
#                        of each routine; a routine keeps its class), baseline_ids.tsv.gz
#                        (every mismatch) and coverage.tsv (the calls of each routine)
#
# Recording needs gfortran 16 and fpm on PATH and takes 15 to 40 minutes; it writes
# about 2 GB.
#
# A run without FAMILY is the full check. It needs a complete recording: every family of
# DIR/logs/manifest.tsv recorded and passed. It fails (exit status 1) if
#   - a family fails to record or to replay;
#   - a mismatch (any kind but `match` and `skipped`) is not in the baseline: each is
#     identified by family, routine, kind and a hash of its inputs, so a new mismatch
#     fails even when another one went away, and a call that moves from `match` to
#     `near` or `cspice` fails too;
#   - a routine has more calls skipped, `near` or `spread` than coverage.tsv records, or
#     fewer calls compared.
# A run with FAMILY checks the mismatches of those families only. The exit status is 2
# for an incomplete recording or a failure of the harness. The results are
# DIR/results/tally.tsv (counts per routine), DIR/results/replay.tsv (the first
# mismatches of each routine) and DIR/results/ids.tsv (every mismatch).

using SHA: sha256

const HERE = @__DIR__
const COUNTS = ["match", "cspice", "near", "spread", "value", "error", "adapter", "skipped"]
const GATED = ["value", "error", "adapter"]

"One replay process per CPU thread, at most one per 3 GB of free memory (a process of a
large family holds about 2 GB)."
default_jobs() = max(1, min(Sys.CPU_THREADS, floor(Int, Sys.free_memory() / 3.0e9)))

function parse_args(args)
    o = Dict{String, Any}(
        "dir" => joinpath(homedir(), ".cache", "almanac-tspice"),
        "record" => false, "jobs" => default_jobs(), "oracle" => true, "update" => false,
        "families" => String[]
    )
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--dir"
            o["dir"] = abspath(expanduser(args[i + 1])); i += 1
        elseif a == "--jobs"
            o["jobs"] = parse(Int, args[i + 1]); i += 1
        elseif a == "--record"
            o["record"] = true
        elseif a == "--no-oracle"
            o["oracle"] = false
        elseif a == "--update-baseline"
            o["update"] = true
        elseif startswith(a, "--")
            error("unknown option $a")
        else
            push!(o["families"], a)
        end
        i += 1
    end
    return o
end

"Read a table of counts: routine => Dict(kind => count)."
function read_counts(path)
    t = Dict{String, Dict{String, Int}}()
    isfile(path) || return t
    lines = readlines(path)
    head = split(lines[1], '\t')
    for l in lines[2:end]
        f = split(l, '\t')
        t[f[1]] = Dict(String(head[k]) => parse(Int, f[k]) for k in 2:length(f) if head[k] in COUNTS)
    end
    return t
end

"""
The class and note of each routine of the baseline: `documented` (a difference from
SPICELIB that Almanac states), `gap` (a limit of the replay, not of Almanac) or
`unexplained`, with what it is.
"""
function read_notes(path)
    t = Dict{String, Tuple{String, String}}()
    isfile(path) || return t
    lines = readlines(path)
    head = split(lines[1], '\t')
    ic, inote = findfirst(==("class"), head), findfirst(==("note"), head)
    for l in lines[2:end]
        f = split(l, '\t')
        t[f[1]] = (ic === nothing ? "unexplained" : String(f[ic]), inote === nothing ? "" : String(f[inote]))
    end
    return t
end

function write_counts(path, t, kinds; notes = nothing)
    return open(path, "w") do io
        println(io, join(["routine"; kinds; notes === nothing ? String[] : ["class", "note"]], '\t'))
        for r in sort(collect(keys(t)))
            extra = notes === nothing ? String[] : collect(get(notes, r, ("unexplained", "")))
            println(io, join([r; [string(get(t[r], k, 0)) for k in kinds]; extra], '\t'))
        end
    end
end

"The logs split into `n` groups of about equal size, the largest first."
function split_logs(logs, n)
    groups = [String[] for _ in 1:n]
    load = zeros(Int, n)
    for f in sort(logs; by = filesize, rev = true)
        k = argmin(load)
        push!(groups[k], f)
        load[k] += filesize(f)
    end
    return filter(!isempty, groups)
end

function main(args)
    o = parse_args(args)
    dir, jobs = o["dir"], o["jobs"]
    logdir = joinpath(dir, "logs")
    full = isempty(o["families"])
    have = isfile(joinpath(logdir, "manifest.tsv"))
    if o["record"] || !have
        println("recording into $dir with $jobs processes ...")
        ok = withenv("JOBS" => string(jobs)) do
            success(`bash $(joinpath(HERE, "record.sh")) $dir $(o["families"])`)
        end
        ok || (println("the recording failed: see $(joinpath(logdir, "status.tsv"))"); return 2)
    end
    problems = check_recording(logdir, full ? nothing : o["families"])
    if !isempty(problems)
        println("the recording is not complete:")
        foreach(p -> println("  ", p), first(problems, 20))
        length(problems) > 20 && println("  ... $(length(problems) - 20) more")
        println("record again with --record")
        return 2
    end
    logs = filter(f -> endswith(f, ".log.gz"), readdir(logdir; join = true))
    if !full
        wanted = Set("F_" * uppercase(replace(f, r"^F_"i => "")) * ".LOG.GZ" for f in o["families"])
        filter!(f -> uppercase(basename(f)) in wanted, logs)
    else
        expected = Set(manifest_families(logdir))
        filter!(f -> replace(basename(f), ".log.gz" => "") in expected, logs)
    end
    isempty(logs) && error("no logs in $logdir")
    results = joinpath(dir, "results")
    rm(results; force = true, recursive = true)
    groups = split_logs(logs, jobs)
    println("replaying $(length(logs)) families in $(length(groups)) processes ...")
    project = dirname(HERE)                # test/
    procs = map(enumerate(groups)) do (k, g)
        out = joinpath(results, "worker$k")
        mkpath(out)
        names = [replace(basename(f), ".log.gz" => "") for f in g]
        env = copy(ENV)
        env["TSPICE_OUT"] = out
        env["TSPICE_ORACLE"] = o["oracle"] ? "1" : "0"
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$project $(joinpath(HERE, "replay.jl")) $logdir $names`
        run(pipeline(setenv(cmd, env); stdout = joinpath(out, "out.txt"), stderr = joinpath(out, "err.txt")); wait = false)
    end
    foreach(wait, procs)
    # every family of every worker must be replayed to the end
    failed = String[]
    for (k, (p, g)) in enumerate(zip(procs, groups))
        w = joinpath(results, "worker$k")
        st = isfile(joinpath(w, "families.tsv")) ?
            Dict(first(split(l, '\t')) => l for l in readlines(joinpath(w, "families.tsv"))) : Dict()
        for f in g
            fam = replace(basename(f), ".log.gz" => "")
            l = get(st, fam, nothing)
            (l !== nothing && endswith(l, "\tdone")) ||
                push!(failed, l === nothing ? "$fam: not replayed (worker $k, exit $(p.exitcode))" : l)
        end
        success(p) || isempty(failed) && push!(failed, "worker $k: exit $(p.exitcode), see $(joinpath(w, "err.txt"))")
    end
    if !isempty(failed)
        println("the replay of some families failed:")
        foreach(l -> println("  ", l), failed)
        return 2
    end
    # merge
    total = Dict{String, Dict{String, Int}}()
    ids = Dict{String, Int}()
    open(joinpath(results, "replay.tsv"), "w") do io
        for k in 1:length(groups)
            w = joinpath(results, "worker$k")
            for l in eachline(joinpath(w, "ids.tsv"))
                ids[l] = get(ids, l, 0) + 1
            end
            for (r, c) in read_counts(joinpath(w, "tally.tsv"))
                d = get!(total, r, Dict{String, Int}())
                for (kind, n) in c
                    d[kind] = get(d, kind, 0) + n
                end
            end
            write(io, read(joinpath(w, "replay.tsv")))
        end
    end
    write_counts(joinpath(results, "tally.tsv"), total, COUNTS)
    write_ids(joinpath(results, "ids.tsv"), ids)
    sums = Dict(k => sum(get(c, k, 0) for c in values(total)) for k in COUNTS)
    println(join(["$k $(sums[k])" for k in COUNTS], ", "))
    baseline_path = joinpath(HERE, "baseline.tsv")
    if o["update"]
        (full && o["oracle"]) || (println("--update-baseline needs a full run with the oracle"); return 2)
        # a routine keeps its class and note; a new one is `unexplained` until classed
        notes = read_notes(baseline_path)
        write_counts(
            baseline_path, Dict(r => c for (r, c) in total if any(get(c, k, 0) > 0 for k in GATED)),
            GATED; notes
        )
        write_ids_gz(joinpath(HERE, "baseline_ids.tsv.gz"), ids)
        write_counts(joinpath(HERE, "coverage.tsv"), total, COVERED)
        println("baseline written to $baseline_path, baseline_ids.tsv.gz and coverage.tsv")
        return 0
    end
    if !o["oracle"]
        # without CSPICE a difference that CSPICE shares is a mismatch: the baseline,
        # made with the oracle, does not apply
        println("--no-oracle: the mismatches are not compared with the baseline")
        return 0
    end
    base = read_counts(baseline_path)
    notes = read_notes(baseline_path)
    byclass = Dict{String, Int}()
    for (r, c) in total, k in GATED
        cl = first(get(notes, r, ("unexplained", "")))
        byclass[cl] = get(byclass, cl, 0) + get(c, k, 0)
    end
    println("gated mismatches by class: ", join(["$k $(byclass[k])" for k in sort(collect(keys(byclass)))], ", "))
    worse = String[]
    # 1. every gated mismatch must be one of the baseline, by identity
    baseids = read_ids_gz(joinpath(HERE, "baseline_ids.tsv.gz"))
    fams = Set(replace(basename(f), ".log.gz" => "") for f in logs)
    gone = 0
    for (id, n) in sort(collect(ids); by = first)
        b = get(baseids, id, 0)
        n > b && push!(worse, "new mismatch: $(replace(id, '\t' => ' '))" * (b > 0 ? " ($n, baseline $b)" : ""))
    end
    for (id, b) in baseids
        first(split(id, '\t')) in fams || continue
        gone += max(0, b - get(ids, id, 0))
    end
    gone > 0 && println("$gone mismatches of the baseline are gone: run --update-baseline to record it")
    # 2. the coverage of each routine, for a full run
    if full
        cov = read_counts(joinpath(HERE, "coverage.tsv"))
        for r in sort(collect(union(keys(cov), keys(total))))
            c, b = get(total, r, Dict{String, Int}()), get(cov, r, Dict{String, Int}())
            for k in ("skipped", "near", "spread")
                get(c, k, 0) > get(b, k, 0) && push!(worse, "$r: $(get(c, k, 0)) calls $k (baseline $(get(b, k, 0)))")
            end
            compared(x) = sum(get(x, k, 0) for k in ("match", "cspice", "near", "spread", "value", "error", "adapter"))
            compared(c) < compared(b) && push!(worse, "$r: $(compared(c)) calls compared (baseline $(compared(b)))")
        end
    else
        println("a run of some families: the coverage is not checked")
    end
    if isempty(worse)
        println("every mismatch is in the baseline, and no routine lost coverage")
        return 0
    end
    println("the run differs from the baseline:")
    foreach(l -> println("  ", l), first(worse, 50))
    length(worse) > 50 && println("  ... $(length(worse) - 50) more")
    println("see $(joinpath(results, "replay.tsv")) and $(joinpath(results, "ids.tsv"))")
    return 1
end

"The families of the manifest of a full recording."
manifest_families(logdir) = [
    String(split(l, '\t')[2]) for l in readlines(joinpath(logdir, "manifest.tsv"))
        if startswith(l, "family\t")
]

"""
The problems of the recording in `logdir` for the families (all those of the manifest if
`nothing`): a family without a log or whose last status is not `passed`, and a manifest
made with another instrument.py or routines.txt (reported, not a problem).
"""
function check_recording(logdir, families)
    problems = String[]
    mf = joinpath(logdir, "manifest.tsv")
    isfile(mf) || return ["no manifest.tsv: the recording was not made by a full record.sh run"]
    m = Dict(String(a) => String(b) for (a, b) in (split(l, '\t') for l in readlines(mf)) if a != "family")
    for f in ("instrument.py", "routines.txt")
        h = bytes2hex(open(io -> sha256(io), joinpath(HERE, f)))
        get(m, f, "") == h || println("note: the recording was made with another $f")
    end
    status = Dict{String, String}()
    sf = joinpath(logdir, "status.tsv")
    isfile(sf) && for l in readlines(sf)
        p = split(l, '\t')
        status[p[1]] = join(p[2:end], ' ')
    end
    fams = families === nothing ? manifest_families(logdir) :
        ["F_" * uppercase(replace(f, r"^F_"i => "")) for f in families]
    for fam in fams
        isfile(joinpath(logdir, fam * ".log.gz")) || (push!(problems, "$fam: no log"); continue)
        st = get(status, fam, "no status")
        st == "passed" || push!(problems, "$fam: $st")
    end
    return problems
end

const COVERED = ["match", "cspice", "near", "spread", "value", "error", "adapter", "skipped"]

function write_ids(path, ids)
    return open(path, "w") do io
        for (id, n) in sort(collect(ids); by = first)
            println(io, id, '\t', n)
        end
    end
end

"The identities as a gzip file (`gzip` of the system, as for the logs)."
function write_ids_gz(path, ids)
    return open(pipeline(`gzip -9n`; stdout = path), "w") do io
        for (id, n) in sort(collect(ids); by = first)
            println(io, id, '\t', n)
        end
    end
end

function read_ids_gz(path)
    ids = Dict{String, Int}()
    isfile(path) || return ids
    for l in eachline(`gzip -dc $path`)
        p = rsplit(l, '\t'; limit = 2)
        ids[String(p[1])] = parse(Int, p[2])
    end
    return ids
end

exit(main(ARGS))
