# The tspice test: record the calls of NAIF's test suite to SPICELIB (record.sh) if no
# recording is there yet, replay them against Almanac.jl in parallel processes
# (replay.jl), and compare the mismatches with the baseline of known ones.
#
#     julia --project=test test/tspice/run.jl [options] [FAMILY ...]
#
#     --dir DIR          the recording (DIR/logs) and the results (DIR/results);
#                        default ~/.cache/almanac-tspice
#     --record           record again, even if DIR/logs has a recording
#     --jobs N           parallel processes for recording and replay (default: CPU threads)
#     --no-oracle        do not check the mismatches with CSPICE
#     --update-baseline  write the counts of this run to baseline.tsv
#
# Recording needs gfortran 16 and fpm on PATH and takes 15 to 40 minutes; it writes
# about 1 GB. The exit status is 1 if a routine has more `value`, `error` or `adapter`
# mismatches than baseline.tsv allows, 0 otherwise. The results are DIR/results/tally.tsv
# (counts per routine) and DIR/results/replay.tsv (the mismatches).

const HERE = @__DIR__
const COUNTS = ["match", "cspice", "near", "spread", "value", "error", "adapter", "skipped"]
const GATED = ["value", "error", "adapter"]

function parse_args(args)
    o = Dict{String,Any}("dir" => joinpath(homedir(), ".cache", "almanac-tspice"),
        "record" => false, "jobs" => Sys.CPU_THREADS, "oracle" => true, "update" => false,
        "families" => String[])
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--dir"
            o["dir"] = abspath(expanduser(args[i+1])); i += 1
        elseif a == "--jobs"
            o["jobs"] = parse(Int, args[i+1]); i += 1
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
    t = Dict{String,Dict{String,Int}}()
    isfile(path) || return t
    lines = readlines(path)
    head = split(lines[1], '\t')
    for l in lines[2:end]
        f = split(l, '\t')
        t[f[1]] = Dict(String(head[k]) => parse(Int, f[k]) for k in 2:length(f))
    end
    return t
end

function write_counts(path, t, kinds)
    open(path, "w") do io
        println(io, join(["routine"; kinds], '\t'))
        for r in sort(collect(keys(t)))
            println(io, join([r; [string(get(t[r], k, 0)) for k in kinds]], '\t'))
        end
    end
end

"The logs split into `n` groups of about equal size, the largest first."
function split_logs(logs, n)
    groups = [String[] for _ in 1:n]
    load = zeros(Int, n)
    for f in sort(logs; by=filesize, rev=true)
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
    have = isdir(logdir) && any(f -> endswith(f, ".log.gz"), readdir(logdir))
    if o["record"] || !have
        println("recording into $dir ...")
        withenv("JOBS" => string(jobs)) do
            run(`bash $(joinpath(HERE, "record.sh")) $dir $(o["families"])`)
        end
    end
    logs = filter(f -> endswith(f, ".log.gz"), readdir(logdir; join=true))
    if !isempty(o["families"])
        wanted = Set("F_" * uppercase(replace(f, r"^F_"i => "")) * ".LOG.GZ" for f in o["families"])
        filter!(f -> uppercase(basename(f)) in wanted, logs)
    end
    isempty(logs) && error("no logs in $logdir")
    results = joinpath(dir, "results")
    rm(results; force=true, recursive=true)
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
        run(pipeline(setenv(cmd, env); stdout=joinpath(out, "out.txt"), stderr=joinpath(out, "err.txt")); wait=false)
    end
    foreach(wait, procs)
    # merge
    total = Dict{String,Dict{String,Int}}()
    open(joinpath(results, "replay.tsv"), "w") do io
        for k in 1:length(groups)
            w = joinpath(results, "worker$k")
            isfile(joinpath(w, "tally.tsv")) || error("worker $k failed, see $(joinpath(w, "err.txt"))")
            for (r, c) in read_counts(joinpath(w, "tally.tsv"))
                d = get!(total, r, Dict{String,Int}())
                for (kind, n) in c
                    d[kind] = get(d, kind, 0) + n
                end
            end
            write(io, read(joinpath(w, "replay.tsv")))
        end
    end
    write_counts(joinpath(results, "tally.tsv"), total, COUNTS)
    sums = Dict(k => sum(get(c, k, 0) for c in values(total)) for k in COUNTS)
    println(join(["$k $(sums[k])" for k in COUNTS], ", "))
    baseline_path = joinpath(HERE, "baseline.tsv")
    if o["update"]
        write_counts(baseline_path, Dict(r => c for (r, c) in total if any(get(c, k, 0) > 0 for k in GATED)), GATED)
        println("baseline written to $baseline_path")
        return 0
    end
    base = read_counts(baseline_path)
    worse = String[]
    for (r, c) in sort(collect(total); by=first), k in GATED
        n, b = get(c, k, 0), get(get(base, r, Dict{String,Int}()), k, 0)
        n > b && push!(worse, "$r $k: $n (baseline $b)")
    end
    if isempty(worse)
        println("no routine has more mismatches than the baseline")
        return 0
    end
    println("more mismatches than the baseline:")
    foreach(l -> println("  ", l), worse)
    println("see $(joinpath(results, "replay.tsv"))")
    return 1
end

exit(main(ARGS))
