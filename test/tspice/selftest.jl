# Self-test of the gate of run.jl: it must reject what it is meant to reject.
#
#     julia --project=test test/tspice/selftest.jl [DIR]
#
# DIR is a recording of run.jl (default ~/.cache/almanac-tspice) with the families
# F_LTIME and F_POOL. Each case runs run.jl on those families in a scratch copy of the
# recording and checks its exit status: 0 for the recording as it is, 1 for a mismatch
# that the baseline lacks and for a fault injected into Almanac's results (TSPICE_FAULT),
# 2 for a family without a log, a family that failed to record, and a family whose log
# cannot be replayed. The baseline files are copied, changed and
# restored; the scratch directory is removed at the end.

const HERE = @__DIR__
const FAMS = ["F_LTIME", "F_POOL"]

run_status(dir, fams) = run(ignorestatus(`$(Base.julia_cmd()) --startup-file=no --project=$(dirname(HERE)) $(joinpath(HERE, "run.jl")) --dir $dir --jobs 1 $fams`)).exitcode

"A scratch recording with the logs of `FAMS` linked from `src`, its manifest and status."
function scratch(src)
    dir = mktempdir()
    logs = mkpath(joinpath(dir, "logs"))
    for f in FAMS, ext in (".log.gz", ".log.files")
        p = joinpath(src, "logs", f * ext)
        ispath(p) && symlink(p, joinpath(logs, f * ext))
    end
    cp(joinpath(src, "logs", "manifest.tsv"), joinpath(logs, "manifest.tsv"))
    write(joinpath(logs, "status.tsv"), join(["$f\tpassed\n" for f in FAMS]))
    return dir
end

function with_baseline_edit(f, edit)
    path = joinpath(HERE, "baseline_ids.tsv.gz")
    saved = read(path)
    try
        lines = readlines(`gzip -dc $path`)
        open(pipeline(`gzip -9n`; stdout = path), "w") do io
            foreach(l -> println(io, l), edit(lines))
        end
        return f()
    finally
        write(path, saved)
    end
end

function main(args)
    src = isempty(args) ? joinpath(homedir(), ".cache", "almanac-tspice") : abspath(args[1])
    results = Pair{String, Bool}[]
    check(name, got, want) = (
        push!(results, name => got == want);
        println(got == want ? "ok    " : "FAIL  ", name, ": exit ", got, ", expected ", want)
    )

    dir = scratch(src)
    check("the recording as it is", run_status(dir, FAMS), 0)

    # a mismatch of F_POOL removed from the baseline is a new mismatch
    got = with_baseline_edit(
        ls -> (
            i = findfirst(l -> startswith(l, "F_POOL\t"), ls);
            i === nothing ? ls : deleteat!(ls, i)
        )
    ) do
        run_status(dir, FAMS)
    end
    check("a mismatch that the baseline lacks", got, 1)

    # faults of Almanac that the comparison must find: an output off by 1e-10 relative,
    # which CSPICE does not share, and an error of another kind
    for (fault, what) in ("value" => "an output off by 1e-10", "kind" => "an error of another kind")
        got = withenv(() -> run_status(dir, FAMS), "TSPICE_FAULT" => fault)
        check(what, got, 1)
    end

    rm(joinpath(dir, "logs", "F_LTIME.log.gz"))
    check("a family without a log", run_status(dir, FAMS), 2)

    write(joinpath(dir, "logs", "F_LTIME.log.gz"), "not a gzip file")
    check("a family that cannot be replayed", run_status(dir, FAMS), 2)

    rm(joinpath(dir, "logs", "F_LTIME.log.gz"))
    symlink(joinpath(src, "logs", "F_LTIME.log.gz"), joinpath(dir, "logs", "F_LTIME.log.gz"))
    write(joinpath(dir, "logs", "status.tsv"), "F_LTIME\tfailed\ttimeout\nF_POOL\tpassed\n")
    check("a family that failed to record", run_status(dir, FAMS), 2)

    rm(dir; recursive = true)
    nfail = count(!last, results)
    println(nfail == 0 ? "the gate rejects every case" : "$nfail case(s) not rejected")
    return nfail == 0 ? 0 : 1
end

exit(main(ARGS))
