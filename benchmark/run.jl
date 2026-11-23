# Run the benchmark suite and summarise the speed-ups over CSPICE.
#
#     julia --startup-file=no --project=benchmark benchmark/run.jl [script ...]
#
# Each script runs in a new Julia process, one after the other, and writes
# `benchmark/results/<script>.md`. Then `benchmark/results/summary.md` lists every Almanac row
# that has a CSPICE row in its group, slowest first, and marks the rows that are slower
# than CSPICE. Without arguments the suite is `states aberration rotations api latency`;
# `states` gets `--threads=8` (ALMANAC_BENCH_THREADS overrides it) for its throughput
# table, every other script one thread.
#
# The scripts that are not in the suite need more than the test kernels:
#
#   pepper.jl     the comparison of PEPPER with ANISE; takes a directory of its kernels
#   table_gpu.jl  ephemeris tables on an NVIDIA GPU; needs CUDA.jl
#
# Environment: ALMANAC_BENCH_SECONDS (seconds of samples per row, default 2),
# ALMANAC_BENCH_REPS (repetitions of `latency`, default 3), ALMANAC_KERNEL_DIR (the
# kernels of the tests, as in `test/kernels.jl`).

using Printf, Dates

const SUITE = ["states", "aberration", "rotations", "api", "latency"]
const RESULTS = joinpath(@__DIR__, "results")

function run_script(name)
    threads = name == "states" ? get(ENV, "ALMANAC_BENCH_THREADS", "8") : "1"
    script = joinpath(@__DIR__, name * ".jl")
    @assert isfile(script) "no benchmark script $script"
    println("\n=== ", name, " (", threads, " thread(s)) ===")
    cmd = `$(Base.julia_cmd()) --startup-file=no --threads=$threads --project=$(@__DIR__) $script`
    t = @elapsed run(cmd)
    @printf("=== %s done in %.0f s\n", name, t)
    return nothing
end

"Rows of `results/<name>.tsv` as named tuples; empty if the script wrote no rows."
function read_rows(name)
    path = joinpath(RESULTS, name * ".tsv")
    isfile(path) || return NamedTuple[]
    lines = readlines(path)
    return map(lines[2:end]) do l
        f = split(l, '\t')
        (; script = name, group = f[1], name = f[2], min_ns = parse(Float64, f[3]), speedup = parse(Float64, f[7]))
    end
end

function summary(names)
    rows = [r for n in names for r in read_rows(n) if startswith(r.name, "Almanac") && !isnan(r.speedup)]
    sort!(rows; by = r -> r.speedup)
    io = IOBuffer()
    println(io, "# Speed-up of Almanac over CSPICE\n")
    println(io, "- date: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "; scripts: ", join(names, ", "))
    println(io, "- speed-up: the fastest CSPICE (SPICE.jl) call of the same group divided by the Almanac call,")
    println(io, "  minimum times; above 1 is faster than CSPICE. The tables of each script are in `results/<script>.md`.\n")
    slow = count(r -> r.speedup < 1, rows)
    println(io, slow == 0 ? "Every Almanac call is faster than CSPICE.\n" : "**$slow Almanac call(s) slower than CSPICE**, marked below.\n")
    println(io, "| Script | Case | Almanac call | min ns/query | speed-up vs CSPICE |")
    println(io, "|---|---|---|---:|---:|")
    for r in rows
        @printf(
            io, "| %s | %s | %s | %.1f | %s%.2f%s |\n", r.script, r.group, r.name, r.min_ns,
            r.speedup < 1 ? "**" : "", r.speedup, r.speedup < 1 ? "** slower" : ""
        )
    end
    text = String(take!(io))
    mkpath(RESULTS)
    write(joinpath(RESULTS, "summary.md"), text)
    println("\n", text)
    println("summary written to ", joinpath(RESULTS, "summary.md"))
    return nothing
end

function main(args)
    names = isempty(args) ? SUITE : [replace(a, r"\.jl$" => "") for a in args]
    # Stale tables of an earlier run would enter the summary.
    for n in names
        rm(joinpath(RESULTS, n * ".tsv"); force = true)
    end
    foreach(run_script, names)
    return summary(names)
end

main(ARGS)
