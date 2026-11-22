# Aberration benchmark: Almanac against CSPICE (SPICE.jl).
#
#     julia --startup-file=no --project=benchmark benchmark/aberration.jl
#
# The case of `states.jl`: 100 epochs on de440s.bsp, the Earth relative to the Moon,
# J2000, with the corrections NONE, LT, LT+S and CN+S; then `spkder` with a correction.
# Both engines load de440s.bsp only. The CSPICE calls take the names "EARTH" and
# "MOON", the Almanac calls the NAIF codes.

include("common.jl")
using Almanac
import SPICE

require_kernels(Kernels.DE440S)

const EARTH, MOON, J2000 = 399, 301, 1

# Every loop sums one component, so the calls are not removed.
function loop_neo_ezr(ctx, epochs, abcorr)
    acc = 0.0
    for et in epochs
        s, _ = spkezr(ctx, EARTH, et, J2000, abcorr, MOON)
        acc += s[1]
    end
    return acc
end
function loop_neo_pos(ctx, epochs, abcorr)
    acc = 0.0
    for et in epochs
        p, _ = spkpos(ctx, EARTH, et, J2000, abcorr, MOON)
        acc += p[1]
    end
    return acc
end
function loop_neo_der(ctx, epochs, abcorr, v::Val{N}) where {N}
    acc = 0.0
    for et in epochs
        jet = spkder(ctx, EARTH, et, J2000, abcorr, MOON, v)
        acc += jet[N + 1][1]
    end
    return acc
end
function loop_cspice_ezr(epochs, abcorr)
    acc = 0.0
    for et in epochs
        s, _ = SPICE.spkezr("EARTH", et, "J2000", abcorr, "MOON")
        acc += s[1]
    end
    return acc
end
function loop_cspice_pos(epochs, abcorr)
    acc = 0.0
    for et in epochs
        p, _ = SPICE.spkpos("EARTH", et, "J2000", abcorr, "MOON")
        acc += p[1]
    end
    return acc
end

function main()
    epochs = anise_epochs()
    n = length(epochs)
    SPICE.kclear(); SPICE.furnsh(Kernels.DE440S)
    ctx = SpiceContext(); furnsh(ctx, Kernels.DE440S)

    # One group per call and correction, so that each Almanac row has its CSPICE row.
    rows = Row[]
    for ab in ("NONE", "LT", "LT+S", "CN+S")
        g = "spkezr $ab"
        push!(rows, bench_row(g, "Almanac spkezr, string correction", loop_neo_ezr, ctx, epochs, ab; n))
        push!(rows, bench_row(g, "Almanac spkezr, Aberration", loop_neo_ezr, ctx, epochs, Aberration(ab); n))
        push!(rows, bench_row(g, "CSPICE spkezr (SPICE.jl)", loop_cspice_ezr, epochs, ab; n))
    end
    for ab in ("LT+S", "CN+S")
        g = "spkpos $ab"
        push!(rows, bench_row(g, "Almanac spkpos", loop_neo_pos, ctx, epochs, ab; n))
        push!(rows, bench_row(g, "CSPICE spkpos (SPICE.jl)", loop_cspice_pos, epochs, ab; n))
    end
    # Jets have no CSPICE counterpart; order 1 is the state of spkezr.
    for ab in ("LT+S", "CN+S"), N in 1:3
        push!(rows, bench_row("spkder $ab", "Almanac spkder Val($N)", loop_neo_der, ctx, epochs, ab, Val(N); n))
    end

    text = run_header("Aberration benchmark: Almanac, CSPICE") * """
        Kernel: `$(Kernels.DE440S)` only, in both engines.
        Epochs: 100, evenly spaced from 1900-01-01 12:00 ET to 2099-01-01 12:00 ET (end excluded).

        ## Time per query

        """ * markdown_table(rows)
    println()
    println(text)
    write_results("aberration", text, rows)
    SPICE.kclear()
    return rows
end

main()
