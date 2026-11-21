# Group: kernels built in memory (`dafops`, `spkops`, `ckops`, `pckops`, `furnsh` of a
# handle, names of a kernel, `load_samples!`, `load_function!`; src/daf_write.jl,
# src/memkernels.jl).
#
# Invariant layer: the writers write the same bytes to a handle in memory as to a file
# (group `dafwrite` compares the files with CSPICE), so a kernel loaded from a handle
# must give the queries of the same file to the last bit. No kernel file of NAIF is
# needed.

FULL || @eval Base.Experimental.@optlevel 0

const MEM_GM = 398600.4418                  # km^3/s^2, a Kepler orbit about the Earth
const MEM_PV0 = [7000.0, 0.0, 0.0, 0.0, 7.5, 1.0]

"States of the Kepler orbit at `ets` (km, km/s)."
mem_states(ets) = [prop2b(MEM_GM, MEM_PV0, et) for et in ets]

"Write SPK segments of types 2, 9 and 13 into `h`, and the Earth relative to the barycenter."
function mem_spk!(h)
    ets = collect(range(-3600.0, 3600.0; length = 25))
    states = mem_states(ets)
    spkw02(h, 499, 4, "J2000", 0.0, 100.0, "type 2", 100.0, [[1.0, 0.5, 2.0, 0.0, 3.0, 0.0]], 0.0)
    spkw09(h, -1000, 399, "J2000", ets[1], ets[end], "type 9", 7, states, ets)
    spkw13(h, -1001, 399, "J2000", ets[1], ets[end], "type 13", 7, states, ets)
    earth = [[1.5e8 - 30.0 * et, 30.0 * et, 0.0, -30.0, 30.0, 0.0] for et in (ets[1], ets[end])]
    spkw09(h, 399, 0, "J2000", ets[1], ets[end], "earth", 1, earth, [ets[1], ets[end]])
    return h
end

"Write a CK segment of type 3 into `h`."
function mem_ck!(h)
    quats = [[1.0, 0.0, 0.0, 0.0], [cos(0.1), sin(0.1), 0.0, 0.0]]
    ckw03(h, 0.0, 10.0, -77000, "J2000", "attitude", [0.0, 10.0], quats, nothing, [0.0])
    return h
end

"Write a binary PCK segment of type 2 into `h`."
function mem_pck!(h)
    pckw02(h, 3000, "J2000", 0.0, 100.0, "angles", 100.0, [[0.1, 0.01, 0.2, 0.0, 0.3, 0.0]], 0.0)
    return h
end

@testset "in-memory handles" begin
    dir = mktempdir()
    @testset "bytes equal a file" begin
        cases = ((spkopn, spkops, mem_spk!), (ckopn, ckops, mem_ck!), (pckopn, pckops, mem_pck!))
        for (opn, ops, write!) in cases
            path = joinpath(dir, string(ops, ".bin"))
            f, m = opn(path, "same", 2000), ops("same", 2000)
            for h in (f, m)
                dafac(h, ["first line", "second line"])
                write!(h)
            end
            @test Almanac._daf_bytes(m) == read(path)
            close(f)
            close(m)
        end
    end

    @testset "read calls" begin
        path = joinpath(dir, "read.bsp")
        spkcls(mem_spk!(spkopn(path, "x")))
        r, m = dafopr(path), mem_spk!(spkops("x"))
        dafbfs(r)
        dafbfs(m)
        n = 0
        while daffna(r)              # invariant: the first n arrays of both agreed
            @test daffna(m)
            @test dafgs(m) == dafgs(r)
            @test dafgn(m) == dafgn(r)
            _, ic = dafus(dafgs(m), 2, 6)
            @test dafgda(m, ic[5], ic[6]) == dafgda(r, ic[5], ic[6])
            n += 1
        end
        @test n == 4 && !daffna(m)
        @test dafec(m) == dafec(r)
        dafcls(r)
        spkcls(m)
    end

    @testset "close and the table of open files" begin
        a, b = spkops(), spkops()
        @test !haskey(Almanac._OPEN_FILES, a)
        @test isopen(a) && isopen(b)
        @test throws_kind(() -> spkcls(a), :nosegmentsfound)
        mem_spk!(b)
        spkcls(b)
        @test !isopen(b)
        @test throws_kind(() -> dafgda(b, 1, 1), :dafnosuchhandle)
    end
end

@testset "furnsh of a handle" begin
    dir = mktempdir()
    path = joinpath(dir, "same.bsp")
    spkcls(mem_spk!(spkopn(path, "same")))
    oracle_load([path])
    ets = collect(range(-3500.0, 3500.0; length = 50))

    @testset "queries equal the file and CSPICE" begin
        h = mem_spk!(spkops("same"))
        mem, file = SpiceContext(), neo_context([path])
        furnsh(mem, h; name = "same")
        spkcls(h)
        for body in (-1000, -1001), corr in ("NONE", "LT+S"), et in ets
            s, lt = spkezr(mem, body, et, "J2000", corr, 399)
            @test (s, lt) == spkezr(file, body, et, "J2000", corr, 399)
            @test s == SPICE.spkezr(string(body), et, "J2000", corr, "399")[1]
        end
        @test spkpos(mem, 499, 50.0, "J2000", "NONE", 4)[1] == spkpos(file, 499, 50.0, "J2000", "NONE", 4)[1]
    end

    @testset "load list" begin
        h = mem_spk!(spkops())
        ctx = SpiceContext()
        furnsh(ctx, h; name = "mine")
        @test kdata(ctx, 1) == ("mine", "SPK", "", 0)
        # A file loaded later has priority; the kernel in memory has priority over an
        # earlier one. The file `line.bsp` puts body -1000 on a straight line.
        line = joinpath(dir, "line.bsp")
        let f = spkopn(line, "line")
            spkw09(f, -1000, 399, "J2000", -3600.0, 3600.0, "line", 1, [zeros(6), [1.0, 0, 0, 0, 0, 0]], [-3600.0, 3600.0])
            spkcls(f)
        end
        furnsh(ctx, line)
        @test spkpos(ctx, -1000, 0.0, "J2000", "NONE", 399)[1] == [0.5, 0.0, 0.0]
        furnsh(ctx, h; name = "mine")
        @test ktotal(ctx) == 3
        @test spkpos(ctx, -1000, 0.0, "J2000", "NONE", 399)[1] != [0.5, 0.0, 0.0]
        unload(ctx, "mine")
        @test ktotal(ctx) == 2
        @test spkpos(ctx, -1000, 0.0, "J2000", "NONE", 399)[1] == [0.5, 0.0, 0.0]
        kclear(ctx)
        @test ktotal(ctx) == 0
        spkcls(h)
    end

    @testset "a load copies the segments" begin
        h = spkops()
        spkw09(h, -1000, 399, "J2000", -10.0, 10.0, "one", 1, [zeros(6), [1.0, 0, 0, 0, 0, 0]], [-10.0, 10.0])
        ctx = SpiceContext()
        furnsh(ctx, h; name = "k")
        spkw09(h, -1000, 399, "J2000", -10.0, 10.0, "two", 1, [zeros(6), [3.0, 0, 0, 0, 0, 0]], [-10.0, 10.0])
        @test spkpos(ctx, -1000, 0.0, "J2000", "NONE", 399)[1] == [0.5, 0.0, 0.0]
        furnsh(ctx, h; name = "k")
        @test spkpos(ctx, -1000, 0.0, "J2000", "NONE", 399)[1] == [1.5, 0.0, 0.0]
        spkcls(h)
    end

    @testset "errors" begin
        ctx = SpiceContext()
        closed = mem_spk!(spkops())
        spkcls(closed)
        onfile = spkopn(joinpath(dir, "open.bsp"), "x")
        adding = spkops()
        dafbna(adding, dafps([0.0, 1.0], [1, 2, 3, 4, 0, 0]), "open array")
        empty = spkops()
        cases = (
            (closed, "a", :dafnosuchhandle), (onfile, "a", :dafinvalidaccess),
            (adding, "a", :dafnewconflict), (empty, "a", :nosegmentsfound),
            (mem_spk!(spkops()), "", :emptystring),
        )
        for (h, name, kind) in cases
            @test throws_kind(() -> furnsh(ctx, h; name = name), kind)
            @test ktotal(ctx) == 0
        end
        close(onfile)
    end

    @testset "CK and binary PCK" begin
        ck, pck = joinpath(dir, "same.bc"), joinpath(dir, "same.bpc")
        ckcls(mem_ck!(ckopn(ck, "x")))
        pckcls(mem_pck!(pckopn(pck, "x")))
        mem, file = SpiceContext(), neo_context([ck, pck])
        hc, hp = mem_ck!(ckops()), mem_pck!(pckops())
        furnsh(mem, hc; name = "ck")
        furnsh(mem, hp; name = "pck")
        @test kdata(mem, 1)[2] == "CK" && kdata(mem, 2)[2] == "PCK"
        for t in (0.0, 2.5, 7.0)
            @test ckgp(mem, -77000, t, 0.0, "J2000") == ckgp(file, -77000, t, 0.0, "J2000")
            @test sxform(mem, "ITRF93", "J2000", 10t) == sxform(file, "ITRF93", "J2000", 10t)
        end
        ckcls(hc)
        pckcls(hp)
    end
end

@testset "names of a kernel" begin
    # A context with the kernel of `mem_spk!` loaded under `name` with `bodies`.
    function named!(ctx, name, bodies)
        h = mem_spk!(spkops())
        furnsh(ctx, h; name = name, bodies = bodies)
        spkcls(h)
        return ctx
    end

    @testset "the name comes and goes with the kernel" begin
        ctx = named!(SpiceContext(), "k", ["MySat" => -1001])
        @test bodn2c(ctx, "MYSAT") == -1001
        @test bodc2n(ctx, -1001) == "MySat"
        @test spkezr(ctx, "mysat", 0.0, "J2000", "NONE", "EARTH")[1] ==
            spkezr(ctx, -1001, 0.0, "J2000", "NONE", 399)[1]
        unload(ctx, "k")
        @test bodn2c(ctx, "MYSAT") === nothing
        @test throws_kind(() -> furnsh(ctx, spkops(); name = "x", bodies = [" " => 1]), :badarg)
        @test throws_kind(() -> furnsh(ctx, spkops(); name = "x", bodies = ["BIG" => 2^40]), :badarg)
    end

    @testset "order: pool, later kernel, earlier kernel, boddef" begin
        ctx = SpiceContext()
        boddef(ctx, "MYSAT", -3000)
        named!(ctx, "early", ["MYSAT" => -1001])
        @test bodn2c(ctx, "MYSAT") == -1001                  # a kernel over boddef
        named!(ctx, "late", ["MYSAT" => -1002])
        @test bodn2c(ctx, "MYSAT") == -1002                  # the later kernel wins
        pcpool(ctx, "NAIF_BODY_NAME", ["MYSAT"])
        pipool(ctx, "NAIF_BODY_CODE", [-2000])
        @test bodn2c(ctx, "MYSAT") == -2000                  # the pool over the kernels
        dvpool(ctx, "NAIF_BODY_NAME")
        dvpool(ctx, "NAIF_BODY_CODE")
        unload(ctx, "late")
        @test bodn2c(ctx, "MYSAT") == -1001                  # the earlier kernel again
        unload(ctx, "early")
        @test bodn2c(ctx, "MYSAT") == -3000                  # boddef again
    end

    @testset "add_ephemeris! with names" begin
        ctx = SpiceContext()
        add_ephemeris!(et -> [7000.0, 0.001et, 0, 0, 0.001, 0], ctx, "PROBE", "EARTH", "J2000", 0.0, 10.0; name = "p", bodies = ["PROBE" => -1003])
        @test spkpos(ctx, "PROBE", 5.0, "J2000", "NONE", "EARTH")[1] == [7000.0, 0.005, 0.0]
        unload(ctx, "p")
        @test bodn2c(ctx, "PROBE") === nothing
    end

    @testset "a dynamic frame that names the body" begin
        ctx = SpiceContext()
        lmpool(
            ctx, [
                "FRAME_MYSAT_RTN = 1850009", "FRAME_1850009_NAME = 'MYSAT_RTN'",
                "FRAME_1850009_CLASS = 5", "FRAME_1850009_CLASS_ID = 1850009",
                "FRAME_1850009_CENTER = 399", "FRAME_1850009_RELATIVE = 'J2000'",
                "FRAME_1850009_DEF_STYLE = 'PARAMETERIZED'", "FRAME_1850009_FAMILY = 'TWO-VECTOR'",
                "FRAME_1850009_PRI_AXIS = 'X'",
                "FRAME_1850009_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'",
                "FRAME_1850009_PRI_OBSERVER = 'EARTH'", "FRAME_1850009_PRI_TARGET = 'MYSAT'",
                "FRAME_1850009_PRI_ABCORR = 'NONE'", "FRAME_1850009_SEC_AXIS = 'Y'",
                "FRAME_1850009_SEC_VECTOR_DEF = 'OBSERVER_TARGET_VELOCITY'",
                "FRAME_1850009_SEC_OBSERVER = 'EARTH'", "FRAME_1850009_SEC_TARGET = 'MYSAT'",
                "FRAME_1850009_SEC_ABCORR = 'NONE'", "FRAME_1850009_SEC_FRAME = 'J2000'",
            ]
        )
        named!(ctx, "k", ["MYSAT" => -1001])
        s = spkezr(ctx, -1001, 100.0, "J2000", "NONE", 399)[1]
        @test pxform(ctx, "J2000", "MYSAT_RTN", 100.0) ≈ twovec(s[1:3], 1, s[4:6], 2) rtol = 1.0e-14
    end
end

@testset "coverage of a handle, save_kernel" begin
    dir = mktempdir()
    spk, ck, pck = joinpath(dir, "c.bsp"), joinpath(dir, "c.bc"), joinpath(dir, "c.bpc")
    spkcls(mem_spk!(spkopn(spk, "c")))
    ckcls(mem_ck!(ckopn(ck, "c")))
    pckcls(mem_pck!(pckopn(pck, "c")))
    hs, hc, hp = mem_spk!(spkops("c")), mem_ck!(ckops("c")), mem_pck!(pckops("c"))

    @testset "the same as the file" begin
        @test spkobj(hs) == spkobj(spk)
        for id in spkobj(spk)
            @test collect(spkcov(hs, id)) == collect(spkcov(spk, id))
        end
        @test ckobj(hc) == ckobj(ck)
        for level in ("SEGMENT", "INTERVAL")
            @test collect(ckcov(hc, -77000, false, level, 0.0, "SCLK")) ==
                collect(ckcov(ck, -77000, false, level, 0.0, "SCLK"))
        end
        @test pckfrm(hp) == pckfrm(pck)
        @test collect(pckcov(hp, 3000)) == collect(pckcov(pck, 3000))
        let r = dafopr(spk)                 # a handle of a file works too
            @test spkobj(r) == spkobj(spk)
            dafcls(r)
        end
        @test throws_kind(() -> spkobj(hc), :invalidfiletype)
    end

    @testset "a search under way goes on" begin
        dafbfs(hs)
        @test daffna(hs)
        first_name = dafgn(hs)
        spkcov(hs, -1000)
        @test dafgn(hs) == first_name
        @test daffna(hs) && dafgn(hs) != first_name
    end

    @testset "save_kernel" begin
        out = joinpath(dir, "saved.bsp")
        save_kernel(hs, out)
        @test read(out) == read(spk)
        @test throws_kind(() -> save_kernel(hs, out), :fileopenfail)
        r = dafopr(spk)
        @test throws_kind(() -> save_kernel(r, joinpath(dir, "x.bsp")), :dafinvalidaccess)
        dafcls(r)
    end
    spkcls(hs)
    ckcls(hc)
    pckcls(hp)
end

@testset "load_samples!" begin
    dir = mktempdir()
    n = nsample(200, 100_000)
    ets = collect(range(-3600.0, 3600.0; length = n))
    states = mem_states(ets)
    probe = sort!(rand(MersenneTwister(SEED + 91), n ÷ 4) .* 7200.0 .- 3600.0)
    for (method, writer) in ((:hermite, spkw13), (:lagrange, spkw09))
        path = joinpath(dir, "$method.bsp")
        let f = spkopn(path, "x")
            writer(f, -1000, 399, "J2000", ets[1], ets[end], "samples", 7, states, ets)
            spkcls(f)
        end
        file = neo_context([path])
        for form in (states, reduce(hcat, states), reduce(vcat, states))
            ctx = SpiceContext()
            load_samples!(ctx, "samples", -1000, 399, "J2000", ets, form; method = method)
            @test all(et -> spkezr(ctx, -1000, et, "J2000", "NONE", 399) == spkezr(file, -1000, et, "J2000", "NONE", 399), probe)
        end
    end

    @testset "segment name, names, errors" begin
        ctx = SpiceContext()
        load_samples!(ctx, "x"^50, "PROBE", "EARTH", "J2000", ets, states; bodies = ["PROBE" => -1000])
        @test bodn2c(ctx, "PROBE") == -1000
        @test kdata(ctx, 1)[1] == "x"^50
        cases = (
            ((; method = :spline), :invalidargument), ((; degree = 6), :invaliddegree),
            ((; segid = "a\tb"), :nonprintablechars),
        )
        for (kw, kind) in cases
            @test throws_kind(() -> load_samples!(ctx, "bad", -1000, 399, "J2000", ets, states; kw...), kind)
            @test ktotal(ctx) == 1
        end
        @test throws_kind(() -> load_samples!(ctx, "bad", -1000, 399, "J2000", Float64[], states), :toofewstates)
        @test throws_kind(() -> load_samples!(ctx, "bad", "NOBODY", 399, "J2000", ets, states), :idcodenotfound)
        @test ktotal(ctx) == 1
    end
end

@testset "load_function!" begin
    f(et) = prop2b(MEM_GM, MEM_PV0, et)
    probe = rand(MersenneTwister(SEED + 92), nsample(200, 20_000)) .* 86400.0
    # The fit measures its error at the n - 1 inner extremes of each interval, which bound
    # the interpolation error to a small factor: 2 is the bound checked here.
    @testset "position and velocity (type 3)" begin
        ctx = SpiceContext()
        poserr, velerr = load_function!(ctx, "orbit", f, -1000, 399, "J2000", 0.0, 86400.0; tol = 1.0e-6)
        @test poserr <= 1.0e-6 && 0 < velerr < 1.0e-6
        @test spkuds(spksfs(ctx, -1000, 10.0)[2])[4] == 3
        dp = maximum(et -> vnorm(spkezr(ctx, -1000, et, "J2000", "NONE", 399)[1][1:3] - f(et)[1:3]), probe)
        dv = maximum(et -> vnorm(spkezr(ctx, -1000, et, "J2000", "NONE", 399)[1][4:6] - f(et)[4:6]), probe)
        @test dp <= 2poserr && dv <= 2velerr
        measured("load_function! type 3: poserr=$poserr dp=$dp velerr=$velerr dv=$dv")
    end
    @testset "position only (type 2)" begin
        ctx = SpiceContext()
        poserr, velerr = load_function!(ctx, "orbit", et -> f(et)[1:3], -1000, 399, "J2000", 0.0, 86400.0; tol = 1.0e-6)
        @test poserr <= 1.0e-6 && isnan(velerr)
        @test spkuds(spksfs(ctx, -1000, 10.0)[2])[4] == 2
        @test maximum(et -> vnorm(spkpos(ctx, -1000, et, "J2000", "NONE", 399)[1] - f(et)[1:3]), probe) <= 2poserr
    end
    @testset "errors" begin
        ctx = SpiceContext()
        cases = (
            (et -> f(et)[1:5], :badarraysize), (et -> et > 100 ? zeros(2) : f(et), :badarraysize),
            (et -> throw(SpiceError(:notcomputable, "no")), :notcomputable),
        )
        for (g, kind) in cases
            @test throws_kind(() -> load_function!(ctx, "bad", g, -1000, 399, "J2000", 0.0, 86400.0), kind)
        end
        for kw in ((; ncoef = 29), (; ncoef = 1), (; tol = 0.0), (; maxintervals = 0))
            @test throws_kind(() -> load_function!(ctx, "bad", f, -1000, 399, "J2000", 0.0, 86400.0; kw...), :badarg)
        end
        @test throws_kind(() -> load_function!(ctx, "bad", f, -1000, 399, "J2000", 0.0, 86400.0; tol = 1.0e-12, maxintervals = 2), :noconvergence)
        @test ktotal(ctx) == 0
    end
end

@testset "handles of other types, closed handles" begin
    ctx = SpiceContext()
    @test throws_kind(() -> furnsh(ctx, dafops("XYZ", 2, 6, "x", 0); name = "x"), :baddaf)
    @test ktotal(ctx) == 0
    h = mem_spk!(spkops())
    spkcls(h)
    @test throws_kind(() -> spkobj(h), :dafnosuchhandle)
    @test throws_kind(() -> pckcov(h, 3000), :dafnosuchhandle)
end
