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

"Write SPK segments of types 2, 9 and 13 into `h`."
function mem_spk!(h)
    ets = collect(range(-3600.0, 3600.0; length = 25))
    states = mem_states(ets)
    spkw02(h, 499, 4, "J2000", 0.0, 100.0, "type 2", 100.0, [[1.0, 0.5, 2.0, 0.0, 3.0, 0.0]], 0.0)
    spkw09(h, -1000, 399, "J2000", ets[1], ets[end], "type 9", 7, states, ets)
    spkw13(h, -1001, 399, "J2000", ets[1], ets[end], "type 13", 7, states, ets)
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
        @test n == 3 && !daffna(m)
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
