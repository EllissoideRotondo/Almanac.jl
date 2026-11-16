# Writing kernels

Almanac.jl writes SPK, CK, binary PCK, DSK and EK files with the contents that the CSPICE writers produce from the same input in a [new process](@ref "No state between calls"): SPK, CK and binary PCK files byte for byte, DSK and EK files word for word in their data, since CSPICE fills the padding of a DAS record with what its buffers held before. All writers follow open, write, close; this page writes an SPK and reads it back.

## Write an SPK

[`spkopn`](@ref) creates a file (which must not exist) and returns a handle; each `spkwNN` adds a segment of type `NN`; [`spkcls`](@ref) closes the file. A file still open for write can be loaded, but a context sees only the segments written before it loaded the file.

Type 5 interpolates states with two-body motion, so a few states per revolution suffice. [`conics`](@ref) computes them from orbital elements.

```@example write
using Almanac

mu = 3.986004418e5                           # km³/s²
t0, t1 = 0.0, 86400.0                        # seconds past J2000 TDB

# perifocal distance (km), eccentricity, inclination, longitude of the node, argument
# of periapsis, mean anomaly at the epoch (rad), epoch, gravitational parameter
elements = [7000.0, 0.01, 0.9, 0.0, 0.0, 0.0, t0, mu]
epochs = collect(range(t0, t1; length=97))   # one state each 15 minutes
states = [conics(elements, t) for t in epochs]

path = joinpath(mktempdir(), "example.bsp")
handle = spkopn(path, "example", 0)          # path, internal file name, characters reserved for comments
# handle, body, center, frame, coverage, segment name, mu, states, epochs
spkw05(handle, -77, 399, "J2000", t0, t1, "low orbit", mu, states, epochs)
spkcls(handle)
spkobj(path)                                 # the bodies in the file
```

As in CSPICE, the file changes on disk at each write call, except for type 14, whose writer ([`spk14b`](@ref), [`spk14a`](@ref), [`spk14e`](@ref)) collects records and writes the segment at the end.

## Inspect a file

[`spkobj`](@ref) and [`spkcov`](@ref) read the bodies and coverage of a file that need not be loaded.

```@example write
spkcov(path, -77)
```

SPK, CK and binary PCK files are DAF (double precision array) files, and the DAF calls reach the segments directly. An SPK segment summary holds 2 doubles (coverage) and 6 integers (body, center, frame, type, two addresses).

```@example write
daf = dafopr(path)
dafbfs(daf)                                  # begin a forward search
while daffna(daf)                            # move to the next segment
    doubles, integers = dafus(dafgs(daf), 2, 6)
    println(dafgn(daf), ": body ", integers[1], ", center ", integers[2], ", type ", integers[4])
end
dafcls(daf)
```

## Use the file

A body without a name is given by its ID code, as an integer or a string. The epoch is a number, so no leap-seconds kernel is needed.

```@example write
furnsh(path)
state, _ = spkezr("-77", 3600.0, "J2000", "NONE", "399")
kclear()
state[1:3]                                   # km
```

## Other writers

| Kernel | Open | Write | Close |
|---|---|---|---|
| SPK, 16 data types | [`spkopn`](@ref), or [`spkopa`](@ref) to append | `spkwNN` for each type, [`spk14b`](@ref) for type 14 | [`spkcls`](@ref) |
| SPK from two-line elements | [`spkopn`](@ref) | [`getelm`](@ref), [`spkw10`](@ref) | [`spkcls`](@ref) |
| CK, 5 data types | [`ckopn`](@ref) | `ckw01`, `ckw02`, `ckw03`, `ckw05`, `ckw06` | [`ckcls`](@ref) |
| Binary PCK, types 2 and 20 | [`pckopn`](@ref) | [`pckw02`](@ref), [`pckw20`](@ref) | [`pckcls`](@ref) |
| DSK type 2 | [`dskopn`](@ref) | [`dskmi2`](@ref), [`dskw02`](@ref) | [`dskcls`](@ref) |
| EK | [`ekopn`](@ref), [`ekops`](@ref) | [`ekbseg`](@ref), [`ekappr`](@ref), `ekace*`, or whole columns with [`ekifld`](@ref) | [`ekcls`](@ref) |
| Any DAF, or any DAS (the container of DSK and EK files) | [`dafonw`](@ref), [`dasonw`](@ref) | `dafada`, `dasadd` and related calls | [`dafcls`](@ref), [`dascls`](@ref) |

A text kernel is plain text: write it with `write`, or set pool variables directly with [`pdpool`](@ref), [`pcpool`](@ref) and [`pipool`](@ref).
