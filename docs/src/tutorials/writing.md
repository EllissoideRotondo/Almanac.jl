# Writing kernels

Almanac.jl writes SPK, CK, binary PCK, DSK and EK files. From the same input, a file holds the data, word for word, of the file that CSPICE writes in a new process, and its bytes differ only where CSPICE leaves them unset, such as the padding of DSK and EK records. Every writer follows the pattern open, write, close; this page writes an SPK, reads it back, and builds kernels in memory.

## Write an SPK

[`spkopn`](@ref) creates a file, which must not exist, and returns a handle. Each `spkwNN` adds a segment of type `NN`, and [`spkcls`](@ref) closes the file. A file still open for write can be loaded, and a context sees only the segments written before it loaded the file.

Type 5 interpolates states with two-body motion, so a few states per revolution are enough. [`conics`](@ref) computes them from orbital elements.

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

As in CSPICE, the file changes on disk at each write call, except for the writers that collect records and write the segment at the end: SPK type 14 ([`spk14b`](@ref), [`spk14a`](@ref), [`spk14e`](@ref)), CK type 4 ([`ckw04b`](@ref), [`ckw04a`](@ref), [`ckw04e`](@ref)) and binary PCK type 3 ([`pck03b`](@ref), [`pck03a`](@ref), [`pck03e`](@ref)).

## Inspect the file

[`spkobj`](@ref) and [`spkcov`](@ref) read the bodies and the coverage of a file, which need not be loaded.

```@example write
spkcov(path, -77)
```

SPK, CK and binary PCK files are DAF (double precision array) files, and the DAF calls reach their segments directly. The summary of an SPK segment holds 2 doubles (the coverage) and 6 integers (body, center, frame, type and two addresses).

```@example write
daf = dafopr(path)
dafbfs(daf)                                  # begin a forward search
while daffna(daf)                            # move to the next segment
    doubles, integers = dafus(dafgs(daf), 2, 6)
    println(dafgn(daf), ": body ", integers[1], ", center ", integers[2], ", type ", integers[4])
end
dafcls(daf)
```

## Load the file

A body without a name is given by its ID code, as an integer or as a string. The epoch is a number, so no leap-seconds kernel is needed.

```@example write
furnsh(path)
state, _ = spkezr("-77", 3600.0, "J2000", "NONE", "399")
kclear()
state[1:3]                                   # km
```

## Kernels in memory

A kernel can also live in memory, with no file. This suits states that a program computes (a numerical propagator) or that exist as samples: unlike a body of [`add_ephemeris!`](@ref), whose function runs at every query, a kernel in memory costs at each query what an SPK file costs. The examples use a circular orbit of radius 8000 km around the Earth:

```@example write
r = 8000.0                                   # km
n = sqrt(mu / r^3)                           # rad/s
circular(et) = [r * cos(n * et), r * sin(n * et), 0.0,
                -r * n * sin(n * et), r * n * cos(n * et), 0.0]
nothing # hide
```

[`load_samples!`](@ref) loads states at epochs, which queries interpolate with SPK type 13 (Hermite, which uses the velocities) or type 9 (Lagrange):

```@example write
ets = collect(0.0:60.0:3600.0)
load_samples!("samples", -1001, "EARTH", "J2000", ets, circular.(ets))
vnorm(spkpos(-1001, 1234.5, "J2000", "NONE", "EARTH")[1] - circular(1234.5)[1:3])  # km
```

[`load_function!`](@ref) fits any function with Chebyshev polynomials to a tolerance and returns the measured errors in position (km) and velocity (km/s). The function is called only during the fit:

```@example write
poserr, velerr = load_function!("fit", circular, -1002, "EARTH", "J2000", 0.0, 86400.0; tol = 1e-6)
```

An `ODESolution` of DifferentialEquations.jl is such a function, through a closure such as `et -> sol(et - et0)[1:6]`. Check three things in it: the time of the solver is usually counted from a start epoch `et0`, not from J2000; the state may hold more than six components (a mass, a transition matrix); and the units must be km and km/s.

For several arcs (a trajectory with burns) or another SPK type, write the segments yourself into an SPK in memory, from [`spkops`](@ref), with any writer, and load it with [`furnsh`](@ref). A kernel can bring names of bodies, which go away with it; [`save_kernel`](@ref) writes the same bytes to a file:

```@example write
h = spkops()
half = length(ets) ÷ 2
spkw13(h, -1003, 399, "J2000", ets[1], ets[half], "arc 1", 7, circular.(ets[1:half]), ets[1:half])
spkw13(h, -1003, 399, "J2000", ets[half], ets[end], "arc 2", 7, circular.(ets[half:end]), ets[half:end])
furnsh(h; name = "arcs", bodies = ["MY PROBE" => -1003])
save_kernel(h, joinpath(mktempdir(), "arcs.bsp"))
spkcls(h)
bodn2c("MY PROBE")
```

The load takes a copy, so close the handle after `furnsh` or `save_kernel`, not before.

## Writers for each kind of kernel

| Kernel | Open | Write | Close |
|---|---|---|---|
| SPK, 16 data types | [`spkopn`](@ref), or [`spkopa`](@ref) to append | `spkwNN` for each type, [`spk14b`](@ref) for type 14 | [`spkcls`](@ref) |
| SPK from two-line elements | [`spkopn`](@ref) | [`getelm`](@ref), [`spkw10`](@ref) | [`spkcls`](@ref) |
| CK, types 1 to 6 | [`ckopn`](@ref) | `ckw01`, `ckw02`, `ckw03`, `ckw05`, `ckw06`, [`ckw04b`](@ref) for type 4 | [`ckcls`](@ref) |
| Binary PCK, types 2, 3 and 20 | [`pckopn`](@ref) | [`pckw02`](@ref), [`pckw20`](@ref), [`pck03b`](@ref) for type 3 | [`pckcls`](@ref) |
| DSK type 2 | [`dskopn`](@ref) | [`dskmi2`](@ref), [`dskw02`](@ref) | [`dskcls`](@ref) |
| EK | [`ekopn`](@ref), [`ekops`](@ref) | [`ekbseg`](@ref), [`ekappr`](@ref), `ekace*`, or whole columns with [`ekifld`](@ref) | [`ekcls`](@ref) |
| Any DAF, or any DAS (the container of DSK and EK files) | [`dafonw`](@ref), [`dasonw`](@ref) | `dafada`, `dasadd` and related calls | [`dafcls`](@ref), [`dascls`](@ref) |

A text kernel is plain text: write it with `write`, or set the pool variables directly with [`pdpool`](@ref), [`pcpool`](@ref) and [`pipool`](@ref).
