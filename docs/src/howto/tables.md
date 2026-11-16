# Ephemeris and rotation tables

An inner loop or a GPU kernel that needs a body's position many times over a known interval need not touch the kernels. [`ephemeris_table`](@ref) extracts the position into an [`EphemerisTable`](@ref): one flat array of Chebyshev coefficients on equal intervals, evaluated without loaded kernels, segment search or body chain.

```@setup tab
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Build and evaluate

```@example tab
using Almanac

t0 = str2et("2025 MAR 01 12:00:00 UTC")
t1 = t0 + 86400.0
tab = ephemeris_table("PROBE", t0, t1, "J2000", "NONE", "EARTH"; tol=1e-6)
(tab.nint, tab.ncoef, tab.maxerr)           # intervals, coefficients per interval, measured error (km)
```

`tol` is the accepted error in km; the default is 1e-12 times the largest distance in the table. [`table_position`](@ref), [`table_state`](@ref) and [`table_jet`](@ref) return the position, the state, and the position with derivatives of any order.

```@example tab
et = t0 + 3600.0
table_state(tab, et) - spkezr("PROBE", et, "J2000", "NONE", "EARTH")[1]
```

```@example tab
position, velocity, acceleration = table_jet(tab, et, Val(2))
acceleration                                # km/s²
```

## How the table is built

For a geometric position in an inertial frame from SPK type 2 data (the planetary ephemerides), the table copies the [records](@ref "Exact derivatives and the rates of CSPICE") of the kernel and agrees with `spkpos` to rounding. Otherwise (rotating frame, aberration correction, other SPK types) it interpolates, doubling the intervals until the measured error is below `tol`; `maxerr` records that error. [`ephemeris_table`](@ref) gives the exact conditions.

The table velocity is always the derivative of its position, so it can differ from the `spkezr` velocity with an `LT` correction or for SPK types that store velocity separately. Epochs can be `ForwardDiff.Dual`, also on the GPU.

An epoch outside the table throws `SpiceError(:spkinsuffdata)`, not the [`SpkInsufficientDataError`](@ref) of `spkezr`. The check is a bounds check: `@inbounds` removes it, as GPU code requires, and an epoch outside then extrapolates the nearest interval.

## On a GPU

The table is parametric in its array type, so Adapt.jl moves it to a device, and evaluation is plain arithmetic that compiles in a GPU kernel. Epochs and bounds stay `Float64` even with `Float32` coefficients: a `Float32` count of seconds past J2000 resolves only about a minute. Below, `epochs` is a vector of `Float64`.

```julia
using CUDA, Adapt

gpu = adapt(CuArray, tab)                   # `cu(tab)` gives Float32 coefficients
states = map(et -> @inbounds(table_state(gpu, et)), CuArray(epochs))
```

## In a CasADi model

[CasADiCodegen.jl](https://github.com/EllissoideRotondo/CasADiCodegen.jl) generates Julia code from a CasADi model; its `callback` embeds a Julia function with the derivatives CasADi needs. Since [`table_jet`](@ref) returns static vectors, a table gives the position of a body to generated code without allocation, and its derivatives to the Jacobian. A sketch (not a complete program) for an Earth satellite perturbed by the Moon, where `t0`, `t1`, `GM_EARTH` and `GM_MOON` come from the application:

```julia
using Almanac, CasADiCodegen, LinearAlgebra

tab = ephemeris_table("MOON", t0, t1, "J2000", "NONE", "EARTH")
order(k) = let v = Val(k); t -> table_jet(tab, t[1], v)[k+1]; end       # derivative k of the position
moon = callback(order(0), "moon", (1,), (3,); jacobian = (order(1), order(2)))

x, t = MX.sym("x", 6), MX.sym("t")
r, s = x[1:3], moon(t)
a = -GM_EARTH * r / norm(r)^3 + GM_MOON * ((s - r) / norm(s - r)^3 - s / norm(s)^3)
z = vertcat(x, t)
f = CasADiCodegen.Function("f", [x, t], [vertcat(x[4:6], a), jacobian(vertcat(x[4:6], a), z)])
g, g! = build_function(f; expression = Val{false})
```

`t` is in seconds past J2000 TDB and `order(k)` returns the `k`-th derivative in km/sᵏ; CasADi takes first derivatives from `order(1)` and second from `order(2)`.

## Measured accuracy and speed

400 days of de440s.bsp, from `bench/table_gpu.jl`: one thread of an Intel CPU, and an NVIDIA RTX 3090 with 4,000,000 epochs. Error is the largest difference from the kernels; tolerance is the default unless stated. `MOON_PA` is the lunar principal-axis frame.

| Case | Intervals | Error (km) | Table (ns per epoch) | `spkezr` (ns per epoch) |
|---|---:|---:|---:|---:|
| Moon relative to Earth, J2000, `NONE`, CPU | 101 | 2.4e-10 | 36 | 152 |
| Earth relative to Sun, ECLIPJ2000, `NONE`, CPU | 101 | 1.2e-7 | 36 | 277 |
| Moon relative to Earth, J2000, `LT+S`, CPU | 128 | 4.3e-6 | 42 | 653 |
| Earth relative to Moon, MOON_PA, `NONE`, `tol=1e-6`, CPU | 1024 | 4.2e-7 | 42 | 385 |
| Sun relative to Earth, IAU_EARTH, `NONE`, CPU | 1024 | 3.3e-3 | 42 | 341 |
| Earth relative to Sun, J2000, GPU, `Float64` | 101 | 1.2e-7 | 1.2 | |
| Earth relative to Sun, J2000, GPU, `Float32` | 101 | 16 | 0.33 | |

A table is no more accurate than its source; the `IAU_EARTH` and `LT+S` rows reach the limits in [Known limits of accuracy](@ref).

## Save a table as a kernel

[`table_to_spk`](@ref) writes a geometric table to an SPK file as one type 2 segment. `spkezr`, of Almanac.jl or CSPICE, then returns the table's state to rounding (2e-16 relative measured for the Moon from the Earth over 200 days).

```@example tab
path = joinpath(mktempdir(), "table.bsp")
handle = spkopn(path, "from a table", 0)
table_to_spk(handle, tab, "PROBE from a table")
spkcls(handle)
spkobj(path)
```

## Rotations

[`rotation_table`](@ref) does for `pxform` what `ephemeris_table` does for `spkpos`, interpolating the nine matrix elements. [`table_rotation`](@ref) returns the matrix, [`table_rotation_jet`](@ref) the matrix and its derivatives.

```@example tab
rtab = rotation_table("J2000", "IAU_EARTH", t0, t1)
maximum(abs, table_rotation(rtab, et) - pxform("J2000", "IAU_EARTH", et))
```

The matrix is not re-orthogonalized: its distance from orthogonality is of the order of the table error. Errors over six days from 2022-03-08, default arguments, as Frobenius norms:

| From | To | Intervals | Measured error of the table | Difference from CSPICE `pxform` |
|---|---|---:|---:|---:|
| `ECLIPJ2000` | `J2000` | 1 | 0 | 0 |
| `IAU_MOON` | `J2000` | 1 | 5.8e-13 | 7.5e-13 |
| `J2000` | `IAU_EARTH` | 16 | 1.5e-11 | 2.0e-11 |
| `MOON_PA` | `IAU_EARTH` | 16 | 1.6e-11 | 1.9e-11 |

The tables that reach the Earth stop at 1.5e-11 because of the rounding of the Earth's rotation angle (see [Known limits of accuracy](@ref)).
