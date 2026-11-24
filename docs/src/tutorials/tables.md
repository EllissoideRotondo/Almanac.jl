# Ephemeris and rotation tables

An inner loop or a GPU kernel that needs the position of a body many times over a known interval does not have to touch the kernels. [`ephemeris_table`](@ref) extracts the position into an [`EphemerisTable`](@ref): one flat array of Chebyshev coefficients on equal intervals, evaluated without loaded kernels, segment search or chain of bodies.

## Build and evaluate an ephemeris table

```@example tab
using Almanac, SPICEKernels                 # SPICEKernels downloads the NAIF generic kernels

furnsh(naif0012_tls())                      # leap seconds
furnsh(de440s())                            # planets, Moon and Sun, 1849 to 2150
furnsh(pck00011())                          # orientation and radii of the bodies

t0 = str2et("2025 MAR 01 12:00:00 UTC")
t1 = t0 + 86400.0
tab = ephemeris_table("MOON", t0, t1, "J2000", "NONE", "EARTH"; tol=1e-6)
(tab.nint, tab.ncoef, tab.maxerr)           # intervals, coefficients per interval, measured error (km)
```

`tol` is the accepted error in km; its default is 1e-12 times the largest distance in the table. [`table_position`](@ref), [`table_state`](@ref) and [`table_jet`](@ref) return the position, the state, and the jet: the position with its time derivatives of orders 0 to `N`.

```@example tab
et = t0 + 3600.0
table_state(tab, et) - spkezr("MOON", et, "J2000", "NONE", "EARTH")[1]
```

```@example tab
position, velocity, acceleration = table_jet(tab, et, Val(2))
acceleration                                # km/s²
```

Epochs can be `ForwardDiff.Dual` numbers, also on the GPU.

An epoch outside the table throws a [`SpiceKindError`](@ref) of kind `:spkinsuffdata`, without the fields of the [`SpkInsufficientDataError`](@ref) that `spkezr` throws. The check is a bounds check, which `@inbounds` removes, as GPU code requires; an epoch outside the table then extrapolates the nearest interval.

## How the table is built

For a geometric position in an inertial frame from SPK type 2 data, as in the planetary ephemerides, the table copies the [records](@ref "Smoothness at record boundaries") of the kernel and agrees with `spkpos` to rounding. Otherwise (a rotating frame, an aberration correction, other SPK types) it interpolates, doubling the number of intervals until the measured error is below `tol`, and `maxerr` records that error. [`ephemeris_table`](@ref) gives the exact conditions.

The velocity of a table is always the derivative of its position, so it can differ from the velocity of `spkezr` with an `LT` correction, or for SPK types that store the velocity separately.

## Rotation tables

[`rotation_table`](@ref) does for `pxform` what `ephemeris_table` does for `spkpos`, by interpolating the nine elements of the matrix. [`table_rotation`](@ref) returns the matrix, and [`table_rotation_jet`](@ref) the matrix and its derivatives.

```@example tab
rtab = rotation_table("J2000", "IAU_EARTH", t0, t1)
maximum(abs, table_rotation(rtab, et) - pxform("J2000", "IAU_EARTH", et))
```

The matrix is not orthogonalized again, so its distance from orthogonality is of the order of the error of the table. The errors below are Frobenius norms over six days from 2022-03-08, with the default arguments; `MOON_PA` is the principal-axis frame of the Moon:

| From | To | Intervals | Measured error of the table | Difference from CSPICE `pxform` |
|---|---|---:|---:|---:|
| `ECLIPJ2000` | `J2000` | 1 | 0 | 0 |
| `IAU_MOON` | `J2000` | 1 | 5.8e-13 | 7.5e-13 |
| `J2000` | `IAU_EARTH` | 16 | 1.5e-11 | 2.0e-11 |
| `MOON_PA` | `IAU_EARTH` | 16 | 1.6e-11 | 1.9e-11 |

The tables that reach the Earth stop at 1.5e-11 because of the rounding of the rotation angle of the Earth; see [Known limits of accuracy](@ref).

## Measured accuracy and speed

The figures below cover 400 days of `de440s.bsp` and 4,000,000 epochs, measured on one machine, on one CPU thread and on a GPU. The error is the largest difference from the kernels, and the tolerance is the default unless stated.

| Case | Intervals | Error (km) | Table (ns per epoch) | `spkezr` (ns per epoch) |
|---|---:|---:|---:|---:|
| Moon relative to Earth, J2000, `NONE`, CPU | 101 | 2.4e-10 | 36 | 152 |
| Earth relative to Sun, ECLIPJ2000, `NONE`, CPU | 101 | 1.2e-7 | 36 | 277 |
| Moon relative to Earth, J2000, `LT+S`, CPU | 128 | 4.3e-6 | 42 | 653 |
| Earth relative to Moon, MOON_PA, `NONE`, `tol=1e-6`, CPU | 1024 | 4.2e-7 | 42 | 385 |
| Sun relative to Earth, IAU_EARTH, `NONE`, CPU | 1024 | 3.3e-3 | 42 | 341 |
| Earth relative to Sun, J2000, GPU, `Float64` | 101 | 1.2e-7 | 1.2 | |
| Earth relative to Sun, J2000, GPU, `Float32` | 101 | 16 | 0.33 | |

A table is no more accurate than its source: the `IAU_EARTH` and `LT+S` rows reach the limits in [Known limits of accuracy](@ref).

## Save a table as a kernel

[`table_to_spk`](@ref) writes a geometric table to an SPK file as one segment of type 2. The `spkezr` of Almanac.jl or of CSPICE then returns the state of the table to rounding: 2e-16 relative, measured for the Moon seen from the Earth over 200 days.

```@example tab
path = joinpath(mktempdir(), "table.bsp")
handle = spkopn(path, "from a table", 0)
table_to_spk(handle, tab, "Moon from a table")
spkcls(handle)
spkobj(path)
```

## Evaluate on a GPU

The table is parametric in its array type, so Adapt.jl moves it to a device, and its evaluation is plain arithmetic that compiles in a GPU kernel. Epochs and bounds stay `Float64` even with `Float32` coefficients, because a `Float32` count of seconds past J2000 resolves only about a minute. Below, `epochs` is a vector of `Float64`.

```julia
using CUDA, Adapt

gpu = adapt(CuArray, tab)                   # `cu(tab)` gives Float32 coefficients
states = map(et -> @inbounds(table_state(gpu, et)), CuArray(epochs))
```
