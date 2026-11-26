# Time derivatives

CSPICE gives a position and its velocity. Almanac.jl gives derivatives of any order, either from a direct call or by automatic differentiation of the ordinary calls. The two agree to rounding.

## Differentiate with `spkder` and `pxder`

[`spkder`](@ref) returns a position and its first `N` derivatives, and [`pxder`](@ref) a rotation and its first `N` derivatives. The order is passed as `Val(N)`, so the result is a tuple of static arrays and the call does not allocate.

```@example der
using Almanac, SPICEKernels

furnsh(naif0012_tls())
furnsh(de440s())
furnsh(pck00011())

et = str2et("2025 MAR 01 12:00:00 UTC")
position, velocity, acceleration = spkder("MOON", et, "J2000", "EARTH", Val(2))
acceleration                                # geometric, km/s²
```

An aberration correction goes before the observer:

```@example der
spkder("MOON", et, "J2000", "LT+S", "EARTH", Val(2))[3]    # apparent, km/s²
```

```@example der
rotation, rate = pxder("J2000", "IAU_EARTH", et, Val(1))
rate[1, 1]                                  # 1/s
```

## Differentiate with ForwardDiff

With ForwardDiff.jl loaded, the epoch can be a `ForwardDiff.Dual` in these calls:

- states: `spkezr`, `spkpos`, `spkgeo`, `spkssb`, the integer-ID forms, `spkltc`, `spkaps`, `spkapo`, `spkacs`, `spkcpo`, `spkcpt`, `spkcvo`, `spkcvt`;
- rotations: `pxform`, `sxform`, `pxfrm2`;
- `spkder`, `pxder` and the table queries;
- surface calls: `subpnt`, `subslr`, `sincpt`, `ilumin`, `illumg`, `illumf`, `phaseq`, `azlcpo`;
- time conversions: `unitim`, `deltet`;
- ellipsoid and conic calls: `nearpt` with a dual point or semi-axes, `prop2b` and `conics` with a dual time.

The other calls (`limbpt`, `termpt`, `tangpt`, `occult`, `fovtrg`, `et2lst`, ...) throw `SpiceError(:badarg)` for a `Dual`. Event searches take windows of `Float64` epochs.

A function built from these calls can be differentiated as a whole. The derivative below includes the change of the light time:

```@example der
using ForwardDiff

apparent(t) = spkpos("MOON", t, "IAU_EARTH", "LT+S", "EARTH")[1]
ForwardDiff.derivative(apparent, et)        # km/s in the rotating frame
```

## Differentiate with other AD packages

The same calls can be differentiated with Enzyme.jl, Mooncake.jl, Zygote.jl (and other users of ChainRulesCore.jl), ReverseDiff.jl and Tracker.jl. Loading the package is enough: an extension gives it the analytical derivatives, and every package returns the derivative that ForwardDiff returns. The bisections of `nearpt` (and so of `subpnt` and `subslr` with `NEAR POINT`) and of `prop2b` and `conics` run on plain numbers, and every package takes the derivative of their root from the equation that the root solves. [DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl) selects a package by its backend type:

```julia
import DifferentiationInterface as DI
using ADTypes: AutoEnzyme, AutoZygote
import Enzyme, Zygote

range_km(t) = sum(abs2, spkpos("MOON", t, "J2000", "LT+S", "EARTH")[1]) |> sqrt
DI.derivative(range_km, AutoEnzyme(), et)   # km/s
DI.gradient(ts -> sum(range_km, ts), AutoZygote(), [et, et + 3600])
```

Reverse mode pays off when one number depends on many epochs, as in the gradient over two epochs above. For a second derivative, run ForwardDiff over one of Zygote, ReverseDiff, or Enzyme in forward mode: `DI.SecondOrder(AutoForwardDiff(), AutoZygote())`. Mooncake and Tracker do not take dual numbers, and Enzyme in reverse mode refuses them. ReverseDiff must run without a compiled tape (`AutoReverseDiff(compile=false)`, the default), because a compiled tape keeps the segment and the record of its first epoch.

## Differences from `spkezr` and `sxform`

Order 0 of [`spkder`](@ref) and [`pxder`](@ref) is the result of `spkpos` and `pxform`. Each higher order is the derivative of that result with respect to the epoch, with the SPK velocity in place of the derivative of the SPK position. For the SPK types whose velocity is the derivative of the position, such as type 2 of the planetary ephemerides, this is the exact derivative. For the types that store the velocity separately, such as type 9, order 1 of a geometric position is the stored velocity, as in `spkezr`, and so is the derivative that ForwardDiff returns for a dual-number epoch.

### Derivatives of positions

With `s = -1` for reception and `s = +1` for transmission, `spkder` differentiates, in the barycentric inertial frame,

```
p(t) = r_targ(t + s ℓ(t)) - r_obs(t)
P(t) = R(t + s |pc(t)| / c) * A(p(t), s * v_obs(t))
```

where `ℓ` is the light time, `A` the stellar aberration and `R` the rotation from J2000 to the output frame.

- With the `LT` corrections, `ℓ(t) = |r_targ(t) - r_obs(t)| / c` is the geometric light time, because the single iteration of CSPICE evaluates the target there.
- With the `CN` corrections, `ℓ` solves `ℓ(t) = |p(t)| / c`, and its derivatives satisfy that equation order by order.
- `A` is the identity without `+S`. Otherwise `A(p, v) = p + |p| sinφ v̂⊥ + |p| (cosφ - 1) p̂`, where `v⊥` is the part of `v` perpendicular to `p` and `sinφ = |v⊥| / c`.
- `R` is constant for an inertial frame. Otherwise its epoch is corrected for the light time to the center of the frame (the Earth for `IAU_EARTH`), and `pc` is `p` with that center as the target.

The velocity of `spkezr` differs from order 1 for two reasons:

- Even for the `LT` corrections, which do not iterate to convergence, CSPICE takes the rate of the light time from the converged equation, which gives up to 8e-9 relative in inertial frames and 7e-8 in body-fixed frames (measured with `de440s.bsp`).
- For `+S`, CSPICE takes the acceleration of the observer from a difference of velocities ±1 s apart, which gives less than 1e-13 relative.

### Derivatives of rotations

| Frame | `pxder` differentiates | Derivative block of `sxform` |
|---|---|---|
| Inertial, PCK, fixed offset | the rotation | the same |
| CK | the interpolation of the CK segment, composed with the clock conversion of `sce2c`; the stored angular velocity is not used | built from the stored angular velocity; only segments that have one are used |
| Dynamic, mean equator or mean ecliptic and equinox of date | the precession and obliquity series | the same |
| Dynamic, true equator and equinox of date | the precession, nutation and obliquity series, with all terms | the nutation rates of CSPICE, which omit the drift of the coefficients |
| Dynamic, two-vector | each defining vector with its aberration correction, as in `spkder` | velocities of `spkezr`; a central difference for a velocity vector; zero for a constant vector |
| Dynamic, Euler | the polynomials of the angles | the same |
| Dynamic, frozen at its `FREEZE_EPOCH` | a constant: zero relative to the base frame | zero |
| Dynamic, with `ROTATION_STATE = 'INERTIAL'` | the rotation of the frame, as if it were rotating | zero |
| Switch, which selects a base frame by epoch | the base that `pxform` selects at the epoch | the block of that base |

For a body whose text PCK has cubic nutation and precession angles (`BODY<sys>_MAX_PHASE_DEGREE = 3`), the rate block of `tisbod_` misses the cubic term, while [`tisbod`](@ref), [`sxform`](@ref) and states in such frames give the exact derivative. With the Moon constants of the NAIF test suite, CSPICE differs from a central difference of `tipbod` by 1.1e-9 and Almanac.jl by 1.6e-13; degrees 1 and 2 agree to rounding.

### Smoothness at record boundaries

For polynomial data (SPK types 2, 3, 8, 9, 12, 13, 14, 18, 20 and binary PCK), every order is the exact derivative of a polynomial of the segment: of the position polynomial, or, for the types that store the velocity separately, of the velocity polynomial from order 1 on. For other types and for frames defined by formulas, the algorithm is evaluated on truncated Taylor series, which gives its derivatives without finite differences.

Either way, a derivative is only as smooth as the data. At the boundary between two SPK records (the time spans of two sets of coefficients or states), the geometric position is continuous to rounding, and higher derivatives can jump. Aberration-corrected positions also have small steps there; see [Known limits of accuracy](@ref).
