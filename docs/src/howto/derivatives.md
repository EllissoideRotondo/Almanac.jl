# Time derivatives

CSPICE gives a position and its velocity. Almanac.jl gives derivatives of any order, either from a direct call or by automatic differentiation of the ordinary calls; the two agree to rounding.

```@setup der
using Almanac
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Direct calls

[`spkder`](@ref) returns a position and its first `N` derivatives, [`pxder`](@ref) a rotation and its first `N` derivatives. The order is passed as `Val(N)`, so the result is a tuple of static arrays and, with upper-case names, the call does not allocate.

```@example der
using Almanac

et = str2et("2025 MAR 01 12:00:00 UTC")
position, velocity, acceleration = spkder("PROBE", et, "J2000", "EARTH", Val(2))
acceleration                                # geometric, km/s²
```

An aberration correction goes before the observer:

```@example der
spkder("PROBE", et, "J2000", "LT+S", "EARTH", Val(2))[3]    # apparent, km/s²
```

```@example der
rotation, rate = pxder("J2000", "IAU_EARTH", et, Val(1))
rate[1, 1]                                  # 1/s
```

## Dual-number epochs

With ForwardDiff.jl loaded, the epoch can be a `ForwardDiff.Dual` in:

- states: `spkezr`, `spkpos`, `spkgeo`, `spkssb`, the integer-ID forms, `spkltc`, `spkaps`, `spkapo`, `spkacs`, `spkcpo`, `spkcpt`, `spkcvo`, `spkcvt`;
- rotations: `pxform`, `sxform`, `pxfrm2`;
- `spkder`, `pxder` and the table queries;
- surface calls: `subpnt`, `subslr`, `sincpt`, `ilumin`, `illumg`, `illumf`, `phaseq`, `azlcpo`;
- time conversions: `unitim`, `deltet`;
- [orbits](@ref "Orbits and relative motion"): `orbit`, and an `Orbit` with a dual epoch or state through `element`, `keplerian_orbit`, `two_body`, the local frames and `relative_state`.

Other calls (`limbpt`, `termpt`, `tangpt`, `occult`, `fovtrg`, `et2lst`, ...) throw `SpiceError(:badarg)` for a `Dual`. Event searches take `Float64` windows.

A function built from these calls differentiates as a whole. The derivative below includes the change of the light time:

```@example der
using ForwardDiff

apparent(t) = spkpos("PROBE", t, "IAU_EARTH", "LT+S", "EARTH")[1]
ForwardDiff.derivative(apparent, et)        # km/s in the rotating frame
```

## Other AD packages

The same calls differentiate with Enzyme.jl, Mooncake.jl, Zygote.jl (and other users of ChainRulesCore.jl), ReverseDiff.jl and Tracker.jl. Loading the package is enough; an extension gives it the analytical derivatives, so no package differentiates the evaluators, and every package returns the derivative that ForwardDiff returns. [DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl) selects a package by its backend type:

```julia
import DifferentiationInterface as DI
using ADTypes: AutoEnzyme, AutoZygote
import Enzyme, Zygote

range_km(t) = sum(abs2, spkpos("MOON", t, "J2000", "LT+S", "EARTH")[1]) |> sqrt
DI.derivative(range_km, AutoEnzyme(), et)   # km/s
DI.gradient(ts -> sum(range_km, ts), AutoZygote(), [et, et + 3600])
```

The orbit calls differentiate with the same packages, with respect to the epoch through [`orbit`](@ref) and with respect to a state or to Keplerian elements:

```julia
x = [8000.0, 0.1, 0.5, 1.0, 2.0, 2.0]       # sma, ecc, inc, raan, aop, ma
radius(x) = sum(abs2, keplerian_orbit(0.0, 399, 1; sma=x[1], ecc=x[2], inc=x[3], raan=x[4],
    aop=x[5], ma=x[6], gm=398600.4415).position) |> sqrt
DI.gradient(radius, AutoZygote(), x)        # km per element unit
```

Reverse mode pays off when one number depends on many epochs, as the gradient above. For a second derivative, run ForwardDiff over Zygote, ReverseDiff or Enzyme in forward mode: `DI.SecondOrder(AutoForwardDiff(), AutoZygote())`. Mooncake and Tracker do not take dual numbers, and Enzyme in reverse mode refuses them. ReverseDiff must run without a compiled tape (`AutoReverseDiff(compile=false)`, the default), because a compiled tape keeps the segment and the record of its first epoch.

## Order 1 and the velocity of `spkezr`

For a geometric position, order 1 is the velocity of `spkezr`. With the `LT` corrections CSPICE approximates the light-time rate, so the two differ by about 1e-8 relative for the Moon seen from the Earth; for some frames `pxder` likewise differs from the rate block of `sxform`. [Exact derivatives and the rates of CSPICE](@ref) gives each case.
