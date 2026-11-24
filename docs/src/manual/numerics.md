# Numerical accuracy

This page says how closely the results of Almanac.jl follow those of CSPICE, where they
differ, and how the derivatives of states and rotations are defined.

## Agreement with CSPICE

Each routine follows the arithmetic of its CSPICE counterpart operation by operation,
including the order of summation, and uses no fused multiply-add. With the same kernels
and `Float64` inputs, a result equals that of CSPICE bit for bit or differs by
rounding. [Tests and comparisons](@ref) lists the checks, and [Known limits of accuracy](@ref)
the measured differences.

The elementary functions (`sin`, `cos`, `exp`, ...) are those of Julia, while CSPICE
calls the C library of the platform. Julia and glibc differ in the last bit for about
3 % of the arguments of `sin` and `cos` (measured on 10⁶ points between 1 and 51), so
a result can differ from CSPICE in its last bit. Over long spans this grows: a type 5
SPK segment propagated over six years puts a planet barycenter 3e-8 km away from the
position of CSPICE, at 2e8 km from the Sun. Geometry near a tangency amplifies it: in
NAIF's test suite, a limb point on a DSK shape differs from that of CSPICE by up to
2.5e-3 km, and the intercept of a ray that grazes the surface by up to 4.4e-4 km. A
change of the epoch by one unit in the last place moves the limb point of CSPICE itself
by 3e-3 km. With the `sin` and `cos` of glibc in place of Julia's, these results equal
those of CSPICE.

## Derivatives and stored rates

Order 0 of [`spkder`](@ref) and [`pxder`](@ref) is the result of `spkpos` and `pxform`.
Each higher order is the derivative of that result with respect to the epoch, with the
SPK velocity in place of the derivative of the SPK position. For the SPK types whose
velocity is the derivative of the position, such as type 2 of the planetary
ephemerides, this is the exact derivative. For the types that store the velocity
separately, such as type 9, order 1 of a geometric position is the stored velocity, as
in `spkezr`, and so is the derivative that ForwardDiff returns for a dual-number epoch.

### Derivatives of positions

With `s = -1` for reception and `s = +1` for transmission, `spkder` differentiates, in
the barycentric inertial frame,

```
p(t) = r_targ(t + s ℓ(t)) - r_obs(t)
P(t) = R(t + s |pc(t)| / c) * A(p(t), s * v_obs(t))
```

where `ℓ` is the light time, `A` the stellar aberration and `R` the rotation from
J2000 to the output frame.

- With the `LT` corrections, `ℓ(t) = |r_targ(t) - r_obs(t)| / c` is the geometric
  light time, because the single iteration of CSPICE evaluates the target there.
- With the `CN` corrections, `ℓ` solves `ℓ(t) = |p(t)| / c`, and its derivatives
  satisfy that equation order by order.
- `A` is the identity without `+S`. Otherwise
  `A(p, v) = p + |p| sinφ v̂⊥ + |p| (cosφ - 1) p̂`, where `v⊥` is the part of `v`
  perpendicular to `p` and `sinφ = |v⊥| / c`.
- `R` is constant for an inertial frame. Otherwise its epoch is corrected for the light
  time to the center of the frame (the Earth for `IAU_EARTH`), and `pc` is `p` with
  that center as the target.

The velocity of `spkezr` differs from order 1 for two reasons:

- Even for the `LT` corrections, which do not iterate to convergence, CSPICE takes the
  rate of the light time from the converged equation, which gives about 1e-8 relative
  for the Moon seen from the Earth (8.4e-9 measured with `de440s.bsp`).
- For `+S`, CSPICE takes the acceleration of the observer from a difference of
  velocities ±1 s apart, which gives about 1e-13 relative.

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

For a body whose text PCK has cubic nutation and precession angles
(`BODY<sys>_MAX_PHASE_DEGREE = 3`), the rate block of `tisbod_` misses the cubic term
(see [CSPICE defects](@ref)), while [`tisbod`](@ref), [`sxform`](@ref) and
states in such frames give the exact derivative. With the Moon constants of the NAIF
test suite, CSPICE differs from a central difference of `tipbod` by 1.1e-9 and
Almanac.jl by 1.6e-13; degrees 1 and 2 agree to rounding.

### Smoothness at record boundaries

For polynomial data (SPK types 2, 3, 8, 9, 12, 13, 14, 18, 20 and binary PCK), every
order is the exact derivative of a polynomial of the segment: of the position
polynomial, or, for the types that store the velocity separately, of the velocity
polynomial from order 1 on. For other types and for frames defined
by formulas, the algorithm is evaluated on truncated Taylor series, which gives its
derivatives without finite differences.

Either way, a derivative is only as smooth as the data. At the boundary between two SPK
records (the time spans of two sets of coefficients or states), the geometric position
is continuous to rounding, and higher derivatives can jump. Aberration-corrected
positions also have small steps there; see [Known limits of accuracy](@ref).

## Known limits of accuracy

For the states of the Moon and the planets from `de440s.bsp`, in
inertial and body-fixed frames and with all nine aberration corrections (`NONE` to
`XCN+S`), the largest relative difference of position, velocity and light time measured
by the default test run was 7.3e-16.

Three limits hold in CSPICE as well, because they come from the size of a number or
from the method of CSPICE:

- The rotation angle of the Earth in `IAU_EARTH` grows to about 50,000 rad this
  century. Its rounding alone moves the Sun by some meters in that frame.
- An aberration-corrected position has millimeter steps at the boundaries between SPK
  [records](@ref "Smoothness at record boundaries"), because target and
  observer are evaluated at different epochs and cross a boundary at different times.
- `spkezr` and `spkpos` differ by up to 2e-11 relative with the converged corrections
  (`CN`, `CN+S`), because they stop iterating by different rules.
