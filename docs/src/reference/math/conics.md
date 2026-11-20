# Conics, polynomials and constants

The [rules for all mathematics calls](@ref "Rules for all mathematics calls") hold on this page.

## Conic elements and two-body propagation

```@docs
conics
oscelt
oscltx
eqncpv
kpsolv
prop2b
stmp03
```

## Polynomials and interpolation

A Chebyshev expansion, a sum of Chebyshev polynomials of the first kind, is given by its
coefficients and by `x2s = [midpoint, radius]` of its interval, as in an SPK segment of
type 2 or 3. A table for interpolation is given by its abscissas, or by the first
abscissa and the step of an equally spaced table; a Hermite table holds the value and
the derivative at each abscissa, in this order.

```@docs
chbval
chbint
chbder
chbigr
lgrint
lgrind
lgresp
hrmint
hrmesp
polyds
qderiv
rquad
diags2
```

## Doubles as hexadecimal strings

A double is written as hexadecimal digits of a mantissa in [1/16, 1) and an exponent of
16.

```@docs
dp2hx
hx2dp
```

## Mathematical and physical constants

```@docs
clight
halfpi
twopi
rpd
dpr
```
