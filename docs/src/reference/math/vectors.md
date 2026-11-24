# Vectors and matrices

## Rules for all mathematics calls

These rules hold for every call on the Mathematics pages:

- No kernel is needed, except by [`recpgr`](@ref), [`pgrrec`](@ref), [`drdpgr`](@ref),
  [`dpgrdr`](@ref) and [`xfmsta`](@ref), which read the radii of a body.
- Integer components are converted to `Float64`.
- A sum of products is evaluated term by term in index order, as in CSPICE, so a
  `Float64` input gives the result of CSPICE. The `*` of StaticArrays can round
  differently.
- An input of the wrong size throws `SpiceError(:badarg)`.

## Vectors

```@docs
vadd
vsub
vscl
vminus
vdot
vcrss
vnorm
vhat
vdist
vsep
vrel
vproj
vperp
vlcom
vlcom3
vrotv
vzero
vequ
vpack
vupack
ucrss
unorm
vtmv
frame
```

## Matrices

```@docs
mxv
mtxv
mxm
mtxm
mxmt
xpose
xpose6
invstm
invert
invort
trace
Almanac.det
isrot
ident
```

## Vectors and matrices of any dimension

The `g` calls take no dimension argument: the arguments give the dimension. A call
whose array arguments are all static gives a static result, and does not allocate; a
call with any other array argument gives a `Vector` or a `Matrix` that counts from 1,
with the same elements. Two vectors must have equal length.

```@docs
vaddg
vsubg
vsclg
vminug
vequg
vdotg
vnormg
vhatg
unormg
vdistg
vsepg
vrelg
vprojg
vlcomg
vzerog
vtmvg
mxvg
mtxvg
mxmg
mtxmg
mxmtg
xposeg
```

## Derivatives of vector functions

```@docs
dvhat
dvdot
dvcrss
ducrss
dvnorm
dvsep
```
