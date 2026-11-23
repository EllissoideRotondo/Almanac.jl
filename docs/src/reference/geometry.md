# Geometry

A reception correction (`"LT"`, `"CN"`, with or without `"+S"`) gives the target as the
observer sees it at `et`, and a transmission correction (`"XLT"`, `"XCN"`) gives it
where a signal that the observer sends at `et` reaches it; [`Aberration`](@ref) lists
the corrections.

## Points on a body

```@docs
subpnt
subslr
sincpt
srfrec
latsrf
srfnrm
```

## Illumination

```@docs
ilumin
illumg
illumf
phaseq
trgsep
```

## Occultation, limb and terminator

```@docs
occult
limbpt
termpt
edterm
tangpt
```

## Fields of view

```@docs
fovtrg
fovray
getfov
getfvn
```

## Azimuth and elevation

```@docs
azlcpo
azlrec
recazl
dazldr
drdazl
```
