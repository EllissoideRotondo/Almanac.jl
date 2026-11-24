# Kernels, contexts and the kernel pool

## Contexts and errors

```@docs
SpiceContext
default_context
Almanac.kernelset
Almanac.KernelSet
SpiceError
SpiceKindError
SpkInsufficientDataError
Almanac.ERROR_KINDS
```

## Loading and inspecting kernels

```@docs
furnsh
unload
kclear
ktotal
kinfo
kdata
tkvrsn
getfat
spkobj
spkobj!
spkcov
spkcov!
```

## Reading the kernel pool

```@docs
gdpool
gcpool
gipool
bodvrd
bodvcd
bodfnd
expool
dtpool
gnpool
szpool
stpool
badkpv
```

## Changing the kernel pool

```@docs
pdpool
pcpool
pipool
dvpool
lmpool
ldpool
clpool
swpool
cvpool
poolversion
```

## Body and frame names

```@docs
bodn2c
bodc2n
bodc2s
bods2c
scn2id
scid2n
boddef
namfrm
frmnam
```
