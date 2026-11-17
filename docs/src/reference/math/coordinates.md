# Coordinates

The [rules for all mathematics calls](@ref "Rules for all mathematics calls") hold on this page.

## Coordinates

```@docs
reclat
latrec
recrad
radrec
recsph
sphrec
reccyl
cylrec
cyllat
latcyl
cylsph
sphcyl
sphlat
latsph
recgeo
georec
recpgr
pgrrec
xfmsta
convrt
```

## Jacobians of coordinate conversions

Element `[i, j]` of each Jacobian is the derivative of output coordinate `i` with respect
to input coordinate `j`. Angles are in radians.

```@docs
drdlat
dlatdr
drdsph
dsphdr
drdcyl
dcyldr
drdgeo
dgeodr
drdpgr
dpgrdr
```
