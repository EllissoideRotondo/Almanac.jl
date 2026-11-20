# Frames and orientation

## Rotations between frames

```@docs
pxform
pxfrm2
sxform
```

## Prepared rotation queries

```@docs
rotation_query
RotationQuery
```

## Frames defined by Julia functions

```@docs
add_frame!
```

## Frame information

```@docs
frinfo
cidfrm
cnmfrm
ccifrm
tkfram
tipbod
tisbod
pckfrm
pckfrm!
pckcov
pckcov!
bltfrm
kplfrm
```

## Orientation from CK files

```@docs
ckgp
ckgpav
ckcov
ckcov!
ckobj
ckobj!
ckmeta
ckfrot
ckfxfm
```

## Records of CK segments

These calls read one segment of an open CK file, given by its descriptor from a DAF
search ([`dafbfs`](@ref), [`daffna`](@ref), [`dafgs`](@ref)).

```@docs
cknr01
ckgr01
cknr02
ckgr02
cknr03
ckgr03
cknr04
ckgr04
cknr05
ckgr05
cknm06
ckmp06
ckgr06
```
