# DAS and DSK files

DSK and EK files are DAS (direct access segregated) files. A file handle is a
[`DASHandle`](@ref) object, not an integer.

`[ctx,]` and `[ctx | ks,]` in the signatures are the optional context of
[Contexts and threads](@ref); [`dskxv!`](@ref) takes it after its two
output arguments.

## DAS and DLA files

```@docs
DASHandle
dasopr
dasopw
dasonw
dasops
dascls
dasllc
dashfn
dasrfr
daslla
dashfs
dasrdd
dasrdi
dasrdc
dasadd
dasadi
dasadc
dasudd
dasudi
dasudc
dasac
dasec
dasdc
DLADescr
dlaopn
dlabns
dlaens
dlabfs
dlafns
dlabbs
dlafps
```

## DSK files

The spatial index of a type 2 segment divides a box around the plates into cubes of
equal size, the fine voxels, which form groups of `corscl^3`, the coarse voxels;
`corscl`, the coarse voxel scale, is the number of fine voxels along one edge of a
coarse voxel. The voxel-plate list gives the plates that meet each fine voxel, and the
vertex-plate list the plates that share each vertex. The DSK Required Reading of NAIF
([https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/dsk.html](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/dsk.html))
describes the index.

```@docs
DSKDescr
DSKSegment
dskobj
dskobj!
dsksrf
dsksrf!
dskgd
dskb02
dskz02
dskn02
dskp02
dskv02
dski02
dskd02
dskx02
dskxv
dskxv!
dskxsi
dskgtl
dskstl
dskopn
dskw02
dskmi2
dskrb2
dskcls
pltexp
pltnp
pltnrm
pltar
pltvol
```

## Surface names

```@docs
srfc2s
srfs2c
srfcss
srfscc
```
