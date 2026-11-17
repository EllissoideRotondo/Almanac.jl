# Coverage

Almanac.jl covers the library of CSPICE N0067. The routines it leaves out are replaced
by facilities of Julia or by the conventions of this package; they are listed in
[What is absent](@ref).

## What is available

- **Kernels.** SPK, CK, binary PCK, DSK and EK files; text kernels (leap seconds, PCK,
  frames, spacecraft clock); meta-kernels.
- **SPK data types.** All the types of CSPICE N0067: 1, 2, 3, 5, 8, 9, 10, 12, 13, 14,
  15, 17, 18, 19, 20 and 21. Two-line elements through `getelm` and `evsgp4`.
- **States.** `spkezr`, `spkpos`, `spkgeo`, `spkssb`; the integer-ID forms `spkez`,
  `spkezp` and `spkgps`; `spkltc`, `spkaps`, `spkacs`, `spkapo`; every aberration
  correction of CSPICE. Observers and targets of constant position or velocity
  (`spkcpo`, `spkcpt`, `spkcvo`, `spkcvt`).
- **Frames.** Inertial, body-fixed (text and binary PCK), fixed-offset, CK-based,
  dynamic (all families) and switch frames, through `pxform` and `sxform`. Frame
  information: `frinfo`, `cidfrm`, `cnmfrm`, `ccifrm`, `tkfram`, `tipbod`, `tisbod`,
  `pckfrm`, `pckcov`, `bltfrm`, `kplfrm`.
- **Orientation and clocks.** CK types 1 to 6; `ckgp`, `ckgpav`, `ckcov`, `ckobj`,
  `ckfrot`, `ckfxfm`; the records of CK types 2 and 3 (`cknr02`, `ckgr02`, `cknr03`,
  `ckgr03`); spacecraft clocks of type 1.
- **Time.** `str2et`, `utc2et`, `timout`, `et2utc`, `etcal`, `unitim`, `tparse`,
  `tpictr`, `timdef`, `tsetyr`, `tparch`, `deltet`, `et2lst`, `lspcn`.
- **Kernel pool.** Reads (`gdpool`, `gcpool`, `gipool`, `bodvrd`, `bodvcd`, `bodfnd`,
  `expool`, `dtpool`, `gnpool`, `szpool`, `stpool`, `badkpv`), changes (`pdpool`,
  `pcpool`, `pipool`, `dvpool`, `lmpool`, `ldpool`, `clpool`, `boddef`), the body and
  frame name calls, and the watchers `swpool` and `cvpool`, one set per context.
- **Geometry.** `subpnt`, `subslr`, `sincpt`, `ilumin`, `illumg`, `illumf`, `phaseq`,
  `occult`, `limbpt`, `termpt`, `tangpt`, `getfov`, `getfvn`, `fovtrg`, `fovray`,
  `azlcpo`, `azlrec`, `recazl`, `dazldr`, `drdazl`, `srfrec`, on ellipsoids and on DSK
  shapes.
- **Geometry finder.** `gfdist`, `gfsep`, `gfposc`, `gfrr`, `gfpa`, `gfsubc`,
  `gfsntc`, `gfilum`, `gfoclt`, `gftfov`, `gfrfov`, `gfuds`, `gfudb`; `gfevnt` for the
  six quantities of CSPICE, `gfocce`, `gffove`; `gfstol`, `gfrefn`, `uddf`, `uddc`.
- **Mathematics.** Vectors, matrices, rotations and quaternions; coordinate
  conversions and their Jacobians; ellipsoids, planes and ellipses; conic elements and
  two-body propagation; equinoctial elements (`eqncpv`); Chebyshev, Lagrange and
  Hermite interpolation; doubles as base 16 strings (`dp2hx`, `hx2dp`).
- **DAF files and writers.** The DAF read and write calls and the comment area. SPK
  writers for all 16 types (type 14 through `spk14b`, `spk14a`, `spk14e`), CK writers
  for 5 types, binary PCK writers of types 2 and 20. SPK segments by descriptor
  (`spksfs`, `spkpvn`, `spkuds`, and `spksub` for every type). `spkobj`, `spkcov`,
  `kinfo`, `kdata`.
- **Shape kernels.** DAS file access; DSK type 2 (the only DSK type of CSPICE N0067),
  read and write; the plate access calls; ray-surface intercepts (`dskxv`, `dskxsi`,
  `dskx02`); `latsrf`, `srfnrm`; the surface name calls.
- **Event kernels.** Read access (`ekopr`, `ekssum`, `ekrcec`, `ekrced`, `ekrcei`);
  loading one file (`eklef`, `ekuef`); the table catalog (`ekntab`, `ektnam`,
  `ekccnt`, `ekcii`); queries (`ekpsel`, `ekfind`, `eknelt`, `ekgc`, `ekgd`, `ekgi`),
  each with its own result object; writers by record (`ekopn`, `ekops`, `ekbseg`,
  `ekappr`, `ekinsr`, `ekdelr`, `ekace*`, `ekuce*`) and the fast writers, which write
  whole columns (`ekifld`, `ekacl*`, `ekffld`).

The features beyond CSPICE (derivatives, tables, user-defined bodies and frames,
orbits) are in [What Almanac.jl adds](@ref).

## What is absent

Of the 651 routines declared in the CSPICE N0067 header, 162 have no counterpart of the
same name. Each group below names its replacement in Julia or in Almanac.jl.

- **Error handling** (`chkin`, `chkout`, `erract`, `errch`, `errdp`, `errint`,
  `errprt`, `failed`, `getmsg`, `qcktrc`, `reset`, `return`, `setmsg`, `sigerr`,
  `trcdep`, `trcnam`, `trcoff` and others): a call throws a [`SpiceError`](@ref);
  catch it with `try`.
- **Cells and sets** (`insrtd`, `removd`, `elemi`, `union`, `inter`, `diff`, `sdiff`,
  `scard`, `size`, `valid` and others): `Vector` and `Set`; windows are a
  [`Window`](@ref).
- **Arrays** (`bsrchd`, `isrchi`, `lstled`, `orderd`, `reordd`, `shelld`, `maxd`,
  `sumad`, `filld`, `moved`, `mequ` and others): `searchsortedfirst`, `findfirst`,
  `sortperm`, `sort!`, `maximum`, `sum`, `fill!`, `copy`.
- **Strings and parsing** (`ucase`, `lcase`, `cmprss`, `lparse`, `nextwd`, `prsdp`,
  `prsint`, `repmc`, `matchw`, `lx4num` and others): `uppercase`, `lowercase`,
  `split`, `parse`, `replace`, `occursin`, `Regex`.
- **Constants** (`pi`, `dpmax`, `dpmin`, `intmax`, `intmin`): `π`, `floatmax`,
  `-floatmax`, `typemax(Int32)`, `typemin(Int32)`.
- **Files and the command line** (`exists`, `rdtext`, `prompt`, `getcml`, `putcml`,
  `ftncls`): `isfile`, `eachline`, `readline`, `ARGS`, `close`.
- **Routines deprecated by NAIF** (`bodvar`, `illum`, `spkapp`, `srfxpt`, `subpt`,
  `subsol`, `dafrda`): [`bodvrd`](@ref) or [`bodvcd`](@ref), [`ilumin`](@ref),
  [`spkaps`](@ref), [`sincpt`](@ref), [`subpnt`](@ref), [`subslr`](@ref),
  [`dafgda`](@ref).
- **Loading one file by handle** (`spklef`, `spkuef`, `cklpf`, `ckupf`, `pcklof`,
  `pckuof`): [`furnsh`](@ref) and [`unload`](@ref). [`eklef`](@ref) and
  [`ekuef`](@ref) remain, for EK files only.
- **Settings of the geometry finder** (`gfsstp`, `gfstep`, `gfrepi`, `gfrepu`,
  `gfrepf`, `gfbail`, `gfinth`, `gfclrh`, `gfseth`): the step argument of each search
  (`udstep` in [`gfevnt`](@ref), [`gffove`](@ref) and [`gfocce`](@ref)), and the
  keywords `report` and `bail`; see [`gfdist`](@ref).
- **The current DAF or DAS file** (`dafcs`, `dafgh`, `daswbr`): the handle passed to
  each call.

The private `zz` routines are left out as well. The toolkit programs (`brief`,
`chronos`, `mkspk`, `spkmerge`, ...) are not included;
[SPICEApplications.jl](https://github.com/jdc-pub/SPICEApplications.jl) runs the NAIF
executables from Julia.
