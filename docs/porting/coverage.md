# Supported routines

Almanac.jl covers the library of CSPICE N0067. The API reference lists every call
(starting at [Kernels, contexts and the kernel pool](@ref)); this page says what is in
scope and what replaces the routines that are not.

## What is covered

- **Kernels:** SPK, CK, binary PCK, DSK and EK files, text kernels and meta-kernels.
- **SPK:** all 16 data types of N0067 (1, 2, 3, 5, 8, 9, 10, 12 to 15 and 17 to 21),
  every aberration correction, and
  two-line elements (`getelm`, `evsgp4`).
- **Frames:** inertial, body-fixed, fixed-offset, CK-based, dynamic and switch frames.
- **Orientation and clocks:** CK types 1 to 6, binary PCK types 2, 3 and 20, and
  spacecraft clocks of type 1.
- **Time:** parsing, formatting and conversion between time systems, local solar time.
- **Kernel pool:** reads, changes and watchers (one set of watchers per context).
- **Geometry:** sub-points, intercepts, illumination, occultation, limb and terminator,
  fields of view, on ellipsoids and on DSK shapes.
- **Geometry finder:** all `gf` searches, including the user-defined ones.
- **Mathematics:** vectors, matrices, rotations, quaternions, coordinate systems,
  ellipsoids, conics and interpolation. `det` is public but not exported, because
  `LinearAlgebra` exports a `det` too: call it as `Almanac.det`.
- **Writers:** SPK (all 16 types), CK (types 1 to 6), binary PCK (types 2, 3 and 20), DSK
  type 2 and EK, plus low-level DAF and DAS access.
- **SPICELIB routines without a C wrapper**, which a caller of CSPICE reaches through
  their f2c names (CSPICE is SPICELIB translated to C by f2c): among them the writers of CK type 4 and binary PCK type 3, the record readers of
  CK types 1 to 6, [`kpsolv`](@ref), [`npsgpt`](@ref), [`incnsg`](@ref),
  [`scn2id`](@ref) and [`scid2n`](@ref). Their docstrings say so.

The features beyond CSPICE (derivatives, tables, user-defined bodies and frames) are in [Beyond CSPICE](@ref).

## What is left out

Routines that Julia or the design of the package makes unnecessary:

- **Error handling** (`chkin`, `sigerr`, `failed`, `reset`, ...): a call throws a
  [`SpiceError`](@ref); catch it with `try`.
- **Cells, sets and arrays** (`insrtd`, `union`, `bsrchd`, `orderd`, `sumad`, ...):
  `Vector`, `Set` and Base functions; windows are a [`Window`](@ref).
- **Strings and parsing** (`ucase`, `lparse`, `prsdp`, `matchw`, ...): `uppercase`,
  `split`, `parse`, `Regex`.
- **Constants** (`pi`, `dpmax`, `intmax`, ...): `π`, `floatmax`, `typemax(Int32)`.
- **Copies and placeholders** (`mequ`, `mequg`, `udf`): `copy`, and an ordinary function
  where a user-defined search needs one.
- **Files and the command line** (`exists`, `rdtext`, `prompt`, `getcml`, ...):
  `isfile`, `eachline`, `readline`, `ARGS`.
- **Routines deprecated by NAIF** (`bodvar`, `illum`, `spkapp`, `srfxpt`, `subpt`,
  `subsol`, `dafrda`): [`bodvrd`](@ref), [`ilumin`](@ref), [`spkaps`](@ref),
  [`sincpt`](@ref), [`subpnt`](@ref), [`subslr`](@ref), [`dafgda`](@ref).
- **Loading one file by handle** (`spklef`, `cklpf`, `pcklof`, and their unload
  forms): [`furnsh`](@ref) and [`unload`](@ref). [`eklef`](@ref) remains for EK files.
- **Settings of the geometry finder** (`gfsstp`, `gfrepi`, `gfbail`, ...): the step
  argument and the keywords of each search: `report`, a progress callback, and `bail`,
  a function that stops the search when it returns `true`; see [`gfdist`](@ref).
- **The current DAF or DAS file** (`dafcs`, `dafgh`, `daswbr`): the handle passed to
  each call.
- **Private `zz` routines and the toolkit programs** (`brief`, `mkspk`, ...).
  [SPICEApplications.jl](https://github.com/jdc-pub/SPICEApplications.jl) runs the
  NAIF executables from Julia.
