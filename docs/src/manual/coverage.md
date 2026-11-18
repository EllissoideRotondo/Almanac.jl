# Coverage of CSPICE

Almanac.jl covers the library of CSPICE N0067. The API reference lists every call
(starting at [Kernels, contexts and the kernel pool](@ref)); this page says what is in
scope and what replaces the routines that are not.

## What is covered

- **Kernels:** SPK, CK, binary PCK, DSK and EK files, text kernels and meta-kernels.
- **SPK:** every data type of N0067 (1 to 21), every aberration correction, and
  two-line elements (`getelm`, `evsgp4`).
- **Frames:** inertial, body-fixed, fixed-offset, CK-based, dynamic and switch frames.
- **Orientation and clocks:** CK types 1 to 6 and spacecraft clocks of type 1.
- **Time:** parsing, formatting and conversion between time systems, local solar time.
- **Kernel pool:** reads, changes and watchers (one set of watchers per context).
- **Geometry:** sub-points, intercepts, illumination, occultation, limb and terminator,
  fields of view, on ellipsoids and on DSK shapes.
- **Geometry finder:** all `gf` searches, including the user-defined ones.
- **Mathematics:** vectors, matrices, rotations, quaternions, coordinate systems,
  ellipsoids, conics and interpolation.
- **Writers:** SPK (all types), CK (5 types), binary PCK (types 2 and 20), DSK type 2
  and EK, plus low-level DAF and DAS access.

The features beyond CSPICE (derivatives, tables, user-defined bodies and frames,
orbits) are in [What Almanac.jl adds](@ref).

## What is left out

Routines that Julia or the design of the package makes unnecessary:

- **Error handling** (`chkin`, `sigerr`, `failed`, `reset`, ...): a call throws a
  [`SpiceError`](@ref); catch it with `try`.
- **Cells, sets and arrays** (`insrtd`, `union`, `bsrchd`, `orderd`, `sumad`, ...):
  `Vector`, `Set` and Base functions; windows are a [`Window`](@ref).
- **Strings and parsing** (`ucase`, `lparse`, `prsdp`, `matchw`, ...): `uppercase`,
  `split`, `parse`, `Regex`.
- **Constants** (`pi`, `dpmax`, `intmax`, ...): `π`, `floatmax`, `typemax(Int32)`.
- **Files and the command line** (`exists`, `rdtext`, `prompt`, `getcml`, ...):
  `isfile`, `eachline`, `readline`, `ARGS`.
- **Routines deprecated by NAIF** (`bodvar`, `illum`, `spkapp`, `srfxpt`, `subpt`,
  `subsol`, `dafrda`): [`bodvrd`](@ref), [`ilumin`](@ref), [`spkaps`](@ref),
  [`sincpt`](@ref), [`subpnt`](@ref), [`subslr`](@ref), [`dafgda`](@ref).
- **Loading one file by handle** (`spklef`, `cklpf`, `pcklof`, and their unload
  forms): [`furnsh`](@ref) and [`unload`](@ref). [`eklef`](@ref) remains for EK files.
- **Settings of the geometry finder** (`gfsstp`, `gfrepi`, `gfbail`, ...): the step
  argument and the keywords `report` and `bail` of each search; see [`gfdist`](@ref).
- **The current DAF or DAS file** (`dafcs`, `dafgh`, `daswbr`): the handle passed to
  each call.
- **Private `zz` routines and the toolkit programs** (`brief`, `mkspk`, ...).
  [SPICEApplications.jl](https://github.com/jdc-pub/SPICEApplications.jl) runs the
  NAIF executables from Julia.
