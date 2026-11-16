# Differences from CSPICE

Almanac.jl keeps the names, the argument order (with the exceptions below), units and numerical results of CSPICE N0067. It departs from CSPICE where the C interface has hidden state, output arguments or arguments for C memory, and in the argument order of a few calls.

Each docstring links to the NAIF page of the routine it ports: the `_c` page of CSPICE, or the SPICELIB page for a Fortran routine that CSPICE does not wrap. The NAIF page is the full specification; the docstring states the Julia signature, units, and what differs.

## Call conventions

| CSPICE | Almanac.jl |
|---|---|
| Names end in `_c` | No suffix: `spkezr_c` is `spkezr` |
| A call fills or changes an argument | A call that changes an argument ends in `!`: `wninsd!` |
| Results come back through output arguments | Results are returned, several as a tuple: `state, lt = spkezr(...)` |
| Vectors and matrices are C arrays | Results are `SVector` and `SMatrix` of StaticArrays.jl; inputs can be any `AbstractVector` or `AbstractMatrix` |
| An error sets a flag that the program must test | An error throws a [`SpiceError`](@ref): `kind` is the CSPICE short message in lower case (`SPICE(INVALIDDEGREE)` is `:invaliddegree`), or a kind of Almanac.jl such as `:badarg` for a condition CSPICE lacks; `msg` is the text; missing SPK data throw a subtype holding the bodies and the epoch |
| A `found` flag reports a failed lookup | A lookup returns `nothing` ([`bodn2c`](@ref) for an unknown name, [`ckgp`](@ref) without pointing); a call that must return a value throws ([`gdpool`](@ref) for a missing variable) |
| Indices count from 0 | Indices count from 1: window intervals, DSK plates, EK rows |
| A file handle is an integer, and some calls act on a hidden "current" file | A handle is an object ([`DAFHandle`](@ref), [`DASHandle`](@ref)) passed to each call |
| A file open for write cannot be opened for read (`RWCONFLICT`), and a loaded kernel cannot be opened for write (`FILEOPENCONFLICT`) | A file open for write can be read by [`dafopr`](@ref), [`dasopr`](@ref), [`ekopr`](@ref) or [`furnsh`](@ref): a DAF array is complete before its summary, and a DAS file changes on disk only at [`dascls`](@ref). A loaded kernel can be opened for write (e.g. [`spkopa`](@ref)); the context sees new segments only after loading the file again. A second write open still throws `SpiceError(:fileopenconflict)` |
| Loaded kernels are one table per process | Loaded kernels belong to a [`SpiceContext`](@ref); see [Contexts and threads](@ref) |
| Output strings need a length argument | A call returns a `String` |
| Some arguments serve only C memory: counts, sizes, `room`, workspaces | They are absent: [`gdpool`](@ref) takes only the name, [`gfdist`](@ref) has no `nintvls` |
| `spkobj_c`, `spkcov_c`, `ckobj_c`, `ckcov_c`, `pckfrm_c`, `pckcov_c`, `dskobj_c`, `dsksrf_c` take the output cell last | [`spkobj!`](@ref), [`spkcov!`](@ref), [`ckobj!`](@ref), [`ckcov!`](@ref), [`pckfrm!`](@ref), [`pckcov!`](@ref), [`dskobj!`](@ref), [`dsksrf!`](@ref) take it first, as SPICE.jl does for all but `ckobj!` and `ckcov!` |
| A cell or window has a declared size | A [`Window`](@ref) grows as needed |

### Porting from SPICE.jl

SPICE.jl already returns results, uses `!`, returns `nothing` for failed lookups and counts from 1. What changes:

- results are static arrays, which cannot be mutated in place;
- the error type is `Almanac.SpiceError` with its `kind` field, not `SPICE.SpiceError`;
- windows are a [`Window`](@ref), not a cell;
- the optional context argument;
- the C-memory arguments that SPICE.jl keeps are gone (`start` and `room` of `gdpool`, `nintvls` of `gfdist`, `maxn` of `bodvrd`);
- a call that must return a value throws where SPICE.jl returns `nothing` (`gdpool` for a missing variable).

## No state between calls

CSPICE keeps buffers and saved search results between calls for speed. Almanac.jl keeps none: a query is a function of the loaded kernels and its arguments. Consequences:

- **Threads.** Queries need no lock, and kernels can be loaded while other threads query; see [Contexts and threads](@ref).
- **Settings.** [`timdef`](@ref), [`tsetyr`](@ref), [`tparch`](@ref), [`gfstol`](@ref) and [`dskstl`](@ref) set values of a context, not of the process. `gfstol` is the default of each search's `tol` keyword.
- **Repeatable answers.** Where a CSPICE result depends on an earlier call (`ckgp` right after `ckgpav` on overlapping segments, dates in a text kernel loaded after certain `str2et` calls, the four padding bytes of a binary PCK segment summary, which CSPICE fills from the last SPK or CK summary it handled), Almanac.jl returns what CSPICE returns in a new process.

## Numerical results

Each routine follows the arithmetic of its CSPICE counterpart operation by operation, including the order of summation, without fused multiply-add. For the same kernels and `Float64` inputs, results equal CSPICE bit for bit or differ by rounding; [Coverage and validation](@ref) gives the checks and the measured differences.

The elementary functions (`sin`, `cos`, `exp`, ...) are Julia's, while CSPICE uses the platform C library. Julia and glibc differ in the last bit for about 3 % of `sin` and `cos` arguments (10⁶ points between 1 and 51), so a result can differ from CSPICE in the last bit: a type 5 SPK segment propagated over six years puts a planet barycenter 3e-8 km from CSPICE's, at 2e8 km from the Sun.

### Exact derivatives and the rates of CSPICE

Order 0 of [`spkder`](@ref) and [`pxder`](@ref) is the result of `spkpos` and `pxform`; each higher order is the derivative of that result with respect to the epoch, with the SPK velocity in place of the derivative of the SPK position. For SPK types whose velocity is the derivative of the position (type 2 of the planetary ephemerides), this is the exact derivative. For types that store velocity separately (type 9), order 1 of a geometric position is the stored velocity, as in `spkezr`, and so is the partial of a dual-number epoch.

For a position, with `s = -1` for reception and `+1` for transmission, `spkder` differentiates, in the barycentric inertial frame,

```
p(t) = r_targ(t + s ℓ(t)) - r_obs(t)
P(t) = R(t + s |pc(t)| / c) * A(p(t), s * v_obs(t))
```

where `ℓ` is the light time, `A` the stellar aberration and `R` the rotation from J2000 to the output frame.

- `LT` corrections: `ℓ(t) = |r_targ(t) - r_obs(t)| / c`, the geometric light time, since the single CSPICE iteration evaluates the target there.
- `CN` corrections: `ℓ` solves `ℓ(t) = |p(t)| / c`, and its derivatives satisfy that equation order by order.
- `A` is the identity without `+S`; otherwise `A(p, v) = p + |p| sinφ v̂⊥ + |p| (cosφ - 1) p̂`, with `v⊥` the part of `v` perpendicular to `p` and `sinφ = |v⊥| / c`.
- `R` is constant for an inertial frame; otherwise its epoch is corrected for the light time to the frame center (the Earth for `IAU_EARTH`), and `pc` is `p` with that center as target.

The velocity of `spkezr` differs from order 1 for two reasons. CSPICE takes the light-time rate from the converged equation also for `LT` corrections: about 1e-8 relative for the Moon seen from the Earth (8.4e-9 measured with `de440s.bsp`). CSPICE also takes the observer acceleration for `+S` from a ±1 s difference of velocities: about 1e-13 relative.

For a body whose text PCK has cubic nutation and precession angles (`BODY<sys>_MAX_PHASE_DEGREE = 3`), `tisbod_` of CSPICE N0067 divides the cubic term of an angle rate by the square of the seconds per century instead of by the seconds per century, so its rate block misses that term. [`tisbod`](@ref), [`sxform`](@ref) and states in such frames give the exact derivative: with the Moon constants of NAIF's test suite, CSPICE differs from a central difference of `tipbod` by 1.1e-9 and Almanac.jl by 1.6e-13. Degrees 1 and 2 agree to rounding.

For a rotation:

| Frame | `pxder` differentiates | Derivative block of `sxform` |
|---|---|---|
| Inertial, PCK, fixed offset | the rotation | the same |
| CK | the interpolation of the CK segment, composed with the clock conversion of `sce2c`; the stored angular velocity is not used | built from the stored angular velocity; only segments that have one are used |
| Dynamic, mean equator or mean ecliptic and equinox of date | the precession and obliquity series | the same |
| Dynamic, true equator and equinox of date | the precession, nutation and obliquity series, with all terms | the nutation rates of CSPICE, which omit the drift of the coefficients |
| Dynamic, two-vector | each defining vector with its aberration correction, as in `spkder` | velocities of `spkezr`; a central difference for a velocity vector; zero for a constant vector |
| Dynamic, Euler | the polynomials of the angles | the same |
| Dynamic, frozen at its `FREEZE_EPOCH` | a constant: zero relative to the base frame | zero |
| Dynamic, with `ROTATION_STATE = 'INERTIAL'` | the rotating frame | zero |
| Switch, which selects a base frame by epoch | the base that `pxform` selects at the epoch | the block of that base |

For polynomial data (SPK types 2, 3, 8, 9, 12, 13, 14, 18, 20 and binary PCK), any order is the exact derivative of the polynomial. For other types and for frames defined by formulas, the algorithm is evaluated on truncated Taylor series, which gives its derivatives without finite differences. Either way a derivative is only as smooth as the data: at the boundary between two SPK records (the time span of one set of coefficients or states), the geometric position is continuous to rounding and higher derivatives can jump. Aberration-corrected positions also have small steps there; see [Known limits of accuracy](@ref).

## What Almanac.jl adds

| Addition | Calls |
|---|---|
| Derivatives of any order | [`spkder`](@ref), [`pxder`](@ref) |
| Dual-number epochs | the calls that [Dual-number epochs](@ref) lists |
| Tables for inner loops and GPUs | [`ephemeris_table`](@ref), [`rotation_table`](@ref) |
| Aberration correction parsed once | [`Aberration`](@ref): `ab = Aberration("LT+S")`, then `spkezr("MOON", et, "J2000", ab, "EARTH")` |
| Scratch files in memory | [`dasops`](@ref), [`ekops`](@ref) |
| States and rotations resolved once | [`state_query`](@ref), [`rotation_query`](@ref) |
| Bodies and frames from Julia functions | [`add_ephemeris!`](@ref), [`add_frame!`](@ref) |
| Text kernels in JSON | [`furnsh`](@ref) and [`ldpool`](@ref) of a `*.json` file; see [Text kernels in JSON](@ref) |

### Text kernels in JSON

A file named `*.json` is read as a text kernel written as one JSON object, in the format of [jsonspice](https://github.com/jacobwilliams/jsonspice), which [PEPPER](https://github.com/jacobwilliams/pepper) also reads. [`furnsh`](@ref), [`ldpool`](@ref), [`unload`](@ref), [`kinfo`](@ref) and meta-kernels accept it; CSPICE does not. The extension decides: a JSON file with another name is read as a classic text kernel and loads nothing.

```json
// gm.json
{
    "BODY399_RADII": [6378.1366, 6378.1366, 6356.7519],
    "+NAIF_BODY_NAME": "SPOCK",
    "+NAIF_BODY_CODE": -1000,
    "EPOCH": "@2026-NOV-11/12:00",
}
```

Each member is one assignment; a key `"+NAME"` appends, like `+=`. Numbers, `true` (1) and `false` (0) are numeric values; strings are string values, and a string starting with `@` is a time. An array holds numbers and times, or strings. `//` line comments and trailing commas are allowed; other JSON5 syntax is an error.

JSON numbers are read correctly rounded, while text-kernel numbers are read as CSPICE reads them, which can be a few ulps off: the same value in the two formats can give pools that differ in the last bits.

## What is absent

Of the 651 routines declared in the CSPICE N0067 header, 162 have no counterpart of the same name: facilities Julia provides, routines NAIF has deprecated, or C-interface details the conventions above make unnecessary.

| CSPICE routines | In Almanac.jl |
|---|---|
| Error handling: `chkin`, `chkout`, `erract`, `errch`, `errdp`, `errint`, `errprt`, `failed`, `getmsg`, `qcktrc`, `reset`, `return`, `setmsg`, `sigerr`, `trcdep`, `trcnam`, `trcoff` and others | A call throws [`SpiceError`](@ref); catch it with `try` |
| Cells and sets: `insrtd`, `removd`, `elemi`, `union`, `inter`, `diff`, `sdiff`, `scard`, `size`, `valid` and others | `Vector` and `Set`; windows are a [`Window`](@ref) |
| Arrays: `bsrchd`, `isrchi`, `lstled`, `orderd`, `reordd`, `shelld`, `maxd`, `sumad`, `filld`, `moved`, `mequ` and others | `searchsortedfirst`, `findfirst`, `sortperm`, `sort!`, `maximum`, `sum`, `fill!`, `copy` |
| Strings and parsing: `ucase`, `lcase`, `cmprss`, `lparse`, `nextwd`, `prsdp`, `prsint`, `repmc`, `matchw`, `lx4num` and others | `uppercase`, `lowercase`, `split`, `parse`, `replace`, `occursin`, `Regex` |
| Constants: `pi`, `dpmax`, `dpmin`, `intmax`, `intmin` | `π`, `floatmax`, `-floatmax`, `typemax(Int32)`, `typemin(Int32)` |
| Files and the command line: `exists`, `rdtext`, `prompt`, `getcml`, `putcml`, `ftncls` | `isfile`, `eachline`, `readline`, `ARGS`, `close` |
| Deprecated by NAIF: `bodvar`, `illum`, `spkapp`, `srfxpt`, `subpt`, `subsol`, `dafrda` | [`bodvrd`](@ref) or [`bodvcd`](@ref), [`ilumin`](@ref), [`spkaps`](@ref), [`sincpt`](@ref), [`subpnt`](@ref), [`subslr`](@ref), [`dafgda`](@ref) |
| Loading one file by handle: `spklef`, `spkuef`, `cklpf`, `ckupf`, `pcklof`, `pckuof` | [`furnsh`](@ref) and [`unload`](@ref); [`eklef`](@ref) and [`ekuef`](@ref) remain, for EK files only |
| Geometry finder settings: `gfsstp`, `gfstep`, `gfrepi`, `gfrepu`, `gfrepf`, `gfbail`, `gfinth`, `gfclrh`, `gfseth` | The step argument of each search (`udstep` in [`gfevnt`](@ref), [`gffove`](@ref), [`gfocce`](@ref)) and the keywords `report` and `bail`; see [`gfdist`](@ref) |
| The current DAF or DAS file: `dafcs`, `dafgh`, `daswbr` | The handle passed to each call |

The private `zz` routines are left out as well.

The toolkit programs (`brief`, `chronos`, `mkspk`, `spkmerge`, ...) are not included; [SPICEApplications.jl](https://github.com/jdc-pub/SPICEApplications.jl) runs the NAIF executables from Julia.
