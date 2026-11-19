# Agreement with CSPICE: state and plan

Working document for the tspice replay work on branch `tspice-harness-review` of
EllissoideRotondo/Almanac.jl. An agent that picks this up reads it first, updates it
as items close, and commits it with the work.

## Rules for this work

- CSPICE N0067 is the reference. SPICELIB (PEPPER) is what tspice recorded; where
  CSPICE and SPICELIB differ, Almanac follows CSPICE.
- Almanac keeps Julia's `sin` and `cos` (decision of the maintainer). Where a function
  needs both of one argument, it calls `sincos(x)` (bit-identical to `sin` and `cos`).
  A difference that only Julia's libm causes is accepted, but only with evidence: the
  replay with glibc's functions swapped in (below) must remove it.
- Follow CSPICE's arithmetic (formula and order of operations) where that costs no
  accuracy, no generality and no support for dual numbers.
- No tolerance is loosened to make a check pass. A tolerance is a measured rounding
  difference with a stated reason.
- No public issues on the fork or on Orbitology/Almanac.jl. Commits end with the
  `Co-Authored-By` and `Claude-Session` lines of the session.

## How to measure

```sh
julia --project=test test/tspice/run.jl [--update-baseline]      # full replay and gates
julia --project=test test/tspice/selftest.jl                     # the gates reject faults
TSPICE_ORACLE=1 TSPICE_OUT=out julia --project=test test/tspice/replay.jl LOGDIR F_XXX ...
TSPICE_ORACLE_ALL=1 ...    # also: cspice_agreement.tsv, every call against CSPICE
TSPICE_AUDIT=1 ...         # also: audit.tsv, matches that compared few or no outputs
TSPICE_AGREE_ROWS=sxform,pxform ...   # with ORACLE_ALL: agreement_rows.tsv, each call of
                           # these routines not equal to CSPICE, outputs and ulps
```

The glibc check: run `replay.jl` through a script that first redefines, for `Float64`,
`Base.Math.sin`, `cos`, `tan`, `asin`, `acos`, `atan` (one and two arguments), `sinh`,
`cosh`, `tanh`, `exp`, `log`, `log10`, `expm1`, `log1p`, `^` as `ccall`s into
`libm.so.6`, and `sincos(x) = (sin(x), cos(x))`, then `include`s `replay.jl`. A
mismatch that remains with glibc is not a libm effect.

## State at commit 5c208be (third recording, new baseline)

Full recording with the recorder of this branch (399 families, all passed), full replay
with the oracle and `--update-baseline`, then the GF bound fixed (2 rows added):

| kind | count | meaning |
|---|---|---|
| match | 16,019,710 | equal to SPICELIB within 1e-12 (GF windows: within 2 tol) |
| cspice | 163,595 | differs from SPICELIB, equal to CSPICE, CSPICE differs from SPICELIB |
| near | 2,893 | within 1e-9 of SPICELIB; CSPICE agrees with SPICELIB |
| spread | 137 | SPICELIB and CSPICE differ; Almanac within twice their distance |
| value | 823 | larger difference (all classed: item 7) |
| error | 14 | other error, or an error where SPICELIB had none (classed) |
| adapter | 0 | |
| skipped | 186,903 | 177,888 unmapped handles (F_ZZDSKBSR 69,728: files past the copy cap and loaded by its utilities), 24 without counterpart, 31 truncated, 14 in an error state |

Written files compared at close: 38,323, all equal to SPICELIB's; 9,685 not compared
(9,454 past the copy cap, mostly F_ZZDSKBSR; 215 opened outside the log; 16 others).
GF searches compared: 712 (rest of the list in `gf_agreement.tsv`), 2 outside 2 tol
(`gfpa`, `gfilum` of F_GFILUM, libm, documented). Classes of the gated rows: documented
833, cspice 3, gap 1, unexplained 0. `selftest.jl` rejects every case. Unit tests: 11,463
pass (full suite on d9fab11; the groups changed since pass).

Earlier measurements, still valid:
- glibc check on F_DYN01, F_XFMSTA, F_GFDIST, F_LS: near 1,352 -> 0, value 18 -> 4 (the
  4 are `tkvrsn`). The large value rows of F_DYN01 (to 3e-7 of the vector magnitude)
  are libm, amplified by a cancellation in a frame that co-rotates with Phobos. The
  remaining near rows (`tangpt` 712, `pxform` 489, `sxform` 405, `spkezr` 287, `spkcpt`
  242, `conics` 216, `spkpos` 105, `xfmsta` 95, `sincpt` 94) are of this kind; `tisbod`
  129 is the CSPICE defect of cubic rates (Almanac exact).
- Agreement with CSPICE, every call of 34 families (plan-rates, glibc check): equal to
  the last bit 878,267 of 960,248 calls; `spkezr`, `spkpos`, `spkez`, `azlcpo` on
  every call.

## Fixed since the last baseline: the 38,164 `dasrdc` gaps

They were a recorder bug, not a limit of the replay. `F_DASA2L` writes strings that
hold every character code, line feeds included, and the line-based log split them, so
Almanac wrote different data and every read differed. `instrument.py` now logs a
string with a character outside printable ASCII as `%` and hex codes; `replay.jl`
decodes it. F_DASA2L re-recorded and replayed: `dasrdc` 87,500 match, 0 value. Only
F_DASA2L had such strings (scan of all logs).

To do first on a new machine: record again (instrument.py changed: items 6 and 8), then
`run.jl --update-baseline` and `selftest.jl`.

## Open items

Each item: what to do, how to know it is done.

1. **Port the CSPICE rate paths that Almanac computes another way.** Done (branch
   plan-rates). The state mode of `_dyn_ofdate` builds the of-date frames as
   `zzdynfrm_`: `_dyn_prec_state` (`zzeprc76_`), `_dyn_enut80_state` (`zzenut80_`),
   `_dyn_mobliq_state` (`zzmobliq_`) through `eul2xf_`, products `mxmg_`, inverse
   `invstm_`. Near-point vectors take `dnearp` in the state mode (`_dyn_dnearp`). A
   switch frame over a PCK base takes the `tisbod_` rates (`_ck_pck_step`), not the
   `Taylor` derivative. `DT_TWOVEC_IN_DYN` followed from the of-date fix. Found on the
   way and fixed: `zzstelab_` takes scaled norms, projections and unit vectors
   (`vnorm_`, `vperp_`, `dvhat_`), which the old `_zzstelab` did not (one ulp in the
   `LT+S` velocity in ECLIPJ2000, seen in `DT_GSE_LTS`); every light time is `vnorm_`
   of the position over `clight_` (the plain norm differed by one ulp in about a third
   of all `spkezr`/`spkpos` calls; `_ltc`, `_apo`, `spkezr`/`spkpos`/`spkgeo`/`spkgps`
   geometric, `geometry/constobs.jl`, `_dyn_apo`/`_dyn_spkez`; and on branch plan-rest
   the rest: `subpnt` and `subslr` (`vdist_`/`vnorm_` of the sub-point vector),
   `zzsfxcor_` of `sincpt`/`ilumin`/`illum*` (light times, `range`, `relerr`, `srflen`
   and `udir` with `vhat_`), `spkezp_`/`spkez_` of the geometry routines
   (`_geo_spkezp`, `_geo_spkez`) and `ltime`; `limbpt`, `termpt`, `tangpt` and the GF
   state functions had them already; the light-time series of `abcorr_jet` are not
   ports and keep the norm of their number type. Glibc check, `TSPICE_ORACLE_ALL=1`, 33
   families (F_SUBPNT, F_SUBSLR, F_SINCPT, F_ILUMIN, F_ILLUMG, F_ILLUMF, F_LIMBPT,
   F_TERMPT, F_TERM, F_TANGPT0..3, F_TANGPTW, F_OCCULT, F_GFDIST, F_GFSUBC, F_GFSNTC,
   F_GFILUM, F_GFOCLT, F_LTIME, and the F_ZZGF* state functions): calls equal to CSPICE
   999,647 -> 1,001,223 of 1,064,791 (`sincpt` 3,124 -> 4,690 of 4,690, with
   `_stelab` also `spkpos` 31,362 -> 31,367, `tangpt` 33,704 -> 33,707, `stelab` 202 ->
   204: every call of the geometry routines now exact); ids identical); `recazl` takes the
   azimuth of the Fortran `reclat_` (unscaled `atan2`), as `recazl_c` and `azlcpo_` do.
   Measured with the glibc check: all 25 probe frames of `test/dynfixtures.jl` equal
   to `sxform_c` on 300 epochs each (before: 9 frames differed); `azlcpo` equal to
   `azlcpo_c` in every output on 400 random cases over IAU_EARTH, IAU_MOON, IAU_VENUS,
   MOON_PA_DE440, MOON_ME_DE440_ME421 (before: light time about 30 percent, azimuth a
   few, the binary PCK rates nearly all). Not covered by a regression test that is
   independent of libm: `dnearp` in two-vector frames (its frames need `IAU_MOON`).
2. **Binary PCK rates.** Done (branch plan-rates). `_binary_pck_state`: angles and
   rates of `chbint_` (`pck_eval` of order 1, type 2) or `spke20_` (type 20), W reduced
   as `pcke02_`, then `eul2xf_` and the right product with the inertial rotation, as
   `pckmat_` and `tisbod_`. Glibc check: `sxform` from J2000, ECLIPJ2000 and IAU_EARTH
   to MOON_PA_DE440 equal to `sxform_c` on 200 epochs each (before: every derivative
   block differed). Families with binary PCKs: F_PCK20 (type 20), F_PCKBSR,
   F_DDHFTSIZE (files only).

   Replay of 34 families with the glibc check, `TSPICE_ORACLE_ALL=1`, before and after
   (F_DYN01..08, F_XFMSTA, F_ZZSWFXFM, F_CKFXFM, F_AZL, F_AZLCPV, F_AB, F_SPKEZ,
   F_SPKEZP, F_SPKAPO, F_SPKCOR, F_SPKCPV, F_SPKGEO, F_SPKGPS, F_ZZSTELAB, F_LTIME,
   F_STLABX, F_SPKAPP, F_ZZSPKFUN, F_SUBPNT, F_ILLUMG, F_PCK20, F_RDPCK, F_RDPCKD,
   F_PCKBUF, F_FRMCHG, F_SXFORM): the tspice classes are the same row for row (no new
   near or value). Calls equal to CSPICE to the last bit: 745,267 -> 878,267 of
   960,248 (`sxform` 239,272 -> 349,575 of 349,632; `spkez` 20,789 -> 32,265 of 32,265;
   `azlcpo`, `spkezr`, `spkpos`, `spkltc`, `spkacs`, `spkcpo`, `spkcpt`, `spkcvo`,
   `spkcvt` now all). The 57 `sxform`, 137 `pxform`, 1 `spkgeo` and 1 `spkgps` calls
   left within 1e-12 but not exact: done (branch plan-rest), all exact now. Found with
   `TSPICE_AGREE_ROWS` (below); four port differences:
   - `_stelab` (`stelab_` with `vrotv_`) took plain norms and a plain projection where
     CSPICE takes `vhat_`, `vnorm_` and `vproj_`: the dynamic frames with stellar
     aberration (pxform 140 -> 25 not exact, sxform 59 -> 26);
   - frame chains (`_chain_rotation` in `src/ck.jl`, used by `frame_rotation` and
     `_ck_state_rotation`): a CK or switch frame was linked to J2000 in one step, so
     the two chains met at J2000 although the CK segment's frame was the target or on
     its chain (CK_-9999 to GALACTIC through J2000 and back, up to 24 ulp). Now as
     `zzrefch0_`/`zzfrmch0_`: the chain of `from` is evaluated link by link to `to` or
     J2000 (a CK frame without data ends it, as `zzrotgt0_`), then the chain of `to`
     until it meets it;
   - `_ck_axisar` (CK types 2 and 3) projected with a plain dot product, `axisar_` uses
     `vrotv_` (`vproj_`);
   - TK frames: ANGLES were scaled by pi/180, `tkfram_` calls `convrt_` (degrees over
     degrees per radian, also from RADIANS); MATRIX and QUATERNION took plain norms,
     `sharpr_` and `vhatg_` scaled ones (TOPOCENTRIC to TST-PHOENIX 240 ulp, GOLDSTONE_TOPO).
   Glibc check of the 16 families (F_DYN01..08, F_XFMSTA, F_ZZSWFXFM, F_CKFXFM, F_FRMCHG,
   F_SXFORM, F_PXFORM, F_SPKGEO, F_SPKGPS), `TSPICE_ORACLE_ALL=1`: calls equal to CSPICE
   350,921 -> 351,176 of 382,496 (all calls of `sxform` 57,109, `pxform` 36,490,
   `spkgeo`, `spkgps`, `stelab`, `ckgp`, `ckgpav`; the rest are routines where SPICELIB
   and CSPICE differ, `mxv` and the like); ids identical. Regression tests (equal to
   CSPICE to the last bit, without libm): `test_ck.jl` "chains meet at the reference
   frame of the CK segment", `test_rotations.jl` "TK frames: MATRIX and QUATERNION as
   tkfram_"; both fail with the old code.
3. **Order of operations in the GF state functions.** Done. `_gf_xf` (`src/gf/surface.jl`)
   now calls `rotate_jet` of order 1, `_gf_xfi` sums the six products of each velocity
   component from left to right with the transposed derivative block first (`mxvg_` with
   the matrix of `invstm`), and the normal state of `_gf_ilurate` (`zzilusta`) goes
   through `rotate_jet` (same values but for the sign of a zero). The other state
   transformations of `src/geometry` and `src/aberration.jl` already use `rotate_jet`
   or `_geo_corsxf`; `zzgfcprx`, `zzgffvu`, `zzgfspu` and the other GF routines apply
   3 by 3 matrices only (`mxv_`). Probe against `zzgfssob_`, `zzgfssin_`, `zzilusta_`
   (50 random epochs, Moon and Earth as target, nine corrections, both methods, rays in
   J2000 and in the observer's body frame), glibc check: bit-identical before 11 to 34 %
   of the calls of `zzgfssob` and `zzgfssin`, after 100 % but one velocity component of
   one `CN+S` intercept (one ulp), which comes from the stellar aberration helpers (see
   below). GF families (51, F_GF*, F_ZZGF*): ids identical to the baseline. The search
   routines are now replayed too (item 11). Unit groups gf_solver gf_quantities
   gf_posc gf_surface gf_illum gf_occult gf_fov: 1725 pass. No regression test: the
   change is one ulp, below what a test without glibc's libm can see.
   Found here and fixed on branch plan-rates (item 1): `_zzstelab` (`src/aberration.jl`) uses `_norm3`, `_vperp` and
   `_dvhat`, written for every number type, which round differently from `vnorm`
   (scaled by the largest component), `vperp` and `dvhat`; with glibc they differ from
   CSPICE in most random calls by an ulp, and `spkaps` of F_SPKCOR stays 176 of 432
   calls within 1e-12 (not exact) with glibc. `recazl` of F_AZL: 32 of 684 not exact
   with glibc.
4. **`spkpos` of F_GFDIST (46) and F_LS (7).** Checked: libm (0 with glibc). Closed.
5. **Oracle wrappers that fail.** Done. `CSPICE_EXTRA` of `replay.jl` has wrappers for
   `mxvg`, `mtxvg`, `mxmg`, `mtxmg`, `mxmtg`, `xposeg` (the f2c routines, column major
   as logged), `vtmvg`, `vequ`, `vequg`, `spkaps` (acceleration of 3), `latsrf`,
   `srfnrm`, `dskrb2`, `pltar`, `pltvol`, `polyds` (degree from the log), `azlrec`,
   `drdazl`, `recazl`, `dazldr` (keywords of SPICE.jl) and `swpool`, which is now
   mirrored into CSPICE as the loads are (its watchers decide `cvpool`).
   `TSPICE_ORACLE_ALL` over all 399 families in four runs: `oracle_failures.tsv` empty
   after the last fixes (F_POLYDS, F_POOL, F_VECTORG, F_PLT, F_AZL, F_AZLCPV run again).
   Each new wrapper agrees with Almanac and SPICELIB on every call (`cspice_agreement.tsv`).
   Comparison logic: `cspice_call` returns `NoOracle()` where CSPICE is not asked, so a
   `nothing` of SPICE.jl (not found) is an answer; it agrees only with a result of
   Almanac that is the output `found` alone and false. For a call that SPICELIB failed
   and CSPICE did not, Almanac's outputs are compared with CSPICE's alone (SPICELIB's
   are the caller's values). `within_spread` takes the distance test for floating-point
   outputs only: integers (handles) must agree. Effect on the ids: `vprjp` 1 and
   `vprjpi` 2 of F_PRJP and `inedpl` 1 of F_SLICE (CSPICE's `inedpl_c` is a C routine
   without the zero-normal check of INEDPL) go from `error` to `cspice`; `kdata` 13 and
   `kinfo` 6 of F_KEEPER and `kinfo` 1 of F_KPBUG go from `spread` to `value` (handles 0 against 9: the spread
   was a false pass; the other 137 spread rows of the baseline, of `limbpt`, `sincpt`,
   `subpnt`, `termpt`, stay), and 5 `kdata` of F_DSKKPR stay `value`. The baseline needs an
   update for these.
6. **Writers are checked only by reading back.** Done (branch plan-writers). The
   recorder copies the file of each close from the test program (`dafcls`, `dascls`,
   `spkcls`, `ckcls`, `pckcls`, `ekcls`, `dskcls` at depth 1) whose handle is open for
   write, after the close (`tlog_close_name`, `tlog_close_snapshot`: the dedupe of
   `tlog_snapshot_file`, and past the 256 MB cap nothing is copied or linked). At the
   close `replay.jl` compares Almanac's file with that copy through Almanac's readers
   (`compare_files`): DAF: ND, NI, internal name, comments, and per array the name, DC,
   IC without the two addresses, and the data; DAS: the file record, comments, the last
   addresses and the data of each type; EK: segments, table, columns, attributes and
   every entry. The close is a match or a value mismatch whose name says what differs
   (`array 3 data[85 of 3027]`). `files.tsv` lists each written file: compared, or why
   not. Changes that this needed:
   - writers that are not recorded are logged with a `W` mark and their handles
     (`WRITE_MARKS` of instrument.py: `ekifld`/`ekacl*`/`ekffld`, `spkw19`, `ckw06`,
     `ckw04*`, `spk14a`, `dafws`, `dafrs`, `dafwda`, the `das*` record writers,
     `daswfr`, `dascud`, the `sg*` generic segment writers, the `spks*` subset writers);
     the replay gives up the file they write, so a file that misses a write is never
     compared;
   - an input array longer than the log keeps goes to a gzip binary file named by its
     (byte-wise FNV) hash, within `TSPICE_MAXBIN` (256 MB) per family, named by a `B`
     line that an older replay ignores. `dskw02` logs `spaixi` with its length from
     dsk02.inc; before, every `dskw02` call was truncated and no DSK was compared;
   - `dafbna` is recorded, with its summary unpacked as `dc` and `ic` (packed integers do
     not survive the decimal log); `dafcad` and `pckw20` are recorded;
   - a truncated `dafada` gives up the current DAF (F_CK06 left a short file unflagged);
     an unload of a handle that no loader loaded keeps the handle (F_SGMETA); the
     internal name `replay` that the replay gives for a blank one is not a difference;
     `pckw20`/`spkw20` of degree -1 and `spkw21` of no record are not mapped (Almanac
     takes the degree or line size from the records).

   Measured on the 75 writer families (F_CK0*, F_CKW0*, F_DAF*, F_DAS*, F_DLA, F_DSK*,
   F_DYN05, F_EK*, F_NNEK*, F_PCK*, F_SGMETA, F_SPK*, F_LTIME, F_POOL), recorded again
   into a scratch directory: 11,891 written files closed, 11,573 compared, all equal;
   no real difference in Almanac's written files was found. 318 not compared: 62 a writer
   outside the log (ekifld 15, spkw19 33, ckw06 8, spk14a 2, dafrs 2, dafwda 1, daswfr
   1), 62 a file opened outside the log (tspice utilities: TSTSPK, ZZDDHOPN in F_DAFAH),
   181 past the copy cap or deleted before the close (F_DSKOBJ 177, F_SGMETA 4: KILFIL
   deletes the files before it closes them), 13 a write the replay could not make (3
   still truncated: 2 `dskw02` past the binary budget, 1 `dasadc`, whose strings have no
   binary file; 10 record writers whose log holds fewer values than N records need). The close mismatches
   found on the way were all the recorder's or the replay's (`dafcad` not recorded,
   truncated `dafada`, a weak hash, the `replay` internal name). The other mismatches of
   these families are those of the old baseline, but: new `cspice` rows where calls that
   were truncated are now replayed (F_PCKBUF `tisbod`/`tipbod` 13 -> 26, F_SPKCPV
   `spkacs` 128 -> 164, `spkez` 747 -> 1,135, `sxform` 325 -> 336), and the DSK families
   now replay `dskw02`. `selftest.jl` passes on F_LTIME and F_POOL of the new recording.

   Second round (the writers that were only marked). Recorded and mapped now:
   `spkw19` and `ckw06` (packet sizes per subtype, `tlog_psum`), the EK fast load
   (`ekifld` notes the rows and, from the declarations, each column's entry size and
   null flag; `ekaclc`/`ekacld`/`ekacli` log their values with those lengths; the replay
   keeps Almanac's `EKFastLoad` where SPICELIB passes record pointers), `spk14a` (the
   degree that `spk14b` noted), `dafws`/`dafrs` (the summary unpacked for the DAF of
   the current search, `tlog_dafsum_cur`). Strings and outputs too long for the log now
   also go to binary files, and the budget counts the compressed size (the DSKMI2 index,
   20 million mostly-zero integers, took 80 MB of it per call before). The workspace of
   `dskmi2` is not logged; its index is compared on the header, the coarse grid and the
   lists, which is what is in use (SPICELIB leaves the rest as the caller's array had
   it). Stay marked, because Almanac has no such routine: `dafwda`, `daswfr`, `dascud`,
   the other `das*` record writers, `ckw04*`, `pck03*`, the `sg*` generic segment
   writers and the `spks*` subset writers.

   Measured on the 50 families affected (the writer families of the routines above,
   the DSK families, F_SRFTRN, F_LTIME, F_POOL), recorded into /tmp/rec-writers2, against
   the replay of the same families from /tmp/tsp2: files compared at close 17,020 ->
   17,140, all equal (F_SPKS19 2 -> 27 of 27, F_SPK19 0 -> 8 of 8, F_CK06 8 -> 15 of 16,
   F_ZZDSKBUX 58 -> 110 of 110, F_EK01/03/04 and F_NNEK01/03/04 0 -> all fast-loaded
   files, F_SPK14 2 -> 4 of 4, F_DYN05 0 -> 2 of 6). Not compared: 9,450 past the copy
   cap (F_ZZDSKBSR 9,272, F_DSKOBJ 178), 30 files opened outside the log (tspice
   utilities), 3 still truncated (1 `ckw06`, 1 `dasadc` with BPOS 0, 2 `dskw02` past the
   binary limits), 1 `ckw04b` (marked). One real difference found and fixed in src:
   `dafws` on a DAF open for read threw `:dafinvalidaccess`; SPICELIB and CSPICE's
   `dafws_` signal `DAFILLEGWRITE` after finding the current array (now so, with a
   regression test against CSPICE in test_dafwrite.jl; `dafrs` keeps
   `:dafinvalidaccess`, as `dafrs_` of CSPICE). No other new mismatch in these families
   (the replay faults found on the way, untyped empty vectors and the room of the
   DSKMI2 index, were the replay's). Remaining note: a public `ekacl*` call with an
   untyped empty `entszs` or `nlflgs` (`[]`) throws a MethodError inside Almanac; the
   docstring example uses typed vectors.
7. **Review the `documented` mismatches one by one.** Done. Each was checked against
   CSPICE (probe or oracle) or the glibc check; wrong ones were fixed in Almanac:
   - fixed: `spkpvn` (39) evaluates types 15 and 17 outside the coverage, as
     `spkpvn_c` does (bit-identical); `tisbod`/`tipbod` (6) signal
     COMPETINGEPOCHSPEC, COMPETINGFRAMESPEC and DEGREEOUTOFRANGE as `tisbod_c`; a null
     EK entry (6) compares the null flag only (the values are undefined); a switch
     frame base ID outside Int32 throws `:intoutofrange` (it threw InexactError);
     unequal SCLK partition lists throw `:numpartsunequal`; bad ABCORR/UNITS of
     dynamic frames throw the kinds of CSPICE.
   - reclassed: `dashfs` (765) came from writes outside the log (DASCUD, DASUFS, the
     EK fast load), now marked and skipped; `kinfo`/`kdata` handle 0 is documented.
   - confirmed: `tkvrsn` (399, by design); `tisbod` (134, CSPICE's cubic rate defect,
     Almanac exact); `tangpt`, `sincpt`, `pxform`, `sxform`, `termpt`, `limbpt`,
     `spkcpt`, `subpnt`, `xfmsta`, `dskxsi` (libm); `pipool`/`pdpool`/`pcpool`/
     `lmpool`/`furnsh` `:badvarname` past 32 characters and `dtpool` not found, as
     CSPICE (probe); F_DYN07 `pxform` `:spkinvalidoption`, as `pxform_c`.
   - kept by the maintainer's decision: `gdpool`/`gipool`/`gcpool` throw
     `:typemismatch` for a variable of the other type, where CSPICE returns not found.
8. **The other gaps (about 50 calls).** Done (branch plan-writers), recorded again and
   replayed: 0 `gap` rows in these families but one.
   - `dafopr`/`dafopw` of F_DAFAH (21), F_DAFNN (4), F_DYN05 (2): the logs predated the
     snapshot at the open; with it the replay reads SPICELIB's file: 0.
   - `daffna`/`daffpa` of F_SGMETA (10): SGCTDF writes with `dafbna`, which was not
     recorded, and the test unloads its write handles with `spkuef`, which made the
     replay forget them. Both fixed; `dafwda` is marked: 0.
   - `eknseg` of F_EK01/F_EK03 (2): TSTEK fills the EK with `ekifld`, `ekacl*`,
     `ekffld`; first marked, now recorded and replayed (item 6, second round), so the
     later `ekopr` reads the file that Almanac wrote, equal to SPICELIB's at the close: 0.
     The `dashfs` of F_NNEK01/03/04 (item 7) match too.
   - `dasrfr` of F_DASFR (1): `daswfr` is marked (Almanac has no `daswfr`), and the
     `dasrfr` on that handle is skipped: 0. `dascud` is marked too, so the 760 `dashfs`
     of F_DASCUD (item 7) are skipped.
   - `pcklof` of F_PCKBSR (1) remains: Almanac has no `pcklof`, the replay loads with
     `furnsh`, whose error for a missing file is `NOSUCHFILE` (as `furnsh_c`), where
     SPICELIB's `pcklof` signals `FILENOTFOUND`.
9. **The skipped calls.** Done: `replay.jl` writes `skips.tsv`, the skipped calls per
   routine by reason, and `run.jl` prints the totals. New recording: 222,936 skipped,
   of which 180,444 name a handle the replay does not have (files of writers that are
   not recorded, see item 6), 22,018 follow a load whose array the log could not hold,
   332 hold a truncated array, 14 were made in an error state, 0 have no mapping.
   Second round (item 6): in the 50 affected families, truncated 328 -> 29 (left: failed
   calls with garbage counts in F_GETFOV, `dasadc`/`dasudc` with positions outside the
   string, 6 `dskmi2` and 3 `dskw02` past the binary limits, 1 `ckw06`),
   `after_truncated_load` 22,018 -> 0 (all of F_SRFTRN, a `pcpool` of long strings now
   in a binary file), `unmapped_handle` 74,457 -> 71,901. The unmapped handles that
   remain are mostly F_ZZDSKBSR (69,728: `dskgd`, `dlafps`, `dlabbs` on files that its
   utilities load, and the 9,272 closes past the copy cap, which count here too).
   `skips_family.tsv` gives the counts per family.
10. **`test/test_derivatives.jl`.** Done, approved by the maintainer: `pxder` against
    `sxform` is checked to 2e-15 relative (measured 6e-16; `sxform` rounds as `tisbod_`,
    `pxder` is the derivative of the evaluated matrix). The check of `pxder` against
    CSPICE's `sxform`, which allowed 1e-7, now allows 2e-15 too (measured 3.4e-16).
11. **The GF searches were not replayed.** Done (branch plan-rest). The recorder
    skipped them: their windows are cells `(lbcell:*)`, for which `length_expr` has no
    length, and the workspace `(lbcell:mw, nw)`. `instrument.py` now has `cell_args`
    (and `CELL_ROUTINES`, two hooks in `instrument_routine`): for the eleven searches
    `gfdist gfilum gfoclt gfpa gfposc gfrfov gfrr gfsep gfsntc gfsubc gftfov` a window
    is logged as an array of doubles from CELL(-1): size, cardinality and the CARD
    endpoints (at most 10**6), with the existing `tlog_d`; the workspace is not logged.
    `gfevnt`, `gffove`, `gfocce`, `gfudb`, `gfuds` take callbacks and stay unrecorded.
    `replay.jl` maps the eleven (`gf_search_mapping`): Almanac takes SPICELIB's
    arguments in the same order without MW, NW, WORK and RESULT, and the result window
    is compared by its number of endpoints and each endpoint within the convergence
    tolerance `tol` of the call (`gfstol`, 1e-6 s by default): both solvers end with the
    midpoint of a bracket no wider than `tol` that holds the root, so two endpoints of
    the same root are within `tol`. A SPICELIB error about the sizes that Almanac does
    not take (INVALIDDIMENSION, OUTOFROOM, WINDOWEXCESS, WINDOWTOOSMALL) is skipped as
    `no_counterpart` (new skip reason; run.jl knows it). No CSPICE oracle for them yet
    (`NoOracle`). `gf_agreement.tsv` counts per routine the results equal to SPICELIB's,
    within the tolerance and outside, with the largest difference within.
    Measured on the 22 GF search families recorded again (F_GFBAIL, F_GFBIRP,
    F_GFDIRP, F_GFDIST, F_GFILUM, F_GFLORP, F_GFOCLT, F_GFPA, F_GFPCRP, F_GFPOSC,
    F_GFREFN, F_GFRFOV, F_GFRPRT, F_GFRR, F_GFSCRP, F_GFSEP, F_GFSNTC, F_GFSPRP,
    F_GFSTEP, F_GFSUBC, F_GFTFOV, F_GFXCRP), glibc check: before, 0 search calls
    compared; after, 906 calls: 882 match (175 errors of SPICELIB's kind and 707
    searches), 24 skipped (no counterpart), no other kind; of the 707 searches 687 equal to SPICELIB's window to the last bit, 20 within the
    tolerance (gfilum 10, gfposc 8, gfoclt 2; largest 8.9e-7 s, a last bisection step
    that went the other way), 0 outside; no new mismatch rows (the ids are the old `conics`, `ilumin` cspice and `tkvrsn` rows).
    Remains: re-record everything (instrument.py changed) and update the baseline; a
    CSPICE oracle for the searches (the `gf*_c` wrappers need windows as SpiceCells);
    the same cell logging would record `wn*`, `spkcov`, `ckcov`, `pckcov`, `spkobj`,
    `ckobj`, `dskobj`, `bltfrm`, `kplfrm` if replay.jl mapped them.
