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
```

The glibc check: run `replay.jl` through a script that first redefines, for `Float64`,
`Base.Math.sin`, `cos`, `tan`, `asin`, `acos`, `atan` (one and two arguments), `sinh`,
`cosh`, `tanh`, `exp`, `log`, `log10`, `expm1`, `log1p`, `^` as `ccall`s into
`libm.so.6`, and `sincos(x) = (sin(x), cos(x))`, then `include`s `replay.jl`. A
mismatch that remains with glibc is not a libm effect.

## State at commit 8d5eb68 plus the recorder fix

Full replay (399 families, CSPICE oracle on):

| kind | count | meaning |
|---|---|---|
| match | 9,607,829 | equal to SPICELIB within 1e-12 |
| cspice | 163,581 | differs from SPICELIB, equal to CSPICE, CSPICE differs from SPICELIB |
| near | 2,895 | within 1e-9 of SPICELIB; CSPICE agrees with SPICELIB |
| spread | 178 | SPICELIB and CSPICE differ; Almanac within twice their distance |
| value | 39,728 | larger difference |
| error | 103 | other error, or an error where SPICELIB had none |
| skipped | 6,443,171 | not replayed (nested calls, calls in an error state, unmapped) |

Verified on F_DYN01, F_XFMSTA, F_GFDIST, F_LS with the glibc check: near 1,352 -> 0,
value 18 -> 4 (the 4 are `tkvrsn`). The large value rows of F_DYN01 (to 3e-7 of the
vector magnitude) are libm: the velocity of Mars from Phobos in PHOBOS_RADIAL is a
small residue of a cancellation, which amplifies a one-ulp difference of `sin`.
The agents found the same for `conics`, `tangpt`, `sincpt`, `spkcpt`, `subpnt`,
`subslr`, `evsgp4` and the type 5 SPK (`prop2b`). Note: on the arguments checked,
glibc's `sin` and `cos` were correctly rounded where Julia's were not (one ulp).

Agreement with CSPICE on every call of 20 families (882,315 calls):

| | exact | within 1e-12 | 1e-12 to 1e-9 | beyond 1e-9 |
|---|---|---|---|---|
| Almanac | 721,956 | 151,808 | 6,514 | 2,037 |
| SPICELIB (PEPPER) | 659,811 | 155,548 | 64,245 | 2,711 |

## Fixed since the last baseline: the 38,164 `dasrdc` gaps

They were a recorder bug, not a limit of the replay. `F_DASA2L` writes strings that
hold every character code, line feeds included, and the line-based log split them, so
Almanac wrote different data and every read differed. `instrument.py` now logs a
string with a character outside printable ASCII as `%` and hex codes; `replay.jl`
decodes it. F_DASA2L re-recorded and replayed: `dasrdc` 87,500 match, 0 value. Only
F_DASA2L had such strings (scan of all logs).

To do first on a new machine: record again (instrument.py changed; at least F_DASA2L),
then `run.jl --update-baseline` and `selftest.jl`.

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
   geometric, `geometry/constobs.jl`, `_dyn_apo`/`_dyn_spkez`); `recazl` takes the
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
   `spkcvt` now all). Left: 57 `sxform` and 137 `pxform` calls within 1e-12 but not
   exact, 1 `spkgeo` and 1 `spkgps`; not looked at.
3. **Order of operations in the GF state functions.** `src/gf/surface.jl` `_gf_xf` and
   `src/gf/illumination.jl` apply `R*v + dR*r`; CSPICE applies the 6x6 matrix as
   `mxvg_` (derivative block first, sums from left to right), as `rotate_jet` of order 1
   now does. Done: the GF families show no new near or value, and the changed functions
   equal their CSPICE counterparts on random input.
4. **`spkpos` of F_GFDIST (46) and F_LS (7).** Checked: libm (0 with glibc). Closed.
5. **Oracle wrappers that fail.** SPICE.jl's `spkaps` (`accobs` of 6), `latsrf`,
   `srfnrm`, `vequ`, `mxvg`, `mxmg` do not take the arguments as the replay passes
   them, so a mismatch of these routines is never checked against CSPICE. Add wrappers
   in `CSPICE_EXTRA`. Done: `oracle_failures.tsv` of a `TSPICE_ORACLE_ALL` run is empty.
   Also `vprjpi` of F_PRJP (2 errors): SPICELIB signals NONUNITNORMAL, CSPICE takes
   the normal; the comparison of a call that SPICELIB failed with outputs CSPICE
   returned is not handled, so it stays `error`.
6. **Writers are checked only by reading back.** A writer call compares nothing; its
   file is checked only if a later recorded call reads it. Add a comparison of the file
   that Almanac writes with SPICELIB's at close (`dafcls`, `dascls`, `spkcls`, `ckcls`,
   `pckcls`, `ekcls`, `dskcls`): snapshot SPICELIB's file at its close in the recorder,
   and compare the arrays and summaries (not the bytes: Almanac may lay out records
   differently) in the replay. Done: each writer family compares every segment it
   wrote; the audit shows no writer whose output is never checked.
7. **Review the 1,607 `documented` value and error mismatches one by one.** The class
   says a docstring or the manual states the difference; that has to be true and the
   difference has to be right. Known so far:
   - `dashfs` 765: 760 are F_DASCUD, which calls `DASCUD` and `DASRRI`, internal
     SPICELIB routines that Almanac does not have, so its files differ; class `gap`,
     not `documented`. 3 are F_NNEK01/03/04 (`free` 212 vs 7) and 1 is F_DASA2L
     (`lastla` of characters 1,000,000 vs 10 after a write the log cannot hold):
     not explained yet.
   - `tkvrsn` 399: the version string names Almanac. By design.
   - `tisbod` 142, `tipbod` 4: SPICELIB's rate of a cubic angle term is wrong; Almanac
     gives the exact derivative (manual, CSPICE defects). Also "the degree check":
     confirm which calls that covers.
   - `spkpvn` 39: Almanac refuses an epoch outside the segment. Check what CSPICE
     does; if `spkpvn_c` evaluates, decide with the maintainer.
   - `tangpt` 186, `sincpt` 19, `pxform` 9, `sxform` 5, `termpt` 7, `limbpt` 5,
     `spkcpt` 3, `subpnt` 2, `xfmsta` 1: libm. Confirm each with the glibc check.
   - `dtpool` 5, `pipool` 1, `gdpool` 2, `ekrcec`/`ekrced`/`ekrcei` 6: confirm against
     CSPICE (the oracle) and the docstrings.
   Done: every `documented` routine of `baseline.tsv` has a note that is checked, and
   any that is wrong is fixed in Almanac or reclassed.
8. **The other gaps (about 50 calls).** `dafopr`/`dafopw` of F_DAFAH, F_DAFNN, F_DYN05
   (files made with Fortran I/O outside the log), `daffna`/`daffpa` of F_SGMETA,
   `eknseg` of F_EK01/F_EK03, `dasrfr` of F_DASFR (`DASWFR` not recorded), `pcklof` of
   F_PCKBSR. For each: record the missing routine, or snapshot the file at the open
   (the recorder snapshots opens at depth 1 only), so that the replay reads
   SPICELIB's file. Done: 0 `gap` rows, or each remaining one has a stated reason
   that the replay cannot reproduce it.
9. **The 6.4 million skipped calls.** Count them by reason (nested call, error state,
   truncated array, no mapping, unmapped handle) in `tally.tsv`, and list the routines
   with no mapping. Nested calls are expected; any routine that tspice calls at the top
   and the replay never compares is a coverage gap to close.
10. **`test/test_derivatives.jl`.** An exact check became a tolerance of 2e-15 (measured
    6e-16) because `sxform` now rounds as `tisbod_`. Needs the maintainer's approval.
