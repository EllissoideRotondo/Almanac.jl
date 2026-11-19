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

1. **Port the CSPICE rate paths that Almanac computes another way.** `zzeprc76`,
   `zzenut80`, `zzmobliq` (of-date frames and the switch frames over them),
   `dnearp` (TARGET_NEAR_POINT vectors), the DT_TWOVEC_IN_DYN case, and the analytic
   rate path of `tisbod` used by `azlcpo` (one ulp in the range rate). Done: the probe
   frames of `test/dynfixtures.jl` are bit-identical to `sxform_c` with the glibc check,
   and `azlcpo` equals `azlcpo_c` with it.
2. **Binary PCK rates.** They come from the Taylor derivative of the evaluation; CSPICE
   evaluates the Chebyshev derivative (`pckeul`, `spke02`). Compare with `sxform_c` on
   a binary PCK, port the CSPICE arithmetic if they differ. Done: bit-identical with
   the glibc check on the binary PCKs of the test kernels.
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
   below). GF families (51, F_GF*, F_ZZGF*): ids identical to the baseline; note that the
   replay compares none of the search routines (`gfsubc`, `gfsntc`, `gfilum` are not
   replayed), so this gate is weak for this item. Unit groups gf_solver gf_quantities
   gf_posc gf_surface gf_illum gf_occult gf_fov: 1725 pass. No regression test: the
   change is one ulp, below what a test without glibc's libm can see.
   Found, not fixed: `_zzstelab` (`src/aberration.jl`) uses `_norm3`, `_vperp` and
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
   `kinfo` 6 of F_KEEPER go from `spread` to `value` (handles 0 against 9: the spread
   was a false pass), and 5 `kdata` of F_DSKKPR stay `value`. The baseline needs an
   update for these.
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
