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

To do first on a new machine: record again (instrument.py changed: items 6 and 8), then
`run.jl --update-baseline` and `selftest.jl`.

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
   Remains: the fast-load EK writers, `spkw19` and `ckw06` need lengths that depend on
   the column declarations or on the packet sizes per subtype; record them to compare
   those files. The coordinator re-records everything (instrument.py changed) and
   updates the baseline (closes now compare, so their counts change).
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
8. **The other gaps (about 50 calls).** Done (branch plan-writers), recorded again and
   replayed: 0 `gap` rows in these families but one.
   - `dafopr`/`dafopw` of F_DAFAH (21), F_DAFNN (4), F_DYN05 (2): the logs predated the
     snapshot at the open; with it the replay reads SPICELIB's file: 0.
   - `daffna`/`daffpa` of F_SGMETA (10): SGCTDF writes with `dafbna`, which was not
     recorded, and the test unloads its write handles with `spkuef`, which made the
     replay forget them. Both fixed; `dafwda` is marked: 0.
   - `eknseg` of F_EK01/F_EK03 (2): TSTEK fills the EK with `ekifld`, `ekacl*`,
     `ekffld`, now marked; the later `ekopr` reads SPICELIB's file: 0. The `dashfs` of
     F_NNEK01/03/04 (item 7) went the same way.
   - `dasrfr` of F_DASFR (1): `daswfr` is marked (Almanac has no `daswfr`), and the
     `dasrfr` on that handle is skipped: 0. `dascud` is marked too, so the 760 `dashfs`
     of F_DASCUD (item 7) are skipped.
   - `pcklof` of F_PCKBSR (1) remains: Almanac has no `pcklof`, the replay loads with
     `furnsh`, whose error for a missing file is `NOSUCHFILE` (as `furnsh_c`), where
     SPICELIB's `pcklof` signals `FILENOTFOUND`.
9. **The 6.4 million skipped calls.** Count them by reason (nested call, error state,
   truncated array, no mapping, unmapped handle) in `tally.tsv`, and list the routines
   with no mapping. Nested calls are expected; any routine that tspice calls at the top
   and the replay never compares is a coverage gap to close.
10. **`test/test_derivatives.jl`.** An exact check became a tolerance of 2e-15 (measured
    6e-16) because `sxform` now rounds as `tisbod_`. Needs the maintainer's approval.
