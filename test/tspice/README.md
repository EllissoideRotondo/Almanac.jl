# NAIF's test suite (tspice)

This directory replays the calls of tspice, NAIF's test suite, against Almanac.jl, and a
difference from the recorded library is checked against CSPICE. The run is not part of
the default test run.

```sh
julia --project=test test/tspice/run.jl
```

## Recording

`TSPICE_SOURCE` selects what the recording runs:

- `cspice`: the 400 families of tspice translated to C by f2c, from NAIF's test suite of
  CSPICE N0067 ([`tspice_c.tar`](https://naif.jpl.nasa.gov/pub/naif/misc/tspice/N0067/PC_Linux_64bit/tspice_c.tar)),
  against the CSPICE of the same archive. It needs gcc only. `cspice/build_cspice.sh`
  downloads the archive (checked by its SHA-256), writes the wrappers with
  `cspice/instrument_c.py` and builds the program; only the calls of the test program
  (the families and NAIF's test and support utilities) go through the wrappers.
  `cspice/signatures.tsv` holds the kind, direction and dimensions of the arguments of
  each recorded routine, which the C prototypes do not give: after a change of
  `routines.txt`, write it again with `cspice/signatures.py` from the Fortran sources of
  SPICELIB (`tspice.tar` in the same directory of NAIF).
- `pepper` (the default): the families of tspice through
  [PEPPER](https://github.com/jacobwilliams/pepper)'s SPICELIB interface, built with
  gfortran 16 and fpm, which must be on `PATH`.

The first run builds the program, runs the families and records their calls into
`~/.cache/almanac-tspice` (about 2 GB, 15 to 40 minutes), with a manifest of the
families and of what made the recording. Later runs only replay the recording, in
parallel processes.

To keep the recording small, a call of a routine that has no state is kept once for
each distinct set of inputs, and a query is kept once between two changes of the loaded
kernels or of the settings. `TSPICE_ALL=1` keeps every call.

The recording also copies each file that a test closes after writing it, so that the
replay compares the file that Almanac.jl wrote with SPICELIB's, by content (arrays and
summaries, DAS data, EK entries), at the close. `files.tsv` (in `results/worker*`) lists
each such file and whether it was compared, or why not. Past 256 MB of copies in a
family, only the copy of the last close of each file is kept: a test that appends to
one file and closes it after each segment has its file compared once, with every
segment, at the last close. An input array too long for the log goes to a compressed
binary file beside it (at most `TSPICE_MAXBIN` bytes, default 256 MB, per family).

The log also names the file of each handle of a file open for reading. A test gets some
handles from routines outside the recording (the DSK segment search of ZZDSKBSR), and
the replay finds them by their file, if Almanac.jl has it loaded.

Windows and sets (cells) are logged as their size, their cardinality and their elements.
The GF searches that take procedures (GFEVNT, GFFOVE, GFOCCE) are replayed with the step
that GFSSTP or the last search set and with Almanac.jl's `gfrefn`, which are the
procedures that every tspice call passes. Not replayed:

- GFUDB, GFUDS, UDDC and UDDF, whose functions are written in the test program;
- writers without a counterpart in Almanac.jl (the SPKSnn subset writers, DAFWDR and
  other record-level writers). A file that one of them wrote is read from
  the copy of SPICELIB's file, and its later writes are skipped;
- calls on a file that the test opened with routines outside the recording;
- DASCUD and the DAS calls of F_DASCUD after it: Almanac.jl keeps a DAS file in memory
  and writes its directories at the close, so it has no step that updates them;
- calls of the writers that Almanac.jl's arguments cannot express, such as a degree
  that disagrees with the length of the data (Almanac.jl takes the degree from the
  records).

`TSPICE_WHY=1` writes `unmapped.tsv` (in `results/worker*`): for each family and routine,
the calls skipped for want of a file, and why the replay has none.

## When a run stops or fails

A run without a list of families is the full check.

It **stops** (exit code 2), without a comparison with the baseline, if a family is missing from the
recording, failed to record, or fails to replay to its end.

It **fails** (exit code 1) if:

- a mismatch is not in the baseline `baseline_ids.tsv.gz`. A mismatch is any call that
  is not an exact match: a value within a tolerance (`near`), a value that only CSPICE
  shares (`cspice`), a different value, a different error kind, or an error of the
  replayed test script itself. Each mismatch is identified by family, routine, kind and a
  hash of its inputs, so a new mismatch fails even when another one went away, and a call
  that moves from a match to `near` or `cspice` fails too. A call counts as `cspice` only
  where CSPICE itself differs from SPICELIB;
- a routine has more calls skipped, or compared with a tolerance, than `coverage.tsv`
  records, or fewer calls compared.

The mismatches are written to `results/replay.tsv` and `results/ids.tsv`.

`selftest.jl` checks the gate itself. On a copy of the recording it removes a mismatch
from the baseline, scales Almanac.jl's outputs by 1 + 1e-10, changes the kind of its
errors, deletes a log, damages a log and marks a family as failed, and requires the run
to fail each time.

## The baseline

Each routine of `baseline.tsv` has a class and a note:

- `documented`: a difference that a docstring or the manual states;
- `cspice`: Almanac.jl follows CSPICE where SPICELIB differs;
- `gap`: a call that the replay cannot reproduce, such as a read of a file that the
  test wrote with routines outside the recording;
- `unexplained`.

## Options

- `--record` records again.
- `--update-baseline` makes the full run the new baseline and keeps the classes.
- A list of families (`F_SPK08 F_CK06`) restricts the run to them, without the coverage
  check.
