# NAIF's test suite (tspice)

This directory replays the calls of tspice, the test suite of SPICELIB, against
Almanac.jl. The calls are recorded through [PEPPER](https://github.com/jacobwilliams/pepper),
and a difference from SPICELIB is checked against CSPICE. The run is not part of the
default test run.

```sh
julia --project=test test/tspice/run.jl
```

## Recording

The first run builds PEPPER with gfortran 16 and fpm, which must be on `PATH`. It then
runs the about 400 tspice families and records their calls into
`~/.cache/almanac-tspice` (about 2 GB, 15 to 40 minutes), with a manifest of the
families and of the PEPPER commit. Later runs only replay the recording, in parallel
processes.

To keep the recording small, a call of a routine that has no state is kept once for
each distinct set of inputs, and a query is kept once between two changes of the loaded
kernels or of the settings. `TSPICE_ALL=1` keeps every call.

The recording also copies each file that a test closes after writing it, so that the
replay compares the file that Almanac.jl wrote with SPICELIB's, by content (arrays and
summaries, DAS data, EK entries), at the close. `files.tsv` (in `results/worker*`) lists
each such file and whether it was compared, or why not. An input array too long for the
log goes to a compressed binary file beside it (at most `TSPICE_MAXBIN` bytes, default
256 MB, per family).

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
