# Almanac.jl

A native Julia port of CSPICE N0067. The C source of CSPICE N0067 is the reference, and
CSPICE (called through SPICE.jl) is the test oracle: match its conventions (km, km/s,
radians) and its arithmetic, down to the last bit where a test compares bits.

## Commands

```sh
julia --startup-file=no --threads=8 test/parallel.jl [group ...]                 # tests, several processes
julia --startup-file=no --threads=8 --project=test test/runtests.jl [group ...]  # tests, one process
julia --project=docs docs/make.jl                                               # manual and doctests
julia --project=@runic -m Runic --inplace .                                     # format (CI checks it)
```

Run the tests with `--threads=8`; the groups `threads`, `gf_threads` and `ek_threads`
skip with one thread. The code is formatted with [Runic](https://github.com/fredrikekre/Runic.jl)
(`julia --project=@runic -e 'using Pkg; Pkg.add("Runic")'` once); CI fails on a file that
Runic would change.

## Skills

Use the skills whenever you change an exported name, a docstring, the README, the
manual or a test. The global skills hold the general rules, and the project skills in
`.claude/skills/` add the rules of this repository:

- `julia-docstrings`, `julia-manual`, `julia-readme`, with `almanac-docs`: docstrings,
  the Documenter manual, the README.
- `julia-doc-reader-review` (agent `julia-doc-reader`): independent reader review after
  a documentation change; `almanac-docs` says what the reader gets.
- `julia-tests`, with `almanac-tests`: test design, regression tests, trimming.
- `julia-dependencies`: every environment (root, `test`, `docs`) on the latest versions
  that `[compat]` allows, checked before a test run or a measurement.

## Conventions that hold everywhere

- Epochs are ET (TDB seconds past J2000). Evaluators are generic over `T<:Real`, so
  `Float64`, `ForwardDiff.Dual` and the internal `Taylor{N,T}` all pass through. Every
  discrete decision (segment search, record index, branch) is made on
  `primal(et)::Float64`. `Taylor` defines no comparisons, so code that branches on a
  number must call `primal`.
- A "jet of order N" is an `NTuple{N+1}` of time derivatives, orders 0 to N.
- A `KernelSet` is an immutable snapshot of everything that is loaded. A query reads
  the snapshot of its context once and takes no lock. A change (`furnsh`, `unload`,
  `kclear`, a pool edit) builds a new snapshot under the lock of the context and swaps
  it in atomically. A stored pool vector is never mutated, because older snapshots
  share it.
- Every public call that needs kernel data has three forms: `f(ks::KernelSet, ...)`,
  `f(ctx::SpiceContext, ...)`, and `f(...)` on `default_context()`.
- Hot paths allocate nothing and return StaticArrays values.
- Errors: `throw(SpiceError(kind, msg))`, where `kind` is the CSPICE short message in
  lower case (`:invaliddegree` for SPICE(INVALIDDEGREE)). A new kind must be added to
  `ERROR_KINDS` in `errors.jl`: the constructor refuses an unknown kind, and a test
  compares the registry with the source. Public calls make the input checks of CSPICE
  in the order of CSPICE.
- No hidden state. The only global state is the table of open files in
  `daf_access.jl`. State that CSPICE hides in saved variables is either an object the
  caller owns (`DAFHandle`, `DASHandle`, `EKResult`, `EKFastLoad`, `SPK14Writer`,
  `CKW04Writer`, `PCK03Writer`) or an atomic field of the `SpiceContext` (`ks`,
  `timedefs`, `dsktols`, `gftol`).
- Where CSPICE computes with 32-bit integers, use `_cint`, `_wrap32`, `_dnint` and
  `_nint` from `math/kernels.jl`, so the results are the same. Integer pool values are
  rounded as `gipool_` rounds them, through `_gipool_value` and `_gipool_first`.
- A call that has a "found" flag in CSPICE returns `nothing` when the flag is false.
  Results are returned, not written to output arguments.

## Source layout

`src/Almanac.jl` sets the include order, and a file uses only files that come before
it. The shared types are in `types.jl`, `dyntypes.jl`, `dsk_types.jl` and
`ek/types.jl`. They are immutable. The snapshot-like types are `mutable` with only
`const` fields, so that an atomic reference to one of them is a single word.

| Area | Files |
|---|---|
| Binary kernels (DAF) | `daf.jl` (read), `daf_access.jl`, `daf_write.jl`, `memkernels.jl` (kernels built in memory: `furnsh` of a handle, `load_samples!`, `load_function!`) |
| SPK evaluation | `spk_eval.jl` (types 2, 3), `spk_more.jl`, `spk_types/`, `conic.jl`, `interp.jl`, `mda.jl` |
| Text and JSON kernels, pool | `textkernel.jl`, `jsonkernel.jl`, `pool_api.jl`, `pool_write.jl` |
| Names, frames, rotations | `names.jl`, `frames.jl`, `rotations.jl`, `dynframes.jl`, `switchframes.jl`, `frame_info.jl` |
| Time, clocks | `time.jl`, `timeparse.jl`, `timeout.jl`, `time_api.jl`, `sclk.jl`, `lst.jl` |
| CK | `ck.jl`, `ck_api.jl`, `ck_records.jl`, `ck_write.jl` |
| Context, queries | `context.jl`, `spk_api.jl`, `frame_api.jl`, `aberration.jl`, `queries.jl`, `providers.jl` |
| Writers | `spk_write.jl`, `ck_write.jl`, `pck_write.jl`, `dsk_write.jl`, `ek/write*.jl` |
| DAS, DSK, EK | `das.jl`, `dla.jl`, `dsk.jl`, `dsk_api.jl`, `surfaces.jl`, `ek/` |
| Math, geometry, events | `math/`, `geometry/`, `gf/`, `windows.jl` |
| Tables, AD, extensions | `table.jl`, `ad.jl`, `ext/` (Adapt, ForwardDiff, ChainRulesCore, EnzymeCore, Mooncake, ReverseDiff, Tracker) |
| Orbit layer (ANISE-style analysis) | `analysis/`: `orbit.jl` (`Orbit`, `orbit`), `elements.jl` (`element`, `elements`), `brouwer.jl` (Brouwer mean elements, Kepler solver), `build.jl` (`keplerian_orbit`, `with_element`, `with_apsides`, `two_body`), `local_frames.jl` |

## Invariants that are easy to break

- **Priority order.** `KernelSet.spk`, `pck`, `ck` and `dsk` are in the CSPICE search
  order: the later file first, and within a file the later segment first. The lists of
  `spk`, `pck` and `ck` are `SegmentList`s: from 16 segments on, `find_segment` answers
  from an index of the segment bounds, built once with the snapshot, and must give the
  index of the linear scan. A CK search with a tolerance scans all the segments.
- **3 by 3 products.** Use `_m_mxv`, `_m_mtxv`, `_m_mxm`, `_m_mtxm` or `_m_mxmt`, never
  `*`. The StaticArrays `*` lets the compiler fuse operations, so its last bit changed
  between Julia versions.
- **Chebyshev sums.** `_cheb3` must match `chbval_`/`chbint_` to the last bit for a
  plain epoch, which matters for angles of many turns.
- **Remainders.** Where CSPICE calls `d_mod`, use `f2c_dmod` (or `_dyn_dmod`), never
  `rem`: the exact remainder moved the state of SPK type 15 by 3e-12 relative after 3e5
  revolutions.
- **Providers.** An `EphemerisProvider` segment has `dtype == PROVIDER_DTYPE` (-1). Any
  code that reads `ks.daf[seg.file]` must test for this first. A `FrameProvider` frame
  has class 7, so dispatch on dynamic frames with `class == 5 || class == 6`, never
  `class >= 5`.
- **Frame chains.** A chain is followed to its first common frame, as in `refchg_`,
  for at most `MAX_FRAME_CHAIN` (10) links. Dynamic frames recurse through
  `ks.inner`, which is `_DYN_LEVELS` (2) deep, and throw `:recursiontoodeep` below
  that.
- **CK rates.** `_CKState()` takes its derivative from the angular velocity of the
  segment, as `sxform` does, which is not the time derivative of the rotation.
  `_CKExact()` and `_exact_frame_jet` give the exact derivative. `frame_jet` picks the
  first for chains that contain class 3, 5 or 6 frames.
- **PCK type 3 rates.** A segment holds series of the angles and of the angular
  velocity. `sxform` at a plain epoch takes the derivative from the angular velocity
  series (`_pck03_state`, as `pcke03`), through the state path of class 2 frames; an AD
  or `Taylor` epoch differentiates the angle series (`_pck03_angles`). W is not reduced
  by whole turns, as in `pcke03`.
- **SPK type 3** (separate velocity series): order k >= 1 of the jet is derivative k-1
  of the velocity series, as in `spkder`.
- **Kernels in memory.** A handle of `dafops` (`spkops`, `ckops`, `pckops`) has an
  `IOBuffer` and `path == ""`, and is never in the table of open files. `furnsh` of a
  handle copies its bytes and parses them with `_daf_parse`, as `open_daf` does a file,
  so writing to the handle later never changes a loaded kernel. The body names of a
  kernel (`bodies`, on a `DAF` or an `EphemerisProvider`) rank below the pool and above
  `boddef`, the later kernel first.
- **Memory maps on Windows.** Windows cannot delete or replace a mapped file, so a map
  that is needed only temporarily is released with `_unmap!` (or `_das_release!`), not
  left to the garbage collector. The maps of a `KernelSet` stay mapped.
- **Writers.** Each writer makes the SPICELIB checks in the same order, and the data it
  writes are equal word for word to what CSPICE writes from the same input. Bytes go
  in before their summary, so a reader always sees a valid file.
- **Parsers of dynamic and switch frames** never throw. They record `errkind`, and the
  error is thrown when the frame is used.
- **Geometry.** Each call is split into a public function, which makes the CSPICE
  checks, and a private function, which does the computation; the geometry finder
  calls the private one. `occultation.jl`, `limbterm.jl`, `fov.jl` and `lst.jl` take
  plain epochs only, through `_plain_epoch`.
- **Geometry finder.** `gf/solver.jl` knows no geometry. A new quantity is a new
  subtype of `GFQuantity` with three methods: `_gf_value`, `_gf_decreasing` and
  `_gf_label`. The "is decreasing" test is the CSPICE formula, not the exact
  derivative, so that the brackets match those of CSPICE.
- **Automatic differentiation.** No AD package differentiates the evaluators. The
  derivative points are `spk_jet`, `frame_rotation`, `frame_jet`, `abcorr_jet`,
  `abcorr_position` and `abcorr_state`; the derivative of each is the next order of its
  jet. The Orbit layer has one more, the Kepler solver `_ma_to_ta(M, e)` of
  `analysis/brouwer.jl`, whose rules take the partials of the implicit function at the
  converged anomaly (ForwardDiff and the tracked numbers get them from the carried
  Newton step `_carry`). `AlmanacForwardDiffExt` assembles `Dual` results from the jet
  at `primal(et)`; `src/ad.jl` gives `_ad_rule` (value and derivative, for the rules of
  ChainRulesCore, Enzyme and Mooncake) and `_ad_lift` (for the tracked numbers of
  ReverseDiff and Tracker). A new derivative point, or a new path to the epoch that
  bypasses them, needs a method in every extension. Lookups and string parsing on the
  way to a derivative point go into `_AD_CONSTANT`. Enzyme compiles the bodies of the
  derivative points too, and fails when two calls of different result types share a
  storage slot (see `_ssbpos_still` in `aberration.jl`); the groups `ad_enzyme`,
  `ad_mooncake`, `ad_zygote` and `ad_tracked` (fixture `test/adfixtures.jl`) cover every
  backend. `element` and the properties of an `Orbit` take the name at run time:
  `getproperty` stays `@inline` (so that Mooncake and Enzyme see one branch), and the
  ChainRulesCore extension has rules for both (so that Zygote compiles one branch, not
  all of `element`).

## Deliberate deviations from CSPICE

Every known defect of CSPICE N0067 in a ported routine is listed in
`docs/src/manual/cspice_defects.md`, with what Almanac does (copies or fixes) and the
test group that covers it. When you find a defect, or change how Almanac treats one,
update that page in the same change.

- `SpiceContext` objects are independent. Each one has its own pool watchers
  (`swpool`/`cvpool`), time defaults, DSK tolerances and GF tolerance.
- Windows have no size limit, and searches return a new `Window`, so there are no
  workspaces or result cells.
- `ekfind` returns an `EKResult`, which the fetch calls take as their first argument.
- `kdata`/`kinfo` return handle 0, because loaded kernels have no DAF handle. `spksfs`
  returns a path in place of a handle.
- CK type 6 uses the stored clock rate of the next mini-segment, which `ckr06_` rounds
  to a whole number.
- `spkpvn` checks the time bounds for every SPK type but 15 and 17, whose single set of
  elements it propagates to any epoch, as CSPICE does.
- `_spks19` throws `:spkstructureerror` where `spks19` reads `data[0]`.
- `ckw02`/`ckw03` take the full record count, and `pckw02` takes the degree from the
  record length. Both fix defects of SPICE.jl 0.2.5.
- The EK query engine fixes three defects of N0067 (`zzekweed`, `zzekillt`/`zzekille`,
  `zzekjsqz`). `compat = true` reproduces them for the tests.
- `dafopw` refuses a write open of a file that another handle has open
  (`:fileopenconflict`).

## Tests

- One group per file `test/test_<group>.jl`, listed in `GROUPS` of `test/runtests.jl`.
  Each group must run on its own, so it includes its fixtures with guards of the form
  `isdefined(@__MODULE__, :NAME) || include(...)`.
- CSPICE has global state and is not thread safe. Each group calls `oracle_load` once,
  and SPICE.jl is never called from a spawned task.
- Time, measured on CI (ubuntu, Julia 1, one process): every group takes about 17 min,
  of which the AD groups take about 12 min (compiling the backends; `ad_enzyme` alone
  8 min) and `static` 40 s; the other groups, which the CI test jobs run, take about
  4.5 min. A group that grows states its new time in its header. Large samples go in the
  exhaustive tier (`ALMANAC_TEST_FULL=1`, `nsample(fast, full)`). `ALMANAC_TEST_SEED`
  sets the random seed. `ALMANAC_TEST_EXCLUDE` (group names separated by spaces) leaves
  groups out of a run without group names; CI runs the AD groups and `static` on the
  coverage job (and two AD groups on Julia 1.10) only, and every group on a tag or a
  manual run. The job `ad-lower-bounds` runs `ad_mooncake` and `ad_enzyme` on Julia
  1.10 with EnzymeCore and Mooncake pinned to the lower bounds of `[compat]`.
- Kernels are downloaded into `test/data/` by `test/kernels.jl` (`SOURCES` lists the
  URLs and SHA-256 sums). `ALMANAC_KERNEL_DIR` overrides the location,
  `ALMANAC_DOWNLOAD_KERNELS=0` turns the download off, and `ALMANAC_REQUIRE_KERNELS=1`
  makes missing kernels a failure.
- Coverage: run with `--code-coverage=@src`, then `julia --startup-file=no
  test/coverage_summary.jl`.
