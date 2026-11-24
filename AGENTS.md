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
julia --startup-file=no --project=benchmark benchmark/run.jl                     # benchmarks against CSPICE
```

Run the tests with `--threads=8`; the groups `threads`, `gf_threads` and `ek_threads`
skip with one thread. CI fails on a file that Runic would change. `benchmark/benchmarks.jl`
is the AirspeedVelocity suite that CI runs on a pull request; the other scripts in
`benchmark/` compare Almanac with CSPICE call by call.

## Skills

Use `julia-docstrings`, `julia-manual` and `julia-readme` for documentation,
`julia-doc-reader-review` after a documentation change, `julia-tests` for tests, and
`julia-dependencies` before a test run or a measurement.

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
  `ERROR_KINDS` in `core/errors.jl`: the constructor refuses an unknown kind, and a test
  compares the registry with the source. Public calls make the input checks of CSPICE
  in the order of CSPICE.
- No hidden state. The only global state is the table of open files in
  `files/daf_access.jl`. State that CSPICE hides in saved variables is either an object the
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
it; the folders group the files by area and do not change that order. The shared types
are in `core/types.jl`, `frames/dyntypes.jl`, `dsk/dsk_types.jl` and `ek/types.jl`. They
are immutable. The snapshot-like types are `mutable` with only `const` fields, so that
an atomic reference to one of them is a single word.

| Folder | Files |
|---|---|
| `core/` | `types.jl`, `errors.jl`, `strings.jl`, `taylor.jl` (`Taylor`, `primal`), `bytebuf.jl` (buffers of the number and time formats), `hex.jl` |
| `files/` (DAF, DAS) | `daf.jl` (read), `daf_access.jl`, `daf_write.jl`, `das.jl`, `dla.jl` |
| `kernels/` (contexts, pool) | `context.jl`, `textkernel.jl`, `jsonkernel.jl`, `pool_api.jl`, `pool_write.jl`, `names.jl`, `kernel_info.jl`, `memkernels.jl` (kernels built in memory: `furnsh` of a handle, `load_samples!`, `load_function!`) |
| `spk/` | `spk_layout.jl` (layouts of types 10, 14, 19 and PCK 3, checked at load), `spk_eval.jl` (types 2, 3; PCK segments, `_cheb3`, `find_segment`), `spk_more.jl`, `types/` (`conic.jl`, `interp.jl`, `mda.jl`, `sgp4.jl`, ...), `spk_api.jl`, `aberration.jl`, `spk_segments.jl`, `spk_write.jl`, `tle.jl` |
| `ck/` | `ck.jl`, `ck_api.jl`, `ck_records.jl`, `ck_write.jl` |
| `frames/` | `dyntypes.jl`, `frames.jl`, `rotations.jl`, `dynframes.jl`, `switchframes.jl`, `frame_api.jl`, `frame_info.jl`, `pck_write.jl` |
| `time/` | `time.jl`, `timeparse.jl`, `timeout.jl`, `time_api.jl`, `sclk.jl`, `lst.jl` |
| `dsk/` | `dsk_types.jl`, `dsk.jl`, `dsk_api.jl`, `dsk_write.jl`, `surfaces.jl` |
| `ek/` | event kernels: file, columns, query parser and engine, writers (`write*.jl`) |
| `queries/` | `queries.jl` (prepared queries), `providers.jl`, `table.jl` |
| `math/`, `geometry/`, `gf/` | SPICELIB math (`kernels.jl` holds the `_m_` products; `windows.jl`), geometry, the geometry finder |
| `analysis/` (Orbit layer, ANISE-style) | `orbit.jl` (`Orbit`, `orbit`), `elements.jl` (`element`, `elements`), `brouwer.jl` (Brouwer mean elements, Kepler solver), `build.jl` (`keplerian_orbit`, `with_element`, `with_apsides`, `two_body`), `local_frames.jl` |
| top level | `Almanac.jl`, `ad.jl` (AD rules), `precompile.jl`; the package extensions are in `ext/` (Adapt, ForwardDiff, ChainRulesCore, EnzymeCore, Mooncake, ReverseDiff, Tracker) |

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
- **Remainders.** Where CSPICE calls `d_mod`, use `f2c_dmod`, never `rem`: the exact
  remainder moved the state of SPK type 15 by 3e-12 relative after 3e5 revolutions.
- **Providers.** An `EphemerisProvider` segment has `dtype == PROVIDER_DTYPE` (-1). Any
  code that reads `ks.daf[seg.file]` must test for this first. A `FrameProvider` frame
  has class 7, so dispatch on dynamic frames with `class == 5 || class == 6`, never
  `class >= 5`.
- **Frame chains.** A chain is followed to its first common frame, as in `zzrefch0_`:
  the chain of the first frame keeps `MAX_FRAME_CHAIN` (10) frames and folds later
  links into the last, and a chain of `_FRAME_CHAIN_LIMIT` (1000) links is taken to be
  a cycle (`:framechain`; CSPICE loops). Dynamic frames recurse through
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
- **Segment layouts.** `_daf_parse` checks the meta data of SPK type 10, 14 and 19
  and binary PCK type 3 segments once (`_segment_layout`) and stores the layout in
  `Segment.layout`. A segment that fails the checks loads with a zero layout, as CSPICE
  loads it, and its evaluator makes the checks at each call and throws there. Any other
  code that builds a `Segment` over data that it did not check leaves `layout` zero, and
  a `Segment` must not be evaluated over data other than the data it was loaded from.
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
  Newton step `_carry`). The bisections of `nearpt` and `prop2b` (`_m_nearpt_plain`,
  `_prop2b_root`) are in `_AD_CONSTANT`; their callers take the derivative of the root
  by Newton steps in the number type, one for `Float64` through `_ad_pin(a, p)` (value
  `p`, derivative of `a`, with a rule in the ChainRulesCore, Enzyme and Mooncake
  extensions), so the plain value stays bit for bit. The angles `:ta`, `:aop` and
  `:raan` use `_ad_pin` too: the arc cosine of ANISE for the value
  (`_AD_CONSTANT_ORBIT`), the arc tangent for the derivative. `AlmanacForwardDiffExt` assembles `Dual` results from the jet
  at `primal(et)`; `src/ad.jl` gives `_ad_rule` (value and derivative, for the rules of
  ChainRulesCore, Enzyme and Mooncake) and `_ad_lift` (for the tracked numbers of
  ReverseDiff and Tracker). A new derivative point, or a new path to the epoch that
  bypasses them, needs a method in every extension. Lookups and string parsing on the
  way to a derivative point go into `_AD_CONSTANT`. Enzyme compiles the bodies of the
  derivative points too, and fails when two calls of different result types share a
  storage slot (see `_ssbpos_still` in `aberration.jl`), and applies a rule only to a
  call that is not inlined (`_ad_pin` is `@noinline`); the groups `ad_enzyme`,
  `ad_mooncake`, `ad_zygote` and `ad_tracked` (fixture `test/adfixtures.jl`) cover every
  backend. `element` and the properties of an `Orbit` take the name at run time:
  `getproperty` stays `@inline` (so that Mooncake and Enzyme see one branch), and the
  ChainRulesCore extension has rules for both (so that Zygote compiles one branch, not
  all of `element`).

## Deviations from CSPICE

Every known defect of CSPICE N0067 in a ported routine is listed in
`docs/src/manual/cspice_defects.md`, with what Almanac does and the test group that
covers it. The other differences of the API are in `docs/src/manual/conventions.md`
(windows without a size limit, pool names without blanks, `dafopw` conflicts, the
writers that fix SPICE.jl), in the contexts tutorial (independent contexts) and in the
docstrings of the calls concerned (`ekfind`, `kdata`, `spksfs`). A change that adds or
alters a deviation updates that text in the same commit.

## Tests

- One group per file `test/test_<group>.jl`, listed in `GROUPS` of `test/runtests.jl`.
  Each group must run on its own, so it includes its fixtures with guards of the form
  `isdefined(@__MODULE__, :NAME) || include(...)`.
- CSPICE has global state and is not thread safe. Each group calls `oracle_load` once,
  and SPICE.jl is never called from a spawned task.
- The AD groups take most of the suite time (`ad_enzyme` about 8 min on CI). CI runs
  them and `static` on the coverage job and on a tag; the job `ad-lower-bounds` runs
  `ad_mooncake` and `ad_enzyme` on Julia 1.10 at the lower bounds of `[compat]`.
  `ALMANAC_TEST_EXCLUDE` leaves groups out, `ALMANAC_TEST_SEED` sets the seed, and large
  samples go in the exhaustive tier (`ALMANAC_TEST_FULL=1`, `nsample(fast, full)`).
- CI tests Julia 1.10, which can allocate where later versions do not: after changing a
  hot path, run `quality` on Julia 1.10 too.
- Two checks are outside `test/runtests.jl`. NAIF's test suite replays with
  `julia --project=test test/tspice/run.jl`: the first run downloads NAIF's archive,
  builds it with gcc and records the calls (about 2 GB), and later runs only replay
  them; `test/tspice/README.md` has the details. The `tspice` workflow runs it, in about
  an hour, on each push to `main` that changes the package. The comparison with ANISE is
  `julia --startup-file=no --project=test/anise test/anise/runtests.jl`, which installs
  Python and the `anise` wheel through CondaPkg.
- Kernels are downloaded into `test/data/` by `test/kernels.jl` (`SOURCES` lists the
  URLs and SHA-256 sums). `ALMANAC_KERNEL_DIR` overrides the location,
  `ALMANAC_DOWNLOAD_KERNELS=0` turns the download off, and `ALMANAC_REQUIRE_KERNELS=1`
  makes missing kernels a failure.
- Coverage: run with `--code-coverage=@src`, then `julia --startup-file=no
  test/coverage_summary.jl`.
