# Analysis layer part 1 (`Orbit`, elements, local frames) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** An `Orbit` type in Almanac.jl with every element and local frame of ANISE's `Orbit`, checked against CSPICE, GMAT values and identities in the default run, and against the ANISE Python package in an opt-in group.

**Architecture:** New code in `src/analysis/`, included after `src/table.jl`. `Orbit{T}` holds a state, epoch, center, frame, GM and radii; `element(o, name)` is one type-stable `if` chain over private `_el_<name>(o)` functions, and `getproperty` forwards to it. Formulas are ported from ANISE's Rust source at commit `b3aaadc`, function by function, including ANISE's error and fallback cases.

**Tech Stack:** Julia ≥ 1.10, StaticArrays, Test; SPICE.jl (oracle); AllocCheck; DifferentiationInterface + six AD backends (last task); PythonCall + CondaPkg (opt-in `anise` group only).

**Spec:** `docs/superpowers/specs/2026-11-15-analysis-orbit-design.md`

**ANISE reference source:** `https://github.com/nyx-space/anise` at `b3aaadc` (`git clone https://github.com/nyx-space/anise && git -C anise checkout b3aaadc`). Paths below like `astro/orbit.rs` are under `anise/src/`; `tests/astro/orbit.rs` is under `anise/tests/`.

## Global Constraints

- Julia 1.10 is the oldest supported version; no new runtime dependency of Almanac (`[deps]` unchanged).
- Units km, s, rad (and products). Never degrees in the API.
- `[ctx,]` optional first argument on every call that reads kernels; one docstring on the generic function (`almanac-docs` skill).
- Errors are `SpiceError(kind)` with kinds from `Almanac.ERROR_KINDS`: `:degeneratecase` where ANISE fails on an undefined quantity, `:badarg` for bad arguments, the kernel kinds that `bodvcd`/`spkezr` already throw.
- Every new function generic in the number type `T` of state and epoch; GM and radii stay `Float64`.
- Default suite under 120 s in one process; each new group states its time (`almanac-tests` skill). Tolerances are measured differences with a stated reason, recorded with `measured(...)`.
- Random input from the seeded generator (`SEED` in `test/helpers.jl`); sample sizes through `nsample(fast, full)`.
- Docstrings of new exports say "This call has no CSPICE counterpart." and name the ANISE function it matches.
- New exports exactly: `Orbit, orbit, element, elements, keplerian, keplerian_orbit, with_element, with_apsides, two_body, vinf_periapsis, vinf_turn_angle, local_frame, local_frame_state, relative_state`.

## Review Focus

1. Inbound hyperbolic state (`r·v < 0`, `ecc > 1`): accepted; `:ta` in (π, 2π) wrapped as ANISE main does, `:ma` the negative hyperbolic mean anomaly `e sinh H − H`. Test in Task 3.
2. Exactly circular (`ecc == 0`) and exactly equatorial (`inc == 0` and `inc == π`): no NaN ever leaves `element`; each quantity either returns ANISE's fallback or throws `:degeneratecase`. Test in Task 2.
3. Center without `RADII` in the pool (barycenter 3, a spacecraft): construction and Keplerian elements work; `:altitude`, `:latitude`, `:height` and the Brouwer elements throw `SpiceError(:badarg)` naming the missing radii (radii are NaN in the struct; the check is in the element). Test in Task 4.
4. Near-parabolic `|ecc − 1| < 1e-11`: anomalies and `:sma` throw `:degeneratecase`, not `Inf`. Test in Task 2.
5. Unknown name: `element(o, :foo)` and `o.foo` throw `SpiceError(:badarg)` whose message lists `elements()`. Test in Task 1.

---

## File structure

| File | Responsibility |
|---|---|
| `src/analysis/orbit.jl` | `Orbit` type, `orbit`, `Orbit` constructors, `element`, `elements`, `getproperty`, `propertynames`, `show` |
| `src/analysis/elements.jl` | `_el_*` functions: Cartesian, Keplerian, geodetic, hyperbolic/B-plane, equinoctial; `keplerian`, `vinf_periapsis`, `vinf_turn_angle` |
| `src/analysis/brouwer.jl` | Brouwer mean short elements (port of `astro/orbit_mean_elements.rs`) |
| `src/analysis/build.jl` | `keplerian_orbit`, `with_element`, `with_apsides`, `two_body` |
| `src/analysis/local_frames.jl` | `local_frame`, `local_frame_state`, `relative_state` |
| `src/Almanac.jl` | includes after `table.jl`; the export line |
| `test/test_analysis_orbit.jl` | default group `analysis_orbit` (Tasks 1–6) |
| `test/gmat_cases.jl` | GMAT/Vallado/Davis reference cases ported from ANISE's `tests/astro/orbit.rs` |
| `test/test_analysis_derivatives.jl` | default group `analysis_derivatives` (Task 9) |
| `test/anise/Project.toml`, `test/anise/CondaPkg.toml`, `test/anise/runtests.jl` | opt-in ANISE comparison (Task 8) |
| `docs/src/howto/orbits.md`, `docs/src/reference/orbits.md`, `docs/make.jl`, `docs/src/explanation/coverage.md` | manual (Task 7) |

---

### Task 1: `Orbit`, constructors, Cartesian elements, name dispatch

**Files:**
- Create: `src/analysis/orbit.jl`, `src/analysis/elements.jl`, `test/test_analysis_orbit.jl`
- Modify: `src/Almanac.jl` (include both after `include("table.jl")`; add the export line), `test/runtests.jl` (`"analysis_orbit"` after `"queries"` in `GROUPS`)

**Interfaces:**
- Produces:
  - `struct Orbit{T<:Real}` with fields `state::SVector{6,T}, et::T, center::Int32, frame::Int32, gm::Float64, radii::SVector{3,Float64}` (radii NaN when unknown).
  - `Orbit(state, et, center::Integer, frame::Integer, gm::Real, radii=SVector(NaN,NaN,NaN))` — no kernels.
  - `Orbit([ctx,] state, et, center::NameOrID, frame::NameOrID)` — `gm = bodvcd(ks, center, "GM")[1]`, radii from `bodvcd(ks, center, "RADII")` when `bodfnd`, else NaN.
  - `orbit([ctx,] targ, et, ref, abcorr, obs)` — `spkezr` state, center `obs`, frame `ref`.
  - `element(o::Orbit, name::Symbol) -> T` (scalars only); `elements() -> NTuple{N,Symbol}` in the order of the spec's table.
  - `Base.getproperty(o::Orbit, s::Symbol)`: fields, then `:position, :velocity, :hvec, :evec` (`SVector{3,T}`), then `element`.
  - `Base.propertynames(o::Orbit)`: fields, vector properties, `elements()`.
  - Private `_el_<name>(o::Orbit)` for each element; `_degenerate(name, what)` throwing `SpiceError(:degeneratecase, ...)`.

- [ ] **Step 1: Write the failing tests** (group header comment as in `test/test_math_conics.jl`; Earth GM `398600.4415`, the GMAT value of ANISE's tests)

```julia
@testset "Orbit: construction and Cartesian elements" begin
    s = SVector(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0)
    o = Orbit(s, 0.0, 399, 1, 398600.4415)
    @test o.state === s && o.center === Int32(399) && o.frame === Int32(1)
    @test (o.x, o.vz, o.rmag) == (s[1], s[6], norm(s[1:3]))
    @test o.position == s[1:3] && o.velocity == s[4:6]
    @test o.hvec == cross(s[1:3], s[4:6])
    @test o.hx ≈ 35065.806679607005 && o.hz ≈ 24796.2925419       # GMAT, val_state_def_circ_inc
    @test o.energy ≈ -25.842247282849137
    @test o.light_time == o.rmag / clight()
    @test element(o, :hz) === o.hz
    @test all(n -> n in propertynames(o), elements())
    @test throws_kind(() -> element(o, :foo), :badarg)
    @test throws_kind(() -> o.foo, :badarg)
    @test throws_kind(() -> element(o, :position), :badarg)       # vectors are properties only
    @test all(isnan, o.radii)
end

@testset "Orbit from kernels" begin
    ctx = neo_context([Kernels.LSK, Kernels.DE440S, Kernels.GM_TPC, Kernels.PCK11])
    et = 1.0e8
    o = orbit(ctx, "MOON", et, "J2000", "NONE", "EARTH")
    @test o.state == spkezr(ctx, "MOON", et, "J2000", "NONE", "EARTH")[1]
    @test o.gm == bodvcd(ctx, 399, "GM")[1] && o.radii == bodvcd(ctx, 399, "RADII")
    @test Orbit(ctx, o.state, et, "EARTH", "J2000") == o
    @test throws_kind(() -> Orbit(ctx, o.state, et, -999, "J2000"), :kernelvarnotfound)
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `julia --startup-file=no --project=test test/runtests.jl analysis_orbit`
Expected: FAIL, `UndefVarError: Orbit`.

- [ ] **Step 3: Implement** the interfaces above. `element` is one `if name === :x ... elseif ...` chain (generated by a macro over the name list is fine) ending in `SpiceError(:badarg, "element: unknown name :foo; the names are ...")`. `getproperty` checks `hasfield` first with `getfield`. Measure the kind `Orbit(ctx, ...)` throws for a body without GM and use it in the test (expected `:kernelvarnotfound`; adjust the test only to the measured kind).

- [ ] **Step 4: Run to verify they pass** — same command, expected PASS.

- [ ] **Step 5: Commit**

```bash
git add src/analysis test/test_analysis_orbit.jl src/Almanac.jl test/runtests.jl
git commit -m "Orbit type, constructors and Cartesian elements"
```

### Task 2: Keplerian elements and their degenerate cases

**Files:**
- Modify: `src/analysis/elements.jl`, `test/test_analysis_orbit.jl`
- Create: `test/gmat_cases.jl`

**Interfaces:**
- Consumes: `Orbit`, `_degenerate` (Task 1).
- Produces: `_el_` for `:sma :ecc :inc :raan :aop :ta :ma :ea :ha :aol :tlong :period :mean_motion :ta_dot :c3 :fpa :semi_parameter :semi_minor_axis :periapsis :apoapsis :right_ascension :declination :velocity_declination`; property `:evec`; `keplerian(o) -> NamedTuple{(:sma,:ecc,:inc,:raan,:aop,:ta)}`.

Port each from `astro/orbit.rs` (and `orbit_hyperbolic.rs` for `:ha`), converting degrees to radians at the end, not inside comparisons with ANISE's thresholds (`ECC_EPSILON = 1e-11`; `INC_EPSILON_DEG = 1e-6` used as `|ẑ×ĥ| < deg2rad(1e-6)` in `aol_deg`). Where ANISE returns `Err`, throw `:degeneratecase`; where it substitutes a value (`aop_deg`, `raan_deg` turn a NaN `acos` into 0 or 180°), return the same value in radians.

- [ ] **Step 1: Port the reference cases** into `test/gmat_cases.jl` as a vector of named tuples `(name, state, gm, expected::Dict{Symbol,Float64})`, from ANISE `tests/astro/orbit.rs`: `val_state_def_circ_inc`, `val_state_def_elliptical`, `val_state_def_circ_eq`, `val_state_def_equatorial`, `hyperbolic_inbound_true_anomaly_is_accepted`, `val_fpa_matches_geometric_definition`. Convert degree values to radians with `deg2rad` in the file; keep ANISE's digits.

- [ ] **Step 2: Write the failing tests**

```julia
include("gmat_cases.jl")
@testset "Keplerian elements: GMAT cases ($(c.name))" for c in GMAT_CASES
    o = Orbit(c.state, 0.0, 399, 1, c.gm)
    for (k, v) in c.expected
        @test isapprox(element(o, k), v; rtol=RTOL_GMAT, atol=ATOL_GMAT)   # measured, stated
    end
end

@testset "Keplerian elements against oscltx" begin
    rng = Xoshiro(SEED + 31)
    for _ in 1:nsample(200, 5000)
        o = random_orbit(rng)            # elliptic and hyperbolic, prograde and retrograde
        el = oscltx(o.state, o.et, o.gm)
        @test o.periapsis ≈ el[1] && o.ecc ≈ el[2] && o.inc ≈ el[3]
        @test o.raan ≈ el[4] && o.aop ≈ el[5] && o.ta ≈ el[9]
        o.ecc < 1 && @test o.ma ≈ el[6] && o.sma ≈ el[10] && o.period ≈ el[11]
    end
end

@testset "Keplerian elements: degenerate cases" begin
    gm = 398600.4415
    orb(v) = Orbit(SVector{6,Float64}(v), 0.0, 399, 1, gm)
    circ = orb((7000.0, 0, 0, 0, 0, sqrt(gm / 7000)))         # polar, ecc ≈ 0
    pro = orb((7000.0, 0, 0, 0, 8.0, 0))                       # equatorial, inc == 0
    retro = orb((7000.0, 0, 0, 0, -8.0, 0))                    # equatorial, inc == π
    para = orb((7000.0, 0, 0, 0, 0, sqrt(2gm / 7000)))         # |ecc - 1| < 1e-11
    @test circ.ecc < 1e-11 && throws_kind(() -> circ.ta, :degeneratecase)
    @test circ.aop in (0.0, π)                                 # ANISE's fallback for NaN acos
    for o in (pro, retro)
        @test o.raan in (0.0, π) && throws_kind(() -> o.aol, :degeneratecase)
    end
    @test retro.inc == π
    @test throws_kind(() -> para.sma, :degeneratecase) && throws_kind(() -> para.ma, :degeneratecase)
    for o in (circ, pro, retro, para), n in elements()        # Review Focus 2: no NaN escapes
        y = outcome(() -> element(o, n))
        @test y isa Symbol || !isnan(y)
    end
    @test keplerian(pro) == (sma=pro.sma, ecc=pro.ecc, inc=pro.inc, raan=pro.raan, aop=pro.aop, ta=pro.ta)
end
```

`random_orbit(rng)` is a helper in the group file, reused by later tasks: periapsis 6500–50000 km, ecc in [0, 0.95] ∪ [1.05, 5], angles uniform, `gm = 398600.4415`, built with `keplerian_orbit` once Task 5 exists and from `conics` until then. `outcome` is the existing helper of `test/mathoracle.jl` (include it as `test_math_conics.jl` does). Before writing each expectation of this testset, read the ANISE function for that case; if ANISE's behaviour differs from the expectation written here (for example `:raan` throws instead of falling back), follow ANISE and note it in the test.

- [ ] **Step 3: Run to verify they fail** — `julia --startup-file=no --project=test test/runtests.jl analysis_orbit`; expected FAIL on unknown names.
- [ ] **Step 4: Implement** the `_el_` functions and `keplerian`.
- [ ] **Step 5: Run to verify they pass**; record the largest differences with `measured(...)` and set `RTOL_GMAT`/`ATOL_GMAT` from them with a comment giving the reason (GMAT prints 15–17 digits).
- [ ] **Step 6: Commit** — `git commit -m "Keplerian elements of Orbit, with ANISE's degenerate cases"`

### Task 3: Hyperbolic, B-plane and v∞ quantities

**Files:** Modify `src/analysis/elements.jl`, `test/gmat_cases.jl`, `test/test_analysis_orbit.jl`

**Interfaces:**
- Produces: `_el_` for `:bplane_bt :bplane_br :bplane_mag :bplane_angle :bplane_ltof`; `vinf_periapsis(o::Orbit, turn_angle::Real) -> T` (km), `vinf_turn_angle(o::Orbit, periapsis::Real) -> T` (rad).

Port from `astro/orbit_hyperbolic.rs` (`BPlane::new`, `b_ltof`) and `astro/orbit.rs` (`vinf_periapsis_km`, `vinf_turn_angle_deg`, main's C3 form). B-plane of an elliptic orbit throws `:degeneratecase`.

- [ ] **Step 1: Add the cases** `b_plane_davis`, `val_b_plane_gmat`, `vinf_turn_angle_uses_excess_velocity` from ANISE `tests/astro/orbit.rs` to `GMAT_CASES` (B-plane values) and a separate `VINF_CASES`.
- [ ] **Step 2: Write the failing tests**: GMAT cases as in Task 2; identities on random hyperbolas: `bplane_mag ≈ hypot(bplane_bt, bplane_br)`, `bplane_angle ≈ atan(bplane_br, bplane_bt)`, `vinf_turn_angle(o, o.periapsis) ≈ 2asin(1/o.ecc)`, `vinf_periapsis(o, vinf_turn_angle(o, rp)) ≈ rp`; an inbound hyperbola (`r·v < 0`) gives `π < o.ta < 2π` and `o.ma < 0` (Review Focus 1); elliptic orbit `throws_kind(() -> o.bplane_bt, :degeneratecase)`.
- [ ] **Step 3: Run, fail. Step 4: implement. Step 5: run, pass; record measured differences.**
- [ ] **Step 6: Commit** — `git commit -m "B-plane and v-infinity quantities of Orbit"`

### Task 4: Geodetic, altitude and equinoctial elements; Brouwer mean elements

**Files:** Modify `src/analysis/elements.jl`, `test/test_analysis_orbit.jl`; create `src/analysis/brouwer.jl` (include after `elements.jl`)

**Interfaces:**
- Produces: `_el_` for `:altitude :periapsis_altitude :apoapsis_altitude :sma_altitude :latitude :longitude :height :equinoctial_a :equinoctial_h :equinoctial_k :equinoctial_p :equinoctial_q :equinoctial_lambda :brouwer_sma :brouwer_ecc :brouwer_inc :brouwer_raan :brouwer_aop :brouwer_ma`.

Mean equatorial radius `(radii[1] + radii[2]) / 2`, polar radius `radii[3]`. Geodetic: port `latlongalt` of `astro/orbit_geodetic.rs` (closed form). Equinoctial: `astro/orbit_equinoctial.rs`. Brouwer: `astro/orbit_mean_elements.rs`, constant `J2_EARTH = 1.082626925638815e-3`, validity check `is_brouwer_short_valid` → `:degeneratecase`. A NaN radius → `SpiceError(:badarg, "<name> needs the radii of body <id>, which the pool does not hold")` (Review Focus 3).

- [ ] **Step 1: Write the failing tests**
  - against `recgeo(o.position, re, f)` with `re = radii[1]`, `f = (radii[1] - radii[3]) / radii[1]` on random body-fixed positions with `radii[1] == radii[2]`: latitude, longitude, height agree (measured tolerance; closed form vs. iteration);
  - `verif_geodetic_vallado` from ANISE `tests/astro/orbit.rs` (added to `gmat_cases.jl`);
  - equinoctial identities: `h == ecc*sin(aop+raan)`, `k == ecc*cos(aop+raan)`, `p == tan(inc/2)*sin(raan)`, `q == tan(inc/2)*cos(raan)` (or the exact forms of the ANISE source if it differs — read first), `equinoctial_a == sma`;
  - Brouwer: ANISE's own `assert_brouwer_roundtrip` cases (`orbit_mean_elements.rs` tests) reproduced: mean → osculating → mean returns within ANISE's tolerances (4 km, 1e-2, 0.7°) and our measured tighter value; invalid orbit (periapsis < 3000 km) throws `:degeneratecase`;
  - orbit about body 3 (no radii in the pool): `o.sma` works, `throws_kind(() -> o.height, :badarg)`.
- [ ] **Step 2: Run, fail. Step 3: implement. Step 4: run, pass; record measured differences.**
- [ ] **Step 5: Commit** — `git commit -m "Geodetic, equinoctial and Brouwer mean elements of Orbit"`

### Task 5: Building and editing orbits; two-body propagation

**Files:** Create `src/analysis/build.jl`; modify `test/test_analysis_orbit.jl`

**Interfaces:**
- Consumes: Tasks 1–2.
- Produces:
  - `keplerian_orbit([ctx,] et, center, frame; sma=nothing, ecc=nothing, periapsis=nothing, apoapsis=nothing, inc, raan, aop, ta=nothing, ma=nothing, gm=nothing) -> Orbit{Float64}`. Exactly one of {`sma`+`ecc`, `periapsis`+`apoapsis`} and one of {`ta`, `ma`}, else `:badarg`. With `gm` given and integer `center`/`frame`, no kernels are read.
  - `with_element(o, name::Symbol, value) -> Orbit` for `:sma :ecc :inc :raan :aop :ta`; other names `:badarg`.
  - `with_apsides(o, periapsis, apoapsis) -> Orbit`.
  - `two_body(o, et2) -> Orbit` (port `at_epoch`, equinoctial mean-longitude advance).

Port from `astro/orbit.rs`: `try_keplerian`, `try_keplerian_apsis_radii`, `try_keplerian_mean_anomaly`, `with_*`/`set_*`, `at_epoch`.

- [ ] **Step 1: Write the failing tests**
  - `keplerian_orbit` against `conics([rp, ecc, inc, raan, aop, m0, et, gm], et)` on random elliptic and hyperbolic elements (state agrees, measured);
  - round trip `keplerian(keplerian_orbit(...; keplerian(o)...))` ≈ `keplerian(o)`;
  - `with_element(o, :sma, a).sma ≈ a` and the other five elements unchanged; same for each of the six names; `with_apsides` gives the asked apsides;
  - `two_body(o, o.et + dt).state` against `prop2b(o.gm, o.state, dt)` for random `dt` within ±10 periods (measured); `verif_orbit_at_epoch` of ANISE `tests/astro/orbit.rs`; circular equatorial orbit propagates without error;
  - argument errors: both `ta` and `ma`, neither, `sma` without `ecc` → `:badarg`.
- [ ] **Step 2: Run, fail. Step 3: implement. Step 4: run, pass.**
- [ ] **Step 5: Commit** — `git commit -m "Orbit constructors from elements, element editing, two-body propagation"`

### Task 6: Local frames and relative states; allocation proof

**Files:** Create `src/analysis/local_frames.jl`; modify `test/test_analysis_orbit.jl`, `test/test_alloccheck.jl`

**Interfaces:**
- Produces:
  - `local_frame(o, kind::Symbol) -> SMatrix{3,3,T}` (local → orbit frame), `kind ∈ (:ric, :vnc, :rcn, :tnw)`, else `:badarg`. Columns: RIC `(r̂, ĥ×r̂, ĥ)`; VNC `(v̂, ĥ, v̂×ĥ)`; TNW `(v̂, ĥ×v̂, ĥ)`; RCN = RIC.
  - `local_frame_state(o, kind) -> SMatrix{6,6,T}`: `[R 0; Ṙ R]` with `Ṙ` exact (the analytic derivative of the columns along the two-body motion: `ṙ = v`, `v̇ = -gm r/|r|³`).
  - `relative_state(o, chief, kind) -> SVector{6,T}`: `transpose(local_frame_state(chief, kind)) * (o.state - chief.state)`, as ANISE's `ric_difference_unchecked` (the transpose of the 6×6 matrix is its inverse, because `RᵀṘ` is skew); `:badarg` unless `o` and `chief` have the same `et`, `center` and `frame`.

- [ ] **Step 1: Write the failing tests**
  - each kind: `R' * R ≈ I`, `det(R) ≈ 1`, first column as listed;
  - `local_frame(o, :rcn) == local_frame(o, :ric)`;
  - `Ṙ` of `local_frame_state` against `ForwardDiff.derivative(t -> local_frame(two_body(o, t), kind), o.et)` (measured);
  - `relative_state(o, o, k) == zeros(6)`; a radial offset of 1 km gives `[1, 0, 0, …]` in RIC;
  - mismatched epoch, center, frame → `:badarg`; unknown kind → `:badarg`.
- [ ] **Step 2: Add to `cases` of `test/test_alloccheck.jl`**: `(element, (Orbit{Float64}, Symbol))`, `(keplerian, (Orbit{Float64},))`, `(local_frame, (Orbit{Float64}, Symbol))`, `(local_frame_state, (Orbit{Float64}, Symbol))`, `(relative_state, (Orbit{Float64}, Orbit{Float64}, Symbol))`, `(two_body, (Orbit{Float64}, Float64))`, `(with_element, (Orbit{Float64}, Symbol, Float64))`.
- [ ] **Step 3: Run, fail. Step 4: implement; fix any allocation AllocCheck reports.**
- [ ] **Step 5: Run** `julia --startup-file=no --project=test test/runtests.jl analysis_orbit alloccheck` — PASS; note the time of `analysis_orbit` (budget: state it in the group header; target under 5 s).
- [ ] **Step 6: Commit** — `git commit -m "Local orbit frames and relative states; allocation proof of the Orbit calls"`

### Task 7: Manual

**Files:** Create `docs/src/howto/orbits.md`, `docs/src/reference/orbits.md`; modify `docs/make.jl`, `docs/src/explanation/coverage.md`; docstrings of Tasks 1–6 as needed.

Follow the `almanac-docs` skill.
- `docs/make.jl`: new part `"Beyond CSPICE" => ["howto/orbits.md", "reference/orbits.md"]` after `"Explanation"`.
- `reference/orbits.md`: `@docs` blocks for all 14 exports, and a table of `elements()` with units and the ANISE name (generated in an `@example` from `elements()` so it cannot drift).
- `howto/orbits.md`: one task from start to end on `docs/fixture.jl` (`PROBE` around `EARTH`): build the orbit, read elements, change the inclination, propagate, RIC difference of two orbits.
- `coverage.md`: section "Beyond CSPICE" with its own validation table (GMAT/Vallado/Davis values, CSPICE equivalents, ANISE group, derivatives) and the statement that the CSPICE accuracy claims above do not cover it.

- [ ] **Step 1: Write the pages; run each example in a scratch script first.**
- [ ] **Step 2: Build** `julia --project=docs docs/make.jl` — PASS, no `checkdocs` missing-export error; read the rendered outputs in `docs/build/` once.
- [ ] **Step 3: Reader review** with the `julia-doc-reader-review` skill (agent `julia-doc-reader`) on the two pages and the reference text from `docs/reference_text.jl`; apply findings.
- [ ] **Step 4: Commit** — `git commit -m "Manual: Orbit how-to and reference, Beyond CSPICE coverage"`

### Task 8: The opt-in `anise` group

**Files:** Create `test/anise/Project.toml`, `test/anise/CondaPkg.toml`, `test/anise/runtests.jl`

**Interfaces:**
- Consumes: all exports of Tasks 1–6.

`Project.toml`: `Almanac` (via `[sources] Almanac = {path = "../.."}`), `PythonCall`, `CondaPkg`, `StaticArrays`, `Random`, `Test`. `CondaPkg.toml`: `[pip.deps] anise = "==1.0.0b0"` now, with a comment that the reference is main at `b3aaadc` and the pin moves to the first release containing it.

`runtests.jl`:
- Build ANISE frames with the same GM and radii as Almanac (`Frame` with `mu_km3_s2` and an `Ellipsoid`), so no ANISE data file is needed; use the stubs `astro.pyi` of the installed package for the exact names.
- Sample `nsample`-like counts (default 2000, `ALMANAC_TEST_FULL=1` 20000) per class: elliptic, near-circular (`ecc < 1e-6`), near-equatorial, retrograde, hyperbolic outbound and inbound; around Earth (399) and Mars (499).
- For each element of `elements()`, each local frame kind and `relative_state`, compare with ANISE's counterpart (degrees converted), collect relative differences, print the largest per quantity with the state that produced it, and `@test` against a per-quantity tolerance set from the first measured run, with the reason in a comment.
- Constant `KNOWN_DIFFERENT_IN_PIN`: quantities where 1.0.0b0 differs from main by definition (`:fpa`, `:ma`/`:ha`/`:ta` of hyperbolas, inbound hyperbolas, TNW, `vinf_*`). These are reported but not tested, with a printed notice; the set is emptied when the pin moves.
- Pinned regression values: a printed Julia literal of 20 states with ANISE's results, for pasting into `test/test_analysis_orbit.jl` once the pin is a release containing main's definitions.

- [ ] **Step 1: Write the files.**
- [ ] **Step 2: Run** `julia --startup-file=no --project=test/anise test/anise/runtests.jl` — the first run installs Python and anise through CondaPkg; expected PASS for every quantity outside `KNOWN_DIFFERENT_IN_PIN`. A difference outside the set is a bug in Almanac's port or a definition change between 1.0.0b0 and main: read both sources before changing anything.
- [ ] **Step 3: Commit** — `git commit -m "Opt-in comparison of the Orbit calls with ANISE"`

### Task 9: Derivatives through six AD systems (after the merge)

**Gate:** starts only after the AD rule extensions are committed on `main` (owned by session almanac-1d; the user decides when) and `main` is merged into `analysis`. Before that, skip this task.

**Files:** Create `test/test_analysis_derivatives.jl`; modify `test/Project.toml` (add DifferentiationInterface, Enzyme, Mooncake, Zygote, ReverseDiff, Tracker), `test/runtests.jl` (`"analysis_derivatives"` after `"analysis_orbit"`).

- [ ] **Step 1: Write the tests**: for backends `AutoForwardDiff(), AutoEnzyme(), AutoMooncake(), AutoZygote(), AutoReverseDiff(), AutoTracker()`:
  - gradient with respect to the state of `s -> element(Orbit(s, 0.0, 399, 1, gm), q)` for each `q` in a representative subset (`:sma :ecc :inc :raan :aop :ta :ma :latitude :bplane_bt :brouwer_sma`) and of `s -> sum(local_frame(Orbit(s,…), :ric))`, `s -> two_body(Orbit(s,…), 3600.0).state[1]`;
  - derivative with respect to the epoch of `t -> element(orbit(ctx, "MOON", t, "J2000", "NONE", "EARTH"), :sma)`;
  - each against ForwardDiff and a central finite difference (step `1e-6 * max(1, |x|)`), tolerance measured.
  A backend that fails to load or to differentiate on the running Julia version: `skip_notice(...)` with the error, not a silent pass.
- [ ] **Step 2: Run** `julia --startup-file=no --project=test test/runtests.jl analysis_derivatives`; measure the time. If the group exceeds 15 s, keep ForwardDiff and one reverse backend (Mooncake) in the default run and move the rest under `FULL`, stating the times in the group header.
- [ ] **Step 3: Commit** — `git commit -m "Derivatives of the Orbit calls through six AD systems"`

### Task 10: Whole-branch check

- [ ] Run `julia --startup-file=no --threads=8 test/parallel.jl` — all groups PASS, total under the 120 s budget (single process) as the skill states.
- [ ] Run `julia --project=docs docs/make.jl` — PASS.
- [ ] Commit nothing new unless a fix was needed.
