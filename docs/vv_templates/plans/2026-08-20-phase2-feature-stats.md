# Phase 2 Feature-Statistics Batch — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** V&V five feature-statistics filters — `ComputeBiasedFeatures`, `ComputeSurfaceAreaToVolume`, `ComputeEuclideanDistMap` (SimplnxCore) and `ComputeShapes`, `ComputeSchmids` (OrientationAnalysis) — against hand-derived Class 1 (+ Class 4) oracles plus full binary A/B vs DREAM3D 6.5.171, fully retiring `6_6_find_biased_features` and removing four of the six remaining live consumers of `6_6_stats_test_v2`.

**Architecture:** Every filter currently compares its output against arrays that live *inside* the same 6.6-produced `.dream3d` it reads as input — certifying "NX reproduces 6.6," not "NX is correct." Each task replaces that with small hand-built fixtures whose expected output is derived from the algorithm contract before any execution, then runs the same fixtures through both binaries as migration evidence. Research (2026-08-20, five parallel source reads) pre-identified a defect slate (§Pre-identified findings) that the fixtures are designed to catch rather than mask.

**Tech Stack:** C++20, Catch2, CMake, DREAM3D 6.5.171 `PipelineRunner`, `nxrunner`, h5py (conda env `dream3d`), EbsdLib 3.1.0, Eigen.

---

## Why these five

From `docs/vv_templates/vv_program_plan.md` §3.3, these are five of the eight remaining opening-batch dual-purpose filters — the whole Feature-statistics group across both plugins. Archive leverage:

| Filter | Plugin | Algorithm LOC | Archive it pins | Archive fate |
|---|---|---|---|---|
| `ComputeBiasedFeatures` | SimplnxCore | 343 | `6_6_find_biased_features` | **fully retires** (sole consumer) |
| `ComputeSurfaceAreaToVolume` | SimplnxCore | 168 | `6_6_stats_test_v2` | consumption removed; archive stays |
| `ComputeEuclideanDistMap` | SimplnxCore | 473 | `6_6_stats_test_v2` | consumption removed; archive stays |
| `ComputeShapes` | OrientationAnalysis | 690 | `6_6_stats_test_v2` | consumption removed; archive stays |
| `ComputeSchmids` | OrientationAnalysis | 113 (+EbsdLib) | `6_6_stats_test_v2` | consumption removed; archive stays |

`6_6_stats_test_v2` does **not** retire: `AlignSectionsMutualInformation` (out of scope, Phase 3) and `ComputeFeatureNeighbors` (`ComputeFeatureNeighborsTest.cpp:859`) remain consumers after this batch. Its `download_test_data()` lines in BOTH plugins' test CMakeLists stay untouched.

Style models on `develop`: `RequireMinNumNeighborsTest.cpp` (`DiscriminatingFixture`), `ErodeDilateBadDataTest.cpp` (inline legacy-verified constants + provenance tripwires), `vv/ComputeKernelAvgMisorientationsFilter.md` and `vv/ComputeFeatureNeighborsFilter.md` (report style; the latter's D1 is the anisotropic-defect-proved-by-patched-legacy precedent Task 2 follows). If the Phase-2 cleanup PRs (#1713–#1717) have merged by execution time, their tests are additional style models.

---

## Standing requester decisions (carried from the cleanup batch)

- **OOC build runs: waived.** Record the dated waiver in each report's outstanding list.
- **Second-engineer sign-off: delegated to the PR reviewer.** Same wording as the cleanup batch reports, header row AND body.
- **One PR per filter**, branch `vv/<FilterName>` off `upstream/develop`, PR titles styled `VV: <Human Name> Fully V&V'ed`.
- **Preflight validation guards are wanted** where legacy had them or where NX silently accepts values that produce garbage/UB (the erode/dilate family precedent, 2026-08-19). Messages include the offending value; unique negative error codes per file's series; failing test first.
- **A/B evidence** goes in `/Users/mjackson/Workspace9/ww_work/<FilterName>/` (never committed) and is uploaded at the end via `rclone copy /Users/mjackson/Workspace9/ww_work/<FilterName> vv_work:<FilterName> --progress` (exclude `__pycache__`).

## Global Constraints

- **"Legacy produced this output" is never a valid oracle for correctness.** Oracle first; the A/B is diff-explanation and migration evidence, run only after SIMPLNX matches the oracle.
- **Derive expected values by hand (or by an independent script for Class 2 cases) BEFORE running.** Never paste observed output as "expected."
- **Mutation verification is a required step, not a review nicety:** for every discriminating fixture, temporarily apply the target mutation, prove ONLY the new test fails (record the failing index), revert, prove the diff is empty. The cleanup batch's five-lens reviews showed "executed" ≠ "discriminated."
- **Absolute claims need absolute proof.** If a behaviour is data-dependent (termination, tie outcomes, truncation flips), the docs must say so with the evidence class labeled (executed vs source-derived) — the CN wording episode is the cautionary tale.
- Every commit uses the `VV:` prefix (per `docs/vv_templates/commit_template.md`); `Signed-off-by: Michael Jackson <mike.jackson@bluequartz.net>` last; **never** any AI attribution; never commit `CLAUDE.md`, `.claude/`, or anything under `ww_work`.
- C++20, Allman, 2-space, 200-col; `k_CamelCase` constants, `m_CamelCase` members, local pointers `camelBack`+`Ptr`, `Geom`/`Ref` suffixes.
- Every `TEST_CASE` ends with `UnitTest::CheckArraysInheritTupleDims(dataStructure);`; all `getDataRefAs<T>()` in tests wrapped in `REQUIRE_NOTHROW()`; `CAPTURE` per index on element-wise loops.
- Run tests only via `ctest -R "<pattern>" --verbose`; preserve every SIMPL `BackwardsCompatibility` `TEST_CASE` untouched.
- clang-format touched `.cpp`/`.hpp`: `/opt/local/clang+llvm-15.0.4-arm64-apple-darwin21.0/bin/clang-format -i -style=file <files>` from repo root.
- Float fixtures use exactly-representable values (dyadic spacings/origins) wherever an exact comparison is asserted.

**Repo:** `/Users/mjackson/Workspace9/simplnx` (planning branch `vv/forward_plan`; task branches off `upstream/develop`).
**Build:** `/Users/mjackson/Workspace9/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel` — target `SimplnxCoreUnitTest` or `OrientationAnalysisUnitTest` per plugin; `nxrunner` at `Bin/nxrunner`.
**Legacy binary:** `/Users/mjackson/Applications/DREAM3D.app/Contents/Bin/PipelineRunner` (6.5.171).
**Legacy sources:** `/Users/mjackson/Workspace/D3D_v6.5.171/DREAM3D/Source/...` (per-task paths below; ignore `DREAM3D-Build/` dirs).
**6.5.172 surgical-patch branch (legacy-shared bugs only):** `/Users/mjackson/Desktop/ShellScripts/init_6_5_172_build.sh 9` → `/Users/mjackson/Workspace9/6.5.172/DREAM3D`, `cmake --preset D3D-Rel-Qt515-6_5_171`.
**A/B tooling to reuse:** `/Users/mjackson/Workspace9/ww_work/ErodeDilateMask/make_fixtures.py` (proven legacy-format writer), `run_ab.sh`/`compare.py` patterns in the three cleanup-batch `ww_work` folders; venv `/Users/mjackson/Workspace9/ww_work/.vvenv` or conda env `dream3d` (h5py fixed).

## THE ONE RULE THAT MATTERS

> **Derive the expected values by hand from the algorithm source. Then run the test.
> If they differ, adjudicate (oracle error vs filter bug) and write down which.
> Never run the filter first and paste its output as "expected."**

## Legacy A/B protocol

Identical to the cleanup batch (`2026-08-18-phase2-cleanup-trio.md` §Legacy A/B protocol) with these batch-specific rules:

1. The unit-test fixtures ARE the A/B fixtures. Write them as legacy-format `.dream3d` with the proven writer; author legacy pipeline JSON with SIMPL **property names**; NX `.d3dpipeline` mirrors.
2. **Expected divergences are pre-identified** (§Pre-identified findings). Plan each A/B run to land on the right side of them; an A/B difference matching a predicted divergence confirms the analysis, an unpredicted difference enters the Bug adjudication protocol.
3. **Legacy-shared output bugs (requester decision 2026-08-20):** document the deviation **against 6.5.171** (the A/B binary is always `~/Applications/DREAM3D.app`'s 6.5.171 `PipelineRunner`); fix NX in-PR; **then** create the surgical patch on the 6.5.172 branch whose sole purpose is *alignment validation* — the patched legacy build must reproduce the fixed NX output on the same fixtures, proving the fix brings simplnx and SIMPL into output alignment. Smallest possible diff, one filter per patch; applies to every confirmed shared bug in this batch (SA-1/SA-2, and BF-2 / SH-2 / SH-4 if they confirm).
4. Float outputs: state the tolerance and its cause per array. Int outputs derived from float truncation (`Poles`) get an explicit off-by-one budget on real data, exact assertion on hand fixtures.
5. Every A/B folder gets `ReadMe.md` + results files; upload via rclone at batch end.

## Bug adjudication protocol

As the cleanup batch, plus one scope clause:

- **Verdict A** — oracle error: fix the expectation, note what you misread.
- **Verdict B** — SIMPLNX bug: fix in-PR; read legacy source (paths per task) to classify; legacy-shared → Deviation entry + 6.5.172 surgical patch; NX-only port regression → report Summary, no deviation entry.
- **Verdict C (new) — EbsdLib bug. Requester decision 2026-08-20: FIX them.** Defects inside EbsdLib are fixed on the locally checked-out `topic/3_1_1_staging` branch of `/Users/mjackson/Workspace9/EbsdLib` (verified clean, tip `fcaa1de`), each as its own commit in that repo (github-workflow conventions: `BUG:` prefix, sign-off, no AI attribution; failing-then-passing evidence recorded in the report since EbsdLib's own test coverage may not exist for these). simplnx work that depends on the fix builds against the LOCAL EbsdLib via CMake preset **`NX-Com-Qt69-Vtk96-Rel-EbsdLib`** (DREAM3DNX `CMakeUserPresets.json:278`, binaryDir `/Users/mjackson/Workspace9/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib` — configure from `/Users/mjackson/Workspace9/DREAM3DNX` if absent). Oracles then assert the EXACT post-fix values. **Merge dependency (accepted):** the affected simplnx PR (Schmids) cannot pass stock CI until EbsdLib 3.1.1 is released and simplnx's `vcpkg.json` pin is bumped — note this prominently in that PR's body and carry the `vcpkg.json` `"version>="` bump in the branch. NX-side mitigations (zero-initializing output arrays, per-iteration local reinit) remain Verdict B fixes in the simplnx repo.

## Pre-identified findings (verify each — designed to be caught, not assumed)

From the 2026-08-20 source reads. Each task's fixtures must discriminate its rows. "Shared" = present in both 6.5.171 and NX.

| ID | Filter | Finding | Class → disposition |
|---|---|---|---|
| BF-1 | BiasedFeatures | 2D path uses `spacing[0]/[1]` for the box regardless of which axis is flat (`Algorithms/ComputeBiasedFeatures.cpp:263`); legacy remaps per axis | **NX port regression** → fix, Summary |
| BF-2 | BiasedFeatures | 2D classification compares raw centroid X/Y against the axis-shifted box (`:318-341`) — wrong for X-/Y-normal slabs | **Shared bug** → fix + Deviation + 6.5.172 patch |
| BF-3 | BiasedFeatures | Legacy 2D hardcodes origin to 0 (`FindBoundingBoxFeatures.cpp:331-333`); NX uses the real origin | Legacy-only bug (NX correct) → Deviation, Trust SIMPLNX |
| SA-1 | SurfaceAreaToVolume | ±X/±Y face areas swapped (`Algorithms/ComputeSurfaceAreaToVolume.cpp:128-139`: ±Y gets `dy·dz`, ±X gets `dz·dx`); masked by isotropic spacing | **Shared bug** → fix + Deviation + 6.5.172 patch |
| SA-2 | SurfaceAreaToVolume | Sphericity exponents `0.333333`/`0.66666` instead of 1/3, 2/3 (`:155-165`) → ~1e-4 relative error | **Shared bug** → fix + Deviation (+ ride the SA-1 patch) |
| SA-3 | SurfaceAreaToVolume | User doc claims id-0 neighbors don't count and never mentions the outer-boundary face skip — both false vs `:99-126` | NX doc defect → fix docs |
| SA-4 | SurfaceAreaToVolume | No preflight tuple/geometry cross-check FeatureIds↔ImageGeom → OOB read | NX gap → guard (standing decision) |
| ED-1 | EuclideanDistMap | Float "Euclidean" output is BFS-seed distance under tie-break +z>+y>+x>−x>−y>−z, NOT the true EDT; doc implies nearest-boundary | Shared semantics + NX doc defect → document, fix doc |
| ED-2 | EuclideanDistMap | Bad-data voxels in float mode: legacy `0.0`, NX `-1` (`ComputeEuclideanDistMap.cpp:301` vs legacy init 0) | Legacy/NX divergence → Deviation, recommend Trust SIMPLNX (`-1` consistent with Manhattan mode) |
| ED-3 | EuclideanDistMap | `nearestNeighbors` sized from FeatureIds tuples (`:277`) but phase 2 iterates geometry cells (`:58`); mismatched selections → OOB | NX gap → guard (standing decision) |
| ED-4 | EuclideanDistMap | `m_ShouldCancel` stored, never read; no progress messages | NX defect → add cancel check (CN precedent); progress optional |
| ED-5 | EuclideanDistMap | Legacy accepts all-toggles-false silently; NX errors `-12802` (whose fmt message drops its argument — no `{}` placeholder, `Filter.cpp:134-137`) | Deviation (guard is better) + fix the message |
| SH-1 | Shapes | 2D degenerate Euler branch writes `pi/180` (`Algorithms/ComputeShapes.cpp:655`); legacy writes `pi/2` — 90× mis-transcription, corroborated by legacy's own unit test | **NX port regression** → fix, Summary |
| SH-2 | Shapes | 3D voxel-center fix (#1124) never applied to the 2D branch (`:423-424`) — half-voxel moment bias | **Shared bug** → fix + Deviation + 6.5.172 patch (or documented-accepted; adjudicate with evidence) |
| SH-3 | Shapes | Legacy 3D uses voxel corners; NX uses centers → 3D A/B **must not match**; legacy is wrong | Deviation, Trust SIMPLNX; prove via patched legacy if cheap |
| SH-4 | Shapes | 2D in-plane axis mapping always uses X/Y components regardless of flat axis (`:401-436`) | **Shared bug** → fix + Deviation + patch (bundle with SH-2) |
| SH-5 | Shapes | Handedness correction is dead code — `isValid()` never returns 0 (`OrientationMatrix.hpp:166-192`) | Shared latent → document (left unfixed: value-changing) |
| SH-6 | Shapes | 2D path never writes `omega3s`; NX preflight fill value must actually zero it — verify | Verify; fix fill if garbage |
| SC-1 | Schmids | Preflight creates outputs with empty fill → uninitialized memory; skip path (`laueClass >= LaueGroupEnd`) leaves garbage; legacy zero-initialized | **NX port regression** → fix (explicit zero fill), Summary |
| SC-2 | Schmids | EbsdLib `1.732f`/`1.414f` constants → Schmid factor can exceed 0.5 (measured 0.500090176 at the maximizing direction); uniform +0.018% bias | **Verdict C (EbsdLib)** → document + upstream recommendation; oracle asserts actual values with bias quantified |
| SC-3 | Schmids | `angleComps` units flip with `OverrideSystem` (cosines vs radians); docs call both "angles" | Shared + doc defect → fix docs, Deviation-style note |
| SC-4 | Schmids | `crystalStructures[featurePhases[i]]` unbounded (`Algorithms/ComputeSchmids.cpp:82`) → OOB for bad phase ids | Shared gap → NX guard (standing decision) |
| SC-5 | Schmids | EbsdLib stub Laue classes leak stale `angleComps`; `HexagonalLowOps` never initializes `schmidfactor`/`slipsys` | Verdict C → document; NX mitigation = per-iteration reinit of the locals (Verdict B, in scope) |

Do not take any row on faith: each task re-verifies its rows against source before designing the fixture that catches it.

---

## File Structure

Per filter `<F>` in its plugin `<P>` ∈ {SimplnxCore, OrientationAnalysis}:

| File | Responsibility |
|---|---|
| `src/Plugins/<P>/test/<F>Test.cpp` (Shapes: `ComputeShapesFilterTest.cpp`) | Replace exemplar tests with inline oracles; keep BC TEST_CASE |
| `src/Plugins/<P>/src/<P>/Filters/Algorithms/<F>.cpp` | Fixes per the findings table |
| `src/Plugins/<P>/src/<P>/Filters/<F>Filter.cpp` | Preflight guards / message fixes |
| `src/Plugins/<P>/docs/<F>Filter.md` | Doc-defect fixes + new constraints |
| `src/Plugins/<P>/vv/<F>Filter.md`, `vv/deviations/<F>Filter.md` | Report + deviations (`python scripts/vv_init.py <F>Filter`) |
| `src/Plugins/SimplnxCore/test/CMakeLists.txt` | Task 1 only: drop the `6_6_find_biased_features` block |
| `/Users/mjackson/Workspace9/ww_work/<F>/` | A/B artifacts + ReadMe (never committed) |
| `/Users/mjackson/Workspace9/6.5.172/DREAM3D` | Surgical patches for shared bugs (SA-1/SA-2 mandatory; BF-2, SH-2/SH-4 per adjudication) |

Reports for Shapes/Schmids follow the OrientationAnalysis `vv/` house style (`ComputeKernelAvgMisorientationsFilter.md`).

---

## Task 0: Sync and setup

- [ ] **Step 1:** `git fetch upstream` and check whether #1713–#1717 merged; note the state in the ledger. Task branches come off current `upstream/develop` either way.
- [ ] **Step 2:** Build both test targets once on a fresh `upstream/develop` checkout branch and run the baseline: `ctest -R "ComputeBiasedFeatures|ComputeSurfaceAreaToVolume|ComputeEuclideanDistMap|ComputeShapes|ComputeSchmids"` — record the pre-change pass state (expect all green).
- [ ] **Step 3:** Verify `PipelineRunner`, `nxrunner`, the legacy writer script, and `conda run -n dream3d python -c "import h5py"` all work. Create the five `ww_work/<F>/` folders.
- [ ] **Step 4:** EbsdLib readiness (controller-verified 2026-08-20: `/Users/mjackson/Workspace9/EbsdLib` on `topic/3_1_1_staging` @ fcaa1de, clean; preset `NX-Com-Qt69-Vtk96-Rel-EbsdLib` defined in DREAM3DNX `CMakeUserPresets.json:278`). Configure that preset before Task 5 if its build dir is absent.

## Task 1: ComputeBiasedFeatures

**Branch `vv/ComputeBiasedFeatures`. Legacy:** `Plugins/Generic/GenericFilters/FindBoundingBoxFeatures.cpp` (NOT Statistics; NOT named FindBiasedFeatures).

### Oracle derivation (verify against source, then freeze in comments)

Greedy per-phase box shrink: full geometry bounds; per surface feature (index order matters — document it) pull the single nearest face onto the centroid; classify biased = centroid `<=`/`>=` any box face. Feature 0 skipped. 2D drops phases entirely. Full trace rules in the 2026-08-20 research; the fixtures:

- **A (3D baseline, CalcByPhase=false):** 4×4×4, origin 0, spacing 1, box `[0,4]³`. f1 `(0.5,2,2)` surf, f2 `(2,2,2)` non-surf, f3 `(3.8,2,2)` surf, f4 `(2,0.5,2)` non-surf. Trace: f1 pulls min-x→0.5; f3 pulls max-x→3.8. Expected `{F,T,F,T,F}` (index 0 first).
- **B (no-shrink + order dependence):** a surface feature exactly on the min-x face (`move=0`, still biased); reorder A's f1/f3 to pin the greedy order.
- **C (origin+anisotropy 3D):** dims (2,3,4), origin (10,20,30), spacing (2,0.5,4) → box `[10,14]×[20,21.5]×[30,46]`.
- **D (CalcByPhase=true):** two phases; a phase-2 feature positioned to NOT shrink phase 1's box; a phase-0 feature stays false; `phases[0]` set large to pin (and prove benign) the max-element-includes-0 quirk.
- **E (2D Z-normal):** 5×5×1 — the self-consistent 2D case.
- **F (2D X-normal, anisotropic):** dims (1,5,6), spacing (1,2,3). Geometrically correct box (Y,Z) = `[0,10]×[0,18]`. Current NX computes `[0,5]×[0,12]` (BF-1) and both implementations classify unshifted components (BF-2). Fixture asserts the geometrically correct answer **after** the fixes; pre-fix it is the RED evidence for both findings. Add the Y-normal sibling.

**Class 4:** `Biased[0]==false`; every surface feature biased within its phase (3D + 2D-Z; this invariant is the automatic BF-2 detector on X/Y-normal); box monotonically shrinks.

- [ ] **Step 1:** Write fixtures A–F + invariants as failing tests (hand-derived comments first; exact dyadic floats).
- [ ] **Step 2:** Run — expect F (and Y-sibling) to fail against current code; A–E per derivation. Adjudicate any surprise.
- [ ] **Step 3:** Fix BF-1 (adopt legacy's per-axis spacing remap) and BF-2 (classify the shifted components). BF-2 is shared → full Deviation + 6.5.172 surgical patch (smallest diff in `FindBoundingBoxFeatures.cpp`'s 2D classify loop; patched build must reproduce fixture F's correct answer). Re-run green; mutation-verify: comparison-flip, no-shrink, phase-gate-drop, feature-0-unskipped each killed by a named fixture.
- [ ] **Step 4:** Legacy A/B: fixtures A–F × relevant toggles through both binaries. Predictions: A–E match; F diverges three ways (BF-1 NX-only pre-fix, BF-3 legacy-only origin bug, BF-2 shared). Legacy needs `SurfaceFeatures` as H5 **bool** (not uint8) — the writer must emit legacy bool layout; Phases only when CalcByPhase.
- [ ] **Step 5:** Deliverables: report (Class 1+4; BF-1 in Summary; BF-2 + BF-3 as Deviations; patch cited), deviations file, doc touch-ups. Remove the `6_6_find_biased_features` `download_test_data()` block (`SimplnxCore/test/CMakeLists.txt:227`) — grep guard first (research confirms sole consumer).
- [ ] **Step 6:** Commit `VV: Compute Biased Features fully V&V'ed` (body: findings, legacy-shared status, archive retirement). Waivers recorded.

## Task 2: ComputeSurfaceAreaToVolume

**Branch `vv/ComputeSurfaceAreaToVolume`. Legacy:** `Plugins/Statistics/StatisticsFilters/FindSurfaceAreaToVolume.cpp`.

### Oracle derivation

Surface = Σ interior faces between differing FeatureIds (id 0 **counts** as differing; outer-boundary faces **never** count — the trap). Volume = `NumCells[i]·dx·dy·dz` (user-supplied NumCells, not recounted). Sphericity `π^(1/3)(6V)^(2/3)/A` — fix SA-2 first so the oracle can assert to 1e-6. Fixtures (research-validated):

| Fixture | Expected (post SA-1/SA-2 fix) | Kills |
|---|---|---|
| F1 single interior voxel (3×3×3, spacing a) | A=6a², SAVR=6/a, Ψ=(π/6)^⅓≈0.80600 | baseline, sphericity constant |
| F2 2×2×2 interior cube (4×4×4) | SAVR=3/a, Ψ≈0.80600 (scale invariance) | volume/area confusion |
| F3 **1×4 rod along X, spacing (1,2,4)** in 6×3×3 | correct A=64 (current code gives 104) + Y-rod mirror (correct 64, code 40) | **SA-1 swap** — the RED fixture |
| F4 corner voxel touching boundary | A=3a², Ψ≈1.612 (>1 — the boundary-rule signature, comment it loudly) | boundary-face rule |
| F5 feature fills volume | A=0, SAVR=0, Ψ=+inf | degenerate divide (document) |
| F6 id-0 ring around center voxel | A=6a² — id 0 IS surface | falsifies the doc claim (SA-3) |
| F7 2D slab | A=4a² (no ±Z faces) | 2D path |
| F8 sphericity off | no Sphericity array created | untested toggle |
| F9 index-0 semantics | `SAVR[0]` never written | loop start |
| F10 error paths | `-5355` negative ids, `-5351` range, `-12802` | validation |

- [ ] **Step 1:** Write the failing tests (F3 asserts the CORRECT 64/40 — it is the RED evidence for SA-1). Fixtures keep `max(id)==numFeatures-1` so they run through legacy (its `-5555` exact-count guard).
- [ ] **Step 2:** Run; F3 must fail pre-fix. Fix **SA-1** (swap the two face-area terms) and **SA-2** (`1.0f/3.0f`, `2.0f/3.0f`); both shared → Deviations + one 6.5.172 surgical patch covering both; patched legacy must reproduce F3's 64/40 and F1's Ψ to 1e-6. Add the **SA-4 guard** (preflight FeatureIds-tuple-count vs geometry-cell-count, failing test first, unique error code). Mutation-verify: face-term swap-back, boundary-skip removal, id-0-skip addition, exponent re-truncation.
- [ ] **Step 3:** Legacy A/B: isotropic fixture set (predict: identical pre-fix behavior on F1/F2/F4–F7 — shared code; post-fix NX differs from stock legacy exactly on F3 + sphericity 4th decimal; patched legacy matches NX). Sphericity-off cannot be A/B'd (legacy never reads the flag — L2; say so). Sparse-id case: legacy `-5555`, NX accepts → Deviation.
- [ ] **Step 4:** Deliverables: report; Deviations for SA-1, SA-2 (legacy-shared, patch cited), the L1 legacy-only sphericity-off crash (NX fixed — Trust SIMPLNX), L4 strictness difference; **rewrite the user doc** (SA-3: boundary-face rule stated, id-0 claim removed, the copy-pasted SurfaceFeatures paragraph deleted). Also fix the missing `REQUIRE_NOTHROW`s the old test had. `6_6_stats_test_v2` CMake lines untouched.
- [ ] **Step 5:** Commit `VV: Compute Surface Area to Volume fully V&V'ed`. Waivers recorded.

## Task 3: ComputeEuclideanDistMap

**Branch `vv/ComputeEuclideanDistMap`. Legacy:** `Plugins/Statistics/StatisticsFilters/FindEuclideanDistMap.cpp`.

### Oracle derivation

Manhattan (int32) = layer-synchronous 6-connected BFS graph distance — deterministic, order-independent, fully hand-derivable. Float mode = straight-line distance **to the BFS-selected seed** under tie-break +z>+y>+x>−x>−y>−z (ED-1: not a true EDT — the Class 1 float oracle must reimplement the propagation + tie-break; a scipy EDT will disagree by design). Seeds: GB ≥1 differing face neighbor (id 0 counts, `>= 0`); TJ ≥2 distinct; QP >2. Fixtures:

- **A — the legacy 10×6×1 fixture** (`FindEuclideanDistMapTest.cpp:120-127`, spacing (1,2,1)) — adopt it: it pins seeds, bad-data column, anisotropy, tie-break (`GBEuclidean[0]=4.0` where true EDT is 2.0 — assert 4.0 with the ED-1 comment), and is independently corroborated by legacy's own published arrays. Re-derive several rows by hand anyway (Verdict A protection).
- **B** 6×1×1 two-feature split: GB `{2,1,0,0,1,2}` both modes.
- **C** 3×3×1 triple junction (TJ seeds). **D** 3×3×1 pinwheel QP (quad points do NOT need 3D — state it) + a 3×3×3 variant with non-unit `spacing[2]` for the z-decode.
- **E** single feature: everything stays −1. **F** all id 0: everything −1 (the ED-2 pin — legacy float gives 0.0 here).
- Toggle matrix per fixture: 3 single-toggle + 1 all-three, × both distance modes; the all-three run doubles as the `ParallelTaskAlgorithm` determinism check (run 3×).
- **Class 4:** `{QP=0} ⊆ {TJ=0} ⊆ {GB=0}`; `-1` exactly where `featureId<=0` (both modes for NX); adjacent-in-feature Manhattan deltas ≤1; `euclidean >= trueEDT` pointwise.

- [ ] **Step 1:** Failing tests per above, plus a small Python reference generator for the float mode checked into `ww_work` (Class 1 by independent reimplementation; document the tie-break in it).
- [ ] **Step 2:** Run; adjudicate. Fix **ED-4** (cancel check at the sweep level, CN precedent) and the **ED-5 message** (`{}` placeholder). Add the **ED-3 guard** (FeatureIds tuples == geometry cells, failing test first). Mutation-verify: seed-threshold shifts, `>=0`→`>0` bad-data admission, fill(−1)→fill(0), tie-break reversal (kills via A's `GBEuclidean[0]`), distance-commit moved inside the sweep.
- [ ] **Step 3:** Legacy A/B: all fixtures × both modes × toggle sets. Predictions: Manhattan identical everywhere; float identical except bad-data voxels (ED-2: legacy 0.0 vs NX −1) → Deviation, Trust SIMPLNX; all-toggles-false: legacy silent no-op vs NX `-12802` → Deviation (guard better). Legacy pipeline sets `SaveNearestNeighbors=false` (one extra run with `true` captured as tie-break evidence). Property-name JSON keys (the `readFilterParameters` key mismatch is dead code — verified).
- [ ] **Step 4:** Deliverables: report; Deviations ED-2, ED-5, plus the ED-1 semantics note **and a user-doc fix** (the "distance to nearest boundary cell" implication); document the reciprocal-decode latent hazard as shared, with recommended integer division; note the concurrent-DataStore-read hazard for OOC (waived runs, but record it). `6_6_stats_test_v2` lines untouched.
- [ ] **Step 5:** Commit `VV: Compute Euclidean Distance Map fully V&V'ed`. Waivers recorded.

## Task 4: ComputeShapes

**Branch `vv/ComputeShapes`. Legacy:** `Plugins/Statistics/StatisticsFilters/FindShapes.cpp` (+ its `FindShapesTest.cpp` — a ready-made validated fixture source).

### Oracle derivation (validated closed form)

8-point octant quadrature per voxel; for an axis-aligned box of `N_x×N_y×N_z` voxels, spacing `d`: `P_i = V_tot·d_i²·[(N_i²−1)/12 + 1/16]` (the `/16` self-moment is the implementation's quadrature, NOT the exact `/12` integral — a deliberate discriminator); `Ixx=P_y+P_z` etc.; `Q_i = 15P_i/(4π)`; `m=(Q_xQ_yQ_z)^{1/5}`, `a=√(Q_x/m)`, `b=√(Q_y/m)`, `c=√(Q_z/m)`; `Ω3 = V⁵/(P_xP_yP_z·2000π²/9)`. Validated against legacy's own test numbers (256×128×64 box → (121.0, 40.3, 10.1) vs asserted (120,40,10)±1.5; Ω3 0.789 vs 0.78715; 2D case likewise). Fixtures:

- **F1** small anisotropic box (8×4×2, spacing (0.75,0.5,0.25)) — the full 3D chain, closed-form.
- **F2** two boxes + one empty feature id — per-feature accumulation + the NaN-on-empty-feature hole. **Requester decision 2026-08-20: FIX it** — zero the axis outputs (lengths, ratios, Eulers, Ω3) for features with no voxels, failing test first; consistent with legacy's zero-initialized arrays. Verdict B, Summary.
- **F3** single voxel; **F4** isotropic cube (a=b=c, ratios 1.0, TripletSort ties).
- **F5** Z-flat 2D box (port legacy's) — 2D chain + **SH-1 directly**: degenerate branch expects π/2; current NX gives π/180. RED evidence.
- **F6** X-flat and Y-flat 2D boxes — SH-4 (+SH-2): assert the geometrically correct values post-fix.
- **F7** 45°-about-Z rotated box — eigen solve; expected via a small numpy reference (Class 2, script + versions recorded); assert axis DIRECTIONS up to sign (|dot| ≈ 1), **never element-wise Euler angles** (sign/gimbal arbitrariness — SH-5).
- **F8** invariants: a≥b≥c; ratios ∈(0,1]; Ω3∈[0,1]; Σvolumes = V_geom on a fully-labelled grid; no NaN/Inf (excluding the documented F2 hole until fixed).

- [ ] **Step 1:** Failing tests (F5 pins SH-1 red; F6 pins SH-2/SH-4 red). Verify SH-6 (2D omega3s fill) while building F5.
- [ ] **Step 2:** Fix **SH-1** (π/2 — port regression, Summary). Adjudicate **SH-2/SH-4** with the binary: both shared → fix both + Deviations + one 6.5.172 surgical patch for the 2D branch (voxel-center + axis remap), patched legacy reproducing F5/F6. Fix the wrong `[SimplnxCore]` test tag. Mutation-verify: quadrature /16→/12, konst swaps, TripletSort inversion, 15/(4π) chain breaks.
- [ ] **Step 3:** Legacy A/B: F1–F6 + Small-IN100-scale optional. Predictions: **3D must NOT match** (SH-3, legacy corner bias — Deviation, Trust SIMPLNX; quantify on F1 where the closed form arbitrates); `Volumes` bit-exact everywhere (the sanity anchor); 2D matches stock legacy only pre-fix; Euler angles excluded from element-wise diff (compare axes up to sign); residual float tolerances ~1e-5 with causes named.
- [ ] **Step 4:** Deliverables: report; Deviations SH-2/SH-3/SH-4 (+SH-5 documented latent, D9 verification result); doc updates. Archive lines untouched.
- [ ] **Step 5:** Commit `VV: Compute Feature Shapes fully V&V'ed`. Waivers recorded.

## Task 5: ComputeSchmids

**Branch `vv/ComputeSchmids`. Legacy:** `Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindSchmids.cpp`. EbsdLib 3.1.0 at `/Users/mjackson/Workspace9/EbsdLib` (conventions verified identical to legacy: quat layout (x,y,z,w), qu2om line-identical, no transpose; the only numeric delta is float→double).

### Step 0 for this task: fix EbsdLib on `topic/3_1_1_staging` (Verdict C, per the requester decision)

In `/Users/mjackson/Workspace9/EbsdLib`: (a) replace the `1.732f`/`1.414f` normalizers in both `CubicOps::getSchmidFactorAndSS` overloads (and any sibling ops using the same literals — grep) with exact `std::sqrt(3.0)`/`std::sqrt(2.0)`-derived constants (SC-2); (b) initialize `schmidfactor`/`slipsys` in `HexagonalLowOps` and make the stub Laue classes write defined `angleComps` (SC-5). One `BUG:` commit per defect, sign-off, no AI attribution, do NOT push without requester confirmation. Then configure/build the `NX-Com-Qt69-Vtk96-Rel-EbsdLib` preset and use ITS build dir for every build/test/`nxrunner` in this task.

### Oracle derivation (exact values, post-EbsdLib-fix)

Identity quaternion, `Cubic_High`, auto-select path. With SC-2 fixed the oracle asserts the EXACT-arithmetic values (left column); the pre-fix EbsdLib values (right) are recorded in the deviations file as the ≤3.1.0 bias (+0.0180%) for migration readers:

| Loading | exact m (ASSERT post-fix) | pre-fix EbsdLib m (deviations record) | SlipSystems | Poles |
|---|---|---|---|---|
| [0,0,1] | 1/√6 = 0.408248290 | 0.408321919 | 0 | (0,0,100) |
| [1,1,1] | 2/(3√6) = 0.272165527 | 0.272214613 | 4 | (57,57,57) |
| [0,1,1] | 1/√6 = 0.408248290 | 0.408321919 | 1 | (0,70,70) |
| [1,2,3] | 0.466569475 | 0.466653622 | 10 | (26,53,80) |

Phis/Lambdas: re-derive the exact cosines post-fix (the research table's values carry the pre-fix normalizer bias — recompute, don't reuse). Poles values must be re-verified post-fix too (truncation boundaries can flip).

Re-derive at least [0,0,1] and [1,2,3] independently before freezing (the executor recomputes; the table is the cross-check, not the source). [1,2,3] has a UNIQUE max (no ties) — the enumeration-order pin; [0,1,1] separates Phis from Lambdas; ties resolve to lowest index via strict `>` (assert SlipSystems==0 on [0,0,1]). Add: a 90°-about-Z quaternion case (equivariance: must equal the rotated-loading answer — Class 4); an `OverrideSystem=true` case (other overload: sym-op-index numbering, RADIAN angleComps — assert and document the unit flip SC-3); `StoreAngleComponents=false`; a `laueClass >= LaueGroupEnd` feature (post-SC-1 fix: outputs zero, not garbage); scale-invariance (L vs 3L). Class 4 (post-fix): the physical bound `0 ≤ m ≤ 0.5` is now assertable (pre-fix EbsdLib violated it at 0.500090176 — record in deviations); `Poles[k]==trunc(100·(om·L̂)[k])`; `m == Phis·Lambdas` (auto path only).

- [ ] **Step 1:** Failing tests per above (hand-derivation comments first).
- [ ] **Step 2:** Fix **SC-1** (explicit zero fill on all five created arrays — failing test = the LaueGroupEnd-skip case asserting zeros). Add **SC-4 guard** (phase id bounds vs ensemble count, failing test first). NX-side **SC-5 mitigation**: reinit `schmid/slipSystem/angleComps` per loop iteration (kills the stale-leak; EbsdLib defect still documented as Verdict C). Mutation-verify: OM transpose, quat component swap, Poles round-vs-trunc, tie-break `>`→`>=`.
- [ ] **Step 3:** Legacy A/B: the 4 loading directions + override case + equivariance quat on one hand-built input, vs stock 6.5.171. Predictions: **floats now differ by the +0.018% normalizer bias** (legacy carries the old constants) — a predicted, documented deviation (Trust SIMPLNX, fixed in EbsdLib 3.1.1); `SlipSystems` unchanged (bias is uniform, argmax stable); `Poles` may flip at trunc boundaries — budget it; legacy pipeline uses `SchmidPhis`/`SchmidLambdas` names and error `-1001` vs NX `-13500`. Alignment validation: rebuild the 6.5.172 branch against... NOT applicable here — legacy's Schmid math lives in OrientationLib, so the alignment proof is the patched-EbsdLib NX values matching the exact-arithmetic oracle; state that explicitly instead of patching legacy.
- [ ] **Step 4:** Deliverables: report (cite Dieter §4-3 / Schmid & Boas for the 0.408 anchor — Class 3 garnish on the Class 1 oracle); Deviations: SC-2 + SC-5 as documented EbsdLib defects with upstream recommendation (include the HexagonalLowOps uninitialized locals), SC-3 doc/unit-flip note, `Poles`-is-not-a-Miller-index doc fix, the `(05)`→`(04)` pipeline cite fix. Archive lines untouched.
- [ ] **Step 5:** Commit `VV: Compute Schmids fully V&V'ed`. Waivers recorded.

## Task 6: Hit rate + program plan

- [ ] Full `ctest` on both plugins' suites at each branch (esp. Task 1's archive removal): no new failures.
- [ ] On `vv/forward_plan`: append this batch's hit-rate line to `vv_program_plan.md` §3.3 (Verdict B count/of 5, legacy-shared count → deviations + 6.5.172 patches shipped, Verdict C EbsdLib count, port regressions) and update the `6_6_stats_test_v2` consumer count (6→2 live) + `6_6_find_biased_features` retired. Commit `VV: Record the feature-statistics batch defect hit rate`.

## Task 7: Archive + PRs

- [ ] `rclone copy` each of the five `ww_work/<F>` folders to `vv_work:<F>` (exclude `__pycache__`); verify with `rclone size`.
- [ ] Push branches; open five PRs titled `VV: <Human Name> Fully V&V'ed`, Summary/Test-Plan bodies, no AI attribution; note in each stats-archive consumer's PR that the archive intentionally stays.

## Definition of done

- [ ] Five test suites assert inline Class 1 (+Class 4) oracles; zero consumption of `6_6_stats_test_v2` by these five tests; `6_6_find_biased_features` fully retired; `grep -rn "6_6_find_biased_features" src/` empty (doc self-references excepted).
- [ ] Every pre-identified finding adjudicated with evidence (fixed / deviation / Verdict C documented) — none silently dropped; every fix has a test that fails without it; every discriminating fixture mutation-verified with the blind-suite proof.
- [ ] Legacy-shared fixes carry Deviation entries AND 6.5.172 surgical patches (SA-1/SA-2 mandatory; BF-2, SH-2/SH-4 per adjudication), patched-legacy reproduction recorded.
- [ ] Full binary A/B per filter with per-array tolerance stories; predicted divergences confirmed and documented; artifacts in ww_work with ReadMes, uploaded via rclone.
- [ ] Reports/deviations per templates; `| Status | READY FOR REVIEW |` unbolded; OOC waiver + PR-reviewer sign-off delegation recorded (header + body); BC TEST_CASEs untouched; every commit `VV:`-prefixed with sign-off, no AI attribution.

## Requester decisions (2026-08-20, all questions resolved)

1. **Shared bugs:** deviations documented against 6.5.171; A/B always against `~/Applications/DREAM3D.app` (6.5.171); the 6.5.172 surgical patch is the *alignment-validation* step for every confirmed shared bug (patched legacy must reproduce fixed-NX output).
2. **Preflight guards** (SA-4, ED-3, SC-4): in scope per the 2026-08-19 family-wide decision.
3. **Shapes F2 NaN:** fix (zero empty features).
4. **EbsdLib:** fix SC-2/SC-5 on the local `topic/3_1_1_staging` branch; build via `NX-Com-Qt69-Vtk96-Rel-EbsdLib`; Schmids PR carries the vcpkg bump and is merge-blocked on the EbsdLib 3.1.1 release (noted in its PR body). Do not push EbsdLib without requester confirmation.
5. **Execute now**, subagent-driven, one PR per filter, sequential tasks.
