# V&V Report: ComputeAvgOrientationsFilter

|           |                          |
|-----------|--------------------------|
| Plugin    | OrientationAnalysis      |
| SIMPLNX UUID               | `086ddb9a-928f-46ab-bad6-b1498270d71e`          |
| SIMPLNX Human Name         | Compute Feature Average Orientations            |
| DREAM3D 6.5.171 equivalent | `FindAvgOrientations` (SIMPL UUID `bf7036d8-25bd-540e-b6de-3a5ab0e42c5f`) — `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindAvgOrientations.{h,cpp}` |
| Verified commit            | `a307946e7` (v7.4.2 release)         |
| Status | COMPLETE     |
| Sign-off  | Michael Jackson <mike.jackson@bluequartz.net> — 2026-07-16 |
| Second-engineer sign-off | Michael Jackson (technical authority) — 2026-07-16 |

## At a glance

| Aspect                 | Current state            |
|------------------------|--------------------------|
| Algorithm Relationship | **Minor changes** (Rodrigues method) **+ New** (vMF/Watson). The Rodrigues average is a faithful port of legacy `FindAvgOrientations::execute()` with 4 deliberate deltas (positive-orientation canonicalization, FeatureId==0 inclusion, divide-by-zero cleanup, EbsdLib library swap). The von Mises-Fisher and Watson methods are **new to DREAM3D-NX** (added PR #1577), have **no legacy equivalent**, are off by default, and delegate their EM math to EbsdLib `DirectionalStats` — itself a C++ port of the EMsoft `dictmod.f90` directional-statistics module (Marc De Graef Research Group, Carnegie Mellon University; original EM implementation by Yu-Hui Chen). The math is published and peer-reviewed [R1] [R2], not novel to this project. |
| Oracle (confirmed)     | **Confirmed.** Rodrigues → **Class 1 (analytical)** + **Class 4 (invariant)**; vMF/Watson → **Class 2 (EbsdLib reference, trusted, not re-tested)** + **Class 4 (invariant)**, scoped to the filter's value-add per the "don't re-test upstream" rule. Encoded as 2 oracle TEST_CASEs (+1 error test) in `test/ComputeAvgOrientationsTest.cpp`; all pass. SIMPLNX matches the oracle on every fixture. |
| Code paths enumerated  | 15 of 16 paths exercised (see Code path coverage); the one gap is the Watson-only dispatch branch (low-value array-selection path).              |
| Tests today            | 9 test cases: Rodrigues analytical (Class 1+4), Rodrigues cubic-symmetry invariant (Class 4, #1660), Rodrigues voxel-ordering independence (Class 4), vMF/Watson EbsdLib-reference (Class 2+4), vMF/Watson phase-0 exclusion regression (#1659), crystal-structure/phase guard warnings (#1661), no-method-enabled error (preflight -54673 + runtime backstop -54670), cell-array tuple-mismatch error (-651, code asserted), and SIMPL 6.4/6.5 backwards-compat (DYNAMIC_SECTION). All inline hand-built data — no exemplar archive. Closes the GCOV vMF/Watson coverage gap.       |
| Exemplar archive       | **None — retired `7_ComputeAvgOrientation_v2.tar.gz`** (it was a circular oracle regenerated from post-fix SIMPLNX output in PR #1577). Oracle is now inline in the test source. `download_test_data()` entry to be removed in Phase 10.               |
| Legacy comparison      | **Run (2026-06-30) vs the official DREAM3D 6.5.171 release.** Two fixtures. On the realistic 480k-cell / 408-feature input the two are **numerically equivalent on all real features** (`AvgQuats` ≤1e-6, zero sign flips; `AvgEulerAngles` max 4.77e-7). Divergences are confined to **feature-0 / unindexed-voxel handling**: **D3** (empty feature 0) and **D2** (FeatureId-0 voxels with phase>0 — demonstrated on a forcing fixture: NX computes the average, legacy writes `(0,0,0,0)`). **D4** sub-epsilon. **D1 downgraded** (no demonstrable divergence). vMF/Watson have no legacy equivalent (N/A). |
| Bug flags              | One **legacy** bug, empirically confirmed: `ComputeAvgOrientationsFilter-D3` (6.5.171 leaves the zero-voxel feature 0 as `(0,0,0,0)` / divides zero-count features by zero; SIMPLNX writes a clean identity). D2 is a deliberate algorithmic-choice divergence (not a bug). One **SIMPLNX** bug found by the post-merge adversarial review and fixed 2026-07-08: **issue #1659** — the vMF/Watson gather loop ignored `Phases`, so a phase-0 voxel inside a feature contributed a garbage quaternion to the EM average. Fixed by gating the gather identically to the counting pass; pinned by the `vMF/Watson Ignores Phase-0 Voxels` regression test. |
| V&V phase | **COMPLETE.** |

## Summary

`ComputeAvgOrientationsFilter` ("Compute Feature Average Orientations") computes a per-Feature average crystallographic orientation from per-Element (voxel) quaternions, using one or more of three independently-toggled methods: the original **Rodrigues** running-quaternion average (a port of legacy `FindAvgOrientations`), and two **new** Expectation-Maximization distribution fits — **von Mises-Fisher** and **Watson** — that delegate their sampling math to EbsdLib `DirectionalStats`, a C++ port of the EMsoft `dictmod.f90` module written by Dr. Marc De Graef's group at Carnegie Mellon University [R1] [R2] [R3]. Verification establishes correctness independently of legacy: Rodrigues via a Class 1 analytical oracle (hand-built single- and multi-voxel features with closed-form averages) plus Class 4 invariants; vMF/Watson via Class 2 (EbsdLib as trusted reference, not re-tested) plus Class 4 invariants scoped to the filter's value-add. Headline result: SIMPLNX matches the oracle on every fixture and is numerically equivalent to the official DREAM3D 6.5.171 on all real features; the only legacy divergences are the deliberate feature-0 / unindexed-voxel handling (D2, D3) and a sub-epsilon precision difference (D4). All 9 unit tests pass.

## Algorithm Relationship

*Classification:* **Minor changes** (Rodrigues method, the legacy-equivalent output path) **+ New filter capability** (vMF/Watson methods, no legacy equivalent).

*Evidence:* Same SIMPL UUID retained (`bf7036d8-…` → SIMPLNX `086ddb9a-…`); SIMPL 6.4/6.5 conversion fixtures at `test/simpl_conversion/{6_4,6_5}/ComputeAvgOrientationsFilter.json` convert only the six legacy parameters (FeatureIds, Phases, Quats, CrystalStructures, AvgQuats, AvgEulerAngles) — i.e. the Rodrigues method. The legacy `FindAvgOrientations::execute()` (DREAM3D 6.5.171, 347 lines) computes **only** the Rodrigues average; `computeRodriguesAverage()` (`ComputeAvgOrientations.cpp` lines 367–457) preserves its control flow line-for-line (per-voxel running-average accumulation with `getNearestQuat`, then per-feature divide+normalize). The vMF/Watson methods were added in **PR #1577** (`parametersVersion` bumped 1→2) and are off by default.

### vMF / Watson provenance — EMsoft port chain

The von Mises-Fisher and Watson methods have no DREAM3D 6.5.171 counterpart, so there is nothing to A/B against. Their verification instead rests on their **provenance**: they are a two-stage port of an established, published, third-party reference implementation, not new mathematics authored for this project.

| Stage | Artifact | Author / owner | Evidence |
|---|---|---|---|
| 1. Original Fortran | EMsoft `Source/EMsoftLib/dictmod.f90`, module `dictmod` — `DI_EMforDD`, `DD_Estep`, `DD_Mstep`, `DD_getQandL`, `DD_Density`, `logCp`, `DI_SampleDD`, `VMFMeanDirDensity`, `WatsonMeanDirDensity` plus the Bessel and PRNG helpers | Marc De Graef Research Group / Carnegie Mellon University; EM routine attributed in-source to "Yu-Hui Chen, U. Michigan / Marc De Graef, Carnegie Mellon University". BSD 3-Clause, Copyright (c) 2014-2022. | `dictmod.f90` header block and the `DI_EMforDD` doc comment, which cite [R1] and [R2] as the derivation. |
| 2. C++ port | EbsdLib `ebsdlib::DirectionalStats` (`EbsdLib/Core/DirectionalStats.{hpp,cpp}`) — `EMforDS`, `Estep_`, `Mstep_`, `getQandL_`, `Density_`, `logCp_` | Michael Jackson, BlueQuartz Software | EbsdLib commit `006a023` (2026-04-03), PR #41 "ENH: Implement Watson and von-Mises Fisher Orientation Averaging". Port fidelity is visible in the source: in-line comments read "mirroring the Fortran routines", "matches the Fortran `r8_uniform_01` exactly", "faithful port of the Fortran `r8vec_normal_01`", and "Fortran 1-based -> 0-based". Both files carry a `THIRD-PARTY ATTRIBUTION: EMsoft` header block giving the full routine correspondence, the paper citations, and the retained EMsoft BSD-3 notice (see *Attribution / license compliance* below). |
| 3. Filter integration | `ComputeAvgOrientations::computeVmfWatsonAverage()` / `VmfWatsonSamplingImpl` | BlueQuartz Software, PR #1577 | Per-feature gather → `getFZQuat` → `DirectionalStats::EMforDS(43514, …)` → `positiveOrientation` → write quat / Euler / kappa at the feature tuple. |

Two consequences for this V&V:

1. **The Class 2 oracle is well-founded.** Treating EbsdLib as a trusted reference is not merely "trust an upstream library" — the library reproduces a published algorithm [R1] [R2] whose Fortran original is publicly available, and EbsdLib's own `DirectionalStatsTest` hard-asserts the result at this filter's exact configuration. See *The EbsdLib boundary* under Oracle.
2. **The statistical model is the *modified* (symmetry-group-invariant) vMF / axial Watson distribution**, not the textbook form. Both `dictmod.f90` and `DirectionalStats::Estep_` build the density as an equal-weight mixture over the `Pmdims = LaueOps::getNumSymOps()` proper rotation operators of the phase's Laue class. This is the formulation of [R2]. The filter docs (`docs/ComputeAvgOrientationsFilter.md`) state this explicitly.

**Attribution / license compliance.** `EbsdLib/Core/DirectionalStats.{hpp,cpp}` each carry a two-part header block: the standard BlueQuartz Software BSD-3 notice for the port, and a `THIRD-PARTY ATTRIBUTION: EMsoft` block naming the Marc De Graef Research Group at Carnegie Mellon University, crediting Yu-Hui Chen (University of Michigan) with the original EM implementation as the EMsoft source does, giving the routine-by-routine correspondence (`EMforDS` ← `DI_EMforDD`, `Estep_` ← `DD_Estep`, `Mstep_` ← `DD_Mstep`, `getQandL_` ← `DD_getQandL`, `Density_` ← `DD_Density`, `logCp_` ← `logCp`, plus the `math.f90` Bessel and PRNG helpers), citing [R1] and [R2], linking the EMsoft repository, and reproducing the EMsoft BSD 3-Clause copyright notice, conditions, and disclaimer verbatim as clause 1 of that license requires. The provenance chain in the table above is therefore evidenced by the source files themselves, not only by this report.

### Rodrigues port-time deltas (each a candidate Deviation — see `vv/deviations/ComputeAvgOrientationsFilter.md`)

1. **Positive-orientation canonicalization** (D1 — **downgraded, not a deviation**). SIMPLNX appends `.getPositiveOrientation()` after normalize (line 448), forcing the average quaternion into the northern hemisphere; legacy `FindAvgOrientations` calls only `QuaternionMathF::UnitQuaternion` (normalize, line 270). Empirically (2026-06-30) this produces **no divergence** from 6.5.171 — symmetry reduction keeps the average at `w ≥ 0`, and the pure-triclinic edge case is the same physical orientation (`q ≡ −q`). Defensive normalization only. See deviations file.
2. **FeatureId == 0 inclusion** (candidate D2). Legacy gates accumulation on `m_FeatureIds[i] > 0 && m_CellPhases[i] > 0` (line 246). SIMPLNX gates on `currentPhase > 0` only (line 412) — deliberately allowing voxels with FeatureId 0 to be averaged, per the documented use-case of averaging an unlabeled bag of orientations (algorithm comment, lines 401–411). Root cause: algorithmic choice. **Observable only** when a dataset has FeatureId 0 with phase > 0 (atypical).
3. **Zero-count finalize cleanup** (candidate D3 — legacy bug). Legacy finalize loop (lines 263–274) starts at feature 1 and, for `counts[i] == 0`, sets identity but then still executes `ScalarDivide(avgQuats[i], counts[i])` — a divide by zero. SIMPLNX (lines 440–444) sets identity and `continue`s, skipping the divide. Root cause: bug in 6.5.171. **Observable** for any zero-voxel feature; legacy also never finalizes feature 0 (loop starts at 1), leaving it `(0,0,0,0)`.
4. **Library swap** (candidate D4 — library + precision). Legacy uses `QuaternionMathF` + `OrientationTransforms::qu2eu`; SIMPLNX uses EbsdLib `ebsdlib::QuatF` + `QuaternionFType::toEuler()`. Same `float32` precision. **Possible** sub-epsilon Euler differences from differing intermediate math/quaternion→Euler conventions.

*Material PRs since the filter's creation:* **#1577** (added vMF/Watson + version 2 + regenerated the `_v2` exemplar — the circular-oracle source; commit note also says "Fix phaseIndex calculation bug"), **#1588** (SIMPL backwards-compat test redesign), **#1535** (removed redundant preflight checks now done by parameters), **#1476** (backwards pipeline compatibility fix), **#1472** (EbsdLib 2.0.0 API — cross-cutting EbsdLib-drift hotspot), **#1458** (EbsdLib 1.0.40), **#1438** (microtexture-related cleanup — cross-cutting hotspot).

### Algorithm review (V&V cycle quality pass)

Line-by-line review performed via the `review-algorithm` skill on the already-oracle-verified code. Findings applied (verified output-neutral — all tests still pass after rebuild):

- **Thread safety:** `computeVmfWatsonAverage` now calls `dataAlg.setParallelizationEnabled(false)`. The worker read/wrote shared `DataArray`/`DataStore` objects concurrently, which the simplnx thread-safety policy forbids even at distinct indices. Serial execution matches the `ComputeFeatureFaceMisorientation` precedent.
- **Cancel checking:** added `m_Filter->getCancel()` (new accessor) to the vMF/Watson per-feature loop, which previously could not be cancelled during EM.
- **Progress messaging:** `computeRodriguesAverage` now emits throttled progress via `MessageHelper`/`ThrottledMessenger` over the per-cell loop.
- **Robustness:** guarded `ops[laueClass]` / `orientationOps[xtal]` against out-of-range `CrystalStructures` values (e.g. `999`/Unknown) to avoid undefined behavior; the feature/voxel is skipped (vMF/Watson outputs stay NaN).
- **Cosmetic:** removed commented-out dead code, fixed the `dataStruture` parameter typo, removed the unused `m_LastProgressInt` member.

**Deferred (documented rationale):** the vMF/Watson voxel gathering is O(numFeatures × numVoxels) — each feature rescans all voxels rather than bucketing once. This is a *performance* optimization, not correctness; with parallelization now disabled it runs serially but remains acceptable for typical feature counts. Left as a follow-up to avoid restructuring verified code in this V&V pass.

## Oracle

*Status: confirmed. Boundary with EbsdLib confirmed; second-engineer oracle review signed off by Michael Jackson (technical authority) 2026-07-16.*

*Class:* **1 (Analytical) + 4 (Invariant)** for Rodrigues; **2 (Reference — EbsdLib, trusted & not re-tested) + 4 (Invariant)** for vMF/Watson.

### The EbsdLib boundary (what we do NOT re-test)

The vMF/Watson EM math lives in EbsdLib `DirectionalStats::EMforDS` — a C++ port of EMsoft's `DI_EMforDD` (`dictmod.f90`), implementing the modified vMF / axial Watson estimators of [R1] [R2] (see *vMF / Watson provenance* above) — and is already covered by `EbsdLib/Source/Test/DirectionalStatsTest.cpp`. Critically, `DirectionalStatsTest:VMF` and `:Watson` run **at this filter's exact configuration** — `numEM=5, numIter=10, seed=43514` — on a fixed set of 22 hand-coded quaternions, and hard-assert the resulting `muhat` (±1e-6) and `kappahat`:
- VMF: `muhat = (w 0.88937, x 0.33220, y −0.19646, z 0.24507)`, `kappahat = 88.9943`
- WAT: `muhat = (w 0.90119, x 0.29483, y −0.21060, z 0.23787)`, `kappahat = 30.5730`

EbsdLib owns: EM convergence/multi-restart (`numEM` = EMsoft `Num_of_init` random restarts, best-likelihood wins; `numIter` = EMsoft `Num_of_iterations` EM iterations per restart, with early exit when the Q-function changes by <0.01), symmetry-aware E/M steps, kappa estimation, FZ mapping of the result, seeded-PRNG reproducibility (the port reproduces EMsoft's Park-Miller `r8_uniform_01` and `r8vec_normal_01` bit-for-bit in double precision). We trust those (Class 2 reference), on the strength of both the published derivation [R1] [R2] and EbsdLib's own hard-asserted unit test. We do **not** reproduce the EM algorithm, and we do **not** re-derive it from the papers.

### Applied

- **Rodrigues (Class 1 — analytical):** Hand-built features with closed-form averages — a single-voxel feature averages to that voxel's FZ-reduced, northern-hemisphere quaternion; N identical-orientation voxels average to that same orientation; a symmetric pair about a known axis averages to the bisector. Expected quats/Eulers derived by hand.
- **Rodrigues (Class 4 — invariant):** unit-norm + northern-hemisphere (`w ≥ 0`) output quats; zero-voxel feature → identity `(0,0,0,1)` + zero Euler; result invariant to voxel ordering within a feature.
- **vMF/Watson (Class 2 — EbsdLib reference, reused not reinvented):** Construct a single SIMPLNX feature whose voxels are exactly EbsdLib's 22 reference quaternions (same single Laue class), run the filter, and assert the filter's per-feature `muhat`/`kappa` match EbsdLib's **already-asserted** VMF/WAT values above. This proves the filter feeds EbsdLib correctly and lands the result in the right tuple — without re-deriving any EM math. (If the 22 quats are not pre-FZ-reduced, account for the filter's `getFZQuat` step in the input or expected value.)
- **vMF/Watson (Class 4 — the filter's value-add):** single-voxel feature → `muhat` equals that voxel's FZ quat and `kappa == 0` (filter shortcut, algorithm lines 103–106 / 139–141, EM skipped); zero-voxel feature → all outputs NaN (fill at lines 339–355); correct per-feature voxel gathering + FZ reduction + phase→crystal-structure lookup; results written at the feature index (not voxel index); unit + northern-hemisphere quats.

*Encoded:*
- `test/ComputeAvgOrientationsTest.cpp::"OrientationAnalysis::ComputeAvgOrientations: Rodrigues Analytical Oracle"` — Class 1 (6 Rodrigues fixtures, exact quats at 1e-6) + Class 4 (unit/northern, zero→identity, vMF/Watson single→FZ-quat & kappa 0, zero→NaN).
- `test/ComputeAvgOrientationsTest.cpp::"OrientationAnalysis::ComputeAvgOrientations: vMF/Watson EbsdLib Reference Oracle"` — Class 2 (22-quat feature reproduces EbsdLib's documented VMF/WAT muhat+kappa; quat margin 5e-3, kappa 2% to absorb the float32 input round-trip) + Class 4 (zero-voxel→NaN).
- `test/ComputeAvgOrientationsTest.cpp::"OrientationAnalysis::ComputeAvgOrientations: No Method Enabled Error"` — error path `-54670`.

- `test/ComputeAvgOrientationsTest.cpp::"OrientationAnalysis::ComputeAvgOrientations: Rodrigues Cubic Symmetry Invariant"` — Class 4 under non-trivial (Cubic_High) symmetry: five representations of the same orientation (Rz(30°+90°k) and −Rz(30°)) must average to Rz(30°) (issue #1660).
- `test/ComputeAvgOrientationsTest.cpp::"OrientationAnalysis::ComputeAvgOrientations: Rodrigues Voxel Ordering Independence"` — Class 4: F3 fixture in both voxel orderings → identical Rz(45°).
- `test/ComputeAvgOrientationsTest.cpp::"OrientationAnalysis::ComputeAvgOrientations: vMF/Watson Ignores Phase-0 Voxels"` — regression pin for issue #1659.
- `test/ComputeAvgOrientationsTest.cpp::"OrientationAnalysis::ComputeAvgOrientations: Unknown Crystal Structure and Out-Of-Range Phase Guards"` — guard/warning coverage for issue #1661 (-54671/-54672).

All 9 test cases pass in the `NX-Com-Qt69-Vtk95-Rel` build.

*Second-engineer review:* **Signed off by Michael Jackson (technical authority), 2026-07-16.**

## Code path coverage

*15 of 16 paths exercised. Source: `src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/Algorithms/ComputeAvgOrientations.cpp`.* Logical phases: **(a)** `operator()` dispatch + output-array selection/validation, **(b)** Rodrigues two-pass average, **(c)** vMF/Watson per-feature EM (serial over features). Paths 12–16 were added by the #1659/#1661 follow-up work (2026-07-08).

| #  | Phase            | Path       | Test case            |
|----|------------------|----|----------------------|
| 1  | (a) Dispatch     | `avgEulerAnglesArrayPath` exists → use it to set `m_NumberOfFeatures` (Rodrigues enabled)       | `Rodrigues Analytical Oracle` (Rodrigues on) |
| 2  | (a) Dispatch     | else `VMFEulerAnglesArrayPath` exists → use it (vMF enabled, Rodrigues off)    | `vMF/Watson EbsdLib Reference Oracle` (Rodrigues off, vMF on) |
| 3  | (a) Dispatch     | else `WatsonEulerAnglesArrayPath` exists → use it (Watson enabled, Rodrigues & vMF off)         | *Not directly tested. Low-value array-selection branch; only sets `m_NumberOfFeatures` source, equivalent to paths 1/2.* |
| 4  | (a) Dispatch     | none exist → error `-54670` "A valid Feature level array … was not found"      | `No Method Enabled Error` |
| 5  | (b) Rodrigues    | `currentPhase > 0` → accumulate running nearest-quat average (incl. FeatureId 0)                | `Rodrigues Analytical Oracle` — F1–F4 |
| 6  | (b) Rodrigues    | `counts[fid] == 1` → reset running average to identity before getNearestQuat   | `Rodrigues Analytical Oracle` — F1, F2 (single-voxel) |
| 7  | (b) Rodrigues    | finalize `counts[fid] == 0` → write identity quat, `continue` (skip divide)    | `Rodrigues Analytical Oracle` — F0, F5 (zero-voxel) |
| 8  | (b) Rodrigues    | finalize `counts[fid] > 0` → divide, normalize, positiveOrientation, qu2eu     | `Rodrigues Analytical Oracle` — F2, F3, F4 |
| 9  | (c) vMF/Watson   | `featureNumVoxels[fid] == 0` → skip (output stays NaN from fill)               | both oracle tests — F0/F5 and background feature → NaN |
| 10 | (c) vMF/Watson   | `fzQuats.size() == 1` → muhat = single FZ quat, kappa = 0 (EM skipped)         | `Rodrigues Analytical Oracle` — F1, F2 (single-voxel) |
| 11 | (c) vMF/Watson   | `fzQuats.size() > 1` → `DirectionalStats::EMforDS` (VMF or WAT) → muhat, kappa; positiveOrientation; write quat/euler/kappa       | `vMF/Watson EbsdLib Reference Oracle` (22-quat feature) + `Rodrigues Analytical Oracle` F3/F4 |
| 12 | (a) Preflight    | all three methods disabled → preflight error `-54673` (runtime `-54670` kept as backstop, pinned by direct algorithm invocation)  | `No Method Enabled Error` |
| 13 | (b)+(c) Guards   | `Phases` value ≥ CrystalStructures tuple count → voxel excluded (no OOB read) + warning `-54672`| `Unknown Crystal Structure and Out-Of-Range Phase Guards` |
| 14 | (b) Rodrigues    | crystal structure ≥ ops count (e.g. 999/Unknown) → voxel excluded + warning `-54671`; fully-excluded feature → identity           | `Unknown Crystal Structure and Out-Of-Range Phase Guards` |
| 15 | (c) vMF/Watson   | crystal structure ≥ ops count → feature skipped (outputs stay NaN) + warning `-54671`           | `Unknown Crystal Structure and Out-Of-Range Phase Guards` |
| 16 | (c) vMF/Watson   | gather loop excludes phase-0 / out-of-range-phase voxels (gate identical to counting pass; issue #1659)          | `vMF/Watson Ignores Phase-0 Voxels` + guards test |

## Test inventory

| Test case | Status | Notes |
|-----------|--------|-------|
| `OrientationAnalysis::ComputeAvgOrientations: Rodrigues Analytical Oracle` | new-for-V&V | Class 1 + Class 4. 6 Triclinic Rodrigues fixtures with closed-form averages (exact quats @1e-6) + invariants; also asserts vMF/Watson value-add invariants (single→FZ-quat & kappa 0, zero→NaN) on the same small data. Inline hand-built data. |
| `OrientationAnalysis::ComputeAvgOrientations: vMF/Watson EbsdLib Reference Oracle` | new-for-V&V | Class 2 + Class 4. One Cubic feature of EbsdLib's 22 reference quaternions reproduces EbsdLib `DirectionalStatsTest`'s documented VMF/WAT muhat+kappa (quat margin 5e-3, kappa 2% to absorb float32 input round-trip). Inline hand-built data. |
| `OrientationAnalysis::ComputeAvgOrientations: Rodrigues Cubic Symmetry Invariant` | new (issue #1660) | Class 4 under Cubic_High: five representations of one orientation (Rz(30+90k), −Rz(30)) must average to Rz(30). Exercises the symmetry-reduction branch under 24 operators. |
| `OrientationAnalysis::ComputeAvgOrientations: Rodrigues Voxel Ordering Independence` | new (issue #1661) | Class 4: F3 fixture in both voxel orderings → identical Rz(45). Encodes the invariant the provenance previously claimed without a test. |
| `OrientationAnalysis::ComputeAvgOrientations: vMF/Watson Ignores Phase-0 Voxels` | new (issue #1659) | Regression pin for the phase-gate fix: a phase-0 voxel inside a feature must not enter the gather (single-voxel shortcut preserved; EM path lands on the valid orientation; cross-checked vs Rodrigues). |
| `OrientationAnalysis::ComputeAvgOrientations: Unknown Crystal Structure and Out-Of-Range Phase Guards` | new (issue #1661) | Guards: unknown crystal structure (999) and out-of-range `Phases` are excluded from both paths with warnings `-54671`/`-54672` (asserted); dropped features → identity/NaN. |
| `OrientationAnalysis::ComputeAvgOrientations: No Method Enabled Error` | new-for-V&V (updated for #1661) | No method enabled → preflight error `-54673` (asserted); runtime backstop `-54670` pinned via direct algorithm invocation. |
| `OrientationAnalysis::ComputeAvgOrientations: Cell Array Tuple Mismatch Error (-651)` | new-for-V&V (updated for #1661) | Preflight error path: mismatched cell-array tuple counts → `-651` (code asserted). The duplicate TEST_CASE that PR #1644 added to this file was removed in favor of this one. |
| `OrientationAnalysis::ComputeAvgOrientationsFilter: SIMPL Backwards Compatibility` | kept | `DYNAMIC_SECTION` over SIMPL 6.4 + 6.5 conversion fixtures; validates UUID + argument-key conversion only. Unaffected by V&V. |

| *(retired)* `OrientationAnalysis::ComputeAvgOrientations` | retired | Exemplar-comparison "valid execution" against `7_ComputeAvgOrientation_v2` (tol 5e-7). **Circular oracle** (exemplar regenerated from SIMPLNX output, PR #1577). Replaced by the Class 1/2/4 oracle tests above. |

**Coverage note (GCOV handoff, `.claude/gcov_ComputeAvgOrientations.md`):** the prior capture showed algorithm `.cpp` at 90.7% line / 46.5% branch, with the entire vMF/Watson path, the `-54670` branch, zero/single-voxel handling, and the VMF/Watson feature-array-detection branches **uncovered**. The V&V suite above now exercises all of them (Rodrigues, vMF, Watson, all-three-combined, both error paths, zero- and single-voxel features), closing the gap. The `-651` tuple branch (filter `.cpp`) is covered here; PR #1644's duplicate `-651` TEST_CASE was removed from this file (issue #1661).

## Exemplar archive

- **Archive:** `7_ComputeAvgOrientation_v2.tar.gz`
- **SHA512:** `2c2a691f1da301c449c20bafec65512d5134db38384ac7cb4c910880ccd87a260a5f011e905f35b97abff3952309f109c737c63ec3c833708926827a62a92efc`
- **Provenance:** `src/Plugins/OrientationAnalysis/vv/provenance/ComputeAvgOrientationsFilter.md` (documents the circular-oracle retirement and the inline replacement oracle)
- **Circular-oracle flag:** This `_v2` archive was regenerated from post-fix SIMPLNX output during PR #1577 (it bundles the filter's own `AvgQuats`, `Watson Avg Quats`, `vMF Avg Quats`, etc. as the comparison target). It is named explicitly in the audit's cross-cutting circular-oracle list. The V&V replaces this with oracle-derived expected values (Phase 10).

## Deviations from DREAM3D 6.5.171

*Comparison run 2026-06-30 against the **official DREAM3D 6.5.171 release** (`~/Applications/DREAM3D.app/Contents/bin/PipelineRunner`), Rodrigues method only (vMF/Watson have no legacy equivalent). Two byte-identical-input fixtures shared with NX `nxrunner`: **(A)** realistic — built from the `ASCIIData` CSVs (480k cells, 408 real features, single cubic phase); **(B)** a hand-built 3-feature fixture forcing the edge cases. **Result:** on fixture A the two are numerically equivalent on all 408 real features (`AvgQuats` ≤1e-6, zero sign flips; `AvgEulerAngles` max 4.77e-7); all divergences are confined to feature-0 / unindexed-voxel handling (D3 on A, D2 demonstrated on B). Full write-up: `vv/deviations/ComputeAvgOrientationsFilter.md` and the comparison summary in the archive (Phase 12).*

- `ComputeAvgOrientationsFilter-D1` — **Downgraded — NOT a deviation.** Northern-hemisphere canonicalization (`getPositiveOrientation`) could not be made to diverge (symmetry reduction keeps `w≥0`; the pure-triclinic edge is `q≡−q`, same orientation). Documented so the source difference isn't re-flagged. See deviations file.
- `ComputeAvgOrientationsFilter-D2` — **Demonstrated:** FeatureId-0 voxels (phase>0) — SIMPLNX averages them, 6.5.171 writes `(0,0,0,0)`. Algorithmic choice (not a bug). See deviations file.
- `ComputeAvgOrientationsFilter-D3` — **Confirmed legacy bug:** empty feature 0 is `(0,0,0,0)` in 6.5.171 vs clean identity `(0,0,0,1)` in SIMPLNX. See deviations file.
- `ComputeAvgOrientationsFilter-D4` — **Confirmed sub-epsilon:** library/precision, `AvgEulerAngles` max diff 4.77e-7. See deviations file.

*Post-merge adversarial review of PR #1645 filed three follow-up issues, all resolved 2026-07-08 and reflected
in the Test inventory and Code path coverage above: **#1659** (SIMPLNX bug — vMF/Watson gather ignored `Phases`;
fixed and pinned by `vMF/Watson Ignores Phase-0 Voxels`), **#1660** (oracle gap — added `Rodrigues Cubic
Symmetry Invariant`), and **#1661** (range-checked phase index with warnings `-54671`/`-54672`; all-methods-
disabled promoted to preflight `-54673`; error tests assert their codes; removed PR #1644's duplicate `-651`
test).*

## References

- **[R1]** Y.H. Chen, S.U. Park, D. Wei, G. Newstadt, M.A. Jackson, J.P. Simmons, M. De Graef, and A.O. Hero. "A Dictionary Approach to Electron Backscatter Diffraction Indexing." *Microscopy and Microanalysis* **21**(3), 739-752 (2015). DOI: [10.1017/S1431927615000756](https://doi.org/10.1017/S1431927615000756). Cited in the `dictmod.f90` module header as the derivation of the modified vMF / axial Watson indexing approach.
- **[R2]** Y.H. Chen, D. Wei, G. Newstadt, M. De Graef, J.P. Simmons, and A.O. Hero. "Parameter Estimation in Spherical Symmetry Groups." *IEEE Signal Processing Letters* **22**(8), 1152-1155 (2015). DOI: [10.1109/LSP.2014.2387206](https://doi.org/10.1109/LSP.2014.2387206). Derives the group-invariant (symmetry-modified) vMF distribution and its EM estimator — the mixture-over-symmetry-operators form implemented in `Estep_` / `Mstep_`.
- **[R3]** EMsoft, Marc De Graef Research Group, Carnegie Mellon University. Module `dictmod`, `Source/EMsoftLib/dictmod.f90`. BSD 3-Clause, Copyright (c) 2014-2022. <https://github.com/EMsoft-org/EMsoft>
- **[R4]** EbsdLib PR #41, commit `006a023` (2026-04-03), "ENH: Implement Watson and von-Mises Fisher Orientation Averaging" — the C++ port that introduced `ebsdlib::DirectionalStats`.
- **[R5]** K.V. Mardia and P.E. Jupp. *Directional Statistics*. Wiley, Chichester (2000). General background for the von Mises-Fisher and Watson distributions.
