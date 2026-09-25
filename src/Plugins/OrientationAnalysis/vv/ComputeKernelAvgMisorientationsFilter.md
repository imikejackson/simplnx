# V&V Report: ComputeKernelAvgMisorientationsFilter

|           |                          |
|-----------|--------------------------|
| Plugin    | OrientationAnalysis      |
| SIMPLNX UUID               | `61cfc9c1-aa0e-452b-b9ef-d3b9e6268035`   |
| SIMPLNX Human Name         | Compute Kernel Average Misorientations   |
| DREAM3D 6.5.171 equivalent | `FindKernelAvgMisorientations` — `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindKernelAvgMisorientations.{h,cpp}` (UUID `88d332c1-cf6c-52d3-a38d-22f6eae19fa6`) |
| Verified commit            | `a307946e7` (v7.4.2 release)  |
| Status | COMPLETE     |
| Sign-off  | Michael Jackson <mike.jackson@bluequartz.net> — 2026-07-29. *Reopened 2026-07-15 for the `use_feature_ids` feature (issue #1613, branch `topic/kam_ignore_feature_ids`), delivered by PR #1674; the prior 2026-06-03 sign-off is superseded.* |
| Second-engineer sign-off   | Jared Duffey — 2026-07-29 (approving reviewer, PR #1674, which supersedes the PR #1631 review by Nathan Young of 2026-06-11) |

## At a glance

| Aspect                 | Current state            |
|------------------------|--------------------------|
| Algorithm Relationship | **Port** with a corrected legacy kernel-bound bug (D2), precision/library changes (D1), and the NX-only `use_feature_ids=false` mode (D3). |
| Oracle (confirmed)     | **Class 1** uses 6 hand-derived KAM fixtures; **Class 4** verifies range, background, uniform-data, and mode-equivalence invariants. All 9 tests pass. |
| Code paths enumerated  | 8 of 8 paths exercised after removal of the unreachable `numVoxel==0` fallback. |
| Tests today            | 9 test cases: 6 analytical fixtures, 2 invariant tests, and 1 SIMPL conversion test. |
| Exemplar archive       | None; fixtures are inline. The circular `6_6_stats_test_v2.tar.gz` KAM exemplar was retired for this filter but remains shared by other tests. |
| Legacy comparison      | **Run** — six cases on synthetic and Small IN100 inputs establish precision D1 and asymmetric-kernel bug D2; the corrected local legacy build agrees within 9.54e-7° in every case. |
| Bug flags              | `ComputeKernelAvgMisorientationsFilter-D2` — the legacy x-loop used the z-radius as its upper bound; SIMPLNX corrected it at port time. |
| V&V phase | **COMPLETE.** |

## Summary

`ComputeKernelAvgMisorientationsFilter` computes the per-cell **Kernel Average Misorientation (KAM)**: for each valid cell (featureId > 0, phase > 0), the algorithm iterates over an axis-aligned kernel of user-specified radius `KernelSize = (x, y, z)`, averages the misorientation between the focal cell's orientation quaternion and every same-feature neighbor's quaternion (including the focal cell itself, which contributes a self-misorientation of 0°), and stores the result. Cells with `featureId == 0` or `phase == 0` are treated as background and assigned KAM = 0 directly.

The filter is the cell-level analog of the feature-level `ComputeFeatureNeighborMisorientationsFilter`. Like that filter, it consumes `Quats` (cell-level avg-orientations) and `CrystalStructures` (per-phase Laue class index), and delegates the actual cubic/hex/etc. symmetry-reduced disorientation calculation to `ebsdlib::LaueOps::calculateMisorientation()`. Unlike the feature-level filter, the kernel iteration ALWAYS visits the focal cell as part of its neighbor list (via the `j=k=l=0` inner iteration), so the per-cell divisor `numVoxel` is always ≥ 1.

The output is `KernelAverageMisorientations`, a `Float32Array` co-located in the same `AttributeMatrix` as the input `Cell Data` arrays, sized one-tuple-per-cell with one component per tuple, in degrees.

**`use_feature_ids` feature (issue #1613):** a new `BoolParameter` (default `true`) selects the neighbor-inclusion rule. With `use_feature_ids = true` (the default, behavior-identical to the legacy filter and to every prior release), a kernel neighbor contributes only if it shares the focal cell's `featureId` — the classic per-grain KAM. With `use_feature_ids = false`, the filter computes a **per-voxel** KAM: a neighbor contributes iff it is in-bounds, has `featureId > 0`, and shares the focal cell's phase — feature boundaries are ignored, so the kernel averages across grain boundaries within the same phase. The focal-validity gate (`featureIds[point] > 0 && cellPhases[point] > 0`) is identical in both modes.

## Algorithm Relationship

*Classification:* **Port (with UUID reassignment + name rename + one inherent legacy bug corrected at port time), plus one NX-only feature addition (`use_feature_ids`, issue #1613).**

*Evidence:* Cross-checked SIMPLNX algorithm against legacy `FindKernelAvgMisorientations.cpp::execute()`. Same per-voxel outer loop, same per-kernel inner loop, same `KAM = totalMisorientation / numVoxel`. On the **default per-grain path (`use_feature_ids = true`)** the neighbor gate is identical: legacy `if(good && m_FeatureIds[point] == m_FeatureIds[neighbor])` (line 292) ≡ SIMPLNX `featureIds[point] == featureIds[neighborIdx]`. Port-time deltas:

- `QuatF` → `QuatD`; `getMisoQuat` → `calculateMisorientation` (EbsdLib 2.4.1+ `2·atan2(|v|, w)` precision form — see D1). The `float32`→`float64` promotion of the misorientation math is the second component of the precision family observed in the runtime A/B (below).
- `setParallelizationEnabled` removed (now always parallel via `ParallelData3DAlgorithm`).
- Iteration order `col→row→plane` → `plane→row→col` (cache-friendlier for x-fastest-varying storage; mathematically identical since all writes go to the same `point` index).
- D2 corrected at port.
- UUID reassigned; `Find` → `Compute` rename.

*NX-only feature addition (issue #1613):* the `use_feature_ids` BoolParameter. When `true` (default) the code path is bit-for-bit the pre-feature per-grain algorithm — the feature is opt-in and does not alter the default output. When `false`, the neighbor gate becomes `featureIds[neighborIdx] > 0 && cellPhases[neighborIdx] == cellPhases[point]` (per-voxel, phase-gated, feature-agnostic). This is a **new capability with no DREAM3D 6.5.171 counterpart** — validated by the Class 1 / Class 4 oracle only, never by legacy comparison (see D3).

*Material PRs since baseline:* commits `25959c1f2` through `3e863da2f` (inclusive) on branch `topic/kam_ignore_feature_ids` add the `use_feature_ids` parameter, the per-mode neighbor gate, three unit tests, and the user-doc update, followed by a small follow-up commit (`3e863da2f`) correcting the A/B tallies in this V&V report. No other change to the core algorithm.

*Review-driven cleanup (2026-07-16, commit `7f9cddc7d`):* a review-flagged, behavior-preserving cleanup pass refactored the kernel boundary-index handling to signed comparisons (the neighbor index is now computed directly from the clamped `zIdx`/`yIdx`/`xIdx`, removing the vacuous unsigned `< 0` checks and the separate `neighbor<0` guard), hoisted `LaueOps` construction out of the per-cell hot path into a worker member built once per run (constructor-time, shared read-only across parallel ranges — same pattern as `ComputeFeatureFaceMisorientationPerTriangleImpl`), tightened cancellation checking from per-plane to per-row granularity, floored the progress-message increment at 1, and removed the dead `numVoxel==0` fallback (folding the invalid-cell reset into an `else` branch, logically equivalent by De Morgan's law since the two branch conditions are complementary for these non-negative `int32` arrays). No algorithmic or numeric change — verified by the unchanged 9-test suite passing.

## Oracle

*Class:* **1 (Analytical) primary, 4 (Invariant) companion.**

*Applied:* Six hand-derived fixtures use pure φ1 rotations to calculate expected KAM values. Companion invariants check output bounds, background values, uniform data, and equivalence of the two modes when their neighbor sets are equal.

### Class 1 (Analytical)

Class 1 oracle derived by hand for each fixture in terms of `|Δφ1|` between cell pairs, justified by the cubic FZ analysis below. All quaternions in the test fixtures use the helper `QuatFromPhi1Deg(phi1)` which returns the quaternion form of a pure Bunge ZXZ Euler rotation `(phi1, 0, 0)` with `Phi = phi2 = 0`. This collapses to a single rotation about the z-axis.

**Cubic FZ argument:** For two cells with pure φ1 rotations differing by Δφ1, the disorientation between them in the cubic FZ equals `|Δφ1|` whenever `|Δφ1| ≤ 45°`. The reasoning: the cubic group's 4-fold rotation about the z-axis reduces φ1 differences modulo 90°; for |Δφ1| ≤ 45°, the reduction is the identity operator (no reduction needed). For all 4 Class 1 fixtures, the maximum φ1 difference between any pair is ≤ 30°, well within the 45° bound. Other cubic symmetry operators (3-fold about [111], 2-fold about [110], etc.) only produce smaller-angle equivalents for misorientations not aligned with a pure z-axis rotation; for pure z-rotations of small magnitude, the identity is the global minimum.

**Per-fixture expected KAM derivation:**

| Fixture                | Geometry | Kernel    | Per-cell expected KAM (degrees)        |
|--------------------------------------------|----------|-----------|------------------------------------------------------------|
| `Class 1 - Uniform 2D Single Feature`      | 3x3x1    | {1,1,0}   | All cells = 0.0 (all in-kernel neighbors share orientation) |
| `Class 1 - 1D x-axis Gradient`             | 5x1x1    | {1,0,0}   | [2.5, 10/3, 10/3, 10/3, 2.5] (see derivation in test comments) |
| `Class 1 - 1D z-axis Gradient (3D path)`   | 1x1x3    | {0,0,1}   | [5.0, 20/3, 5.0]    |
| `Class 1 - Multi-Feature Multi-Voxel + BG` | 6x1x1    | {1,0,0}   | [5.0, 5.0, 10.0, 10.0, 0.0, 0.0]        |
| `Class 1 - Per-Voxel Mode` (`use_feature_ids=false`) | 6x1x1 | {1,0,0} | [5.0, 20/3, 10.0, 10.0, 0.0, 0.0] |
| `Class 1 - Per-Voxel Mode Two-Phase Gates` (`use_feature_ids=false`) | 5x1x1 | {1,0,0} | [5.0, 5.0, 0.0, 0.0, 0.0] |

Detailed per-cell hand-derivations are in the test file's TEST_CASE comments and in `vv/provenance/ComputeKernelAvgMisorientationsFilter.md`.

**Per-voxel-mode fixture derivations (issue #1613):**

*`Class 1 - Per-Voxel Mode`* — the same 6×1×1 layout as the Multi-Feature fixture (featureIds `[1,1,2,2,0,1]`, phases `[1,1,1,1,0,1]`, φ1 `[0,10,0,20,–,30]°`), kernel `{1,0,0}`, run with `use_feature_ids=false`. A neighbor now contributes iff `featureId>0 && phase==focalPhase`, so feature boundaries are crossed but `featureId=0` and phase-mismatched cells are still excluded:
- cell 0: {self=0, x1(F1,P1)=|10−0|} → 10/2 = **5.0** (same as per-grain).
- cell 1: {x0(F1,P1)=|0−10|, self=0, x2(F2,P1)=|0−10|} → 20/3 ≈ **6.667** — the *diagnostic* cell: per-grain gives 5.0 (x2 skipped as different feature); per-voxel includes x2 because it is a valid same-phase cell. This one value proves the mode actually changed behavior.
- cell 2: {x1(F1,P1)=|10−0|, self=0, x3(F2,P1)=|20−0|} → 30/3 = **10.0**.
- cell 3: {x2(F2,P1)=|0−20|, self=0}; x4 excluded (`featureId=0`) → 20/2 = **10.0**.
- cell 4: focal-invalid (`featureId=0, phase=0`) → **0.0** exactly.
- cell 5: {self=0}; x4 excluded (`featureId=0`) → 0/1 = **0.0**.

*`Class 1 - Per-Voxel Mode Two-Phase Gates`* — 5×1×1, every cell its own feature, two cubic phases (ensemble `[999,1,1]`); featureIds `[1,2,3,4,0]`, phases `[1,1,2,1,1]`, φ1 `[0,10,20,30,40]°`, kernel `{1,0,0}`, `use_feature_ids=false`:
- cell 0: {self=0, x1(P1)=|0−10|} → 10/2 = **5.0** (per-grain would give 0.0 — every feature is a single cell).
- cell 1: {x0(P1)=|10−0|, self=0}; x2 skipped (phase 2 ≠ 1) → 10/2 = **5.0**.
- cell 2: {self=0}; x1, x3 skipped (phase 1 ≠ 2) → 0/1 = **0.0**.
- cell 3: {self=0}; x2 skipped (phase), x4 skipped (`featureId=0`) → 0/1 = **0.0**.
- cell 4: focal `featureId=0` (invalid) → **0.0** exactly, even though its phase>0 — proves the focal gate is unchanged in per-voxel mode.

### Class 4 (Invariant)

Class 4 invariants asserted in the `Class 4 - Invariants` TEST_CASE across 3 sub-sections:

1. **Uniform-orientation single-feature → KAM == 0 everywhere.** Asserted on a 3x3x3 uniform-identity-quaternion fixture with kernel `{1,1,1}` (3D path coverage).
2. **Background cell → KAM == 0 exactly.** Asserted on a 3x1x1 fixture with the middle cell flagged as `(featureId=0, phase=0)`.
3. **Range and non-triviality on the x-axis gradient fixture:** (i) `KAM[i] >= 0` for all cells, (ii) `KAM[i] <= 62.8°` (Mackenzie cubic upper bound), (iii) at least one cell has `KAM > 0` (sanity check that the algorithm actually computed something).

**Mode-equivalence invariant (`Class 4 - Mode Equivalence on Single Feature`):** on single-feature single-phase data, the per-grain gate (`featureId` match) and the per-voxel gate (`featureId>0 && phase match`) admit *exactly the same* neighbor set for every focal cell, so the two modes must produce **bit-for-bit identical** output. Asserted on a 3×3×3 gradient fixture (φ1 = 2x+3y+4z°, max 18° < 45° FZ bound) with kernel `{1,1,1}`: `REQUIRE(kamPerGrain[i] == kamPerVoxel[i])` for every cell, plus a non-triviality guard that at least one cell is non-zero. This is a derived property (not a hand-computed value), so it holds regardless of the EbsdLib precision class and pins the invariant "the feature is opt-in and does not change output where the two gates coincide."

The Class 4 invariants are oracle-agnostic — they hold for any input, so they catch regressions even if specific Class 1 expected values were edited away.

*Encoded:* `test/ComputeKernelAvgMisorientationsTest.cpp` — 6 Class 1 fixtures and 2 Class 4 tests; all 9 active test cases pass.

### Class 2, 3, 5

N/A — no reference-library invocation (Class 2), no published-paper figure reproduction (Class 3), no expert-visual sign-off (Class 5) needed. Class 1 + Class 4 are sufficient.

*Second-engineer review:* **Jared Duffey — 2026-07-29** (approving reviewer of PR #1674, which delivered the reopened `use_feature_ids` work). This supersedes the earlier review of PR #1631 by Nathan Young (2026-06-11).

The area flagged for the second pair of eyes during the cycle was the multi-feature multi-voxel fixture, which carries 6 hand-derived per-cell expected values resting on symmetry-reduced cubic misorientation reasoning.

## Bugs found and fixed

| Deviation | Defect | Affected released versions | Resolution in this branch |
|-----------|--------|----------------------------|---------------------------|
| `ComputeKernelAvgMisorientationsFilter-D2` | The legacy x-loop used `KernelSize.z` for its upper bound and produced the wrong kernel shape when the x- and z-radii differed. | DREAM.3D 6.5.171 only. DREAM3D-NX was not affected. | SIMPLNX uses `KernelSize.x` for both x-loop bounds and has done so since the port. |

## Code path coverage

*8 of 8 paths exercised.*

Source: `src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/Algorithms/ComputeKernelAvgMisorientations.cpp` (209 lines).

Re-enumerated for the per-mode neighbor gate (`ComputeKernelAvgMisorientations.cpp:124`). Paths 3–4 are the per-grain branch (`use_feature_ids = true`); paths 5–7 are the per-voxel branch (`use_feature_ids = false`). **Updated 2026-07-16 (commit `7f9cddc7d`):** the former path 9 (`numVoxel==0` fallback) was removed as unreachable dead code — see path 8's note. Path 2's boundary clamp is now a signed-index comparison (`zIdx`/`yIdx`/`xIdx` at `.cpp:100,107,114`) with the neighbor index computed directly from the clamped indices (`.cpp:120`); the separate `neighbor < 0` guard no longer exists because it is now provably unreachable by construction rather than checked at runtime.

| Path | Description      | Exercised by |
|------|------------------------------------------------------------------------------|--------------|
| 1    | Focal-valid gate (`featureIds[point] > 0 && cellPhases[point] > 0`) → enter kernel               | All 6 Class 1 fixtures; both Class 4 tests |
| 2    | Kernel cell out-of-bounds (signed boundary clamp on `zIdx`/`yIdx`/`xIdx`) → `continue`           | `Class 1 - 1D x-axis Gradient` cells 0/4 (x-boundary); `Class 1 - 1D z-axis Gradient` planes 0/2 (z-boundary); `Class 1 - Uniform 2D` corners; both 3×3×3 fixtures (all faces) |
| 3    | **Per-grain** (`use_feature_ids=true`): in-bounds neighbor + `featureId` match → accumulate miso + numVoxel++        | All default-mode Class 1 fixtures; Class 4 Invariants (i)/(iii); Class 4 Mode-Equivalence (per-grain run) |
| 4    | **Per-grain**: in-bounds neighbor + `featureId` mismatch → skip (no accumulate) | `Class 1 - Multi-Feature Multi-Voxel + Background` (cells 1, 2 see different-feature in-bounds neighbors) |
| 5    | **Per-voxel** (`use_feature_ids=false`): in-bounds neighbor + `featureId>0` + phase match → accumulate               | `Class 1 - Per-Voxel Mode` cells 1, 2 (cross-feature same-phase include); `Class 1 - Per-Voxel Two-Phase Gates` cells 0, 1; Class 4 Mode-Equivalence (per-voxel run) |
| 6    | **Per-voxel**: neighbor `featureId == 0` → skip          | `Class 1 - Per-Voxel Mode` cells 3, 5 (x=4 background excluded); `Class 1 - Per-Voxel Two-Phase Gates` cell 3 (x=4 excluded) |
| 7    | **Per-voxel**: neighbor phase mismatch (`cellPhases[neighbor] != cellPhases[point]`) → skip      | `Class 1 - Per-Voxel Two-Phase Gates` cell 2 (excludes phase-1 x=1/x=3), cells 1 & 3 (exclude phase-2 x=2) |
| 8    | Focal-invalid (`featureIds == 0 \|\| cellPhases == 0`, reached via `else` since commit `7f9cddc7d`) → KAM = 0 directly. The former path 9 (`numVoxel == 0` fallback) was removed by this commit: in both modes the focal voxel always self-contributes when the focal is valid (per-grain: `featureId==featureId`; per-voxel: focal-valid ⇒ `featureId>0` and `phase==phase`), guaranteeing numVoxel ≥ 1 whenever path 1 is entered, so the fallback was unreachable by construction. | `Class 1 - Multi-Feature` cell 4; `Class 1 - Per-Voxel Mode` cell 4; `Class 1 - Per-Voxel Two-Phase Gates` cell 4 (`featureId=0` yet phase>0 — proves focal gate unchanged); Class 4 Invariants (ii) |

8 of 8 paths exercised by the V&V suite. The former path 9 (`numVoxel==0` fallback, previously "unreachable by construction") was removed from SIMPLNX by the review-driven cleanup commit `7f9cddc7d` (2026-07-16); it was purely dead code and its removal is behavior-preserving, confirmed by the unchanged 9-test suite passing.

## Test inventory

| TEST_CASE | Category | Status | ctest entry                |
|--------|----------|--------|--------|
| `: SIMPL Backwards Compatibility`              | Compat   | kept   | Yes (2 dynamic sections: 6.4 + 6.5)            |
| `: Class 1 - Uniform 2D Single Feature`        | Class 1  | kept   | Yes    |
| `: Class 1 - 1D x-axis Gradient`               | Class 1  | kept   | Yes    |
| `: Class 1 - 1D z-axis Gradient (3D path)`     | Class 1  | kept   | Yes    |
| `: Class 1 - Multi-Feature Multi-Voxel + Background`               | Class 1  | kept   | Yes    |
| `: Class 1 - Per-Voxel Mode (use_feature_ids = false)`             | Class 1  | **new (#1613)** | Yes — per-voxel cross-feature include + `featureId=0` exclude; expected `{5.0, 20/3, 10.0, 10.0, 0, 0}` |
| `: Class 1 - Per-Voxel Mode Two-Phase Gates`   | Class 1  | **new (#1613)** | Yes — per-voxel phase-mismatch exclude + `featureId=0` focal gate; expected `{5.0, 5.0, 0, 0, 0}` |
| `: Class 4 - Mode Equivalence on Single Feature`  | Class 4  | **new (#1613)** | Yes — per-grain ≡ per-voxel bit-for-bit on single-feature single-phase 3×3×3 |
| `: Class 4 - Invariants` (3 sub-sections)      | Class 4  | kept   | Yes (3 SECTIONs)           |
| ~~`: ComputeKernelAvgMisorientationsFilter` (legacy exemplar test)~~  | RETIRED  | —      | Retired 2026-06-03 (circular oracle from pre-EbsdLib-2.4.1 SIMPLNX output)             |

9 active TEST_CASEs / 9 ctest entries (confirmed via `ctest -N -R "ComputeKernelAvgMisorientations"`); 100% pass (full regression sweep, 2026-07-15).

## Exemplar archive

**None** — inline-constructed. The pre-V&V test (now retired) consumed `6_6_stats_test_v2.tar.gz` (SHA512 `e84999...089723`, downloaded from the BlueQuartz Data_Archive release). The archive contains exemplar `KernelAverageMisorientations` arrays generated from a pre-2.4.1 SIMPLNX build, which embeds the spurious self-misorientation precision noise described in D1. Comparing the post-2.4.1 SIMPLNX output against those exemplars fails by ~0.01-0.05° per cell on real Small_IN100 data. The exemplar arrays cannot be re-generated against the post-2.4.1 build (circular oracle pattern), so the test was retired and replaced with the analytical / invariant suite above.

The shared archive remains referenced in `src/Plugins/OrientationAnalysis/test/CMakeLists.txt` (line 130) for use by `AlignSectionsMutualInformation`, `ComputeShapes`, and `ComputeSchmids` tests. Only F#5's consumption line was removed.

- **Provenance:** `vv/provenance/ComputeKernelAvgMisorientationsFilter.md` — the canonical record of how the inlined data fixtures were designed and how the expected values were derived.

## Deviations from DREAM3D 6.5.171

See `vv/deviations/ComputeKernelAvgMisorientationsFilter.md` for the canonical, ID-stable list:

- **`ComputeKernelAvgMisorientationsFilter-D1`** — Orientation-library precision difference. On Small IN100 with the symmetric kernel, DREAM3D 6.5.171 differs from SIMPLNX by up to 0.056°; 2,281 of 256,000 cells exceed 0.01°. The local legacy build with the surgical precision correction reduces the maximum residual to 9.54e-7°.
- **`ComputeKernelAvgMisorientationsFilter-D2`** — Legacy x-loop upper-bound bug. The 2026-09-17 asymmetric-kernel rerun exercises both over-reach (`{1,1,2}`) and truncation (`{2,1,1}`); differences reach 1.16° on Small IN100. Applying the one-character correction to a local legacy build reduces every asymmetric case to at most 9.54e-7° with zero cells above 0.01°, proving the root cause.
- **`ComputeKernelAvgMisorientationsFilter-D3`** — `use_feature_ids = false` (per-voxel KAM) is an **NX-only capability** added for issue #1613. DREAM3D 6.5.171 `FindKernelAvgMisorientations` has no equivalent (it is per-grain only), so there is nothing to A/B against; the mode is validated by the Class 1 per-voxel fixtures and the Class 4 mode-equivalence invariant. The default (`use_feature_ids = true`) is unchanged and remains legacy-comparable.
