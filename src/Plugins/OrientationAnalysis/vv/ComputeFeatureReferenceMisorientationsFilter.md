# V&V Report: ComputeFeatureReferenceMisorientationsFilter

|           |                          |
|-----------|--------------------------|
| Plugin    | OrientationAnalysis      |
| SIMPLNX UUID               | `24b54daf-3bf5-4331-93f6-03a49f719bf1`  |
| SIMPLNX Human Name         | Compute Feature Reference Misorientations                |
| DREAM3D 6.5.171 equivalent | `FindFeatureReferenceMisorientations` — `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindFeatureReferenceMisorientations.{h,cpp}` (UUID `428e1f5b-e6d8-5e8b-ad68-56ff14ee0e8c`) |
| Verified commit            | `a307946e7` (v7.4.2 release) |
| Status | COMPLETE     |
| Sign-off  | Michael Jackson <mike.jackson@bluequartz.net> — 2026-06-03 (PR #1629 author) |
| Second-engineer sign-off   | Nathan Young — 2026-06-03 (approving reviewer, PR #1629)   |

## At a glance

| Aspect                 | Current state            |
|------------------------|--------------------------|
| Algorithm Relationship | **Port** — The two reference modes and per-voxel math are preserved. SIMPLNX adds API updates, linear iteration, cancel checks, and the optional `EuclideanCenters` output. |
| Oracle (confirmed)     | **Class 1 (Analytical)** uses 6 hand-derived 2D, 3D, multi-feature, and edge-case fixtures; **Class 4 (Invariant)** checks bounds, skip behavior, and feature averages. |
| Code paths enumerated  | 7 of 8 paths exercised directly; only the cancellation path is not directly tested. |
| Tests today            | 8 test cases cover 6 Class 1 fixtures, 1 Class 4 invariant sweep, and SIMPL conversion. No exemplar archive is used. |
| Exemplar archive       | **None** — the circular `compute_feature_reference_misorientation.tar.gz` test is retired and replaced by inline fixtures. |
| Legacy comparison      | **Run** — Six analytical fixtures and Small IN100 confirmed D1's precision difference and D2's Mode 1 crash; the local legacy proof build matched SIMPLNX at the stated precision. |
| Bug flags              | `ComputeFeatureReferenceMisorientationsFilter-D2` is a legacy Mode 1 null-pointer defect; SIMPLNX is not affected. |
| V&V phase | **COMPLETE.** |

## Summary

`ComputeFeatureReferenceMisorientationsFilter` compares each cell orientation with either its feature average or its feature's farthest interior voxel and also computes feature averages. Six Class 1 fixtures and a Class 4 invariant sweep verify both modes, including 2D, 3D, multi-feature, and edge-case data. The comparison confirms D1's precision difference and D2's legacy Mode 1 crash; the local legacy proof build reproduces SIMPLNX on the analytical cases and to float32 precision on production data.

## Algorithm Relationship

*Classification:* **Port (with UUID reassignment and API modernization)** ~~| Minor changes | Rewrite | New filter~~

*Evidence:* The SIMPLNX algorithm at `Algorithms/ComputeFeatureReferenceMisorientations.cpp` (~175 lines) is a near line-by-line translation of legacy `FindFeatureReferenceMisorientations::execute()` (DREAM3D 6.5.171, ~110 lines). Same two-mode dispatch (`ReferenceOrientation` parameter), same per-voxel main loop computing misorientation via `LaueOps`, same per-feature averaging finalization. The SIMPLNX filter was assigned a **new UUID** (`24b54daf-3bf5-4331-93f6-03a49f719bf1` vs legacy `428e1f5b-e6d8-5e8b-ad68-56ff14ee0e8c`) for the `Find` → `Compute` rename; SIMPL 6.4/6.5 pipelines still open correctly via the conversion fixtures at `test/simpl_conversion/6_*/`.

*Port-time deltas (non-deviation — preserve algorithmic equivalence at hand-built-fixture precision):*

1. **EbsdLib API**: `getMisoQuat` → `calculateMisorientation` (axis + angle returned together; underlying `LaueOps` math identical).
2. **Quaternion precision**: `QuatF` → `QuatD` inside the algorithm (float32 inputs promoted, output cast back to float32). Equivalent for non-precision-sensitive inputs; the EbsdLib 2.4.1 CubicOps fix is visible for cubic sym-op-aligned boundaries (see D1).
3. **Voxel iteration**: `for col/row/plane` triple loop → `for voxelIdx` linear loop. Equivalent iteration order; raster-order tie-break (`>=`, later-voxel-wins) preserved in `m_Centers` selection.
4. **Auxiliary storage**: legacy `avgMisoPtr` (interleaved count + sum) → two separate vectors `avgMisorientationSums` + `avgMisorientationCounts`. Equivalent semantics.
5. **Cancel checks** added (UX only; no algorithmic effect on completed runs).
6. **New optional `EuclideanCenters` array** (Mode 1 only) — SIMPLNX feature, not a port artifact; does not affect pre-existing output.
7. **EbsdLib 2.4.1 CubicOps precision improvement** — non-observable on V&V data fixtures (no sym-op-aligned features); ~ULP-scale per-feature drift on real EBSD data. See D1 + `BadDataNeighborOrientationCheckFilter`'s V&V cycle.

*Material PRs since baseline (2025-10-01):* None identified that materially change this filter's algorithm. PR #1472 (EbsdLib 2.0.0 API bump) is the closest, and that just affects the `getMisoQuat` → `calculateMisorientation` API delta noted above.

## Oracle

*Class:* **1 (Analytical)** primary + **4 (Invariant)** companion. Class 3 (Paper-based) N/A — this filter delegates misorientation math to `ebsdlib::LaueOps::calculateMisorientation`; the Rowenhorst 2015 paper-based verification of that math is part of EbsdLib's own V&V, not this filter's.

*Applied:* The Class 1 and Class 4 oracles are applied as described below.

### Class 1 — Analytical

Expected per-voxel `FRM` and per-feature `avgRefMis` outputs are derived in closed form from the input `Quats` + `Phases` + `FeatureIds` + reference-quaternion source (Mode 0: `AvgQuats[fid]`; Mode 1: `Quats[centerVoxelIdx]`) by hand-tracing the algorithm. The fixtures use pure φ1-rotation quaternions (Bunge ZXZ Euler `(φ1, 0, 0)`) so that misorientation between any two voxels equals `|Δφ1|` modulo the cubic c-axis 4-fold symmetry. For Δφ1 ∈ {0°, 5°, 10°} the symmetry reduction is the identity (no fold below the input), so expected FRM values are exactly `|Δφ1|`. Per-feature averages are `sum(FRM[v ∈ feature fid, phase>0]) / count(v ∈ feature fid, phase>0)` for `fid > 0` with non-empty count, or `0` when count is empty (path 7 in the code-path table below).

Mode 1 hand-picks `GBEuclideanDistances` values so that `m_Centers[fid]` selection has a unique closed-form answer (or, for the tied-distance multi-feature fixture, a deterministic later-voxel-wins tie-break per the `>=` comparison semantics that both legacy and SIMPLNX share).

### Class 4 — Invariant

Five invariants every filter run must satisfy regardless of input configuration, asserted via `namespace AnalyticalFixtures::AssertClass4Invariants()` in the test source:

- **Non-negativity**: `FRM[i] >= 0` ∀ voxel
- **Cubic max-angle bound**: `FRM[i] <= 62.8°` for cubic phases (the maximum symmetry-reduced misorientation under m-3m symmetry)
- **Skip-path correctness**: `FRM[i] == 0` when `featureIds[i] == 0` OR `cellPhases[i] == 0`
- **Background-feature zero**: `avgRefMis[0] == 0` (background)
- **Per-feature averaging formula**: `avgRefMis[fid] == sum(FRM[v ∈ feature fid, phase>0]) / count(v ∈ feature fid, phase>0)` for `fid > 0` with count > 0; `avgRefMis[fid] == 0` when count == 0

*Encoded:* The tests below encode the oracle.

- **Class 1 (Analytical)**: `test/ComputeFeatureReferenceMisorientationsTest.cpp` — 6 `TEST_CASE` blocks under the `Class 1 - …` family. Per-voxel and per-feature expected values asserted via `AnalyticalFixtures::RequireFRMClose()` / `RequireAvgClose()` with 1e-3° tolerance (degrees) and `Approx().margin(1e-5f)` for `EuclideanCenters` coord assertions.
- **Class 4 (Invariant)**: `ComputeFeatureReferenceMisorientationsFilter: Class 4 - Invariants Sweep` — two configurations (Mode 0 mixed 3×3×1 and Mode 1 3×3×1) each asserting all five invariants via the `AssertClass4Invariants()` helper.
- *(kept)* `ComputeFeatureReferenceMisorientationsFilter: SIMPL Backwards Compatibility` — SIMPL 6.4 + 6.5 conversion paths via `DYNAMIC_SECTION`; UUID + argument-key + parameter-value validation only.

*Second-engineer review:* Nathan Young — 2026-06-03 (approving reviewer, PR #1629).

Recorded review topics:

- *The Class 1 hand-derivations in the 6 data fixtures + 1 invariants sweep for plausibility (the fixtures are small enough to walk through in ~30 minutes).*
- *The Class 4 invariant set for completeness — are there other properties this algorithm must satisfy?*
- *The decision to retire the `compute_feature_reference_misorientation.tar.gz` Small-IN100 exemplar archive in favor of inline data fixtures.*

## Bugs found and fixed

| Deviation | Defect | Affected released versions | Resolution in this branch |
|-----------|--------|----------------------------|---------------------------|
| `ComputeFeatureReferenceMisorientationsFilter-D2` | DREAM3D 6.5.171 dereferences a Mode-0-only quaternion pointer in Mode 1 and terminates before calculation. | DREAM.3D 6.5.171 only. DREAM3D-NX was not affected. | SIMPLNX obtains the feature count without the Mode-0-only input and completes both modes. |

## Code path coverage

*7 of 8 paths exercised directly; 1 (cancel) implicitly via the unconditional cancel-check-at-loop-top instrumentation.*

Source: `src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/Algorithms/ComputeFeatureReferenceMisorientations.cpp` (185 lines).

The algorithm has three logical phases: (a) Mode 1 pre-loop (populate `centers[fid]` from `gbEuclideanDistances` + write `EuclideanCenters` coords); (b) main per-voxel loop (compute per-voxel misorientation, accumulate per-feature sums + counts); (c) per-feature finalize (compute per-feature average from sum/count).

| #  | Phase             | Path        | Test case       |
|----|-------------------|-------------------------------------------------------------------------------------------------------------------|----------------------------------------|
| 1  | (b) Main loop     | Mode 0 — `q2 = avgQuatsPtr[featureId * 4 .. ]`               | Fixtures A, B, C (`Class 1 - Mode 0 SingleGrainIdentity` / `KnownAngle5deg` / `MultiGrain EdgeCases`)               |
| 2  | (a) Mode 1 pre-1  | Mode 1 — first sweep populates `centers[fid]` via `>=` tie-break against `gbEuclideanDistances` | Fixtures D, E, F (`Mode 1 KnownCenter` / `MultiGrain CenterIsolation` / `3D Volume`)              |
| 3  | (a) Mode 1 pre-2  | Mode 1 — second sweep writes `EuclideanCenters[fid] = imageGeom.getCoordsf(centers[fid])`      | Fixtures D, E, F               |
| 4  | (b) Main loop     | `featureIds[i] > 0 && cellPhases[i] > 0` → compute misorientation, accumulate `sums[fid]++` and `counts[fid]+=`  | All Class 1 fixtures            |
| 5  | (b) Main loop     | Skip (`featureIds[i] == 0` OR `cellPhases[i] == 0`) → FRM stays 0 (initial `fill(0.0f)` value) | Fixture C (`Mode 0 MultiGrain EdgeCases`) — background voxel 0, mid-feature unphased voxels 6-7, entire feature 4 unphased voxels 8-11|
| 6  | (b) Main loop     | Cancel check at loop top (`m_ShouldCancel.load()` → early return)              | *Not directly tested.* Unconditional check at loop top in both passes; failure mode (silent cancel-disregard) would manifest as test hang in any test, but is not specifically exercised. Low-value gap. |
| 7  | (c) Finalize      | `avgMisorientationCounts[fid] == 0` → `avgRefMis[fid] = 0`    | Fixture C — feature 4 has all-unphased voxels, so `count[4] == 0` |
| 8  | (c) Finalize      | Otherwise → `avgRefMis[fid] = sums[fid] / counts[fid]`         | All Class 1 fixtures with non-empty features     |

## Test inventory

| Test case         | Status      | Notes            |
|--------|-------------|--|
| `ComputeFeatureReferenceMisorientationsFilter: Class 1 - Mode 0 SingleGrainIdentity`  | new-for-V&V | 2×2×2; single feature; all identity quats; expected FRM=0, avg=0. Class 4 invariants asserted.         |
| `ComputeFeatureReferenceMisorientationsFilter: Class 1 - Mode 0 KnownAngle5deg`       | new-for-V&V | 2×2×2; single feature; all voxel quats = 5° about z, AvgQuats[1] = identity; expected FRM=5° everywhere, avg=5°. Verifies non-zero misorientation magnitude. Class 4 invariants asserted.      |
| `ComputeFeatureReferenceMisorientationsFilter: Class 1 - Mode 0 MultiGrain EdgeCases` | new-for-V&V | 4×3×1; 5 features (background + 4 grains); covers skip paths (background, mid-feature unphased), zero-count finalize path (feature 4 all-unphased → avg=0), and normal accumulation. Class 4 invariants asserted.             |
| `ComputeFeatureReferenceMisorientationsFilter: Class 1 - Mode 1 KnownCenter`          | new-for-V&V | 3×3×1; single feature; center voxel hand-picked via unique max `GBEuclideanDistances`; verifies `centers[]` selection + `EuclideanCenters` coord writing. Class 4 invariants asserted.          |
| `ComputeFeatureReferenceMisorientationsFilter: Class 1 - Mode 1 MultiGrain CenterIsolation`            | new-for-V&V | 2×3×1; 2 features; verifies `centers[fid]` isolation per feature + tied-distance `>=` later-voxel-wins tie-break. Class 4 invariants asserted.                |
| `ComputeFeatureReferenceMisorientationsFilter: Class 1 - Mode 1 3D Volume`            | new-for-V&V | 3×3×2; single feature; verifies linear `voxelIdx → (x,y,z)` arithmetic when `dimZ > 1`. Class 4 invariants asserted.          |
| `ComputeFeatureReferenceMisorientationsFilter: Class 4 - Invariants Sweep`            | new-for-V&V | Runs Mode 0 and Mode 1 configurations distinct from the value-specific fixtures; asserts only the Class 4 invariants. Catches future regressions where specific values shift but invariants still hold.       |
| `ComputeFeatureReferenceMisorientationsFilter: SIMPL Backwards Compatibility`         | retained    | `DYNAMIC_SECTION` over SIMPL 6.4 + 6.5 conversion fixtures (`test/simpl_conversion/6_*/ComputeFeatureReferenceMisorientationsFilter.json`); validates UUID + argument-key + parameter-value decoding.            |
| *(retired)* `ComputeFeatureReferenceMisorientationsFilter_AverageMisorientation`      | retired     | Removed 2026-06-01. Regression-against-archive test consuming `compute_feature_reference_misorientation.tar.gz` exemplar arrays; archive's exemplar values were a circular oracle (regenerated from pre-EbsdLib-2.4.1 SIMPLNX output). Test failure surfaced when EbsdLib 2.4.1 CubicOps precision fix shifted exemplar values by 2× epsilon. Replaced by inline Class 1 + Class 4 fixtures above. |
| *(retired)* `ComputeFeatureReferenceMisorientationsFilter_EuclideanDistance`          | retired     | Same as above for Mode 1; archive exemplar shifted by 10× epsilon post-EbsdLib-2.4.1.|

All 8 active TEST_CASEs pass at the verified commit (`100% tests passed, 0 tests failed out of 8` in ~0.7s).

## Exemplar archive

**None — data inlined in `test/ComputeFeatureReferenceMisorientationsTest.cpp` namespace `AnalyticalFixtures`.**

The pre-existing `compute_feature_reference_misorientation.tar.gz` archive was retired during this V&V cycle. See `src/Plugins/OrientationAnalysis/vv/provenance/ComputeFeatureReferenceMisorientationsFilter.md` for the retirement rationale and the methodology used to construct the replacement data fixtures.

- **Archive:** None
- **SHA512:** N/A
- **Provenance:** `src/Plugins/OrientationAnalysis/vv/provenance/ComputeFeatureReferenceMisorientationsFilter.md`

## Deviations from DREAM3D 6.5.171

Two active deviations are documented: D1 (orientation-library precision) and D2 (a legacy Mode 1 crash).

### ComputeFeatureReferenceMisorientationsFilter-D1

- **Symptom:** Per-cell and per-feature misorientations differ at the precision level. On the analytical fixtures the baseline difference reaches 1.86e-4°; on Small IN100 it reaches 0.0738° per cell and 0.01761° per feature average. A local legacy build with the shared orientation-precision correction reduces the per-cell residual to at most one float32 ULP.
- **Root cause:** Precision — propagation of the EbsdLib 2.4.1 `CubicOps::calculateMisorientationInternal` precision fix. This filter is a clean Port of the legacy algorithm; no algorithmic deviation. See `vv/deviations/ComputeFeatureReferenceMisorientationsFilter.md` for full root-cause walkthrough.

### ComputeFeatureReferenceMisorientationsFilter-D2

- **Symptom:** DREAM3D 6.5.171 terminates with a segmentation fault whenever Mode 1 (Euclidean-distance reference) is selected; SIMPLNX completes and returns the analytical result.
- **Root cause:** **Bug** in DREAM3D 6.5.171. The legacy execution path reads the Mode-0-only `AvgQuats` pointer to obtain the feature count before entering the Mode 1 center-selection path. A surgical local correction uses the tuple count of the always-created feature-average output instead; the corrected build matches SIMPLNX bit-for-bit on all three analytical Mode 1 fixtures.
- See `vv/deviations/ComputeFeatureReferenceMisorientationsFilter.md` for the debugger evidence and patch proof.
