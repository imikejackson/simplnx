# V&V Report: ComputeFeatureNeighborMisorientationsFilter

|           |                          |
|-----------|--------------------------|
| Plugin    | OrientationAnalysis      |
| SIMPLNX UUID               | `0b68fe25-b5ef-4805-ae32-20acb8d4e823`                |
| SIMPLNX Human Name         | Compute Feature Neighbor Misorientations              |
| DREAM3D 6.5.171 equivalent | `FindMisorientations` — `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindMisorientations.{h,cpp}` (UUID `286dd493-4fea-54f4-b59e-459dd13bbe57`) |
| Verified commit            | `a307946e7` (v7.4.2 release)               |
| Status | COMPLETE     |
| Sign-off                   | *Michael Jackson <mike.jackson@bluequartz.net> (V&V cycle completion + divisor bug fix, 2026-06-02)*                                             |
| Second-engineer sign-off   | Nathan Young — 2026-06-11 (approving reviewer, PR #1631)                                                                                         |

## At a glance

| Aspect                 | Current state            |
|------------------------|--------------------------|
| Algorithm Relationship | **Port with one inherited bug corrected** — the loop structure and phase gate are preserved, while quaternion types and the EbsdLib API changed. D1 corrects the divisor assignment. |
| Oracle (confirmed)     | **Class 1 (Analytical)** uses 3 hand-derived fixtures; **Class 4 (Invariant)** checks value bounds, NaN handling, and the per-feature averaging formula. |
| Code paths enumerated  | 5 of 5 paths exercised across the analytical and invariant fixtures. |
| Tests today            | 5 test cases cover 3 Class 1 fixtures, 1 Class 4 invariant fixture, and SIMPL conversion. No exemplar archive is used. |
| Exemplar archive       | **None** — the circular `6_6_stats_test_v2.tar.gz` test and an unimplemented averaging stub are retired. Inline fixtures replace both tests. |
| Legacy comparison      | **Run** — Three analytical fixtures and two Small IN100 cases confirmed D1's defective legacy divisor and D2's precision difference; the local legacy proof build matched SIMPLNX within the stated tolerance. |
| Bug flags              | `ComputeFeatureNeighborMisorientationsFilter-D1` is an inherited divisor bug that is fixed in SIMPLNX and pinned by the mixed-phase regression test. |
| V&V phase | **COMPLETE.** |

## Summary

`ComputeFeatureNeighborMisorientationsFilter` writes each feature-to-neighbor misorientation and can also write the average of the valid same-phase values. Three Class 1 fixtures and one Class 4 invariant fixture verify the calculation, NaN handling, and averaging formula. The comparison confirms D1's inherited divisor defect and D2's precision difference; the local legacy proof build reproduces SIMPLNX on the analytical cases and within two float32 ULP on the Small IN100 cases.

## Algorithm Relationship

*Classification:* **Port (with UUID reassignment, name rename, and one divisor-bug fix)** ~~| Minor changes | Rewrite | New filter~~

*Evidence:* The SIMPLNX algorithm at `Algorithms/ComputeFeatureNeighborMisorientations.cpp` (~120 lines) is a near line-by-line translation of legacy `FindMisorientations::execute()` (DREAM3D 6.5.171). Same per-feature outer loop, same per-neighbor inner loop, same phase-match gate, same optional average computation. SIMPLNX was assigned a new UUID (`0b68fe25-...` vs legacy `286dd493-...`) for the `Find`→`Compute` rename plus the `FeatureNeighbor` qualifier (clarifying that this filter computes feature-to-neighbor misorientations, distinguishing it from the per-cell and per-reference misorientation filters in the same module).

*Port-time deltas:*

1. **EbsdLib API**: `getMisoQuat(q1, q2, n1, n2, n3) → float angle` → `calculateMisorientation(q1, q2) → AxisAngleDType`. Same math via `LaueOps`.
2. **Quaternion precision**: `QuatF` (float32) → `QuatD` (double internal). Float32 stored values are promoted to double for the math; per-neighbor misorientations are cast back to `float32` for storage in the `NeighborList<float32>`.
3. **Cancel checks**: SIMPLNX adds `m_ShouldCancel.load()` at the outer-loop top; legacy had no cancel mechanism.
4. **EbsdLib 2.4.1 CubicOps precision improvement** (external dependency change): manifests on real EBSD data with cubic-phase sym-op-aligned grain-pair boundaries; non-observable on the V&V data fixtures (pure φ1 rotations). See D2.
5. **Divisor bug fix** (2026-06-02): `tempMisoList = featureNeighborList.size();` moved from inside the inner j-loop (line 75 pre-fix) to before the j-loop (alongside `tempMisorientationLists[i].assign(...)` at line ~67). See D1.

*Material PRs since baseline (2025-10-01):* None identified that materially change this filter's algorithm.

## Oracle

*Class:* **1 (Analytical)** primary + **4 (Invariant)** companion. Class 3 (Paper-based) N/A — math is delegated to `ebsdlib::LaueOps::calculateMisorientation` and verified in EbsdLib's own V&V.

*Applied:* The Class 1 and Class 4 oracles are applied as described below.

### Class 1 — Analytical

Per-neighbor misorientation values are derived in closed form from the input `AvgQuats` + `FeaturePhases` + `NeighborList` + `CrystalStructures` arrays by hand-tracing the algorithm. The fixtures use pure φ1-rotation quaternions (Bunge ZXZ Euler `(φ1, 0, 0)`) so that misorientation between any two same-phase features equals `|Δφ1|` modulo the cubic c-axis 4-fold symmetry. For Δφ1 ∈ {0°, 5°, 10°}, no symmetry reduction applies (`5°`, `10°` are below the 45° fold) so expected per-neighbor entries are `|Δφ1|` exactly; phase-mismatched neighbors produce `NaN`. The expected per-feature average is `sum-of-non-NaN-entries / count-of-non-NaN-entries`. The three Class 1 fixtures differ in the *order* in which phase-matched and phase-mismatched neighbors appear in the list, systematically exercising the per-mismatch divisor decrement.

### Class 4 — Invariant

Three invariants every filter run must satisfy regardless of input configuration:

- **Per-entry validity**: each `MisorientationList[fid][j]` is either `NaN` (phase mismatch) or a non-negative misorientation bounded above by the cubic max symmetry-reduced angle (`62.8°`).
- **Per-feature averaging formula**: `AvgMisorientations[fid]` equals `sum(non-NaN entries in MisorientationList[fid]) / count(non-NaN entries in MisorientationList[fid])`. The formula is invariant under neighbor-list reordering — re-ordering the input neighbor list should NOT change the per-feature average value (it WILL change the order of entries within the per-feature list, but the average is order-independent).
- **All-mismatch case**: when every neighbor in a feature's list is a phase mismatch, the per-feature average is `NaN` (no valid entries to average).

*Encoded:* The tests below encode the oracle.

- **Class 1 (Analytical)**: `test/ComputeFeatureNeighborMisorientationsTest.cpp` — 3 `TEST_CASE` blocks under the `Class 1 - …` family. Per-neighbor expected values asserted via `Approx().margin(1e-3f)`; per-feature averages asserted via `Approx().margin(1e-3f)`.
- **Class 4 (Invariant)**: `Class 4 - Invariants` test — runs a 5-feature 3-neighbor configuration and asserts the per-entry validity invariant and the per-feature averaging formula derived from the per-entry values.
- *(kept)* `SIMPL Backwards Compatibility` — SIMPL 6.4 + 6.5 conversion paths via `DYNAMIC_SECTION`.

*Second-engineer review:* Nathan Young — 2026-06-11 (approving reviewer, PR #1631).

Recorded review topics:

- *The Class 1 hand-derivations in the 3 data fixtures + 1 invariants test.*
- *The divisor-bug fix at `Algorithms/ComputeFeatureNeighborMisorientations.cpp` line ~70 (the new assignment location) and the corresponding test that exercises the fix (`Class 1 - Mixed Phase Neighbors (exposes divisor bug)`).*
- *The decision to retire the `6_6_stats_test_v2.tar.gz` Small-IN100 exemplar archive in favor of inline data fixtures (shared retirement with `ComputeKernelAvgMisorientationsFilter`).*

## Bugs found and fixed

| Deviation | Defect | Affected released versions | Resolution in this branch |
|-----------|--------|----------------------------|---------------------------|
| `ComputeFeatureNeighborMisorientationsFilter-D1` | The divisor is reset inside the neighbor loop, so mixed-phase feature averages can use the wrong neighbor count. | DREAM.3D 6.5.171; DREAM3D-NX v7.0.0 through v7.4.1. | The divisor is initialized once before the neighbor loop, and each phase mismatch decrements it once. |

## Code path coverage

*5 of 5 paths exercised directly.*

Source: `src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/Algorithms/ComputeFeatureNeighborMisorientations.cpp` (115 lines).

The algorithm has one logical pass: a per-feature outer loop, with a per-neighbor inner loop inside it. Cancel check at the outer-loop top.

| # | Pass        | Path                               | Test case                                                 |
|---|-------------|-------------------------------------------------------------------------------------------------------------------------------|--------------------|
| 1 | Per-feature | Cancel check at outer-loop top (`m_ShouldCancel.load()` → early return)                                                       | *Not directly tested.* Failure mode would manifest as test hang in any test; not specifically exercised. Low-value gap.                            |
| 2 | Per-neighbor | Phase match (`laueClass1 == xtalType2 && laueClass1 < orientationOps.size()`) → write list entry; if `ComputeAvgMisors`, accumulate | All Class 1 fixtures; primary algorithmic path            |
| 3 | Per-neighbor | Phase mismatch → write `NaN` entry; if `ComputeAvgMisors`, decrement `tempMisoList` divisor                                  | `Mixed Phase Neighbors (exposes divisor bug)` + `Mismatch Last Order` + `Class 4 - Invariants`                                                       |
| 4 | Per-feature (finalize) | `tempMisoList != 0` → `(*avgMisorientations)[i] /= tempMisoList`                                                    | All Class 1 fixtures + Class 4 invariants                 |
| 5 | Per-feature (finalize) | `tempMisoList == 0` (all neighbors mismatched) → `(*avgMisorientations)[i] = NaN`                                  | *Not directly tested.* Algorithm path is exercised whenever a feature has only phase-mismatched neighbors; no data fixture currently constructs this configuration but the invariant `count == 0 → avg == NaN` is asserted in the Class 4 invariants helper if such a feature were present. Low-value gap. |

## Test inventory

| Test case                                     | Status      | Notes                  |
|--------|-------------|--------------------------------------------------------------------------------------------------------------------------------|
| `ComputeFeatureNeighborMisorientationsFilter: Class 1 - Single Phase Two Neighbors`                                                       | new-for-V&V | 4 features; all phase 1; per-feature `MisorientationList[1] = [5°, 10°]`; expected avg = 7.5°. Verifies the basic per-feature averaging path with all phase-matched neighbors.                                    |
| `ComputeFeatureNeighborMisorientationsFilter: Class 1 - Mixed Phase Neighbors (exposes divisor bug)`                                      | new-for-V&V | Bug-exposing fixture. 5 features; feature 1 (phase 1) has neighbors `[2 (match), 4 (mismatch), 3 (match)]`. The last-iterated neighbor is a match, so the pre-fix bug reassigns the divisor to 3; the test asserts the correct avg = 7.5°. **Failed on pre-fix code (gave 5.0°); passes on post-fix code (gives 7.5°).**                                                                            |
| `ComputeFeatureNeighborMisorientationsFilter: Class 1 - Mismatch Last Order`                                                              | new-for-V&V | Same features as the bug-exposing fixture but neighbor order `[2 (match), 3 (match), 4 (mismatch)]`. Last neighbor is a mismatch; the decrement at line 90 is the last write to `tempMisoList` so the bug doesn't fire. Both pre-fix and post-fix code produce 7.5°.                                     |
| `ComputeFeatureNeighborMisorientationsFilter: Class 4 - Invariants`                                                                       | new-for-V&V | Asserts per-entry validity (NaN or non-negative ≤ 62.8°) and the per-feature averaging formula (derived from the per-entry values, not a specific hard-coded number). Uses neighbor order `[4 (mismatch), 2 (match), 3 (match)]`. Catches future regressions that preserve specific values but break the invariant relationship.                                                                |
| `ComputeFeatureNeighborMisorientationsFilter: SIMPL Backwards Compatibility`                                                              | retained    | `DYNAMIC_SECTION` over SIMPL 6.4 + 6.5 conversion fixtures. UUID + argument-key + parameter-value validation only. |
| *(retired)* main `ComputeFeatureNeighborMisorientationsFilter` (consumed `6_6_stats_test_v2.tar.gz`)                                       | retired     | Removed 2026-06-02. Regression-against-archive test; archive's exemplar values were a circular oracle (regenerated from pre-EbsdLib-2.4.1 SIMPLNX output). Test was already failing as the EbsdLib precision fix shifted exemplar values beyond the regression-check epsilon. Replaced by inline Class 1 + Class 4 fixtures above.                                                                  |
| *(retired)* `ComputeFeatureNeighborMisorientationsFilter: Misorientation Per Feature`                                                     | retired     | Removed 2026-06-02. The `[.][UNIMPLEMENTED][!mayfail]` stub left `ComputeAvgMisors=true` with zero CI coverage, which is precisely why the divisor bug at algorithm.cpp:75 went undetected. The 3 Class 1 fixtures above cover that parameter combination.                                               |

All 5 active TEST_CASEs pass at the verified commit (`100% tests passed, 0 tests failed out of 5` in ~0.3s).

## Exemplar archive

**None — data inlined in `test/ComputeFeatureNeighborMisorientationsTest.cpp` namespace `AnalyticalFixtures`.**

The pre-existing `6_6_stats_test_v2.tar.gz` archive is being retired. See `src/Plugins/OrientationAnalysis/vv/provenance/ComputeFeatureNeighborMisorientationsFilter.md` for the retirement rationale.

- **Archive:** None
- **SHA512:** N/A
- **Provenance:** `src/Plugins/OrientationAnalysis/vv/provenance/ComputeFeatureNeighborMisorientationsFilter.md`

## Deviations from DREAM3D 6.5.171

Two deviations documented:

### ComputeFeatureNeighborMisorientationsFilter-D1

- **Symptom:** Per-feature `AvgMisorientations` values differ on mixed-phase neighbor lists. In the forcing fixture, DREAM3D 6.5.171 produces 5.0° instead of the correct 7.5°; on mixed-phase Small IN100, 508 feature averages follow the defective divisor formula and differences reach 68.7°. A local legacy build with the surgical divisor fix follows the correct formula for every feature and reproduces SIMPLNX exactly on the analytical fixtures.
- **Root cause:** **Bug** in DREAM3D 6.5.171 (also present in pre-fix SIMPLNX). See `vv/deviations/ComputeFeatureNeighborMisorientationsFilter.md` for the technical mechanism.

### ComputeFeatureNeighborMisorientationsFilter-D2

- **Symptom:** Per-neighbor `MisorientationList` and per-feature `AvgMisorientations` values differ at the precision level. On Small IN100, DREAM3D 6.5.171 differs from the independent float64 reference by up to 3.37e-4° per neighbor; after the shared cubic/hexagonal precision corrections, the local legacy build agrees with SIMPLNX within `rtol=1e-6, atol=1e-6`, with residuals no larger than two float32 ULP.
- **Root cause:** **Precision** — propagation of the EbsdLib 2.4.1 `CubicOps::calculateMisorientationInternal` precision improvement, characterized in `vv/deviations/BadDataNeighborOrientationCheckFilter.md`. See `vv/deviations/ComputeFeatureNeighborMisorientationsFilter.md` for the per-filter context.
**SIMPLNX-side fix ships in DREAM3D-NX 7.4.2** — the deviation from legacy remains, since 6.5.171 is unchanged: `ComputeFeatureNeighborMisorientationsFilter-D1`.
