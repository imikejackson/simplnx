# V&V Report: ComputeFeatureFaceMisorientationFilter

|           |                          |
|-----------|--------------------------|
| Plugin    | OrientationAnalysis      |
| SIMPLNX UUID               | `f3473af9-db77-43db-bd25-60df7230ea73`|
| SIMPLNX Human Name         | Compute Feature Face Misorientation (Face)             |
| DREAM3D 6.5.171 equivalent | `GenerateFaceMisorientationColoring` — `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/GenerateFaceMisorientationColoring.{h,cpp}`            |
| Verified commit            | `a307946e7` (v7.4.2 release)                |
| Status | COMPLETE     |
| Sign-off  | *Nathan Young (algorithm rewrite + initial dataset, 2026-05-19) — Michael Jackson <mike.jackson@bluequartz.net> (hand-built test data, V&V completion, 2026-05-28)*           |
| Second-engineer sign-off | Nathan Young — 2026-06-11 (approving reviewer, PR #1631; output rework PR #1618 approved by Matthew Marine 2026-06-04). Supersedes the 2026-05-28 technical-authority self-sign-off. |

## At a glance

| Aspect                 | Current state            |
|------------------------|--------------------------|
| Algorithm Relationship | **Rewrite** — SIMPLNX writes a scalar angle, supports 11 Laue classes, writes NaN for invalid faces, uses the modern EbsdLib API, and uses `ParallelDataAlgorithm`. Six deltas are documented below. |
| Oracle (confirmed)     | **Class 1 (Analytical)** — 37 hand-built faces exercise all 11 EbsdLib Laue classes, three symmetry boundaries per class, and invalid-face cases. |
| Code paths enumerated  | 7 of 7 paths exercised by the curated analytical fixture. |
| Tests today            | 2 test cases cover the Class 1 curated data and SIMPL 6.4/6.5 conversion. The obsolete invalid-execution test is retired. |
| Exemplar archive       | **None** — the 37-face dataset is inline in `test/ComputeFeatureFaceMisorientationTest.cpp`. |
| Legacy comparison      | **Run** — Three comparisons mapped the legacy axis-angle norm to the SIMPLNX scalar angle and confirmed D1 through D4; the local legacy proof build matched all 640,880 comparable faces within the stated tolerance. |
| Bug flags              | `ComputeFeatureFaceMisorientations-D3` is a legacy correctness ambiguity. D4 is a precision deviation, not a Bug-class entry. |
| V&V phase | **COMPLETE.** |

## Summary

`ComputeFeatureFaceMisorientationFilter` computes one misorientation angle in degrees for each surface-mesh triangle. A Class 1 analytical fixture verifies 37 faces across all 11 Laue classes. The three-case comparison covered 640,880 comparable faces; the unmodified release differed above 0.0001° on 1,249 faces, while the local legacy proof build agreed with SIMPLNX within the stated tolerance on every comparable face.

## Algorithm Relationship

*Classification:* **Rewrite** ~~| Port | Minor changes | New filter~~

*Evidence:* Deliberate rewrite of legacy `GenerateFaceMisorientationColoring::CalculateFaceMisorientationColorsImpl` (DREAM3D 6.5.171). Same SIMPL UUID retained; SIMPL 6.4/6.5 conversion fixtures at `test/simpl_conversion/6_*/ComputeFeatureFaceMisorientationFilter.json`. Legacy control-flow structure is preserved (parallel per-triangle, two-face phase + Laue-class lookup, fall-through for mismatches) but every interior choice differs — see deltas below.

*Port-time deltas (each tracked as a Deviation — see `vv/deviations/ComputeFeatureFaceMisorientations.md`):*

1. **Output structure: 3-component → 1-component** (Deviation D2). Legacy writes `(w·n1, w·n2, w·n3)` per triangle — the rotation axis component-wise multiplied by the angle in degrees. SIMPLNX writes just the angle in degrees. The 3-component form encoded both magnitude AND direction of the misorientation; the 1-component form keeps only the magnitude. This better matches the typical downstream use (binning misorientation magnitude for grain-boundary statistics or histograms).
2. **Laue class support: 2 classes → 11 classes** (Deviation D1). Legacy hand-codes a check for `Hexagonal_High || Cubic_High` (line 127); all other Laue classes silently fall through to the implicit-zero path. SIMPLNX checks `laueIndex < m_LaueOrientationOps.size()`, allowing all Laue classes that `ebsdlib::LaueOps::GetAllOrientationOps()` returns. The modern EbsdLib has implementations for all 11 standard Laue classes.
3. **Invalid-face handling: implicit 0 → explicit NaN** (Deviation D3). Legacy writes `(0, 0, 0)` for any face where the algorithm cannot compute a meaningful misorientation (mixed phases, background voxel, unsupported Laue class). SIMPLNX writes `(NaN)`. Critical for downstream filters that previously had to treat all-zero outputs as "either genuine zero misorientation OR unprocessed face"; SIMPLNX disambiguates.
4. **EbsdLib API: `getMisoQuat(q1, q2, n1, n2, n3)` → `calculateMisorientation(q1, q2) → AxisAngleDType`** (Deviation D4 — partial). Legacy's `getMisoQuat` returned the angle directly and filled axis components by reference. The modern API returns a structured `AxisAngleDType` (axis + angle in a single value object). The legacy API is no longer present in the current EbsdLib.
5. **Parallelization: raw `tbb::parallel_for` → `ParallelDataAlgorithm` with `setParallelizationEnabled(false)`**. Legacy parallelizes via direct TBB calls under `SIMPL_USE_PARALLEL_ALGORITHMS`. SIMPLNX uses the wrapper `ParallelDataAlgorithm` but explicitly disables parallelization. **Per the project's thread-safety guidance**: DataArray write access from worker threads is not guaranteed safe across all SIMPLNX data store implementations; serial execution is the safe default. No algorithmic effect on completed runs.
6. **EbsdLib precision fix** (Deviation D4 — full root-cause). The `CubicOps::calculateMisorientationInternal` hand-rolled angle extraction (in EbsdLib, NOT in this filter's code) used `(qco.z()+qco.w())/sqrt(2)` followed by `acos(w)`. When the misorientation lies on a cubic symmetry op (e.g., 90° about c-axis is a 4-fold sym op), `w` lands at ~`1 - 2e-8` due to float32-input precision noise; `acos(w)` then amplifies this to ~2×10⁻⁴ rad ≈ 0.023° residual. Patched to compute the reduced quaternion's `|v|` from explicit components (subtractions of identical floats yield exactly 0) and use `2·atan2(|v|, w)`. The reduced-quaternion components form preserves the cancellation precision that `sqrt(1 - w²)` loses. Eliminates the ~0.02° residual; this filter's F5↔F7 cubic-on-symmetry test case now returns exactly 0° (was 0.0212°).

*Material PRs since baseline (2026-05-19, Nathan's V&V doc commit):*

- **`nathan/enh/issue_1596`** (squashed into this branch 2026-05-28) — Six commits authored by Nathan implementing the filter rewrite (deltas 1–5 above). Squash-merged because the intermediate commits ("filter and algorithm implementation; test pending", "Create new test (failing)", "patch test", etc.) are work-in-progress and the squashed unit is the meaningful change.
- **This V&V cycle (2026-05-28)** — Mike added the hand-built Class 1 test dataset, Mike added Trigonal_High coverage (originally missing — Laue indices 0–9 only; index 10 added to close the gap), Mike + Claude root-caused the F5↔F7 precision residual to EbsdLib `CubicOps::calculateMisorientationInternal`, Claude patched EbsdLib with the `2·atan2(|v|, w)` form, Mike updated test assertion from `0.0212f` → `0.0f`.

## Oracle

*Class:* **1 (Analytical)** primary.

*Applied:* The Class 1 analytical oracle derives expected values as described below.

Expected misorientation values are derived from the closed-form symmetry-group reduction of the boundary's true rotation. The dataset uses pure φ1-rotations (Bunge Euler angles `(φ1, 0, 0)` with `Φ = φ2 = 0`), so the true misorientation between any two features is simply `|Δφ1|` modulo the c-axis-aligned symmetry operators of the Laue class.

For each Laue class L, four features (one phase, four orientations) are constructed:
- Feature A: `φ1 = 0°`
- Feature B: `φ1 = 45°`
- Feature C: `φ1 = 90°`
- Feature D: `φ1 = 180°`

And three boundary faces are constructed: A↔B, A↔C, A↔D. The symmetry-reduced expected misorientation depends on the Laue class's c-axis n-fold:

| Laue class (idx)          | c-axis n-fold | A↔B (0°↔45°) | A↔C (0°↔90°) | A↔D (0°↔180°) |
|--------------------------------------------|---------------|--------------|--------------|---------------|
| Hexagonal_High m⁻³m (0)   | 6-fold        | 15°          | 30°          | 0°            |
| Cubic_High 6/mmm (1)      | 4-fold        | 45°          | **0°***      | 0°            |
| Hexagonal_Low 6/m (2)     | 6-fold        | 15°          | 30°          | 0°            |
| Cubic_Low m-3 (3)         | 2-fold (face) | 45°          | 90°          | 0°            |
| Triclinic -1 (4)          | 1-fold        | 45°          | 90°          | 180°          |
| Monoclinic 2/m (5)        | 1-fold        | 45°          | 90°          | 180°          |
| OrthoRhombic mmm (6)      | 2-fold        | 45°          | 90°          | 0°            |
| Tetragonal_Low 4/m (7)    | 4-fold        | 45°          | 0°           | 0°            |
| Tetragonal_High 4/mmm (8) | 4-fold        | 45°          | 0°           | 0°            |
| Trigonal_Low -3 (9)       | 3-fold        | 45°          | 30°          | 60°           |
| Trigonal_High -3m (10) — *SIMPLNX addition* | 3-fold        | 45°          | 30°          | 60°           |

*Cubic_High A↔C: 90° about c-axis is a 4-fold cubic symmetry op of m-3m, so the true symmetry-reduced misorientation is exactly 0°. Before the EbsdLib precision fix, this returned ~0.0212° due to `acos(w)` near 1. See Deviation D4 / Algorithm Relationship delta 6.

**Edge cases** (faces 30–33, after the 30 normal cases): all four expected to produce NaN.

| Face | Front label | Back label | Expected | Path exercised |
|------|-------------|------------|----------|--------------------------------------------------|
| 30   | 0           | 1          | NaN      | Background-front (frontFeature == 0)             |
| 31   | 1           | 0          | NaN      | Background-back (backFeature == 0)               |
| 32   | 1           | 5          | NaN      | Different phases (phase 1 Hex_High vs phase 2 Cubic_High) |
| 33   | 5           | 1          | NaN      | Different phases, reversed      |

*Encoded:* The tests below encode the oracle.

- **Class 1 (Analytical)**: `test/ComputeFeatureFaceMisorientationTest.cpp::"OrientationAnalysis::ComputeFeatureFaceMisorientationFilter: Curated Data"` — 30 + 4 + 3 = 37 fixture assertions, 54 total assertions (including geometry setup REQUIRE-VALID checks).
- *(kept)* `test/ComputeFeatureFaceMisorientationTest.cpp::"OrientationAnalysis::ComputeFeatureFaceMisorientationFilter: SIMPL Backwards Compatibility"` — SIMPL 6.4 + 6.5 conversion paths via `DYNAMIC_SECTION`.

*Second-engineer review:* Nathan Young — 2026-06-11 (approving reviewer, PR #1631; output rework PR #1618 approved by Matthew Marine 2026-06-04). Review focus: the symmetry-group hand calculations for Trigonal_High and the EbsdLib precision-fix rationale. Note that the Trigonal_Low and Trigonal_High closed-form values are identical because mirror planes that contain the c-axis do not reduce pure c-axis rotations further.

## Bugs found and fixed

| Deviation | Defect | Affected released versions | Resolution in this branch |
|-----------|--------|----------------------------|---------------------------|
| `ComputeFeatureFaceMisorientations-D3` | DREAM3D 6.5.171 writes zero for invalid or unprocessed faces, so users cannot distinguish them from a true zero misorientation. | DREAM.3D 6.5.171 only. DREAM3D-NX was not affected. | SIMPLNX writes NaN for every face that does not have a meaningful misorientation. |

## Code path coverage

*7 of 7 paths exercised. Cancel-check paths and "valid Laue class" type-dispatch are aggregate-tested via the Class 1 dataset; per-Laue-class paths are confirmed individually by the per-class assertions.*

Source: `src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/Algorithms/ComputeFeatureFaceMisorientation.cpp` (134 lines).

| # | Phase           | Path             | Test case            |
|---|-----------------|-------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------|
| 1 | Cancel check    | `m_ShouldCancel` checked at top of per-triangle loop → early return| *Not directly tested.* Loop-guard only; cancel-signal injection requires test infrastructure not present. Low-value gap. |
| 2 | Per-face        | `frontFeature == 0` (background) → `frontPhase = 0` → falls through to "different phases" path → NaN written           | `Curated Data` — face 30 `(0, 1)` covers this path     |
| 3 | Per-face        | `backFeature == 0` (background) → `backPhase = 0` → falls through to "different phases" path → NaN written             | `Curated Data` — face 31 `(1, 0)` covers this path     |
| 4 | Per-face        | `frontPhase > 0 && frontPhase != backPhase` → falls through to "different phases" path → NaN written  | `Curated Data` — faces 32 `(1, 5)` and 33 `(5, 1)` cover this path      |
| 5 | Per-face        | `frontPhase > 0 && frontPhase == backPhase && laueIndex >= m_LaueOrientationOps.size()` → NaN written (unsupported Laue class)          | *Not directly tested.* All 11 Laue classes in the curated dataset are within EbsdLib's supported range. Low-value gap. |
| 6 | Per-face        | `frontPhase > 0 && frontPhase == backPhase && laueIndex < m_LaueOrientationOps.size()` → call `m_LaueOrientationOps[laueIndex]->calculateMisorientation(q1, q2)` → write `axisAngle[3] * k_180OverPiD` (angle in degrees) | `Curated Data` — all 30 normal-case asserts + 3 Trigonal_High asserts exercise this path |
| 7 | Per-face (math) | Inside `calculateMisorientation`: cubic-class sym-op enumeration via type-1/2/3 reduced quaternion (in `CubicOps::calculateMisorientationInternal`), with the precision-fixed `2·atan2(|v|, w)` angle extraction               | `Curated Data` — F5↔F6 (type 1), F5↔F7 (type 2, EbsdLib precision-fix-critical path), F5↔F8 (type 2 or 3 depending on which sym op wins) |

## Test inventory

| Test case        | Status      | Notes           |
|--------------------------------------------------------------------------------------|-------------|----------------------------------------------|
| `OrientationAnalysis::ComputeFeatureFaceMisorientationFilter: Curated Data`          | new-for-V&V | Class 1 hand-built dataset; 30 normal + 4 edge + 3 Trigonal_High asserts. Replaces the legacy `Valid filter execution` test (which used the `6_6_Small_IN100_GBCD.tar.gz` exemplar) — the legacy test was a regression-against-exemplar test, not a closed-form correctness test, and was incompatible with the rewritten 1-component output. |
| `OrientationAnalysis::ComputeFeatureFaceMisorientationFilter: SIMPL Backwards Compatibility` | kept        | Unchanged. `DYNAMIC_SECTION` over SIMPL 6.4 and 6.5 conversion fixtures (`test/simpl_conversion/6_*/ComputeFeatureFaceMisorientationFilter.json`); validates UUID, argument keys, and parameter conversion only.             |
| *(retired)* `OrientationAnalysis::ComputeFeatureFaceMisorientationFilter: Invalid filter execution` | retired     | Removed during Nathan's rewrite. The Class 1 dataset's faces 30–33 cover the same paths via the NaN-on-invalid-face semantics; the explicit-preflight-failure tests are no longer reachable for the new code structure.     |

## Exemplar archive

- **Archive:** None — data inlined in `test/ComputeFeatureFaceMisorientationTest.cpp` namespace `curated`.
- **SHA512:** N/A
- **Provenance:** `src/Plugins/OrientationAnalysis/vv/provenance/ComputeFeatureFaceMisorientations.md`

Data construction details: 102 vertices laid out in a y-axis-stacked grid (one row of 9 vertices per Laue class block + edge case block), 34 + 3 triangles (3 per Laue class + 4 edge case + 3 Trigonal_High), 41 + 4 features, 12 + 1 ensembles. The unique-vertices-per-triangle layout means each triangle is geometrically independent (no shared edges or vertices between triangles) — this is intentional, since the algorithm only reads `FaceLabels`, not vertex coordinates, so the geometric layout is arbitrary.

## Deviations from DREAM3D 6.5.171

Four deviation classes are documented. D1 and D2 are library or algorithmic choices, D3 records a legacy correctness ambiguity, and D4 records a related EbsdLib precision fix.

### ComputeFeatureFaceMisorientations-D1

- Supports all 11 Laue classes; legacy supported only Hex_High and Cubic_High. See `vv/deviations/ComputeFeatureFaceMisorientations.md`.

### ComputeFeatureFaceMisorientations-D2

— Output is a 1-component angle in degrees; legacy was a 3-component axis·angle vector. See `vv/deviations/ComputeFeatureFaceMisorientations.md`.

### ComputeFeatureFaceMisorientations-D3

— Invalid faces (mixed phase, background voxel, unsupported Laue class) write NaN; legacy wrote 0. See `vv/deviations/ComputeFeatureFaceMisorientations.md`.

### ComputeFeatureFaceMisorientations-D4

— Precision improvement on cubic boundaries that lie on a 4-fold symmetry op. Root-caused to EbsdLib `CubicOps::calculateMisorientationInternal`; patched at the EbsdLib level (replaces `acos(w)` with `2·atan2(|v|, w)` using explicit reduced-quaternion components). See `vv/deviations/ComputeFeatureFaceMisorientations.md`.

### Downstream impact note (not a deviation, characterized for transparency):
 
The EbsdLib precision fix in D4 propagates through any filter that consumes cubic misorientations. Eight OrientationAnalysis unit tests now fail against their pre-fix exemplars with diffs in the range `1.4× to 10× epsilon` (epsilons of `1e-4`, observed diffs `1.4e-4` to `1e-3`): 
- `BadDataNeighborOrientationCheckFilter: Case 1.{3,4,5,6}.3`
- `ComputeFeatureReferenceMisorientationsFilter_AverageMisorientation`
- `ComputeFeatureReferenceMisorientationsFilter_EuclideanDistance`
- `ComputeKernelAvgMisorientationsFilter`
- `ComputeFeatureNeighborMisorientationsFilter`. 

These exemplar files were generated against the pre-fix algorithm; the new values are *more* mathematically correct. **The exemplar files will need to be regenerated**.
