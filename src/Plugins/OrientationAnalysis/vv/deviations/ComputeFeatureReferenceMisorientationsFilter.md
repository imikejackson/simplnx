# Deviations from DREAM3D 6.5.171: ComputeFeatureReferenceMisorientationsFilter

This file lists every documented behavioral difference between this SIMPLNX filter and its DREAM3D 6.5.171 equivalent (`FindFeatureReferenceMisorientations`, source at `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindFeatureReferenceMisorientations.{h,cpp}` in DREAM3D 6.5.171).

Entries are referenced by stable ID (`ComputeFeatureReferenceMisorientationsFilter-D<N>`) from the V&V report and from public migration guidance. The ID is stable across renames; the Filter UUID field is the permanent cross-reference anchor.

## Comparison summary

The legacy comparison was rerun empirically on 2026-09-17 using the six analytical fixtures and a 748,800-cell Small IN100 case in both reference modes. Every application read byte-identical scientific inputs. DREAM3D 6.5.171 completed Mode 0 but terminated with a segmentation fault in all four direct Mode 1 cases. A local legacy build with the surgical Mode 1 and orientation-precision corrections completed every case, matched SIMPLNX bit-for-bit on the three analytical Mode 1 fixtures, and agreed within one per-cell float32 ULP on Small IN100.

Mode 0 quantified D1 rather than merely inferring it. On the analytical fixtures, baseline differences reached 1.86e-4°; on Small IN100 they reached 0.0738° per cell and 0.01761° per feature average. The corrected local legacy build reduced the per-cell production residual to at most 1.907e-6°. Mode 1 exposed the independent legacy crash documented as D2.

---

## ComputeFeatureReferenceMisorientationsFilter-D1

| Field            | Value                                                       |
|------------------|-------------------------------------------------------------|
| **Deviation ID** | `ComputeFeatureReferenceMisorientationsFilter-D1`           |
| **Filter UUID**  | `24b54daf-3bf5-4331-93f6-03a49f719bf1`                      |
| **Status**       | active |

**Symptom:** Per-cell and per-feature misorientations differ between SIMPLNX and DREAM3D 6.5.171 at the precision level. The empirical rerun measured up to 1.86e-4° on the analytical Mode 0 fixtures and, on Small IN100, up to 0.0738° per cell and 0.01761° per feature average. All six analytical fixtures remain within `1e-3°` of the independent expected values.

**Root cause:** **Precision** — not an algorithm change in either implementation.

The deviation traces to the EbsdLib 2.4.1 release commit `5c8c993` (BlueQuartz Software, 2026-05-29), which replaces a precision-fragile `acos(w)` form in `CubicOps::calculateMisorientationInternal` with a numerically-stable `2·atan2(|v|, w)` form using the explicit reduced-quaternion `v` components. The precision improvement is real and mathematically more correct; it manifests for cubic misorientations whose minimum-rotation-axis representation lies on or near a cubic symmetry operator (e.g., 90° about the cubic c-axis is a 4-fold sym op of m-3m; pre-fix `acos`-form yielded `~0.02°` residual due to float32-input ULP noise, post-fix yields the mathematically correct value).

This filter is a clean Port of `FindFeatureReferenceMisorientations` from DREAM3D 6.5.171; the SIMPLNX algorithm reproduces the legacy two-mode dispatch + per-voxel misorientation accumulation + per-feature averaging structure exactly. The legacy filter consumes `OrientationLib::CubicOps::getMisoQuat` (pre-fix `acos`-form, float32); the SIMPLNX filter consumes `ebsdlib::CubicOps::calculateMisorientation` (post-fix `2·atan2`-form, QuatD). The difference is entirely in the EbsdLib precision improvement, NOT in this filter.

For the full root-cause walkthrough of the EbsdLib precision improvement, see the precedent characterization in `vv/deviations/BadDataNeighborOrientationCheckFilter.md` §"Non-deviations" → "EbsdLib 2.4.1 CubicOps precision improvement". The characterization there applies equally to this filter; the only difference is that this filter's per-feature averaging amplifies the per-voxel precision shift across the feature's voxels (typically hundreds to thousands for real EBSD data), making the deviation more visible at the per-feature output level than at the per-voxel level.

**Affected users:** Anyone running this filter in DREAM3D 6.5.171 on EBSD data with cubic-phase grains that have grain boundaries near 4-fold (90° c-axis), 3-fold (120° [111]), or 2-fold (180° face-diagonal) cubic symmetry operators, and comparing per-feature `Feature Avg Misorientations` output across the version boundary. On non-cubic-phase data, no deviation. On cubic data without sym-op-aligned boundaries, no observable deviation.

**Recommendation:** **Trust SIMPLNX.** The 6.5.171 result was limited by float32-input ULP noise amplified by `acos`-near-1 catastrophic cancellation; SIMPLNX returns the mathematically correct value. The `~0.02°` shift is well below typical EBSD measurement resolution (per the BadDataNeighborOrientationCheckFilter V&V cycle's precedent characterization) and will not materially affect downstream microstructural analyses for users migrating from DREAM3D 6.5.171.

---

## ComputeFeatureReferenceMisorientationsFilter-D2

| Field            | Value                                                       |
|------------------|-------------------------------------------------------------|
| **Deviation ID** | `ComputeFeatureReferenceMisorientationsFilter-D2`           |
| **Filter UUID**  | `24b54daf-3bf5-4331-93f6-03a49f719bf1`                      |
| **Status**       | active |

**Symptom:** DREAM3D 6.5.171 terminates with a segmentation fault whenever Mode 1 (Euclidean-distance reference) is selected. It produces no output for any direct Mode 1 fixture. SIMPLNX completes all Mode 1 cases and matches the independent analytical expectations.

**Root cause:** **Bug** in DREAM3D 6.5.171. The legacy `execute()` method obtains `totalFeatures` by locking the Mode-0-only `AvgQuats` input pointer before dispatching to the Mode 1 center-selection path. `dataCheck()` does not initialize that pointer in Mode 1, so the dereference is invalid. A surgical one-line correction uses the tuple count of the always-created feature-average output instead, which is valid in both modes.

**Affected users:** Every DREAM3D 6.5.171 user who selected the Euclidean-distance reference mode. The filter terminates before producing reference-misorientation output.

**Recommendation:** **Trust SIMPLNX.** The SIMPLNX port does not depend on the Mode-0-only input for its feature count. A local legacy build with the surgical correction completed all four Mode 1 cases, matched SIMPLNX bit-for-bit on the three analytical fixtures, and agreed within one per-cell float32 ULP on Small IN100.

---

## Non-deviations (algorithm characteristics common to both filters)

The following behaviors are NOT deviations — SIMPLNX and 6.5.171 agree on them. Captured here so future engineers don't re-discover them and propose them as deviations.

### Raster-order tie-break in `centers[]` selection (Mode 1)

Both implementations use `if(distance >= centerDistances[featureId])` in the Mode 1 pre-loop that selects each feature's reference voxel. The `>=` (rather than `>`) means that when two or more voxels within a feature have identical `GBEuclideanDistances` values, the LATER voxel (in linear iteration order) overwrites earlier candidates and is selected as the feature's reference. The choice is therefore raster-order dependent — different DataStructure layouts that expose the same logical voxels in a different iteration order would yield different `centers[]` and different `EuclideanCenters`. **Both filters share this behavior** — verified by source inspection of the legacy `FindFeatureReferenceMisorientations::execute()` lines 320-325 vs SIMPLNX `ComputeFeatureReferenceMisorientations.cpp` lines 89-103.

### Background voxel and unphased voxel handling

Both implementations skip voxels where `featureIds[i] == 0` (background) or `cellPhases[i] == 0` (unphased). In both, the per-voxel `FRM` array is initialized to 0 (or `fill(0.0f)` in SIMPLNX) and skipped voxels retain that zero value. Per-feature `avgRefMis` is computed only over the non-skipped voxels in each feature; if a feature consists entirely of skipped voxels, its `count == 0` and `avgRefMis[fid]` is set to `0`. **Both filters share this behavior** — algorithm characteristic, not a defect.

### Background feature (featureId = 0) → `avgRefMis[0] == 0`

Both implementations leave `avgRefMis[0]` at its initialized `0.0f` value (since the main per-voxel loop only computes misorientations for `featureIds[i] > 0`, and the per-feature finalize loop iterates `for(i = 1; i < totalFeatures; i++)`, skipping index 0 entirely). Algorithm characteristic, not a defect.

---

## Comparison artifacts

The archived comparison record contains six analytical fixtures, the Small IN100 production case, direct Mode 0 and Mode 1 pipelines, an independent Mode 1 emulation through Mode 0, debugger evidence for the release crash, outputs from the local surgically corrected legacy build, and machine-readable comparisons. Thirty-six pipelines were executed. All twelve SIMPLNX and twelve corrected-legacy executions succeeded; DREAM3D 6.5.171 succeeded on the eight Mode 0/emulation cases and terminated with exit 139 on all four direct Mode 1 cases.
