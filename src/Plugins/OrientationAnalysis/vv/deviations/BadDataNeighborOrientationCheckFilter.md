# Deviations from DREAM3D 6.5.171: BadDataNeighborOrientationCheckFilter

This file lists every documented behavioral difference between this SIMPLNX filter and its DREAM3D 6.5.171 equivalent (`BadDataNeighborOrientationCheck`, source at `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/BadDataNeighborOrientationCheck.{h,cpp}` in DREAM3D 6.5.171).

Entries are referenced by stable ID (`BadDataNeighborOrientationCheckFilter-D<N>`) from the V&V report and from public migration guidance. The ID is stable across renames; the Filter UUID field is the permanent cross-reference anchor.

## Comparison summary

A direct A/B comparison was rerun on 2026-09-17 across all 27 algorithmic fixtures and one 4,444,713-cell Small IN100 case. Inputs were identical for all compared applications: the same `Quats`, `Phases`, `Mask`, `CrystalStructures`, `MisorientationTolerance`, and `NumberOfNeighbors` arrays and parameters were read from each shared input file.

| Comparison | Identical cases | Differing mask values |
|---|---|---|
| SIMPLNX vs DREAM3D 6.5.171 | 12 of 28 | 45,229, all SIMPLNX good / legacy bad |
| SIMPLNX vs local legacy build with the surgical fixes | 28 of 28 | 0 |

The patch-isolation sequence separated the causes. After the D1 loop-bound and shared orientation-precision corrections, 41 fixture differences remained: 37 from D2's stale cross-phase misorientation and four from the legacy tolerance conversion at the exact 5° boundary. Applying the D2 same-phase gate and matching the corrected tolerance conversion removed those differences. The final local legacy build reproduced SIMPLNX exactly on all four compared arrays in all 28 cases.

---

## BadDataNeighborOrientationCheckFilter-D1

| Field | Value |
|---|---|
| **Deviation ID** | `BadDataNeighborOrientationCheckFilter-D1` |
| **Filter UUID** | `3f342977-aea1-49e1-a9c2-f73760eba0d3` |
| **Status** | active |

**Symptom:** DREAM3D 6.5.171 fails to flip a bad voxel whose good-neighbor count is exactly equal to the user-supplied `NumberOfNeighbors`. SIMPLNX correctly flips such voxels. The 27 analytical fixtures contain 288 such mask differences; the production case adds 44,941, for 45,229 across the complete comparison.

**Root cause:** Bug in DREAM3D 6.5.171.

The legacy iterative-decay loop is `while(currentLevel > m_NumberOfNeighbors)` (`Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/BadDataNeighborOrientationCheck.cpp:299`). With user-supplied `NumberOfNeighbors = N`, this walks `currentLevel` from 6 down through `N + 1` and never executes the `currentLevel == N` iteration. A bad voxel with exactly N good neighbors can therefore never be flipped — contradicting the parameter description ("Required Number of Neighbors") which implies that count to be sufficient.

SIMPLNX corrects this to `while(currentLevel >= m_InputValues->NumberOfNeighbors)` (`src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/Algorithms/BadDataNeighborOrientationCheck.cpp:129`). The fix is also explicitly documented as "BUG: Fix only checking values greater than the supplied min number of neighbors" in the merge commit of PR #1499 and in the engineer's V&V test archive README at `bad_data_neighbor_orientation_check_v2/README.md` §"Issue 1".

**Affected users:** Anyone running the filter in DREAM3D 6.5.171 with `NumberOfNeighbors < 6` on a dataset where the bottom level matters (i.e., where any bad voxel's eligible-neighbor count equals the user's `NumberOfNeighbors`). In practice this is the typical usage — the Small IN100 reconstruction pipeline (the canonical DREAM3D example) uses `NumberOfNeighbors = 4`. The 6.5.171 output left bad voxels unflipped that should have been flipped, manifesting downstream as smaller-than-expected grain reconstructions, more "rough" grain boundaries, and lower fraction of good voxels.

**Recommendation:** Trust SIMPLNX. The 6.5.171 result was mathematically incorrect for the stated parameter semantics. The minimal legacy patch is a one-line change from `>` to `>=`; the root cause was proven by applying this fix (bundled with D2) to a local build of the legacy source — contact the DREAM3D team for the legacy-parity patch.

---

## BadDataNeighborOrientationCheckFilter-D2

| Field | Value |
|---|---|
| **Deviation ID** | `BadDataNeighborOrientationCheckFilter-D2` |
| **Filter UUID** | `3f342977-aea1-49e1-a9c2-f73760eba0d3` |
| **Status** | active |

**Symptom:** DREAM3D 6.5.171 can count a different-phase neighbor's misorientation as within tolerance if a *previous* same-phase neighbor's `w` was small, because the misorientation-threshold check sits outside the same-phase conditional and inherits the stale `w` from the prior iteration. SIMPLNX prevents this by moving the threshold check inside the same-phase conditional.

D1 masks D2 in the unmodified release comparison. The isolation rerun first corrected D1 while retaining the legacy phase gate; 37 fixture mask differences then appeared. Applying the same-phase gate used by SIMPLNX closed all 37. Four additional exact-tolerance differences were separately removed by matching the corrected degree-to-radian conversion. This sequence demonstrates D2 independently of D1.

**Root cause:** Bug in DREAM3D 6.5.171.

The legacy per-neighbor loop body is (`BadDataNeighborOrientationCheck.cpp:283-291`):

```cpp
if(m_CellPhases[i] == m_CellPhases[neighbor] && m_CellPhases[i] > 0)
{
  w = m_OrientationOps[phase1]->getMisoQuat(q1, q2, n1, n2, n3);
}
if(w < misorientationTolerance)  // <-- outside the same-phase conditional!
{
  neighborCount[i]++;
}
```

When the current neighbor has a different phase than the voxel, the `w = getMisoQuat(...)` assignment is skipped, and the subsequent `if(w < misorientationTolerance)` reads `w` from the *previous neighbor iteration* (or from `w`'s initial value `10000.0f` if no previous iteration matched). The previous iteration's `w` may be small (e.g., from a same-phase good neighbor with an identical orientation), in which case the comparison succeeds and the count is incorrectly bumped.

SIMPLNX moves both the misorientation computation AND the increment inside the same-phase conditional (`Algorithms/BadDataNeighborOrientationCheck.cpp:105-117`):

```cpp
if(cellPhases[voxelIndex] == cellPhases[neighborPoint] && cellPhases[voxelIndex] > 0)
{
  ebsdlib::QuatD quat2(quats[neighborPoint * 4], ...);
  quat2.positiveOrientation();
  ebsdlib::AxisAngleDType axisAngle = orientationOps[laueClass1]->calculateMisorientation(quat1, quat2);
  if(axisAngle[3] < misorientationTolerance)
  {
    neighborCount[voxelIndex]++;
  }
}
```

The bug is documented as "Issue 2" in the engineer's V&V test archive README at `bad_data_neighbor_orientation_check_v2/README.md`, and was bundled into PR #1499's REV cleanup.

**Affected users:** Anyone running the filter in DREAM3D 6.5.171 on a dataset with mixed phases adjacent to grain boundaries. The bug would manifest as voxels at phase boundaries being incorrectly flipped to "good" because they appear to have more within-tolerance neighbors than they actually do.

**Recommendation:** Trust SIMPLNX. The 6.5.171 result was mathematically incorrect. D1 and D2 were isolated in sequence on a local build of the legacy source, and the final corrected build reproduced SIMPLNX on every comparison case.

---

### EbsdLib 2.4.1 CubicOps precision improvement (precision improvement; not a behavioral deviation in this filter's test data)

SIMPLNX delegates misorientation math to `ebsdlib::LaueOps::calculateMisorientation` (EbsdLib 2.4.1+); legacy 6.5.171 delegates to `OrientationLib::CubicOps::getMisoQuat`. The modern API recovers ~0.02° of precision for cubic misorientations on 4-fold, 3-fold, or 2-fold symmetry operators by replacing the precision-fragile `acos(w)` near 1 with the stable `2·atan2(|v|, w)` form. The improvement is documented in the EbsdLib 2.4.1 release notes.

**Not observed as a deviation in this filter** because the engineer's test fixtures do not include any voxel pair whose misorientation lands on a cubic sym op. The improvement is real and affects other downstream filters (see `ComputeFeatureFaceMisorientationFilter` V&V cycle's D4); for `BadDataNeighborOrientationCheck` specifically, this is a transparent dependency upgrade.

---

## Comparison artifacts

The archived record contains the 27 analytical fixtures, the Small IN100 production input, matching pipelines for all compared applications, all outputs, the original engineer test design, source snapshots, execution logs, and machine-readable array comparisons. The archive is reproducible from its included scripts and records the input and binary hashes used for the 2026-09-17 rerun.
