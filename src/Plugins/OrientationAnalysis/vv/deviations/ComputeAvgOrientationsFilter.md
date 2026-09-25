# Deviations from DREAM3D 6.5.171: ComputeAvgOrientationsFilter

This file lists every documented behavioral difference between this SIMPLNX filter and its DREAM3D 6.5.171 equivalent, `FindAvgOrientations` (SIMPL UUID `bf7036d8-25bd-540e-b6de-3a5ab0e42c5f`).

Entries are referenced by stable ID (`ComputeAvgOrientationsFilter-D<N>`) from the V&V report and from public migration guidance. The deviations apply **only to the Rodrigues (original) averaging method**; the von Mises-Fisher and Watson methods are new in SIMPLNX and have no 6.5.171 equivalent.

> **Status:** Empirically validated against the official DREAM3D 6.5.171 release and refreshed 2026-09-17 with a local legacy build carrying the surgical D3 correction. Every application read the same legacy-format input bytes.
>
> Four fixtures were run: **(A)** a realistic 480,000-cell / 409-tuple crop, **(B)** a five-tuple fixture forcing D2 and both D3 mechanisms, **(C1)** the Class 1 analytical unit-test fixture, and **(C2)** the cubic-symmetry invariant fixture.
>
> **Headline:** SIMPLNX and 6.5.171 agree within `1e-6` on every real feature (`AvgQuats` max 8.94e-8; `AvgEulerAngles` max 4.77e-7; zero sign flips). D2 remains an intentional FeatureId-0 policy difference. D3 is now empirically confirmed at tuple zero and for zero-count tuples at index ≥1; the surgical D3 correction makes those empty-tuple outputs agree with SIMPLNX. D4 is sub-epsilon. D1 remains downgraded.

---

## ComputeAvgOrientationsFilter-D1 — NOT A DEVIATION (defensive normalization)

| Field | Value |
|---|---|
| **Deviation ID** | `ComputeAvgOrientationsFilter-D1` |
| **Filter UUID** | `086ddb9a-928f-46ab-bad6-b1498270d71e` |
| **Status** | retired 2026-06-30 — not a deviation (defensive normalization) |

**Symptom:** Withdrawn. The originally reported symptom was a possible `q` versus `−q` difference when SIMPLNX appends `.getPositiveOrientation()` after normalization and the average lands with `w < 0`.

**Root cause:** Algorithmic choice. This record is not a deviation because `.getPositiveOrientation()` is a defensive canonicalization and no behavioral difference was demonstrated.

**Evidence (empirical + structural, 2026-06-30):** It could not be made to diverge. The Rodrigues algorithm builds each average by symmetry-reducing every voxel toward the running average (`getNearestQuat`), so the result lands in the fundamental zone — rotation angle ≤ 180° ⇒ `w ≥ 0`. A deliberate fixture (feature with `Rz(170°)` + `Rz(200°)`, intended to push the sum to `w < 0`) was reduced identically by **both** 6.5.171 and SIMPLNX to `(0, 0, 0.0436, 0.999)` (`w > 0`) — no sign flip. Even in the only theoretical corner where `w < 0` could survive (pure-triclinic, the incremental average overshooting 180°), the two quaternions represent the **same physical orientation** (`q ≡ −q`), so there is no correctness or downstream effect. `getPositiveOrientation()` is therefore a harmless defensive canonicalization, not a behavioral difference from legacy.

**Affected users:** None.

**Recommendation:** Do not count as a deviation. Keep this record so the extra `getPositiveOrientation()` call is not mistaken for an unverified deviation in a future audit.

---

## ComputeAvgOrientationsFilter-D2

| Field | Value |
|---|---|
| **Deviation ID** | `ComputeAvgOrientationsFilter-D2` |
| **Filter UUID** | `086ddb9a-928f-46ab-bad6-b1498270d71e` |
| **Status** | active |

**Symptom:** Voxels labeled `FeatureId == 0` (with `Phase > 0`) contribute to averaging in SIMPLNX and produce a computed average for feature 0; in 6.5.171 they are skipped and feature 0 is left at `(0,0,0,0)`.

**Empirical (2026-06-30):** Forcing fixture B — feature 0 given two phase-1 cells (identity + `Rz(90°)`). SIMPLNX computed `AvgQuats[0] = (0, 0, 0.382683, 0.92388)` (the `Rz(45°)` average); the official 6.5.171 wrote `AvgQuats[0] = (0, 0, 0, 0)` (skipped). A clear, large divergence — not sub-epsilon. (The realistic fixture A had no `FeatureId == 0` cells, so this gate is dormant on conventional EBSD data.)

**Root cause:** Algorithmic choice. Legacy gates accumulation on `m_FeatureIds[i] > 0 && m_CellPhases[i] > 0` (`FindAvgOrientations.cpp:246`). SIMPLNX gates on `currentPhase > 0` only (`ComputeAvgOrientations.cpp:412`), deliberately allowing the documented use-case of averaging an unlabeled bag of orientations all tagged `FeatureId 0` (algorithm comment lines 401–411).

**Affected users:** Only datasets where `FeatureId 0` legitimately carries `Phase > 0` data (atypical — `FeatureId 0` is conventionally background/unindexed with `Phase 0`). For conventional data there is no observable difference.

**Recommendation:** Either acceptable within tolerance. For conventional data the outputs match; the SIMPLNX behavior is a strict superset supporting an additional use-case (computing the average of an unlabeled orientation set).

---

## ComputeAvgOrientationsFilter-D3

| Field | Value |
|---|---|
| **Deviation ID** | `ComputeAvgOrientationsFilter-D3` |
| **Filter UUID** | `086ddb9a-928f-46ab-bad6-b1498270d71e` |
| **Status** | active |

**Symptom:** For a feature with zero contributing voxels, SIMPLNX writes identity quaternion `(0,0,0,1)` and zero Euler angles. DREAM3D 6.5.171 leaves tuple zero at `(0,0,0,0)` and produces NaN quaternion and Euler values for zero-count tuples at index ≥1. Both mechanisms were reproduced directly in the forcing and analytical fixtures.

**Root cause:** Bug in 6.5.171, via two mechanisms. (a) The legacy init and finalize loops both run `for(i = 1; i < totalFeatures)` (`FindAvgOrientations.cpp:239,263`), so feature 0 is never finalized and keeps `(0,0,0,0)`. (b) For a zero-count feature at index ≥1, the legacy finalize sets `Identity` and then divides it by zero. SIMPLNX (`ComputeAvgOrientations.cpp:440–444`) includes tuple zero and, when `counts == 0`, writes identity plus zero Euler angles and continues. Applying those same changes to a local build of the legacy source reproduced the SIMPLNX empty-tuple outputs on all forcing and analytical cases.

**Affected users:** Anyone whose Feature Attribute Matrix contains a feature index with no contributing voxels (feature 0 always; index ≥ 1 after feature removal/renumbering gaps), and any consumer of feature 0's `AvgQuats` in legacy.

**Recommendation:** Trust SIMPLNX. The legacy zero-count result was uninitialized/undefined (`(0,0,0,0)` or division by zero); SIMPLNX's identity quaternion is well-defined.

---

## ComputeAvgOrientationsFilter-D4

| Field | Value |
|---|---|
| **Deviation ID** | `ComputeAvgOrientationsFilter-D4` |
| **Filter UUID** | `086ddb9a-928f-46ab-bad6-b1498270d71e` |
| **Status** | active |

**Symptom:** `AvgEulerAngles` (and possibly `AvgQuats` at the last ULPs) may differ from 6.5.171 at the sub-epsilon level.

**Root cause:** Library + precision. Legacy uses `QuaternionMathF` arithmetic and `OrientationTransforms::qu2eu` for the quaternion→Euler conversion; SIMPLNX uses EbsdLib `ebsdlib::QuatF` and `QuaternionFType::toEuler()`. Both operate in `float32`, but the differing intermediate-math implementations and quaternion→Euler routines can produce last-bit differences.

**Empirical (refreshed 2026-09-17):** Across all 408 real features, `AvgEulerAngles` differed by at most **4.77e-7** (86 of 1224 components in the 1e-7–1e-6 band, the rest below 1e-7; mean 2.29e-8). `AvgQuats` differed by at most 8.94e-8. This confirms the divergence is float32 round-off in the two independent library code paths.

**Affected users:** Anyone doing bit-exact comparison of `AvgEulerAngles` between versions. Differences are at the floating-point-noise level and not materially significant for any downstream calculation.

**Recommendation:** Either acceptable within tolerance (≈1e-6). Neither implementation is more correct; the difference is float round-off in independent library code paths.
