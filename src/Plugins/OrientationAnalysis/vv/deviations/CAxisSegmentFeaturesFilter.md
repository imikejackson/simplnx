# Deviations from DREAM3D 6.5.171: CAxisSegmentFeaturesFilter

This file lists every documented behavioral difference between this SIMPLNX filter and its DREAM3D 6.5.171 equivalent.

Entries are referenced by stable ID (`CAxisSegmentFeaturesFilter-D<N>`) from the V&V report and from public migration guidance. The ID is stable across renames; the Filter UUID field is the permanent cross-reference anchor.

## Comparison summary

A ten-case analytical release matrix was run on 2026-09-17 against DREAM3D 6.5.171, a deterministic local legacy build, DREAM3D-NX 7.4.0, source-rebuilt DREAM3D-NX 7.4.1 with EbsdLib 2.2.0, and DREAM3D-NX 7.4.2 with EbsdLib 3.1.2. Every successful output matches the independent Class 1 partition. Controlled release-before/release-after execution proves D1, D4, D5, and D6; repeated legacy runs prove D2 is a label permutation; D3 remains a deliberate domain guard.

---

## CAxisSegmentFeaturesFilter-D1

| Field | Value |
|---|---|
| **Deviation ID** | `CAxisSegmentFeaturesFilter-D1` |
| **Filter UUID** | `9fe07e17-aef1-4bf1-834c-d3a73dafc27d` |
| **Status** | active |

**Symptom:** PR #1466 (2025-11-14) introduced this, so **7.4.1 (2026-03-23) is the only affected release** — 7.4.0 (2025-10-27) and earlier predate it. In affected builds the filter could report one extra (empty) feature, shift every FeatureId up by one, or grow a feature from a masked-out / unindexed voxel — whenever voxel 0 of the image was not a legitimate seed. 6.5.171 never exhibits this.

**Root cause:** Bug (SIMPLNX). The shared driver `src/simplnx/Utilities/SegmentFeatures.cpp::execute()` started the flood fill from the raw index `seed = 0` without calling `getSeed()`, so the first seed was neither validated against the mask/phase requirements nor stamped with its FeatureId. Legacy `SegmentFeatures::execute()` obtains every seed — including the first — from `getSeed()`. Restored in this V&V cycle (`seed = getSeed(gnum, nextSeed)` before the loop); the fix also applies to `EBSDSegmentFeatures` and `ScalarSegmentFeatures`, which share the driver. Pinned by the `Class 1 Analytical (Mask Excludes Voxel 0)` and `Execute Error - No Features Found (-87000)` test cases.

**Empirical patch proof:** On the masked fixture, 7.4.1 writes real ids 2 and 3 with an empty id-1 tuple; 7.4.2 writes ids 1 and 2 with no empty tuple, matching the oracle and legacy feature count. On the all-masked fixture, 7.4.1 incorrectly succeeds with zero real features and a two-tuple Active array; 7.4.2 and DREAM3D 6.5.171 both report `-87000`. Source ancestry confirms 7.4.0 predates the raw-seed regression and 7.4.2 contains the correction.

**Affected users:** DREAM3D-NX 7.4.1 users whose datasets have a masked-out, unindexed, or already-owned cell at linear index 0 — common in EBSD scans with a mask. Legacy 6.5.171 users are unaffected.

**Recommendation:** Trust SIMPLNX at or after this V&V commit (which agrees with 6.5.171). Results from affected intermediate SIMPLNX builds on masked data should be regenerated.

---

## CAxisSegmentFeaturesFilter-D2

| Field | Value |
|---|---|
| **Deviation ID** | `CAxisSegmentFeaturesFilter-D2` |
| **Filter UUID** | `9fe07e17-aef1-4bf1-834c-d3a73dafc27d` |
| **Status** | active |

**Symptom:** FeatureIds from 6.5.171 are a different (random) labeling on every run, while SIMPLNX produces the same FeatureIds on every run; the two versions never produce bit-identical FeatureIds arrays.

**Root cause:** Algorithmic choice. 6.5.171 hard-codes `m_RandomizeFeatureIds = true` and seeds its RNG from the wall clock (`CAxisSegmentFeatures.cpp::initializeVoxelSeedGenerator`), and the option is not exposed as a pipeline parameter — legacy output labeling is irreproducible by construction. SIMPLNX exposes `Randomize Feature Ids` as a parameter (default `false`) and, when enabled, uses a fixed-seed `std::mt19937_64` (`ClusterUtilities::RandomizeFeatureIds`), so output is deterministic either way. The segmentation *partition* (which cells share a feature) is unaffected; the A/B runs (2026-07-22, rerun with the 3-D fixture 2026-07-24, `vv/comparisons/CAxisSegmentFeaturesFilter/`) matched partitions exactly on all five fixtures.

**Empirical proof:** Three independent 6.5.171 runs on the same chain input produced three different raw FeatureId permutations and one repeated permutation, while all four canonicalized to the same oracle partition. A local legacy build with randomization disabled reproduced the 7.4.2 raw FeatureIds on every successful shared ImageGeom case.

**Affected users:** Anyone diffing raw FeatureIds arrays between versions or between two 6.5.171 runs; downstream statistics keyed by feature (sizes, misorientations) are invariant to the labeling.

**Recommendation:** Trust SIMPLNX. Deterministic labeling is strictly more reproducible; compare segmentations at the partition level when validating against legacy runs.

---

## CAxisSegmentFeaturesFilter-D3

| Field | Value |
|---|---|
| **Deviation ID** | `CAxisSegmentFeaturesFilter-D3` |
| **Filter UUID** | `9fe07e17-aef1-4bf1-834c-d3a73dafc27d` |
| **Status** | active |

**Symptom:** On data containing non-hexagonal phases that participate in segmentation, 6.5.171 silently produces a segmentation; SIMPLNX fails with error `-8363` (and `-8364` for phase values with no CrystalStructures entry).

**Root cause:** Algorithmic choice (deliberate SIMPLNX guard). The c-axis is only a physically meaningful unique axis for hexagonal (6/m, 6/mmm) Laue classes, but the c-axis math itself never consults the crystal structure — legacy computes the [001] misalignment for cubic/other phases and returns scientifically meaningless groupings. SIMPLNX validates that every cell that can participate in segmentation (phase > 0, not masked out) has a Hexagonal_High or Hexagonal_Low crystal structure (`CAxisSegmentFeatures.cpp::operator()`). Unindexed (phase 0) cells and masked-out cells are exempt, since they can never seed or join a feature.

**Empirical evidence:** DREAM3D 6.5.171 segments the participating cubic fixture into the partition `[1,1,2]` after label canonicalization. DREAM3D-NX 7.4.0, 7.4.1, and 7.4.2 all reject it with `-8363`. This run documents the migration difference; the Class 1/domain argument, not legacy output, establishes that the NX rejection is correct.

**Affected users:** Anyone who ran the legacy filter on multi-phase data with non-hexagonal phases — their legacy results for those phases were never meaningful. Pure-hexagonal workflows are unaffected.

**Recommendation:** Trust SIMPLNX. The error is a correctness guard; mask out non-hexagonal phases (supported) to segment only the hexagonal cells.

---

## CAxisSegmentFeaturesFilter-D4

| Field | Value |
|---|---|
| **Deviation ID** | `CAxisSegmentFeaturesFilter-D4` |
| **Filter UUID** | `9fe07e17-aef1-4bf1-834c-d3a73dafc27d` |
| **Status** | active |

**Symptom:** In every released version up to and including 7.4.1, SIMPLNX rejected (error `-8363`) any dataset containing unindexed (phase 0) cells — whose `CrystalStructures[0]` entry is the conventional `999` sentinel — or masked-out non-hexagonal cells; 6.5.171 processes such datasets normally.

**Root cause:** Bug (SIMPLNX). The D3 validation loop checked every cell's crystal structure, including phase-0 cells and cells excluded by the mask, neither of which can ever participate in segmentation (`getSeed` requires phase > 0 and a set mask bit; `determineGrouping` requires equal phases and a set mask bit). It also indexed `CrystalStructures[phase]` without a bounds check, so an out-of-range phase value read out of bounds instead of producing an error. Fixed by exempting phase ≤ 0 and masked-out cells and adding the `-8364` bounds error. Pinned by the `Phase 0 (Unindexed) Cells Tolerated`, `Masked Non-Hexagonal Cells Tolerated`, and `Execute Error - Phase Out of Ensemble Bounds (-8364)` test cases; the 2026-07-22 A/B run confirms post-fix parity with 6.5.171 on phase-0 data (TC4).

**Empirical patch proof:** Releases 7.4.0 and 7.4.1 reject the phase-zero and masked-cubic fixtures with `-8363`; release 7.4.2 accepts both and matches the oracle and legacy partitions. The older releases also execute an out-of-range phase fixture after the unchecked access; release 7.4.2 stops with `-8364` and reports the cell, phase value, selected array, and tuple count.

**Affected users:** Users of any released version through 7.4.1 with EBSD scans containing unindexed points — a very common case — or deliberately masked non-hexagonal phases.

**Recommendation:** Trust SIMPLNX at or after this V&V commit; upgrade if the filter spuriously rejects data with unindexed points.

---

## CAxisSegmentFeaturesFilter-D5

| Field | Value |
|---|---|
| **Deviation ID** | `CAxisSegmentFeaturesFilter-D5` |
| **Filter UUID** | `9fe07e17-aef1-4bf1-834c-d3a73dafc27d` |
| **Status** | active |

**Symptom:** SIMPLNX accepts a RectGrid geometry as input (6.5.171 accepts only Image geometry); in every released version up to and including 7.4.1, selecting a RectGrid geometry passed preflight and then crashed at execute.

**Root cause:** Bug (SIMPLNX). The algorithm fetched the geometry with a stale `getDataAs<ImageGeom>()` cast, returning a null pointer for RectGrid input, which the shared segmentation driver dereferenced. Fixed to `getDataAs<IGridGeometry>()`, matching the parameter's allowed types and the sibling EBSD/Scalar segmentation algorithms. Pinned by the `Class 1 Analytical (RectGrid Geometry)` test case.

**Empirical patch proof:** Releases 7.4.0 and 7.4.1 terminate with signal 11 on the same three-cell RectGrid input. Release 7.4.2 completes and writes FeatureIds `[1,1,2]`, exactly matching the analytical expectation. Legacy comparison is not applicable because DREAM3D 6.5.171 does not expose RectGrid input for this filter.

**Affected users:** Users of any released version through 7.4.1 segmenting RectGrid data — e.g., regularized serial-sectioning data. Legacy users are unaffected (the capability does not exist in 6.5.171).

**Recommendation:** Trust SIMPLNX at or after this V&V commit.

---

## CAxisSegmentFeaturesFilter-D6

| Field | Value |
|---|---|
| **Deviation ID** | `CAxisSegmentFeaturesFilter-D6` |
| **Filter UUID** | `9fe07e17-aef1-4bf1-834c-d3a73dafc27d` |
| **Status** | active |

**Symptom:** In DREAM3D-NX releases through 7.4.0, the cell partition and FeatureIds are correct, but the Cell Feature AttributeMatrix contains one extra empty tuple. DREAM3D-NX 7.4.1 resolves this defect. For example, a four-feature chain has six Active tuples instead of the required five (background plus four real features).

**Root cause:** Bug (SIMPLNX). Before PR #1466, the shared segmentation driver stored `gnum` directly as `m_FoundFeatures` after the loop, although `gnum` was the next unused FeatureId. The filter then resized the feature AttributeMatrix to `m_FoundFeatures + 1`, producing one extra tuple. PR #1466 corrected the count to `gnum - 1`; that correction resolved D6 in 7.4.1 but simultaneously introduced D1 by replacing the validated first `getSeed()` call with raw index zero.

**Empirical proof:** DREAM3D-NX 7.4.0 produces exactly one extra Active tuple in every successful matrix case: real feature count plus background plus one. The 7.4.1 ordinary-seed cases and all 7.4.2 cases have the correct tuple count. DREAM3D 6.5.171 also has the correct count.

**Affected users:** Users of released DREAM3D-NX versions through 7.4.0 whose downstream processing assumes every non-background feature tuple is populated. The extra tuple can appear as an empty or zero-valued feature in feature-level statistics and exports.

**Recommendation:** Trust DREAM3D-NX 7.4.2. Release 7.4.1 corrects this count but contains D1; 7.4.2 is the first release that has neither defect.
