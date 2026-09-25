# Deviations from DREAM3D 6.5.171: ComputeKernelAvgMisorientationsFilter

This file lists every documented behavioral difference between this SIMPLNX filter and its DREAM3D 6.5.171 equivalent (`FindKernelAvgMisorientations`, source at `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindKernelAvgMisorientations.{h,cpp}` in DREAM3D 6.5.171).

Entries are referenced by stable ID (`ComputeKernelAvgMisorientationsFilter-D<N>`) from the V&V report and from public migration guidance. The ID is stable across renames; the Filter UUID field is the permanent cross-reference anchor.

## Comparison summary

The legacy comparison was refreshed empirically on 2026-09-17 with two shared inputs and three kernels per input: symmetric `{1,1,1}` plus asymmetric `{1,1,2}` and `{2,1,1}`. All compared applications read bit-identical `FeatureIds`, `Phases`, and `Quats`. The symmetric cases quantify D1 while D2 is dormant; the asymmetric cases exercise D2 in both the over-reach and truncation directions. A local legacy build with the surgical D1/D2 corrections agrees with SIMPLNX to at most 9.54e-7° in every case, with zero cells above 0.01°.

**Runtime A/B setup.** A seeded 12×12×12 synthetic volume and a 256,000-cell Small IN100 crop were each run with kernels `{1,1,1}`, `{1,1,2}`, and `{2,1,1}`. All six cases use the legacy-equivalent per-grain mode. The archived generators, matching pipelines, inputs, outputs, logs, and comparison scripts make the experiment reproducible without relying on ephemeral scratch files.

---

## ComputeKernelAvgMisorientationsFilter-D1

| Field            | Value                                                                            |
|------------------|----------------------------------------------------------------------------------|
| **Deviation ID** | `ComputeKernelAvgMisorientationsFilter-D1`                                       |
| **Filter UUID**  | `61cfc9c1-aa0e-452b-b9ef-d3b9e6268035`                                           |
| **Status**       | active                     |

**Symptom:** Per-cell `KernelAverageMisorientations` values differ between SIMPLNX and DREAM3D 6.5.171 at the orientation-math precision level even when D2 is dormant. On Small IN100 with kernel `{1,1,1}`, the measured mean absolute difference is 0.00154° and the maximum is 0.05595°; 2,281 of 256,000 cells exceed 0.01°. A local legacy build with the surgical precision correction reduces the maximum residual to 9.54e-7°.

**Dependency:** the fix lives in EbsdLib commit `5c8c993`, contained in the `v2.4.1` tag. As of this PR the SIMPLNX `vcpkg.json` pins `ebsdlib version>=2.4.1`, so the **standard vcpkg build now links the fixed EbsdLib** and the artifact no longer appears in any supported configuration: both the standard build (`NX-Com-Qt69-Vtk95-Rel`) and the local-source build (`NX-Com-Qt69-Vtk95-Rel-EbsdLib`, `SIMPLNX_USE_LOCAL_EBSD_LIB=ON`) produce the correct, self-miso-free result. The V&V data-fixture unit tests assert the exact analytical oracle (margin `1e-3`) and pass in both configurations. The artifact reappears only if EbsdLib is pinned below `2.4.1` (e.g. an older vcpkg baseline); the *Empirical confirmation* below was captured against the pre-fix `2.4.0` — the version that shipped before this PR — to characterize the symptom and verify the fix.

**Root cause:** **Precision** — not an algorithm change in either implementation.

The deviation traces to the EbsdLib 2.4.1 release commit `5c8c993` (BlueQuartz Software, 2026-05-29), which replaces a precision-fragile `acos(w)` form in `CubicOps::calculateMisorientationInternal` with a numerically-stable `2·atan2(|v|, w)` form using the explicit reduced-quaternion `v` components. The precision improvement is real and mathematically more correct; for `ComputeKernelAvgMisorientationsFilter` specifically it manifests *more strongly than for the per-pair misorientation filters in this cycle* because the kernel inclusion of the focal voxel triggers a per-cell self-misorientation call. For the pre-fix `acos(w)`-form:

- `q_self_miso = q_focal * q_focal.conjugate() = (0, 0, 0, 1)` mathematically (identity quaternion), so the *true* self-misorientation is exactly 0°.
- The error is **not** introduced by the raw `q * q.conjugate()` product (its `w` component is `|q|² ≈ 1`); it is introduced by the **symmetry reduction**. `calculateMisorientationInternal` maximizes `wmin` over the 24 cubic sym-op candidates, evaluating three candidate forms per sym op: `qco.w()`, `(qco.z() + qco.w())/√2` (the 4-fold-about-c form), and `(qco.x()+qco.y()+qco.z()+qco.w())/2`. For the identity misorientation, several candidates equal 1.0 mathematically — but on float32-sourced quaternions a non-trivial candidate such as `(qco.z() + qco.w())/√2` evaluates to `1 − ε` (with `ε ≈ 1.7e-8`) and can be selected as the maximum. **Which focal quaternions trigger this is candidate-dependent: most reduce on the trivial `qco.w() == 1.0` branch and yield exactly 0; only those whose maximizing candidate is a non-trivial sym-op form land at `1 − ε`.**
- `acos(1 − ε)` near 1 is precision-fragile: the derivative of `acos` at 1 is `-1/√(1-x²) → -∞`, so a `1-ULP` error in `wmin` propagates to a `√(2ε)`-scale error in the angle (then doubled by the `2 * acos(wmin)` step). For `ε` of order `1e-8` this puts the spurious self-miso in the `~0.02–0.03°` range — the fix commit message cites `~0.02°`; the value **measured empirically on this branch is `0.0326°`** (see below). The exact constant depends on the winning sym-op candidate and the platform's float32 quantization, so treat the magnitude as order-of-`0.03°`, not a fixed number.
- The post-fix `2 * atan2(|v|, w)` form, using the **explicit** reduced-quaternion vector components, is numerically stable: components like `(qco.z() - qco.w())` evaluate to *exactly* 0 in IEEE-754 when `qco.z() == qco.w()` regardless of upstream float32 truncation, so `|v| = 0` and the result is exactly 0 for every identity self-misorientation.

The KAM filter is *more sensitive* than `ComputeFeatureNeighborMisorientations` and `BadDataNeighborOrientationCheck` to this precision improvement because:

1. **Self-misorientation contribution.** The KAM kernel includes the focal cell (via the `j=k=l=0` inner iteration). For each focal cell, the algorithm therefore makes one call to `calculateMisorientation` with `q1 == q2`. With pre-fix EbsdLib this call returns a spurious ~0.03° **for the subset of focal orientations whose symmetry reduction lands on a non-trivial sym-op candidate** (exactly 0 for the rest), which gets added to `totalMisorientation` and shifts that cell's average up by `(spurious_self / numVoxel)`. With fixed EbsdLib the call returns 0° for *every* focal orientation and contributes nothing. This is the *entire* KAM-specific deviation — see point 3.

2. **Same-feature large-N averaging.** For a cell in the middle of a large grain with kernel `{1,1,1}`, numVoxel = 27 (all same-feature). The cumulative effect of 27 small precision noises averages out somewhat, but the systematic self-miso contribution is always present.

3. **The focal self term is the dominant filter-specific amplifier, but not the only last-bit effect.** The KAM kernel always includes the focal cell, so the legacy self-misorientation artifact enters every affected average. On general 3D orientations, `QuatF` versus `QuatD` and other intermediate-rounding differences also affect distinct pairs at the `~1e-3°` scale in either direction. The six-case rerun quantifies the combined precision family; no symmetric-kernel cell differs by more than 0.056° on Small IN100, and the surgically corrected local legacy build reduces the residual to float32 ULP scale.

**Affected users:** Anyone migrating from DREAM3D 6.5.171 to SIMPLNX on cubic-phase EBSD data with this filter, or anyone running SIMPLNX with an obsolete EbsdLib baseline. On Small IN100 with the symmetric kernel, 2,276 of the 2,281 cells above 0.01° had legacy above SIMPLNX and five had the opposite sign; the effect is therefore strongly, but not absolutely, one-sided. Its magnitude depends on focal orientation, the number of contributing neighbors, and the distinct-pair rounding terms.

**Recommendation:** **Trust SIMPLNX (EbsdLib 2.4.1+).** The 6.5.171 result was limited by the well-understood `acos(w near 1)` precision pathology amplified by float32-sourced quaternion inputs; SIMPLNX returns the mathematically correct value. The shift is well below typical EBSD measurement resolution and will not materially affect downstream microstructural analyses, but the cumulative effect on KAM-based maps will be visibly smoother in the post-2.4.1 output. Users requiring exact 6.5.171 reproduction can compile against EbsdLib < 2.4.1 (not recommended).

For the full root-cause walkthrough of the EbsdLib precision improvement, see the precedent characterization in `vv/deviations/BadDataNeighborOrientationCheckFilter.md` §"Non-deviations" → "EbsdLib 2.4.1 CubicOps precision improvement". The characterization there applies equally to this filter, with the additional amplification factor described above.

**Empirical confirmation (V&V cycle, branch `topic/vv/ComputeFeatureNeighborMisorientationsFilter`, 2026-06-04):** The Class 1 / Class 4 data fixtures were run on Apple Silicon against both EbsdLib builds, with per-pair `calculateMisorientation` results instrumented:

- **Distinct-orientation pairs are exact on both builds.** Pairs of `5°`, `10°`, and `15°` apart returned `4.99991°`, `4.99988°`, etc. (`< 0.0002°` from the analytical value) on both vcpkg `2.4.0` and the fixed local EbsdLib. This rules out per-pair precision noise as a contributor and confirms point 3 above.
- **Self-misorientations are 0 for most focal orientations even pre-fix.** `q1 == q2` for the identity and for the `5°`, `10°`, `15°`-about-c focal cells returned *exactly* `0.0°` on vcpkg `2.4.0`.
- **Only the `20°`-about-c focal cells triggered the artifact pre-fix.** Their self-misorientation returned `0.0325663°` on vcpkg `2.4.0` (matching the `~0.033°` derived above), inflating those cells' KAM by `0.0326°/numVoxel`. Example: the 1D x-axis gradient fixture's last cell (`numVoxel = 2`) read `2.51628°` against an analytical `2.5°` — exceeding the test's `1e-3` margin.
- **The fixed EbsdLib zeroes every self-misorientation.** Rebuilding the same fixtures with the `NX-Com-Qt69-Vtk95-Rel-EbsdLib` preset (local EbsdLib at `5c8c993`), all three misorientation suites pass exactly: KAM `134/134` assertions, `ComputeFeatureNeighborMisorientations` `56/56`, `ComputeFeatureReferenceMisorientations` `238/238`.

The data-fixture unit tests assert the analytical oracle directly (margin `1e-3`, no tolerance for the pre-fix artifact). With EbsdLib pinned `≥ 2.4.1` in `vcpkg.json` this is the correct, regression-sensitive choice: it holds in every supported build and would immediately flag any future regression of the EbsdLib precision fix, rather than silently absorbing it under a loose tolerance.

**Runtime A/B confirmation (refreshed 2026-09-17):** identical synthetic and Small IN100 inputs were run through DREAM3D 6.5.171, SIMPLNX, and a local legacy build with the surgical precision correction. The default symmetric kernel isolates D1 because D2 cannot change the neighbor set.

**Input recipe and seed (for reproducibility):** the archived generator uses `np.random.default_rng(1613)` on a 12×12×12 grid partitioned into eight 6×6×6 octant features. Each feature receives a random-axis 5–25° base rotation; each cell receives an independent ≤3° perturbation. The archive also records a Small IN100 crop and includes the scripts, inputs, and hashes required to reproduce both datasets.

| Metric | Value |
|---|---|
| legacy KAM range (min/mean/max) | 0.951784 / 2.141960 / 3.549188° |
| nx KAM range (min/mean/max) | 0.952303 / 2.141357 / 3.549206° |
| \|Δ\| min / mean / max | 0 / 7.3522e-4 / **7.1352e-3°** |
| cells \|Δ\| > 0.001° | 462 / 1728 |
| cells \|Δ\| > 0.01° | **0** / 1728 |
| signed (legacy − nx): cells legacy>nx / legacy<nx / equal | 907 / 820 / 1 |

Interpretation: the delta is entirely precision-class and is fully explained by D1's family. **Gating is provably identical** on this path — legacy line 292 (`m_FeatureIds[point] == m_FeatureIds[neighbor]`) and SIMPLNX's `use_feature_ids=true` branch admit the same neighbor set for every focal cell, and both include the focal self, so `numVoxel` (the divisor) is identical per cell in both builds. The remaining difference is therefore purely in the per-pair `calculateMisorientation` values, from two combined precision effects: (a) the EbsdLib 2.4.1 symmetry-reduction fix on the focal self-misorientation term (the effect characterized above), and (b) the `QuatF`→`QuatD` port delta — legacy does the misorientation math in `float32`, SIMPLNX in `float64`. Effect (b) is why the delta is **bidirectional** here whereas the earlier pure-φ1 empirical confirmation (2026-06-04) saw legacy ≥ nx: those fixtures used high-symmetry pure-z-axis rotations for which distinct-pair misorientations happen to agree between the two forms to `<1e-4°`, isolating the one-directional self-miso term; on **general 3D orientations** the `float32`-vs-`float64` distinct-pair difference surfaces at the `~1e-3°` scale and takes either sign. Both effects are precision, not algorithmic — no cell exceeds `0.01°` (well below EBSD angular resolution) and there is no structural/gating pattern (a gating difference would show as `O(degrees)` jumps on specific cells, not uniform sub-`0.01°` noise). **Recommendation stands: trust SIMPLNX.**

The complete A/B record is retained in the filter verification archive: reproducible generators, six NX pipelines, twelve legacy pipelines, both shared inputs, all outputs, execution logs, environment details, and per-case comparisons. The archive supersedes the earlier ephemeral scratch record.

---

## ComputeKernelAvgMisorientationsFilter-D2

| Field            | Value                                                                              |
|------------------|------------------------------------------------------------------------------------|
| **Deviation ID** | `ComputeKernelAvgMisorientationsFilter-D2`                                         |
| **Filter UUID**  | `61cfc9c1-aa0e-452b-b9ef-d3b9e6268035`                                             |
| **Status**       | active              |

**Symptom:** Per-cell `KernelAverageMisorientations` values differ between SIMPLNX and DREAM3D 6.5.171 whenever the user-supplied `KernelSize` has `KernelSize.x != KernelSize.z`. For symmetric kernels (`{1,1,1}`, `{2,2,2}`, etc. — the default and the most common use), the deviation is **dormant**. For asymmetric kernels (e.g., `{1, 1, 2}` — common when the user is processing serial-section data with non-isotropic voxel spacing), the legacy code iterates the x-direction inner loop with the WRONG bound, producing a kernel of incorrect shape and an incorrect KAM.

Concrete examples: with `KernelSize = {1,1,2}`, legacy iterates `l = -1..2` (four columns) instead of `-1..1` (three), adding an unintended `x+2` column. With `{2,1,1}`, it iterates `l = -2..1` (four columns) instead of `-2..2` (five), dropping the intended `x+2` column. The same one-character upper-bound correction fixes both directions.

**Root cause:** **Bug** in legacy DREAM3D 6.5.171 only.

The legacy code at `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindKernelAvgMisorientations.cpp:264` is:

```cpp
for(int32_t l = -m_KernelSize.x; l < m_KernelSize.z + 1; l++)
//                                            ^ should be .x
```

The two surrounding outer loops use the correct axis: `m_KernelSize.z` for `j` (line 258) and `m_KernelSize.y` for `k` (line 261). Line 264 is a copy-paste typo where the upper bound `m_KernelSize.z + 1` was carried over from the z-loop instead of being changed to `m_KernelSize.x + 1`.

The SIMPLNX algorithm at `Algorithms/ComputeKernelAvgMisorientations.cpp:111` is correct:

```cpp
for(int32_t l = -kernelSize[0]; l < kernelSize[0] + 1; l++)
```

where `kernelSize[0]` is X. The port from legacy to SIMPLNX silently corrected the bug — most likely the porter manually wrote the loop bound instead of mechanically copy-pasting the legacy line, and used `kernelSize[0]` consistently for both the lower and upper bounds.

**Why this bug went undetected in 6.5.171:** Default and most shipping pipelines use symmetric kernels (`{1,1,1}` is the parameter default; the Small_IN100 reference pipelines all use `{1,1,1}`). The bug is dormant for any symmetric kernel and produces correct output. Asymmetric kernels are uncommon in published DREAM3D workflows but are a real use case for non-isotropic-voxel-spacing serial-section EBSD data.

**Affected users:** DREAM3D 6.5.171 users who ran `FindKernelAvgMisorientations` with an asymmetric `KernelSize`. The output is silently wrong: cells near the upper x-boundary may also see different in-kernel neighbor counts than expected due to the boundary clamp now interacting with the wider-than-requested x-iteration range.

**Recommendation:** **Trust SIMPLNX.** The bug was fixed at port time and SIMPLNX has produced the correct kernel shape for all kernel parameters since the OrientationAnalysis plugin was first ported. Users migrating from DREAM3D 6.5.171 with asymmetric kernels should expect KAM values to change toward the mathematically correct (intended-kernel) value.

The D2 root cause was proven empirically by applying the one-character upper-bound correction to a local build of the legacy source. On both asymmetric kernels and both inputs, the corrected build reduces the SIMPLNX difference to at most 9.54e-7° and zero cells above 0.01°.

This bug is documented in the internal V&V triage record as a known DREAM3D 6.5.171 issue with no SIMPLNX-side action required.

---

## ComputeKernelAvgMisorientationsFilter-D3

| Field            | Value                                                                              |
|------------------|------------------------------------------------------------------------------------|
| **Deviation ID** | `ComputeKernelAvgMisorientationsFilter-D3`                                         |
| **Filter UUID**  | `61cfc9c1-aa0e-452b-b9ef-d3b9e6268035`                                             |
| **Status**       | active                                 |

**Symptom:** SIMPLNX exposes a `use_feature_ids` boolean parameter (default `true`) that DREAM3D 6.5.171 `FindKernelAvgMisorientations` does not have. With `use_feature_ids = false`, SIMPLNX computes a **per-voxel** Kernel Average Misorientation in which a kernel neighbor contributes whenever it is in-bounds, has `featureId > 0`, and shares the focal cell's phase — regardless of whether it belongs to the same feature. There is no way to produce this output with DREAM3D 6.5.171, which only ever computes the per-grain KAM (neighbor must share the focal cell's `featureId`).

**Root cause:** **Algorithmic choice** — a deliberate feature addition (issue #1613), not a bug, precision effect, or library difference. SIMPLNX adds a second neighbor-inclusion mode; the legacy filter has only the per-grain mode.

The two modes differ only in the neighbor gate at `Algorithms/ComputeKernelAvgMisorientations.cpp:124`:

```cpp
const bool neighborContributes = useFeatureIds
    ? (featureIds[point] == featureIds[neighborIdx])                                  // per-grain (legacy-equivalent)
    : (featureIds[neighborIdx] > 0 && cellPhases[neighborIdx] == cellPhases[point]);  // per-voxel (NX-only, #1613)
```

The focal-validity gate (`featureIds[point] > 0 && cellPhases[point] > 0`), the boundary clamps, the divisor semantics (focal self always included), and the background short-circuit are all identical between the two modes and unchanged from the legacy behavior.

**Relationship to the default path:** `use_feature_ids = true` is the default and is behavior-identical to every prior SIMPLNX release and to DREAM3D 6.5.171 (up to the D1/D2 precision/bug notes). The per-voxel mode is strictly opt-in; enabling it cannot change the default output. On single-feature single-phase data the two modes are provably equivalent (see the Class 4 mode-equivalence invariant in the V&V report), because every neighbor that passes the per-grain gate also passes the per-voxel gate and vice-versa.

**Affected users:** none in the migration sense — this is additive. Users who want per-voxel KAM (e.g. to visualize sub-grain orientation gradients without feature segmentation, or to include grain-boundary-adjacent lattice curvature) now have it in SIMPLNX with no DREAM3D 6.5.171 equivalent. Users reproducing legacy pipelines leave the parameter at its `true` default and see no change.

**Validation:** Class 1 analytical fixtures (`Class 1 - Per-Voxel Mode`, expected `{5.0, 20/3, 10.0, 10.0, 0, 0}`; `Class 1 - Per-Voxel Mode Two-Phase Gates`, expected `{5.0, 5.0, 0, 0, 0}`) plus the Class 4 `Mode Equivalence on Single Feature` invariant. Because there is no legacy counterpart, this mode is **never** validated by an A/B numeric comparison — the oracle is the sole authority, per V&V policy.

**Recommendation:** **Trust SIMPLNX.** This is a new, oracle-validated capability. No legacy-parity concern applies.

---

## Non-deviations (algorithm characteristics common to both filters)

The following behaviors are NOT deviations — SIMPLNX (post-EbsdLib 2.4.1) and DREAM3D 6.5.171 (with D2 dormant on symmetric kernels) agree on them where D1 precision noise is below the user's tolerance. Captured here so future engineers don't re-discover them and propose them as deviations.

### Focal voxel always included in the kernel sum

Both implementations have the focal cell as a same-feature neighbor of itself (the `j=k=l=0` inner iteration produces `neighbor = point`). The focal cell's self-misorientation contributes 0° to `totalMisorientation` and 1 to `numVoxel`. This is intentional algorithm characteristic — it provides a non-zero divisor for cells with no in-kernel same-feature neighbors (a single-voxel isolated grain). **Both filters share this behavior** — algorithm characteristic, not a defect.

### Background cell short-circuit to KAM = 0

Both implementations short-circuit background cells to `KAM = 0` at the *end* of the per-cell processing. Legacy explicitly checks `featureIds[point] == 0 || cellPhases[point] == 0` (line 311). **Since the review-driven cleanup (commit `7f9cddc7d`, 2026-07-16), SIMPLNX expresses the identical condition as an `else` branch** (`.cpp:144-146`) on the preceding `if(featureIds[point] > 0 && cellPhases[point] > 0)` gate (`.cpp:86`) rather than a second explicit `if` (previously at the pre-cleanup line 146) — logically identical by De Morgan's law for these `int32` arrays whose valid domain is non-negative (the only way to fail `> 0` is to equal `0`). For hypothetical out-of-domain negative ids the pre-cleanup code left the output cell holding uninitialized memory (the `DataStore` allocates without a fill value); the `else` now writes a deterministic `0.0f` there — a strict improvement over undefined behavior, not a preserved-behavior change. The logic is AFTER the kernel loop but only fires when the focal validity check at the top failed (so the kernel loop didn't run). **Both filters share this behavior**; the SIMPLNX code shape changed, the observable output did not (confirmed by the unchanged 9-test suite).

### `numVoxel == 0` fallback (dead code — removed from SIMPLNX 2026-07-16)

Legacy `FindKernelAvgMisorientations` still includes an `if(numVoxel == 0) { KAM[point] = 0; }` guard immediately after the `KAM[point] = totalMiso / numVoxel` divide. Prior to the review-driven cleanup, SIMPLNX carried the identical dead-code guard (the former path 9 in the V&V report's code path coverage). In practice the focal voxel always self-matches (former path 6), so `numVoxel >= 1` whenever the focal-validity gate passed — the fallback was provably unreachable in both implementations. **Commit `7f9cddc7d` (2026-07-16) removed the SIMPLNX-side guard as pure dead code** (behavior-preserving — confirmed by the unchanged 9-test suite passing); legacy retains its copy unchanged. This is now a code-shape difference only, not a behavioral one: both implementations compute the identical KAM for every reachable input.

### Multi-threading model

Both implementations parallelize over the outer cell loop. Legacy uses `tbb::parallel_for` over a single dimension after marshalling; SIMPLNX uses `ParallelData3DAlgorithm` with a 3D `Range3D`. Both make concurrent reads of the shared input `DataArray`s (FeatureIds, CellPhases, Quats, CrystalStructures). Per the SIMPLNX project policy (`CLAUDE.md`), DataArray subscript access is not formally thread-safe for concurrent reads, but in practice this works for read-only access on contiguous in-memory DataStores. The algorithm has been stable under parallel execution on shipping pipelines; no thread-safety issue surfaced during the V&V cycle's test suite. This filter explicitly calls `parallelAlgorithm.requireArraysInMemory(algArrays)` at line 205 (pre-cleanup: line 196), which disables parallelization when any input array is not resident in memory (the algorithm then runs serially rather than concurrently accessing a store that is not thread-safe for that access pattern); such inputs are still processed correctly — just single-threaded. **Since the review-driven cleanup (commit `7f9cddc7d`), `FindKernelAvgMisorientationsImpl` also builds its `LaueOps` list once in its constructor (a worker member, `.cpp:29,166`) rather than once per `convert()` call — the object is constructed once per `parallelAlgorithm.execute(...)` call and its `m_OrientationOps` member is read-only for the remainder of execution, then shared by const-reference across all parallel `convert()` invocations. This follows the same constructor-built/shared-across-parallel-ranges precedent as `ComputeFeatureFaceMisorientationPerTriangleImpl`, and `calculateMisorientation` remains a const, stateless call on the shared `LaueOps::Pointer` objects — no new thread-safety exposure. Confirmed by the unchanged 9-test suite passing.**
