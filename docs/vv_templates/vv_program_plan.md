# DREAM3D-NX Fleet V&V Program Plan

Companion to [`vv_policy.md`](./vv_policy.md). The policy says *how* to V&V one filter.
This says *which* filters, *in what order*, and *at what rigor*.

**Status:** DRAFT — updated 2026-08-10
**Baseline measured at:** `upstream/develop` @ `165f8713a`

---

## 1. The situation

| | count |
|---|---|
| Filters declaring `SIMPLNX_DEF_FILTER_TRAITS` | 302 |
| — ITKImageProcessing (slated for removal) | 88 |
| — **real surface** | **214** |
| V&V reports authored (all branches) | 40 |
| — of those, present on `develop` | **34** |
| Filters used in the 68 shipped example pipelines | 115 (114 non-ITK) |
| Filters pinned to a declared `6_6_*.tar.gz` archive | 36 (29 not yet authored) |
| Distinct `6_6_*` archives referenced by tests | 21 |

**Measured throughput:** 25 reports between 2026-05-27 and 2026-07-24 at full v2 rigor
≈ **3 filters/week**. Verifying all 214 at that rigor is ~14 months of continuous effort.
That is not the plan.

### 1.1 The FAA obligation is closed on paper, not on `develop`

The regulated artifact is the MTR analysis pipeline, not the filter fleet. The template
routine `PW_ANG_routine_v75.d3dpipeline` is **51 steps / 28 unique enabled filters**.
Including the CTF variant of the routine, the full closure is **29 filters**.

Both remaining gaps were authored 2026-07-28 by @JDuffeyBQ and are now merged —
`ITKImageWriterFilter` (#1693) and `RequireMinNumNeighbors` (#1694).

**All 29 closure reports are authored. `develop` carries 25/29 (86%).**
The remaining four are a review problem, not an engineering one — see §6.

#### Live risk: the `ITKImageWriterFilter` report is a decaying asset

PR #1693 merged as-is, which settles the immediate question — `PW_ANG_routine_v75` calls
`ITKImageWriterFilter` today, so it is a genuine closure filter and the closure counts it.
The underlying tension is unresolved, only deferred:

- **If the submission tags before ITK removal**, nothing further is needed.
- **If ITK is removed first**, the report is deleted along with the plugin,
  `src/Plugins/ITKImageProcessing/vv/` (currently one file) disappears, and
  `WriteImageFilter` reopens as a closure gap that nobody is tracking.

Whoever executes the ITK removal must check the closure before deleting the plugin.
Encode this in the Phase 0 CI gate (§3.1): removing a plugin that contains a closure
filter's V&V report should fail the build, not silently reduce coverage.

Mitigation if it goes that way: the oracle design transfers. The `x + 10y + 100z` fixture
over XY/XZ/YZ planes is format-agnostic and carries over to `WriteImageFilter` despite the
different backends (stb/libtiff vs ITK) — the re-spend is plumbing, not analysis.

The legacy CTF routine's `ConvertData` and `ChangeAngleRepresentation` steps are **not**
in the closure — both are now handled inside the CTF reader:

- `ChangeAngleRepresentation` (degrees → radians) is subsumed by `ReadCtfDataFilter`'s
  `k_DegreesToRadians_Key` parameter (default `true`); conversion happens in
  `Algorithms/ReadCtfData.cpp:143`. Confirmed in source.
- `ConvertData` converted `BC` → `BC_Float`. `ReadCtfDataFilter` types each column from
  EbsdLib's native CTF column type (`ReadCtfDataFilter.cpp:147-152`), so `BC` still lands
  as `int32`. **Caveat:** this is subsumed only if the NX CTF routine has no downstream
  need for `BC` as `float32` — i.e. if the legacy step was a 6.5-era workaround. Confirm
  when the NX CTF routine is authored (§3.1).

Everything beyond those two is **discretionary product quality**, not compliance.
It is therefore optimized for risk reduction per engineer-hour at a sustainable pace,
not for coverage against a deadline.

### 1.2 `6_6_*` exemplars are a defect backlog, not merely a gap

A test pinned to a `6_6_*` archive is green *because* DREAM3D-NX reproduces
DREAM3D 6.6 output. No release of 6.6 ever shipped and it was never rigorously tested.
If 6.6 was wrong, CI is actively certifying the wrong answer — strictly worse than
having no test at all.

This is not hypothetical. The EbsdLib 2.4.1 `CubicOps` `2·atan2(|v|,w)` precision fix
broke five filters precisely by exposing exemplars that encoded the pre-fix wrong values.

36 filters are pinned to a declared `6_6_*.tar.gz` archive; 29 of them have no V&V report
anywhere. A further ~16 test files reference 6.6-era data without a declared archive.
Hit rate unknown.

---

## 2. Rigor tiers

| Tier | Oracle | Legacy 6.5.171 A/B | Deviation file | Applies to |
|---|---|---|---|---|
| **A** | Class 1–4, second-engineer reviewed | Required | Required | MTR closure (§1.1) |
| **B** | Class 1–4 | **Omitted** unless the filter is on a real 6.5 migration path | Only if A/B was run | MTR-adjacent (§3.4) |
| **C** | Class 4 invariants via shared harness, or re-derived archive | Omitted | No | Long tail (§3.5) |
| **D** | None | — | — | Explicitly registered as unverified (§4) |

**The single largest throughput lever is dropping the legacy A/B below Tier A.**
`vv_policy.md` already states that matching 6.5.171 is a *diff-explanation* check and
never establishes correctness. For a filter no customer migrates from 6.5, the legacy
pipeline authoring, dual runs, discrepancy investigation, and deviation write-up cost
roughly half the per-filter effort and buy zero verification value.

Tier B is expected to run at 5–8 filters/week against Tier A's 3.

---

## 3. Phases

### 3.1 Phase 0 — pin the deliverable

The authoritative definition of the FAA surface currently lives only at
`/Users/Shared/Data/MTR_Data/dot_78565_DS1/Templates/` — outside version control.
That is a single point of failure for the entire SBIR deliverable.

1. Commit the ANG and CTF template routines (paths sanitized) under
   `docs/vv_templates/mtr_closure/`.
2. Add `scripts/vv_status.py`: resolves each pipeline step's UUID against
   `SIMPLNX_DEF_FILTER_TRAITS`, joins to `src/Plugins/*/vv/*.md`, prints coverage
   for the closure, the MTR-adjacent set, and the fleet.
3. CI check: fail if any filter in the MTR closure lacks a V&V report or changes UUID.
4. Normalize the `| Status |` line to a fixed vocabulary —
   `DRAFT` / `READY FOR REVIEW` / `COMPLETE`.
   **Largely done by PR #1701** (approved, pending merge), which flips 26 reports to
   `COMPLETE`. The leading token is now consistent even though free-text suffixes vary
   (`COMPLETE — 2026-07-16`, `COMPLETE — all V&V phases complete; …`), so `vv_status.py`
   should anchor on the leading keyword and treat the tail as human annotation rather
   than forcing the suffixes into a schema.

5. Author the NX CTF routine (no `.d3dpipeline` equivalent of `PW_CTF_routine_v65.json`
   exists yet) and confirm the `ConvertData` caveat in §1.1 while doing so.

`vv_status.py` must report against a **named git ref**, not the working tree, and default
to `develop`. The gap between "all branches" and "`develop`" is currently 11 points of
closure coverage; a tool that hides that gap is worse than no tool.

**Exit:** `vv_status.py` reports the closure at a ref, and CI enforces it.

### 3.2 Phase 1 — close the FAA gap (Tier A)

**Authoring is done — all 29 closure reports exist.** What remains is review and merge:

1. **Clear the four `CHANGES_REQUESTED` closure PRs** — #1688, #1687, #1683, #1672.
   These are now the *only* thing between `develop` and a complete deliverable. Until they
   land, `develop` reports 25/29 and any release tag captures an incomplete evidence
   package. Three of the four belong to one engineer. See §6.
2. **Merge #1701** (`VV: Mark filters 'COMPLETE'`) — approved 2026-08-07 and still
   unmerged. It flips 26 reports from `READY FOR REVIEW` to `COMPLETE` and largely
   satisfies Phase 0 item 4 (§3.1). It also grows a conflict surface against every V&V PR
   that touches a report header, so the longer it sits the more it costs.
3. **Re-examine PR #1640** (`VV/BUG: Identify Sample Full V&V`), merged with
   `CHANGES_REQUESTED` and no approving review. `IdentifySampleFilter` is in the closure,
   so its report is the one merged closure entry without a clean second-engineer sign-off.
4. **Backfill the oracle attestation** (§3.6) across all 29 closure reports. This is the
   audit trail an external reviewer will ask for and it does not exist yet.
5. **Guard the ITK removal** against silently reopening the `ITKImageWriterFilter` gap
   (§1.1).

**Estimate:** review-bound, not effort-bound. **Exit:** `develop` reports 29/29 and every
closure report carries a recorded attestation — at which point the SBIR filter work is
genuinely, verifiably closed.

### 3.3 Phase 2 — `6_6_` defect sweep (Tier C, cheapest-first)

Ordered by **archive**, not by filter — the leverage is uneven. 36 filters are pinned to
21 declared `6_6_*.tar.gz` archives; 29 of those filters have no V&V report anywhere:

| Archive | Filters pinned, not yet authored |
|---|---|
| `6_6_stats_test_v2` | **5** of 9 — AlignSectionsMutualInformation, ComputeEuclideanDistMap, ComputeSchmids, ComputeShapes, ComputeSurfaceAreaToVolume |
| `6_6_Small_IN100_GBCD` | **4** of 5 — ComputeGBCD, ComputeGBCDPoleFigure, WriteGBCDGMTFile, WriteGBCDTriangleData |
| `6_6_combine_stl_files_v2` | 2 — CombineStlFiles, WriteStlFile |
| `6_6_erode_dilate_test` | 2 of 3 — ErodeDilateCoordinationNumber, ErodeDilateMask |
| `6_6_find_feature_centroids` | 2 — ExtractComponentAsArray, WriteAbaqusHexahedron |
| `6_6_avizo_writers` | 2 — WriteAvizoRectilinearCoordinate, WriteAvizoUniformCoordinate |
| remaining 15 archives | 1 each |

**2a — the writers.** `WriteAvizoUniformCoordinate`, `WriteAvizoRectilinearCoordinate`,
`WriteINLFile`, `WriteSPParksSites`, `WriteStlFile`, `WriteLosAlamosFFT`,
`WriteAbaqusHexahedron`, `WriteGBCDGMTFile`, `WriteGBCDTriangleData`.

Every one has a free Class 1 oracle — the published file-format specification. Build a
2×2×2 geometry, hand-write the expected bytes, assert equality. No legacy run, no
exemplar archive, no `6_6_` dependency afterwards. ~9 filters, purges ~12 archive refs.

**2b — simple element-wise.** `ExtractComponentAsArray`, `ComputeVectorColors`,
`ComputeNumFeatures`, `ComputeVolumeFractions`, `ComputeFeaturePhasesBinary`,
`ComputeBoundaryCells`, `ReverseTriangleWinding`. Class 1/4 on hand-built input.

**Deferred as expensive** (revisit in Phase 3 or 4): the GBCD family, the two metric-based
filters, `MapPointCloudToRegularGrid`, `AlignSectionsMutualInformation`,
the surface-meshing group.

**Deliverable is per archive, not per filter:** a re-derived archive plus a provenance
sidecar. Filters pinned to it get their V&V report in Phase 3. The point of this phase
is to stop CI from certifying wrong answers, not to move the report counter.

**Report the hit rate.** Every archive that turns out to encode a wrong value is
evidence for how much the rest of this program is worth. If the rate is near zero,
Phase 3 can be scaled back; if it is high, Phase 3 should be scaled up.

### 3.4 Phase 3 — MTR-adjacent coverage (Tier B)

Filters a Ti supplier could plausibly substitute into or bolt onto the PW routine.
**63 filters total = 29 closure + 34 additions. 29 authored, 34 remaining.**

| Group | Filters |
|---|---|
| Alternate EBSD readers | `ReadH5Ebsd`, `ReadH5OimData`, `ReadH5EspritData`, `ReadH5OinaData`, `EbsdToH5Ebsd`, `ConvertHexGridToSquareGrid` |
| Alignment | `AlignSectionsList`, `AlignSectionsMisorientation`, `AlignSectionsFeatureCentroid`, `AlignSectionsMutualInformation` |
| Segmentation | `EBSDSegmentFeatures`, `ScalarSegmentFeatures`, `MergeTwins` |
| Cleanup | `ErodeDilateMask`, `ErodeDilateCoordinationNumber`, `RequireMinimumSizeFeatures` |
| Feature statistics | `ComputeShapes`, `ComputeSurfaceFeatures`, `ComputeBiasedFeatures`, `ComputeEuclideanDistMap`, `ComputeLargestCrossSections` |
| Data plumbing | `CreateDataArray`, `CreateAttributeMatrix`, `DeleteData`, `MoveData`, `CombineAttributeArrays`, `CropImageGeometry`, `ConditionalSetValue`, `ArrayCalculator`, `ComputeArrayHistogram` |
| Output / reporting | `WriteFeatureDataCSV`, `WriteASCIIData`, `ReadDREAM3D`, `CreateColorMap` |

All 34 additions are name-verified against the source tree. This list is a judgment
call and is expected to be edited; it is committed so the judgment is reviewable.

Order within the phase: filters already touched by Phase 2 first (the oracle work is
half done), then by shipped-pipeline usage count descending.

**First divergence to reconcile:** V&V of `DBSCANFilter` is complete and submitted as
#1702 — the first post-closure work finished, and it is *not* in the list
above. `DBSCAN` appears 6 times in the shipped example pipelines, so it is defensible
Phase 3 work under the exposure criterion, but it entered the queue by engineer choice
rather than by this plan. Either add it to the table or record why it was picked ahead of
the listed filters. If organic selection keeps outrunning the list, that is evidence the
list is wrong, not that the engineers are.

**Estimate:** ~34 filters at 5/week ≈ 7 weeks. **Exit:** MTR-adjacent 63/63.

### 3.5 Phase 4 — long tail (Tier C, open-ended)

The remaining ~151 non-ITK filters. No deadline, no coverage target — this runs
alongside feature work indefinitely.

Build a shared Class 4 invariant harness (`test/UnitTestCommon` additions) asserting
derivable properties across many filters at once:

- FeatureIds are contiguous and start at 1
- Feature-array tuple count equals max(FeatureId) + 1
- Sum of per-feature cell counts equals the number of non-zero cells
- Neighbor lists are symmetric
- Phase fractions sum to 1
- Geometry element counts are consistent after any topology-modifying filter

This is weak per-filter evidence but it is the only mechanism that scales to 151,
and it catches whole classes of breakage that exemplar diffing misses entirely.

### 3.6 Cross-cutting — make the oracle attestation auditable

`vv_policy.md` asks that the oracle and its design be reviewed by a second engineer,
on the grounds that "a wrong oracle silently confirms buggy filters, and the filter
author is the least likely person to notice."

In practice this is delegated to the PR reviewer, and **the mechanism largely works**:
of 43 V&V PRs, the merged ones nearly all carry an approving review from a different
engineer — `nyoungbq` (10), `imikejackson` (5), `mmarineBlueQuartz` (3), `JDuffeyBQ` (3).

Two problems remain:

1. **A PR approval is not an oracle review.** A reviewer approving a diff is reviewing
   code, tests, and documentation. Nothing in the process obliges them to attest that
   they examined the oracle and found it independent of the implementation, and nothing
   records that they did. The policy's central control is therefore unauditable.
2. **At least one merge bypassed review entirely** — PR #1640, merged with
   `CHANGES_REQUESTED` and no approval, for a filter inside the MTR closure.

PR #1701 moves in the right direction — one report now reads
`COMPLETE — … second-engineer sign-off recorded at PR review` — but it names no reviewer
and no PR number, so it asserts the review happened without providing evidence. Prose is
not an audit trail. Now is the cheapest time to fix this, while #1701 is flipping 26
reports at once.

Fix — an attestation, not a new process:

- Add to `report_template.md`: `| Oracle reviewed by | @<handle>, PR #<n> |`
- `vv_status.py` refuses to report `COMPLETE` for any filter whose attestation is empty.
- Add a matching checkbox to the PR template: *"I reviewed the oracle design and confirm
  it is independent of the implementation under test."*
- Require an approving review on branch protection for paths under `src/Plugins/*/vv/`.

This converts an informal practice already being followed into something an external
reviewer can verify, at close to zero ongoing cost.

---

## 4. The unverified register

`docs/vv_templates/unverified.md`, generated by `vv_status.py`: every non-ITK filter
with no V&V report, its tier assignment, and why it has not been done.

For a regulator, an explicit and accurate statement of what is *not* verified is more
credible than an unqualified claim of blanket coverage. This file is a deliverable,
not an embarrassment.

---

## 5. Decisions recorded

| Decision | Rationale |
|---|---|
| ITK's 88 filters are out of scope entirely | Plugin is being removed; the cheapest V&V is a filter that does not ship |
| Legacy 6.5.171 A/B is Tier A only | `vv_policy.md` already holds that legacy match never establishes correctness; it is migration guidance, and only migrating filters need it |
| Phase 2 unit of work is the archive, not the filter | `6_6_stats_test_v2` alone backs 11 filters; re-deriving one archive clears all of them |
| Phase 3 priority set is a judgment call, committed as a list | A mechanical rule (e.g. "all filters in shipped example pipelines") measures what BlueQuartz chose to demo, not what MTR work needs |
| No fleet-wide coverage target | 214 filters at full rigor is ~14 months; a target nobody can hit produces either burnout or a diluted definition of "verified" |
| Second-engineer oracle review stays delegated to the PR reviewer, made auditable via a recorded attestation | The practice already exists and works; what is missing is evidence of it, not the review itself |

## 6. The merge backlog is the critical path

**40 V&V reports have been authored. 34 are on `develop`.** The other 6 exist only on
unmerged branches, and **4 of those 6 are MTR closure filters — every one stalled in a
review round-trip**:

| Filter | PR | Blocked on |
|---|---|---|
| `MultiThresholdObjects` | #1688 | `CHANGES_REQUESTED` (@mmarineBlueQuartz) |
| `ErodeDilateBadData` | #1687 | `CHANGES_REQUESTED` (@mmarineBlueQuartz) |
| `WriteDREAM3D` | #1683 | `CHANGES_REQUESTED` (@mmarineBlueQuartz) |
| `ComputeFeaturePhases` | #1672 | `CHANGES_REQUESTED` (@nyoungbq) |

The remaining two are non-closure: `ComputeTriangleGeomCentroids` and `DBSCANFilter`
(#1702, awaiting first review).

**Progress:** `develop` went 18/29 → 23/29 → 24/29 → **25/29** across 2026-07-28 to
2026-08-10. Merged in that window: #1693, #1692, #1689, #1685, #1638, #1694, #1695, #1679.

**The backlog is now one engineer deep.** Every approved or awaiting-review closure PR has
landed; all four that remain are waiting on an author to answer review comments, and
**three of the four belong to @mmarineBlueQuartz** (#1688, #1687, #1683). The critical
path to the SBIR deliverable is one person's revision queue. That is worth knowing
explicitly — it is a scheduling risk, not a process failure, and it is addressed by
reassignment or by prioritizing those three above other work, not by more V&V.

The consequence is unchanged in kind, only in size: **the SBIR filter work is finished as
engineering and 86% complete as evidence.** A release tag cut today would pin
`(commit hash, archive SHA512)` pairs — the mechanism `vv_policy.md` relies on — against a
`develop` missing four of the twenty-nine closure reports.

No amount of new V&V improves this number. Clearing the backlog is the highest-value work
available until `develop` reports 29/29.

## 7. Open questions

- **Does ITK removal land before the submission tag?** (§1.1) — #1693 merged, so this no
  longer gates anything today, but removing the plugin silently reopens a closure gap.
  The Phase 0 CI gate must catch it.
- Are there customer-authored MTR pipelines (beyond the PW templates) that should
  expand the §1.1 closure? If so they must be committed alongside the templates.
- Does the FAA submission pin a DREAM3D-NX version? If so, the closure's
  (commit hash, archive SHA512) pairs need to be captured at that tag — and §6 must be
  resolved before that tag is cut.
- Does the NX CTF routine need `BC` as `float32`? Determines whether `ConvertDataFilter`
  re-enters the closure (§1.1).
