# DREAM3D-NX Fleet V&V Program Plan

Companion to [`vv_policy.md`](./vv_policy.md). The policy says *how* to V&V one filter.
This says *which* filters, *in what order*, and *at what rigor*.

**Status:** DRAFT — updated 2026-08-12
**Baseline measured at:** `upstream/develop` @ `bbc470e6d`

---

## 1. The situation

| | count |
|---|---|
| Filters declaring `SIMPLNX_DEF_FILTER_TRAITS` | 302 |
| — ITKImageProcessing (slated for removal) | 88 |
| — **real surface** | **214** |
| V&V reports authored (all branches) | 40 |
| — of those, present on `develop` | **35** |
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

**All 29 closure reports are authored. `develop` carries 27/29 (93%).**
The remaining two are a review problem, not an engineering one — and they both belong
to one engineer. See §6.

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
4. **Normalize the `| Status |` line — and then *enforce* it.** PR #1701 (merged
   2026-08-10) swept all reports to `DRAFT` / `READY FOR REVIEW` / `COMPLETE`.
   **It regressed on the very next merge.** `ComputeFeaturePhasesFilter.md` (#1672,
   2026-08-12) writes `| **Status** | COMPLETE |` with the key bold-wrapped, which does
   not match `^\| Status \|` — so 35 reports are on `develop` but only 34 are countable.

   The lesson is that a one-time sweep does not hold. `vv_status.py` must (a) tolerate
   optional `**` around the key, and (b) **fail** on any `vv/*.md` it cannot parse, rather
   than silently omitting it. An unparseable report is a reporting bug, not a zero.
   Free-text suffixes after the keyword remain fine as human annotation.

5. Author the NX CTF routine (no `.d3dpipeline` equivalent of `PW_CTF_routine_v65.json`
   exists yet) and confirm the `ConvertData` caveat in §1.1 while doing so.

`vv_status.py` must report against a **named git ref**, not the working tree, and default
to `develop`. The gap between "all branches" and "`develop`" is currently 3 points of
closure coverage — down from 11 — and a tool that hides that gap is worse than no tool.

**Exit:** `vv_status.py` reports the closure at a ref, and CI enforces it.

### 3.2 Phase 1 — close the FAA gap (Tier A)

**Authoring is done — all 29 closure reports exist.** What remains is review and merge:

1. **Clear the two remaining `CHANGES_REQUESTED` closure PRs** — #1688, #1683.
   (#1687, `ErodeDilateBadData`, merged 2026-08-18.) These are the *only* thing between
   `develop` and a complete deliverable, and **both belong to @mmarineBlueQuartz.** Until
   they land, `develop` reports 27/29 and any release tag captures an incomplete evidence
   package. See §6.
2. **Land the two open non-closure PRs** — #1702 (`DBSCAN`) and #1703
   (`GroupMicroTextureRegions`), both awaiting first review. Neither affects the closure,
   but #1703 closes the last outstanding `DRAFT` gate items from #1637.
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

### 3.3 Phase 2 — MTR-adjacent coverage (Tier B)

Filters a Ti supplier could plausibly substitute into or bolt onto the PW routine.
**64 filters total = 29 closure + 35 additions. 29 authored, 35 remaining.**

> **Reordered 2026-08-18.** This was Phase 3 and now runs ahead of the `6_6_` sweep.
> Rationale: 11 of these filters are themselves `6_6_`-pinned, so this phase delivers the
> defect sweep's measurement as a by-product of work that also carries MTR exposure. The
> `6_6_`-pinned filters left behind are almost entirely niche exporters and specialty
> analyses no Ti supplier runs. The writers plan at
> `docs/vv_templates/plans/2026-08-18-phase2a-writer-oracles.md` targets what is now Phase 3.

| Group | Filters |
|---|---|
| Alternate EBSD readers | `ReadH5Ebsd`, `ReadH5OimData`, `ReadH5EspritData`, `ReadH5OinaData`, `EbsdToH5Ebsd`, `ConvertHexGridToSquareGrid` |
| Alignment | `AlignSectionsList`, `AlignSectionsMisorientation`, `AlignSectionsFeatureCentroid` |
| Segmentation | `EBSDSegmentFeatures`, `ScalarSegmentFeatures`, `MergeTwins` |
| Cleanup | `ErodeDilateMask`, `ErodeDilateCoordinationNumber`, `RequireMinimumSizeFeatures` |
| Feature statistics | `ComputeShapes`, `ComputeSurfaceFeatures`, `ComputeBiasedFeatures`, `ComputeEuclideanDistMap`, `ComputeLargestCrossSections`, `ComputeSchmids`, `ComputeSurfaceAreaToVolume` |
| Data plumbing | `CreateDataArray`, `CreateAttributeMatrix`, `DeleteData`, `MoveData`, `CombineAttributeArrays`, `CropImageGeometry`, `ConditionalSetValue`, `ArrayCalculator`, `ComputeArrayHistogram` |
| Output / reporting | `WriteFeatureDataCSV`, `WriteASCIIData`, `ReadDREAM3D`, `CreateColorMap` |

All 35 additions are name-verified against the source tree. This list is a judgment
call and is expected to be edited; it is committed so the judgment is reviewable.

**Roster changes 2026-08-18:**
- **Removed `AlignSectionsMutualInformation`** — expensive to verify and not part of any
  plausible Ti-supplier workflow.
- **Added `ComputeSchmids` and `ComputeSurfaceAreaToVolume`** — both plausible for MTR
  reporting, and both pinned to `6_6_stats_test_v2`, the highest-leverage archive in the
  repo (9 filters). **Note:** that archive does *not* fully retire, because
  `AlignSectionsMutualInformation` also pins it and is now out of scope. Retiring
  `6_6_stats_test_v2` requires either verifying that filter later or converting its test.

### Opening batch — the 11 dual-purpose filters

Sequence these first. Each carries MTR exposure *and* is pinned to a `6_6_` archive, so
each one advances the SBIR story and drains the defect backlog at the same time (the
Cleanup group's three filters completed 2026-08-18, branches pending merge; the
`ErodeDilateBadData` conversion completed the same day and rides its own branch —
the filter's V&V itself merged as #1687):

`AlignSectionsFeatureCentroid`, `AlignSectionsMisorientation`, `ComputeBiasedFeatures`,
`ComputeEuclideanDistMap`, `ComputeShapes`, `ComputeSchmids`, `ComputeSurfaceAreaToVolume`,
`ErodeDilateCoordinationNumber`, `ErodeDilateMask`, `ReadH5OinaData`,
`RequireMinimumSizeFeatures`

Archives fully retired by this batch (target outcome): `6_6_ImportH5Data`,
`6_6_erode_dilate_test`, `6_6_find_biased_features`, `6_6_min_size_input`,
`6_6_min_size_output`. (Status 2026-08-18: the two `6_6_min_size_*` archives are retired
on the `RequireMinimumSizeFeatures` branch; `6_6_erode_dilate_test` retires once the last
of its three consumer PRs merges — see the hit-rate note below; the rest await their
filters.)

**Report the hit rate from this batch.** It is the measurement Phase 3 was originally
meant to produce, and it now arrives earlier and attached to filters that matter.

**Hit rate, first three of the batch (2026-08-18).** `ErodeDilateMask`,
`ErodeDilateCoordinationNumber`, and `RequireMinimumSizeFeatures` are the first three of
the eleven to finish: Class 1 (+Class 4 for `RequireMinimumSizeFeatures`) oracle V&V plus a
full 6.5.171 binary A/B. Verdict B (SIMPLNX defect): **1 of 3** — `ErodeDilateMask`
accepted its `X`/`Y`/`Z Direction` flags but never read them, a port regression (6.5.171
honoured them), fixed in-PR. That is the same defect class as
`ErodeDilateBadDataFilter-D1` from #1687 — two independent instances of one
copy-paste-era defect family across the erode/dilate group, worth flagging even though
neither gets a deviation entry (port regressions get none, per policy). Legacy-shared
output bugs: **0 of 3** — zero deviation entries, no 6.5.172 surgical patch triggered by
this batch. Post-fix A/B parity was exact: `ErodeDilateMask` 28/28 pipeline pairs,
`ErodeDilateCoordinationNumber` 44/44 array pairs, `RequireMinimumSizeFeatures` 32/32
array pairs, all element-wise identical. **A truthful low number is a useful result**: one
port regression and zero legacy-shared bugs across three filters is the measurement, and
here the A/B's marginal contribution was catching that port regression and independently
confirming all five hand oracles, not finding output bugs. The work also surfaced findings
outside the hit-rate definition proper: one behavioural deviation filed
(`ErodeDilateMaskFilter-D1` — legacy errors out with `-5555` when `NumIterations <= 0`, NX
silently no-ops; recommendation pending a product decision), one NX-only defect fixed in
`ErodeDilateCoordinationNumber` (the cancel flag was plumbed but never read, leaving the
filter uncancellable), and two legacy-shared hazards documented but not fixed
(`ErodeDilateCoordinationNumber`'s `Loop=true, CN<=1` non-termination, and
`RequireMinimumSizeFeatures`'s unguarded feature-id indexing / no-progress loop — both
verified present in 6.5.171 source). Archives: `6_6_min_size_input` and
`6_6_min_size_output` retire with the `RequireMinimumSizeFeatures` commit;
`6_6_erode_dilate_test` has all three consumer conversions complete on their respective
branches, all pending merge (BadData's underlying V&V already merged as #1687; its
archive-test conversion is a separate follow-on branch), and retires via a
`CMakeLists.txt` removal once the last of the three PRs merges.

**Hit rate, feature-statistics five (2026-08-20).** `ComputeBiasedFeatures`,
`ComputeSurfaceAreaToVolume`, `ComputeEuclideanDistMap`, `ComputeShapes` and
`ComputeSchmids` finished the Feature-statistics group: Class 1 (+Class 4) oracles plus
full 6.5.171 binary A/B, with every predicted divergence derived from source before any
run (11/11, 10/10, 87/87, and 93/93 A/B rows matched prediction across the batch — zero
unpredicted differences). Verdict B (SIMPLNX output defect, fixed in-PR): **3 of 5** —
`ComputeBiasedFeatures` (2D spacing remap port regression, plus an uninitialized-output
fill), `ComputeSurfaceAreaToVolume` (shared ±X/±Y face-area swap and truncated sphericity
exponents), `ComputeShapes` (2D degenerate-Euler `π/180`-for-`π/2` mis-transcription —
unreachable until the shared 2D voxel-center bug was also fixed — plus the 2D axis remap
and an empty-feature NaN); `ComputeEuclideanDistMap` needed **no output change**, as
predicted, and `ComputeSchmids`' suspected uninitialized-output defect was **overturned by
experiment** into a latent-only regression (the in-core store factory hard-codes zero
init). Legacy-shared bugs: **three alignment patches** on the 6.5.172 branch
(`FindBoundingBoxFeatures`, `FindSurfaceAreaToVolume`, `FindShapes`), each proving
patched-legacy ≡ fixed-NX ≡ oracle. Legacy-only bugs newly documented: six — including
`FindEuclideanDistMap`'s zero-init corrupting TJ-only/QP-only maps in both distance modes
(release-note impact) and `FindShapes`' voxel-corner moment bias. **EbsdLib itself was
fixed** on `topic/3_1_1_staging` (Schmid normalizer constants that let m exceed 0.5;
uninitialized Hexagonal-Low outputs), with twelve further hexagonal sqrt-divisor literals
recorded as known-open follow-up; the `ComputeSchmids` PR is merge-blocked until EbsdLib
3.1.1 ships. Archives: `6_6_find_biased_features` retired outright;
`6_6_stats_test_v2`'s live consumers drop from six to two.

**Hit rate, align-sections pair (2026-08-21).** `AlignSectionsFeatureCentroid` and
`AlignSectionsMisorientation`: Class 1 oracles plus full 6.5.171 binary A/B, every divergence
predicted from source before any run (59 comparables against 6.5.171 with 21 predicted
divergences for FeatureCentroid; 95 checks all exact-value equal for Misorientation — zero
unpredicted differences across the pair). Verdict B (SIMPLNX output defect, fixed in-PR):
**1 of 2** — `AlignSectionsFeatureCentroid` corrected six defects, headlined by the
*Reference Slice* semantics bug **shared with 6.5.171** (the parameter indexed the internal
top-down iteration order, so "slice k" anchored physical slice `Z-1-k`; per product ruling
slice 0 is now the slice at the Z origin — a user-visible behavior change for any pipeline
with *Use Reference Slice* on, release-note required), plus the shared all-masked-slice
NaN→int64 UB, the shared/NX bounds guard, uninitialized shift-array tuple 0, demoted-and-broken
diagnostics, and four further missing preflight guards. One 6.5.172 alignment patch proves
patched-legacy ≡ fixed-NX ≡ oracle on the divergence fixtures. `AlignSectionsMisorientation`
needed **no output change**, as predicted — its yield was five new guards closing OOB/UB
acceptance paths, with legacy parity exact. Legacy-only bugs newly documented: three
(FeatureCentroid's always-zero "New X/Y Shift" CSV columns, its `>` reference-bound
off-by-one, and iteration-index diagnostics). New shared-code findings for the follow-up
queue: the `AlignSections` base silently discards `findShifts` warnings (all four consumers)
and still crashes the two out-of-batch filters on non-`IDataArray` cell children; guard
parity for `AlignSectionsList`/`AlignSectionsMutualInformation` is assigned to the engineer
track. Archives: the orphaned `6_6_align_sections_feature_centroids.tar.gz` download and the
vestigial `align_sections.tar.gz` retire with these PRs; both filters' circular
`output_*.dream3d` shift-array comparisons were replaced by hand-derived assertions
(`AlignSectionsListTest` still consumes the Misorientation output exemplar as a golden input —
repo-wide circularity closure is a recorded follow-up).

**Hit rate, H5OINA reader (2026-08-24).** `ReadH5OinaDataFilter` is the program's first
new-filter V&V (no 6.5.171 equivalent exists — verified against the whole legacy tree), so
the legacy binary A/B is replaced by a Class 2 independent readback: an h5py script re-reads
every fixture and the real AZtec file and independently derives the expected NX arrays.
Verdict B: **11 defects corrected — 6 in the filter, 5 in EbsdLib's `H5OinaReader`** (on
`topic/3_1_1_staging`). Filter side: hexagonal alignment added the literal `30.0` (thirty
*radians*) to radian-valued φ2 with the option on by default — garbage orientations for
every hexagonal point in every release that shipped the filter, invisible because the only
test file was cubic; two multi-scan stacking corruptions (Euler slab offset off by 3×;
hex alignment applied S times to scan 0 and never to the rest); a process crash reachable
from a well-formed multi-scan file (ensemble arrays sized from the first scan only); pattern
import descoped behind an honest unsupported error; stacking order actually implemented
(previously accepted and ignored). EbsdLib side: the lattice-constant γ slot echoed β; lattice
angles were imported in radians where every other importer reports degrees (**breaking change
to a published output, release-noted**); a phase group missing `Lattice Angles` crashed the
process. The circular NX-generated exemplar comparison was replaced by readback-derived
assertions (2 → 17 test cases). **Release-gating note:** the EbsdLib 3.1.1 release now gates
two PRs (`ComputeSchmids`, this one); measured at this batch's head, staging EbsdLib breaks
exactly **one** test on develop (`ComputeSchmidsFilter`, the known CubicOps drift, corrected
by the open `ComputeSchmids` PR) — an earlier 29-failure count was a stale-binary artifact
(`PIPELINE::`/`PY::` tests need `--target all`, not just the unit-test target). Sibling
exposure (`ReadH5OimData`, `ReadH5EspritData` share several defect shapes through
`IEbsdOemReader`) is enumerated in the report's follow-ups for the engineer track.

After the opening batch, order by shipped-pipeline usage count descending.

**Throughput caveat:** these are harder than the deferred writers. Do not schedule the
opening batch against the 5–8 filters/week Tier B figure — the alignment and feature-stats
filters need realistic EBSD-shaped fixtures and will run closer to Tier A rates.

**First divergence to reconcile:** V&V of `DBSCANFilter` is complete and submitted as
#1702 — the first post-closure work finished, and it is *not* in the list
above. `DBSCAN` appears 6 times in the shipped example pipelines, so it is defensible
Phase 3 work under the exposure criterion, but it entered the queue by engineer choice
rather than by this plan. Either add it to the table or record why it was picked ahead of
the listed filters. If organic selection keeps outrunning the list, that is evidence the
list is wrong, not that the engineers are.

**Estimate:** ~35 filters; the 11-filter opening batch at Tier A-ish rates, the remainder
nearer Tier B. **Exit:** MTR-adjacent 64/64.

### 3.4 Phase 3 — `6_6_` defect sweep (Tier C, cheapest-first)

Ordered by **archive**, not by filter — the leverage is uneven. 36 filters are pinned to
21 declared `6_6_*.tar.gz` archives; 29 of those filters have no V&V report anywhere:

| Archive | Filters pinned, not yet authored |
|---|---|
| `6_6_stats_test_v2` | **1** of 9 — AlignSectionsMutualInformation (the other four were authored 2026-08-20, branches pending merge; live test consumers after that batch: AlignSectionsMutualInformation and ComputeFeatureNeighbors) |
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

**Deferred as expensive** (revisit in Phase 4): the GBCD family, the two metric-based
filters, `MapPointCloudToRegularGrid`, `AlignSectionsMutualInformation`,
the surface-meshing group.

`AlignSectionsMutualInformation` was dropped from Phase 2's roster on 2026-08-18 and lands
here. It is the sole remaining blocker on retiring `6_6_stats_test_v2` once Phase 2's
opening batch clears the other four consumers.

**Deliverable is per archive, not per filter:** a re-derived archive plus a provenance
sidecar. Filters pinned to it get their V&V report in Phase 4. The point of this phase
is to stop CI from certifying wrong answers, not to move the report counter.

**Report the hit rate.** Every archive that turns out to encode a wrong value is
evidence for how much the rest of this program is worth. Phase 2's opening batch produces
the first such measurement; this phase extends it. If the rate is near zero, Phase 4 can be
scaled back; if it is high, Phase 4 should be scaled up.

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
of 44 V&V PRs, the merged ones nearly all carry an approving review from a different
engineer — `nyoungbq` (10), `imikejackson` (5), `mmarineBlueQuartz` (3), `JDuffeyBQ` (3).

Three problems remain:

1. **A PR approval is not an oracle review.** A reviewer approving a diff is reviewing
   code, tests, and documentation. Nothing in the process obliges them to attest that
   they examined the oracle and found it independent of the implementation, and nothing
   records that they did. The policy's central control is therefore unauditable.
2. **At least one merge bypassed review entirely** — PR #1640, merged with
   `CHANGES_REQUESTED` and no approval, for a filter inside the MTR closure.
3. **Reports can pass review and still be internally inconsistent.** PR #1703 reworks the
   `GroupMicroTextureRegions` report from the #1637 cycle and, by its own summary, fixes a
   code-path count that claimed "9 of 9 exercised" while its own table showed 7, a legacy
   comparison marked "not run" for the wrong reason, and an overstated 6.5.171 migration
   impact. That report was `DRAFT`, so nothing was mis-certified — but the defects are the
   kind a diff-focused review does not catch, and there is no reason to think the
   `COMPLETE` population is immune. A cheap mitigation: have `vv_status.py` cross-check
   the dashboard's code-path count against the section table and fail on mismatch.

Two near-misses show the idea is wanted but the shape is not settled:

- PR #1701 made one report read `COMPLETE — … second-engineer sign-off recorded at PR
  review`, which names no reviewer and no PR number — it asserts the review happened
  without evidencing it.
- `ComputeFeaturePhasesFilter.md` (#1672) adds a
  `| **Sign-off** | Nathan Young - 7/15/2026 |` row. That is closer, but it names the
  **author** of the PR, not a second engineer, so it records self-certification rather
  than the independent check the policy asks for.

Both are reaching for the same missing field. Settling it now — while reports are actively
being edited — is cheaper than retrofitting 35 of them later.

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

**40 V&V reports have been authored. 36 are on `develop`.** The other 4 exist only on
unmerged branches, and **2 of those 4 are MTR closure filters — both owned by the
same engineer, both stalled in a review round-trip** (`ErodeDilateBadData`, #1687, merged
2026-08-18, is no longer in this backlog):

| Filter | PR | Blocked on |
|---|---|---|
| `MultiThresholdObjects` | #1688 | `CHANGES_REQUESTED` (@mmarineBlueQuartz) |
| `WriteDREAM3D` | #1683 | `CHANGES_REQUESTED` (@mmarineBlueQuartz) |

The remaining two are non-closure: `ComputeTriangleGeomCentroids` and `DBSCANFilter`
(#1702, awaiting first review). `GroupMicroTextureRegions` (#1703, awaiting first review)
is already on `develop` but in `DRAFT`; #1703 promotes it to `READY FOR REVIEW`.

**Progress:** `develop` went 18/29 → 23/29 → 24/29 → 25/29 → 26/29 across 2026-07-28
to 2026-08-12, then → **27/29** on 2026-08-18 (#1687, `ErodeDilateBadData`, merged).
Merged in the earlier window: #1693, #1692, #1689, #1685, #1638, #1694, #1695,
#1679, #1701, #1672.

**The SBIR deliverable is now gated on one engineer's revision queue.** Every other
contributor's closure work has landed. @mmarineBlueQuartz owes revisions on #1688 and
#1683 — two PRs standing between `develop` and 29/29. Nothing else on this plan
changes that, and no amount of new V&V substitutes for it. If those two cannot be
prioritized, reassigning the review responses is the only other lever.

**The backlog is now one engineer deep.** Every approved or awaiting-review closure PR has
landed; all three that remain are waiting on an author to answer review comments, and
**two of the three belong to @mmarineBlueQuartz** (#1688, #1683). The critical
path to the SBIR deliverable is one person's revision queue. That is worth knowing
explicitly — it is a scheduling risk, not a process failure, and it is addressed by
reassignment or by prioritizing those two above other work, not by more V&V.

The consequence is unchanged in kind, only in size: **the SBIR filter work is finished as
engineering and 93% complete as evidence.** A release tag cut today would pin
`(commit hash, archive SHA512)` pairs — the mechanism `vv_policy.md` relies on — against a
`develop` missing two of the twenty-nine closure reports.

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
