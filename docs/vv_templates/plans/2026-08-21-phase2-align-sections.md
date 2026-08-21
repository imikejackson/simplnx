# Phase 2 Align-Sections Batch Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Status:** READY — research dossiers complete (`.superpowers/sdd/2026-08-21-phase2-align-sections/research-{afc,am}.md`, working files, not committed); controller decisions ratified by Mike 2026-08-21.

**Goal:** Full oracle-first V&V of `AlignSectionsFeatureCentroidFilter` (SimplnxCore) and `AlignSectionsMisorientationFilter` (OrientationAnalysis), including their shared `AlignSections` base algorithm, with binary A/B against DREAM3D 6.5.171 and 6.5.172 alignment patches for any confirmed shared bugs.

**Architecture:** Same batch machinery as the feature-stats batch (`2026-08-20-phase2-feature-stats.md`): one branch + one PR per filter off `upstream/develop`, subagent-driven (implementer → adversarial review gate → fix waves), evidence folders under `/Users/mjackson/Workspace9/ww_work/<Filter>/`, rclone archival at close.

**Tech Stack:** C++20/Catch2/CMake (NX), h5py via conda env `dream3d` (legacy-format fixtures), PipelineRunner 6.5.171 at `~/Applications/DREAM3D.app`, nxrunner from `NX-Com-Qt69-Vtk96-Rel`.

## Global Constraints

- Base commit: `upstream/develop` = `9bcf794e9`. Branches: `vv/AlignSectionsFeatureCentroid`, `vv/AlignSectionsMisorientation`.
- **Out of scope:** `AlignSectionsMutualInformationFilter` (reserved for human engineers), `AlignSectionsListFilter`, `ReadH5OinaData` (later batch). Do not modify them beyond what a shared-base change mechanically requires — and if a shared-base change would alter MutualInformation behavior, STOP and escalate to the controller.
- **Shared base caution:** `AlignSections` base algorithm is consumed by FeatureCentroid, Misorientation, MutualInformation, and List. Any fix in the base must be verified against ALL consumers' existing tests, and the batch's two PRs must not both edit the base (second merge would conflict) — if a base fix is needed, it lands in exactly one PR and the other branch stacks on it or waits.
- Oracle-first discipline: derive expected values by hand from algorithm source BEFORE running the filter; never paste observed output as "expected". Oracle classes per `docs/vv_templates/oracle_classes_quick_reference.md`.
- Mutation verification: apply mutation → prove ONLY the killing test fails (record assertion index) → revert → empty diff → blind-suite proof.
- Evidence classes labeled (executed vs source-derived); absolute claims need absolute proof; data-dependent behavior stated as such.
- A/B protocol: toy fixtures written as legacy-format .dream3d (h5py, conda env `dream3d`); legacy SIMPL pipeline JSON + NX `.d3dpipeline`; run 6.5.171 PipelineRunner vs nxrunner; element-wise diff; ALL divergences predicted from source BEFORE the runs. Deviations documented vs 6.5.171. For every confirmed shared bug, a surgical patch in `/Users/mjackson/Workspace9/6.5.172/DREAM3D` proving patched-legacy ≡ fixed-NX ≡ oracle. The shift .txt file both filters can write is an additional A/B comparable.
- EbsdLib: Misorientation branch builds against local `/Users/mjackson/Workspace9/EbsdLib` @ `topic/3_1_1_staging` via preset build dir `NX-Com-Qt69-Vtk96-Rel-EbsdLib` **only if** the staging fixes are load-bearing for this filter (dossier to determine); otherwise the standard `NX-Com-Qt69-Vtk96-Rel` build is used and no vcpkg bump occurs. Never push EbsdLib without user confirmation.
- House rules: `VV:` commit prefix; `Signed-off-by: Michael Jackson <mike.jackson@bluequartz.net>` last; NEVER any AI attribution anywhere; never commit CLAUDE.md/.claude/ww_work; every TEST_CASE ends with `UnitTest::CheckArraysInheritTupleDims(dataStructure);`; `getDataRefAs` wrapped in `REQUIRE_NOTHROW`; CAPTURE-per-index; ctest only; BackwardsCompatibility TEST_CASEs untouched; clang-format `/opt/local/clang+llvm-15.0.4-arm64-apple-darwin21.0/bin/clang-format -i -style=file`; error codes unique per-file series; report gate `| Status | READY FOR REVIEW |`.
- V&V report format-of-record: `report_template.md` style (At a glance / Summary / Algorithm Relationship / Oracle / Code path coverage / Test inventory / Exemplar archive / Deviations), co-located at `src/Plugins/<Plugin>/vv/<FilterName>.md`; provenance sidecars under `src/Plugins/<Plugin>/vv/provenance/`.
- OOC runs waived (standing decision 2026-08-19). PR reviewer = second-engineer sign-off. PR titles: `VV: <Human Name> Fully V&V'ed`.
- Implementers RE-DERIVE every number in their brief independently; trust their derivation and flag discrepancies (brief numbers were wrong in 4 of 5 feature-stats tasks).
- Exclude `PY::` chained example-pipeline tests from filtered ctest expectations (known FS Task 0 finding).

## Ratified controller decisions (Mike, 2026-08-21)

1. **AFC-3 (ReferenceSlice semantics): FIX.** "Slice 0 is the slice at the Z origin (most negative Z)." `ReferenceSlice = k` must mean the physical slice at z-index k. Today (NX and legacy alike) it indexes the inverted top-down iteration array, so `ReferenceSlice = 0` means the TOP slice. Full fix requires two parts: (a) map the user's k to the internal inverted index; (b) extend the shared base `AlignSections.cpp` transfer loop to start at `i = 0` so the top slice's shift (nonzero only in FeatureCentroid reference mode after this fix) is applied — **with a mandatory proof that this is a no-op for AlignSectionsList / Misorientation / MutualInformation / FC-consecutive** (their `xShifts[0]/yShifts[0]` are always 0; add a cheap skip when both are 0). Shared bug → Deviation vs 6.5.171 + 6.5.172 surgical alignment patch. Equivalence mapping for legacy-parity runs: old `ReferenceSlice=0` ≡ new `ReferenceSlice = zDim-1`.
2. **AFC-4 (all-masked slice): zero shift + Warning.** An empty (zero in-mask voxels) non-reference slice contributes relative shift 0 and emits a Warning naming the slice; its centroid carries forward from the last valid slice so consecutive-mode chaining stays NaN-free. An empty **reference** slice is an execute error (alignment target is undefined). No more `static_cast<int64>(NaN)` UB on any path.
3. **Guard scope: per-filter preflight only.** No guard code in the shared base. The identical gaps in AlignSectionsList / AlignSectionsMutualInformation are logged as a follow-up for the human engineers, not fixed here. (The base transfer-loop bound change in decision 1 is the sole base edit in this batch, and it lands only in the FeatureCentroid PR.)
4. **AFC-1 fix shape (controller call): fill value.** Pass fillValue `"0"` to the four `CreateArrayAction`s so tuple 0 is deterministic zeros — consistent with Misorientation's existing tuple-0 = {0,0} convention (AM-10). No anchor-row write.
5. **AM-9 findShifts store/no-store duplication: NOT refactored in the V&V PR** (behavior-preserving refactors ride separately, per the adjustValidNeighbors precedent). Logged as follow-up. Guards added in this batch must be placed OUTSIDE the duplicated branches wherever possible.
6. **Archives:** remove the orphaned `6_6_align_sections_feature_centroids.tar.gz` download (FC PR) and the vestigial `align_sections.tar.gz` download+sentinel (AM PR). Keep both `6_6_` exemplars (legacy-6.6 A/B-grade regression value). Replace the two circular `output_*.dream3d` shift-array comparisons with Class 1 hand-derived assertions; the archives themselves are not re-cut.
7. **A/B binary:** 6.5.171 official `~/Applications/DREAM3D.app` (standing decision; no EbsdLib interaction on the legacy side for either filter).
8. **Docs fixes in-scope:** AFC-5 truncation-vs-round wording, AFC-12/AM-6 phantom "Linear Background Subtraction" paragraphs, AM-10 tuple-0 documentation.

## Pre-identified findings

Implementers MUST re-derive all of this from source; the dossiers were wrong in 0 known places so far, but the feature-stats batch corrected its briefs in 4 of 5 tasks. Line numbers cite `upstream/develop = 9bcf794e9`.

### AlignSectionsFeatureCentroid (dossier: research-afc.md)

| ID | Finding | Suspected verdict | Action |
|---|---|---|---|
| AFC-1 | Tuple 0 of all 4 shift arrays uninitialized (`CreateArrayAction` gets no fill value; findShifts loop starts at 1) | NX-only bug | FIX per decision 4 |
| AFC-2 | No upper-bound preflight on ReferenceSlice → OOB vector read; legacy had `-5556` but with `>` off-by-one | NX-only guard gap + legacy off-by-one | FIX (guard -68071) |
| AFC-3 | ReferenceSlice indexes inverted iteration order | Shared semantics bug | FIX per decision 1 + Deviation + 6.5.172 patch |
| AFC-4 | All-masked slice → NaN centroid → `static_cast<int64>(NaN)` UB, platform-divergent | Shared bug | FIX per decision 2 + Deviation |
| AFC-5 | Shifts truncate toward zero; NX doc claims rounding; consecutive mode sums truncated relatives | Shared design wart + NX doc bug | Document + doc fix |
| AFC-6 | `relativexshift/relativeyshift` typed `size_t`, negative via wraparound | NX-only latent (benign today) | FIX type (int64) |
| AFC-7 | Base throws on non-IDataArray cell-AM children (StringArray/NeighborList); legacy tolerated | NX-only crash | Per-filter preflight guard (-68073); base untouched; follow-up for List/MI |
| AFC-8 | No 3D-dims guard (legacy `-3010`) | NX-only guard gap | FIX (guard -68072) |
| AFC-9 | Warnings demoted to Info; NaN messages missing `{}` placeholder (4 sites) | NX-only bugs | FIX |
| AFC-10 | Range/NaN warning flags cross-suppress | Shared benign | Document |
| AFC-11 | ReferenceSlice validated (`<0`) even when unused; not validated against zDim when used | Benign + AFC-2 | Make -68064/-68071 conditional on UseReferenceSlice |
| AFC-12 | Dead codes -68001..-68004; phantom doc paragraph | NX cruft | FIX |

### AlignSectionsMisorientation (dossier: research-am.md)

| ID | Finding | Suspected verdict | Action |
|---|---|---|---|
| AM-1 | Preflight never ties Quats/Phases/Mask tuple count to selected geometry cell count → OOB | NX-only guard gap | FIX (guard -68006) |
| AM-2 | No 3D-dims guard; dims of 1 → silent no-op or negative-index UB on memo vector | NX-only guard gap | FIX (guard -68005) |
| AM-3 | `misorients[idx]` read before bounds check; OOB near half-volume shifts | Shared latent | Document; SIZE FIXTURES to avoid (≥32×32) |
| AM-4 | `count==0` → NaN disorientation — unreachable per bounds argument | Shared benign | Verify unreachability claim, document |
| AM-5 | `crystalStructures[cellPhases[pos]]` OOB when phase ≥ ensemble tuple count | Shared UB on hostile data | FIX execute-time max-phase guard (-68008) |
| AM-6 | Doc describes Linear Background Subtraction this filter doesn't have | NX doc bug | FIX |
| AM-7 | Error-code collision: two constants = -68063; unused constants | NX cruft | FIX |
| AM-8 | Tolerance-constant ULP drift + float→double angle path can flip only within ~1e-6 rad of threshold | Shared benign, quantified | Document as A/B deviation class; fixtures ≥0.5° from tolerance |
| AM-9 | findShifts duplicated ~100 lines store/no-store; redundant double cancel check | NX maintainability | Follow-up (decision 5); dedup cancel check OK |
| AM-10 | Shift arrays tuple 0 never written (zero by default-init luck) | NX benign quirk | Make deterministic explicit (fillValue "0") + document |
| AM-11 | Asymmetric OR tie-break (shared with legacy) order-sensitive on ties | Shared benign | Document; oracle fixtures need unique strict minimum |

Cross-checks required of implementers: AFC/AM dossiers disagree on nothing structural, but AM found EbsdLib CubicOps atan2 fix already ships in the pinned v3.1.0 (`vcpkg.json`) — so **the Misorientation branch builds in the standard `NX-Com-Qt69-Vtk96-Rel` dir, no EbsdLib preset, no vcpkg bump** (verify vcpkg pin before relying on this).

---

### Task 0: Batch sync + baseline

- [ ] Verify all prior-batch PRs' branches untouched; fetch upstream; confirm base `9bcf794e9`.
- [ ] Rebuild `NX-Com-Qt69-Vtk96-Rel` on a scratch checkout of `upstream/develop`; run `ctest -R "SimplnxCore::AlignSections|OrientationAnalysis::AlignSections"` baseline; record pass/fail.
- [ ] Verify legacy PipelineRunner runs (`~/Applications/DREAM3D.app/Contents/Bin/PipelineRunner`), conda env `dream3d` importable h5py.
- [ ] Confirm archive consumers for `align_sections.tar.gz` (shared with MutualInformation?) — retirement decisions must not break the reserved filter's test.

### Task 1: AlignSectionsFeatureCentroid V&V

**Branch:** `vv/AlignSectionsFeatureCentroid` off `9bcf794e9`. Evidence: `/Users/mjackson/Workspace9/ww_work/AlignSectionsFeatureCentroid/`.

**Files:**
- Modify: `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/AlignSectionsFeatureCentroidFilter.cpp` (guards -68071/-68072/-68073 conditional per AFC-11; fillValue "0" on the 4 CreateArrayActions; delete dead -68001..-68004)
- Modify: `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/AlignSectionsFeatureCentroid.cpp` (AFC-3 physical→internal index mapping; AFC-4 empty-slice policy; AFC-6 int64 locals; AFC-9 Warning severity + fmt placeholders)
- Modify: `src/simplnx/Utilities/AlignSections.cpp` (transfer loop `i=0` extension with zero-shift skip — SOLE base edit, no-op proof required against List/Misorientation/MutualInformation existing tests)
- Modify: `src/Plugins/SimplnxCore/docs/AlignSectionsFeatureCentroidFilter.md` (AFC-5 truncation wording, AFC-12 phantom paragraph, new semantics + guards)
- Modify: `src/Plugins/SimplnxCore/test/AlignSectionsFeatureCentroidTest.cpp` (oracle TEST_CASEs added; "output test" circular exemplar comparison replaced by hand-derived assertions; "Algorithm Test" updated to `ReferenceSlice = zDim-1` equivalence mapping so the 6_6 exemplar stays a valid legacy-parity pin — document the mapping in-test)
- Modify: `src/Plugins/SimplnxCore/test/CMakeLists.txt` (remove orphaned `6_6_align_sections_feature_centroids.tar.gz` download, line ~222)
- Create: `src/Plugins/SimplnxCore/vv/AlignSectionsFeatureCentroidFilter.md` (report, format-of-record) + provenance sidecar for the retained archive.

**Oracle fixtures (Class 1, dossier §6 F1–F10 — re-derive every number):** F1 consecutive-mode integer offsets; F2 reference-slice physical-semantics probe (post-fix: `ReferenceSlice=2` on a 3-slice volume anchors physical slice 2); F3 fractional-centroid truncation pins (±0.6, ±1.5); F4 off-edge push / zero-fill incl. the Mask array itself; F5 StoreAlignmentShifts full contract incl. deterministic tuple 0; F6 all-masked middle slice → zero-shift+Warning; F6b all-masked REFERENCE slice → execute error; F7 degenerate dims → -68072; F8 ReferenceSlice bounds → -68071 (and -68064/-68071 skipped when UseReferenceSlice=false); F9 uint8-mask parity; F10 non-unit spacing invariance. RED-first for every behavior change; mutation verification per policy.

**A/B matrix:** legacy-format toy .dream3d fixtures (bool mask only — legacy requires DataArray<bool>); modes: consecutive, reference(old-0≡new-top), reference(non-trivial post-fix slice) — the last is the AFC-3 divergence run, predicted BEFORE execution; with/without shift output. Compare all cell arrays + legacy shift .txt vs NX shift arrays (row iter ↔ tuple iter; pre-document legacy "New X/Y Shift"-always-0 bug, NX extra tuple 0, text precision). Exclude NaN/UB fixtures (F6 family) from A/B — platform-divergent. 6.5.172 surgical patches: AFC-3 (RS mapping + top-slice application) and, if graded shared-and-fixed, AFC-4; prove patched-legacy ≡ fixed-NX ≡ oracle.

**Deviations expected (grade in execution):** AFC-3 fix (D, alignment-patched), AFC-4 fix (D), legacy New-Shift-columns-always-0 (legacy-only bug), legacy -5556 `>` off-by-one (legacy-only), warning-severity/message deltas, shift output file→arrays (#1237, structural).

### Task 2: AlignSectionsMisorientation V&V

**Branch:** `vv/AlignSectionsMisorientation` off `9bcf794e9` — **independent of Task 1's base edit**; its findShifts never writes shift[0], so it is provably unaffected. Standard build dir `NX-Com-Qt69-Vtk96-Rel` (no EbsdLib preset, no vcpkg change — the CubicOps atan2 fix already ships in the pinned EbsdLib ≥ 3.1.0; verify the pin first). Evidence: `/Users/mjackson/Workspace9/ww_work/AlignSectionsMisorientation/`.

**Files:**
- Modify: `src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/AlignSectionsMisorientationFilter.cpp` (guards -68005/-68006/-68007; fix -68063 collision; delete unused constants; fillValue "0" on the 3 shift CreateArrayActions)
- Modify: `src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/Algorithms/AlignSectionsMisorientation.cpp` (execute-time max-phase guard -68008 + unknown-structure warning -68009 placed BEFORE/OUTSIDE the duplicated branches; dedup the redundant double cancel check only)
- Modify: `src/Plugins/OrientationAnalysis/docs/AlignSectionsMisorientationFilter.md` (AM-6 phantom paragraph; tuple-0 semantics; guards)
- Modify: `src/Plugins/OrientationAnalysis/test/AlignSectionsMisorientationTest.cpp` (oracle TEST_CASEs; "output test" circular shift-array comparison replaced by hand-derived assertions; Small IN100 6_6 exemplar test KEPT as legacy-parity pin; drop the `align_sections.tar.gz` sentinel)
- Modify: `src/Plugins/OrientationAnalysis/test/CMakeLists.txt` (remove `align_sections.tar.gz` download, line ~136)
- Create: `src/Plugins/OrientationAnalysis/vv/AlignSectionsMisorientationFilter.md` + provenance sidecar.

**Oracle fixtures (Class 1, dossier §7 O1–O6 — re-derive, ESPECIALLY the shift sign convention of O1):** 32×32×3 ImageGeom (sized to stay out of AM-3's OOB zone), stride-4-aligned block features, cubic pair identity vs 30°@[001] (6× margin over 5° tolerance). O1 multi-hop convergence + accumulation + array wiring + zero edge-fill; O2 tolerance bracket 29°/31° (deg→rad kill) — never assert exactly at 30°; O3 mask XOR + count-denominator + all-false-mask tie-walk; O4 multi-phase/phase-0/cross-Laue mismatch semantics; O5 guard fixtures (-68005/-68006/-68008 + structure-999 warning); O6 hexagonal path (20°@[0001]) covering the generic 2·acos EbsdLib path. Unique strict minima everywhere (AM-11). RED-first; mutation verification.

**A/B matrix:** legacy-format toy .dream3d (Quats float32×4 (x,y,z,w), Phases int32, Mask bool, CrystalStructures uint32 with tuple0=999); legacy pipeline uses `MisorientationTolerance`/`UseGoodVoxels`/`WriteAlignmentShifts`/`AlignmentShiftFileName` + path objects; absolute shift-file path (legacy preflight requires it). Fixture disorientations ≥0.5° from tolerance → predicted bit-identical shifts and cell arrays; shift .txt rows ↔ NX arrays via D11 mapping. All divergences (expected: none beyond structural file-vs-array) predicted pre-run. No 6.5.172 patch expected (no confirmed shared output bug); if one emerges, standard alignment-patch protocol applies.

**Deviations expected:** shift file→arrays (#1237 structural), float→double angle path (AM-8, quantified benign), mask type widening (NX superset), guard additions (NX errors where legacy accepted or vice versa — enumerate).

### Task 3: Hit-rate recording

- [ ] Append batch hit-rate paragraph to `docs/vv_templates/vv_program_plan.md` on `vv/forward_plan`; update archive-consumer tables.
- [ ] Log follow-ups for the human engineers / ticket queue: (a) guard parity for AlignSectionsList + AlignSectionsMutualInformation (3D-dims, tuple-vs-geometry, non-IDataArray children — per decision 3); (b) AM-9 findShifts store/no-store dedup refactor (separate ENH PR, adjustValidNeighbors precedent).

### Task 4: Archival + PRs

- [ ] `rclone copy /Users/mjackson/Workspace9/ww_work/<Filter> vv_work:<Filter> --exclude "__pycache__/**"` per filter; verify with `rclone size`.
- [ ] Push branches to origin (imikejackson fork); open PRs against `upstream develop`, titles `VV: Align Sections Feature Centroid Fully V&V'ed` / `VV: Align Sections Misorientation Fully V&V'ed`; Summary/Test Plan bodies; cross-PR notes for any shared-base stacking.
