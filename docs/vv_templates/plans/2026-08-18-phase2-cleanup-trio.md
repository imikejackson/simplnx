# Phase 2 Opening Batch — Cleanup Filter Oracles Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** V&V three MTR-adjacent cleanup filters — `ErodeDilateMask`, `ErodeDilateCoordinationNumber`, `RequireMinimumSizeFeatures` — against Class 1 analytical oracles **plus a full legacy 6.5.171 A/B per filter**, and convert `ErodeDilateBadData`'s remaining archive-pinned test, fully retiring the `6_6_erode_dilate_test`, `6_6_min_size_input`, and `6_6_min_size_output` archives.

**Architecture:** Each of the three unverified filters currently loads a DREAM3D 6.6-produced `.dream3d` and diffs its output against another 6.6-produced `.dream3d`. That certifies "NX reproduces 6.6," not "NX is correct." Each task replaces it with a small hand-built voxel grid whose expected output is derived from the algorithm's stated contract — and then, **after** the oracle is green, runs the same toy fixtures through the DREAM3D 6.5.171 `PipelineRunner` binary and `nxrunner` and diffs the outputs (migration evidence, not correctness evidence). All three are face-neighbour operations on `FeatureIds`/mask arrays, so a 5×5×1 or 6×6×6 grid is enough to exercise every branch. `ErodeDilateBadData` was fully V&V'ed and merged today (PR #1687); its task is only the archive-retirement conversion its own report lists as outstanding.

**Tech Stack:** C++20, Catch2, CMake, `simplnx` `UnitTestCommon` helpers, DREAM3D 6.5.171 `PipelineRunner`, `nxrunner`, `h5py` (conda env `dream3d`).

---

## Revision history

- **2026-08-18 (rev 2):** (1) Added `ErodeDilateBadData` archive-retirement task — its full V&V merged as PR #1687 on 2026-08-18, so only the `6_6_erode_dilate_test`-pinned `(Erode)` test remains. (2) All paths retargeted to **Workspace9** (Workspace2 does not exist on this machine). (3) **Legacy A/B upgraded from Tier B source-read to full binary A/B** per requester decision — rationale: the process is largely LLM-automated so the A/B cost argument no longer dominates, and binary-verified A/B both strengthens the migration case for 6.5.171 holdouts and catches what source reading cannot (see the #1687 precedent below). (4) Added Task 0: the working branch is 17 commits behind `upstream/develop` and must sync first.

---

## Why these four, in this order

From `docs/vv_templates/vv_program_plan.md` §3.3, the first three are dual-purpose filters — each carries MTR exposure *and* is pinned to a `6_6_` archive, so each advances the SBIR coverage story and drains the defect backlog simultaneously. The fourth is a small conversion that unlocks an archive the first two cannot retire alone.

| Filter | Plugin | Scope | Archive it pins | Archive fate |
|---|---|---|---|---|
| `ErodeDilateMask` | SimplnxCore | Full V&V (101 LOC) | `6_6_erode_dilate_test` | shared — see below |
| `ErodeDilateCoordinationNumber` | SimplnxCore | Full V&V (181 LOC) | `6_6_erode_dilate_test` | shared — see below |
| `RequireMinimumSizeFeatures` | SimplnxCore | Full V&V (332 LOC) | `6_6_min_size_input`, `6_6_min_size_output` | **both fully retire** |
| `ErodeDilateBadData` | SimplnxCore | Archive conversion only — V&V already merged (#1687) | `6_6_erode_dilate_test` | **fully retires** once all three consumers are converted |

`6_6_erode_dilate_test` has **three** consumers (`ErodeDilateMaskTest`, `ErodeDilateCoordinationNumberTest`, `ErodeDilateBadDataTest`), not two. It retires only when Task 4 lands alongside Tasks 1 and 2.

`RequireMinNumNeighbors` (merged as #1694) is a near-sibling of `RequireMinimumSizeFeatures` — same feature-removal-and-reassignment shape. Read its post-V&V test at `src/Plugins/SimplnxCore/test/RequireMinNumNeighborsTest.cpp` **after the Task 0 sync** (the pre-sync copy on this branch predates #1694); its `DiscriminatingFixture` namespace is the house style for this kind of fixture and Task 3 should follow it closely.

---

## Global Constraints

- **"Legacy produced this output" is never a valid oracle for correctness.** Neither 6.5.171 nor 6.6 output is. The oracle comes first; the legacy A/B is a *diff-explanation and migration-evidence* exercise run only after SIMPLNX matches the oracle.
- **Pick the oracle before running any comparison.** Derive expected values from the algorithm source and write them into the test *before* executing the filter.
- **Every commit uses the `VV:` prefix — always.** Follow `docs/vv_templates/commit_template.md`. A V&V is *assumed* to have possibly found and fixed bugs, so a fix riding along does **not** change the prefix. Never invent `VV/BUG:` or `VV/PERF:` variants — `git log --grep='^VV:'` must return every V&V commit.
- Commit body: summarise any defect found and whether 6.5.171 shared it.
- `Signed-off-by: <name from git config>`. **Never** mention AI assistance or add an AI co-author line. Never commit `CLAUDE.md` or anything under `.claude/`.
- **A/B artifacts are never committed to git.** All legacy pipelines, toy legacy `.dream3d` inputs, run outputs, diff scripts, and notes live in `/Users/mjackson/Workspace9/ww_work/<FilterName>/` and are uploaded to OneDrive per policy at the end.
- C++20, Allman braces, 2-space indent, 200-column limit.
- Private members `m_CamelCase`; constants `k_CamelCase`; methods `camelBack`.
- Geometry variables take a `Geom` suffix; DataStore references take a `Ref` suffix.
- Every `TEST_CASE` ends with `UnitTest::CheckArraysInheritTupleDims(dataStructure);`.
- All `getDataRefAs<T>()` calls in tests wrapped in `REQUIRE_NOTHROW()`.
- Run tests with `ctest -R "<name>" --verbose`, never the test binary directly.
- Preserve every existing SIMPL `BackwardsCompatibility` `TEST_CASE` untouched.
- No filter parameters are added, removed, or re-keyed by this plan, so no parameter-version bumps are expected. If a task turns out to need one, bump the version and add a terse comment in `parameters()` stating what changed.
- clang-format the files you touched before committing: `/opt/local/clang+llvm-15.0.4-arm64-apple-darwin21.0/bin/clang-format -i -style=file <files>` run from inside the repo. Never run the full-tree sweep script for a focused change.

**Repo:** `/Users/mjackson/Workspace9/simplnx` (branch `vv/forward_plan`)
**Build directory:** `/Users/mjackson/Workspace9/DREAM3D-Build/simplnx-Rel`
(configure from repo root with `cmake --preset simplnx-Rel` if absent — the `CMakeUserPresets.json` preset resolves `binaryDir` to `${sourceDir}/../DREAM3D-Build/simplnx-Rel`).
**Legacy binary:** `/Users/mjackson/Applications/DREAM3D.app/Contents/Bin/PipelineRunner` (6.5.171)
**Legacy source:** `/Users/mjackson/Workspace/D3D_v6.5.171/DREAM3D/Source/Plugins/Processing/ProcessingFilters/` (ignore anything under `DREAM3D-Build/` — generated `moc_*` files)
**6.5.172 patch branch (only if a legacy bug is confirmed):** `/Users/mjackson/Desktop/ShellScripts/init_6_5_172_build.sh 9`, then work in `/Users/mjackson/Workspace9/6.5.172/DREAM3D`, configure with `cmake --preset D3D-Rel-Qt515-6_5_171`.

---

## THE ONE RULE THAT MATTERS

> **Derive the expected values by hand from the algorithm source. Then run the test.
> If they differ, you have found either a filter bug or an oracle error — adjudicate
> and write down which. Never run the filter first and paste its output as "expected."**

Pasting observed output is how the 6.6 exemplars were made. A test built that way is green by construction and proves nothing. The legacy A/B in each task runs **after** the oracle is green, never before — its job is to explain differences to migrating users, not to define correct.

---

## Legacy A/B protocol (applies to Tasks 1–3)

Upgraded from the earlier Tier B source-read scoping — see Revision history. The `#1687` precedent justifies the cost: during the `ErodeDilateBadData` V&V, a plausible-looking tie-break bug hypothesis survived source-level comparison (legacy's source has the identical line) and was only falsified by running both binaries on the same input and diffing an array that was not blind to the tie-break. Source reading alone provably misses things in exactly this family of filters.

Per filter, after the oracle tests are green:

1. **Write the toy input as a legacy-format `.dream3d`** using the `compare-legacy-dream3d` skill's `legacy_dream3d/writer.py` (h5py, conda env `dream3d`). The same 5×5×1 / 6×6×6 grids the unit tests build — do not invent different data for the A/B.
2. **Author the legacy pipeline JSON** per the `write-simpl-pipeline` skill: `DataContainerReader` → legacy filter → `DataContainerWriter`, one pipeline per parameter combination the unit tests exercise (for the erode/dilate filters, mirror #1687's sweep: 7 direction combinations × erode/dilate; for `RequireMinimumSizeFeatures`, the threshold/at-threshold/single-phase combinations).
3. **Author the matching NX `.d3dpipeline`** per the `write-d3dpipeline` skill and run it with `nxrunner` from the build directory.
4. **Run both**: `PipelineRunner -p <legacy.json>` and `nxrunner <nx.d3dpipeline>`; diff the outputs array-by-array with an h5py script kept in `ww_work`.
5. **Every array-level difference becomes a structured Deviation entry** per `docs/vv_templates/deviation_template.md` (Deviation ID, Filter UUID, Symptom, Root cause, Affected users, Recommendation). Because the A/B **was run**, the deviation file is now **required** for these filters (tier rules: "Only if A/B was run" — it was).
6. **If the A/B confirms a bug in legacy** (not just a difference): build the 6.5.172 branch (`init_6_5_172_build.sh 9`) and create a **surgical patch** there — smallest possible diff, one filter, no drive-by cleanup — so 6.5-line users have a corrected reference. Record the patch location in the deviation entry.
7. Everything from steps 1–6 lives in `/Users/mjackson/Workspace9/ww_work/<FilterName>/` with a `ReadMe.md` describing how to re-run it. Never committed; uploaded to OneDrive at the end.

---

## Bug adjudication protocol

Every mismatch (oracle vs. SIMPLNX, or later legacy vs. SIMPLNX) gets exactly one of two verdicts.

**Verdict A — oracle error.** Your derivation was wrong. Correct the expectation and note what you misread.

**Verdict B — filter bug.** SIMPLNX is wrong. Then:

1. **Fix the algorithm in this same PR.** Do not defer, do not open a separate issue. The test that caught it and the fix that resolves it belong in one reviewable change.
2. **Determine whether 6.5.171 had the same bug** — first by reading the legacy source, then confirmed by the task's binary A/B run.

   | SIMPLNX filter | DREAM3D 6.5.171 equivalent |
   |---|---|
   | `ErodeDilateMask` | `ErodeDilateMask.cpp` |
   | `ErodeDilateCoordinationNumber` | `ErodeDilateCoordinationNumber.cpp` |
   | `RequireMinimumSizeFeatures` | `MinSize.cpp` |

3. **Write the outcome down**, per `docs/vv_templates/deviation_template.md`:
   - **6.5.171 had it too** → a **legacy bug now corrected**. Full deviation entry: `Root cause: Bug in 6.5.171`, cite legacy file and line range, describe affected users, `Recommendation: Trust SIMPLNX`. Copy the template's `SegmentFeatures-D2` example shape. This is a user-visible behavioural change between versions and must be documented as one. This is also the trigger for the **6.5.172 surgical patch** (Legacy A/B protocol step 6).
   - **6.5.171 was correct** → a **port regression**. Fix it and describe it in the report's Summary. It produces **no deviation entry of its own** (once fixed there is no remaining difference) — but the A/B run that confirms post-fix parity still happens and is still recorded in the report's Legacy comparison section.
4. Re-run and confirm green with the fix in place.

Keep a scratch note throughout — index, expected, actual, verdict, and for Verdict B whether legacy shared the defect. It is the raw material for the deviations file and the Task 5 hit rate.

---

## Confirmed findings to start from

**1. `ErodeDilateMask` ignores its `X Direction`, `Y Direction`, and `Z Direction` parameters entirely.**

- The filter declares all three (`ErodeDilateMaskFilter.cpp:59-61`) and copies them into `InputValues` (`:107-109`).
- `ErodeDilateMaskInputValues` declares them (`Algorithms/ErodeDilateMask.hpp:29-31`).
- **`Algorithms/ErodeDilateMask.cpp` never reads them.** The neighbour loop visits all six face neighbours unconditionally.
- DREAM3D 6.5.171 *does* honour them — `ErodeDilateMask.cpp:227-247` gates each neighbour on `m_XDirOn` / `m_YDirOn` / `m_ZDirOn`.

That makes it a **port regression** — Verdict B, second branch. Three user-facing parameters silently do nothing. Task 1 must reproduce it with a test, fix it, and describe it in the report Summary with **no** deviation entry.

**2. The identical defect existed in `ErodeDilateBadData` and was already fixed and merged (PR #1687, 2026-08-18, deviation `ErodeDilateBadDataFilter-D1`).** The fix pattern to mirror: a file-local `adjustValidNeighbors(isValidFaceNeighbor, xDir, yDir, zDir)` helper that ANDs each of the six `isValidFaceNeighbor` entries against the correct axis flag using the named `VoxelNeighbors<Image3D>` constants, called immediately after `computeValidFaceNeighbors` — see `Algorithms/ErodeDilateBadData.cpp:64-80,162-163` on post-sync `develop`. Read that deviation entry's *branch-history note* before writing the Task 1 fix: an earlier attempt at this exact fix bitwise-ANDed face indices against booleans and swapped axes. Do not repeat it.

Do not take finding 1 on faith — verify each claim before writing the fix. It is stated here so the fixture is designed to catch it rather than accidentally masking it. Also check `ErodeDilateCoordinationNumber` for the same pattern while in Task 2 (it declares no direction parameters, so it is expected clean — confirm rather than assume).

---

## File Structure

| File | Responsibility |
|---|---|
| `src/Plugins/SimplnxCore/test/ErodeDilateMaskTest.cpp` | Replace exemplar test with Class 1 oracles incl. direction-flag coverage |
| `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/ErodeDilateMask.cpp` | Honour the direction flags (mirror the #1687 `adjustValidNeighbors` pattern) |
| `src/Plugins/SimplnxCore/test/ErodeDilateCoordinationNumberTest.cpp` | Replace exemplar test with Class 1 oracles |
| `src/Plugins/SimplnxCore/test/RequireMinimumSizeFeaturesTest.cpp` | Replace exemplar test with Class 1 + Class 4 oracles, incl. single-phase path |
| `src/Plugins/SimplnxCore/test/ErodeDilateBadDataTest.cpp` | Drop the archive-pinned `(Erode)` exemplar test (inline Expanded tests already cover the paths) |
| `src/Plugins/SimplnxCore/test/CMakeLists.txt` | Drop `download_test_data()` for the three retired archives |
| `src/Plugins/SimplnxCore/vv/<Filter>.md` | V&V report per filter (`python scripts/vv_init.py <FilterName>`); for `ErodeDilateBadDataFilter.md` an *update*, not a new report |
| `src/Plugins/SimplnxCore/vv/deviations/<Filter>.md` | Deviation entries (now required — A/B is run) |
| `src/Plugins/SimplnxCore/vv/provenance/6_6_erode_dilate_test.md` | Mark superseded/retired once the archive is dropped |
| `docs/vv_templates/vv_program_plan.md` | Record the hit rate; refresh the stale #1687 rows |
| `/Users/mjackson/Workspace9/ww_work/<FilterName>/` | A/B pipelines, toy legacy inputs, diffs, ReadMe (never committed) |

---

## Task 0: Sync the branch and stand up the toolchain

**Files:** none committed (branch maintenance + build setup).

The working branch is **17 commits behind `upstream/develop`**, missing at least: #1687 (`ErodeDilateBadData` V&V — the `adjustValidNeighbors` fix pattern and the current state of its test file), #1694 (`RequireMinNumNeighbors` — the `DiscriminatingFixture` house style Task 3 copies), and #1701/#1672 (report status conventions). Executing against the stale tree would rewrite files #1687 already changed.

- [ ] **Step 1: Clear the dirty working tree.** `git status` currently shows unstaged deletions under `Code_Review/vv/WritePoleFigure/`. Confirm with the requester whether these deletions are intentional; commit them separately or restore them. Do not rebase over a dirty tree.
- [ ] **Step 2: Sync.** `git fetch upstream && git rebase upstream/develop` on `vv/forward_plan`. Expected conflicts: none (local commits are docs-only).
- [ ] **Step 3: Re-verify the post-sync landscape.** Confirm: (a) `Algorithms/ErodeDilateBadData.cpp` now contains `adjustValidNeighbors` wired in at `:162-163`; (b) `RequireMinNumNeighborsTest.cpp` contains the `DiscriminatingFixture` namespace; (c) `Algorithms/ErodeDilateMask.cpp` still does **not** read the direction flags (finding 1 stands); (d) `grep -rn "6_6_erode_dilate_test" src/` still shows all three test consumers.
- [ ] **Step 4: Configure and build.** From repo root: `cmake --preset simplnx-Rel`, then `cd /Users/mjackson/Workspace9/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest`. Also run `cmake --build . --target Fetch_Remote_Data_Files` once so the baseline suite can run before archives are removed.
- [ ] **Step 5: Verify the legacy toolchain.** `ls /Users/mjackson/Applications/DREAM3D.app/Contents/Bin/PipelineRunner` and run it with no args to confirm it executes. Create `/Users/mjackson/Workspace9/ww_work/`.
- [ ] **Step 6: Baseline test run.** `ctest -R "SimplnxCore::ErodeDilate|SimplnxCore::RequireMinimum" --verbose` — record the pre-change pass/fail state in the scratch note.

---

## Task 1: ErodeDilateMask

**Files:**
- Modify: `src/Plugins/SimplnxCore/test/ErodeDilateMaskTest.cpp`
- Modify: `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/ErodeDilateMask.cpp`
- Read: `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/ErodeDilateMaskFilter.cpp` for the real `k_*_Key` names, and `Algorithms/ErodeDilateBadData.cpp:64-80` (post-sync) for the fix pattern

**Interfaces:**
- Produces: `::MaskFixture` — a 5×5×1 bool-mask builder reused conceptually (not shared across TUs) by Task 2.

### Oracle derivation

`ErodeDilateMask::operator()` double-buffers: it copies `mask` into `maskCopy`, walks every voxel, and for voxels where `!mask[voxelIndex]` inspects the six face neighbours. Then:

- **Dilate** (`Operation == k_DilateIndex`): if any face neighbour is `true`, set `maskCopy[voxelIndex] = true`. The false voxel becomes true — the true region grows.
- **Erode** (`Operation == k_ErodeIndex`): if a face neighbour is `true`, set `maskCopy[neighpoint] = false`. The *neighbour* is cleared — the true region shrinks from its boundary.

Note the asymmetry: dilate writes the centre voxel, erode writes the neighbour. Satisfy yourself that both are correct formulations before accepting them.

Use a 5×5×1 grid with a single true voxel at the centre `(2,2,0)`:

```
. . . . .
. . . . .
. . X . .
. . . . .
. . . . .
```

- **Dilate, 1 iteration** → the centre plus its 4 in-plane face neighbours become true (a plus shape). Z neighbours do not exist in a 1-thick grid.
- **Erode, 1 iteration** → the centre's false neighbours each clear it, so the grid goes all-false.
- **Dilate with `XDirOn = false`** → only the Y neighbours are added, giving a vertical 3-cell bar. **This is the assertion that catches confirmed finding 1.**
- **Erode with `YDirOn = false`, seed = horizontal 3-cell bar `{(1,2),(2,2),(3,2)}`** → only the X-axis neighbours can clear, so after 1 iteration the bar's ends are cleared from the sides leaving `{(2,2)}`; with all directions on, the Y-neighbours also participate — same surviving set here, so instead assert the flag via the *dilate* asymmetry above **and** add this erode-with-flags case with `ZDirOn = false` on a 3×3×3 grid seeded true everywhere except the centre: with Z off, the centre's ±Z true neighbours are *not* cleared; with Z on, they are. Derive both grids by hand in a comment. (The fix touches both operation branches; both need a flag-sensitive assertion.)

- [ ] **Step 1: Write the failing tests**

Replace the exemplar `TEST_CASE`s (`(Dilate)`, `(Erode)`) with a fixture and four `TEST_CASE`s (dilate, erode, direction-gated dilate, direction-gated erode). Keep the `SIMPL Backwards Compatibility` `TEST_CASE` untouched.

```cpp
namespace
{
namespace MaskFixture
{
constexpr usize k_XDim = 5;
constexpr usize k_YDim = 5;
constexpr usize k_ZDim = 1;
constexpr usize k_CellCount = k_XDim * k_YDim * k_ZDim;

constexpr usize GetIndex(usize x, usize y)
{
  return y * k_XDim + x;
}

// 5x5x1 grid, single true voxel at (2, 2).
inline DataStructure CreateFixture()
{
  DataStructure dataStructure;
  auto* imageGeomPtr = ImageGeom::Create(dataStructure, "Image");
  imageGeomPtr->setDimensions({k_XDim, k_YDim, k_ZDim});
  imageGeomPtr->setOrigin({0.0F, 0.0F, 0.0F});
  imageGeomPtr->setSpacing({1.0F, 1.0F, 1.0F});

  auto* cellAmPtr = AttributeMatrix::Create(dataStructure, "CellData", {k_ZDim, k_YDim, k_XDim}, imageGeomPtr->getId());
  imageGeomPtr->setCellData(*cellAmPtr);

  auto* maskPtr = BoolArray::CreateWithStore<DataStore<bool>>(dataStructure, "Mask", {k_ZDim, k_YDim, k_XDim}, {1}, cellAmPtr->getId());
  auto& maskRef = maskPtr->getDataStoreRef();
  for(usize i = 0; i < k_CellCount; i++)
  {
    maskRef[i] = false;
  }
  maskRef[GetIndex(2, 2)] = true;

  return dataStructure;
}

inline DataPath GeomPath()
{
  return DataPath({"Image"});
}

inline DataPath MaskPath()
{
  return DataPath({"Image", "CellData", "Mask"});
}

inline void RunFilter(DataStructure& dataStructure, uint64 operation, int32 numIterations, bool xDirOn, bool yDirOn, bool zDirOn)
{
  ErodeDilateMaskFilter filter;
  Arguments args;
  args.insertOrAssign(ErodeDilateMaskFilter::k_Operation_Key, std::make_any<ChoicesParameter::ValueType>(operation));
  args.insertOrAssign(ErodeDilateMaskFilter::k_NumIterations_Key, std::make_any<int32>(numIterations));
  args.insertOrAssign(ErodeDilateMaskFilter::k_XDirOn_Key, std::make_any<bool>(xDirOn));
  args.insertOrAssign(ErodeDilateMaskFilter::k_YDirOn_Key, std::make_any<bool>(yDirOn));
  args.insertOrAssign(ErodeDilateMaskFilter::k_ZDirOn_Key, std::make_any<bool>(zDirOn));
  args.insertOrAssign(ErodeDilateMaskFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(GeomPath()));
  args.insertOrAssign(ErodeDilateMaskFilter::k_MaskArrayPath_Key, std::make_any<DataPath>(MaskPath()));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
}

inline void CheckMask(const DataStructure& dataStructure, const std::vector<usize>& expectedTrueIndices)
{
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<BoolArray>(MaskPath()));
  const auto& maskRef = dataStructure.getDataRefAs<BoolArray>(MaskPath()).getDataStoreRef();
  for(usize i = 0; i < k_CellCount; i++)
  {
    CAPTURE(i);
    const bool expected = std::find(expectedTrueIndices.cbegin(), expectedTrueIndices.cend(), i) != expectedTrueIndices.cend();
    REQUIRE(maskRef[i] == expected);
  }
}
} // namespace MaskFixture
} // namespace

TEST_CASE("SimplnxCore::ErodeDilateMaskFilter: Class 1 dilate, one iteration", "[SimplnxCore][ErodeDilateMaskFilter]")
{
  UnitTest::LoadPlugins();
  DataStructure dataStructure = MaskFixture::CreateFixture();

  MaskFixture::RunFilter(dataStructure, ::detail::k_DilateIndex, 1, true, true, true);

  // Centre plus its four in-plane face neighbours.
  MaskFixture::CheckMask(dataStructure, {MaskFixture::GetIndex(2, 2), MaskFixture::GetIndex(1, 2), MaskFixture::GetIndex(3, 2), MaskFixture::GetIndex(2, 1), MaskFixture::GetIndex(2, 3)});
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ErodeDilateMaskFilter: Class 1 erode, one iteration", "[SimplnxCore][ErodeDilateMaskFilter]")
{
  UnitTest::LoadPlugins();
  DataStructure dataStructure = MaskFixture::CreateFixture();

  MaskFixture::RunFilter(dataStructure, ::detail::k_ErodeIndex, 1, true, true, true);

  // The lone true voxel is cleared by its false neighbours.
  MaskFixture::CheckMask(dataStructure, {});
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ErodeDilateMaskFilter: Class 1 dilate honours direction flags", "[SimplnxCore][ErodeDilateMaskFilter]")
{
  UnitTest::LoadPlugins();
  DataStructure dataStructure = MaskFixture::CreateFixture();

  // X disabled: growth is confined to Y, producing a vertical 3-cell bar.
  MaskFixture::RunFilter(dataStructure, ::detail::k_DilateIndex, 1, false, true, true);

  MaskFixture::CheckMask(dataStructure, {MaskFixture::GetIndex(2, 2), MaskFixture::GetIndex(2, 1), MaskFixture::GetIndex(2, 3)});
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}
```

Add the fourth `TEST_CASE` (erode honours direction flags) per the 3×3×3 derivation in the Oracle derivation section — build a second small fixture function (`CreateHollowCubeFixture()`) in the same namespace, seeded true everywhere except the centre `(1,1,1)`, and assert the ±Z-neighbour survival difference between `ZDirOn = true` and `ZDirOn = false`. Derive the two expected index sets by hand in a comment before running.

Resolve `detail::k_DilateIndex` / `detail::k_ErodeIndex` against their real namespace — they are declared at `Algorithms/ErodeDilateMask.hpp:21-22` inside `namespace nx::core::detail`; include that header and use whatever qualification compiles. If qualification fights you, pass the literal choice indices `0` (dilate) and `1` (erode) as documented at `ErodeDilateMaskFilter.cpp:57`.

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /Users/mjackson/Workspace9/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest \
  && ctest -R "SimplnxCore::ErodeDilateMask" --verbose
```
Expected: the first two may pass; **both direction-flag tests must fail**, because the flags are ignored. If either unexpectedly passes, stop — your fixture is not exercising the flag and needs redesigning before you go further.

- [ ] **Step 3: Fix the algorithm**

Mirror the merged #1687 pattern: add a file-local `adjustValidNeighbors(std::array<bool, k_NumFaceNeighbors>& isValidFaceNeighbor, bool xDir, bool yDir, bool zDir)` to `Algorithms/ErodeDilateMask.cpp` that ANDs each entry against its axis flag using the named `VoxelNeighbors<Image3D>` face constants (copy the axis mapping from `Algorithms/ErodeDilateBadData.cpp:64-80` — do **not** re-derive it; a previous re-derivation swapped axes, see deviation `ErodeDilateBadDataFilter-D1`'s branch-history note). Call it immediately after `computeValidFaceNeighbors(...)` at `Algorithms/ErodeDilateMask.cpp:71`, alongside the existing bounds check, not instead of it. Hoisting the helper into `NeighborUtilities.hpp` to de-duplicate with `ErodeDilateBadData` is acceptable if the reviewer prefers it — note the choice in the report either way.

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd /Users/mjackson/Workspace9/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest \
  && ctest -R "SimplnxCore::ErodeDilateMask" --verbose
```
Expected: all four PASS.

- [ ] **Step 5: Legacy A/B**

Follow the **Legacy A/B protocol**: the 5×5×1 fixture written as a legacy `.dream3d`, 7 direction combinations × {erode, dilate} × 1 iteration through both `PipelineRunner` (legacy `ErodeDilateMask`) and `nxrunner`, mask arrays diffed element-wise. Artifacts in `/Users/mjackson/Workspace9/ww_work/ErodeDilateMask/`. Expected outcome: **post-fix parity on every combination** (the regression is NX-only). Any surviving difference → Bug adjudication protocol.

- [ ] **Step 6: Write the V&V deliverables**

```bash
cd /Users/mjackson/Workspace9/simplnx && python scripts/vv_init.py ErodeDilateMaskFilter
```

Fill in per `docs/vv_templates/report_template.md`, gated by `report_gates.md`:

- **Oracle (confirmed):** `Class 1 analytical` — 5×5×1 single-seed mask (+ 3×3×3 hollow cube), expected dilate/erode neighbourhoods derived from the face-neighbour contract.
- **Summary:** describe the direction-flag regression — three parameters accepted and ignored, honoured by 6.5.171, fixed here mirroring the #1687 pattern. **Port regression → no deviation entry for it.**
- **Legacy comparison:** full binary A/B run (14 combinations); state the result per combination and cite the `ww_work` ReadMe. Cite `ErodeDilateMask.cpp:227-247` for the source-level classification.
- **Deviations file:** created; contains any A/B differences that are not the fixed regression (expected: none — say so explicitly if so).
- **Exemplar archive:** none — oracle is inline; A/B artifacts on OneDrive.
- **Status line:** `| Status | READY FOR REVIEW |` — key **not** bold-wrapped.

- [ ] **Step 7: Commit**

```bash
git add src/Plugins/SimplnxCore/test/ErodeDilateMaskTest.cpp \
        src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/ErodeDilateMask.cpp \
        src/Plugins/SimplnxCore/vv/ErodeDilateMaskFilter.md \
        src/Plugins/SimplnxCore/vv/deviations/ErodeDilateMaskFilter.md
git commit -m "$(cat <<'EOF'
VV: Erode Dilate Mask fully V&V'ed

Replaces the 6.6-derived exemplar comparison with Class 1 analytical
oracles on 5x5x1 and 3x3x3 seed masks, and records a 14-combination
binary A/B against DREAM3D 6.5.171 PipelineRunner.

Fixes a port regression: the X/Y/Z Direction parameters were accepted
and copied into the algorithm's input values but never read, so all six
face neighbours were always visited. DREAM3D 6.5.171 honoured them
(ErodeDilateMask.cpp:227-247), so this is a regression introduced by the
port rather than a deviation from legacy behaviour. Same defect class as
ErodeDilateBadDataFilter-D1 (PR #1687); fix mirrors that pattern.
Post-fix A/B shows parity on all combinations.

Signed-off-by: YOUR_NAME <YOUR_EMAIL>
EOF
)"
```

---

## Task 2: ErodeDilateCoordinationNumber

**Files:**
- Modify: `src/Plugins/SimplnxCore/test/ErodeDilateCoordinationNumberTest.cpp`
- Read: `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/ErodeDilateCoordinationNumber.cpp`, `.../ErodeDilateCoordinationNumberFilter.cpp`

### Oracle derivation

The algorithm counts, for each voxel, how many of its six face neighbours sit across a "bad/good" boundary — the condition is `(featureName > 0 && feature == 0) || (featureName == 0 && feature > 0)`, i.e. a good voxel adjacent to bad, or a bad voxel adjacent to good. That count is the voxel's *coordination number*. Voxels whose coordination number meets the `CoordinationNumber` threshold get reassigned to the most common qualifying neighbour feature. With `Loop == true` the sweep repeats until no voxel qualifies.

Two things to verify rather than assume:

1. **`featureCount` is reset selectively.** It is allocated once outside the loop (`:80`), incremented at `:120`, and zeroed at `:156` only for the features that were touched. Confirm that every increment is matched by a reset on every path — an unmatched path would leak counts between voxels and corrupt the "most common neighbour" choice. This is the single most likely defect in this filter. (The merged `ErodeDilateBadData` V&V examined its own copy of this accumulator pattern — read `vv/ErodeDilateBadDataFilter.md`'s coverage table entry for it before starting.)
2. **Tie-breaking** when two neighbour features have equal counts. Determine what the code does, then decide whether it is deterministic. A non-deterministic tie-break is a finding. **Warning from #1687:** a tie-break hypothesis in the sibling filter looked like a bug from the source alone and was falsified only by the binary A/B — do not fix a tie-break "bug" before Step 4's A/B confirms it is real.

Use a 5×5×1 grid with a single bad voxel (`FeatureIds == 0`) at the centre, surrounded by feature 1:

```
1 1 1 1 1
1 1 1 1 1
1 1 0 1 1
1 1 1 1 1
1 1 1 1 1
```

The centre has 4 good in-plane neighbours, so coordination number 4. With `CoordinationNumber = 4` it should be reassigned to feature 1; with `CoordinationNumber = 5` it should be left alone. Derive both expected grids by hand. Add a third fixture with **two** feature regions (left half feature 1, right half feature 2, one bad voxel on the seam) so the "most common neighbour" vote and its tie-break are actually exercised — derive its expectation by hand, and if the tie-break turns out unspecified, assert the invariant (reassigned to *a* qualifying neighbour feature) and record the determinism question for the A/B.

- [ ] **Step 1: Write the failing tests**

Build a `CoordinationFixture` namespace following the shape of `MaskFixture` in Task 1 — same 5×5×1 geometry, but an `Int32Array` named `FeatureIds` initialised to 1 everywhere with index `GetIndex(2, 2)` set to 0. Provide `GeomPath()`, `FeatureIdsPath()`, a `RunFilter(dataStructure, coordinationNumber, loop)` helper reading the real `k_*_Key` names from `ErodeDilateCoordinationNumberFilter.cpp:56-66`, and a `CheckFeatureIds(dataStructure, expectedIds)` helper that `CAPTURE(i)`s and `REQUIRE`s each of the 25 values. Keep the `SIMPL Backwards Compatibility` `TEST_CASE` untouched.

Write three `TEST_CASE`s:
- `CoordinationNumber = 4`, `Loop = false` → all 25 values are 1.
- `CoordinationNumber = 5`, `Loop = false` → the grid is unchanged, centre still 0.
- Two-feature seam fixture per the derivation above.

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /Users/mjackson/Workspace9/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest \
  && ctest -R "SimplnxCore::ErodeDilateCoordinationNumber" --verbose
```
Expected: FAIL until the fixture and expectations are correct.

- [ ] **Step 3: Reconcile the oracle**

Apply the **Bug adjudication protocol**. Legacy source:
`/Users/mjackson/Workspace/D3D_v6.5.171/DREAM3D/Source/Plugins/Processing/ProcessingFilters/ErodeDilateCoordinationNumber.cpp`

Pay particular attention to the `featureCount` reset question above — compare the NX reset placement against how legacy managed the same accumulator.

Re-run until green, with any algorithm fix included.

- [ ] **Step 4: Legacy A/B**

Follow the **Legacy A/B protocol**: all three fixtures, `CoordinationNumber ∈ {4, 5}`, `Loop ∈ {false, true}`, through both binaries; diff `FeatureIds` element-wise. This is also where any tie-break question from the oracle derivation gets settled with evidence instead of source reasoning. Artifacts in `/Users/mjackson/Workspace9/ww_work/ErodeDilateCoordinationNumber/`. Every unexplained difference → Bug adjudication protocol; legacy-shared bugs additionally trigger the 6.5.172 surgical patch.

- [ ] **Step 5: Write the V&V deliverables and retire the shared archive**

```bash
python scripts/vv_init.py ErodeDilateCoordinationNumberFilter
grep -rn "6_6_erode_dilate_test" src/ || echo "NO REMAINING REFERENCES"
```

Report contents mirror Task 1 Step 6 (Class 1 oracle; full A/B section; deviations file required).

The archive has **three** consumers; it is removable only after Tasks 1, 2, **and 4** are all converted. If Task 4 has already landed when you get here, and the grep is empty, remove the `6_6_erode_dilate_test.tar.gz` `download_test_data()` block from `src/Plugins/SimplnxCore/test/CMakeLists.txt` and mark `src/Plugins/SimplnxCore/vv/provenance/6_6_erode_dilate_test.md` as retired (do not delete it — it documents where the old evidence came from). Otherwise leave the block and note in the report that retirement completes in Task 4.

- [ ] **Step 6: Commit**

```bash
git add src/Plugins/SimplnxCore/test/ErodeDilateCoordinationNumberTest.cpp \
        src/Plugins/SimplnxCore/test/CMakeLists.txt \
        src/Plugins/SimplnxCore/vv/ErodeDilateCoordinationNumberFilter.md \
        src/Plugins/SimplnxCore/vv/deviations/ErodeDilateCoordinationNumberFilter.md
git commit -m "$(cat <<'EOF'
VV: Erode Dilate Coordination Number fully V&V'ed

Replaces the 6.6-derived exemplar comparison with Class 1 analytical
oracles on 5x5x1 single-bad-voxel and two-feature seam grids, and
records a binary A/B against DREAM3D 6.5.171 PipelineRunner across the
CoordinationNumber and Loop parameter grid.

Signed-off-by: YOUR_NAME <YOUR_EMAIL>
EOF
)"
```

---

## Task 3: RequireMinimumSizeFeatures

**Files:**
- Modify: `src/Plugins/SimplnxCore/test/RequireMinimumSizeFeaturesTest.cpp`
- Read: `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/RequireMinimumSizeFeatures.cpp`, and `src/Plugins/SimplnxCore/test/RequireMinNumNeighborsTest.cpp` (post-sync) for the house fixture style

### Oracle derivation

This is the largest of the three (332 lines) and the closest sibling to `RequireMinNumNeighbors`, which was V&V'ed in #1694. Read that test first — its `DiscriminatingFixture` namespace, `constexpr` index helpers, and `GetInputFeatureId` / `GetExpectedSourceIndex` pattern are the model to follow.

The filter removes features whose cell count falls below `MinAllowedFeaturesSize`, reassigns their cells to neighbouring features, and compacts the feature-level arrays. **Confirmed:** the single-phase path exists — `k_ApplySinglePhase_Key` (`"apply_single_phase"`) and `k_SinglePhaseNumber_Key` (`"phase_number"`) at `RequireMinimumSizeFeaturesFilter.hpp:32-33`, linked to `k_FeaturePhasesPath_Key` at `.cpp:81-82` — so the fixture **must** carry a second phase. Design a 6×6×6 fixture containing:

- one large feature comfortably above the threshold,
- one feature exactly *at* the threshold (must survive — verify the comparison is `<` and not `<=`),
- one feature one cell below the threshold (must be removed),
- a feature touching the volume boundary, so boundary-neighbour handling is exercised,
- a second phase with its own below-threshold feature, so `ApplySinglePhase = true, PhaseNumber = 1` removes only phase-1 offenders and leaves the phase-2 offender intact — and the all-phases run removes both.

Assert both **Class 1** (specific cell values after reassignment) and **Class 4 invariants**:

- `FeatureIds` are contiguous starting at 1 with no gaps after compaction,
- the feature-level array tuple count equals `max(FeatureId) + 1`,
- total cell count is unchanged,
- no cell retains a removed feature's id.

- [ ] **Step 1: Write the failing tests**

Build the fixture and `TEST_CASE`s as described (all-phases run; single-phase run; at-threshold survival). Follow `RequireMinNumNeighborsTest.cpp`'s structure closely — same anonymous-namespace fixture, same `CAPTURE`-per-index assertion style. Read the real `k_*_Key` names from `RequireMinimumSizeFeaturesFilter.cpp` rather than guessing. Keep the `SIMPL Backwards Compatibility` `TEST_CASE` untouched.

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /Users/mjackson/Workspace9/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest \
  && ctest -R "SimplnxCore::RequireMinimumSizeFeatures" --verbose
```
Expected: FAIL.

- [ ] **Step 3: Reconcile the oracle**

Apply the **Bug adjudication protocol**. Legacy source:
`/Users/mjackson/Workspace/D3D_v6.5.171/DREAM3D/Source/Plugins/Processing/ProcessingFilters/MinSize.cpp`

The at-threshold boundary case (`<` vs `<=`) and the volume-boundary feature are the two most likely places for a discrepancy. Re-run until green, with any algorithm fix included.

- [ ] **Step 4: Legacy A/B**

Follow the **Legacy A/B protocol**: the 6×6×6 fixture as a legacy `.dream3d` (cell `FeatureIds` + feature-level `Phases` + `NumCells`, matching what legacy `MinSize` requires — read its `dataCheck` for the exact required arrays), runs for {all-phases, single-phase} × {at-threshold, above, below} through both binaries; diff `FeatureIds` and the compacted feature arrays. Artifacts in `/Users/mjackson/Workspace9/ww_work/RequireMinimumSizeFeatures/`. Legacy-shared bugs additionally trigger the 6.5.172 surgical patch.

- [ ] **Step 5: Write the V&V deliverables and retire both archives**

```bash
python scripts/vv_init.py RequireMinimumSizeFeaturesFilter
grep -rn "6_6_min_size_input\|6_6_min_size_output" src/ || echo "NO REMAINING REFERENCES"
```

Report contents mirror Task 1 Step 6 (Class 1 + Class 4 oracles; full A/B section; deviations file required). Remove both `download_test_data()` blocks from `src/Plugins/SimplnxCore/test/CMakeLists.txt` if the grep is empty.

- [ ] **Step 6: Commit**

```bash
git add src/Plugins/SimplnxCore/test/RequireMinimumSizeFeaturesTest.cpp \
        src/Plugins/SimplnxCore/test/CMakeLists.txt \
        src/Plugins/SimplnxCore/vv/RequireMinimumSizeFeaturesFilter.md \
        src/Plugins/SimplnxCore/vv/deviations/RequireMinimumSizeFeaturesFilter.md
git commit -m "$(cat <<'EOF'
VV: Require Minimum Size Features fully V&V'ed

Replaces the 6.6-derived exemplar comparison with Class 1 analytical and
Class 4 invariant oracles on a two-phase 6x6x6 fixture, records a binary
A/B against DREAM3D 6.5.171 PipelineRunner, and retires the
6_6_min_size_input and 6_6_min_size_output archives.

Signed-off-by: YOUR_NAME <YOUR_EMAIL>
EOF
)"
```

---

## Task 4: ErodeDilateBadData — convert the last archive-pinned test

**Files:**
- Modify: `src/Plugins/SimplnxCore/test/ErodeDilateBadDataTest.cpp`
- Modify: `src/Plugins/SimplnxCore/vv/ErodeDilateBadDataFilter.md`, `src/Plugins/SimplnxCore/vv/deviations/ErodeDilateBadDataFilter.md`
- Modify: `src/Plugins/SimplnxCore/vv/provenance/6_6_erode_dilate_test.md` (mark retired)
- Modify: `src/Plugins/SimplnxCore/test/CMakeLists.txt` (if this task lands last among 1, 2, 4)

**This is NOT a fresh V&V.** `ErodeDilateBadData` was fully V&V'ed and merged as PR #1687 on 2026-08-18 (status `READY FOR REVIEW`, Class 2 oracle: 28 direction/operation/iteration combinations of genuine 6.5.171 output compiled into the test as constants, direction-flag bug fixed as deviation `ErodeDilateBadDataFilter-D1`). What remains — and what its own report lists as outstanding — is that the `(Erode)` `TEST_CASE` still loads the `6_6_erode_dilate_test` archive as a production-scale exemplar regression test.

- [ ] **Step 1: Read the merged state.** Post-sync, read `ErodeDilateBadDataTest.cpp` (the `(Erode)` archive test vs. the `(Erode) Expanded` / `(Dilate) Expanded` inline tests), the report's "Outstanding before promotion to COMPLETE" list, and the deviations doc's follow-up items 1 and 3.

- [ ] **Step 2: Decide the replacement, with the requester.** Two defensible options:
  - **(a) Drop the `(Erode)` archive test outright** *(recommended)*. The Expanded tests already cover every code path inline with legacy-binary-verified constants across 28 combinations, which is strictly stronger evidence than the 6.6 exemplar diff. The archive's production-scale value is regression breadth, not correctness — and per VV_Notes policy, exemplar archives are to be retired in favour of inline toy data.
  - **(b) Replace with a compact, provenanced legacy archive** — package the legacy 6.5.171 input/output pairs produced during this batch's A/B runs as a new (non-`6_6_`) archive with a provenance sidecar, per the deviations doc's follow-up item 1. Keeps an executable-scale A/B in CI at the cost of a new archive upload.
  Record the decision and rationale in the report.
- [ ] **Step 3: Implement it.** For (a): delete the `(Erode)` `TEST_CASE` (lines ~353–397 post-sync) and its sentinel; keep every other `TEST_CASE` including `SIMPL Backwards Compatibility` untouched. For (b): follow the `publish-exemplars` skill for packaging/upload, then rewrite the test against the new archive.
- [ ] **Step 4: Retire the archive if last.** `grep -rn "6_6_erode_dilate_test" src/ || echo "NO REMAINING REFERENCES"` — if empty (Tasks 1 and 2 already landed), remove the `download_test_data()` block from `src/Plugins/SimplnxCore/test/CMakeLists.txt` and add a "Retired 2026-08-18, superseded by inline Class 1/Class 2 oracles" banner to `vv/provenance/6_6_erode_dilate_test.md`.
- [ ] **Step 5: Update the V&V docs.** In `vv/ErodeDilateBadDataFilter.md`, strike the corresponding item from "Outstanding before promotion to COMPLETE" and describe the conversion. Do **not** change its Status line — second-engineer sign-off and the cancel path are still outstanding, and promotion to COMPLETE is that reviewer's call, not this task's.
- [ ] **Step 6: Run and commit.**

```bash
cd /Users/mjackson/Workspace9/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest \
  && ctest -R "SimplnxCore::ErodeDilateBadData" --verbose
```
Expected: all remaining `TEST_CASE`s PASS.

```bash
git add src/Plugins/SimplnxCore/test/ErodeDilateBadDataTest.cpp \
        src/Plugins/SimplnxCore/test/CMakeLists.txt \
        src/Plugins/SimplnxCore/vv/ErodeDilateBadDataFilter.md \
        src/Plugins/SimplnxCore/vv/deviations/ErodeDilateBadDataFilter.md \
        src/Plugins/SimplnxCore/vv/provenance/6_6_erode_dilate_test.md
git commit -m "$(cat <<'EOF'
VV: Retire the archive-pinned Erode Dilate Bad Data exemplar test

The 28-combination inline oracle from PR #1687 supersedes the
6_6-derived production-scale exemplar diff. Completes the conversion of
all three 6_6_erode_dilate_test consumers so the archive is retired.

Signed-off-by: YOUR_NAME <YOUR_EMAIL>
EOF
)"
```

---

## Task 5: Report the hit rate

**Files:**
- Modify: `docs/vv_templates/vv_program_plan.md` §3.3 and §6

- [ ] **Step 1: Run the full SimplnxCore suite**

```bash
cd /Users/mjackson/Workspace9/DREAM3D-Build/simplnx-Rel && ctest -R "SimplnxCore::" --verbose 2>&1 | tail -40
```
Expected: no new failures. If removing a `download_test_data()` entry broke an unrelated test, restore it — a shared archive is not retirable until every consumer is converted.

- [ ] **Step 2: Record the hit rate**

§3.3 asks the opening batch to report its defect hit rate. Add a line under the **Opening batch** subsection giving:

- how many of the three fresh V&V filters produced a Verdict B finding,
- how many of those 6.5.171 also exhibited (legacy bugs → deviation entries + 6.5.172 surgical patches),
- how many were port regressions (→ report Summary only).

The `ErodeDilateMask` direction-flag regression is one confirmed port regression before you start, so the floor is 1 of 3. Note alongside it (without double-counting) that the same defect class was independently found and fixed in `ErodeDilateBadData` by #1687 — two instances of one copy-paste-era defect family is itself a finding worth a sentence.

**This number is the point of the phase.** It is the evidence for how much Phases 3 and 4 are worth. A truthful low number is a useful result.

- [ ] **Step 3: Refresh the stale program-plan rows**

While in `vv_program_plan.md`: §3.2/§6 still list PR #1687 as `CHANGES_REQUESTED`; it merged 2026-08-18. Update the merge-backlog table and the `develop` closure count accordingly (verify the other two, #1688 and #1683, with `gh pr view` before touching their rows).

- [ ] **Step 4: Commit**

```bash
git add docs/vv_templates/vv_program_plan.md
git commit -m "$(cat <<'EOF'
VV: Record the Phase 2 opening-batch defect hit rate

Also refreshes the merge-backlog table for the #1687 merge.

Signed-off-by: YOUR_NAME <YOUR_EMAIL>
EOF
)"
```

---

## Task 6: Archive the A/B evidence

- [ ] **Step 1:** Ensure each `/Users/mjackson/Workspace9/ww_work/<FilterName>/` folder has a `ReadMe.md` (inputs, pipelines, binary versions, how to re-run, results summary) per the `archive-filter-verification` skill's structure.
- [ ] **Step 2:** Hand the folders to the requester for manual OneDrive upload per policy. **Nothing from `ww_work` is ever committed to git.**

---

## Definition of done

- [ ] Three filter tests assert against inline Class 1 (and for Task 3, Class 4) oracles with zero `.tar.gz` dependencies; the `ErodeDilateBadData` archive test is converted per the Task 4 decision.
- [ ] `grep -rn "6_6_erode_dilate_test\|6_6_min_size_input\|6_6_min_size_output" src/` returns nothing (provenance sidecars under `vv/provenance/` excepted — they document history).
- [ ] The `ErodeDilateMask` direction flags are honoured, with tests (dilate *and* erode) that fail without the fix.
- [ ] Every fresh-V&V filter has a full binary A/B against 6.5.171 recorded in its report, with artifacts + ReadMe in `ww_work` (never committed).
- [ ] Each fresh-V&V filter has `vv/<Filter>.md` and `vv/deviations/<Filter>.md`, `| Status |` un-bolded, status `READY FOR REVIEW`; the `ErodeDilateBadData` docs are updated, not recreated, and its Status line is untouched.
- [ ] Every Verdict B finding is fixed in its own task's PR, not deferred; legacy-shared defects have full deviation entries citing legacy file and lines **and** a 6.5.172 surgical patch; port regressions are in the report Summary with no deviation entry.
- [ ] `ctest -R "SimplnxCore::" --verbose` shows no new failures.
- [ ] `vv_program_plan.md` §3.3 records the hit rate and §6 reflects the #1687 merge.
- [ ] Every SIMPL `BackwardsCompatibility` `TEST_CASE` is preserved untouched.
- [ ] Every commit uses the `VV:` prefix.

## Open questions for the requester

1. ~~Does `RequireMinimumSizeFeatures` have an `ApplyToSinglePhase` parameter?~~ **Resolved 2026-08-18: yes** (`k_ApplySinglePhase_Key` / `k_SinglePhaseNumber_Key`); Task 3 carries the two-phase fixture.
2. **Are these four PRs, or one?** The plan writes separate commits and assumes separate PRs so each stays reviewable (Task 4 is small enough to ride with Task 1 or 2 if preferred). If you want them batched, say so.
3. `6_6_stats_test_v2` is **not** retired by this batch — `AlignSectionsMutualInformation` still pins it and is out of scope per the 2026-08-18 roster change. Confirm that is intended before someone tries to delete the archive.
4. **Task 4 Step 2's replacement choice** — drop the production-scale archive test (recommended) vs. package a fresh provenanced legacy archive. Needs your call before Task 4 executes.
5. **The dirty working tree** — the unstaged deletions under `Code_Review/vv/WritePoleFigure/` need a decision (commit or restore) before the Task 0 rebase.
