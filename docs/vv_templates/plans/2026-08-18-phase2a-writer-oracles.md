# Phase 2a — Writer Filter Class 1 Oracles Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the circular 6.6-era exemplar comparisons in four file-writer filters with Class 1 analytical oracles derived from the file format, eliminating their dependence on untrusted golden data.

**Architecture:** Each writer currently loads a `.dream3d` produced by DREAM3D 6.6, runs the filter, and diffs the output against a `.am`/`.txt` file also produced by DREAM3D 6.6. That certifies "NX reproduces 6.6," not "NX is correct." Each task replaces that with an in-memory geometry of a dozen cells and a hand-derived expected file, asserted byte-for-byte. No archive, no legacy run, no network fetch.

**Tech Stack:** C++20, Catch2, CMake, `simplnx` `UnitTestCommon` helpers.

---

## Global Constraints

Copied verbatim from `docs/vv_templates/vv_policy.md` and `CLAUDE.md`:

- **"Legacy 6.5.171 produced this output" is never a valid oracle for correctness.** Neither is 6.6 output. This plan exists because those exemplars are untrusted.
- **Pick the oracle before running any comparison.** The expected file contents must be derived by reading the algorithm source and the file format, and written into the test, *before* the filter is executed even once.
- C++20, Allman braces, 2-space indent, 200-column limit, `.hpp`/`.cpp`.
- Private members `m_CamelCase`; globals/constants `k_CamelCase`; methods/functions `camelBack`.
- Geometry variables take a `Geom` suffix (`imageGeom`); DataStore references take a `Ref` suffix (`featureIdsRef`).
- Every `TEST_CASE` ends with `UnitTest::CheckArraysInheritTupleDims(dataStructure);`.
- All `getDataRefAs<T>()` calls in tests are wrapped in `REQUIRE_NOTHROW()`.
- Run tests with `ctest -R "<name>" --verbose`, never the test binary directly.
- **Every commit in this plan that lands a V&V cycle uses the `VV:` prefix — always.** Follow `docs/vv_templates/commit_template.md`. A V&V is *assumed* to have possibly found and fixed bugs, so a bug fix riding along does **not** change the prefix. Do not invent `VV/BUG:` or `VV/PERF:` variants: `git log --grep='^VV:'` must return every V&V commit, and a slashed variant breaks that anchor.
- Commit body: summarise the defect and whether 6.5.171 shared it, per the Bug adjudication protocol.
- `Signed-off-by: <name from git config>`. **Never** mention AI assistance or add an AI co-author line. Never commit `CLAUDE.md` or anything under `.claude/`.

**Build directory:** `/Users/mjackson/Workspace2/DREAM3D-Build/simplnx-Rel`
(configure from repo root with `cmake --preset simplnx-Rel` if it does not exist).

---

## THE ONE RULE THAT MATTERS

> **Derive the expected bytes by hand from the algorithm source. Then run the test.
> If they differ, you have found either a filter bug or an oracle error — adjudicate
> and write down which. Never run the filter first and paste its output as "expected."**

Pasting observed output is how the 6.6 exemplars were created in the first place. A test built that way is green by construction and proves nothing. If you find yourself tempted, stop and re-read the algorithm.

A bug found this way is the highest-value output of this entire phase. Handle it via the protocol below.

---

## Bug adjudication protocol

Every mismatch between your derived expectation and the filter's actual output gets exactly one of two verdicts.

**Verdict A — oracle error.** Your derivation was wrong. Correct the expectation. Note what you misread; if the source was genuinely ambiguous, that is worth a sentence in the report's Oracle section.

**Verdict B — filter bug.** SIMPLNX is wrong. Then:

1. **Fix the algorithm in this same PR.** Do not defer it, do not open a separate issue. The test that caught it and the fix that resolves it belong in one reviewable change.
2. **Determine whether 6.5.171 had the same bug** by reading the legacy source. This is a source read, not a pipeline run — you are not performing the Tier A legacy A/B, which stays out of scope.

   | SIMPLNX filter | DREAM3D 6.5.171 equivalent |
   |---|---|
   | `WriteAvizoUniformCoordinate` | `AvizoUniformCoordinateWriter.cpp` |
   | `WriteAvizoRectilinearCoordinate` | `AvizoRectilinearCoordinateWriter.cpp` |
   | `WriteSPParksSites` | `SPParksSitesWriter.cpp` |
   | `WriteLosAlamosFFT` | `LosAlamosFFTWriter.cpp` |

   All four live in:
   ```
   /Users/mjackson/Workspace/D3D_v6.5.171/DREAM3D/Source/Plugins/ImportExport/ImportExportFilters/
   ```
   Ignore anything under `DREAM3D-Build/` — those are generated `moc_*` files.

3. **Write the outcome into the deviations file**, following `docs/vv_templates/deviation_template.md`. Which case you are in determines the entry:

   - **Bug present in 6.5.171 too** → this is a **legacy bug now corrected**. Write a full deviation entry: `Root cause: Bug in 6.5.171`, cite the legacy file and line range, describe who was affected, and set `Recommendation: Trust SIMPLNX`. The template's "Legacy-bug example" (`SegmentFeatures-D2`) is the exact shape to copy. **This is a user-visible behavioural change between versions and must be documented as one** — a user comparing old and new output will see the difference and needs to know why.
   - **Bug introduced in the SIMPLNX port** (6.5.171 was correct) → this is a regression, not a deviation. Fix it, and note it in the report's Summary as a port defect found by V&V. Do **not** write a deviation entry: once fixed, there is no remaining difference from 6.5.171 to document.

4. Re-run the test and confirm green with the fix in place.

Keep a scratch note as you go — line index, expected, actual, verdict, and for Verdict B whether legacy shared the defect. That note is the raw material for both the deviations file and the hit-rate line in Task 6.

---

## Scope

**In scope — four writers whose output is short, plain-text, and fully determined by their inputs:**

| Filter | Plugin | Algorithm LOC | Current data dependency |
|---|---|---|---|
| `WriteAvizoUniformCoordinate` | SimplnxCore | 117 | `6_6_avizo_writers.tar.gz` |
| `WriteAvizoRectilinearCoordinate` | SimplnxCore | 143 | `6_6_avizo_writers.tar.gz` |
| `WriteSPParksSites` | SimplnxCore | 141 | `6_5_spparks_sites_writer.tar.gz` |
| `WriteLosAlamosFFT` | SimplnxCore | 93 | `bin_feature_phases.tar.gz` → `6_6_find_feature_phases_binary.dream3d` |

**Explicitly out of scope**, and why — `docs/vv_templates/vv_program_plan.md` §3.3 lists nine "cheap" writers. That was too optimistic. These five are not cheap and each needs its own plan:

| Filter | Algorithm LOC | Why deferred |
|---|---|---|
| `WriteStlFile` | 692 | Binary STL, per-feature file splitting, multiple grouping modes |
| `WriteAbaqusHexahedron` | 460 | Emits four interdependent files |
| `WriteGBCDGMTFile` | 316 | Requires a GBCD computation upstream; oracle is not the file format |
| `WriteINLFile` | 205 | Needs crystallography/ensemble data to build a meaningful fixture |
| `WriteGBCDTriangleData` | 84 | Small, but coupled to the GBCD archive and its upstream filter |

**Update `vv_program_plan.md` §3.3 to reflect this split** as the final task.

---

## File Structure

| File | Responsibility |
|---|---|
| `test/UnitTestCommon/include/simplnx/UnitTest/UnitTestCommon.hpp` | Add `CompareTextFilesExact` — byte-exact comparison with named volatile-line predicates. Replaces blind index-based line skipping. |
| `src/Plugins/SimplnxCore/test/WriteAvizoUniformCoordinateTest.cpp` | Replace exemplar test with Class 1 oracle; keep the SIMPL back-compat test. |
| `src/Plugins/SimplnxCore/test/WriteAvizoRectilinearCoordinateTest.cpp` | Same. |
| `src/Plugins/SimplnxCore/test/WriteSPParksSitesTest.cpp` | Same. |
| `src/Plugins/SimplnxCore/test/WriteLosAlamosFFTTest.cpp` | Same. |
| `src/Plugins/SimplnxCore/test/CMakeLists.txt` | Remove `download_test_data()` entries once no test references them. |
| `src/Plugins/SimplnxCore/vv/Write*.md` | V&V report per filter (scaffold with `python scripts/vv_init.py <FilterName>`). |
| `src/Plugins/SimplnxCore/vv/deviations/Write*.md` | Deviation entries; "none observed" is a valid content if true. |
| `docs/vv_templates/vv_program_plan.md` | Correct the §3.3 writer roster. |

---

## Task 1: Byte-exact comparison helper with volatile-line predicates

The existing `UnitTest::CompareAsciiFiles(computed, exemplar, linesToSkip)` takes a vector of line *indices* to ignore. Both Avizo tests pass `{6, 7}` to skip the `Author` and `DateTime` lines. Index-based skipping is brittle — it silently stops checking the right line the moment a header line is added — and it discards the ability to assert anything at all about those lines.

This task adds a helper that compares against an expected *string* and matches designated lines by regex instead of ignoring them.

**Files:**
- Modify: `test/UnitTestCommon/include/simplnx/UnitTest/UnitTestCommon.hpp` (add alongside `CompareAsciiFiles` at line ~1550)

**Interfaces:**
- Produces:
  ```cpp
  namespace nx::core::UnitTest
  {
  // Compares `actualFileContents` line-by-line against `expectedLines`.
  // An expected line of the form "<REGEX>...</REGEX>" is matched as a std::regex
  // against the actual line instead of compared literally.
  // Fails via CAPTURE + REQUIRE on the first differing line, reporting the index.
  void CompareTextFilesExact(const std::filesystem::path& actualFilePath, const std::vector<std::string>& expectedLines);
  }
  ```

- [ ] **Step 1: Write the failing test**

Create `test/UnitTestCommon/CompareTextFilesExactTest.cpp` — if the UnitTestCommon target has no test executable, put this `TEST_CASE` temporarily at the top of `src/Plugins/SimplnxCore/test/WriteAvizoUniformCoordinateTest.cpp` and move it later.

```cpp
TEST_CASE("UnitTest::CompareTextFilesExact: literal and regex lines", "[UnitTestCommon]")
{
  const fs::path tempFile = fs::path(fmt::format("{}/compare_text_files_exact.txt", unit_test::k_BinaryTestOutputDir));
  fs::create_directories(tempFile.parent_path());
  {
    std::ofstream out(tempFile);
    out << "header line\n";
    out << "DateTime \"Mon Aug 18 09:15:00 2026\"\n";
    out << "trailing value 42 \n";
  }

  const std::vector<std::string> expectedLines = {
      "header line",
      R"(<REGEX>DateTime "[A-Z][a-z]{2} [A-Z][a-z]{2} .*"</REGEX>)",
      "trailing value 42 ",
  };

  UnitTest::CompareTextFilesExact(tempFile, expectedLines);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/mjackson/Workspace2/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest
```
Expected: FAIL to compile — `'CompareTextFilesExact' is not a member of 'nx::core::UnitTest'`.

- [ ] **Step 3: Write minimal implementation**

Add to `UnitTestCommon.hpp`:

```cpp
inline void CompareTextFilesExact(const std::filesystem::path& actualFilePath, const std::vector<std::string>& expectedLines)
{
  std::ifstream actualFile(actualFilePath);
  REQUIRE(actualFile.is_open());

  std::vector<std::string> actualLines;
  for(std::string line; std::getline(actualFile, line);)
  {
    actualLines.push_back(line);
  }

  const std::string regexOpen = "<REGEX>";
  const std::string regexClose = "</REGEX>";

  for(usize i = 0; i < expectedLines.size(); i++)
  {
    CAPTURE(i);
    REQUIRE(i < actualLines.size());

    const std::string& expected = expectedLines[i];
    const std::string& actual = actualLines[i];
    CAPTURE(expected);
    CAPTURE(actual);

    if(expected.rfind(regexOpen, 0) == 0 && expected.size() >= regexOpen.size() + regexClose.size())
    {
      const std::string pattern = expected.substr(regexOpen.size(), expected.size() - regexOpen.size() - regexClose.size());
      REQUIRE(std::regex_match(actual, std::regex(pattern)));
    }
    else
    {
      REQUIRE(actual == expected);
    }
  }

  // No extra trailing content beyond what was specified.
  REQUIRE(actualLines.size() == expectedLines.size());
}
```

Add `#include <regex>` to the header's include block if absent.

- [ ] **Step 4: Run test to verify it passes**

```bash
cd /Users/mjackson/Workspace2/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest \
  && ctest -R "UnitTestCommon" --verbose
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add test/UnitTestCommon/include/simplnx/UnitTest/UnitTestCommon.hpp
git commit -m "$(cat <<'EOF'
TEST: Add CompareTextFilesExact with regex-matched volatile lines

Index-based line skipping in CompareAsciiFiles silently stops checking
the intended line whenever a header line is added, and discards any
assertion about the skipped content.

Signed-off-by: YOUR_NAME <YOUR_EMAIL>   # from `git config user.name` / `user.email`
EOF
)"
```

---

## Task 2: WriteAvizoUniformCoordinate Class 1 oracle

**Files:**
- Modify: `src/Plugins/SimplnxCore/test/WriteAvizoUniformCoordinateTest.cpp` (replace the first `TEST_CASE`; leave the `BackwardsCompatibility` `TEST_CASE` untouched)
- Read for the oracle: `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/WriteAvizoUniformCoordinate.cpp`, `src/Plugins/SimplnxCore/src/SimplnxCore/utils/AvizoWriter.cpp`

**Interfaces:**
- Consumes: `UnitTest::CompareTextFilesExact` from Task 1.
- Produces: the `::AvizoFixture` namespace (fixture geometry builder) reused by Task 3.

### Oracle derivation — do this before writing any code

Read `WriteAvizoUniformCoordinate::generateHeader` and `writeData` and derive, on paper, what a 3×2×2 ImageGeom at origin (0,0,0) with spacing (1,1,1) and `FeatureIds = {1,2,3,4,5,6,7,8,9,10,11,12}` must produce. Points to reason through:

1. Line 1 is `# AmiraMesh 3D ASCII 2.0` when `WriteBinaryFile == false`.
2. `Content` uses `%llux%llux%llu` → `3x2x2 int, uniform coordinates`.
3. `BoundingBox` uses `%f` (six decimals) and computes `origin + spacing * dims` — **not** `origin + spacing * (dims - 1)`. Decide for yourself whether that is correct for a cell-centred grid before you accept it; if you think it is wrong, that is a deviation to record, not a test to bend.
4. `FeatureIds Path` prints `DataPath::toString()`. Confirm the separator character by reading `DataPath::toString()` — do not guess.
5. The ASCII data loop is:
   ```cpp
   int count = 0;
   for(size_t i = 0; i < totalPoints; ++i)
   {
     fprintf(outputFile, "%d", featureIds[i]);
     if(count < 20) { fprintf(outputFile, " "); count++; }
     else           { fprintf(outputFile, "\n"); count = 0; }
   }
   fprintf(outputFile, "\n");
   ```
   Work out exactly how many values land on the first line, and whether a trailing space precedes the newline. With 12 values none of them trips the `count == 20` branch, so all 12 are emitted followed by a space, then the loop's final `fprintf` adds the newline. Write down the resulting line character-for-character.
6. `Author` and `DateTime` are volatile — `DateTime` embeds `std::ctime`. These are the two `<REGEX>` lines.

- [ ] **Step 1: Write the failing test**

Replace the `Valid Filter Execution` `TEST_CASE` with:

```cpp
namespace
{
namespace AvizoFixture
{
constexpr usize k_XDim = 3;
constexpr usize k_YDim = 2;
constexpr usize k_ZDim = 2;
constexpr usize k_CellCount = k_XDim * k_YDim * k_ZDim;

const std::string k_GeomName = "Image";
const std::string k_CellDataName = "CellData";
const std::string k_FeatureIdsName = "FeatureIds";

// Builds a 3x2x2 ImageGeom, origin (0,0,0), spacing (1,1,1),
// FeatureIds = 1..12 in x-fastest order.
inline DataStructure CreateFixture()
{
  DataStructure dataStructure;
  auto* imageGeomPtr = ImageGeom::Create(dataStructure, k_GeomName);
  imageGeomPtr->setDimensions({k_XDim, k_YDim, k_ZDim});
  imageGeomPtr->setOrigin({0.0F, 0.0F, 0.0F});
  imageGeomPtr->setSpacing({1.0F, 1.0F, 1.0F});

  auto* cellAmPtr = AttributeMatrix::Create(dataStructure, k_CellDataName, {k_ZDim, k_YDim, k_XDim}, imageGeomPtr->getId());
  imageGeomPtr->setCellData(*cellAmPtr);

  auto* featureIdsPtr = Int32Array::CreateWithStore<DataStore<int32>>(dataStructure, k_FeatureIdsName, {k_ZDim, k_YDim, k_XDim}, {1}, cellAmPtr->getId());
  auto& featureIdsRef = featureIdsPtr->getDataStoreRef();
  for(usize i = 0; i < k_CellCount; i++)
  {
    featureIdsRef[i] = static_cast<int32>(i + 1);
  }

  return dataStructure;
}

inline DataPath GeomPath()
{
  return DataPath({k_GeomName});
}

inline DataPath FeatureIdsPath()
{
  return DataPath({k_GeomName, k_CellDataName, k_FeatureIdsName});
}
} // namespace AvizoFixture
} // namespace

TEST_CASE("SimplnxCore::WriteAvizoUniformCoordinateFilter: Class 1 ASCII oracle", "[SimplnxCore][WriteAvizoUniformCoordinateFilter]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure = AvizoFixture::CreateFixture();
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(AvizoFixture::GeomPath()));

  const fs::path outputPath(fmt::format("{}/avizo_uniform_class1.am", unit_test::k_BinaryTestOutputDir));

  WriteAvizoUniformCoordinateFilter filter;
  Arguments args;
  args.insertOrAssign(WriteAvizoUniformCoordinateFilter::k_OutputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(outputPath));
  args.insertOrAssign(WriteAvizoUniformCoordinateFilter::k_WriteBinaryFile_Key, std::make_any<bool>(false));
  args.insertOrAssign(WriteAvizoUniformCoordinateFilter::k_GeometryPath_Key, std::make_any<DataPath>(AvizoFixture::GeomPath()));
  args.insertOrAssign(WriteAvizoUniformCoordinateFilter::k_FeatureIdsArrayPath_Key, std::make_any<DataPath>(AvizoFixture::FeatureIdsPath()));
  args.insertOrAssign(WriteAvizoUniformCoordinateFilter::k_Units_Key, std::make_any<StringParameter::ValueType>("microns"));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)

  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)

  // Derived by hand from WriteAvizoUniformCoordinate::generateHeader/writeData.
  // This listing is one derivation; check every line against the source yourself
  // before running, and treat any disagreement as a finding to adjudicate.
  const std::vector<std::string> expectedLines = {
      "# AmiraMesh 3D ASCII 2.0",
      "",
      "# Dimensions in x-, y-, and z-direction",
      "define Lattice 3 2 2",
      "Parameters {",
      "     DREAM3DParams {",
      R"(<REGEX>         Author "DREAM3D-NX SimplnxCore Version .*",</REGEX>)",
      R"(<REGEX>         DateTime ".*"</REGEX>)",
      R"(         FeatureIds Path "Image/CellData/FeatureIds")",
      "     }",
      "     Units {",
      R"(         Coordinates "microns")",
      "     }",
      R"(     Content "3x2x2 int, uniform coordinates",)",
      "     # Bounding Box is xmin xmax ymin ymax zmin zmax",
      "     BoundingBox 0.000000 3.000000 0.000000 2.000000 0.000000 2.000000",
      R"(     CoordType "uniform")",
      "}",
      "",
      "Lattice { int FeatureIds } = @1",
      "# Data section follows",
      "@1",
      "1 2 3 4 5 6 7 8 9 10 11 12 ",
  };

  UnitTest::CompareTextFilesExact(outputPath, expectedLines);
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/mjackson/Workspace2/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest \
  && ctest -R "SimplnxCore::WriteAvizoUniformCoordinate" --verbose
```
Expected: FAIL. **Read every reported difference and classify it before changing anything:**
- expected string wrong → your derivation erred; fix the expectation and note why
- filter output wrong → a genuine finding; record it (Step 4) and raise it

- [ ] **Step 3: Reconcile the oracle**

Apply the **Bug adjudication protocol** above to every mismatch. Legacy source for this filter:
`/Users/mjackson/Workspace/D3D_v6.5.171/DREAM3D/Source/Plugins/ImportExport/ImportExportFilters/AvizoUniformCoordinateWriter.cpp`

Pay particular attention to the `BoundingBox` line. It computes `origin + spacing * dims`, which places the upper bound one full cell beyond the last cell's origin. Decide whether that is correct for the AmiraMesh uniform-coordinate convention, and check whether the legacy writer did the same. If both are wrong in the same way, that is a legacy bug to fix and document; if only SIMPLNX is wrong, it is a port regression.

Re-run until green, with any algorithm fix included.

- [ ] **Step 4: Write the V&V deliverables**

```bash
cd /Users/mjackson/Workspace9/simplnx && python scripts/vv_init.py WriteAvizoUniformCoordinateFilter
```

Fill in the scaffolded `src/Plugins/SimplnxCore/vv/WriteAvizoUniformCoordinateFilter.md` following `docs/vv_templates/report_template.md`, gated by `docs/vv_templates/report_gates.md`. Specifics for this filter:

- **Oracle (confirmed):** `Class 1 analytical` — expected file derived from the AmiraMesh ASCII layout emitted by `generateHeader`/`writeData`, on a 3×2×2 fixture.
- **Legacy comparison:** the executable A/B against 6.5.171 was **not run** — per §2 of `vv_program_plan.md` this filter is Tier B and the A/B is omitted. State that explicitly rather than leaving the section blank. If Step 3 produced any Verdict B finding, also state that the legacy source *was* read to classify it, and cite the file.
- **Exemplar archive:** none — the oracle is inline. Say so; do not create a provenance sidecar for an archive that no longer exists.
- **Status:** `READY FOR REVIEW` on the `| Status |` line, key **not** bold-wrapped — `| Status | READY FOR REVIEW |`. A bold-wrapped key breaks the fleet-wide status grep.
- Record every `filter finding` from Step 3 in `vv/deviations/WriteAvizoUniformCoordinateFilter.md`. If there were none, say "No deviations observed" — do not delete the file.

- [ ] **Step 5: Commit**

```bash
# Include the algorithm/filter source if Step 3 produced a fix:
#   src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/WriteAvizoUniformCoordinate.cpp
#   src/Plugins/SimplnxCore/src/SimplnxCore/utils/AvizoWriter.cpp
git add src/Plugins/SimplnxCore/test/WriteAvizoUniformCoordinateTest.cpp \
        src/Plugins/SimplnxCore/vv/WriteAvizoUniformCoordinateFilter.md \
        src/Plugins/SimplnxCore/vv/deviations/WriteAvizoUniformCoordinateFilter.md
git commit -m "$(cat <<'EOF'
VV: Write Avizo Uniform Coordinate Class 1 oracle

Replaces the 6.6-derived exemplar comparison with an inline analytical
oracle on a 3x2x2 fixture, removing the dependence on golden data that
was never itself verified.


Signed-off-by: YOUR_NAME <YOUR_EMAIL>   # from `git config user.name` / `user.email`
EOF
)"
```

---

## Task 3: WriteAvizoRectilinearCoordinate Class 1 oracle

**Files:**
- Modify: `src/Plugins/SimplnxCore/test/WriteAvizoRectilinearCoordinateTest.cpp`
- Read for the oracle: `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/WriteAvizoRectilinearCoordinate.cpp`

**Interfaces:**
- Consumes: `UnitTest::CompareTextFilesExact` (Task 1); the `AvizoFixture` pattern from Task 2 — **copy it into this test file**, do not share it across translation units.

### Oracle derivation

The rectilinear header differs from the uniform one in exactly these ways — verify each against the source rather than trusting this list:

1. Adds `define Coordinates %llu` where the count is `dims[0] + dims[1] + dims[2]` (for 3×2×2 → `7`), followed by a blank line.
2. `CoordType` is `"rectilinear"`, not `"uniform"`.
3. There is **no** `Content` line and **no** `BoundingBox` line.
4. Adds a second data marker: `Coordinates { float xyz } = @2`.
5. The `@1` marker line carries a comment: `@1 # FeatureIds in z, y, x with X moving fastest, then Y, then Z`.
6. There is a `@2` section holding the coordinate values. Read `writeData` to determine the exact ordering and formatting, and derive the 7 expected values for the fixture.

- [ ] **Step 1: Write the failing test**

Copy the `AvizoFixture` namespace from Task 2 Step 1 into this file verbatim (3×2×2, origin 0, spacing 1, FeatureIds 1..12), then:

```cpp
TEST_CASE("SimplnxCore::WriteAvizoRectilinearCoordinateFilter: Class 1 ASCII oracle", "[SimplnxCore][WriteAvizoRectilinearCoordinateFilter]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure = AvizoFixture::CreateFixture();
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(AvizoFixture::GeomPath()));

  const fs::path outputPath(fmt::format("{}/avizo_rectilinear_class1.am", unit_test::k_BinaryTestOutputDir));

  WriteAvizoRectilinearCoordinateFilter filter;
  Arguments args;
  args.insertOrAssign(WriteAvizoRectilinearCoordinateFilter::k_OutputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(outputPath));
  args.insertOrAssign(WriteAvizoRectilinearCoordinateFilter::k_WriteBinaryFile_Key, std::make_any<bool>(false));
  args.insertOrAssign(WriteAvizoRectilinearCoordinateFilter::k_GeometryPath_Key, std::make_any<DataPath>(AvizoFixture::GeomPath()));
  args.insertOrAssign(WriteAvizoRectilinearCoordinateFilter::k_FeatureIdsArrayPath_Key, std::make_any<DataPath>(AvizoFixture::FeatureIdsPath()));
  args.insertOrAssign(WriteAvizoRectilinearCoordinateFilter::k_Units_Key, std::make_any<StringParameter::ValueType>("microns"));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)

  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)

  // Derived by hand from WriteAvizoRectilinearCoordinate::generateHeader/writeData.
  // `define Coordinates` is dims[0] + dims[1] + dims[2] = 3 + 2 + 2 = 7.
  // The @2 block emits, per axis d, `origin[d] + spacing[d] * i` for i in [0, dims[d]),
  // each formatted "%f " (note the trailing space before each newline).
  const std::vector<std::string> expectedLines = {
      "# AmiraMesh 3D ASCII 2.0",
      "",
      "# Dimensions in x-, y-, and z-direction",
      "define Lattice 3 2 2",
      "define Coordinates 7",
      "",
      "Parameters {",
      "     DREAM3DParams {",
      R"(<REGEX>         Author "DREAM3D-NX SimplnxCore Version .*",</REGEX>)",
      R"(<REGEX>         DateTime ".*"</REGEX>)",
      R"(         FeatureIds Path "Image/CellData/FeatureIds")",
      "     }",
      "     Units {",
      R"(         Coordinates "microns")",
      "     }",
      R"(     CoordType "rectilinear")",
      "}",
      "",
      "Lattice { int FeatureIds } = @1",
      "Coordinates { float xyz } = @2",
      "",
      "# Data section follows",
      "@1 # FeatureIds in z, y, x with X moving fastest, then Y, then Z",
      "1 2 3 4 5 6 7 8 9 10 11 12 ",
      "@2 # x coordinates, then y, then z",
      "0.000000 1.000000 2.000000 ",
      "0.000000 1.000000 ",
      "0.000000 1.000000 ",
  };

  UnitTest::CompareTextFilesExact(outputPath, expectedLines);
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}
```

Note there is **no** `Content` line and **no** `BoundingBox` line in the rectilinear header — confirm that against the source rather than trusting this listing.

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/mjackson/Workspace2/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest \
  && ctest -R "SimplnxCore::WriteAvizoRectilinearCoordinate" --verbose
```
Expected: FAIL.

- [ ] **Step 3: Reconcile the oracle**

Apply the **Bug adjudication protocol** above to every mismatch. Legacy source for this filter:
`/Users/mjackson/Workspace/D3D_v6.5.171/DREAM3D/Source/Plugins/ImportExport/ImportExportFilters/AvizoRectilinearCoordinateWriter.cpp`

The `@2` coordinate block emits `origin[d] + spacing[d] * i` for `i` in `[0, dims[d])` — node positions, not cell centres. Check that against the legacy writer and against what AmiraMesh expects for a rectilinear lattice.

Re-run until green, with any algorithm fix included.

- [ ] **Step 4: Write the V&V deliverables**

```bash
python scripts/vv_init.py WriteAvizoRectilinearCoordinateFilter
```
Fill in as in Task 2 Step 4.

- [ ] **Step 5: Retire the shared archive and commit**

Both Avizo tests are now archive-free. Remove the `6_6_avizo_writers.tar.gz` `download_test_data()` block from `src/Plugins/SimplnxCore/test/CMakeLists.txt`, then confirm nothing else references it:

```bash
cd /Users/mjackson/Workspace9/simplnx && grep -rn "6_6_avizo_writers" src/ || echo "NO REMAINING REFERENCES"
```

```bash
# Include src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/WriteAvizoRectilinearCoordinate.cpp
# if Step 3 produced a fix.
git add src/Plugins/SimplnxCore/test/WriteAvizoRectilinearCoordinateTest.cpp \
        src/Plugins/SimplnxCore/test/CMakeLists.txt \
        src/Plugins/SimplnxCore/vv/WriteAvizoRectilinearCoordinateFilter.md \
        src/Plugins/SimplnxCore/vv/deviations/WriteAvizoRectilinearCoordinateFilter.md
git commit -m "$(cat <<'EOF'
VV: Write Avizo Rectilinear Coordinate Class 1 oracle

Retires the 6_6_avizo_writers archive; both Avizo writers now assert
against inline analytical oracles instead of 6.6-era golden files.

Signed-off-by: YOUR_NAME <YOUR_EMAIL>   # from `git config user.name` / `user.email`
EOF
)"
```

---

## Task 4: WriteSPParksSites Class 1 oracle

**Files:**
- Modify: `src/Plugins/SimplnxCore/test/WriteSPParksSitesTest.cpp`
- Read for the oracle: `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/WriteSPParksSites.cpp`

### Oracle derivation

This is the simplest of the four — no volatile lines at all, so the whole file is literal. From `WriteHeader` and `WriteFile`, for a 3×2×2 fixture with `FeatureIds = 1..12`:

- Header is a `-` line, `3 dimension`, `<totalPoints> sites`, `26 max neighbors`, then `0 <dimX> xlo xhi`, `0 <dimY> ylo yhi`, `0 <dimZ> zlo zhi`, a blank line, `Values`, a blank line.
- Body is one line per cell: `<k+1> <featureIds[k]>`.
- Note the header reports `dims` — the cell counts — as the upper bounds, ignoring origin and spacing entirely. Judge whether that is right for the SPParks site format before accepting it; if not, it is a deviation.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("SimplnxCore::WriteSPParksSitesFilter: Class 1 oracle", "[SimplnxCore][WriteSPParksSitesFilter]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure = SPParksFixture::CreateFixture();
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(SPParksFixture::GeomPath()));

  const fs::path outputPath(fmt::format("{}/spparks_sites_class1.txt", unit_test::k_BinaryTestOutputDir));

  WriteSPParksSitesFilter filter;
  Arguments args;
  args.insertOrAssign(WriteSPParksSitesFilter::k_OutputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(outputPath));
  args.insertOrAssign(WriteSPParksSitesFilter::k_ImageGeomPath_Key, std::make_any<DataPath>(SPParksFixture::GeomPath()));
  args.insertOrAssign(WriteSPParksSitesFilter::k_FeatureIdsArrayPath_Key, std::make_any<DataPath>(SPParksFixture::FeatureIdsPath()));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)

  const std::vector<std::string> expectedLines = {
      "-", "3 dimension", "12 sites", "26 max neighbors",
      "0 3 xlo xhi", "0 2 ylo yhi", "0 2 zlo zhi",
      "", "Values", "",
      "1 1", "2 2", "3 3", "4 4", "5 5", "6 6",
      "7 7", "8 8", "9 9", "10 10", "11 11", "12 12",
  };

  UnitTest::CompareTextFilesExact(outputPath, expectedLines);
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}
```

Prepend this fixture to the file:

```cpp
namespace
{
namespace SPParksFixture
{
constexpr usize k_XDim = 3;
constexpr usize k_YDim = 2;
constexpr usize k_ZDim = 2;
constexpr usize k_CellCount = k_XDim * k_YDim * k_ZDim;

inline DataStructure CreateFixture()
{
  DataStructure dataStructure;
  auto* imageGeomPtr = ImageGeom::Create(dataStructure, "Image");
  imageGeomPtr->setDimensions({k_XDim, k_YDim, k_ZDim});
  imageGeomPtr->setOrigin({0.0F, 0.0F, 0.0F});
  imageGeomPtr->setSpacing({1.0F, 1.0F, 1.0F});

  auto* cellAmPtr = AttributeMatrix::Create(dataStructure, "CellData", {k_ZDim, k_YDim, k_XDim}, imageGeomPtr->getId());
  imageGeomPtr->setCellData(*cellAmPtr);

  auto* featureIdsPtr = Int32Array::CreateWithStore<DataStore<int32>>(dataStructure, "FeatureIds", {k_ZDim, k_YDim, k_XDim}, {1}, cellAmPtr->getId());
  auto& featureIdsRef = featureIdsPtr->getDataStoreRef();
  for(usize i = 0; i < k_CellCount; i++)
  {
    featureIdsRef[i] = static_cast<int32>(i + 1);
  }

  return dataStructure;
}

inline DataPath GeomPath()
{
  return DataPath({"Image"});
}

inline DataPath FeatureIdsPath()
{
  return DataPath({"Image", "CellData", "FeatureIds"});
}
} // namespace SPParksFixture
} // namespace
```

**Verify the actual parameter key names** in `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/WriteSPParksSitesFilter.hpp` — the keys used in the `TEST_CASE` above are inferred from the algorithm's `InputValues` struct and may not match the filter's declared `k_*_Key` constants.

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/mjackson/Workspace2/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest \
  && ctest -R "SimplnxCore::WriteSPParksSites" --verbose
```
Expected: FAIL.

- [ ] **Step 3: Reconcile the oracle**

Apply the **Bug adjudication protocol** above to every mismatch. Legacy source for this filter:
`/Users/mjackson/Workspace/D3D_v6.5.171/DREAM3D/Source/Plugins/ImportExport/ImportExportFilters/SPParksSitesWriter.cpp`

The header reports cell counts as the `xhi`/`yhi`/`zhi` bounds and ignores origin and spacing entirely. For a geometry with non-zero origin or non-unit spacing that is very likely wrong. Build a second fixture with origin `(10, 20, 30)` and spacing `(0.5, 0.5, 0.5)` to find out, then check the legacy writer for the same behaviour.

Re-run until green, with any algorithm fix included.

- [ ] **Step 4: Write the V&V deliverables and retire the archive**

```bash
python scripts/vv_init.py WriteSPParksSitesFilter
grep -rn "6_5_spparks_sites_writer" src/ || echo "NO REMAINING REFERENCES"
```

Remove the `6_5_spparks_sites_writer.tar.gz` entry from `src/Plugins/SimplnxCore/test/CMakeLists.txt` **only if** the grep comes back empty. Note in the report that this filter's prior exemplar was 6.5-derived rather than 6.6 — the same circularity argument applies.

- [ ] **Step 5: Commit**

```bash
# Include src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/WriteSPParksSites.cpp
# if Step 3 produced a fix.
git add src/Plugins/SimplnxCore/test/WriteSPParksSitesTest.cpp \
        src/Plugins/SimplnxCore/test/CMakeLists.txt \
        src/Plugins/SimplnxCore/vv/WriteSPParksSitesFilter.md \
        src/Plugins/SimplnxCore/vv/deviations/WriteSPParksSitesFilter.md
git commit -m "$(cat <<'EOF'
VV: Write SPParks Sites Class 1 oracle

Signed-off-by: YOUR_NAME <YOUR_EMAIL>   # from `git config user.name` / `user.email`
EOF
)"
```

---

## Task 5: WriteLosAlamosFFT Class 1 oracle

**Files:**
- Modify: `src/Plugins/SimplnxCore/test/WriteLosAlamosFFTTest.cpp`
- Read for the oracle: `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/WriteLosAlamosFFT.cpp`

This filter needs a richer fixture: it consumes Euler angles and cell phases in addition to FeatureIds.

### Oracle derivation

The body is a single `fmt::format` per cell:

```cpp
phi1 = cellEulerAngles[index * 3]     * 180.0f * Constants::k_1OverPiF;
phi  = cellEulerAngles[index * 3 + 1] * 180.0f * Constants::k_1OverPiF;
phi2 = cellEulerAngles[index * 3 + 2] * 180.0f * Constants::k_1OverPiF;
file << fmt::format("{:.3f} {:.3f} {:.3f} {} {} {} {} {}\n",
                    phi1, phi, phi2, x + 1, y + 1, z + 1,
                    featureIds.getValue(index), cellPhases.getValue(index));
```

Choose Euler values whose degree conversion is exact at three decimals so the expected strings are unambiguous. Radians of `0`, `π/2`, `π` convert to `0.000`, `90.000`, `180.000`. Use `Constants::k_PiOver2F` and `Constants::k_PiF` in the fixture rather than literal `1.5707964F`, so the test does not encode a rounding assumption.

There is **no header** — the source comments that the original filter left `writeHeader` unimplemented. Confirm that and note it in the report; a downstream consumer expecting a header would be a real compatibility concern worth flagging.

Iteration order is z-outer, y-middle, x-inner, and the emitted indices are 1-based.

- [ ] **Step 1: Write the failing test**

Build a 2×2×1 fixture (4 cells — keeps the expected block short):

```cpp
namespace
{
namespace FftFixture
{
constexpr usize k_XDim = 2;
constexpr usize k_YDim = 2;
constexpr usize k_ZDim = 1;
constexpr usize k_CellCount = k_XDim * k_YDim * k_ZDim;

inline DataStructure CreateFixture()
{
  DataStructure dataStructure;
  auto* imageGeomPtr = ImageGeom::Create(dataStructure, "Image");
  imageGeomPtr->setDimensions({k_XDim, k_YDim, k_ZDim});
  imageGeomPtr->setOrigin({0.0F, 0.0F, 0.0F});
  imageGeomPtr->setSpacing({1.0F, 1.0F, 1.0F});

  auto* cellAmPtr = AttributeMatrix::Create(dataStructure, "CellData", {k_ZDim, k_YDim, k_XDim}, imageGeomPtr->getId());
  imageGeomPtr->setCellData(*cellAmPtr);

  auto* eulersPtr = Float32Array::CreateWithStore<DataStore<float32>>(dataStructure, "EulerAngles", {k_ZDim, k_YDim, k_XDim}, {3}, cellAmPtr->getId());
  auto& eulersRef = eulersPtr->getDataStoreRef();
  const std::array<float32, 12> eulerValues = {
      0.0F,                     0.0F,                     0.0F,
      Constants::k_PiOver2F,    0.0F,                     0.0F,
      0.0F,                     Constants::k_PiOver2F,    0.0F,
      0.0F,                     0.0F,                     Constants::k_PiF,
  };
  for(usize i = 0; i < eulerValues.size(); i++)
  {
    eulersRef[i] = eulerValues[i];
  }

  auto* phasesPtr = Int32Array::CreateWithStore<DataStore<int32>>(dataStructure, "Phases", {k_ZDim, k_YDim, k_XDim}, {1}, cellAmPtr->getId());
  auto* featureIdsPtr = Int32Array::CreateWithStore<DataStore<int32>>(dataStructure, "FeatureIds", {k_ZDim, k_YDim, k_XDim}, {1}, cellAmPtr->getId());
  auto& phasesRef = phasesPtr->getDataStoreRef();
  auto& featureIdsRef = featureIdsPtr->getDataStoreRef();
  for(usize i = 0; i < k_CellCount; i++)
  {
    phasesRef[i] = 1;
    featureIdsRef[i] = static_cast<int32>(i + 1);
  }

  return dataStructure;
}
} // namespace FftFixture
} // namespace
```

Expected block, derived from the format string and the z/y/x iteration order:

```cpp
const std::vector<std::string> expectedLines = {
    "0.000 0.000 0.000 1 1 1 1 1",
    "90.000 0.000 0.000 2 1 1 2 1",
    "0.000 90.000 0.000 1 2 1 3 1",
    "0.000 0.000 180.000 2 2 1 4 1",
};
```

And the test itself:

```cpp
TEST_CASE("SimplnxCore::WriteLosAlamosFFTFilter: Class 1 oracle", "[SimplnxCore][WriteLosAlamosFFTFilter]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure = FftFixture::CreateFixture();
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(DataPath({"Image"})));

  const fs::path outputPath(fmt::format("{}/los_alamos_fft_class1.txt", unit_test::k_BinaryTestOutputDir));

  WriteLosAlamosFFTFilter filter;
  Arguments args;
  args.insertOrAssign(WriteLosAlamosFFTFilter::k_OutputFile_Key, std::make_any<FileSystemPathParameter::ValueType>(outputPath));
  args.insertOrAssign(WriteLosAlamosFFTFilter::k_ImageGeomPath_Key, std::make_any<DataPath>(DataPath({"Image"})));
  args.insertOrAssign(WriteLosAlamosFFTFilter::k_FeatureIdsArrayPath_Key, std::make_any<DataPath>(DataPath({"Image", "CellData", "FeatureIds"})));
  args.insertOrAssign(WriteLosAlamosFFTFilter::k_CellEulerAnglesArrayPath_Key, std::make_any<DataPath>(DataPath({"Image", "CellData", "EulerAngles"})));
  args.insertOrAssign(WriteLosAlamosFFTFilter::k_CellPhasesArrayPath_Key, std::make_any<DataPath>(DataPath({"Image", "CellData", "Phases"})));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)

  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)

  // Derived from the single fmt::format in WriteLosAlamosFFT::operator().
  // Fields: phi1 phi phi2 (degrees, 3dp) x+1 y+1 z+1 featureId phase,
  // iterated z-outer / y-middle / x-inner.
  const std::vector<std::string> expectedLines = {
      "0.000 0.000 0.000 1 1 1 1 1",
      "90.000 0.000 0.000 2 1 1 2 1",
      "0.000 90.000 0.000 1 2 1 3 1",
      "0.000 0.000 180.000 2 2 1 4 1",
  };

  UnitTest::CompareTextFilesExact(outputPath, expectedLines);
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}
```

**Verify the parameter key names** against `src/Plugins/SimplnxCore/src/SimplnxCore/Filters/WriteLosAlamosFFTFilter.hpp` — the `k_*_Key` constants above are inferred from the algorithm's `InputValues` field names.

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/mjackson/Workspace2/DREAM3D-Build/simplnx-Rel && cmake --build . --target SimplnxCoreUnitTest \
  && ctest -R "SimplnxCore::WriteLosAlamosFFT" --verbose
```
Expected: FAIL.

- [ ] **Step 3: Reconcile the oracle**

Apply the **Bug adjudication protocol** above to every mismatch. Legacy source for this filter:
`/Users/mjackson/Workspace/D3D_v6.5.171/DREAM3D/Source/Plugins/ImportExport/ImportExportFilters/LosAlamosFFTWriter.cpp`

Watch whether `180.000` emerges as `180.000` or `179.999`. If the latter, that is a float32 precision property of the `radians * 180 * (1/pi)` conversion — record it in the report's Oracle section rather than loosening the expectation. Compare against how the legacy writer performed the same conversion; if it used `double`, the difference is a precision deviation worth documenting even though neither is a bug.

- [ ] **Step 4: Write the V&V deliverables and check the archives**

```bash
python scripts/vv_init.py WriteLosAlamosFFTFilter
grep -rn "LosAlamosFFTExemplar\|bin_feature_phases" src/ || echo "NO REMAINING REFERENCES"
```

`bin_feature_phases.tar.gz` is also used by `ComputeFeaturePhasesBinaryTest` — **do not remove it** unless the grep shows this test was its only consumer. `LosAlamosFFTExemplar.tar.gz` should become removable.

- [ ] **Step 5: Commit**

```bash
# Include src/Plugins/SimplnxCore/src/SimplnxCore/Filters/Algorithms/WriteLosAlamosFFT.cpp
# if Step 3 produced a fix.
git add src/Plugins/SimplnxCore/test/WriteLosAlamosFFTTest.cpp \
        src/Plugins/SimplnxCore/test/CMakeLists.txt \
        src/Plugins/SimplnxCore/vv/WriteLosAlamosFFTFilter.md \
        src/Plugins/SimplnxCore/vv/deviations/WriteLosAlamosFFTFilter.md
git commit -m "$(cat <<'EOF'
VV: Write Los Alamos FFT Class 1 oracle

Signed-off-by: YOUR_NAME <YOUR_EMAIL>   # from `git config user.name` / `user.email`
EOF
)"
```

---

## Task 6: Report findings and correct the program plan

**Files:**
- Modify: `docs/vv_templates/vv_program_plan.md` §3.3

- [ ] **Step 1: Run the full SimplnxCore suite**

```bash
cd /Users/mjackson/Workspace2/DREAM3D-Build/simplnx-Rel && ctest -R "SimplnxCore::" --verbose 2>&1 | tail -40
```
Expected: no new failures. If removing a `download_test_data()` entry broke an unrelated test, restore that entry — a shared archive is not retirable until every consumer is converted.

- [ ] **Step 2: Update the §3.3 writer roster**

Replace the nine-writer "2a — the writers" list with the four-in / five-deferred split from the **Scope** section of this plan, and state the reason: `WriteStlFile`, `WriteAbaqusHexahedron`, `WriteGBCDGMTFile`, `WriteINLFile`, and `WriteGBCDTriangleData` are not cheap Class 1 targets and each needs its own plan.

- [ ] **Step 3: Record the defect hit rate**

`vv_program_plan.md` §3.3 asks Phase 2 to report its hit rate. Add a line under §3.3 giving the count of `filter finding` verdicts across all four tasks, out of four filters converted. **This number is the point of the phase** — it is the evidence for how much Phases 3 and 4 are worth. If it is zero, say zero; a truthful zero is a useful result.

- [ ] **Step 4: Commit**

```bash
git add docs/vv_templates/vv_program_plan.md
git commit -m "$(cat <<'EOF'
DOC: Record Phase 2a writer results and correct the §3.3 roster

Signed-off-by: YOUR_NAME <YOUR_EMAIL>   # from `git config user.name` / `user.email`
EOF
)"
```

---

## Definition of done

- [ ] Four writer tests assert against inline analytical oracles with zero `.tar.gz` dependencies.
- [ ] `grep -rn "6_6_avizo_writers" src/` returns nothing.
- [ ] Each of the four filters has a `vv/<Filter>.md` and `vv/deviations/<Filter>.md`, `| Status |` un-bolded, status `READY FOR REVIEW`.
- [ ] `ctest -R "SimplnxCore::" --verbose` shows no new failures.
- [ ] `vv_program_plan.md` §3.3 records the corrected roster and the defect hit rate.
- [ ] Every SIMPL `BackwardsCompatibility` `TEST_CASE` is preserved untouched.
- [ ] Every Verdict B finding is **fixed in this PR**, not deferred.
- [ ] Every Verdict B finding that 6.5.171 also exhibited has a full deviation entry with `Root cause: Bug in 6.5.171`, a legacy file/line citation, and `Recommendation: Trust SIMPLNX`.
- [ ] Every Verdict B finding that 6.5.171 did *not* exhibit is described in the report Summary as a port regression, with no deviation entry.

## Open questions for the requester

1. **Is `test/UnitTestCommon` covered by its own ctest target?** Task 1 Step 1 assumes not and offers a fallback. Confirm before starting.
2. `WriteLosAlamosFFT` emits no header at all. If any downstream consumer expects one, that is a compatibility bug rather than a documentation note — worth a decision before the report is written.
