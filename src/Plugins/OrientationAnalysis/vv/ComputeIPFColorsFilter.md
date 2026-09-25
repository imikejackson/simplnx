# V&V Report: ComputeIPFColorsFilter

|           |                          |
|-----------|--------------------------|
| Plugin    | OrientationAnalysis      |
| SIMPLNX UUID | `64cb4f27-6e5e-4dd2-8a03-0c448cb8f5e6` |
| SIMPLNX Human Name | Compute IPF Colors |
| DREAM3D 6.5.171 equivalent | `GenerateIPFColors` (SIMPL UUID `a50e6532-8075-5de5-ab63-945feb0de7f7`) — `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/GenerateIPFColors.{h,cpp}` |
| Verified commit | `a307946e7` (v7.4.2 release) |
| Status | COMPLETE     |
| Sign-off | Michael Jackson <mike.jackson@bluequartz.net> — 2026-07-16 |
| Second-engineer sign-off | Matthew Marine — 2026-07-13 (approving reviewer, PR #1666). Supersedes the 2026-07-16 technical-authority self-sign-off. |

## At a glance

| Aspect                 | Current state            |
|------------------------|--------------------------|
| Algorithm Relationship | **Port** of DREAM3D 6.5.171 `GenerateIPFColors`; material deltas are the EbsdLib color library, three color-key choices, bool-or-uint8 masks, and cancel checks. |
| Oracle (confirmed)     | **Classes 1 and 4** verify orchestration and invariants; **Classes 2 and 3** verify EbsdLib routing and the standard red IPF corner. All 8 tests in `test/ComputeIPFColorsTest.cpp` pass. |
| Code paths enumerated  | 16 of 18 exercised; the mid-loop cancel branch and unreachable color-key default are not exercised. |
| Tests today            | 8 test cases cover the analytical fixture, mask variants, reference-direction normalization, errors, color-key wiring, and SIMPL conversion. |
| Exemplar archive       | None; the oracle is inline. `so3_cubic_high_ipf_001.tar.gz` is retired as an oracle but remains available to another test. |
| Legacy comparison      | **Run** — SIMPLNX and the stored legacy array are byte-identical; a fresh official DREAM3D 6.5.171 run differs in 14 of 343,963 cells by ±1/255 in one channel (D1). |
| Bug flags              | None. |
| V&V phase | **COMPLETE.** |

## Summary

`ComputeIPFColorsFilter` assigns each cell an inverse-pole-figure RGB color by delegating the color computation to EbsdLib `LaueOps::generateIPFColor()` for that cell's phase/crystal structure, with optional masking (bad cells → black) and a selectable color key. Because the color math is EbsdLib's (and is covered by EbsdLib's own `TSLColorKeyTest` / `PUCMColorKeyTest` / `NolzeHielscherColorKeyTest` / `ColorKeyKindTest`), V&V here targets only the SIMPLNX **value-add** — per-cell dispatch/indexing, masking, reference-direction normalization, phase-bounds handling, and output packing — using a Class 1/4 analytical-orchestration oracle, a Class 2 in-process EbsdLib cross-check, and a Class 3 standard-IPF-corner anchor. SIMPLNX matched the oracle on the first run (no filter bugs found) and is byte-identical to the stored legacy reference; the only difference vs a fresh 6.5.171 build is ±1/255 quantization jitter on 0.004% of cells (deviation D1).

## Algorithm Relationship

*Classification:* **Port**

*Evidence:* The SIMPLNX algorithm `Algorithms/ComputeIPFColors.cpp` (203 lines) is a direct translation of `GenerateIPFColors::execute()` / `GenerateIPFColorsImpl::convert()` from DREAM3D 6.5.171. The SIMPL UUID `a50e6532-8075-5de5-ab63-945feb0de7f7` is retained via `OrientationAnalysisLegacyUUIDMapping.hpp` plus SIMPL 6.4/6.5 conversion fixtures. The control flow is preserved verbatim: per-cell init-to-black → mask lookup → phase-bounds sanity → `phase < numPhases && calcIPF && crystalStructures[phase] < LaueGroupEnd` coloring guard → `generateIPFColor` → RGB write, followed by the post-loop `-48000` phase-warning error.

*Port-time deltas (each with its effect on output):*

1. **Color library: OrientationLib → EbsdLib.** The `generateIPFColor` call now resolves to EbsdLib. Same standard conventions; produces the ±1/255 quantization jitter documented in `ComputeIPFColorsFilter-D1`.
2. **`Color Key` choice added (PR #1631).** New `ChoicesParameter` routing to `ebsdlib::ColorKeyKind` {TSL, PUCM, Nolze-Hielscher}. Legacy was TSL-only; the default (index 0 = TSL) reproduces legacy behavior. PUCM/Nolze-Hielscher have no legacy equivalent.
3. **Mask accepts bool *or* uint8.** Legacy required a `bool` mask; SIMPLNX `run()` dispatches `convert<bool>`/`convert<uint8>`. No effect on output values; widens accepted input.
4. **Cancel checks added.** `shouldCancel()` guard at the top of the per-cell loop. UX-only.
5. **Reference direction normalization** is unchanged in intent (`FloatVec3::normalize()` vs legacy `MatrixMath::Normalize3x1`); both make the reference direction unit length before use.

*Material PRs since baseline:* #1631 ("EbsdLib 3.0.0 + V&V of 6 Filters", added the Color Key), #1472 (EbsdLib 2.0.0 API), #1438 (microtexture cleanup), #1501 (Vec3 unification). None alter the coloring logic beyond the deltas above.

## Oracle

*Class:* **1 (Analytical)** + **4 (Invariant)** primary, **2 (Reference — EbsdLib)** and **3 (Paper/standard-IPF)** companions.

Per the "test the value-add, not upstream" principle: EbsdLib is the trusted reference for the color math (it has its own test suite), so the oracle verifies that SIMPLNX *routes data correctly into and out of EbsdLib*, not that the color algorithm is correct.

*Applied:*
- **Class 1 (Analytical, orchestration):** a masked-off cell is exactly `(0,0,0)`; a cell whose crystal structure is `≥ LaueGroupEnd` (999/Unknown) is exactly `(0,0,0)`; a non-unit reference direction `[0,0,5]` yields identical output to `[0,0,1]` (normalization); a cell whose phase index `≥ numPhases` produces error `-48000`. These are pure filter logic, derivable without reference to any DREAM3D implementation.
- **Class 4 (Invariant):** output is a 3-component uint8 array with the parent tuple count; non-colored cells are exactly black; colored cells are non-black.
- **Class 2 (Reference, EbsdLib):** for every colored cell, the filter output equals a direct, independent in-process `ebsdlib::LaueOps::GetAllOrientationOps()[crystalStruct]->generateIPFColor(euler, refDir, false, kind)` call — verifying the phase→crystal-structure→ops indexing, the 3-tuple Euler slicing, the reference-direction hand-off, and the color-key pass-through.
- **Class 3 (Standard IPF):** an identity-orientation (`(0,0,0)`) cubic cell viewed down `[001]` is the red corner of the standard IPF triangle. EbsdLib's own `TSLColorKeyTest` fixes `[001] → r=1, g=0, b=0`; the test asserts R ≥ 250, G = 0, B = 0.

*Toy data:* built inline in C++ (`BuildAnalyticalDataset`) — a 6-cell ImageGeom, ensemble `CrystalStructures = [999 Unknown, 1 Cubic_High, 0 Hexagonal_High]`, cells chosen to hit identity-cubic, arbitrary cubic, arbitrary hex (two Laue classes), a masked cell, and an invalid-crystal-structure cell.

*Encoded:* `test/ComputeIPFColorsTest.cpp` — 8 `TEST_CASE`s, all pass. The main `Class 1/2/3 Oracle (inline analytical dataset)` case carries the Class 2/3/1/4 assertions across four `SECTION`s.

*Second-engineer review:* **Matthew Marine — 2026-07-13 (approving reviewer, PR #1666).** Supersedes the 2026-07-16 technical-authority self-sign-off.

## Bugs found and fixed

None.

## Code path coverage

*16 of 18 code paths exercised. The 2 uncovered paths are the mid-loop cancel branch (requires injecting a cancel signal) and the color-key `default` error `-23510` (unreachable through normal preflight because `ChoicesParameter` already constrains the value to 0–2). Both are low-value.*

Source: `src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/Algorithms/ComputeIPFColors.cpp` (203 lines).

Logical phases: (a) `operator()` setup + parallel dispatch, (b) per-cell `convert()` loop, (c) post-loop phase-warning error; plus (d) filter-level preflight / executeImpl wiring.

| #  | Phase        | Path                                                                                          | Test case                                                                 |
|----|--------------|-----------------------------------------------------------------------------------------------|---------------------------------------------------------------------------|
| 1  | (a) Dispatch | mask present & `boolean` → `convert<bool>`                                                     | `Class 1/2/3 Oracle` (bool `Mask`)                                        |
| 2  | (a) Dispatch | mask present & `uint8` → `convert<uint8>`                                                      | `uint8 mask array drives the black-out path`                              |
| 3  | (a) Dispatch | mask null (`useMask=false`) → `convert<bool>` with null mask                                   | `no-mask path colors every valid cell`                                    |
| 4  | (a) Setup    | reference direction normalized before use                                                     | `reference direction is normalized`                                       |
| 5  | (b) Per-cell | `shouldCancel()` → early return                                                               | *Not directly tested. Requires mid-execution cancel-signal injection; low-value UX guard.* |
| 6  | (b) Per-cell | init each cell to `(0,0,0)`                                                                    | `Class 1/2/3 Oracle` — masked/invalid cells remain black                  |
| 7  | (b) Per-cell | `maskArray != null` → `calcIPF = mask[i]`                                                      | `Class 1/2/3 Oracle`, `uint8 mask ...`                                    |
| 8  | (b) Per-cell | `phase >= numPhases` → `incrementPhaseWarningCount()`                                          | `phase index out of range returns -48000`                                 |
| 9  | (b) Per-cell | coloring branch: `phase<numPhases && calcIPF && CS<LaueGroupEnd` → `generateIPFColor` + write  | `Class 1/2/3 Oracle` — Class 2 (all colored cells) + Class 3 (cell 0)     |
| 10 | (b) Per-cell | not colored — `calcIPF==false` (masked) → stays black                                         | `Class 1/2/3 Oracle` (cell 3), `uint8 mask ...`                           |
| 11 | (b) Per-cell | not colored — `CS >= LaueGroupEnd` (999) → stays black                                         | `Class 1/2/3 Oracle` (cell 4)                                             |
| 12 | (b) Per-cell | not colored — `phase >= numPhases` → stays black                                              | `phase index out of range returns -48000`                                 |
| 13 | (c) Finalize | `m_PhaseWarningCount > 0` → return error `-48000`                                              | `phase index out of range returns -48000`                                 |
| 14 | (c) Finalize | `m_PhaseWarningCount == 0` → return success `{}`                                               | every passing test                                                        |
| 15 | (d) Filter   | `Color Key` switch 0/1/2 → TSL/PUCM/Nolze-Hielscher                                            | `ColorKey choice reaches algorithm`                                       |
| 16 | (d) Filter   | `Color Key` `default` → error `-23510`                                                        | *Not directly tested. Unreachable via preflight — `ChoicesParameter` constrains the value to [0,2].* |
| 17 | (d) Filter   | `validateNumberOfTuples` fails → error `-651`                                                 | `Preflight Error - Cell array tuple count mismatch (-651)`                |
| 18 | (d) Filter   | preflight creates the 3-component uint8 output array                                          | every passing test (preflight is `REQUIRE`-valid)                         |

## Test inventory

| Test case | Status | Notes |
|-----------|--------|-------|
| `OrientationAnalysis::ComputeIPFColorsFilter: Class 1/2/3 Oracle (inline analytical dataset)` | new-for-V&V | Replaces the retired circular-oracle test. 4 SECTIONs: Class 3 red-corner, Class 2 EbsdLib cross-check over 4 colored cells, Class 1 masked→black, Class 1 invalid-CS→black. |
| `OrientationAnalysis::ComputeIPFColorsFilter: uint8 mask array drives the black-out path` | new-for-V&V | Exercises `convert<uint8>` dispatch; masked cell black, good cell colored. |
| `OrientationAnalysis::ComputeIPFColorsFilter: no-mask path colors every valid cell` | new-for-V&V | Exercises null-mask `convert<bool>`; previously-masked cell now matches EbsdLib reference; invalid-CS cell still black. |
| `OrientationAnalysis::ComputeIPFColorsFilter: reference direction is normalized` | new-for-V&V | `[0,0,5]` output byte-identical to `[0,0,1]`. |
| `OrientationAnalysis::ComputeIPFColorsFilter: phase index out of range returns -48000` | new-for-V&V | Corrupts a phase to an out-of-range ensemble index; asserts execute error `-48000`. |
| `OrientationAnalysis::ComputeIPFColorsFilter: ColorKey choice reaches algorithm` | kept (retargeted) | Was on the retired exemplar; now runs on the inline dataset. Asserts TSL ≠ PUCM ≠ Nolze-Hielscher output. |
| `OrientationAnalysis::ComputeIPFColorsFilter: Preflight Error - Cell array tuple count mismatch (-651)` | kept | Synthetic mismatched tuple counts → `-651`. |
| `OrientationAnalysis::ComputeIPFColorsFilter: SIMPL Backwards Compatibility` | kept | `DYNAMIC_SECTION` over SIMPL 6.4 + 6.5 conversion fixtures; validates UUID + argument conversion. |
| *(retired)* `OrientationAnalysis::ComputeIPFColors` | retired | Circular oracle — compared filter output against the legacy-produced `IPF Colors` array inside `so3_cubic_high_ipf_001.dream3d` ("produced by SIMPL/DREAM3D … our results should match theirs"). Replaced by the analytical oracle above. |

## Exemplar archive

- **Archive:** None for this filter. The oracle dataset is constructed inline in C++ (`BuildAnalyticalDataset` in `test/ComputeIPFColorsTest.cpp`); no golden `.dream3d` is downloaded or compared, so no provenance sidecar is required.
- **Retired:** `so3_cubic_high_ipf_001.tar.gz` (SHA512 `dfe4598c…dd616b85`) as an oracle for this filter. Its `download_test_data()` entry in `test/CMakeLists.txt` remains because `CreateEnsembleInfoTest` still consumes the archive.

## Deviations from DREAM3D 6.5.171

Comparison run on `so3_cubic_high_ipf_001` (343,963 single-phase cubic cells, TSL, refDir [001]) — see `vv/comparisons/ComputeIPFColorsFilter/`.

NX 7.4.1 with EbsdLib 2.2.0 and NX 7.4.2 with EbsdLib 3.1.2 are byte-identical across all 343,963 RGB tuples. Both reproduce the stored legacy array exactly. Each differs from the fresh official DREAM3D 6.5.171 run in the same 14 tuples.

- `ComputeIPFColorsFilter-D1` — 14/343,963 cells (0.004%) differ from a fresh 6.5.171 run by ±1/255 in one channel (SIMPLNX is byte-identical to the stored legacy reference); precision + library quantization jitter — see `vv/deviations/ComputeIPFColorsFilter.md`.
