# V&V Report: ConvertOrientationsFilter

|           |                          |
|-----------|--------------------------|
| Plugin    | OrientationAnalysis      |
| SIMPLNX UUID | `501e54e6-a66f-4eeb-ae37-00e649c00d4b` |
| SIMPLNX Human Name | Convert Orientation Representation |
| DREAM3D 6.5.171 equivalent | `ConvertOrientations` (SIMPL UUID `e5629880-98c4-5656-82b8-c9fe2b9744de`) — `Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/ConvertOrientations.{h,cpp}`; mapped in `OrientationAnalysisLegacyUUIDMapping.hpp` |
| Verified commit | `a307946e7` (v7.4.2 release) |
| Status | COMPLETE |
| Sign-off | Michael Jackson <mike.jackson@bluequartz.net> — 2026-07-16 |
| Second-engineer sign-off | Jared Duffey — 2026-07-14 (approving reviewer, PR #1648). Supersedes the 2026-07-16 technical-authority self-sign-off. |

## At a glance

| Aspect                 | Current state            |
|------------------------|--------------------------|
| Algorithm Relationship | **Rewrite** of the filter plumbing under the retained SIMPL UUID; SIMPLNX adds Stereographic, accepts float32 only, and removes destructive Euler sanitization. |
| Oracle (confirmed)     | **Class 3** Rowenhorst transformations, **Class 1** Stereographic closed form, and **Class 4** round-trip invariants are encoded in `test/ConvertOrientationsTest.cpp`; 1032 assertions pass. |
| Code paths enumerated  | 7 of 9 reachable paths exercised; 11 paths are enumerated, including 2 unreachable dispatch arms. |
| Tests today            | 5 test cases cover the 56 conversion pairs, striding, Stereographic, invalid preflight, equal representations, and SIMPL conversion. |
| Exemplar archive       | None; dispatch landmarks are inline. Unknown-provenance values were retired and replaced with EbsdLib-derived landmarks plus independent analytical pins. |
| Legacy comparison      | **Run** — official DREAM3D 6.5.171 is bit-identical on 1,288 interior tuples; D5 proves destructive legacy Euler sanitization, D2–D4 remain active, and provenance-invalid D1 is retired. |
| Bug flags              | `ConvertOrientationsFilter-D5` — legacy mutates valid `Φ = π` and out-of-range Euler input; SIMPLNX preserves the supplied tuple. |
| V&V phase | **COMPLETE.** |

For worked instances see `src/Plugins/OrientationAnalysis/vv/BadDataNeighborOrientationCheckFilter.md` and `src/Plugins/OrientationAnalysis/vv/ComputeAvgCAxesFilter.md`.

## Summary

`ConvertOrientationsFilter` ("Convert Orientation Representation") converts a Float32 orientation array from any one of 8 representations (Euler, Orientation Matrix, Quaternion, Axis-Angle, Rodrigues, Homochoric, Cubochoric, Stereographic) to any other, applying the conversion per-tuple in parallel. The actual transformation equations are delegated to EbsdLib (Rowenhorst 2015) and are verified by EbsdLib's own test suite; this V&V therefore verifies the **filter's value-add only** — that the `(inputType, outputType)` `switch` dispatches to the correct conversion, that components are read/written with the correct per-tuple stride, and that the preflight contract holds — using EbsdLib-3.0.0-derived values as **dispatch landmarks** plus Class 4 round-trip/`isValid` invariants. The refreshed legacy comparison withdraws one provenance-invalid precision claim (D1) and patch-proves a genuine legacy input-mutation bug (D5); four reportable differences remain.

## Algorithm Relationship

*Classification:* **Rewrite** (of the filter plumbing) under the retained SIMPL UUID `e5629880-98c4-5656-82b8-c9fe2b9744de`.

*Evidence:* The SIMPLNX algorithm `Algorithms/ConvertOrientations.cpp` is structurally distinct from legacy `ConvertOrientations::execute()` / `generateRepresentation<T>()` (DREAM3D 6.5.171):

- **Legacy** built a `QVector` of **7** `OrientationConverter<T>` subclass instances (Euler / OM / Quaternion / AxisAngle / Rodrigues / Homochoric / Cubochoric — **no Stereographic**), called `setInputData()` on the input-type converter, then `convertRepresentationTo(outputType)`, which dispatches each requested pair to its **direct** pairwise transform (`OrientationConverter.hpp:492` `eu2om`, `:517` `eu2cu`→`eu2ho→ho2cu`) — **not** through a quaternion intermediate. Every `toX()` also ran `sanityCheckInputData()`, which for Euler input rewrote the input array in place (D5). Supported **both `float` and `double`** input arrays.
- **SIMPLNX** uses an outer `if`-chain on `OutputType` (8 cases) wrapping an inner `switch` on `InputType` (8 cases) that dispatches to a macro-generated `TO_REP##Convertor` functor calling `inputInstance.to##TO_REP()` (the EbsdLib 2.0 `Orientation` member methods, some direct, some via OM/quaternion). Supports **8** types (adds **Stereographic**) and **float32 only**.
- Per V&V policy, **a Rewrite under the same UUID is a claim of functional equivalence** for the 7 shared types — the Deviations file (Step 8) must defend it.

*Port-time deltas / material changes:*

1. **Dispatch rewrite** — legacy converter-class hierarchy → direct `input.toX()` 8×8 switch. Both dispatch to the same direct pairwise transforms. The originally reported D1 precision residual was withdrawn after its archived "legacy" artifact proved to be the local legacy proof build rather than official 6.5.171; the refreshed official matrix is bit-identical away from singular boundaries.
2. **Stereographic added** (SIMPLNX type 7) — no legacy equivalent; out of scope for legacy A/B.
3. **float32 only** (SIMPLNX) vs **float + double** (legacy) — legacy `double` arrays carried out the math in double precision; SIMPLNX always float32. Deviation candidate for any pipeline that fed `double` arrays to legacy.

*Material PRs since baseline:* #1468 ("ConvertOrientationsFilter uses an Algorithm Class"), #1301 ("Add missing algorithm classes"), #1472 ("Update to EbsdLib 2.0.0 API"), #1535 ("Remove redundant preflight checks"). #1472 is the one that swapped the conversion API to EbsdLib 2.0 `input.toX()`.

## Oracle

*Class:* **3 (Paper-based — Rowenhorst 2015)** primary, **1 (Analytical)** for Stereographic, **4 (Invariant)** companion.

*Citation:* D. Rowenhorst, A. D. Rollett, G. S. Rohrer, M. Groeber, M. Jackson, P. J. Konijnenberg, M. De Graef, "Consistent representations of and conversions between 3D rotations," *Modelling and Simulation in Materials Science and Engineering* **23**(8) 083501 (2015), DOI 10.1088/0965-0393/23/8/083501 — cited throughout `EbsdLib/Source/EbsdLib/Core/OrientationTransformation.hpp`.

*Applied:* The transform equations are implemented and verified **inside EbsdLib** (`EbsdLib/Source/Test/OrientationTest.cpp` exercises the full 8×8 conversion matrix incl. Stereographic with round-trip + `isValid()`; `OrientationTransformationTest.cpp` pins analytical landmarks — identity, `ax2om` 90°-about-Z; `OrientationConverterTest.cpp` pins a Rowenhorst-style Euler→Quaternion exemplar). This filter test does **not** re-verify that math. Instead it takes **one general orientation expressed in all 8 representations** — values generated directly from the same EbsdLib 3.0.0 the filter links (reference implementation, independent of the filter's parallel-convertor plumbing) — and uses them as **dispatch landmarks**: for every `(inputType, outputType)` pair the filter must transform `R[inputType]` into `R[outputType]` within tolerance. Wired to the wrong conversion, the output would be a detectably different number. Multi-tuple input arrays additionally pin the per-tuple component striding. Stereographic specifically is cross-checked against its closed form (`st = (qₓ,qᵧ,q_z)/(1+q_w)`; inverse `ω = 4·atan(|st|), n̂ = st/|st|`) — Class 1. Class 4 round-trip (`A→B→A` ≈ identity) and `isValid()` predicates cover the full matrix cheaply.

*Encoded:* `test/ConvertOrientationsTest.cpp`:
- `"Dispatch and striding (8x8 matrix)"` — 56 `DYNAMIC_SECTION`s (every `(in,out)` pair incl. Stereographic), 3 distinct general orientations per multi-tuple input array, exact-value comparison vs `k_Ref` (EbsdLib-3.0.0-derived landmarks) at tol 1e-4, plus output component-count/tuple-count striding assertions.
- `"Stereographic closed form (Class 1)"` — Quaternion→Stereographic, expected `st = (x,y,z)/(1+w)` computed in-test from the closed form (no EbsdLib call), tol 1e-5.
- **1032 assertions, all pass** (`NX-Com-Qt69-Vtk96-Rel`).

**Oracle-independence caveat (honest scope):** the `k_Ref` landmarks are generated from EbsdLib 3.0.0 — the same library the filter links — so with respect to the *transform math* the 8×8 dispatch test is a **consistency check against EbsdLib's reference implementation**, not an EbsdLib-independent one. What it independently verifies is the filter's own value-add: dispatch routing and per-tuple striding (a mis-wired switch or stride bug produces a detectably wrong number regardless of the landmark's provenance). Two elements are genuinely EbsdLib-independent: (1) the Stereographic path, checked against its closed form computed in-test with no EbsdLib call (Class 1); (2) seed-0, whose orientation is the Rowenhorst 2015 worked example — its quaternion matches EbsdLib `OrientationConverterTest`'s exemplar `{-0.2919894…, 0.319372, 0.1502762…, 0.8889099…}`, but the ultimate authority for that value is the paper's Table (Class 3), not the EbsdLib fixture. The transform math itself is verified inside EbsdLib's own `OrientationTest.cpp` / `OrientationTransformationTest.cpp` suite, which this V&V relies on rather than duplicates.

*Second-engineer review:* **Jared Duffey — 2026-07-14 (approving reviewer, PR #1648).** Supersedes the 2026-07-16 technical-authority self-sign-off.

## Bugs found and fixed

| Deviation | Defect | Affected released versions | Resolution in this branch |
|-----------|--------|----------------------------|---------------------------|
| `ConvertOrientationsFilter-D5` | Legacy Euler sanitization rewrites the input array and changes valid `Φ = π` and out-of-range orientations. | DREAM.3D 6.5.171 only. DREAM3D-NX was not affected. | SIMPLNX converts a local copy of the supplied tuple without the destructive sanitizer. |

## Code path coverage

*7 of 9 reachable paths exercised; 11 paths enumerated. Rows 7–8 (the 8 same-type dispatch arms and 8 `Type::Unknown` arms) are unreachable through the filter and are excluded from the coverage ratio. The 2 gaps are the `-67003` multi-dimensional-component-shape guard and the per-tuple cancel branch. Path 3 (`-67004` component-count mismatch) is covered by the `Invalid preflight` test.*

Source: `src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/Algorithms/ConvertOrientations.cpp` (410 lines) + `Filters/ConvertOrientationsFilter.cpp` preflight. Logical phases: (a) filter `preflightImpl` validation, (b) execute dispatch (8×8 output/input `switch`), (c) per-tuple parallel convertor.

| #  | Phase          | Path   | Test case                                                                 |
|----|----------------|---------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------|
| 1  | (a) Preflight  | `inputType == outputType` → `-67005`                                                              | `Equal Representations` (GENERATE over all 8 types)                       |
| 2  | (a) Preflight  | input component shape has >1 dimension → `-67003`                                                 | *Not directly tested. Low-value guard; selection params produce 1-D component shapes in normal use.* |
| 3  | (a) Preflight  | input component count ≠ expected for input type → `-67004`                                         | `Invalid preflight` — 3-component array declared as Quaternion (expects 4) → `-67004` |
| 4  | (a) Preflight  | out-of-range input/output type index → framework `k_Validate_OutOfRange_Error`                    | `Invalid preflight` (input/output type = 8)                               |
| 5  | (a) Preflight  | valid → `CreateArrayAction` with output type's component count                                    | `Dispatch and striding` (preflight of all 56 pairs) + `Invalid preflight` (does-not-exist) |
| 6  | (b) Dispatch   | 56 cross-type arms (8 outputs × 7 inputs, incl. Stereographic)                                    | `Dispatch and striding (8x8 matrix)` — all 56 `DYNAMIC_SECTION`s; `Stereographic closed form` for qu→st |
| 7  | (b) Dispatch   | 8 same-type arms (`case == output`)                                                               | *Unreachable through the filter — blocked at preflight by path 1.*        |
| 8  | (b) Dispatch   | 8 `case Type::Unknown: break;` arms                                                               | *Unreachable — `ChoicesParameter` range-validates the index (path 4).*    |
| 9  | (c) Convertor  | per-tuple read `inNumComps` → `input.toX()` → write `outNumComps` (striding)                      | `Dispatch and striding` — 3 distinct multi-tuple orientations + output component/tuple-count assertions |
| 10 | (c) Convertor  | `m_Filter->shouldCancel()` → early return                                                         | *Not directly tested. Requires injecting a cancel signal mid-execution; low-value coverage gap.* |
| 11 | (c) Convertor  | `sendThreadSafeProgressMessage()` per chunk (mutex + 1s throttle)                                 | Exercised by every dispatch run (message emitted, not asserted).          |

## Test inventory

| Test case | Status | Notes |
|-----------|--------|-------|
| `OrientationAnalysis::ConvertOrientations: Dispatch and striding (8x8 matrix)` | new-for-V&V | Replaces the retired `Valid filter execution`. 56 `DYNAMIC_SECTION`s (every cross-type pair incl. Stereographic), 3 distinct general orientations per multi-tuple input, exact-value vs EbsdLib-3.0.0 landmarks (tol 1e-4) + output component/tuple-count striding checks. ~1020 assertions. |
| `OrientationAnalysis::ConvertOrientations: Stereographic closed form (Class 1)` | new-for-V&V | Quaternion→Stereographic; expected `st=(x,y,z)/(1+w)` computed in-test (no EbsdLib call), tol 1e-5. Independent analytical pin for the type with no legacy equivalent. |
| `OrientationAnalysis::ConvertOrientations: Invalid preflight` | kept | Negative: does-not-exist input path + out-of-range input/output type index (`ChoicesParameter`). |
| `OrientationAnalysis::ConvertOrientations: Equal Representations` | kept | Negative: same input/output type → `-67005`. `GENERATE` over all 8 types. |
| `OrientationAnalysis::ConvertOrientationsFilter: SIMPL Backwards Compatibility` | kept | `DYNAMIC_SECTION` over SIMPL 6.4 + 6.5 conversion fixtures; validates UUID + argument-key conversion. |
| *(retired)* `OrientationAnalysis::ConvertOrientations: Valid filter execution` | retired | Removed: compared against `k_InitValues` of **unknown provenance** (7×7 only, single tuple, no Stereographic) — a circular-oracle risk. Superseded by the 8×8 dispatch test with EbsdLib-derived, cross-validated landmarks. |

## Exemplar archive

- **None.** This filter's oracle is encoded as **inline dispatch landmarks** (`k_Ref` in `test/ConvertOrientationsTest.cpp`), not a cached `.dream3d`. The landmarks are generated from EbsdLib 3.0.0 (the same library the filter links) and cross-validated independently (seed-0 quaternion == EbsdLib `OrientationConverterTest` exemplar; all stereographic values == closed-form projection). No `download_test_data()` entry and no provenance sidecar are required.

## Deviations from DREAM3D 6.5.171

Four reportable deviations (D2–D5) plus one withdrawn entry (D1) are retained under the plumbing **Rewrite** for traceability. SIMPLNX is independently verified-correct against the oracle. The refreshed comparison covers the original toy, 1,288 canonical-interior tuples, a 1,288-tuple boundary matrix, and an out-of-range adversarial tuple; see `vv/deviations/ConvertOrientationsFilter.md` and the OneDrive working archive.

- `ConvertOrientationsFilter-D1` — **withdrawn:** the archived result attributed to 6.5.171 was produced by local version `1.2.832`, not official 6.5.171 (`1.2.828`). Fresh official comparisons are bit-identical on the original toy and all 1,288 interior tuples for all six shared Euler-input conversions.
- `ConvertOrientationsFilter-D2` — legacy accepted float64 orientation arrays (double-precision math); SIMPLNX is float32-only (precision / scope reduction).
- `ConvertOrientationsFilter-D3` — SIMPLNX adds the Stereographic representation; no legacy equivalent (new capability).
- `ConvertOrientationsFilter-D4` — preflight error-code surface changed (range validation delegated to `ChoicesParameter`); invalid configs still rejected.
- `ConvertOrientationsFilter-D5` — **legacy bug, patch-proven:** in-place Euler "sanitization" mutates out-of-range input and valid float32 `Φ = π` input into a different orientation. A one-line bypass preserves the source tuple and makes corrected legacy agree with NX at float32 precision for the physical orientation.
