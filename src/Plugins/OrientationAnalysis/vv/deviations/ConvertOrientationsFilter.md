# Deviations from DREAM3D 6.5.171: ConvertOrientationsFilter

This file lists every documented behavioral difference between this SIMPLNX filter and its DREAM3D 6.5.171 equivalent (`ConvertOrientations`, SIMPL UUID `e5629880-98c4-5656-82b8-c9fe2b9744de`).

Entries are referenced by stable ID (`ConvertOrientationsFilter-D<N>`) from the V&V report and from public migration guidance. The ID is stable across renames; the Filter UUID field is the permanent cross-reference anchor.

The SIMPLNX algorithm is a **Rewrite** of the filter plumbing under the retained UUID (see `../ConvertOrientationsFilter.md`). Four reportable deviations remain below; D1 is retained as a withdrawn entry because the original comparison artifact was misidentified. SIMPLNX is verified-correct independently of 6.5.171 against the Class 3 (Rowenhorst 2015) / Class 1 / Class 4 oracle encoded in `test/ConvertOrientationsTest.cpp` (1032 assertions).

> **Comparison status:** The comparison was refreshed 2026-09-17 using official DREAM3D 6.5.171 (`1.2.828.f45085c83`), NX 7.4.1/EbsdLib 2.2.0, NX 7.4.2/EbsdLib 3.1.2, and a surgically corrected local legacy build. The matrix covers 1,288 canonical-interior tuples, 1,288 tuples including the `Φ = π` boundary, and an out-of-range adversarial tuple. D2–D4 are additionally pinned by source/API inspection.

---

## ConvertOrientationsFilter-D1

| Field | Value |
|---|---|
| **Deviation ID** | `ConvertOrientationsFilter-D1` |
| **Filter UUID** | `501e54e6-a66f-4eeb-ae37-00e649c00d4b` (SIMPL `e5629880-98c4-5656-82b8-c9fe2b9744de`) |
| **Status** | retired 2026-09-17 — comparison artifact was not DREAM3D 6.5.171 |

**Symptom:** Withdrawn — the archived toy comparison reported four of six Euler-input conversions bit-identical and residuals of 1.5e-8 (orientation matrix) and 1.78e-6 (cubochoric), attributed to official DREAM3D 6.5.171.

**Root cause:** Not applicable — retired. The archived artifact came from the local proof build, not the official DREAM3D 6.5.171 release. Its root HDF5 producer string is DREAM3D `1.2.832.f70859912`; the official application reports `1.2.828.f45085c83`. Re-running the embedded pipeline through the official runner made all six Euler-input outputs bit-identical to NX 7.4.2. A broader 1,288-tuple interior matrix was then bit-identical for all components of all six outputs against both NX 7.4.1/EbsdLib 2.2.0 and NX 7.4.2/EbsdLib 3.1.2. The old residual is not evidence of a deviation from the required baseline.

**Affected users:** None. The recorded residual belonged to a mislabeled local-build artifact.

**Recommendation:** Do not count this entry as a deviation. It is retained only to preserve the audit trail.

---

## ConvertOrientationsFilter-D2

| Field | Value |
|---|---|
| **Deviation ID** | `ConvertOrientationsFilter-D2` |
| **Filter UUID** | `501e54e6-a66f-4eeb-ae37-00e649c00d4b` (SIMPL `e5629880-98c4-5656-82b8-c9fe2b9744de`) |
| **Status** | active |

**Symptom:** A pipeline that supplied a **`double`** (float64) orientation array to 6.5.171 cannot supply one to SIMPLNX (the input parameter accepts float32 only), and the output dtype is float32 rather than float64.

**Root cause:** Precision (deliberate scope reduction). 6.5.171 `ConvertOrientations::execute()` branched on the input array type and ran the conversion in `double` for `DoubleArrayType` inputs (`generateRepresentation<double>`). SIMPLNX restricts the input `ArraySelectionParameter` to `DataType::float32` and converts in float32 only. For float64 inputs the legacy intermediate math carried ~16 digits vs SIMPLNX's ~7.

**Affected users:** The small number of legacy pipelines that stored orientations as `double`. Standard EBSD ingest produces float32 orientations, which are unaffected.

**Recommendation:** *trust SIMPLNX for float32 workflows.* Users with float64 orientation arrays who require double-precision conversion should note this scope reduction; for EBSD-scale data float32 is the native precision and the difference is immaterial.

---

## ConvertOrientationsFilter-D3

| Field | Value |
|---|---|
| **Deviation ID** | `ConvertOrientationsFilter-D3` |
| **Filter UUID** | `501e54e6-a66f-4eeb-ae37-00e649c00d4b` (SIMPL `e5629880-98c4-5656-82b8-c9fe2b9744de`) |
| **Status** | active |

**Symptom:** SIMPLNX offers a **Stereographic** representation (input/output type index 7) that 6.5.171 does not.

**Root cause:** Algorithmic choice (new capability). 6.5.171 `generateRepresentation<T>` constructed a 7-element converter vector (Euler, OrientationMatrix, Quaternion, AxisAngle, Rodrigues, Homochoric, Cubochoric); Stereographic did not exist. SIMPLNX adds the 8th type, verified analytically (Class 1, `st = (x,y,z)/(1+w)`) in the unit test.

**Affected users:** None negatively — purely additive. There is no 6.5.171 output to compare against for any conversion involving Stereographic.

**Recommendation:** *trust SIMPLNX.* New capability with an independent analytical oracle; no legacy equivalent exists.

---

## ConvertOrientationsFilter-D4

| Field | Value |
|---|---|
| **Deviation ID** | `ConvertOrientationsFilter-D4` |
| **Filter UUID** | `501e54e6-a66f-4eeb-ae37-00e649c00d4b` (SIMPL `e5629880-98c4-5656-82b8-c9fe2b9744de`) |
| **Status** | active |

**Symptom:** Selecting the same representation for input and output produces an error in SIMPLNX (`-67005`) before any computation; 6.5.171 errored similarly (`-1000`) but with a different code/message, and 6.5.171 additionally emitted distinct error codes for out-of-range type indices (`-1001`/`-1002`).

**Root cause:** Algorithmic choice (preflight refactor). SIMPLNX delegates input/output type-range validation to the `ChoicesParameter` (which emits the framework `k_Validate_OutOfRange_Error`) and keeps only the same-type check (`-67005`) and array-shape checks (`-67003`/`-67004`) in `preflightImpl`. The legacy `-1001`/`-1002`/`-1004` (converter-failure) error codes have no SIMPLNX equivalent.

**Affected users:** Scripts or tests that matched on the legacy numeric error codes `-1000`/`-1001`/`-1002`. No effect on successful conversions.

**Recommendation:** *trust SIMPLNX.* Behavior is equivalent (invalid configurations are still rejected at preflight); only the error-code surface changed.

---

## ConvertOrientationsFilter-D5

| Field | Value |
|---|---|
| **Deviation ID** | `ConvertOrientationsFilter-D5` |
| **Filter UUID** | `501e54e6-a66f-4eeb-ae37-00e649c00d4b` (SIMPL `e5629880-98c4-5656-82b8-c9fe2b9744de`) |
| **Status** | active |
| **Bug flag** | **Legacy bug — valid-input corruption, empirically confirmed and patch-proven** |

**Symptom:** DREAM3D 6.5.171 mutates the user's stored Euler array before every conversion and can produce a genuinely different orientation. This occurs for out-of-range input and for the valid float32 endpoint `Φ = π`: because float32 π is slightly larger than the double constant used by `fmod(Φ, π)`, legacy changes it to approximately `8.74e-8`. In a 1,288-tuple canonical matrix, all 143 tuples at `Φ = π` were mutated; SIMPLNX preserved every input tuple.

**Root cause:** Bug in 6.5.171. Every legacy `toX()` ran `sanityCheckInputData()` (`OC_CONVERT_BODY`, `OrientationConverter.hpp:386,403`); for Euler input, `EulerConverter::sanityCheckInputData()` ran `EulerSanityCheck` (`:425-446`) **in place on the actual stored array**: `fmod(φ1, 2π)`, `fmod(Φ, π)`, `fmod(φ2, 2π)` followed by sign flips for negative values. This is not rotation-preserving, and applying modulo π to a closed-range endpoint is wrong even before the destructive side effect is considered. SIMPLNX copies each tuple into a local `OrientationF` and converts exactly the orientation supplied.

**Controlled proof (2026-09-17):** A local legacy proof build was changed by one line: return before the Euler sanitizer. The exact patched commit is recorded in the filter A/B archive. On adversarial input `(-0.75, 4.0, 7.0)`, official legacy stored `(0.75, 0.8584073, 0.7168147)` and differed from NX by 1.09–1.66 across the six conversions; corrected legacy preserved the original tuple and all six outputs agreed with NX within `1.2e-7`. On the canonical matrix, corrected legacy preserved all 1,288 tuples and its orientation matrices agreed with NX within `5.96e-8` (quaternions bit-identical). Raw Rodrigues/cubochoric coordinates at the 180° singular boundary retain representation-convention differences, but the physical orientation is the same; they do not refute the D5 mechanism.

**Affected users:** Pipelines containing valid `Φ = π` values as well as Euler angles outside the fundamental Bunge ranges (for example, negative angles from upstream arithmetic or Φ ∈ (π, 2π)). Under 6.5.171 these inputs were silently rewritten and converted as a different orientation.

**Recommendation:** Trust SIMPLNX. It preserves the input and the valid `Φ = π` endpoint. Users with out-of-range Euler data should normalize it explicitly upstream rather than rely on a lossy implicit rewrite.
