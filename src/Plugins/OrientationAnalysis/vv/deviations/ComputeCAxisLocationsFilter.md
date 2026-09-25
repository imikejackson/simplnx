# Deviations from DREAM3D 6.5.171: ComputeCAxisLocationsFilter

This file lists every documented behavioral difference between this SIMPLNX filter and its DREAM3D 6.5.171 equivalent.

Entries are referenced by stable ID (`ComputeCAxisLocationsFilter-D<N>`) from the V&V report and from public migration guidance. The ID is stable across renames; the Filter UUID field is the permanent cross-reference anchor.

---

## ComputeCAxisLocationsFilter-D1

| Field | Value |
|---|---|
| **Deviation ID** | `ComputeCAxisLocationsFilter-D1` |
| **Filter UUID** | `a51c257a-ddc1-499a-9b21-f2d25a19d098` |
| **Status** | active |

**Symptom:** For non-hexagonal cells, NX places NaN values whereas 6.5.171 places meaningless-but-finite computed values.

**Root cause:** Algorithmic choice. This intentional improvement tells the user that values cannot be calculated for those phases.

**Comparison evidence (2026-09-17):** On a byte-identical 15-cell mixed fixture, all 8 hexagonal cells are bit-identical between official DREAM3D 6.5.171 and NX 7.4.2. Legacy writes finite vectors for the 7 cubic cells; NX writes `(NaN, NaN, NaN)` for all 7, isolating this domain guard without changing valid-domain output.

**Affected users:** Anyone who has non-hexagonal phases in their input to ComputeCAxisLocationsFilter.

**Recommendation:** Trust SIMPLNX - ComputeCAxisLocationsFilter's calculation is only correct for hexagonal cells. The legacy filter computes meaningless values for non-hexagonal phases.

---

## ComputeCAxisLocationsFilter-D2

| Field | Value |
|---|---|
| **Deviation ID** | `ComputeCAxisLocationsFilter-D2` |
| **Filter UUID** | `a51c257a-ddc1-499a-9b21-f2d25a19d098` |
| **Status** | active |

**Symptom:** When there are no hexagonal phases present, NX emits an error (-3522) whereas 6.5.171 executes.

**Root cause:** Algorithmic choice. This intentional improvement prevents execution on data for which the filter is not valid.

**Comparison evidence (2026-09-17):** On a byte-identical all-cubic fixture, official DREAM3D 6.5.171 executes and writes an output, while NX 7.4.2 returns `-3522` and writes no output.

**Affected users:** Anyone who has no hexagonal phases in their input to ComputeCAxisLocationsFilter.

**Recommendation:** Trust SIMPLNX - ComputeCAxisLocationsFilter's calculation is only correct for hexagonal cells. If there are no hexagonal cells, then the filter does no actual calculation.

---

## ComputeCAxisLocationsFilter-D3

| Field | Value |
|---|---|
| **Deviation ID** | `ComputeCAxisLocationsFilter-D3` |
| **Filter UUID** | `a51c257a-ddc1-499a-9b21-f2d25a19d098` |
| **Status** | active |

**Symptom:** NX emits an unconditional warning (-3521) in preflight which advises the user to make sure their data has hexagonal phases and emits a warning (-3523) if there are non-hexagonal phases. 6.5.171 does not emit any warning in either case.

**Root cause:** Algorithmic choice. This intentional improvement warns the user that the filter output is valid only for hexagonal phases.

**Comparison evidence (2026-09-17):** The mixed NX run emits the unconditional `-3521` reminder and the `-3523` non-hex warning. The all-cubic run emits `-3521` before the D2 `-3522` rejection. Legacy emits none of these warnings.

**Affected users:** Anyone running the filter for the preflight warning, and anyone who has mixed non-hexagonal phases in their input to ComputeCAxisLocationsFilter.

**Recommendation:** Trust SIMPLNX - ComputeCAxisLocationsFilter's calculation is only correct for hexagonal cells.

---
