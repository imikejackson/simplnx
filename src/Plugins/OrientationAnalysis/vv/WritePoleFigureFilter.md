# V&V Report: WritePoleFigureFilter

|           |                          |
|-----------|--------------------------|
| Plugin    | OrientationAnalysis      |
| SIMPLNX UUID | 00cbb97e-a5c2-43e6-9a35-17a0f9ce26ed |
| DREAM3D 6.5.171 equivalent | WritePoleFigure (legacy SIMPL UUID `a10bb78e-fcff-553d-97d6-830a43c85385`) |
| Verified commit | `a307946e7` (v7.4.2 release) |
| Status | COMPLETE |
| Sign-off | Michael Jackson <mike.jackson@bluequartz.net> — 2026-07-16 |
| Second-engineer sign-off | Jared Duffey — 2026-07-10 (approving reviewer, PR #1647). Supersedes the 2026-07-16 technical-authority self-sign-off. |

## At a glance

| Aspect                 | Current state            |
|------------------------|--------------------------|
| Algorithm Relationship | **Rewrite** with the same pole-figure intent; SIMPLNX replaces per-phase PDF output and libharu with image geometry, optional intensity arrays, and EbsdLib raster rendering. |
| Oracle (confirmed)     | **Class 5** expert review covers hexagonal and cubic renders; **Class 4** tests mask, convention, discrete-mode, and marker-radius wiring. All 4 tests pass. |
| Code paths enumerated  | 10 of 13 wrapper paths exercised; 3 defensive or validation paths remain uncovered. EbsdLib owns per-Laue-class rendering tests. |
| Tests today            | 4 test cases cover mask behavior, hex convention, discrete markers, and SIMPL conversion. |
| Exemplar archive       | `Pole_Figure_Exemplars_v6.tar.gz` contains input orientations and a mask only; no rendered image is used as an oracle. |
| Legacy comparison      | **Run** — expert review found visually identical pole data across official DREAM3D 6.5.171, the local legacy proof build, and SIMPLNX; D1–D5 are non-defect presentation or API differences. |
| Bug flags              | None; all 5 deviations are cosmetic, intentional rendering changes, or output-format differences. |
| V&V phase | **COMPLETE.** |

## Summary

`WritePoleFigureFilter` generates `<001>/<011>/<111>` (or the hexagonal/trigonal equivalents) pole figures for each phase from per-cell Euler angles, phases, and crystal structures, optionally masked. It was verified by expert (Class 5) side-by-side comparison of hex and cubic pole figures rendered through official DREAM3D 6.5.171, a locally patched legacy build, and SIMPLNX (EbsdLib 3.1.0), backed by Class 4 invariant unit tests for the simplnx-unique wiring. The pole-figure data is visually identical across all versions; the differences are cosmetic labeling/font, an intentional discrete-marker rendering improvement, and PNG-only output — five documented, non-defect deviations.

## Algorithm Relationship

**Rewrite**

*Evidence:* SIMPLNX inherits the legacy SIMPL UUID `a10bb78e-…`, but the implementation and output are substantially different:

1. **Output medium** — legacy writes one **PDF per phase** to disk (libharu/HPDF); nothing enters the data structure. SIMPLNX creates an **ImageGeometry** with an RGB `Phase_N` array (and optional Float64 intensity arrays) in the DataStructure, and optionally writes a raster image (TIFF/PNG) to disk.
2. **Rendering stack** — legacy draws pole-figure chrome (circle, axes, labels, color bar) directly with libharu primitives; SIMPLNX delegates to the **EbsdLib pole-figure compositor** (`GeneratePoleFigureComposite`), whose rasterized output is byte-tested upstream in `PoleFigureCompositorTest::All_Laue_Classes`.
3. **Discrete mode** — SIMPLNX renders discrete figures with EbsdLib 3.1.0's **vector-marker renderer** (configurable **Discrete Marker Radius**), replacing legacy's single-black-pixel-per-orientation.
4. **New capabilities** — optional intensity Float64 arrays, MRD normalization, in-DataStructure image geometry, and the X‖a / X‖a* Hex/Trig basis convention parameter.

The shared pole-figure projection math (modified Lambert for Color, stereographic for Discrete) descends from the same lineage, which is why the rendered figures agree visually. Being a Rewrite under the same UUID, the Deviations file defends the equivalence claim.

## Oracle

*Class:* **5 (Expert-visual)** primary, **4 (Invariant)** companion.

*Justification for Class 5:* the filter's output is a rasterized pole-figure image whose correctness is inherently visual; legacy DREAM3D emits **only PDFs** (no numeric ground truth to diff), and the pixel-level rendering is owned and byte-tested by EbsdLib upstream. No Class 1–3 oracle fully specifies the rendered image. The analytically-checkable part (pole positions) and the simplnx-unique wiring are covered by the Class 4 invariants below.

*Applied:* the same 502 hex-Ti orientations (and a cubic-relabeled variant) were rendered through official DREAM3D 6.5.171, a locally patched legacy build, and SIMPLNX (EbsdLib 3.1.0) in Color and Discrete modes; a domain expert reviewed the figures side by side. Pole positions, intensity distribution, and color-intensity mapping are visually identical across all three; the differences are the five cosmetic/rendering/API items in the Deviations file.

*Encoded:*
- **Class 4 (Invariant):** `test/WritePoleFigureTest.cpp::"OrientationAnalysis::WritePoleFigureFilter: Mask filter changes the rendered pole figure"` (masked output differs from unmasked by >1% of bytes → mask is wired) and `::"…: HexConvention choice reaches algorithm"` (X‖a vs X‖a* rotates the basal families 30° in both the intensity array and the composite RGB → both plumbing paths honor the convention). Both pass against EbsdLib 3.1.0.
- **Class 5 (Expert-visual):** the hex + cubic renders (legacy PDFs + SIMPLNX PNGs), signed off. Generator scripts + pipelines are committed under `Code_Review/vv/WritePoleFigure/`; the binary renders are archived to OneDrive — see the provenance sidecar and that folder's `README.md`.
- **Class 2 (Reference), cited not duplicated:** EbsdLib `PoleFigureCompositorTest::All_Laue_Classes` pins per-Laue-class pixel reproduction.

*Second-engineer review:* **Jared Duffey — 2026-07-10 (approving reviewer, PR #1647).** Supersedes the 2026-07-16 technical-authority self-sign-off. Review focus: the five documented differences are confirmed as the complete set and are all non-defects.

## Bugs found and fixed

None.

## Code path coverage

10 of 13 simplnx-wrapper paths exercised in CI. The filter is a wrapper around the EbsdLib compositor; per-Laue-class projection and pixel rendering are owned/tested by EbsdLib. Logical phases: (a) preflight validation + array creation, (b) parameter→enum translation, (c) per-phase mask filtering, (d) intensity generation, (e) composite image generation.

Source: `src/Plugins/OrientationAnalysis/src/OrientationAnalysis/Filters/Algorithms/WritePoleFigure.cpp` (787 lines) + `WritePoleFigureFilter.cpp`.

| #  | Phase              | Path                                                                     | Test case                                                        |
|----|--------------------|--------------------------------------------------------------------------|------------------------------------------------------------------|
| 1  | (b) Param→enum     | GenerationAlgorithm = Color (0)                                          | Mask + HexConvention tests; legacy A/B (hex + cubic)             |
| 2  | (b) Param→enum     | GenerationAlgorithm = Discrete (1) → vector markers                     | `Discrete mode and marker radius reach algorithm` (CI) + legacy A/B Discrete renders (hex + cubic), expert sign-off |
| 3  | (b) Param→enum     | HexConvention X‖a (0) vs X‖a* (1)                                        | `HexConvention choice reaches algorithm` (intensity + composite) |
| 4  | (b) Param→enum     | ImageLayout (Horizontal/Vertical/Square)                                | *Not directly asserted; Horizontal exercised in all runs. EbsdLib layout enum — low-value to sweep here.* |
| 5  | (c) Mask filter    | UseMask off / on                                                        | `Mask filter changes the rendered pole figure`                   |
| 6  | (c) Per-phase      | per-phase Euler extraction; empty phase → skip                          | Single-phase fixtures exercise extraction. *Empty-phase skip not directly tested — defensive guard.* |
| 7  | (d) Intensity      | SaveIntensityDataArrays on + NormalizeToMRD                             | `HexConvention choice reaches algorithm` (SaveIntensity=true, MRD=true) |
| 8  | (d) Intensity      | crystal-structure dispatch (Cubic_High, Hexagonal_High)                 | Hex: tests + A/B; Cubic: A/B. Other Laue classes owned by EbsdLib upstream. |
| 9  | (d) Intensity      | unknown crystal structure → warning, skip                              | *Not directly tested — defensive warning branch.*                |
| 10 | (e) Composite      | SaveAsImageGeometry / WriteImageToDisk; DiscreteMarkerRadius            | `Discrete mode and marker radius reach algorithm` (1 px vs 10 px composites differ); A/B renders (write-to-disk); Mask/HexConv tests (image geometry) |
| 11 | (a) Preflight      | mask array wrong type → error `-53900`                                  | *Not directly tested — low-value validation branch.*             |
| 12 | (a) Preflight      | ImageSize ≤ 0 → error `-680002`; Discrete mode with marker radius < 1 → error `-680003` | *Guards added during review; not directly tested — low-value validation branches.* |
| 13 | (h) SIMPL convert  | 6.4 / 6.5 SIMPL JSON → Arguments (legacy ImageFormat dropped, D5)       | `SIMPL Backwards Compatibility`                                  |

## Test inventory

| Test case | Status | Notes |
|-----------|--------|-------|
| `Mask filter changes the rendered pole figure` | kept | Class 4 invariant: masked vs unmasked composite RGB differ by >1% of bytes. Consumes `Pole_Figure_Exemplars_v6`. Passes on EbsdLib 3.1.0. |
| `HexConvention choice reaches algorithm` | kept | Class 4 invariant: X‖a vs X‖a* differ in both the intensity array and the composite RGB (both plumbing paths). Passes on EbsdLib 3.1.0. |
| `Discrete mode and marker radius reach algorithm` | new-for-V&V | Class 4 invariant: Discrete vs Color composites differ; 1 px vs 10 px marker radii differ. Pins the GenerationAlgorithm and DiscreteMarkerRadius plumbing in CI (previously covered only by manual A/B renders). |
| `SIMPL Backwards Compatibility` | kept | 6.4 + 6.5 SIMPL JSON → Arguments round-trip. The legacy `ImageFormat` key is intentionally not converted (D5). |

## Exemplar archive

- **Archive:** `Pole_Figure_Exemplars_v6.tar.gz` (inputs only — 502 hex-Ti orientations + 251/251 mask + crystal structures + phase names).
- **SHA512:** *(see `test/CMakeLists.txt` `download_test_data(... Pole_Figure_Exemplars_v6.tar.gz ...)`)*
- **Provenance:** `src/Plugins/OrientationAnalysis/vv/provenance/WritePoleFigureFilter.md`
- No baked image exemplar: pixel-level reproduction is owned by EbsdLib upstream; duplicating it here couples simplnx CI to EbsdLib byte-identity (the source of prior v5 baseline drift).

## Deviations from DREAM3D 6.5.171

Established by expert (Class 5) visual comparison on hex and cubic; official 6.5.171 and the locally patched legacy build are visually indistinguishable (three raster-identical cases; one 7-pixel residual). Full renders are archived to OneDrive. All five are cosmetic, labeling, intentional-rendering, or output-format/API differences; none is a defect.

- `WritePoleFigureFilter-D1` — axis labels `X`/`Y` (legacy) → `A1`/`A2` (SIMPLNX). Cosmetic.
- `WritePoleFigureFilter-D2` — font/text-metrics differ (libharu Helvetica → EbsdLib canvas_ity). Library, cosmetic.
- `WritePoleFigureFilter-D3` — hex/trig per-figure family labels differ but name symmetry-equivalent families. Cosmetic.
- `WritePoleFigureFilter-D4` — Discrete markers: filled circles (configurable radius) → intentional rendering improvement over legacy single black pixels.
- `WritePoleFigureFilter-D5` — legacy Image Format choice (tif/bmp/png/pdf) dropped; SIMPLNX always writes PNG. The legacy `ImageFormat` key is not converted from SIMPL JSON.

See `vv/deviations/WritePoleFigureFilter.md`.
