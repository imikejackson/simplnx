# V&V Report: ITKImageWriterFilter

| | |
|---|---|
| Plugin | ITKImageProcessing |
| SIMPLNX UUID | `a181ee3e-1678-4133-b9c5-a9dd7bfec62f` |
| DREAM3D 6.5.171 equivalent | `ITKImageWriter` (SIMPL UUID `11473711-f94d-5d96-b749-ec36a81ad338`) - `Source/Plugins/ITKImageProcessing/ITKImageProcessingFilters/ITKImageWriter.{h,cpp}` |
| Verified commit | `a307946e7` (v7.4.2 release) |
| Status | COMPLETE |
| Sign-off | Jared Duffey — 2026-07-31 |
| Second-engineer sign-off | Michael A. Jackson <mike.jackson@bluequartz.net> — 2026-08-07 |

## At a glance

| Aspect | Current state |
|---|---|
| Algorithm Relationship | Minor changes - same XY/XZ/YZ pixel extraction; NX uses SIMPLNX stores, AtomicFile, and current ITK APIs. Plane spacing and origin now follow the selected physical axes (D2). |
| Oracle (confirmed) | Classes 1 + 4 - one 3x2x2 scalar fixture (`value(x,y,z)=x+10y+100z`) has same XY/XZ/YZ pixels for all ten accepted scalar types; it verifies selected-plane spacing for all types and origin for MetaImage output. uint8 uses TIFF and the remaining types use MetaImage. |
| Code paths enumerated | 28 of 48 scoped logical paths exercised; 20 documented gaps cover unsupported inputs, defensive branches, progress and cancellation, untested vector component counts, injected write failures, and SIMPL conversion failure. |
| Tests today | 7 named test cases - 1 Class 1+4 Oracle, 2 preflight error-path tests, 1 single-slice exact-name test, 1 RGBA output test, 1 stack-writing test, and 1 SIMPL backwards-compatibility test |
| Exemplar archive | None - Class 1+4 oracle uses inline data |
| Legacy comparison | Reproduced against DREAM3D 6.5.171 with a pipeline-generated copy of the inline oracle and the archived Small IN100 input. All 577 decoded slices match exactly per implementation; D1 and D2 record filename-formatting and plane-metadata differences. See `vv/provenance/ITKImageWriterFilter.md`. |
| Bug flags | Unsigned uint32/uint64 dispatch was corrected. Current changes add invalid fill-character validation, component-count validation, RGBA dispatch, selected-plane physical metadata, and tuple-copy error propagation. |
| V&V phase | **COMPLETE.** |

## Summary

ITKImageWriterFilter exports ImageGeom cell data as an ITK image or a 2D image stack. A hand-derived non-square fixture verifies decoded pixel orientation and plane metadata independently of legacy, while the legacy comparison checks migration output on the toy and production fixtures.

## Algorithm Relationship

**Minor changes.** Port-time changes are SIMPLNX DataStructure stores, AtomicFile writes, and current ITK APIs. Scalar XY/XZ/YZ decoded pixels match the legacy behavior. The NX writer deliberately preserves the selected physical axes in XZ and YZ output spacing/origin; legacy output used identity 2D metadata (D2). The parameter set is unchanged, so `parametersVersion()` remains `2`; validation has been tightened without a parameter-schema change.

## Oracle

**Class:** 1 (Analytical) + 4 (Invariant).

**Applied:** The 3x2x2 scalar fixture has `value(x,y,z)=x+10y+100z`; exact XY, XZ, and YZ pixel matrices are hand-derived and decoded for all ten accepted types. Its non-uniform origin `(10,20,40)` and spacing `(1,2,4)` verify that each written plane has the correct two physical axes: spacing is decoded for every type and origin for MetaImage output. TIFF does not preserve ITK origin metadata. File creation and decoded dimensions are companion invariants.

**Encoded:** `test/ITKImageWriterTest.cpp::ITKImageProcessing::ITKImageWriterFilter: Analytical Pixel Order` - templated over all ten accepted scalar types.

**Second-engineer review:** **Signed off by Michael A. Jackson, 2026-08-07.** The V&V work was authored by Jared Duffey (PR #1693), so the second-engineer review is independent of the author.

*Second-engineer review:* **Michael A. Jackson <mike.jackson@bluequartz.net> — 2026-08-07** (reviewer of the delivering PR; the PR reviewer is the second engineer under project policy).

## Code path coverage

28 of 48 scoped logical paths exercised. Source: `src/Plugins/ITKImageProcessing/src/ITKImageProcessing/Filters/ITKImageWriterFilter.cpp` (518 lines).

Scope: one row represents one filter-controlled logical behavior or failure exit. Repeated equivalent checks in the XY, XZ, and YZ loops are grouped. Element-type and component-count dispatch are counted as orthogonal families rather than as a Cartesian product. Internal implementation branches of parameter classes, dispatch utilities, AtomicFile, and ITK are excluded; configured parameter outcomes and success or failure results handled at this filter's call boundaries are included.

| # | Phase | Path | Test case |
|---|---|---|---|
| 1 | Preflight | Valid dimensions, fill character, in-memory storage, and supported selections reach the preview calculation. | `Analytical Pixel Order`, `Single Slice Keeps Exact Output Name`, `RGBA Image Output`, and `Write Stack`. |
| 2 | Preflight | Dimension mismatch returns `-25600`. | `Dimension Mismatch Validation` uses a 1x1x1 array with a 1x1x2 geometry and asserts `-25600`. |
| 3 | Preflight | Fill-character length other than one returns `-25601`. | `Fill Character Validation` uses an empty string and asserts `-25601`. |
| 4 | Preflight | A one-character format, path, or ASCII control character returns `-25602`. | `Fill Character Validation` asserts `-25602` for `{`, `/`, and newline. |
| 5 | Preflight | An input array not held in memory returns `ITK::Constants::k_OutOfCoreDataNotSupported`. | *Not directly tested; ITK tests run only with in-memory arrays.* |
| 6 | Parameter validation | Unsupported element type is rejected by `ArraySelectionParameter`. | *Not directly tested; the selected array uses one of the ten accepted numeric element types in every execution test.* |
| 7 | Parameter validation | Unsupported component shape is rejected by `ArraySelectionParameter`. | *Not directly tested; current tests use component shapes `{1}` and `{4}`.* |
| 8 | Preflight | XY plane selects the Z slice count for the example output file. | `Write Stack` asserts the XY preview starts at the configured index offset. |
| 9 | Preflight | XZ plane selects the Y slice count for the example output file. | `Single Slice Keeps Exact Output Name` asserts the one-slice XZ preview. |
| 10 | Preflight | YZ plane selects the X slice count for the example output file. | *Not directly assertion-covered; YZ execution is tested, but its preview value is not asserted.* |
| 11 | Preflight | Defensive plane-switch default leaves `maxSlice == 1`. | *Not directly tested; `ChoicesParameter` rejects values outside XY, XZ, and YZ before `preflightImpl`.* |
| 12 | Preflight | Multi-slice preview includes the formatted index suffix and offset. | `Write Stack` asserts `slice_100.tif`. |
| 13 | Preflight | Single-slice preview suppresses the index suffix. | `Single Slice Keeps Exact Output Name` asserts the exact unsuffixed output path. |
| 14 | Execute | Runtime store creation and copy dispatch cover int8, uint8, int16, uint16, int32, uint32, int64, uint64, float32, and float64. | `Analytical Pixel Order` is templated over all ten element types and asserts decoded pixels. |
| 15 | Execute | XY extraction and selected-axis physical metadata. | `Analytical Pixel Order` asserts hand-derived pixels, spacing `(1,2)`, and MetaImage origin `(10,20)`. |
| 16 | Execute | XZ extraction and selected-axis physical metadata. | `Analytical Pixel Order` asserts hand-derived pixels, spacing `(1,4)`, and MetaImage origin `(10,40)`. |
| 17 | Execute | YZ extraction and selected-axis physical metadata. | `Analytical Pixel Order` asserts hand-derived pixels, spacing `(2,4)`, and MetaImage origin `(20,40)`. |
| 18 | Execute | Cancellation between slices returns before copying or writing the next slice. | *Not directly tested; requires cancel-signal injection during execution.* |
| 19 | Execute | Tuple copy succeeds and the next slice operation proceeds. | `Analytical Pixel Order` asserts every copied pixel for all planes and element types. |
| 20 | Execute | Tuple-copy failure propagates from `copyFrom`. | *Not directly tested. This defensive path is excluded by construction: `sliceData` uses the input data type and component shape and exactly the selected plane's tuple shape; preflight requires the input tuple dimensions to match the selected geometry, so computed source and destination ranges are valid.* |
| 21 | Execute | A per-slice message reports the one-based file number, total file count, and output path. | *Not directly tested; current tests do not capture the message handler.* |
| 22 | Filename | Multi-slice execution adds the configured padding, fill character, and index offset. | `Write Stack` asserts the expected offset filenames and that no extra files remain. |
| 23 | Filename | Single-slice execution keeps the exact output filename. | `Single Slice Keeps Exact Output Name` asserts the unsuffixed file exists and decodes correctly. |
| 24 | Filesystem | Output parent directory already exists. | Multi-slice tests assert files after the first slice, which reuse the directory created for that stack. |
| 25 | Filesystem | Missing output parent directory is created successfully. | Successful output tests use new random nested directories and assert the first output file. |
| 26 | Filesystem | Output parent directory creation fails and returns `-19000`. | *Not directly tested; requires filesystem failure injection.* |
| 27 | Component dispatch | Component shape `{1}` uses the scalar path. | `Analytical Pixel Order` covers `{1}` for all ten accepted element types. |
| 28 | Component dispatch | Component shape `{2}` uses the vector path. | *Not directly tested.* |
| 29 | Component dispatch | Component shape `{3}` uses the vector path. | *Not directly tested.* |
| 30 | Component dispatch | Component shape `{4}` uses the RGBA path. | `RGBA Image Output` decodes one uint8 pixel as `(10,20,30,40)`; other element types are not cross-product tested. |
| 31 | Component dispatch | Component shape `{10}` uses the vector path. | *Not directly tested.* |
| 32 | Component dispatch | Component shape `{11}` uses the vector path. | *Not directly tested.* |
| 33 | Component dispatch | Component shape `{36}` uses the vector path. | *Not directly tested.* |
| 34 | Component dispatch | Defensive `ArraySwitchFunc` failure returns `-21010`. | *Not directly tested; allowed element types and component shapes are constrained by the array-selection parameter.* |
| 35 | Atomic write | `AtomicFile::Create` succeeds and provides a temporary path. | Successful output tests assert the committed output files. |
| 36 | Atomic write | `AtomicFile::Create` failure propagates. | *Not directly tested; requires failure injection.* |
| 37 | Atomic write | ITK writer succeeds. | Successful output tests decode the written images and assert their contents. |
| 38 | Atomic write | ITK throws and the filter returns `-21011`. | *Not directly tested; requires ITK writer failure injection.* |
| 39 | Atomic write | Atomic commit succeeds. | Successful output tests assert the final output paths exist and decode correctly. |
| 40 | Atomic write | Atomic commit failure propagates. | *Not directly tested; requires commit failure injection.* |
| 41 | Execute | Successful `SaveImageData` result continues the plane loop. | `Analytical Pixel Order` and `Write Stack` assert that every expected slice is written. |
| 42 | Execute | Failed `SaveImageData` result propagates from the plane loop. | *Not directly tested; covered failure sources require injection.* |

## Test inventory

| Test case | Status | Notes |
|---|---|---|
| `ITKImageProcessing::ITKImageWriterFilter: Analytical Pixel Order` | new-for-V&V | Tests Class 1+4 Oracle over all accepted scalar types, including decoded XY/XZ/YZ pixels, spacing, and MetaImage origin. |
| `ITKImageProcessing::ITKImageWriterFilter: Fill Character Validation` | new-for-V&V | Empty fill character is rejected with `-25601`; format-control and path-separator characters are rejected with `-25602`. |
| `ITKImageProcessing::ITKImageWriterFilter: Dimension Mismatch Validation` | new-for-V&V | Mismatched array and geometry input is rejected during preflight with error `-25600`. |
| `ITKImageProcessing::ITKImageWriterFilter: Single Slice Keeps Exact Output Name` | new-for-V&V | A 3x1x2 ImageGeom XZ plane previews and writes an unsuffixed `.mha` file with exact decoded pixels. |
| `ITKImageProcessing::ITKImageWriterFilter: RGBA Image Output` | new-for-V&V | A uint8 RGBA pixel is dispatched and decodes as `(10,20,30,40)`. |
| `ITKImageProcessing::ITKImageWriterFilter: Write Stack` | kept | Checks that the preflight preview starts at the configured index offset and that the image is written as a stack of files. |
| `ITKImageProcessing::ITKImageWriterFilter: SIMPL Backwards Compatibility` | kept | Covers SIMPL json backwards compatibility. |

## Exemplar archive

No new exemplar archive was created for this V&V cycle: the Class 1 oracle is encoded entirely as inline expected values in the test source.

The oracle and DREAM3D 6.5.171 comparison was run on 2026-07-30. Its provenance includes the legacy setup pipeline, paired writer pipelines, runner hashes, input hashes, zero-tolerance comparison method, and machine-readable results. The reproducible comparison artifacts were uploaded to OneDrive on 2026-07-31.

## Deviations from DREAM3D 6.5.171

- `ITKImageWriterFilter-D1` - NX defaults to zero-padded slice indices while 6.5.171 does not - see `vv/deviations/ITKImageWriterFilter.md`.
- `ITKImageWriterFilter-D2` - NX writes the selected plane's physical spacing and origin, while 6.5.171 used identity 2D metadata for XY/XZ/YZ output - see `vv/deviations/ITKImageWriterFilter.md`.

All pixels otherwise match exactly with zero tolerance: the analytical fixture has 2 XY, 2 XZ, and 3 YZ slices for each of ten scalar types; Small IN100 has 117 XY, 201 XZ, and 189 YZ float32 TIFF slices.
