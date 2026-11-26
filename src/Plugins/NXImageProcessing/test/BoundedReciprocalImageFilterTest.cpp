#include "UnaryFilterTestUtils.hpp"

#include "NXImageProcessing/Filters/BoundedReciprocalImageFilter.hpp"

using namespace nx::core;

// NOTE: this filter has no parity test. The test validates the filter directly against computed-expected
// values: output = 1 / (1 + x), forced to float64.
TEST_CASE("NXImageProcessing::BoundedReciprocalImageFilter: Correctness", "[NXImageProcessing][BoundedReciprocalImageFilter]")
{
  unit_test::LoadPlugins();
  const unit_test::PreferencesSentinel prefsSentinel(nx::core::DataStorageMode::ForceInCore, 0);

  DataStructure dataStructure;
  const DataPath inputPath = ip_test::BuildRampImage<float32>(dataStructure, 12, 0.0, 1.0); // value(i) = i, i in [0, 1727]

  BoundedReciprocalImageFilter filter;
  Arguments args;
  args.insertOrAssign(BoundedReciprocalImageFilter::k_InputImageGeomPath_Key, std::make_any<DataPath>(DataPath({"Image Geometry"})));
  args.insertOrAssign(BoundedReciprocalImageFilter::k_InputImageDataPath_Key, std::make_any<DataPath>(inputPath));
  args.insertOrAssign(BoundedReciprocalImageFilter::k_OutputImageArrayName_Key, std::make_any<std::string>("Output"));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const DataPath outputPath({"Image Geometry", "CellData", "Output"});
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64ArrayType>(outputPath)); // output is forced to float64
  auto& outStore = dataStructure.getDataRefAs<Float64ArrayType>(outputPath).getDataStoreRef();
  for(usize i = 0; i < outStore.getSize(); ++i)
  {
    const double x = static_cast<double>(i);
    const double expected = 1.0 / (1.0 + x);
    REQUIRE(outStore.getValue(i) == Approx(expected).margin(1.0e-9));
  }
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

// -----------------------------------------------------------------------------
// ITK-sourced real-image golden (plan Sec.3/Sec.6). Mirrors ITKBoundedReciprocalImageTest.cpp(defaults): read
// Ramp-Zero-One-Float.nrrd (float32) through OUR reader, run OUR BoundedReciprocalImageFilter (AlwaysFloat64Type output,
// no params), then:
//   (A) DURABLE golden -- compare our float64 output to ITK's committed baseline .nrrd (type: double, read through
//       the SAME reader) at the ITK test's tolerance 0.01 via ip_golden::CompareImages.
// The ITK `vector` case (VM1111Shrink-RGB.png) is SKIPPED as an RGB/vector input (plan Sec.5).
// -----------------------------------------------------------------------------
TEST_CASE("NXImageProcessing::BoundedReciprocalImageFilter: ITK real-image golden (defaults)", "[NXImageProcessing][ItkGolden][BoundedReciprocalImageFilter]")
{
  ip_test::CheckItkRampBaseline(ip_test::RunPointwiseFilter<BoundedReciprocalImageFilter>, "BasicFilters_BoundedReciprocalImageFilter_defaults.nrrd", 0.01);
}
