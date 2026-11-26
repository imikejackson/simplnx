#include "UnaryFilterTestUtils.hpp"

#include "NXImageProcessing/Filters/AcosImageFilter.hpp"

using namespace nx::core;

// PART A domain edge (documents ITK-faithful behavior; does NOT "fix" it): acos is undefined for |x|>1. The op is
// static_cast<U>(std::acos(static_cast<double>(x))) -- exactly ITK's equally-unguarded Acos functor -- so acos(2)
// is IEEE NaN. This pins that observed output.
TEST_CASE("NXImageProcessing::AcosImageFilter: acos of |x|>1 is NaN (domain edge)", "[NXImageProcessing][AcosImageFilter]")
{
  unit_test::LoadPlugins();
  const unit_test::PreferencesSentinel prefsSentinel(nx::core::DataStorageMode::ForceInCore, 0);

  DataStructure dataStructure;
  const DataPath input = ip_test::BuildScalarValuesImage<float32>(dataStructure, {2.0f});
  const DataPath output = ip_test::RunPointwiseFilter<AcosImageFilter>(dataStructure, input);

  const auto& out = dataStructure.getDataRefAs<DataArray<float32>>(output).getDataStoreRef();
  REQUIRE(std::isnan(out.getValue(0)));
}

// -----------------------------------------------------------------------------
// ITK-sourced real-image golden (plan Sec.3/Sec.6). Mirrors ITKAcosImageTest.cpp(defaults): read
// Ramp-Zero-One-Float.nrrd (float32) through OUR ITK-free reader, run OUR AcosImageFilter with the ITK case's
// params (defaults), then compare our output to ITK's committed baseline .nrrd (read through the SAME reader) at
// the ITK test's tolerance 0.01 via ip_golden::CompareImages. CompareImages dispatches on the runtime type and --
// unlike CompareDataArrays -- tolerates the non-finite tails these functors can emit.
// -----------------------------------------------------------------------------
TEST_CASE("NXImageProcessing::AcosImageFilter: ITK real-image golden (defaults)", "[NXImageProcessing][ItkGolden][AcosImageFilter]")
{
  ip_test::CheckItkRampBaseline(ip_test::RunPointwiseFilter<AcosImageFilter>, "BasicFilters_AcosImageFilter_defaults.nrrd", 0.01);
}
