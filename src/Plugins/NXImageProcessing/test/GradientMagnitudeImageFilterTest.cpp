#include "ItkGoldenTestUtils.hpp"
#include "ReconstructionFilterTestUtils.hpp"

#include "NXImageProcessing/Filters/GradientMagnitudeImageFilter.hpp"

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/Filter/IFilter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Utilities/ImageProcessing/ImageProcessingConstants.hpp"

#include <catch2/catch.hpp>

#include <string>
#include <vector>

using namespace nx::core;
namespace rt = recon_test;

namespace
{
using GradMagFilterType = GradientMagnitudeImageFilter;

// Sets the geometry, input, and output keys plus the UseImageSpacing flag on the Arguments. The helper runs preflight
// and execute, and requires that both steps succeed.
void RunGradMag(IFilter& filter, DataStructure& dataStructure, const DataPath& inputPath, bool useImageSpacing, unit_test::AlgorithmTestScope* algorithmTestScope = nullptr,
                const std::string& outputName = "Output")
{
  const DataPath geomPath = inputPath.getParent().getParent();
  Arguments args;
  args.insertOrAssign(GradMagFilterType::k_InputImageGeomPath_Key, std::make_any<DataPath>(geomPath));
  args.insertOrAssign(GradMagFilterType::k_InputImageDataPath_Key, std::make_any<DataPath>(inputPath));
  args.insertOrAssign(GradMagFilterType::k_OutputImageArrayName_Key, std::make_any<std::string>(outputName));
  args.insertOrAssign(GradMagFilterType::k_UseImageSpacing_Key, std::make_any<bool>(useImageSpacing));
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeFilter = [&]() { return filter.execute(dataStructure, args); };
  auto executeResult = algorithmTestScope == nullptr ? executeFilter() : algorithmTestScope->execute(executeFilter);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);
}

} // namespace

// -----------------------------------------------------------------------------
// (1) Preflight guards: a multi-component input is rejected (requireScalar); a valid scalar input produces a FIXED
//     float32 output array (AlwaysFloat32Type), regardless of the input element type.
// -----------------------------------------------------------------------------
TEST_CASE("NXImageProcessing::GradientMagnitudeImageFilter: preflight guards", "[NXImageProcessing][GradientMagnitudeImageFilter]")
{
  unit_test::LoadPlugins();
  const unit_test::PreferencesSentinel prefsSentinel(nx::core::DataStorageMode::ForceInCore, 0);

  SECTION("multi-component rejected")
  {
    DataStructure dataStructure;
    const DataPath vecPath = rt::BuildThreeComponentImage(dataStructure);
    const DataPath geomPath = vecPath.getParent().getParent();
    GradMagFilterType filter;
    Arguments args;
    args.insertOrAssign(GradMagFilterType::k_InputImageGeomPath_Key, std::make_any<DataPath>(geomPath));
    args.insertOrAssign(GradMagFilterType::k_InputImageDataPath_Key, std::make_any<DataPath>(vecPath));
    args.insertOrAssign(GradMagFilterType::k_OutputImageArrayName_Key, std::make_any<std::string>("Output"));
    args.insertOrAssign(GradMagFilterType::k_UseImageSpacing_Key, std::make_any<bool>(false));
    const auto result = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_INVALID(result.output_actions);
    REQUIRE(result.output_actions.errors().front().code == nx::core::ImageProcessing::k_NonScalarInput);
  }
  SECTION("scalar input creates a float32 output")
  {
    const auto scenario = GENERATE(from_range(unit_test::SelectAlgorithmTestScenariosForInMemoryStores()));
    CAPTURE(scenario);
    unit_test::AlgorithmTestScope algorithmTestScope(scenario);
    DataStructure dataStructure;
    const std::vector<uint8> field = rt::MakeWrappedRamp<uint8>(6, 6, 1);
    const DataPath inputPath = rt::BuildImageFromPattern<uint8>(dataStructure, 6, 6, 1, field);
    GradMagFilterType filter;
    RunGradMag(filter, dataStructure, inputPath, /*useImageSpacing=*/false, &algorithmTestScope);
    const DataPath outputPath({"Image Geometry", "CellData", "Output"});
    REQUIRE(dataStructure.getDataRefAs<IDataArray>(outputPath).getDataType() == DataType::float32);
  }
  SECTION("non-unit spacing is accepted when UseImageSpacing is off")
  {
    DataStructure dataStructure;
    const std::vector<uint8> field = rt::MakeWrappedRamp<uint8>(6, 6, 6);
    const DataPath inputPath = rt::BuildImageFromPattern<uint8>(dataStructure, 6, 6, 6, field);
    REQUIRE(dataStructure.getDataRefAs<ImageGeom>(inputPath.getParent().getParent()).setSpacing(FloatVec3Type{2.0f, 1.0f, 1.0f}).valid());
    GradMagFilterType filter;
    RunGradMag(filter, dataStructure, inputPath, /*useImageSpacing=*/false);
    const DataPath outputPath({"Image Geometry", "CellData", "Output"});
    REQUIRE(dataStructure.getDataRefAs<IDataArray>(outputPath).getDataType() == DataType::float32);
  }
}

// -----------------------------------------------------------------------------
// (2) FromSIMPLJson: the UseImageSpacing bool and the geometry/array/name DataPaths.
// -----------------------------------------------------------------------------
TEST_CASE("NXImageProcessing::GradientMagnitudeImageFilter: FromSIMPLJson", "[NXImageProcessing][GradientMagnitudeImageFilter]")
{
  const nlohmann::json json = {
      {"UseImageSpacing", true},
      {"SelectedCellArrayPath", {{"Data Container Name", "ImageDC"}, {"Attribute Matrix Name", "CellAM"}, {"Data Array Name", "InArray"}}},
      {"NewCellArrayName", "GradMagOut"},
  };
  const Result<Arguments> result = GradientMagnitudeImageFilter::FromSIMPLJson(json);
  SIMPLNX_RESULT_REQUIRE_VALID(result);
  const Arguments& args = result.value();
  REQUIRE(args.value<bool>(GradMagFilterType::k_UseImageSpacing_Key) == true);
  REQUIRE(args.value<DataPath>(GradMagFilterType::k_InputImageGeomPath_Key) == DataPath({"ImageDC"}));
  REQUIRE(args.value<DataPath>(GradMagFilterType::k_InputImageDataPath_Key) == DataPath({"ImageDC", "CellAM", "InArray"}));
  REQUIRE(args.value<std::string>(GradMagFilterType::k_OutputImageArrayName_Key) == "GradMagOut");
}

// -----------------------------------------------------------------------------
// ITK-sourced real-image golden -- "default" case. Duplicates ITKGradientMagnitudeImageTest.cpp(default) on our
// ITK-free filter: the test reads RA-Float.nrrd (float32) through our reader and runs our filter with the ITK case's
// parameters (all defaults; UseImageSpacing defaults to true). The durable golden compares the output with ITK's
// committed baseline BasicFilters_GradientMagnitudeImageFilter_default.nrrd at the ITK test's tolerance of 1e-05 via
// ip_golden::CompareImages. The output is a fixed float32 (AlwaysFloat32Type). The helper pins ForceInCore.
// -----------------------------------------------------------------------------
TEST_CASE("NXImageProcessing::GradientMagnitudeImageFilter: ITK real-image golden (default)", "[NXImageProcessing][ItkGolden][GradientMagnitudeImageFilter]")
{
  rt::RunDistanceMapItkGoldenBaseline<GradientMagnitudeImageFilter>("RA-Float.nrrd", "BasicFilters_GradientMagnitudeImageFilter_default.nrrd", /*tolerance=*/1e-05);
}
