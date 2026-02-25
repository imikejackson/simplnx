#include <catch2/catch.hpp>

#include "SimplnxCore/Filters/CreateScatterChartFilter.hpp"
#include "SimplnxCore/SimplnxCore_test_dirs.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"
#include "simplnx/Parameters/VectorParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

using namespace nx::core;

TEST_CASE("SimplnxCore::CreateScatterChartFilter: Valid Filter Execution", "[SimplnxCore][CreateScatterChartFilter]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure;

  // Create input X and Y data arrays
  const usize numPoints = 50;
  DataPath xDataPath({"XData"});
  DataPath yDataPath({"YData"});

  auto* xArray = DataArray<float32>::CreateWithStore<DataStore<float32>>(dataStructure, "XData", {numPoints}, {1});
  auto* yArray = DataArray<float32>::CreateWithStore<DataStore<float32>>(dataStructure, "YData", {numPoints}, {1});

  REQUIRE(xArray != nullptr);
  REQUIRE(yArray != nullptr);

  // Fill with simple data
  for(usize i = 0; i < numPoints; i++)
  {
    (*xArray)[i] = static_cast<float32>(i);
    (*yArray)[i] = static_cast<float32>(i * 2.0f + 1.0f);
  }

  // Set up filter arguments
  Arguments args;
  args.insertOrAssign(CreateScatterChartFilter::k_XDataPath_Key, std::make_any<DataPath>(xDataPath));
  args.insertOrAssign(CreateScatterChartFilter::k_YDataPath_Key, std::make_any<DataPath>(yDataPath));
  args.insertOrAssign(CreateScatterChartFilter::k_ChartDimensions_Key, std::make_any<VectorUInt32Parameter::ValueType>(std::vector<uint32>{400, 600}));
  args.insertOrAssign(CreateScatterChartFilter::k_ChartTitle_Key, std::make_any<StringParameter::ValueType>("Test Scatter Chart"));
  args.insertOrAssign(CreateScatterChartFilter::k_XAxisTitle_Key, std::make_any<StringParameter::ValueType>("X Values"));
  args.insertOrAssign(CreateScatterChartFilter::k_YAxisTitle_Key, std::make_any<StringParameter::ValueType>("Y Values"));
  args.insertOrAssign(CreateScatterChartFilter::k_OutputImageGeometryPath_Key, std::make_any<DataPath>(DataPath({"Chart Image"})));
  args.insertOrAssign(CreateScatterChartFilter::k_CellDataName_Key, std::make_any<DataObjectNameParameter::ValueType>(ImageGeom::k_CellAttributeMatrixName));
  args.insertOrAssign(CreateScatterChartFilter::k_ImageDataName_Key, std::make_any<DataObjectNameParameter::ValueType>("ImageData"));

  CreateScatterChartFilter filter;

  // Preflight the filter
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)

  // Execute the filter
  auto executeResult = filter.execute(dataStructure, args, nullptr, IFilter::MessageHandler{[](const IFilter::Message& message) { fmt::print("{}\n", message.message); }});
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)

  // Verify the output ImageGeometry was created
  DataPath outputGeomPath({"Chart Image"});
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(outputGeomPath));
  const auto& imageGeom = dataStructure.getDataRefAs<ImageGeom>(outputGeomPath);

  // Verify dimensions: X=600 (width), Y=400 (height), Z=1
  SizeVec3 dims = imageGeom.getDimensions();
  REQUIRE(dims[0] == 600);
  REQUIRE(dims[1] == 400);
  REQUIRE(dims[2] == 1);

  // Verify the image data array exists and has correct size
  DataPath imageDataPath = outputGeomPath.createChildPath(ImageGeom::k_CellAttributeMatrixName).createChildPath("ImageData");
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<DataArray<uint8>>(imageDataPath));
  const auto& imageData = dataStructure.getDataRefAs<DataArray<uint8>>(imageDataPath);

  // 600 * 400 pixels, 3 components (RGB) per pixel
  REQUIRE(imageData.getNumberOfTuples() == 600 * 400);
  REQUIRE(imageData.getNumberOfComponents() == 3);

  // Write the DataStructure out to the file system
#ifdef SIMPLNX_WRITE_TEST_OUTPUT
  UnitTest::WriteTestDataStructure(dataStructure, fs::path(fmt::format("{}/CreateScatterChartFilterTest.dream3d", unit_test::k_BinaryTestOutputDir)));
#endif

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::CreateScatterChartFilter: Mismatched Array Sizes", "[SimplnxCore][CreateScatterChartFilter]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure;

  // Create input arrays with different tuple counts
  DataPath xDataPath({"XData"});
  DataPath yDataPath({"YData"});

  auto* xArray = DataArray<float32>::CreateWithStore<DataStore<float32>>(dataStructure, "XData", {10}, {1});
  auto* yArray = DataArray<float32>::CreateWithStore<DataStore<float32>>(dataStructure, "YData", {20}, {1});

  REQUIRE(xArray != nullptr);
  REQUIRE(yArray != nullptr);

  Arguments args;
  args.insertOrAssign(CreateScatterChartFilter::k_XDataPath_Key, std::make_any<DataPath>(xDataPath));
  args.insertOrAssign(CreateScatterChartFilter::k_YDataPath_Key, std::make_any<DataPath>(yDataPath));
  args.insertOrAssign(CreateScatterChartFilter::k_ChartDimensions_Key, std::make_any<VectorUInt32Parameter::ValueType>(std::vector<uint32>{400, 600}));
  args.insertOrAssign(CreateScatterChartFilter::k_ChartTitle_Key, std::make_any<StringParameter::ValueType>("Test"));
  args.insertOrAssign(CreateScatterChartFilter::k_XAxisTitle_Key, std::make_any<StringParameter::ValueType>("X"));
  args.insertOrAssign(CreateScatterChartFilter::k_YAxisTitle_Key, std::make_any<StringParameter::ValueType>("Y"));
  args.insertOrAssign(CreateScatterChartFilter::k_OutputImageGeometryPath_Key, std::make_any<DataPath>(DataPath({"Chart Image"})));
  args.insertOrAssign(CreateScatterChartFilter::k_CellDataName_Key, std::make_any<DataObjectNameParameter::ValueType>(ImageGeom::k_CellAttributeMatrixName));
  args.insertOrAssign(CreateScatterChartFilter::k_ImageDataName_Key, std::make_any<DataObjectNameParameter::ValueType>("ImageData"));

  CreateScatterChartFilter filter;

  // Preflight should fail due to mismatched array sizes
  auto preflightResult = filter.preflight(dataStructure, args);
  REQUIRE(preflightResult.outputActions.invalid());
}
