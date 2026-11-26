#include "SimplnxCore/Filters/ChangeAngleRepresentationFilter.hpp"
#include "SimplnxCore/SimplnxCore_test_dirs.hpp"

#include "simplnx/Common/Numbers.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Utilities/AlgorithmDispatch.hpp"
#include "simplnx/Utilities/DataStoreUtilities.hpp"

#include <array>
#include <catch2/catch.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <nonstd/span.hpp>

using namespace nx::core;
namespace fs = std::filesystem;

TEST_CASE("SimplnxCore::ChangeAngleRepresentationFilter: Invalid Execution", "[OrientationAnalysis][ChangeAngleRepresentationFilter]")
{
  // Configure the filter arguments.
  ChangeAngleRepresentationFilter filter;
  DataStructure dataStructure;
  Arguments args;

  // This should fail
  args.insertOrAssign(ChangeAngleRepresentationFilter::k_ConversionType_Key, std::make_any<ChoicesParameter::ValueType>(0));
  args.insertOrAssign(ChangeAngleRepresentationFilter::k_AnglesArrayPath_Key, std::make_any<DataPath>(DataPath{}));

  auto preflightResult = filter.preflight(dataStructure, args);
  REQUIRE(preflightResult.output_actions.invalid());

  // This should fail because parameter is out of range
  args.insertOrAssign(ChangeAngleRepresentationFilter::k_ConversionType_Key, std::make_any<ChoicesParameter::ValueType>(2));
  preflightResult = filter.preflight(dataStructure, args);
  REQUIRE(preflightResult.output_actions.invalid());

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ChangeAngleRepresentationFilter: Degrees To Radians")
{
  const auto scenario = GENERATE(from_range(unit_test::SelectAlgorithmTestScenariosForInMemoryStores()));
  CAPTURE(scenario);
  unit_test::AlgorithmTestScope scope(scenario);
  unit_test::LoadPlugins();

  // Configure the filter arguments.
  ChangeAngleRepresentationFilter filter;
  DataStructure dataStructure;
  Arguments args;

  DataGroup* topLevelGroupPtr = DataGroup::Create(dataStructure, constants::k_SmallIN100);
  DataGroup* scanDataPtr = DataGroup::Create(dataStructure, constants::k_EbsdScanData, topLevelGroupPtr->getId());

  std::vector<size_t> tupleShape = {10};
  std::vector<size_t> componentShape = {3};

  Float32ArrayType* anglesPtr = unit_test::CreateTestDataArray<float>(dataStructure, constants::k_EulerAngles, tupleShape, componentShape, scanDataPtr->getId());

  for(size_t t = 0; t < tupleShape[0]; t++)
  {
    for(size_t c = 0; c < componentShape[0]; c++)
    {
      (*anglesPtr)[t * componentShape[0] + c] = static_cast<float>(t * c);
    }
  }

  // This should fail
  args.insertOrAssign(ChangeAngleRepresentationFilter::k_ConversionType_Key, std::make_any<ChoicesParameter::ValueType>(0));
  args.insertOrAssign(ChangeAngleRepresentationFilter::k_AnglesArrayPath_Key, std::make_any<DataPath>(DataPath({constants::k_SmallIN100, constants::k_EbsdScanData, constants::k_EulerAngles})));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);

  auto executeResult = scope.executeFilter(filter, dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  // Check the results
  float d2r = nx::core::numbers::k_Pi / 180.0f;
  for(size_t t = 0; t < tupleShape[0]; t++)
  {
    for(size_t c = 0; c < componentShape[0]; c++)
    {
      REQUIRE((*anglesPtr)[t * componentShape[0] + c] == static_cast<float>(t * c) * d2r);
    }
  }

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ChangeAngleRepresentationFilter: Radians To Degrees")
{
  const auto scenario = GENERATE(from_range(unit_test::SelectAlgorithmTestScenariosForInMemoryStores()));
  CAPTURE(scenario);
  unit_test::AlgorithmTestScope scope(scenario);
  unit_test::LoadPlugins();

  // Configure the filter arguments.
  ChangeAngleRepresentationFilter filter;
  DataStructure dataStructure;
  Arguments args;

  DataGroup* topLevelGroupPtr = DataGroup::Create(dataStructure, constants::k_SmallIN100);
  DataGroup* scanDataPtr = DataGroup::Create(dataStructure, constants::k_EbsdScanData, topLevelGroupPtr->getId());

  std::vector<size_t> tupleShape = {10};
  std::vector<size_t> componentShape = {3};

  Float32ArrayType* anglesPtr = unit_test::CreateTestDataArray<float>(dataStructure, constants::k_EulerAngles, tupleShape, componentShape, scanDataPtr->getId());

  for(size_t t = 0; t < tupleShape[0]; t++)
  {
    for(size_t c = 0; c < componentShape[0]; c++)
    {
      (*anglesPtr)[t * componentShape[0] + c] = static_cast<float>(t * c);
    }
  }

  // This should fail
  args.insertOrAssign(ChangeAngleRepresentationFilter::k_ConversionType_Key, std::make_any<ChoicesParameter::ValueType>(1));
  args.insertOrAssign(ChangeAngleRepresentationFilter::k_AnglesArrayPath_Key, std::make_any<DataPath>(DataPath({constants::k_SmallIN100, constants::k_EbsdScanData, constants::k_EulerAngles})));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);

  auto executeResult = scope.executeFilter(filter, dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  // Check the results
  float r2d = 180.0f / nx::core::numbers::k_Pi;
  for(size_t t = 0; t < tupleShape[0]; t++)
  {
    for(size_t c = 0; c < componentShape[0]; c++)
    {
      REQUIRE((*anglesPtr)[t * componentShape[0] + c] == static_cast<float>(t * c) * r2d);
    }
  }

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ChangeAngleRepresentationFilter: SIMPL Backwards Compatibility", "[SimplnxCore][ChangeAngleRepresentationFilter][BackwardsCompatibility]")
{
  const auto version = GENERATE("6_4", "6_5");
  CAPTURE(version);
  const fs::path fixturePath = fs::path(nx::core::unit_test::k_SourceDir.view()) / "test" / "simpl_conversion" / version / "ChangeAngleRepresentationFilter.json";
  const auto args = unit_test::LoadSIMPLConversionArguments(fixturePath, FilterTraits<ChangeAngleRepresentationFilter>::uuid);
  CHECK(args.value<ChoicesParameter::ValueType>(ChangeAngleRepresentationFilter::k_ConversionType_Key) == 0);
  CHECK(args.value<DataPath>(ChangeAngleRepresentationFilter::k_AnglesArrayPath_Key) == DataPath({"DataContainer", "CellData", "TestArray"}));
}
