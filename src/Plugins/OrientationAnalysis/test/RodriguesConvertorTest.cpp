#include <array>
#include <catch2/catch.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>

#include "simplnx/Core/Application.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Utilities/AlgorithmDispatch.hpp"
#include "simplnx/Utilities/DataArrayUtilities.hpp"
#include "simplnx/Utilities/DataStoreUtilities.hpp"

#include "OrientationAnalysis/Filters/RodriguesConvertorFilter.hpp"
#include "OrientationAnalysis/OrientationAnalysis_test_dirs.hpp"

using namespace nx::core;
namespace fs = std::filesystem;
namespace
{
const std::string k_InputArrayName = "Input";
const std::string k_ConvertedName = "Output";
const std::string k_Exemplar0 = "Exemplar0";

} // namespace

TEST_CASE("OrientationAnalysis::RodriguesConvertorFilter", "[OrientationAnalysis][RodriguesConvertorFilter]")
{
  unit_test::LoadPlugins();
  // AlgorithmTestScope forces the selected path and records its target-call
  // witness.
  const auto scenario = GENERATE(from_range(unit_test::SelectAlgorithmTestScenariosForInMemoryStores()));
  CAPTURE(scenario);
  unit_test::AlgorithmTestScope scope(scenario);

  DataStructure dataStructure;

  Float32ArrayType* quatsPtr = unit_test::CreateTestDataArray<float32>(dataStructure, k_InputArrayName, {4ULL}, {3ULL}, {});
  scope.requireExpectedStore(*quatsPtr);

  for(usize i = 0; i < 12; i++)
  {
    (*quatsPtr)[i] = static_cast<float32>(i);
  }
  Float32ArrayType* exemplarDataPtr = unit_test::CreateTestDataArray<float32>(dataStructure, k_Exemplar0, {4ULL}, {4ULL}, {});
  (*exemplarDataPtr)[0] = 0.0F;
  (*exemplarDataPtr)[1] = 0.447214F;
  (*exemplarDataPtr)[2] = 0.894427F;
  (*exemplarDataPtr)[3] = 2.23607F;
  (*exemplarDataPtr)[4] = 0.424264F;
  (*exemplarDataPtr)[5] = 0.565685F;
  (*exemplarDataPtr)[6] = 0.707107F;
  (*exemplarDataPtr)[7] = 7.07107F;
  (*exemplarDataPtr)[8] = 0.491539F;
  (*exemplarDataPtr)[9] = 0.573462F;
  (*exemplarDataPtr)[10] = 0.655386F;
  (*exemplarDataPtr)[11] = 12.2066F;
  (*exemplarDataPtr)[12] = 0.517893F;
  (*exemplarDataPtr)[13] = 0.575437F;
  (*exemplarDataPtr)[14] = 0.632980F;
  (*exemplarDataPtr)[15] = 17.3781F;
  {
    const RodriguesConvertorFilter filter;

    Arguments args;

    args.insertOrAssign(RodriguesConvertorFilter::k_RodriguesDataArrayPath_Key, std::make_any<DataPath>(DataPath({k_InputArrayName})));
    args.insertOrAssign(RodriguesConvertorFilter::k_OutputDataArrayName_Key, std::make_any<std::string>(k_ConvertedName));
    args.insertOrAssign(RodriguesConvertorFilter::k_DeleteOriginalData_Key, std::make_any<bool>(false));

    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions)

    auto executeResult = scope.executeFilter(filter, dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)

    REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float32ArrayType>(DataPath({k_ConvertedName})));
    auto& outputArray = dataStructure.getDataRefAs<Float32ArrayType>(DataPath({k_ConvertedName}));
    scope.requireExpectedStore(outputArray);

    unit_test::CompareDataArrays<float32>(*exemplarDataPtr, outputArray);
  }

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("OrientationAnalysis::RodriguesConvertorFilter: SIMPL Backwards Compatibility", "[OrientationAnalysis][RodriguesConvertorFilter][BackwardsCompatibility]")
{
  const auto version = GENERATE("6_4", "6_5");
  CAPTURE(version);
  const fs::path fixturePath = fs::path(nx::core::unit_test::k_SourceDir.view()) / "test" / "simpl_conversion" / version / "RodriguesConvertorFilter.json";
  const auto args = unit_test::LoadSIMPLConversionArguments(fixturePath, FilterTraits<RodriguesConvertorFilter>::uuid);
  CHECK(args.value<DataPath>(RodriguesConvertorFilter::k_RodriguesDataArrayPath_Key) == DataPath({"DataContainer", "CellData", "TestArray"}));
  CHECK(args.value<std::string>(RodriguesConvertorFilter::k_OutputDataArrayName_Key) == "TestArray");
  CHECK(args.value<bool>(RodriguesConvertorFilter::k_DeleteOriginalData_Key) == true);
}
