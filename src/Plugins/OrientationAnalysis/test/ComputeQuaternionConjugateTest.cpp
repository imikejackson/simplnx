#include <array>
#include <catch2/catch.hpp>
#include <filesystem>
#include <fstream>
#include <memory>

#include <nonstd/span.hpp>

#include "simplnx/Core/Application.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Utilities/AlgorithmDispatch.hpp"
#include "simplnx/Utilities/DataArrayUtilities.hpp"
#include "simplnx/Utilities/DataStoreUtilities.hpp"

#include "OrientationAnalysis/Filters/ComputeQuaternionConjugateFilter.hpp"
#include "OrientationAnalysis/OrientationAnalysis_test_dirs.hpp"

using namespace nx::core;
namespace fs = std::filesystem;
namespace
{
const std::string k_QuatName = "Quats";
const std::string k_ConvertedName = "Converted";
const std::string k_Exemplar0 = "Exemplar0";

} // namespace

TEST_CASE("OrientationAnalysis::ComputeQuaternionConjugateFilter", "[OrientationAnalysis][ComputeQuaternionConjugateFilter]")
{
  unit_test::LoadPlugins();

  // AlgorithmTestScope forces the selected path and records its target-call
  // witness.
  const auto scenario = GENERATE(from_range(unit_test::SelectAlgorithmTestScenariosForInMemoryStores()));
  CAPTURE(scenario);
  unit_test::AlgorithmTestScope scope(scenario);

  DataStructure dataStructure;

  Float32ArrayType* quatsPtr = unit_test::CreateTestDataArray<float32>(dataStructure, k_QuatName, {4ULL}, {4ULL}, {});

  for(size_t i = 0; i < 16; i++)
  {
    (*quatsPtr)[i] = static_cast<float32>(i);
  }
  Float32ArrayType* exemplarDataPtr = unit_test::CreateTestDataArray<float32>(dataStructure, k_Exemplar0, {4ULL}, {4ULL}, {});

  (*exemplarDataPtr)[0] = -0.0F;
  (*exemplarDataPtr)[1] = -1.0F;
  (*exemplarDataPtr)[2] = -2.0F;
  (*exemplarDataPtr)[3] = 3.0F;

  (*exemplarDataPtr)[4] = -4.0F;
  (*exemplarDataPtr)[5] = -5.0F;
  (*exemplarDataPtr)[6] = -6.0F;
  (*exemplarDataPtr)[7] = 7.0F;

  (*exemplarDataPtr)[8] = -8.0F;
  (*exemplarDataPtr)[9] = -9.0F;
  (*exemplarDataPtr)[10] = -10.0F;
  (*exemplarDataPtr)[11] = 11.0F;

  (*exemplarDataPtr)[12] = -12.0F;
  (*exemplarDataPtr)[13] = -13.0F;
  (*exemplarDataPtr)[14] = -14.0F;
  (*exemplarDataPtr)[15] = 15.0F;
  {
    const ComputeQuaternionConjugateFilter filter;

    Arguments args;

    args.insertOrAssign(ComputeQuaternionConjugateFilter::k_CellQuatsArrayPath_Key, std::make_any<DataPath>(DataPath({k_QuatName})));
    args.insertOrAssign(ComputeQuaternionConjugateFilter::k_OutputDataArrayName_Key, std::make_any<std::string>(k_ConvertedName));
    args.insertOrAssign(ComputeQuaternionConjugateFilter::k_DeleteOriginalData_Key, std::make_any<bool>(false));

    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions)

    auto executeResult = scope.executeFilter(filter, dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)

    auto& outputArray = dataStructure.getDataRefAs<Float32ArrayType>(DataPath({k_ConvertedName}));
    unit_test::CompareDataArrays<float32>(*exemplarDataPtr, outputArray);
  }

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("OrientationAnalysis::ComputeQuaternionConjugateFilter: SIMPL Backwards Compatibility", "[OrientationAnalysis][ComputeQuaternionConjugateFilter][BackwardsCompatibility]")
{
  const auto version = GENERATE("6_4", "6_5");
  CAPTURE(version);
  const fs::path fixturePath = fs::path(nx::core::unit_test::k_SourceDir.view()) / "test" / "simpl_conversion" / version / "ComputeQuaternionConjugateFilter.json";
  const auto args = unit_test::LoadSIMPLConversionArguments(fixturePath, FilterTraits<ComputeQuaternionConjugateFilter>::uuid);
  CHECK(args.value<DataPath>(ComputeQuaternionConjugateFilter::k_CellQuatsArrayPath_Key) == DataPath({"DataContainer", "CellData", "TestArray"}));
  CHECK(args.value<std::string>(ComputeQuaternionConjugateFilter::k_OutputDataArrayName_Key) == "TestArray");
  CHECK(args.value<bool>(ComputeQuaternionConjugateFilter::k_DeleteOriginalData_Key) == true);
}
