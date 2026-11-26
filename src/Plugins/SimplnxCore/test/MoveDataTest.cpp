#include "SimplnxCore/Filters/MoveDataFilter.hpp"
#include "SimplnxCore/SimplnxCore_test_dirs.hpp"

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataGroup.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <catch2/catch.hpp>

#include <filesystem>
#include <fstream>
#include <set>

using namespace nx::core;
namespace fs = std::filesystem;

namespace
{
constexpr StringLiteralType k_Group1Name = "Group1";
constexpr StringLiteralType k_Group2Name = "Group2";
constexpr StringLiteralType k_Group3Name = "Group3";
constexpr StringLiteralType k_Group4Name = "Group4";

const int32 k_TupleDimMismatchWarningCode = -69250;

constexpr StringLiteralType k_AM1Name = "Matrix1";
constexpr StringLiteralType k_DataArray1Name = "Array1";

DataStructure CreateDataStructure()
{
  DataStructure data;
  auto* group1Ptr = DataGroup::Create(data, k_Group1Name);
  auto* group2Ptr = DataGroup::Create(data, k_Group2Name);
  auto* group3Ptr = DataGroup::Create(data, k_Group3Name, group2Ptr->getId());
  auto* group4Ptr = DataGroup::Create(data, k_Group4Name, group2Ptr->getId());

  auto* am1Ptr = AttributeMatrix::Create(data, k_AM1Name, std::vector<usize>{2});

  auto dataStore = std::make_unique<DataStore<uint8>>(std::vector<usize>{20}, std::vector<usize>{3}, 0);
  auto* dataArray1Ptr = DataArray<uint8>::Create(data, k_DataArray1Name, std::move(dataStore));

  data.setAdditionalParent(group4Ptr->getId(), group3Ptr->getId());
  return data;
}
} // namespace

TEST_CASE("SimplnxCore::MoveDataFilter Successful", "[Simplnx::Core][MoveDataFilter]")
{
  MoveDataFilter filter;
  Arguments args;
  DataStructure dataStructure = CreateDataStructure();

  const DataPath kGroup1Path({k_Group1Name});
  const DataPath kGroup3Path({k_Group2Name, k_Group3Name});
  const DataPath kGroup4Path({k_Group2Name, k_Group4Name});

  args.insertOrAssign(MoveDataFilter::k_SourceDataPaths_Key, std::make_any<std::vector<DataPath>>({kGroup3Path, kGroup4Path}));
  args.insertOrAssign(MoveDataFilter::k_DestinationParentPath_Key, std::make_any<DataPath>(kGroup1Path));

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions)

  // Execute the filter and check the result
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)

  const auto* group1Ptr = dataStructure.getDataAs<DataGroup>(kGroup1Path);
  REQUIRE(group1Ptr->getDataMap().getSize() == 2);

  const auto* group2Ptr = dataStructure.getDataAs<DataGroup>(DataPath({k_Group2Name}));
  REQUIRE(group2Ptr->getDataMap().getSize() == 0);

  REQUIRE(dataStructure.getDataAs<DataGroup>(kGroup3Path) == nullptr);

  const DataPath newGroup3Path = kGroup1Path.createChildPath(k_Group3Name);
  REQUIRE(dataStructure.getDataAs<DataGroup>(newGroup3Path) != nullptr);

  REQUIRE(dataStructure.getDataAs<DataGroup>(kGroup4Path) == nullptr);

  const DataPath newGroup4Path = kGroup1Path.createChildPath(k_Group4Name);
  REQUIRE(dataStructure.getDataAs<DataGroup>(newGroup4Path) != nullptr);

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::MoveDataFilter Unsuccessful", "[Simplnx::Core][MoveDataFilter]")
{
  MoveDataFilter filter;
  Arguments args;
  DataStructure dataStructure = CreateDataStructure();

  const DataPath kGroup2Path({k_Group2Name});
  const DataPath kGroup3Path({k_Group3Name});

  SECTION("Object already exists in new data path")
  {
    args.insertOrAssign(MoveDataFilter::k_SourceDataPaths_Key, std::make_any<std::vector<DataPath>>({kGroup3Path}));
    args.insertOrAssign(MoveDataFilter::k_DestinationParentPath_Key, std::make_any<DataPath>(kGroup2Path));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.output_actions);
  }
  SECTION("Cannot reparent object to itself")
  {
    args.insertOrAssign(MoveDataFilter::k_SourceDataPaths_Key, std::make_any<std::vector<DataPath>>({kGroup3Path}));
    args.insertOrAssign(MoveDataFilter::k_DestinationParentPath_Key, std::make_any<DataPath>(kGroup3Path));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.output_actions);
  }
  SECTION("Cannot reparent object to a child object")
  {
    args.insertOrAssign(MoveDataFilter::k_SourceDataPaths_Key, std::make_any<std::vector<DataPath>>({kGroup2Path}));
    args.insertOrAssign(MoveDataFilter::k_DestinationParentPath_Key, std::make_any<DataPath>(kGroup3Path));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.output_actions);
  }

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::MoveDataFilter Tuple Size Mismatches Warning and Failure", "[Simplnx::Core][MoveDataFilter]")
{
  MoveDataFilter filter;
  Arguments args;
  DataStructure dataStructure = CreateDataStructure();

  const DataPath kDataArrayPath({k_DataArray1Name});
  const DataPath kAmPath({k_AM1Name});

  args.insertOrAssign(MoveDataFilter::k_SourceDataPaths_Key, std::make_any<std::vector<DataPath>>({kDataArrayPath}));
  args.insertOrAssign(MoveDataFilter::k_DestinationParentPath_Key, std::make_any<DataPath>(kAmPath));

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(dataStructure, args);

  REQUIRE_FALSE(preflightResult.output_actions.warnings().empty());

  bool found = false;
  for(const auto& warning : preflightResult.output_actions.warnings())
  {
    if(warning.code == ::k_TupleDimMismatchWarningCode)
    {
      found = true;
    }
  }

  REQUIRE(found);

  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);

  auto result = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(result.result);

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::MoveDataFilter: SIMPL Backwards Compatibility", "[SimplnxCore][MoveDataFilter][BackwardsCompatibility]")
{
  const auto version = GENERATE("6_4", "6_5");
  CAPTURE(version);
  const fs::path fixturePath = fs::path(nx::core::unit_test::k_SourceDir.view()) / "test" / "simpl_conversion" / version / "MoveDataFilter.json";
  const auto args = unit_test::LoadSIMPLConversionArguments(fixturePath, FilterTraits<MoveDataFilter>::uuid);
  // Fixture has WhatToMove=0 (move attribute matrix): source is the AttributeMatrixSource path,
  // destination is the DataContainerDestination path.
  CHECK(args.value<DataPath>(MoveDataFilter::k_DestinationParentPath_Key) == DataPath({"DataContainer"}));
  CHECK(args.value<std::vector<DataPath>>(MoveDataFilter::k_SourceDataPaths_Key) == std::vector<DataPath>{DataPath({"DataContainer", "CellData"})});
}
