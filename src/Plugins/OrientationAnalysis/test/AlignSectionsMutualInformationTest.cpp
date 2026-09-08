#include <catch2/catch.hpp>

#include "OrientationAnalysis/Filters/AlignSectionsMutualInformationFilter.hpp"
#include "OrientationAnalysis/OrientationAnalysis_test_dirs.hpp"

#include "simplnx/Core/Application.hpp"
#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Pipeline/Pipeline.hpp"
#include "simplnx/Pipeline/PipelineFilter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Utilities/DataStoreUtilities.hpp"

#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;

using namespace nx::core;

TEST_CASE("OrientationAnalysis::AlignSectionsMutualInformationFilter: Valid filter execution")
{
  UnitTest::LoadPlugins();
  const nx::core::UnitTest::TestFileSentinel testDataSentinel(nx::core::unit_test::k_TestFilesDir, "align_sections_mutual_information.tar.gz", "align_sections_mutual_information");

  // We are just going to generate a big number so that we can use that in the output
  // file path. This tests the creation of intermediate directories that the filter
  // would be responsible to create.
  const uint64_t millisFromEpoch = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  auto* filterList = Application::Instance()->getFilterList();

  const DataPath k_ExemplarShiftsPath = Constants::k_ExemplarDataContainerPath.createChildPath("Exemplar Shifts");

  // Read Exemplar DREAM3D File Filter
  auto exemplarFilePath = fs::path(fmt::format("{}/align_sections_mutual_information/6_5_align_sections_mutual_information.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(exemplarFilePath);

  // Align Sections Mutual Information Filter
  {
    // Instantiate the filter, a DataStructure object and an Arguments Object
    AlignSectionsMutualInformationFilter filter;
    Arguments args;

    // Create valid Parameters for the filter.
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_MisorientationTolerance_Key, std::make_any<float32>(5.0f));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_UseMask_Key, std::make_any<bool>(true));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(Constants::k_DataContainerPath));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_QuatsArrayPath_Key, std::make_any<DataPath>(Constants::k_QuatsArrayPath));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_CellPhasesArrayPath_Key, std::make_any<DataPath>(Constants::k_PhasesArrayPath));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_MaskArrayPath_Key, std::make_any<DataPath>(Constants::k_MaskArrayPath));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_CrystalStructuresArrayPath_Key, std::make_any<DataPath>(Constants::k_CrystalStructuresArrayPath));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)

    // Execute the filter and check the result
    auto executeResult = filter.execute(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
  }

  UnitTest::CompareExemplarToGeneratedData(dataStructure, dataStructure, Constants::k_CellAttributeMatrix, Constants::k_ExemplarDataContainer);

#ifdef SIMPLNX_WRITE_TEST_OUTPUT
  UnitTest::WriteTestDataStructure(dataStructure, fmt::format("{}/align_sections_mutual_information.dream3d", unit_test::k_BinaryTestOutputDir));
#endif

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("OrientationAnalysis::AlignSectionsMutualInformationFilter: InValid filter execution")
{
  UnitTest::LoadPlugins();
  const nx::core::UnitTest::TestFileSentinel testDataSentinel(nx::core::unit_test::k_TestFilesDir, "6_6_stats_test_v2.tar.gz", "6_6_stats_test_v2.dream3d");

  // Read the Small IN100 Data set
  auto baseDataFilePath = fs::path(fmt::format("{}/6_6_stats_test_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);
  // Instantiate the filter and an Arguments Object
  AlignSectionsMutualInformationFilter filter;
  Arguments args;

  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_MisorientationTolerance_Key, std::make_any<float32>(5.0f));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_UseMask_Key, std::make_any<bool>(true));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(Constants::k_DataContainerPath));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_QuatsArrayPath_Key, std::make_any<DataPath>(Constants::k_QuatsArrayPath));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_CellPhasesArrayPath_Key, std::make_any<DataPath>(Constants::k_PhasesArrayPath));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_MaskArrayPath_Key,
                      std::make_any<DataPath>(DataPath({Constants::k_DataContainer, Constants::k_CellFeatureData, Constants::k_ActiveName})));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_CrystalStructuresArrayPath_Key, std::make_any<DataPath>(Constants::k_CrystalStructuresArrayPath));

  SECTION("Mismatching cell data tuples")
  {
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_MaskArrayPath_Key,
                        std::make_any<DataPath>(DataPath({Constants::k_DataContainer, Constants::k_CellFeatureData, Constants::k_ActiveName})));
  }

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.outputActions);
  REQUIRE(preflightResult.outputActions.errors()[0].code == -3542);

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("OrientationAnalysis::AlignSectionsMutualInformationFilter: output test", "[Reconstruction][AlignSectionsMutualInformationFilter]")
{
  const nx::core::UnitTest::TestFileSentinel testDataSentinel(nx::core::unit_test::k_TestFilesDir, "align_sections_mutual_information.tar.gz", "align_sections_mutual_information");

  // Read Exemplar DREAM3D File Filter
  auto baseFilePath = fs::path(fmt::format("{}/align_sections_mutual_information/6_5_align_sections_mutual_information.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseFilePath);

  // Align Sections Mutual Information Filter
  {
    // Instantiate the filter, a DataStructure object and an Arguments Object
    AlignSectionsMutualInformationFilter filter;
    Arguments args;

    // Create valid Parameters for the filter.
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_MisorientationTolerance_Key, std::make_any<float32>(5.0f));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_UseMask_Key, std::make_any<bool>(true));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(Constants::k_DataContainerPath));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_QuatsArrayPath_Key, std::make_any<DataPath>(Constants::k_QuatsArrayPath));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_CellPhasesArrayPath_Key, std::make_any<DataPath>(Constants::k_PhasesArrayPath));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_MaskArrayPath_Key, std::make_any<DataPath>(Constants::k_MaskArrayPath));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_CrystalStructuresArrayPath_Key, std::make_any<DataPath>(Constants::k_CrystalStructuresArrayPath));

    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_StoreAlignmentShifts_Key, std::make_any<bool>(true));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_AlignmentAMName_Key, std::make_any<std::string>(Constants::k_AlignmentAMName));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_SlicesArrayName_Key, std::make_any<std::string>(Constants::k_SlicesArrayName));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_RelativeShiftsArrayName_Key, std::make_any<std::string>(Constants::k_RelativeShiftsArrayName));
    args.insertOrAssign(AlignSectionsMutualInformationFilter::k_CumulativeShiftsArrayName_Key, std::make_any<std::string>(Constants::k_CumulativeShiftsArrayName));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)

    // Execute the filter and check the result
    auto executeResult = filter.execute(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
  }

  // Read Exemplar data structure
  auto exemplarFilePath = fs::path(fmt::format("{}/align_sections_mutual_information/output_align_sections_mutual_information.dream3d", unit_test::k_TestFilesDir));
  DataStructure exemplarDataStructure = UnitTest::LoadDataStructure(exemplarFilePath);

  const DataPath alignmentAMPath = Constants::k_DataContainerPath.createChildPath(Constants::k_AlignmentAMName);

  const DataPath slicesPath = alignmentAMPath.createChildPath(Constants::k_SlicesArrayName);
  UnitTest::CompareDataArrays<uint32>(exemplarDataStructure.getDataRefAs<IDataArray>(slicesPath), dataStructure.getDataRefAs<IDataArray>(slicesPath));

  const DataPath relativeShiftsPath = alignmentAMPath.createChildPath(Constants::k_RelativeShiftsArrayName);
  UnitTest::CompareDataArrays<int64>(exemplarDataStructure.getDataRefAs<IDataArray>(relativeShiftsPath), dataStructure.getDataRefAs<IDataArray>(relativeShiftsPath));

  const DataPath cumulativeShiftsPath = alignmentAMPath.createChildPath(Constants::k_CumulativeShiftsArrayName);
  UnitTest::CompareDataArrays<int64>(exemplarDataStructure.getDataRefAs<IDataArray>(cumulativeShiftsPath), dataStructure.getDataRefAs<IDataArray>(cumulativeShiftsPath));

// Write the output .dream3d file when test-output generation is enabled.
#ifdef SIMPLNX_WRITE_TEST_OUTPUT
  UnitTest::WriteTestDataStructure(dataStructure, fmt::format("{}/output_align_sections_mutual_information.dream3d", unit_test::k_BinaryTestOutputDir));
#endif

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("OrientationAnalysis::AlignSectionsMutualInformationFilter: SIMPL Backwards Compatibility", "[OrientationAnalysis][AlignSectionsMutualInformationFilter][BackwardsCompatibility]")
{
  auto app = Application::GetOrCreateInstance();
  UnitTest::LoadPlugins();
  auto filterList = app->getFilterList();

  const fs::path conversionDir = fs::path(nx::core::unit_test::k_SourceDir.view()) / "test" / "simpl_conversion";

  const std::vector<std::pair<std::string, fs::path>> fixtures = {
      {"SIMPL 6.5 (UUID)", conversionDir / "6_5" / "AlignSectionsMutualInformationFilter.json"},
      {"SIMPL 6.4 (Filter_Name)", conversionDir / "6_4" / "AlignSectionsMutualInformationFilter.json"},
  };

  for(const auto& [label, fixturePath] : fixtures)
  {
    DYNAMIC_SECTION(label)
    {
      auto pipelineResult = Pipeline::FromSIMPLFile(fixturePath, filterList);
      REQUIRE(pipelineResult.valid());

      auto& pipeline = pipelineResult.value();
      REQUIRE(pipeline.size() == 1);

      auto* pipelineFilter = dynamic_cast<PipelineFilter*>(pipeline.at(0));
      REQUIRE(pipelineFilter != nullptr);

      const IFilter* filter = pipelineFilter->getFilter();
      REQUIRE(filter != nullptr);
      REQUIRE(filter->uuid() == FilterTraits<AlignSectionsMutualInformationFilter>::uuid);

      CHECK(pipelineFilter->getComments().empty());

      const Arguments args = pipelineFilter->getArguments();
      CHECK(args.value<bool>(AlignSectionsMutualInformationFilter::k_StoreAlignmentShifts_Key) == true);
      CHECK(args.value<float32>(AlignSectionsMutualInformationFilter::k_MisorientationTolerance_Key) == 2.5f);
      CHECK(args.value<bool>(AlignSectionsMutualInformationFilter::k_UseMask_Key) == true);
      CHECK(args.value<DataPath>(AlignSectionsMutualInformationFilter::k_QuatsArrayPath_Key) == DataPath({"DataContainer", "CellData", "TestArray"}));
      CHECK(args.value<DataPath>(AlignSectionsMutualInformationFilter::k_CellPhasesArrayPath_Key) == DataPath({"DataContainer", "CellData", "TestArray"}));
      CHECK(args.value<DataPath>(AlignSectionsMutualInformationFilter::k_MaskArrayPath_Key) == DataPath({"DataContainer", "CellData", "TestArray"}));
      CHECK(args.value<DataPath>(AlignSectionsMutualInformationFilter::k_CrystalStructuresArrayPath_Key) == DataPath({"DataContainer", "CellData", "TestArray"}));
      CHECK(args.value<DataPath>(AlignSectionsMutualInformationFilter::k_SelectedImageGeometryPath_Key) == DataPath({"DataContainer"}));
    }
  }
}

TEST_CASE("OrientationAnalysis::AlignSectionsMutualInformationFilter: Execute Error - Phase and Laue Index Bounds", "[OrientationAnalysis][AlignSectionsMutualInformationFilter]")
{
  UnitTest::LoadPlugins();

  constexpr usize k_DimX = 4;
  constexpr usize k_DimY = 4;
  constexpr usize k_DimZ = 2;
  constexpr usize k_TotalVoxels = k_DimX * k_DimY * k_DimZ;
  constexpr usize k_NumCrystalStructures = 2;
  const bool invalidPhaseIdx = GENERATE(false, true);
  CAPTURE(invalidPhaseIdx);

  DataStructure dataStructure;
  auto* imageGeomPtr = ImageGeom::Create(dataStructure, "ImageGeometry");
  imageGeomPtr->setDimensions({k_DimX, k_DimY, k_DimZ});
  imageGeomPtr->setSpacing({1.0F, 1.0F, 1.0F});
  imageGeomPtr->setOrigin({0.0F, 0.0F, 0.0F});

  const ShapeType cellTupleShape = {k_DimZ, k_DimY, k_DimX};
  auto* cellDataPtr = AttributeMatrix::Create(dataStructure, "CellData", cellTupleShape, imageGeomPtr->getId());
  imageGeomPtr->setCellData(*cellDataPtr);
  auto* ensembleDataPtr = AttributeMatrix::Create(dataStructure, "CellEnsembleData", ShapeType{k_NumCrystalStructures}, imageGeomPtr->getId());

  auto quatsStore = DataStoreUtilities::CreateDataStore<float32>(cellTupleShape, {4}, IDataAction::Mode::Execute);
  auto* quatsArrayPtr = Float32Array::Create(dataStructure, "Quats", quatsStore, cellDataPtr->getId());
  auto cellPhasesStore = DataStoreUtilities::CreateDataStore<int32>(cellTupleShape, {1}, IDataAction::Mode::Execute);
  auto* cellPhasesArrayPtr = Int32Array::Create(dataStructure, "Phases", cellPhasesStore, cellDataPtr->getId());
  auto maskStore = DataStoreUtilities::CreateDataStore<bool>(cellTupleShape, {1}, IDataAction::Mode::Execute);
  auto* maskArrayPtr = BoolArray::Create(dataStructure, "Mask", maskStore, cellDataPtr->getId());
  auto crystalStructuresStore = DataStoreUtilities::CreateDataStore<uint32>(ShapeType{k_NumCrystalStructures}, {1}, IDataAction::Mode::Execute);
  auto* crystalStructuresArrayPtr = UInt32Array::Create(dataStructure, "CrystalStructures", crystalStructuresStore, ensembleDataPtr->getId());

  quatsArrayPtr->fill(0.0F);
  for(usize voxelIdx = 0; voxelIdx < k_TotalVoxels; ++voxelIdx)
  {
    (*quatsArrayPtr)[voxelIdx * 4 + 3] = 1.0F;
  }
  cellPhasesArrayPtr->fill(invalidPhaseIdx ? static_cast<int32>(k_NumCrystalStructures) : 1);
  maskArrayPtr->fill(true);
  (*crystalStructuresArrayPtr)[0] = 999U;
  (*crystalStructuresArrayPtr)[1] = invalidPhaseIdx ? 1U : 999U;

  const DataPath imageGeomPath({"ImageGeometry"});
  const DataPath cellDataPath = imageGeomPath.createChildPath("CellData");
  const DataPath ensembleDataPath = imageGeomPath.createChildPath("CellEnsembleData");

  AlignSectionsMutualInformationFilter filter;
  Arguments args = filter.getDefaultArguments();
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(imageGeomPath));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_QuatsArrayPath_Key, std::make_any<DataPath>(cellDataPath.createChildPath("Quats")));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_CellPhasesArrayPath_Key, std::make_any<DataPath>(cellDataPath.createChildPath("Phases")));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_MaskArrayPath_Key, std::make_any<DataPath>(cellDataPath.createChildPath("Mask")));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_CrystalStructuresArrayPath_Key, std::make_any<DataPath>(ensembleDataPath.createChildPath("CrystalStructures")));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_UseMask_Key, std::make_any<bool>(false));
  args.insertOrAssign(AlignSectionsMutualInformationFilter::k_StoreAlignmentShifts_Key, std::make_any<bool>(false));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(executeResult.result);
  REQUIRE(executeResult.result.errors()[0].code == (invalidPhaseIdx ? -53703 : -53704));

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}
