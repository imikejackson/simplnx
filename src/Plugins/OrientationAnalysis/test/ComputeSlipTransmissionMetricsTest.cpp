#include <catch2/catch.hpp>

#include "OrientationAnalysis/Filters/ComputeSlipTransmissionMetricsFilter.hpp"
#include "OrientationAnalysis/OrientationAnalysis_test_dirs.hpp"

#include "simplnx/UnitTest/UnitTestCommon.hpp"

using namespace nx::core;
using namespace nx::core::UnitTest;

namespace
{
DataPath grainDataPath = DataPath({Constants::k_SmallIN100, Constants::k_Grain_Data});
DataPath neighborListPath = grainDataPath.createChildPath("NeighborList");
DataPath avgQuatsPath = DataPath({Constants::k_SmallIN100, Constants::k_Grain_Data, Constants::k_AvgQuats});
DataPath featurePhasesPath = DataPath({Constants::k_SmallIN100, Constants::k_Grain_Data, Constants::k_Phases});
DataPath crystalStructuresPath = DataPath({Constants::k_SmallIN100, Constants::k_Phase_Data, Constants::k_CrystalStructures});

const std::string k_f1s = "F1List_COMPUTED";
const std::string k_f1spts = "F1sptList_COMPUTED";
const std::string k_f7s = "F7List_COMPUTED";
const std::string k_mPrimes = "mPrimeList_COMPUTED";
} // namespace

TEST_CASE("OrientationAnalysis::ComputeSlipTransmissionMetricsFilter: Valid Filter Execution", "[OrientationAnalysis][ComputeSlipTransmissionMetricsFilter]")
{
  Application::GetOrCreateInstance()->loadPlugins(unit_test::k_BuildDir.view(), true);

  const nx::core::UnitTest::TestFileSentinel testDataSentinel(nx::core::unit_test::k_CMakeExecutable, nx::core::unit_test::k_TestFilesDir, "feature_boundary_neighbor_slip_transmission_1.tar.gz",
                                                              "feature_boundary_neighbor_slip_transmission_1");

  DataStructure dataStructure =
      UnitTest::LoadDataStructure(fs::path(fmt::format("{}/feature_boundary_neighbor_slip_transmission_1/6_6_feature_boundary_neighbor_slip_transmission.dream3d", unit_test::k_TestFilesDir)));
  {
    // Instantiate the filter, a DataStructure object and an Arguments Object
    ComputeSlipTransmissionMetricsFilter filter;
    Arguments args;

    // Create default Parameters for the filter.
    args.insertOrAssign(ComputeSlipTransmissionMetricsFilter::k_NeighborListArrayPath_Key, std::make_any<DataPath>(neighborListPath));
    args.insertOrAssign(ComputeSlipTransmissionMetricsFilter::k_AvgQuatsArrayPath_Key, std::make_any<DataPath>(avgQuatsPath));
    args.insertOrAssign(ComputeSlipTransmissionMetricsFilter::k_FeaturePhasesArrayPath_Key, std::make_any<DataPath>(featurePhasesPath));
    args.insertOrAssign(ComputeSlipTransmissionMetricsFilter::k_CrystalStructuresArrayPath_Key, std::make_any<DataPath>(crystalStructuresPath));

    args.insertOrAssign(ComputeSlipTransmissionMetricsFilter::k_F1ListArrayName_Key, std::make_any<std::string>(k_f1s));
    args.insertOrAssign(ComputeSlipTransmissionMetricsFilter::k_F1sptListArrayName_Key, std::make_any<std::string>(k_f1spts));
    args.insertOrAssign(ComputeSlipTransmissionMetricsFilter::k_F7ListArrayName_Key, std::make_any<std::string>(k_f7s));
    args.insertOrAssign(ComputeSlipTransmissionMetricsFilter::k_mPrimeListArrayName_Key, std::make_any<std::string>(k_mPrimes));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    REQUIRE(preflightResult.outputActions.valid());

    // Execute the filter and check the result
    auto executeResult = filter.execute(dataStructure, args);
    REQUIRE(executeResult.result.valid());
  }

  // #ifdef SIMPLNX_WRITE_TEST_OUTPUT
  WriteTestDataStructure(dataStructure, fs::path(fmt::format("{}/7_0_ComputeSlipTransmissionMetricsFilter.dream3d", unit_test::k_BinaryTestOutputDir)));
  // #endif

  UnitTest::CompareNeighborLists<float32>(dataStructure, grainDataPath.createChildPath("F1List"), grainDataPath.createChildPath(k_f1s));
  UnitTest::CompareNeighborLists<float32>(dataStructure, grainDataPath.createChildPath("F1sptList"), grainDataPath.createChildPath(k_f1spts));
  UnitTest::CompareNeighborLists<float32>(dataStructure, grainDataPath.createChildPath("F7List"), grainDataPath.createChildPath(k_f7s));
  UnitTest::CompareNeighborLists<float32>(dataStructure, grainDataPath.createChildPath("mPrimeList"), grainDataPath.createChildPath(k_mPrimes));
}
