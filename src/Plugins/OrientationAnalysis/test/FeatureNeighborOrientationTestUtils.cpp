#include "FeatureNeighborOrientationTestUtils.hpp"

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataGroup.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Utilities/DataStoreUtilities.hpp"

#include <EbsdLib/Core/EbsdLibConstants.h>
#include <catch2/catch.hpp>

#include <algorithm>
#include <memory>
#include <vector>

namespace nx::core::unit_test::FeatureNeighbors
{
InputPaths::InputPaths(const DataPath& featureMatrix, const DataPath& ensembleMatrix)
: feature_phases(featureMatrix.createChildPath("FeaturePhases"))
, avg_quats(featureMatrix.createChildPath("AvgQuats"))
, neighbor_list(featureMatrix.createChildPath("NeighborList"))
, crystal_structures(ensembleMatrix.createChildPath("CrystalStructures"))
{
}

void InitializeFeatures(Int32ArrayType& phases, Float32ArrayType& quats, NeighborList<int32>& neighbors)
{
  const usize numFeatures = phases.getNumberOfTuples();
  for(usize i = 0; i < numFeatures; ++i)
  {
    phases[i] = 0;
    quats[i * 4 + 0] = 0.0f;
    quats[i * 4 + 1] = 0.0f;
    quats[i * 4 + 2] = 0.0f;
    quats[i * 4 + 3] = 1.0f;
    neighbors.setList(i, std::make_shared<std::vector<int32>>());
  }
}

void CreateChunkTailInputs(DataStructure& dataStructure, const InputPaths& paths, uint32 crystalStructure, const std::array<float32, 4>& lastFullChunkQuat, const std::array<float32, 4>& tailQuat)
{
  constexpr usize k_FeatureTuples = k_TailTuple + 1;
  const DataPath featureMatrixPath = paths.feature_phases.getParent();
  auto* dataGroup = DataGroup::Create(dataStructure, featureMatrixPath.getParent().getTargetName());
  REQUIRE(dataGroup != nullptr);
  auto* featureData = AttributeMatrix::Create(dataStructure, featureMatrixPath.getTargetName(), {k_FeatureTuples}, dataGroup->getId());
  auto* ensembleDataPtr = AttributeMatrix::Create(dataStructure, paths.crystal_structures.getParent().getTargetName(), {2}, dataGroup->getId());
  REQUIRE(featureData != nullptr);
  REQUIRE(ensembleDataPtr != nullptr);

  auto featurePhasesStore = data_store_utilities::CreateDataStore<int32>(dataStructure, paths.feature_phases, {k_FeatureTuples}, {1});
  auto avgQuatsStore = data_store_utilities::CreateDataStore<float32>(dataStructure, paths.avg_quats, {k_FeatureTuples}, {4});
  auto crystalStructuresStore = data_store_utilities::CreateDataStore<uint32>(dataStructure, paths.crystal_structures, {2}, {1});
  auto* featurePhases = Int32ArrayType::Create(dataStructure, paths.feature_phases.getTargetName(), featurePhasesStore, featureData->getId());
  auto* avgQuats = Float32ArrayType::Create(dataStructure, paths.avg_quats.getTargetName(), avgQuatsStore, featureData->getId());
  auto* neighborList = NeighborList<int32>::Create(dataStructure, paths.neighbor_list.getTargetName(), {k_FeatureTuples}, featureData->getId());
  auto* crystalStructures = UInt32ArrayType::Create(dataStructure, paths.crystal_structures.getTargetName(), crystalStructuresStore, ensembleDataPtr->getId());
  REQUIRE(featurePhases != nullptr);
  REQUIRE(avgQuats != nullptr);
  REQUIRE(neighborList != nullptr);
  REQUIRE(crystalStructures != nullptr);

  // Write the resident fixture buffers in blocks, including the partial final chunk.
  constexpr usize k_BlockTuples = 1024;
  std::array<int32, k_BlockTuples> phases;
  phases.fill(1);
  std::array<float32, 4 * k_BlockTuples> quats{};
  for(usize i = 3; i < quats.size(); i += 4)
  {
    quats[i] = 1.0F;
  }
  for(usize start = 0; start < k_FeatureTuples; start += k_BlockTuples)
  {
    const usize count = std::min(k_BlockTuples, k_FeatureTuples - start);
    phases[0] = start == 0 ? 0 : 1;
    const auto phaseWrite = featurePhases->getDataStoreRef().copyFromBuffer(start, nonstd::span<const int32>(phases.data(), count));
    SIMPLNX_RESULT_REQUIRE_VALID(phaseWrite);
    const auto quatWrite = avgQuats->getDataStoreRef().copyFromBuffer(start * 4, nonstd::span<const float32>(quats.data(), count * 4));
    SIMPLNX_RESULT_REQUIRE_VALID(quatWrite);
  }
  const auto lastFullWrite = avgQuats->getDataStoreRef().copyFromBuffer(k_LastFullChunkTuple * 4, nonstd::span<const float32>(lastFullChunkQuat));
  SIMPLNX_RESULT_REQUIRE_VALID(lastFullWrite);
  const auto tailWrite = avgQuats->getDataStoreRef().copyFromBuffer(k_TailTuple * 4, nonstd::span<const float32>(tailQuat));
  SIMPLNX_RESULT_REQUIRE_VALID(tailWrite);
  const std::array<uint32, 2> crystalStructureValues = {ebsdlib::CrystalStructure::UnknownCrystalStructure, crystalStructure};
  auto crystalStructuresWriteResult = crystalStructures->getDataStoreRef().copyFromBuffer(0, nonstd::span<const uint32>(crystalStructureValues.data(), crystalStructureValues.size()));
  SIMPLNX_RESULT_REQUIRE_VALID(crystalStructuresWriteResult);

  neighborList->setList(k_LastFullChunkTuple, std::vector<int32>{1});
  neighborList->setList(k_TailTuple, std::vector<int32>{1});

  REQUIRE(featurePhases->getDataStoreRef().getDataFormat() == "HDF5-OOC");
  REQUIRE(avgQuats->getDataStoreRef().getDataFormat() == "HDF5-OOC");
  const auto quaternionMetadata = avgQuats->getDataStoreRef().getRecoveryMetadata();
  REQUIRE(quaternionMetadata.at("OocChunkShape") == "65536");
  REQUIRE(crystalStructures->getDataStoreRef().getDataFormat() == "HDF5-OOC");
}
} // namespace nx::core::unit_test::FeatureNeighbors
