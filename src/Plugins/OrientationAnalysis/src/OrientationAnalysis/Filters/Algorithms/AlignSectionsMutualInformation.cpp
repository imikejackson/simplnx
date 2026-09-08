#include "AlignSectionsMutualInformation.hpp"

#include "simplnx/Common/Constants.hpp"
#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/IGridGeometry.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Utilities/FilterUtilities.hpp"
#include "simplnx/Utilities/StringUtilities.hpp"

#include <EbsdLib/LaueOps/LaueOps.h>

#include <vector>

using namespace nx::core;

// -----------------------------------------------------------------------------
AlignSectionsMutualInformation::AlignSectionsMutualInformation(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                                                               AlignSectionsMutualInformationInputValues* inputValues)
: AlignSections(dataStructure, shouldCancel, messageHandler)
, m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(messageHandler)
{
}

// -----------------------------------------------------------------------------
AlignSectionsMutualInformation::~AlignSectionsMutualInformation() noexcept = default;

// -----------------------------------------------------------------------------
Result<> AlignSectionsMutualInformation::operator()()
{
  if(m_ShouldCancel)
  {
    return {};
  }
  const auto& gridGeom = m_DataStructure.getDataRefAs<IGridGeometry>(m_InputValues->ImageGeometryPath);

  return execute(gridGeom.getDimensions(), m_InputValues->ImageGeometryPath);
}

// -----------------------------------------------------------------------------
Result<> AlignSectionsMutualInformation::findShifts(std::vector<int64>& xShifts, std::vector<int64>& yShifts)
{
  const auto& imageGeom = m_DataStructure.getDataRefAs<ImageGeom>(m_InputValues->ImageGeometryPath);
  const AttributeMatrix* cellData = imageGeom.getCellData();
  auto totalPoints = static_cast<int64>(cellData->getNumberOfTuples());

  if(m_InputValues->UseMask)
  {
    try
    {
      m_MaskCompare = MaskCompareUtilities::InstantiateMaskCompare(m_DataStructure, m_InputValues->MaskArrayPath);
    } catch(const std::out_of_range& exception)
    {
      // This really should NOT be happening as the path was verified during preflight BUT we may be calling this from
      // somewhere else that is NOT going through the normal nx::core::IFilter API of Preflight and Execute
      std::string message = fmt::format("Mask Array DataPath does not exist or is not of the correct type (Bool | UInt8) {}", m_InputValues->MaskArrayPath.toString());
      return MakeErrorResult(-53702, message);
    }
  }

  SizeVec3 udims = imageGeom.getDimensions();
  int64 dims[3] = {
      static_cast<int64>(udims[0]),
      static_cast<int64>(udims[1]),
      static_cast<int64>(udims[2]),
  };

  std::vector<int32> sliceFeatureIds(totalPoints, 0);
  std::vector<int32> sliceFeatureCounts(dims[2], 0);

  std::vector<std::vector<float32>> mutualInfo12;
  std::vector<float32> mutualInfo1;
  std::vector<float32> mutualInfo2;

  // Segment each slice
  Result<> formFeaturesResult = formFeaturesSections(sliceFeatureIds, sliceFeatureCounts);
  if(formFeaturesResult.invalid())
  {
    return formFeaturesResult;
  }

  std::vector<std::vector<float32>> misorientations(dims[0]);
  for(int64 i = 0; i < dims[0]; i++)
  {
    misorientations[i].assign(dims[1], 0.0f);
  }

  if(m_InputValues->StoreAlignmentShifts)
  {
    auto& slicesStore = m_DataStructure.getDataAs<UInt32Array>(m_InputValues->SlicesArrayPath)->getDataStoreRef();
    auto& relativeShiftsStore = m_DataStructure.getDataAs<Int64Array>(m_InputValues->RelativeShiftsArrayPath)->getDataStoreRef();
    auto& cumulativeShiftsStore = m_DataStructure.getDataAs<Int64Array>(m_InputValues->CumulativeShiftsArrayPath)->getDataStoreRef();
    for(int64 iter = 1; iter < dims[2]; iter++)
    {
      if(m_ShouldCancel)
      {
        return {};
      }
      m_MessageHandler(IFilter::Message::Type::Info, fmt::format("Determining Shifts: Slice {}/{} complete", iter, dims[2]));

      float32 minDisorientation = std::numeric_limits<float32>::max();
      int64 slice = (dims[2] - 1) - iter;
      int32 featureCount1 = sliceFeatureCounts[slice];
      int32 featureCount2 = sliceFeatureCounts[slice + 1];
      mutualInfo12 = std::vector<std::vector<float32>>(featureCount1, std::vector<float32>(featureCount2, 0.0f));
      mutualInfo1 = std::vector<float32>(featureCount1, 0.0f);
      mutualInfo2 = std::vector<float32>(featureCount2, 0.0f);

      int64 oldXShift = -1;
      int64 oldYShift = -1;
      int64 newXShift = 0;
      int64 newYShift = 0;
      for(int64 i = 0; i < dims[0]; i++)
      {
        for(int64 j = 0; j < dims[1]; j++)
        {
          misorientations[i][j] = 0.0F;
        }
      }
      while(newXShift != oldXShift || newYShift != oldYShift)
      {
        oldXShift = newXShift;
        oldYShift = newYShift;
        for(int32 j = -3; j < 4; j++)
        {
          for(int32 k = -3; k < 4; k++)
          {
            float32 disorientation = 0.0F;
            float32 count = 0.0F;
            if(misorientations[k + oldXShift + dims[0] / 2][j + oldYShift + dims[1] / 2] == 0 && llabs(k + oldXShift) < (dims[0] / 2) && (j + oldYShift) < (dims[1] / 2))
            {
              for(int64 dim1Index = 0; dim1Index < dims[1]; dim1Index = dim1Index + 4)
              {
                for(int64 dim0Index = 0; dim0Index < dims[0]; dim0Index = dim0Index + 4)
                {
                  if((dim1Index + j + oldYShift) >= 0 && (dim1Index + j + oldYShift) < dims[1] && (dim0Index + k + oldXShift) >= 0 && (dim0Index + k + oldXShift) < dims[0])
                  {
                    int64 refPosition = ((slice + 1) * dims[0] * dims[1]) + (dim1Index * dims[0]) + dim0Index;
                    int64 curPosition = (slice * dims[0] * dims[1]) + ((dim1Index + j + oldYShift) * dims[0]) + (dim0Index + k + oldXShift);
                    int32 refGNum = sliceFeatureIds[refPosition];
                    int32 curGNum = sliceFeatureIds[curPosition];
                    if(curGNum >= 0 && refGNum >= 0)
                    {
                      mutualInfo12[curGNum][refGNum]++;
                      mutualInfo1[curGNum]++;
                      mutualInfo2[refGNum]++;
                      count++;
                    }
                  }
                  else
                  {
                    mutualInfo12[0][0]++;
                    mutualInfo1[0]++;
                    mutualInfo2[0]++;
                  }
                }
              }
              for(int32 featureCount1Index = 0; featureCount1Index < featureCount1; featureCount1Index++)
              {
                mutualInfo1[featureCount1Index] = mutualInfo1[featureCount1Index] / count;
              }
              for(int32 featureCount2Index = 0; featureCount2Index < featureCount2; featureCount2Index++)
              {
                mutualInfo2[featureCount2Index] = mutualInfo2[featureCount2Index] / static_cast<float32>(count);
              }
              for(int32 featureCount1Index = 0; featureCount1Index < featureCount1; featureCount1Index++)
              {
                for(int32 featureCount2Index = 0; featureCount2Index < featureCount2; featureCount2Index++)
                {
                  mutualInfo12[featureCount1Index][featureCount2Index] = mutualInfo12[featureCount1Index][featureCount2Index] / count;

                  float32 value = 0.0f;
                  if(mutualInfo1[featureCount1Index] > 0 && mutualInfo2[featureCount2Index] > 0)
                  {
                    value = (mutualInfo12[featureCount1Index][featureCount2Index] / (mutualInfo1[featureCount1Index] * mutualInfo2[featureCount2Index]));
                  }
                  if(value != 0)
                  {
                    disorientation = disorientation + (mutualInfo12[featureCount1Index][featureCount2Index] * logf(value));
                  }
                }
              }
              for(int32 featureCount1Index = 0; featureCount1Index < featureCount1; featureCount1Index++)
              {
                for(int32 featureCount2Index = 0; featureCount2Index < featureCount2; featureCount2Index++)
                {
                  mutualInfo12[featureCount1Index][featureCount2Index] = 0.0f;
                  mutualInfo1[featureCount1Index] = 0.0f;
                  mutualInfo2[featureCount2Index] = 0.0f;
                }
              }
              disorientation = 1.0f / disorientation;
              misorientations[k + oldXShift + dims[0] / 2][j + oldYShift + dims[1] / 2] = disorientation;
              if(disorientation < minDisorientation)
              {
                newXShift = k + oldXShift;
                newYShift = j + oldYShift;
                minDisorientation = disorientation;
              }
            }
          }
        }
      }
      xShifts[iter] = xShifts[iter - 1] + newXShift;
      yShifts[iter] = yShifts[iter - 1] + newYShift;

      usize xIndex = iter * 2;
      usize yIndex = (iter * 2) + 1;
      slicesStore[xIndex] = slice;
      slicesStore[yIndex] = slice + 1;
      relativeShiftsStore[xIndex] = newXShift;
      relativeShiftsStore[yIndex] = newYShift;
      cumulativeShiftsStore[xIndex] = xShifts[iter];
      cumulativeShiftsStore[yIndex] = yShifts[iter];
    }
  }
  else
  {
    for(int64 iter = 1; iter < dims[2]; iter++)
    {
      m_MessageHandler(IFilter::Message::Type::Info, fmt::format("Determining Shifts: Slice {}/{} complete", iter, dims[2]));

      float32 minDisorientation = std::numeric_limits<float32>::max();
      int64 slice = (dims[2] - 1) - iter;
      int32 featureCount1 = sliceFeatureCounts[slice];
      int32 featureCount2 = sliceFeatureCounts[slice + 1];
      mutualInfo12 = std::vector<std::vector<float32>>(featureCount1, std::vector<float32>(featureCount2, 0.0f));
      mutualInfo1 = std::vector<float32>(featureCount1, 0.0f);
      mutualInfo2 = std::vector<float32>(featureCount2, 0.0f);

      int64 oldXShift = -1;
      int64 oldYShift = -1;
      int64 newXShift = 0;
      int64 newYShift = 0;
      for(int64 i = 0; i < dims[0]; i++)
      {
        for(int64 j = 0; j < dims[1]; j++)
        {
          misorientations[i][j] = 0.0F;
        }
      }
      while(newXShift != oldXShift || newYShift != oldYShift)
      {
        oldXShift = newXShift;
        oldYShift = newYShift;
        for(int32 j = -3; j < 4; j++)
        {
          for(int32 k = -3; k < 4; k++)
          {
            float32 disorientation = 0.0F;
            float32 count = 0.0F;
            if(misorientations[k + oldXShift + dims[0] / 2][j + oldYShift + dims[1] / 2] == 0 && llabs(k + oldXShift) < (dims[0] / 2) && (j + oldYShift) < (dims[1] / 2))
            {
              for(int64 dim1Index = 0; dim1Index < dims[1]; dim1Index = dim1Index + 4)
              {
                for(int64 dim0Index = 0; dim0Index < dims[0]; dim0Index = dim0Index + 4)
                {
                  if((dim1Index + j + oldYShift) >= 0 && (dim1Index + j + oldYShift) < dims[1] && (dim0Index + k + oldXShift) >= 0 && (dim0Index + k + oldXShift) < dims[0])
                  {
                    int64 refPosition = ((slice + 1) * dims[0] * dims[1]) + (dim1Index * dims[0]) + dim0Index;
                    int64 curPosition = (slice * dims[0] * dims[1]) + ((dim1Index + j + oldYShift) * dims[0]) + (dim0Index + k + oldXShift);
                    int32 refGNum = sliceFeatureIds[refPosition];
                    int32 curGNum = sliceFeatureIds[curPosition];
                    if(curGNum >= 0 && refGNum >= 0)
                    {
                      mutualInfo12[curGNum][refGNum]++;
                      mutualInfo1[curGNum]++;
                      mutualInfo2[refGNum]++;
                      count++;
                    }
                  }
                  else
                  {
                    mutualInfo12[0][0]++;
                    mutualInfo1[0]++;
                    mutualInfo2[0]++;
                  }
                }
              }
              for(int32 featureCount1Index = 0; featureCount1Index < featureCount1; featureCount1Index++)
              {
                mutualInfo1[featureCount1Index] = mutualInfo1[featureCount1Index] / count;
              }
              for(int32 featureCount2Index = 0; featureCount2Index < featureCount2; featureCount2Index++)
              {
                mutualInfo2[featureCount2Index] = mutualInfo2[featureCount2Index] / static_cast<float32>(count);
              }
              for(int32 featureCount1Index = 0; featureCount1Index < featureCount1; featureCount1Index++)
              {
                for(int32 featureCount2Index = 0; featureCount2Index < featureCount2; featureCount2Index++)
                {
                  mutualInfo12[featureCount1Index][featureCount2Index] = mutualInfo12[featureCount1Index][featureCount2Index] / count;

                  float32 value = 0.0f;
                  if(mutualInfo1[featureCount1Index] > 0 && mutualInfo2[featureCount2Index] > 0)
                  {
                    value = (mutualInfo12[featureCount1Index][featureCount2Index] / (mutualInfo1[featureCount1Index] * mutualInfo2[featureCount2Index]));
                  }
                  if(value != 0)
                  {
                    disorientation = disorientation + (mutualInfo12[featureCount1Index][featureCount2Index] * logf(value));
                  }
                }
              }
              for(int32 featureCount1Index = 0; featureCount1Index < featureCount1; featureCount1Index++)
              {
                for(int32 featureCount2Index = 0; featureCount2Index < featureCount2; featureCount2Index++)
                {
                  mutualInfo12[featureCount1Index][featureCount2Index] = 0.0f;
                  mutualInfo1[featureCount1Index] = 0.0f;
                  mutualInfo2[featureCount2Index] = 0.0f;
                }
              }
              disorientation = 1.0f / disorientation;
              misorientations[k + oldXShift + dims[0] / 2][j + oldYShift + dims[1] / 2] = disorientation;
              if(disorientation < minDisorientation)
              {
                newXShift = k + oldXShift;
                newYShift = j + oldYShift;
                minDisorientation = disorientation;
              }
            }
          }
        }
      }
      xShifts[iter] = xShifts[iter - 1] + newXShift;
      yShifts[iter] = yShifts[iter - 1] + newYShift;
    }
  }

  return {};
}

// -----------------------------------------------------------------------------
Result<> AlignSectionsMutualInformation::formFeaturesSections(std::vector<int32>& sliceFeatureIds, std::vector<int32>& sliceFeatureCounts)
{
  const auto& imageGeom = m_DataStructure.getDataRefAs<ImageGeom>(m_InputValues->ImageGeometryPath);

  const SizeVec3 gridDimensions = imageGeom.getDimensions();
  const std::array<int64, 3> dimensions = {
      static_cast<int64>(gridDimensions[0]),
      static_cast<int64>(gridDimensions[1]),
      static_cast<int64>(gridDimensions[2]),
  };

  const auto orientationOps = ebsdlib::LaueOps::GetAllOrientationOps();

  const auto& quatsArrayRef = m_DataStructure.getDataRefAs<Float32Array>(m_InputValues->QuatsArrayPath);
  const auto& cellPhasesArrayRef = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->CellPhasesArrayPath);
  const auto& crystalStructuresArrayRef = m_DataStructure.getDataRefAs<UInt32Array>(m_InputValues->CrystalStructuresArrayPath);
  const usize numCrystalStructures = crystalStructuresArrayRef.getNumberOfTuples();

  constexpr usize initialVoxelListSize = 1000;

  const float32 misorientationTolerance = m_InputValues->MisorientationTolerance * nx::core::Constants::k_PiOver180F;

  sliceFeatureCounts.resize(dimensions[2]);

  std::vector<int64> voxelList(initialVoxelListSize, -1);
  const std::array<int64, 4> neighborOffsets = {-dimensions[0], -1, 1, dimensions[0]};

  for(int64 sliceIdx = 0; sliceIdx < dimensions[2]; sliceIdx++)
  {
    m_MessageHandler(IFilter::Message::Type::Info, fmt::format("Identifying Features: Slice {}/{} complete", sliceIdx, dimensions[2]));

    const int64 sliceStartVoxelIdx = sliceIdx * dimensions[0] * dimensions[1];
    const int64 sliceEndVoxelIdx = (sliceIdx + 1) * dimensions[0] * dimensions[1];
    int64 nextSeedSearchVoxelIdx = sliceStartVoxelIdx;

    int32 featureCount = 1;
    bool hasNoSeeds = false;
    while(!hasNoSeeds)
    {
      int64 seedVoxelIdx = -1;

      for(int64 voxelIdx = nextSeedSearchVoxelIdx; voxelIdx < sliceEndVoxelIdx; voxelIdx++)
      {
        if((!m_InputValues->UseMask || (m_MaskCompare != nullptr && m_MaskCompare->isTrue(voxelIdx))) && sliceFeatureIds[voxelIdx] == 0 && cellPhasesArrayRef[voxelIdx] > 0)
        {
          seedVoxelIdx = voxelIdx;
          nextSeedSearchVoxelIdx = voxelIdx;
        }
        if(seedVoxelIdx > -1)
        {
          break;
        }
      }

      if(seedVoxelIdx == -1)
      {
        hasNoSeeds = true;
      }
      if(seedVoxelIdx >= 0)
      {
        std::vector<int64>::size_type voxelListSize = 0;
        sliceFeatureIds[seedVoxelIdx] = featureCount;
        voxelList[voxelListSize] = seedVoxelIdx;
        voxelListSize++;
        for(usize voxelListIdx = 0; voxelListIdx < voxelListSize; ++voxelListIdx)
        {
          const int64 currentVoxelIdx = voxelList[voxelListIdx];
          const int64 xIdx = currentVoxelIdx % dimensions[0];
          const int64 yIdx = (currentVoxelIdx / dimensions[0]) % dimensions[1];

          const usize currentQuatOffset = currentVoxelIdx * 4;
          const ebsdlib::QuatD currentQuat(quatsArrayRef[currentQuatOffset], quatsArrayRef[currentQuatOffset + 1], quatsArrayRef[currentQuatOffset + 2], quatsArrayRef[currentQuatOffset + 3]);
          const int32 currentCellPhaseIdx = cellPhasesArrayRef[currentVoxelIdx];
          if(static_cast<usize>(currentCellPhaseIdx) >= numCrystalStructures)
          {
            return MakeErrorResult(-53703, fmt::format("Cell Phases array '{}' has value {} at voxel index {}, but Crystal Structures array '{}' has {} tuples. Valid Phase indices are in [0, {}).",
                                                       m_InputValues->CellPhasesArrayPath.toString(), currentCellPhaseIdx, currentVoxelIdx, m_InputValues->CrystalStructuresArrayPath.toString(),
                                                       numCrystalStructures, numCrystalStructures));
          }
          const uint32 currentLaueIndex = crystalStructuresArrayRef[currentCellPhaseIdx];
          if(currentLaueIndex >= orientationOps.size())
          {
            return MakeErrorResult(-53704, fmt::format("Crystal Structures array '{}' has value {} at Phase index {}, but only {} Laue operations are available. Valid Laue indices are in [0, {}).",
                                                       m_InputValues->CrystalStructuresArrayPath.toString(), currentLaueIndex, currentCellPhaseIdx, orientationOps.size(), orientationOps.size()));
          }
          for(usize faceIdx = 0; faceIdx < neighborOffsets.size(); faceIdx++)
          {
            const int64 neighborVoxelIdx = currentVoxelIdx + neighborOffsets[faceIdx];
            if((faceIdx == 0) && yIdx == 0)
            {
              continue;
            }
            if((faceIdx == 3) && yIdx == (dimensions[1] - 1))
            {
              continue;
            }
            if((faceIdx == 1) && xIdx == 0)
            {
              continue;
            }
            if((faceIdx == 2) && xIdx == (dimensions[0] - 1))
            {
              continue;
            }
            const int32 neighborCellPhaseIdx = cellPhasesArrayRef[neighborVoxelIdx];
            if(sliceFeatureIds[neighborVoxelIdx] <= 0 && neighborCellPhaseIdx > 0)
            {
              if(static_cast<usize>(neighborCellPhaseIdx) >= numCrystalStructures)
              {
                return MakeErrorResult(-53703,
                                       fmt::format("Cell Phases array '{}' has value {} at voxel index {}, but Crystal Structures array '{}' has {} tuples. Valid Phase indices are in [0, {}).",
                                                   m_InputValues->CellPhasesArrayPath.toString(), neighborCellPhaseIdx, neighborVoxelIdx, m_InputValues->CrystalStructuresArrayPath.toString(),
                                                   numCrystalStructures, numCrystalStructures));
              }
              float32 misorientationAngle = std::numeric_limits<float>::max();
              const usize neighborQuatOffset = neighborVoxelIdx * 4;
              const ebsdlib::QuatD neighborQuat(quatsArrayRef[neighborQuatOffset], quatsArrayRef[neighborQuatOffset + 1], quatsArrayRef[neighborQuatOffset + 2], quatsArrayRef[neighborQuatOffset + 3]);
              const uint32 neighborLaueIndex = crystalStructuresArrayRef[neighborCellPhaseIdx];
              if(neighborLaueIndex >= orientationOps.size())
              {
                return MakeErrorResult(-53704,
                                       fmt::format("Crystal Structures array '{}' has value {} at Phase index {}, but only {} Laue operations are available. Valid Laue indices are in [0, {}).",
                                                   m_InputValues->CrystalStructuresArrayPath.toString(), neighborLaueIndex, neighborCellPhaseIdx, orientationOps.size(), orientationOps.size()));
              }

              if(currentLaueIndex == neighborLaueIndex)
              {
                const ebsdlib::AxisAngleDType axisAngle = orientationOps[currentLaueIndex]->calculateMisorientation(currentQuat, neighborQuat);
                misorientationAngle = axisAngle[3];
              }
              if(misorientationAngle < misorientationTolerance)
              {
                sliceFeatureIds[neighborVoxelIdx] = featureCount;
                voxelList[voxelListSize] = neighborVoxelIdx;
                voxelListSize++;
                if(voxelListSize >= voxelList.size())
                {
                  voxelListSize = voxelList.size();
                  voxelList.resize(voxelListSize + initialVoxelListSize);
                  for(usize resetIdx = voxelListSize; resetIdx < voxelList.size(); ++resetIdx)
                  {
                    voxelList[resetIdx] = -1;
                  }
                }
              }
            }
          }
        }
        voxelList.erase(std::remove(voxelList.begin(), voxelList.end(), -1), voxelList.end());
        featureCount++;
        voxelList.assign(initialVoxelListSize, -1);
      }
    }
    sliceFeatureCounts[sliceIdx] = featureCount;
  }
  return {};
}
