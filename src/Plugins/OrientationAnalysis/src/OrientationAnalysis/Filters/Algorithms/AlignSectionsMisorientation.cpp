#include "AlignSectionsMisorientation.hpp"

#include "simplnx/Common/Numbers.hpp"
#include "simplnx/DataStructure/DataGroup.hpp"
#include "simplnx/DataStructure/Geometry/IGridGeometry.hpp"
#include "simplnx/Utilities/FilterUtilities.hpp"
#include "simplnx/Utilities/MaskCompareUtilities.hpp"

#include <EbsdLib/LaueOps/LaueOps.h>

#include <iostream>

using namespace nx::core;

// -----------------------------------------------------------------------------
AlignSectionsMisorientation::AlignSectionsMisorientation(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                                                         AlignSectionsMisorientationInputValues* inputValues)
: AlignSections(dataStructure, shouldCancel, messageHandler)
, m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(messageHandler)
{
}

// -----------------------------------------------------------------------------
AlignSectionsMisorientation::~AlignSectionsMisorientation() noexcept = default;

// -----------------------------------------------------------------------------
Result<> AlignSectionsMisorientation::operator()()
{
  if(m_ShouldCancel)
  {
    return {};
  }
  const auto& gridGeom = m_DataStructure.getDataRefAs<IGridGeometry>(m_InputValues->ImageGeometryPath);

  return execute(gridGeom.getDimensions(), m_InputValues->ImageGeometryPath);
}

// -----------------------------------------------------------------------------
Result<> AlignSectionsMisorientation::findShifts(std::vector<int64_t>& xShifts, std::vector<int64_t>& yShifts)
{
  std::unique_ptr<MaskCompareUtilities::MaskCompare> maskComparePtr = nullptr;
  if(m_InputValues->UseMask)
  {
    try
    {
      maskComparePtr = MaskCompareUtilities::InstantiateMaskCompare(m_DataStructure, m_InputValues->MaskArrayPath);
    } catch(const std::out_of_range& exception)
    {
      // This really should NOT be happening as the path was verified during preflight BUT we may be calling this from
      // somewhere else that is NOT going through the normal nx::core::IFilter API of Preflight and Execute
      std::string message = fmt::format("Mask Array DataPath does not exist or is not of the correct type (Bool | UInt8) {}", m_InputValues->MaskArrayPath.toString());
      return MakeErrorResult(-53900, message);
    }
  }

  const auto& gridGeom = m_DataStructure.getDataRefAs<IGridGeometry>(m_InputValues->ImageGeometryPath);

  const auto& cellPhasesArrayRef = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->CellPhasesArrayPath);
  const auto& quatsArrayRef = m_DataStructure.getDataRefAs<Float32Array>(m_InputValues->QuatsArrayPath);
  const auto& crystalStructuresArrayRef = m_DataStructure.getDataRefAs<UInt32Array>(m_InputValues->CrystalStructuresArrayPath);
  const usize numCrystalStructures = crystalStructuresArrayRef.getNumberOfTuples();

  const SizeVec3 gridDimensions = gridGeom.getDimensions();

  const std::array<int64_t, 3> dimensions = {
      static_cast<int64_t>(gridDimensions[0]),
      static_cast<int64_t>(gridDimensions[1]),
      static_cast<int64_t>(gridDimensions[2]),
  };

  std::vector<ebsdlib::LaueOps::Pointer> orientationOps = ebsdlib::LaueOps::GetAllOrientationOps();

  // Allocate a 2D Array which will be reused from slice to slice
  std::vector<bool> evaluatedShifts(dimensions[0] * dimensions[1], false);

  const auto halfXDimension = static_cast<int64_t>(dimensions[0] * 0.5f);
  const auto halfYDimension = static_cast<int64_t>(dimensions[1] * 0.5f);

  constexpr double degreesToRadians = nx::core::numbers::pi / 180.0;
  ThrottledMessenger throttledMessenger = getMessageHelper().createThrottledMessenger();
  if(m_InputValues->StoreAlignmentShifts)
  {
    auto& slicesStoreRef = m_DataStructure.getDataAs<UInt32Array>(m_InputValues->SlicesArrayPath)->getDataStoreRef();
    auto& relativeShiftsStoreRef = m_DataStructure.getDataAs<Int64Array>(m_InputValues->RelativeShiftsArrayPath)->getDataStoreRef();
    auto& cumulativeShiftsStoreRef = m_DataStructure.getDataAs<Int64Array>(m_InputValues->CumulativeShiftsArrayPath)->getDataStoreRef();
    // Process adjacent sections from the largest Z index to the smallest Z index.
    for(int64_t shiftIdx = 1; shiftIdx < dimensions[2]; shiftIdx++)
    {
      if(m_ShouldCancel)
      {
        return {};
      }
      throttledMessenger.sendThrottledMessage([&]() { return fmt::format("Determining Shifts || {:.2f}% Complete", CalculatePercentComplete(shiftIdx, dimensions[2])); });
      if(getCancel())
      {
        return {};
      }
      float minimumDisorientation = std::numeric_limits<float>::max();
      const int64 sliceIdx = (dimensions[2] - 1) - shiftIdx;
      int64 previousXShift = -1;
      int64 previousYShift = -1;
      int64 currentXShift = 0;
      int64 currentYShift = 0;

      std::fill(evaluatedShifts.begin(), evaluatedShifts.end(), false);

      const float misorientationTolerance = static_cast<float>(m_InputValues->MisorientationTolerance * degreesToRadians);

      while(currentXShift != previousXShift || currentYShift != previousYShift)
      {
        previousXShift = currentXShift;
        previousYShift = currentYShift;
        for(int32 yShiftOffset = -3; yShiftOffset < 4; yShiftOffset++)
        {
          for(int32 xShiftOffset = -3; xShiftOffset < 4; xShiftOffset++)
          {
            float disorientationScore = 0.0F;
            float sampleCount = 0.0F;
            int64 evaluatedShiftXIdx = xShiftOffset + previousXShift + halfXDimension;
            int64 evaluatedShiftYIdx = yShiftOffset + previousYShift + halfYDimension;
            int64 evaluatedShiftIdx = (dimensions[0] * evaluatedShiftYIdx) + evaluatedShiftXIdx;
            if(!evaluatedShifts[evaluatedShiftIdx] && llabs(xShiftOffset + previousXShift) < halfXDimension && llabs(yShiftOffset + previousYShift) < halfYDimension)
            {
              for(int64 sampledYIdx = 0; sampledYIdx < dimensions[1]; sampledYIdx += 4)
              {
                for(int64 sampledXIdx = 0; sampledXIdx < dimensions[0]; sampledXIdx += 4)
                {
                  if((sampledYIdx + yShiftOffset + previousYShift) >= 0 && (sampledYIdx + yShiftOffset + previousYShift) < dimensions[1] && (sampledXIdx + xShiftOffset + previousXShift) >= 0 &&
                     (sampledXIdx + xShiftOffset + previousXShift) < dimensions[0])
                  {
                    sampleCount++;
                    const int64 referenceVoxelIdx = ((sliceIdx + 1) * dimensions[0] * dimensions[1]) + (sampledYIdx * dimensions[0]) + sampledXIdx;
                    const int64 currentVoxelIdx =
                        (sliceIdx * dimensions[0] * dimensions[1]) + ((sampledYIdx + yShiftOffset + previousYShift) * dimensions[0]) + (sampledXIdx + xShiftOffset + previousXShift);
                    if(!m_InputValues->UseMask || maskComparePtr->bothTrue(referenceVoxelIdx, currentVoxelIdx))
                    {
                      float misorientationAngle = std::numeric_limits<float>::max();
                      const int32 referenceCellPhaseIdx = cellPhasesArrayRef[referenceVoxelIdx];
                      const int32 currentCellPhaseIdx = cellPhasesArrayRef[currentVoxelIdx];
                      if(referenceCellPhaseIdx > 0 && currentCellPhaseIdx > 0)
                      {
                        if(static_cast<usize>(referenceCellPhaseIdx) >= numCrystalStructures)
                        {
                          return MakeErrorResult(
                              -53901, fmt::format("Cell Phases array '{}' has value {} at voxel index {}, but Crystal Structures array '{}' has {} tuples. Valid Phase indices are in [0, {}).",
                                                  m_InputValues->CellPhasesArrayPath.toString(), referenceCellPhaseIdx, referenceVoxelIdx, m_InputValues->CrystalStructuresArrayPath.toString(),
                                                  numCrystalStructures, numCrystalStructures));
                        }
                        if(static_cast<usize>(currentCellPhaseIdx) >= numCrystalStructures)
                        {
                          return MakeErrorResult(
                              -53901, fmt::format("Cell Phases array '{}' has value {} at voxel index {}, but Crystal Structures array '{}' has {} tuples. Valid Phase indices are in [0, {}).",
                                                  m_InputValues->CellPhasesArrayPath.toString(), currentCellPhaseIdx, currentVoxelIdx, m_InputValues->CrystalStructuresArrayPath.toString(),
                                                  numCrystalStructures, numCrystalStructures));
                        }
                        const ebsdlib::QuatD referenceQuat(quatsArrayRef[referenceVoxelIdx * 4], quatsArrayRef[referenceVoxelIdx * 4 + 1], quatsArrayRef[referenceVoxelIdx * 4 + 2],
                                                           quatsArrayRef[referenceVoxelIdx * 4 + 3]);
                        const uint32 referenceLaueIndex = crystalStructuresArrayRef[referenceCellPhaseIdx];
                        const ebsdlib::QuatD currentQuat(quatsArrayRef[currentVoxelIdx * 4], quatsArrayRef[currentVoxelIdx * 4 + 1], quatsArrayRef[currentVoxelIdx * 4 + 2],
                                                         quatsArrayRef[currentVoxelIdx * 4 + 3]);
                        const uint32 currentLaueIndex = crystalStructuresArrayRef[currentCellPhaseIdx];
                        if(referenceLaueIndex == currentLaueIndex && referenceLaueIndex < orientationOps.size())
                        {
                          const ebsdlib::AxisAngleDType axisAngle = orientationOps[referenceLaueIndex]->calculateMisorientation(referenceQuat, currentQuat);
                          misorientationAngle = axisAngle[3];
                        }
                      }
                      if(misorientationAngle > misorientationTolerance)
                      {
                        disorientationScore++;
                      }
                    }
                    if(m_InputValues->UseMask)
                    {
                      if(maskComparePtr->isTrue(referenceVoxelIdx) && !maskComparePtr->isTrue(currentVoxelIdx))
                      {
                        disorientationScore++;
                      }
                      if(!maskComparePtr->isTrue(referenceVoxelIdx) && maskComparePtr->isTrue(currentVoxelIdx))
                      {
                        disorientationScore++;
                      }
                    }
                  }
                }
              }
              disorientationScore /= sampleCount;
              evaluatedShiftXIdx = xShiftOffset + previousXShift + halfXDimension;
              evaluatedShiftYIdx = yShiftOffset + previousYShift + halfYDimension;
              evaluatedShiftIdx = (dimensions[0] * evaluatedShiftYIdx) + evaluatedShiftXIdx;
              evaluatedShifts[evaluatedShiftIdx] = true;
              if(disorientationScore < minimumDisorientation ||
                 (disorientationScore == minimumDisorientation && ((llabs(xShiftOffset + previousXShift) < llabs(currentXShift)) || (llabs(yShiftOffset + previousYShift) < llabs(currentYShift)))))
              {
                currentXShift = xShiftOffset + previousXShift;
                currentYShift = yShiftOffset + previousYShift;
                minimumDisorientation = disorientationScore;
              }
            }
          }
        }
      }
      xShifts[shiftIdx] = xShifts[shiftIdx - 1] + currentXShift;
      yShifts[shiftIdx] = yShifts[shiftIdx - 1] + currentYShift;
      const usize xShiftValueIdx = shiftIdx * 2;
      const usize yShiftValueIdx = (shiftIdx * 2) + 1;
      slicesStoreRef[xShiftValueIdx] = sliceIdx;
      slicesStoreRef[yShiftValueIdx] = sliceIdx + 1;
      relativeShiftsStoreRef[xShiftValueIdx] = currentXShift;
      relativeShiftsStoreRef[yShiftValueIdx] = currentYShift;
      cumulativeShiftsStoreRef[xShiftValueIdx] = xShifts[shiftIdx];
      cumulativeShiftsStoreRef[yShiftValueIdx] = yShifts[shiftIdx];
    }
  }
  else
  {
    // Process adjacent sections from the largest Z index to the smallest Z index.
    for(int64_t shiftIdx = 1; shiftIdx < dimensions[2]; shiftIdx++)
    {
      throttledMessenger.sendThrottledMessage([&]() { return fmt::format("Determining Shifts || {:.2f}% Complete", CalculatePercentComplete(shiftIdx, dimensions[2])); });
      if(getCancel())
      {
        return {};
      }
      float minimumDisorientation = std::numeric_limits<float>::max();
      const int64 sliceIdx = (dimensions[2] - 1) - shiftIdx;
      int64 previousXShift = -1;
      int64 previousYShift = -1;
      int64 currentXShift = 0;
      int64 currentYShift = 0;

      std::fill(evaluatedShifts.begin(), evaluatedShifts.end(), false);

      const float misorientationTolerance = static_cast<float>(m_InputValues->MisorientationTolerance * degreesToRadians);

      while(currentXShift != previousXShift || currentYShift != previousYShift)
      {
        previousXShift = currentXShift;
        previousYShift = currentYShift;
        for(int32 yShiftOffset = -3; yShiftOffset < 4; yShiftOffset++)
        {
          for(int32 xShiftOffset = -3; xShiftOffset < 4; xShiftOffset++)
          {
            float disorientationScore = 0.0F;
            float sampleCount = 0.0F;
            int64 evaluatedShiftXIdx = xShiftOffset + previousXShift + halfXDimension;
            int64 evaluatedShiftYIdx = yShiftOffset + previousYShift + halfYDimension;
            int64 evaluatedShiftIdx = (dimensions[0] * evaluatedShiftYIdx) + evaluatedShiftXIdx;
            if(!evaluatedShifts[evaluatedShiftIdx] && llabs(xShiftOffset + previousXShift) < halfXDimension && llabs(yShiftOffset + previousYShift) < halfYDimension)
            {
              for(int64 sampledYIdx = 0; sampledYIdx < dimensions[1]; sampledYIdx += 4)
              {
                for(int64 sampledXIdx = 0; sampledXIdx < dimensions[0]; sampledXIdx += 4)
                {
                  if((sampledYIdx + yShiftOffset + previousYShift) >= 0 && (sampledYIdx + yShiftOffset + previousYShift) < dimensions[1] && (sampledXIdx + xShiftOffset + previousXShift) >= 0 &&
                     (sampledXIdx + xShiftOffset + previousXShift) < dimensions[0])
                  {
                    sampleCount++;
                    const int64 referenceVoxelIdx = ((sliceIdx + 1) * dimensions[0] * dimensions[1]) + (sampledYIdx * dimensions[0]) + sampledXIdx;
                    const int64 currentVoxelIdx =
                        (sliceIdx * dimensions[0] * dimensions[1]) + ((sampledYIdx + yShiftOffset + previousYShift) * dimensions[0]) + (sampledXIdx + xShiftOffset + previousXShift);
                    if(!m_InputValues->UseMask || maskComparePtr->bothTrue(referenceVoxelIdx, currentVoxelIdx))
                    {
                      float misorientationAngle = std::numeric_limits<float>::max();
                      const int32 referenceCellPhaseIdx = cellPhasesArrayRef[referenceVoxelIdx];
                      const int32 currentCellPhaseIdx = cellPhasesArrayRef[currentVoxelIdx];
                      if(referenceCellPhaseIdx > 0 && currentCellPhaseIdx > 0)
                      {
                        if(static_cast<usize>(referenceCellPhaseIdx) >= numCrystalStructures)
                        {
                          return MakeErrorResult(
                              -53901, fmt::format("Cell Phases array '{}' has value {} at voxel index {}, but Crystal Structures array '{}' has {} tuples. Valid Phase indices are in [0, {}).",
                                                  m_InputValues->CellPhasesArrayPath.toString(), referenceCellPhaseIdx, referenceVoxelIdx, m_InputValues->CrystalStructuresArrayPath.toString(),
                                                  numCrystalStructures, numCrystalStructures));
                        }
                        if(static_cast<usize>(currentCellPhaseIdx) >= numCrystalStructures)
                        {
                          return MakeErrorResult(
                              -53901, fmt::format("Cell Phases array '{}' has value {} at voxel index {}, but Crystal Structures array '{}' has {} tuples. Valid Phase indices are in [0, {}).",
                                                  m_InputValues->CellPhasesArrayPath.toString(), currentCellPhaseIdx, currentVoxelIdx, m_InputValues->CrystalStructuresArrayPath.toString(),
                                                  numCrystalStructures, numCrystalStructures));
                        }
                        const ebsdlib::QuatD referenceQuat(quatsArrayRef[referenceVoxelIdx * 4], quatsArrayRef[referenceVoxelIdx * 4 + 1], quatsArrayRef[referenceVoxelIdx * 4 + 2],
                                                           quatsArrayRef[referenceVoxelIdx * 4 + 3]);
                        const uint32 referenceLaueIndex = crystalStructuresArrayRef[referenceCellPhaseIdx];
                        const ebsdlib::QuatD currentQuat(quatsArrayRef[currentVoxelIdx * 4], quatsArrayRef[currentVoxelIdx * 4 + 1], quatsArrayRef[currentVoxelIdx * 4 + 2],
                                                         quatsArrayRef[currentVoxelIdx * 4 + 3]);
                        const uint32 currentLaueIndex = crystalStructuresArrayRef[currentCellPhaseIdx];
                        if(referenceLaueIndex == currentLaueIndex && referenceLaueIndex < orientationOps.size())
                        {
                          const ebsdlib::AxisAngleDType axisAngle = orientationOps[referenceLaueIndex]->calculateMisorientation(referenceQuat, currentQuat);
                          misorientationAngle = axisAngle[3];
                        }
                      }
                      if(misorientationAngle > misorientationTolerance)
                      {
                        disorientationScore++;
                      }
                    }
                    if(m_InputValues->UseMask)
                    {
                      if(maskComparePtr->isTrue(referenceVoxelIdx) && !maskComparePtr->isTrue(currentVoxelIdx))
                      {
                        disorientationScore++;
                      }
                      if(!maskComparePtr->isTrue(referenceVoxelIdx) && maskComparePtr->isTrue(currentVoxelIdx))
                      {
                        disorientationScore++;
                      }
                    }
                  }
                }
              }
              disorientationScore /= sampleCount;
              evaluatedShiftXIdx = xShiftOffset + previousXShift + halfXDimension;
              evaluatedShiftYIdx = yShiftOffset + previousYShift + halfYDimension;
              evaluatedShiftIdx = (dimensions[0] * evaluatedShiftYIdx) + evaluatedShiftXIdx;
              evaluatedShifts[evaluatedShiftIdx] = true;
              if(disorientationScore < minimumDisorientation ||
                 (disorientationScore == minimumDisorientation && ((llabs(xShiftOffset + previousXShift) < llabs(currentXShift)) || (llabs(yShiftOffset + previousYShift) < llabs(currentYShift)))))
              {
                currentXShift = xShiftOffset + previousXShift;
                currentYShift = yShiftOffset + previousYShift;
                minimumDisorientation = disorientationScore;
              }
            }
          }
        }
      }
      xShifts[shiftIdx] = xShifts[shiftIdx - 1] + currentXShift;
      yShifts[shiftIdx] = yShifts[shiftIdx - 1] + currentYShift;
    }
  }

  return {};
}
