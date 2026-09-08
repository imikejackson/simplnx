#include "BadDataNeighborOrientationCheck.hpp"

#include "simplnx/Common/Numbers.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Utilities/MaskCompareUtilities.hpp"
#include "simplnx/Utilities/MessageHelper.hpp"
#include "simplnx/Utilities/NeighborUtilities.hpp"

#include <EbsdLib/LaueOps/LaueOps.h>

using namespace nx::core;

// -----------------------------------------------------------------------------
BadDataNeighborOrientationCheck::BadDataNeighborOrientationCheck(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                                                                 BadDataNeighborOrientationCheckInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(messageHandler)
{
}

// -----------------------------------------------------------------------------
BadDataNeighborOrientationCheck::~BadDataNeighborOrientationCheck() noexcept = default;

// -----------------------------------------------------------------------------
Result<> BadDataNeighborOrientationCheck::operator()()
{
  // Convert the tolerance with double-precision pi. The float value of pi is slightly larger than the exact value.
  // A float conversion can incorrectly include a misorientation that is exactly on the strict tolerance boundary.
  const double misorientationTolerance = static_cast<double>(m_InputValues->MisorientationTolerance) * numbers::pi_v<double> / 180.0;

  const auto& imageGeom = m_DataStructure.getDataRefAs<ImageGeom>(m_InputValues->ImageGeomPath);
  const SizeVec3 gridDimensions = imageGeom.getDimensions();
  const auto& cellPhasesArrayRef = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->CellPhasesArrayPath);
  const auto& quatsArrayRef = m_DataStructure.getDataRefAs<Float32Array>(m_InputValues->QuatsArrayPath);
  const auto& crystalStructuresArrayRef = m_DataStructure.getDataRefAs<UInt32Array>(m_InputValues->CrystalStructuresArrayPath);
  const usize totalVoxels = quatsArrayRef.getNumberOfTuples();
  const usize numCrystalStructures = crystalStructuresArrayRef.getNumberOfTuples();

  std::unique_ptr<MaskCompareUtilities::MaskCompare> maskComparePtr;
  try
  {
    maskComparePtr = MaskCompareUtilities::InstantiateMaskCompare(m_DataStructure, m_InputValues->MaskArrayPath);
  } catch(const std::out_of_range& exception)
  {
    // Defensive: the path was verified during preflight, but this algorithm may be called outside the standard
    // IFilter Preflight/Execute path.
    return MakeErrorResult(-54900,
                           fmt::format("Mask Array at '{}' could not be loaded; expected Bool or UInt8 backing. Underlying error: {}", m_InputValues->MaskArrayPath.toString(), exception.what()));
  }

  const std::array<int64, 3> dimensions = {
      static_cast<int64>(gridDimensions[0]),
      static_cast<int64>(gridDimensions[1]),
      static_cast<int64>(gridDimensions[2]),
  };

  // VoxelNeighbors<Image3D>::k_FaceNeighborCount = 6 is the maximum possible face-neighbor count.
  // computeValidFaceNeighbors() skips +/-Z neighbors when dimensions[2] is 1, so this
  // 3D-typed array correctly handles 2D images without any change here.
  constexpr FaceNeighborType k_NumFaceNeighbors = VoxelNeighbors<Image3D>::k_FaceNeighborCount;
  const std::array<int64, k_NumFaceNeighbors> neighborVoxelOffsets = initializeFaceNeighborOffsets(dimensions);
  constexpr std::array<FaceNeighborType, k_NumFaceNeighbors> faceNeighborIndices = initializeFaceNeighborInternalIdx();

  const std::vector<ebsdlib::LaueOps::Pointer> orientationOps = ebsdlib::LaueOps::GetAllOrientationOps();

  // Validate each Crystal Structures value before the voxel loops. Allow UnknownCrystalStructure as the Phase 0 sentinel.
  // A voxel that resolves to the sentinel does not use a Laue operation.
  const usize numOrientationOps = orientationOps.size();
  for(usize ensembleIdx = 0; ensembleIdx < crystalStructuresArrayRef.getSize(); ++ensembleIdx)
  {
    if(crystalStructuresArrayRef[ensembleIdx] >= numOrientationOps && crystalStructuresArrayRef[ensembleIdx] != ebsdlib::CrystalStructure::UnknownCrystalStructure)
    {
      return MakeErrorResult(-54901, fmt::format("Crystal structure at ensemble index {} has value {}, which is not a valid Laue-group index. Valid range is [0, {}).", ensembleIdx,
                                                 crystalStructuresArrayRef[ensembleIdx], numOrientationOps));
    }
  }

  // Per-voxel running count of within-tolerance face-neighbors. Allocated proportional to the
  // input geometry size: 4 bytes per voxel (~4 GB for a 1B-voxel dataset). Cannot be in-place
  // on the mask array because the algorithm needs to distinguish "newly flipped" from "still bad".
  std::vector<int32> neighborCounts(totalVoxels, 0);

  MessageHelper messageHelper(m_MessageHandler);
  ThrottledMessenger throttledMessenger = messageHelper.createThrottledMessenger();
  // Loop over every point finding the number of neighbors that fall within the
  // user defined angle tolerance.
  for(usize voxelIdx = 0; voxelIdx < totalVoxels; voxelIdx++)
  {
    if(m_ShouldCancel)
    {
      return {};
    }
    throttledMessenger.sendThrottledMessage([&] { return fmt::format("Processing Data {:.2f}% completed", CalculatePercentComplete(voxelIdx, totalVoxels)); });
    // If the mask was set to false, then we check this voxel
    // "Bad" voxels are those whose mask value is false; only these get processed.
    const bool voxelIsBad = !maskComparePtr->isTrue(voxelIdx);
    if(voxelIsBad)
    {
      // We precalculate the positive voxel quaternion and laue class here to prevent reading and recalculating it for each face below
      ebsdlib::QuatD currentQuat(quatsArrayRef[voxelIdx * 4], quatsArrayRef[voxelIdx * 4 + 1], quatsArrayRef[voxelIdx * 4 + 2], quatsArrayRef[voxelIdx * 4 + 3]);
      currentQuat.positiveOrientation();
      const int32 currentCellPhaseIdx = cellPhasesArrayRef[voxelIdx];
      if(currentCellPhaseIdx <= 0)
      {
        continue;
      }
      if(static_cast<usize>(currentCellPhaseIdx) >= numCrystalStructures)
      {
        return MakeErrorResult(-54902, fmt::format("Cell Phases array '{}' has value {} at voxel index {}, but Crystal Structures array '{}' has {} tuples. Valid Phase indices are in [0, {}).",
                                                   m_InputValues->CellPhasesArrayPath.toString(), currentCellPhaseIdx, voxelIdx, m_InputValues->CrystalStructuresArrayPath.toString(),
                                                   numCrystalStructures, numCrystalStructures));
      }
      const uint32 currentLaueIndex = crystalStructuresArrayRef[currentCellPhaseIdx];
      // Defensive: skip voxels whose phase resolves to an out-of-range Laue index (e.g., the
      // UnknownCrystalStructure sentinel allowed by the validation above). Without this, the
      // orientationOps[currentLaueIndex] dereference below would be out-of-bounds.
      if(currentLaueIndex >= numOrientationOps)
      {
        continue;
      }

      const int64 voxelIdxI64 = static_cast<int64>(voxelIdx);
      int64 xIdx = voxelIdxI64 % dimensions[0];
      int64 yIdx = (voxelIdxI64 / dimensions[0]) % dimensions[1];
      int64 zIdx = voxelIdxI64 / (dimensions[0] * dimensions[1]);

      // Loop over the 6 face neighbors of the voxel
      const std::array<bool, k_NumFaceNeighbors> isValidFaceNeighbor = computeValidFaceNeighbors(xIdx, yIdx, zIdx, dimensions);
      for(const auto& faceIdx : faceNeighborIndices)
      {
        if(!isValidFaceNeighbor[faceIdx])
        {
          continue;
        }
        const int64 neighborVoxelIdx = voxelIdxI64 + neighborVoxelOffsets[faceIdx];

        // Compare orientations only when the neighbor mask identifies a good voxel.
        if(maskComparePtr->isTrue(neighborVoxelIdx))
        {
          // Both Cell Phases MUST be the same and be a valid Phase
          if(cellPhasesArrayRef[voxelIdx] == cellPhasesArrayRef[neighborVoxelIdx] && cellPhasesArrayRef[voxelIdx] > 0)
          {
            ebsdlib::QuatD neighborQuat(quatsArrayRef[neighborVoxelIdx * 4], quatsArrayRef[neighborVoxelIdx * 4 + 1], quatsArrayRef[neighborVoxelIdx * 4 + 2], quatsArrayRef[neighborVoxelIdx * 4 + 3]);
            neighborQuat.positiveOrientation();
            // Compute the Axis_Angle misorientation between those 2 quaternions
            ebsdlib::AxisAngleDType axisAngle = orientationOps[currentLaueIndex]->calculateMisorientation(currentQuat, neighborQuat);
            // if the angle is less than our tolerance, then we increment the neighbor count
            // for this voxel
            if(axisAngle[3] < misorientationTolerance)
            {
              neighborCounts[voxelIdx]++;
            }
          }
        }
      }
    }
  }

  // Start at the maximum face-neighbor count. A 2D image cannot reach the first two levels.
  // The constant keeps this limit consistent with VoxelNeighbors.
  constexpr int32 startLevel = static_cast<int32>(k_NumFaceNeighbors);
  int32 currentLevel = startLevel;
  int32 counter = 0;

  // Repeat each level until no additional bad voxels become good.
  while(currentLevel >= m_InputValues->NumberOfNeighbors)
  {
    if(m_ShouldCancel)
    {
      return {};
    }
    counter = 1;
    int32 loopNumber = 0;
    while(counter > 0)
    {
      if(m_ShouldCancel)
      {
        return {};
      }
      counter = 0; // Set this while control variable to zero
      for(usize voxelIdx = 0; voxelIdx < totalVoxels; voxelIdx++)
      {
        if(m_ShouldCancel)
        {
          return {};
        }
        throttledMessenger.sendThrottledMessage([&] {
          return fmt::format("Level '{}' of '{}' || Processing Data ('{}') {:.2f}% completed", (startLevel - currentLevel) + 1, startLevel - m_InputValues->NumberOfNeighbors, loopNumber,
                             CalculatePercentComplete(voxelIdx, totalVoxels));
        });

        // If the current voxel's neighbor count is >= the current level and the mask is FALSE,
        // we flip the voxel to TRUE and recompute its (still-bad) neighbors' counts below.
        const bool voxelIsBad = !maskComparePtr->isTrue(voxelIdx);
        if(neighborCounts[voxelIdx] >= currentLevel && voxelIsBad)
        {
          maskComparePtr->setValue(voxelIdx, true);
          counter++; // Increment the `counter` to force the loop to iterate again

          // We precalculate the positive voxel quaternion and laue class here to prevent reading and recalculating it for each face below
          ebsdlib::QuatD currentQuat(quatsArrayRef[voxelIdx * 4], quatsArrayRef[voxelIdx * 4 + 1], quatsArrayRef[voxelIdx * 4 + 2], quatsArrayRef[voxelIdx * 4 + 3]);
          currentQuat.positiveOrientation();
          const int32 currentCellPhaseIdx = cellPhasesArrayRef[voxelIdx];
          if(currentCellPhaseIdx <= 0)
          {
            continue;
          }
          if(static_cast<usize>(currentCellPhaseIdx) >= numCrystalStructures)
          {
            return MakeErrorResult(-54902, fmt::format("Cell Phases array '{}' has value {} at voxel index {}, but Crystal Structures array '{}' has {} tuples. Valid Phase indices are in [0, {}).",
                                                       m_InputValues->CellPhasesArrayPath.toString(), currentCellPhaseIdx, voxelIdx, m_InputValues->CrystalStructuresArrayPath.toString(),
                                                       numCrystalStructures, numCrystalStructures));
          }
          const uint32 currentLaueIndex = crystalStructuresArrayRef[currentCellPhaseIdx];
          // Defensive: skip voxels with out-of-range Laue index. See matching guard in pass 1.
          if(currentLaueIndex >= numOrientationOps)
          {
            continue;
          }

          // Update each bad neighbor after the current voxel becomes good. This update permits valid cascade changes in later iterations.
          const int64 voxelIdxI64 = static_cast<int64>(voxelIdx);
          int64 xIdx = voxelIdxI64 % dimensions[0];
          int64 yIdx = (voxelIdxI64 / dimensions[0]) % dimensions[1];
          int64 zIdx = voxelIdxI64 / (dimensions[0] * dimensions[1]);

          // Loop over the 6 face neighbors of the voxel
          const std::array<bool, k_NumFaceNeighbors> isValidFaceNeighbor = computeValidFaceNeighbors(xIdx, yIdx, zIdx, dimensions);
          for(const auto& faceIdx : faceNeighborIndices)
          {
            if(!isValidFaceNeighbor[faceIdx])
            {
              continue;
            }

            const int64 neighborVoxelIdx = voxelIdxI64 + neighborVoxelOffsets[faceIdx];

            // If the neighbor voxel's mask is false, then compute misorientation angle
            const bool neighborIsBad = !maskComparePtr->isTrue(neighborVoxelIdx);
            if(neighborIsBad)
            {
              // Make sure both cells phase values are identical and valid
              if(cellPhasesArrayRef[voxelIdx] == cellPhasesArrayRef[neighborVoxelIdx] && cellPhasesArrayRef[voxelIdx] > 0)
              {
                ebsdlib::QuatD neighborQuat(quatsArrayRef[neighborVoxelIdx * 4], quatsArrayRef[neighborVoxelIdx * 4 + 1], quatsArrayRef[neighborVoxelIdx * 4 + 2],
                                            quatsArrayRef[neighborVoxelIdx * 4 + 3]);
                neighborQuat.positiveOrientation();
                // Quaternion Math is not commutative so do not reorder
                ebsdlib::AxisAngleDType axisAngle = orientationOps[currentLaueIndex]->calculateMisorientation(currentQuat, neighborQuat);
                if(axisAngle[3] < misorientationTolerance)
                {
                  neighborCounts[neighborVoxelIdx]++;
                }
              }
            }
          }
        }
      }
      ++loopNumber;
    }
    currentLevel = currentLevel - 1;
  }

  return {};
}
