#include "ComputeAvgCAxes.hpp"

#include "OrientationAnalysis/utilities/OrientationUtilities.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/Utilities/ImageRotationUtilities.hpp"
#include "simplnx/Utilities/Math/GeometryMath.hpp"

#include <EbsdLib/Core/Orientation.hpp>
#include <EbsdLib/Orientation/OrientationFwd.hpp>
#include <EbsdLib/Orientation/OrientationMatrix.hpp>
#include <EbsdLib/Orientation/Quaternion.hpp>

using namespace nx::core;
using namespace nx::core::OrientationUtilities;

// -----------------------------------------------------------------------------
ComputeAvgCAxes::ComputeAvgCAxes(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel, ComputeAvgCAxesInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(messageHandler)
{
}

// -----------------------------------------------------------------------------
ComputeAvgCAxes::~ComputeAvgCAxes() noexcept = default;

// -----------------------------------------------------------------------------
Result<> ComputeAvgCAxes::operator()()
{

  // Figure out if all phases are either Hexagonal-Low 6/m or Hexagonal-High 6/mmm Laue Phases
  const auto& crystalStructuresArrayRef = m_DataStructure.getDataRefAs<UInt32Array>(m_InputValues->CrystalStructuresArrayPath);
  const usize numCrystalStructures = crystalStructuresArrayRef.getNumberOfTuples();
  bool allPhasesHexagonal = true;
  bool noPhasesHexagonal = true;
  for(usize phaseIdx = 1; phaseIdx < numCrystalStructures; ++phaseIdx)
  {
    const uint32 currentLaueIndex = crystalStructuresArrayRef[phaseIdx];
    const bool isHex = currentLaueIndex == ebsdlib::CrystalStructure::Hexagonal_High || currentLaueIndex == ebsdlib::CrystalStructure::Hexagonal_Low;
    allPhasesHexagonal = allPhasesHexagonal && isHex;
    noPhasesHexagonal = noPhasesHexagonal && !isHex;
  }

  // Return an error when no Phase has hexagonal symmetry.
  if(noPhasesHexagonal)
  {
    return MakeErrorResult(-76402, "No phases that have a crystal symmetry of Hexagonal (6/mmm or 6/m) were found.");
  }

  Result<> result;

  // Throw a warning for any NON-Hex Laue Phases
  if(!allPhasesHexagonal)
  {
    result.warnings().push_back({-76403, "Non Hexagonal phases were found. All calculations for non Hexagonal phases will be skipped and a NaN value inserted."});
  }

  const auto& featureIdsArrayRef = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->FeatureIdsArrayPath);
  const auto& quatsArrayRef = m_DataStructure.getDataRefAs<Float32Array>(m_InputValues->QuatsArrayPath);
  const auto& cellPhasesArrayRef = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->CellPhasesArrayPath);
  auto& avgCAxesArrayRef = m_DataStructure.getDataRefAs<Float32Array>(m_InputValues->AvgCAxesArrayPath);
  avgCAxesArrayRef.fill(0.0f); // Initialize all output values to ZERO defensively.

  const usize totalVoxels = featureIdsArrayRef.getNumberOfTuples();
  const usize totalFeatures = avgCAxesArrayRef.getNumberOfTuples();

  const Eigen::Vector3d cAxis{0.0, 0.0, 1.0};

  std::vector<int32> cellCounts(totalFeatures, 0);

  m_MessageHandler({IFilter::Message::Type::Info, "Computing cell contributions"});

  // Loop over each cell
  for(usize voxelIdx = 0; voxelIdx < totalVoxels; voxelIdx++)
  {
    if(m_ShouldCancel)
    {
      return result;
    }

    const int32 currentFeatureIdx = featureIdsArrayRef[voxelIdx];
    // If the featureId for a given cell is valid ( > 0) then analyze that value
    if(currentFeatureIdx > 0)
    {
      const int32 currentCellPhaseIdx = cellPhasesArrayRef[voxelIdx];
      if(currentCellPhaseIdx <= 0)
      {
        continue;
      }
      if(static_cast<usize>(currentCellPhaseIdx) >= numCrystalStructures)
      {
        return MakeErrorResult(-76404, fmt::format("Cell Phases array '{}' has value {} at voxel index {}, but Crystal Structures array '{}' has {} tuples. Valid Phase indices are in [0, {}).",
                                                   m_InputValues->CellPhasesArrayPath.toString(), currentCellPhaseIdx, voxelIdx, m_InputValues->CrystalStructuresArrayPath.toString(),
                                                   numCrystalStructures, numCrystalStructures));
      }
      const uint32 currentLaueIndex = crystalStructuresArrayRef[currentCellPhaseIdx];
      const usize cAxisOffset = 3 * currentFeatureIdx;

      // If the Laue class is not Hexagonal, then continue to the next cell
      if(currentLaueIndex != ebsdlib::CrystalStructure::Hexagonal_High && currentLaueIndex != ebsdlib::CrystalStructure::Hexagonal_Low)
      {
        continue;
      }

      cellCounts[currentFeatureIdx]++;
      const usize quatOffset = voxelIdx * 4;

      // Create the 3x3 Orientation Matrix from the Quaternion. This represents a passive rotation matrix
      const ebsdlib::OrientationMatrixDType orientationMatrix =
          ebsdlib::QuaternionDType(quatsArrayRef[quatOffset], quatsArrayRef[quatOffset + 1], quatsArrayRef[quatOffset + 2], quatsArrayRef[quatOffset + 3]).toOrientationMatrix();

      // Convert the passive rotation matrix to an active rotation matrix by taking the transpose
      // Multiply the active transformation matrix by the C-Axis (as Miller Index). This actively rotates
      // the crystallographic C-Axis (which is along the <0,0,1> direction) into the physical sample
      // reference frame
      Eigen::Vector3d cellCAxis = orientationMatrix.transpose() * cAxis;

      // normalize so that the magnitude is 1
      cellCAxis.normalize();

      // Compute the running average c-axis and normalize the result
      Eigen::Vector3d runningCAxisAvg{avgCAxesArrayRef[cAxisOffset] / static_cast<float32>(cellCounts[currentFeatureIdx]),
                                      avgCAxesArrayRef[cAxisOffset + 1] / static_cast<float32>(cellCounts[currentFeatureIdx]),
                                      avgCAxesArrayRef[cAxisOffset + 2] / static_cast<float32>(cellCounts[currentFeatureIdx])};
      runningCAxisAvg.normalize();

      // Ensure that angle between the current point's sample reference frame C-Axis
      // and the running average sample C-Axis is positive
      float64 cosAngle = ImageRotationUtilities::CosBetweenVectors(cellCAxis, runningCAxisAvg);
      if(cosAngle < 0.0)
      {
        cellCAxis *= -1.0f;
      }

      // Accumulate per-component into the float32 output (Eigen math is double; narrow on store).
      avgCAxesArrayRef[cAxisOffset] = static_cast<float32>(avgCAxesArrayRef[cAxisOffset] + cellCAxis[0]);
      avgCAxesArrayRef[cAxisOffset + 1] = static_cast<float32>(avgCAxesArrayRef[cAxisOffset + 1] + cellCAxis[1]);
      avgCAxesArrayRef[cAxisOffset + 2] = static_cast<float32>(avgCAxesArrayRef[cAxisOffset + 2] + cellCAxis[2]);
    }
  }

  // Compute the final average C-axis for each Feature.
  m_MessageHandler({IFilter::Message::Type::Info, "Computing final feature average C-Axis values"});

  for(usize featureIdx = 0; featureIdx < totalFeatures; featureIdx++)
  {
    if(m_ShouldCancel)
    {
      return result;
    }

    const usize cAxisOffset = 3 * featureIdx;
    if(cellCounts[featureIdx] == 0)
    {
      // Feature is either non-hexagonal or has no assigned voxels; either way, no meaningful average exists.
      avgCAxesArrayRef[cAxisOffset] = NAN;
      avgCAxesArrayRef[cAxisOffset + 1] = NAN;
      avgCAxesArrayRef[cAxisOffset + 2] = NAN;
    }
    else
    {
      // Divide the accumulated sum by the cell count, then normalize so the
      // output is a unit-magnitude C-axis direction. The antipodal-flip rule
      // guarantees |sum| >= sqrt(cellCounts), so the divided vector's magnitude
      // is >= 1/sqrt(cellCounts) > 0 -- no near-zero guard needed.
      Eigen::Vector3d finalAverageCAxis{avgCAxesArrayRef[cAxisOffset] / static_cast<float64>(cellCounts[featureIdx]), avgCAxesArrayRef[cAxisOffset + 1] / static_cast<float64>(cellCounts[featureIdx]),
                                        avgCAxesArrayRef[cAxisOffset + 2] / static_cast<float64>(cellCounts[featureIdx])};
      finalAverageCAxis.normalize();
      avgCAxesArrayRef[cAxisOffset] = static_cast<float32>(finalAverageCAxis[0]);
      avgCAxesArrayRef[cAxisOffset + 1] = static_cast<float32>(finalAverageCAxis[1]);
      avgCAxesArrayRef[cAxisOffset + 2] = static_cast<float32>(finalAverageCAxis[2]);
    }
  }
  return result;
}
