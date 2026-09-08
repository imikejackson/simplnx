#pragma once

#include "OrientationAnalysis/OrientationAnalysis_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"

namespace nx::core
{

/**
 * @struct BadDataNeighborOrientationCheckInputValues
 * @brief Stores the input values for the bad-data neighbor orientation check.
 */
struct ORIENTATIONANALYSIS_EXPORT BadDataNeighborOrientationCheckInputValues
{
  float32 MisorientationTolerance;
  int32 NumberOfNeighbors;
  DataPath ImageGeomPath;
  DataPath QuatsArrayPath;
  DataPath MaskArrayPath;
  DataPath CellPhasesArrayPath;
  DataPath CrystalStructuresArrayPath;
};

/**
 * @class BadDataNeighborOrientationCheck
 * @brief Converts bad voxels to good when enough neighbors have similar orientations.
 */
class ORIENTATIONANALYSIS_EXPORT BadDataNeighborOrientationCheck
{
public:
  /**
   * @brief Constructs the neighbor-orientation algorithm.
   * @param dataStructure Contains the geometry and arrays that the algorithm uses and modifies.
   * @param messageHandler Receives progress messages.
   * @param shouldCancel Stops execution when set.
   * @param inputValues Specifies the input paths, tolerance, and neighbor threshold.
   */
  BadDataNeighborOrientationCheck(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                                  BadDataNeighborOrientationCheckInputValues* inputValues);

  /**
   * @brief Destroys the neighbor-orientation algorithm.
   */
  ~BadDataNeighborOrientationCheck() noexcept;

  BadDataNeighborOrientationCheck(const BadDataNeighborOrientationCheck&) = delete;
  BadDataNeighborOrientationCheck(BadDataNeighborOrientationCheck&&) noexcept = delete;
  BadDataNeighborOrientationCheck& operator=(const BadDataNeighborOrientationCheck&) = delete;
  BadDataNeighborOrientationCheck& operator=(BadDataNeighborOrientationCheck&&) noexcept = delete;

  /**
   * @brief Applies the iterative neighbor-orientation correction.
   * @return An error if a participating Phase or Laue index is out of bounds.
   */
  Result<> operator()();

private:
  DataStructure& m_DataStructure;
  const BadDataNeighborOrientationCheckInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;
};

} // namespace nx::core
