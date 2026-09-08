#pragma once

#include "OrientationAnalysis/OrientationAnalysis_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"

namespace nx::core
{

/**
 * @struct ComputeAvgCAxesInputValues
 * @brief Stores the input values for the average C-axis calculation.
 */
struct ORIENTATIONANALYSIS_EXPORT ComputeAvgCAxesInputValues
{
  DataPath QuatsArrayPath;
  DataPath FeatureIdsArrayPath;
  DataPath CellPhasesArrayPath;
  DataPath CellFeatureDataPath;
  DataPath AvgCAxesArrayPath;
  DataPath CrystalStructuresArrayPath;
};

/**
 * @class ComputeAvgCAxes
 * @brief Calculates the average C-axis direction of each Feature.
 */

class ORIENTATIONANALYSIS_EXPORT ComputeAvgCAxes
{
public:
  /**
   * @brief Constructs the average C-axis algorithm.
   * @param dataStructure Contains the arrays that the algorithm uses and modifies.
   * @param messageHandler Receives progress messages.
   * @param shouldCancel Stops execution when set.
   * @param inputValues Specifies the input and output paths.
   */
  ComputeAvgCAxes(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel, ComputeAvgCAxesInputValues* inputValues);

  /**
   * @brief Destroys the average C-axis algorithm.
   */
  ~ComputeAvgCAxes() noexcept;

  ComputeAvgCAxes(const ComputeAvgCAxes&) = delete;
  ComputeAvgCAxes(ComputeAvgCAxes&&) noexcept = delete;
  ComputeAvgCAxes& operator=(const ComputeAvgCAxes&) = delete;
  ComputeAvgCAxes& operator=(ComputeAvgCAxes&&) noexcept = delete;

  /**
   * @brief Calculates the average C-axis direction for each eligible Feature.
   * @return An error if a participating Phase cannot index the Crystal Structures array.
   */
  Result<> operator()();

private:
  DataStructure& m_DataStructure;
  const ComputeAvgCAxesInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;
};

} // namespace nx::core
