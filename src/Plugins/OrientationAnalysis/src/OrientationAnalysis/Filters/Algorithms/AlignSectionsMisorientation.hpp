#pragma once

#include "OrientationAnalysis/OrientationAnalysis_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Utilities/AlignSections.hpp"

#include <vector>

namespace nx::core
{

/**
 * @struct AlignSectionsMisorientationInputValues
 * @brief Stores the input values for section alignment by misorientation.
 */
struct ORIENTATIONANALYSIS_EXPORT AlignSectionsMisorientationInputValues
{
  DataPath ImageGeometryPath;
  bool UseMask;
  DataPath MaskArrayPath;

  float32 MisorientationTolerance;
  DataPath QuatsArrayPath;
  DataPath CellPhasesArrayPath;
  DataPath CrystalStructuresArrayPath;

  bool StoreAlignmentShifts;
  DataPath AlignmentAMPath;
  DataPath SlicesArrayPath;
  DataPath RelativeShiftsArrayPath;
  DataPath CumulativeShiftsArrayPath;
};

/**
 * @class AlignSectionsMisorientation
 * @brief Aligns adjacent sections by minimizing their crystallographic disorientation.
 */
class ORIENTATIONANALYSIS_EXPORT AlignSectionsMisorientation : public AlignSections
{
public:
  /**
   * @brief Constructs the alignment algorithm.
   * @param dataStructure Contains the geometry and arrays that the algorithm uses and modifies.
   * @param messageHandler Receives progress messages.
   * @param shouldCancel Stops execution when set.
   * @param inputValues Specifies the input paths, tolerance, mask use, and optional shift outputs.
   */
  AlignSectionsMisorientation(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel, AlignSectionsMisorientationInputValues* inputValues);

  /**
   * @brief Destroys the alignment algorithm.
   */
  ~AlignSectionsMisorientation() noexcept override;

  AlignSectionsMisorientation(const AlignSectionsMisorientation&) = delete;
  AlignSectionsMisorientation(AlignSectionsMisorientation&&) noexcept = delete;
  AlignSectionsMisorientation& operator=(const AlignSectionsMisorientation&) = delete;
  AlignSectionsMisorientation& operator=(AlignSectionsMisorientation&&) noexcept = delete;

  /**
   * @brief Calculates and applies the section-alignment shifts.
   * @return An error if a participating Phase cannot index the selected Crystal Structures array.
   */
  Result<> operator()();

protected:
  /**
   * @brief Calculates the relative X and Y shifts between adjacent sections.
   * @param xShifts Receives the cumulative X shift for each section.
   * @param yShifts Receives the cumulative Y shift for each section.
   * @return An error if a participating Phase cannot index the selected Crystal Structures array.
   */
  Result<> findShifts(std::vector<int64_t>& xShifts, std::vector<int64_t>& yShifts) override;

private:
  DataStructure& m_DataStructure;
  const AlignSectionsMisorientationInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;
};
} // namespace nx::core
