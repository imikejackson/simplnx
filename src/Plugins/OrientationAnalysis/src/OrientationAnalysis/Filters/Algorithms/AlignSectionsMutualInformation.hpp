#pragma once

#include "OrientationAnalysis/OrientationAnalysis_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/Utilities/AlignSections.hpp"
#include "simplnx/Utilities/MaskCompareUtilities.hpp"

namespace nx::core
{
/**
 * @struct AlignSectionsMutualInformationInputValues
 * @brief Stores the input values for section alignment by mutual information.
 */
struct ORIENTATIONANALYSIS_EXPORT AlignSectionsMutualInformationInputValues
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
 * @class AlignSectionsMutualInformation
 * @brief Aligns adjacent sections by maximizing mutual information between segmented regions.
 */
class ORIENTATIONANALYSIS_EXPORT AlignSectionsMutualInformation : public AlignSections
{
public:
  /**
   * @brief Constructs the section-alignment algorithm.
   * @param dataStructure Contains the geometry and arrays that the algorithm uses and modifies.
   * @param messageHandler Receives progress messages.
   * @param shouldCancel Stops execution when set.
   * @param inputValues Specifies the input paths, tolerance, mask use, and optional shift outputs.
   */
  AlignSectionsMutualInformation(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                                 AlignSectionsMutualInformationInputValues* inputValues);

  /**
   * @brief Destroys the section-alignment algorithm.
   */
  ~AlignSectionsMutualInformation() noexcept override;

  AlignSectionsMutualInformation(const AlignSectionsMutualInformation&) = delete;
  AlignSectionsMutualInformation(AlignSectionsMutualInformation&&) noexcept = delete;
  AlignSectionsMutualInformation& operator=(const AlignSectionsMutualInformation&) = delete;
  AlignSectionsMutualInformation& operator=(AlignSectionsMutualInformation&&) noexcept = delete;

  /**
   * @brief Calculates and applies the section-alignment shifts.
   * @return An error if a participating Phase or Laue index is out of bounds.
   */
  Result<> operator()();

protected:
  /**
   * @brief Calculates the relative X and Y shifts between adjacent sections.
   * @param xShifts Receives the cumulative X shift for each section.
   * @param yShifts Receives the cumulative Y shift for each section.
   * @return An error if section segmentation fails.
   */
  Result<> findShifts(std::vector<int64>& xShifts, std::vector<int64>& yShifts) override;

  /**
   * @brief Segments each section into regions for the mutual-information calculation.
   * @param sliceFeatureIds Receives the temporary region identifier for each voxel.
   * @param sliceFeatureCounts Receives the number of regions in each section.
   * @return An error if a participating Phase or Laue index is out of bounds.
   */
  Result<> formFeaturesSections(std::vector<int32>& sliceFeatureIds, std::vector<int32>& sliceFeatureCounts);

private:
  DataStructure& m_DataStructure;
  const AlignSectionsMutualInformationInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;

  std::unique_ptr<MaskCompareUtilities::MaskCompare> m_MaskCompare = nullptr;
};

} // namespace nx::core
