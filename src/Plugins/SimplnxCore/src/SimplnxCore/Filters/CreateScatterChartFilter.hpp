#pragma once

#include "SimplnxCore/SimplnxCore_export.hpp"

#include "simplnx/Filter/FilterTraits.hpp"
#include "simplnx/Filter/IFilter.hpp"

namespace nx::core
{
/**
 * @class CreateScatterChartFilter
 * @brief This filter creates a scatter chart image from X and Y data arrays
 * using the McPlotty library and stores the result in an output ImageGeometry.
 */
class SIMPLNXCORE_EXPORT CreateScatterChartFilter : public IFilter
{
public:
  CreateScatterChartFilter() = default;
  ~CreateScatterChartFilter() noexcept override = default;

  CreateScatterChartFilter(const CreateScatterChartFilter&) = delete;
  CreateScatterChartFilter(CreateScatterChartFilter&&) noexcept = delete;

  CreateScatterChartFilter& operator=(const CreateScatterChartFilter&) = delete;
  CreateScatterChartFilter& operator=(CreateScatterChartFilter&&) noexcept = delete;

  // Parameter Keys
  static constexpr StringLiteral k_XDataPath_Key = "x_data_path";
  static constexpr StringLiteral k_YDataPath_Key = "y_data_path";
  static constexpr StringLiteral k_ChartDimensions_Key = "chart_dimensions";
  static constexpr StringLiteral k_ChartTitle_Key = "chart_title";
  static constexpr StringLiteral k_XAxisTitle_Key = "x_axis_title";
  static constexpr StringLiteral k_YAxisTitle_Key = "y_axis_title";
  static constexpr StringLiteral k_OutputImageGeometryPath_Key = "output_image_geometry_path";
  static constexpr StringLiteral k_CellDataName_Key = "cell_data_name";
  static constexpr StringLiteral k_ImageDataName_Key = "image_data_name";

  /**
   * @brief Returns the name of the filter.
   * @return
   */
  std::string name() const override;

  /**
   * @brief Returns the C++ classname of this filter.
   * @return
   */
  std::string className() const override;

  /**
   * @brief Returns the uuid of the filter.
   * @return
   */
  Uuid uuid() const override;

  /**
   * @brief Returns the human readable name of the filter.
   * @return
   */
  std::string humanName() const override;

  /**
   * @brief Returns the default tags for this filter.
   * @return
   */
  std::vector<std::string> defaultTags() const override;

  /**
   * @brief Returns the parameters of the filter (i.e. its inputs)
   * @return
   */
  Parameters parameters() const override;

  /**
   * @brief Returns parameters version integer.
   * Initial version should always be 1.
   * Should be incremented everytime the parameters change.
   * @return VersionType
   */
  VersionType parametersVersion() const override;

  /**
   * @brief Returns a copy of the filter.
   * @return
   */
  UniquePointer clone() const override;

protected:
  /**
   * @brief Takes in a DataStructure and checks that the filter can be run on it with the given arguments.
   * Returns any warnings/errors. Also returns the changes that would be applied to the DataStructure.
   * @param dataStructure The input DataStructure instance
   * @param filterArgs These are the input values for each parameter that is required for the filter
   * @param messageHandler The MessageHandler object
   * @param shouldCancel Atomic boolean value that can be checked to cancel the filter
   * @param executionContext The ExecutionContext that can be used to determine the correct absolute path from a relative path
   * @return Returns a Result object with error or warning values if any of those occurred during execution of this function
   */
  PreflightResult preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                                const ExecutionContext& executionContext) const override;

  /**
   * @brief Applies the filter's algorithm to the DataStructure with the given arguments. Returns any warnings/errors.
   * On failure, there is no guarantee that the DataStructure is in a correct state.
   * @param dataStructure The input DataStructure instance
   * @param filterArgs These are the input values for each parameter that is required for the filter
   * @param pipelineNode The node in the pipeline that is being executed
   * @param messageHandler The MessageHandler object
   * @param shouldCancel Atomic boolean value that can be checked to cancel the filter
   * @param executionContext The ExecutionContext that can be used to determine the correct absolute path from a relative path
   * @return Returns a Result object with error or warning values if any of those occurred during execution of this function
   */
  Result<> executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                       const ExecutionContext& executionContext) const override;
};
} // namespace nx::core

SIMPLNX_DEF_FILTER_TRAITS(nx::core, CreateScatterChartFilter, "b81f8631-10bc-4637-8789-af36ac2e656d");
