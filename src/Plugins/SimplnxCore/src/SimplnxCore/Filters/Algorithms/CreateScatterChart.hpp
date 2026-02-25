#pragma once

#include "SimplnxCore/SimplnxCore_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"

#include <string>
#include <vector>

namespace nx::core
{

struct SIMPLNXCORE_EXPORT CreateScatterChartInputValues
{
  DataPath XDataPath;
  DataPath YDataPath;
  std::vector<uint32> ChartDimensions;
  std::string ChartTitle;
  std::string XAxisTitle;
  std::string YAxisTitle;
  DataPath OutputImageGeometryPath;
  DataPath ImageDataPath;
};

/**
 * @class CreateScatterChart
 * @brief This algorithm uses McPlotty to create a scatter chart from X and Y
 * data arrays and writes the resulting RGBA image into an output ImageGeometry.
 */
class SIMPLNXCORE_EXPORT CreateScatterChart
{
public:
  CreateScatterChart(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, CreateScatterChartInputValues* inputValues);
  ~CreateScatterChart() noexcept;

  CreateScatterChart(const CreateScatterChart&) = delete;
  CreateScatterChart(CreateScatterChart&&) noexcept = delete;
  CreateScatterChart& operator=(const CreateScatterChart&) = delete;
  CreateScatterChart& operator=(CreateScatterChart&&) noexcept = delete;

  Result<> operator()();

private:
  DataStructure& m_DataStructure;
  const CreateScatterChartInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;
};

} // namespace nx::core
