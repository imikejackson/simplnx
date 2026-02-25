#include "CreateScatterChart.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Utilities/FilterUtilities.hpp"

#include <mcplotty/core/Brush.hpp>
#include <mcplotty/core/Color.hpp>
#include <mcplotty/core/Pen.hpp>
#include <mcplotty/plot/Plot.hpp>
#include <mcplotty/plot/PlotScatterChart.hpp>
#include <mcplotty/render/Symbol.hpp>
#include <mcplotty/core/Font.hpp>
#include <mcplotty/core/Text.hpp>
#include <mcplotty/plot/Plot.hpp>
#include <mcplotty/plot/PlotCurve.hpp>

#include <fmt/format.h>

#include <cstring>
#include <memory>
#include <vector>

using namespace nx::core;

namespace
{
struct ConvertToDoubleVectorFunctor
{
  template <typename T>
  std::vector<double> operator()(const IDataArray& dataArray)
  {
    const auto& typedArray = dynamic_cast<const DataArray<T>&>(dataArray);
    usize numTuples = typedArray.getNumberOfTuples();
    std::vector<double> result(numTuples);
    for(usize i = 0; i < numTuples; i++)
    {
      result[i] = static_cast<double>(typedArray[i]);
    }
    return result;
  }
};
} // namespace

// -----------------------------------------------------------------------------
CreateScatterChart::CreateScatterChart(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, CreateScatterChartInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(mesgHandler)
{
}

// -----------------------------------------------------------------------------
CreateScatterChart::~CreateScatterChart() noexcept = default;

// -----------------------------------------------------------------------------
Result<> CreateScatterChart::operator()()
{
  const auto& xDataArray = m_DataStructure.getDataRefAs<IDataArray>(m_InputValues->XDataPath);
  const auto& yDataArray = m_DataStructure.getDataRefAs<IDataArray>(m_InputValues->YDataPath);

  m_MessageHandler(IFilter::Message::Type::Info, "Converting input data to double...");

  // Convert input arrays to std::vector<double> for McPlotty
  std::vector<double> xValues = ExecuteDataFunction(ConvertToDoubleVectorFunctor{}, xDataArray.getDataType(), xDataArray);
  std::vector<double> yValues = ExecuteDataFunction(ConvertToDoubleVectorFunctor{}, yDataArray.getDataType(), yDataArray);

  if(m_ShouldCancel)
  {
    return {};
  }

  uint32 chartHeight = m_InputValues->ChartDimensions[0];
  uint32 chartWidth = m_InputValues->ChartDimensions[1];

  m_MessageHandler(IFilter::Message::Type::Info, fmt::format("Creating scatter chart ({}x{})...", chartWidth, chartHeight));

  // Create the McPlotty plot
  auto plot = std::make_shared<mcplotty::Plot>(static_cast<int>(chartWidth), static_cast<int>(chartHeight));
  plot->setTitle( {m_InputValues->ChartTitle, mcplotty::Font::defaultFont(mcplotty::FontFamily::LatoBold, 32.0) });
  plot->enableAxis(mcplotty::Plot::xBottom, true);
  plot->enableAxis(mcplotty::Plot::yLeft, true);

  mcplotty::Font axisTitleFont = mcplotty::Font::defaultFont(mcplotty::FontFamily::LatoItalic, 24.0);

  plot->setAxisTitle(mcplotty::Plot::xBottom, mcplotty::Text(m_InputValues->XAxisTitle, axisTitleFont));
  plot->setAxisTitle(mcplotty::Plot::yLeft, mcplotty::Text(m_InputValues->YAxisTitle, axisTitleFont));
  plot->setAxisAutoScale(mcplotty::Plot::xBottom, true);
  plot->setAxisAutoScale(mcplotty::Plot::yLeft, true);

  // Create scatter chart series and add data
  auto scatter = std::make_shared<mcplotty::PlotScatterChart>("Data");
  scatter->setSamples(xValues, yValues);

  // Set a default symbol style
  mcplotty::Symbol symbol(mcplotty::SymbolStyle::Ellipse, mcplotty::Brush(mcplotty::Color(0.2f, 0.5f, 0.9f)), mcplotty::Pen(mcplotty::Color(0.1f, 0.2f, 0.5f), 1.0), mcplotty::Size(6.0, 6.0));
  scatter->setSymbol(symbol);

  plot->addItem(scatter);

  if(m_ShouldCancel)
  {
    return {};
  }

  m_MessageHandler(IFilter::Message::Type::Info, "Rendering chart to image...");

  // Render the chart to RGBA pixel data
  std::vector<uint8_t> imageData = plot->renderToImage(static_cast<int>(chartWidth), static_cast<int>(chartHeight));

  if(m_ShouldCancel)
  {
    return {};
  }

  // Copy the rendered image data into the output DataArray, converting RGBA to RGB
  auto& outputArray = m_DataStructure.getDataRefAs<DataArray<uint8>>(m_InputValues->ImageDataPath);

  usize numPixels = static_cast<usize>(chartWidth) * static_cast<usize>(chartHeight);
  usize expectedRgbaBytes = numPixels * 4;
  if(imageData.size() != expectedRgbaBytes)
  {
    return MakeErrorResult(-3000, fmt::format("Rendered image size ({}) does not match expected size ({}).", imageData.size(), expectedRgbaBytes));
  }

  // Convert RGBA to RGB by dropping the alpha channel
  for(usize i = 0; i < numPixels; i++)
  {
    outputArray[i * 3 + 0] = imageData[i * 4 + 0];
    outputArray[i * 3 + 1] = imageData[i * 4 + 1];
    outputArray[i * 3 + 2] = imageData[i * 4 + 2];
  }

  m_MessageHandler(IFilter::Message::Type::Info, "Scatter chart creation complete.");

  return {};
}
