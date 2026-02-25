#include "CreateScatterChartFilter.hpp"

#include "SimplnxCore/Filters/Algorithms/CreateScatterChart.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"
#include "simplnx/Parameters/VectorParameter.hpp"

using namespace std::string_literals;
using namespace nx::core;

namespace nx::core
{
//------------------------------------------------------------------------------
std::string CreateScatterChartFilter::name() const
{
  return FilterTraits<CreateScatterChartFilter>::name.str();
}

//------------------------------------------------------------------------------
std::string CreateScatterChartFilter::className() const
{
  return FilterTraits<CreateScatterChartFilter>::className;
}

//------------------------------------------------------------------------------
Uuid CreateScatterChartFilter::uuid() const
{
  return FilterTraits<CreateScatterChartFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string CreateScatterChartFilter::humanName() const
{
  return "Create Scatter Chart";
}

//------------------------------------------------------------------------------
std::vector<std::string> CreateScatterChartFilter::defaultTags() const
{
  return {className(), "Core", "Visualization", "Chart", "Scatter Plot"};
}

//------------------------------------------------------------------------------
Parameters CreateScatterChartFilter::parameters() const
{
  Parameters params;

  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});

  params.insert(std::make_unique<ArraySelectionParameter>(k_XDataPath_Key, "X Data", "The data array to use for the X axis of the scatter chart", DataPath{},
                                                          ArraySelectionParameter::AllowedTypes{DataType::float32, DataType::float64, DataType::int8, DataType::int16, DataType::int32, DataType::int64,
                                                                                               DataType::uint8, DataType::uint16, DataType::uint32, DataType::uint64},
                                                          ArraySelectionParameter::AllowedComponentShapes{{1}}));

  params.insert(std::make_unique<ArraySelectionParameter>(k_YDataPath_Key, "Y Data", "The data array to use for the Y axis of the scatter chart", DataPath{},
                                                          ArraySelectionParameter::AllowedTypes{DataType::float32, DataType::float64, DataType::int8, DataType::int16, DataType::int32, DataType::int64,
                                                                                               DataType::uint8, DataType::uint16, DataType::uint32, DataType::uint64},
                                                          ArraySelectionParameter::AllowedComponentShapes{{1}}));

  params.insert(std::make_unique<VectorUInt32Parameter>(k_ChartDimensions_Key, "Chart Pixel Dimensions (H x W)", "The height and width of the output chart image in pixels",
                                                        std::vector<uint32>{800, 600}, std::vector<std::string>{"H"s, "W"s}));

  params.insert(std::make_unique<StringParameter>(k_ChartTitle_Key, "Chart Title", "The title displayed at the top of the chart", "Scatter Chart"));
  params.insert(std::make_unique<StringParameter>(k_XAxisTitle_Key, "X Axis Title", "The label for the X axis", "X"));
  params.insert(std::make_unique<StringParameter>(k_YAxisTitle_Key, "Y Axis Title", "The label for the Y axis", "Y"));

  params.insertSeparator(Parameters::Separator{"Output Parameter(s)"});

  params.insert(
      std::make_unique<DataGroupCreationParameter>(k_OutputImageGeometryPath_Key, "Output Image Geometry", "The path to the output ImageGeometry that will hold the chart image", DataPath({"Chart Image"})));
  params.insert(std::make_unique<DataObjectNameParameter>(k_CellDataName_Key, "Cell Data Name", "The name of the cell Attribute Matrix to be created", ImageGeom::k_CellAttributeMatrixName));
  params.insert(std::make_unique<DataObjectNameParameter>(k_ImageDataName_Key, "Image Data Name", "The name of the RGBA image data array", "ImageData"));

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType CreateScatterChartFilter::parametersVersion() const
{
  return 1;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer CreateScatterChartFilter::clone() const
{
  return std::make_unique<CreateScatterChartFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult CreateScatterChartFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler,
                                                                 const std::atomic_bool& shouldCancel, const ExecutionContext& executionContext) const
{
  auto xDataPath = filterArgs.value<DataPath>(k_XDataPath_Key);
  auto yDataPath = filterArgs.value<DataPath>(k_YDataPath_Key);
  auto chartDims = filterArgs.value<VectorUInt32Parameter::ValueType>(k_ChartDimensions_Key);
  auto outputImageGeomPath = filterArgs.value<DataPath>(k_OutputImageGeometryPath_Key);
  auto cellDataName = filterArgs.value<DataObjectNameParameter::ValueType>(k_CellDataName_Key);
  auto imageDataName = filterArgs.value<DataObjectNameParameter::ValueType>(k_ImageDataName_Key);

  const auto& xDataArray = dataStructure.getDataRefAs<IDataArray>(xDataPath);
  const auto& yDataArray = dataStructure.getDataRefAs<IDataArray>(yDataPath);

  if(xDataArray.getNumberOfTuples() != yDataArray.getNumberOfTuples())
  {
    return {MakeErrorResult<OutputActions>(-2000, "X Data and Y Data arrays must have the same number of tuples.")};
  }

  uint32 chartHeight = chartDims[0];
  uint32 chartWidth = chartDims[1];

  if(chartHeight == 0 || chartWidth == 0)
  {
    return {MakeErrorResult<OutputActions>(-2001, "Chart dimensions must be greater than zero.")};
  }

  nx::core::Result<OutputActions> resultOutputActions;

  // Create the output ImageGeometry: dims are X (width), Y (height), Z (1 for 2D image)
  auto createImageGeometryAction = std::make_unique<CreateImageGeometryAction>(outputImageGeomPath,
                                                                               CreateImageGeometryAction::DimensionType({chartWidth, chartHeight, 1}),
                                                                               CreateImageGeometryAction::OriginType({0.0F, 0.0F, 0.0F}),
                                                                               CreateImageGeometryAction::SpacingType({1.0F, 1.0F, 1.0F}), cellDataName);
  resultOutputActions.value().appendAction(std::move(createImageGeometryAction));

  // Create the RGB image data array inside the cell attribute matrix
  // Tuple dimensions must match the attribute matrix shape which is reversed: (Z, Y, X)
  DataPath imageDataPath = outputImageGeomPath.createChildPath(cellDataName).createChildPath(imageDataName);
  auto createArrayAction = std::make_unique<CreateArrayAction>(DataType::uint8, std::vector<usize>{1, static_cast<usize>(chartHeight), static_cast<usize>(chartWidth)},
                                                               std::vector<usize>{3}, imageDataPath);
  resultOutputActions.value().appendAction(std::move(createArrayAction));

  return {std::move(resultOutputActions)};
}

//------------------------------------------------------------------------------
Result<> CreateScatterChartFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                                               const std::atomic_bool& shouldCancel, const ExecutionContext& executionContext) const
{
  CreateScatterChartInputValues inputValues;

  inputValues.XDataPath = filterArgs.value<DataPath>(k_XDataPath_Key);
  inputValues.YDataPath = filterArgs.value<DataPath>(k_YDataPath_Key);
  inputValues.ChartDimensions = filterArgs.value<VectorUInt32Parameter::ValueType>(k_ChartDimensions_Key);
  inputValues.ChartTitle = filterArgs.value<StringParameter::ValueType>(k_ChartTitle_Key);
  inputValues.XAxisTitle = filterArgs.value<StringParameter::ValueType>(k_XAxisTitle_Key);
  inputValues.YAxisTitle = filterArgs.value<StringParameter::ValueType>(k_YAxisTitle_Key);
  inputValues.OutputImageGeometryPath = filterArgs.value<DataPath>(k_OutputImageGeometryPath_Key);

  auto cellDataName = filterArgs.value<DataObjectNameParameter::ValueType>(k_CellDataName_Key);
  auto imageDataName = filterArgs.value<DataObjectNameParameter::ValueType>(k_ImageDataName_Key);
  inputValues.ImageDataPath = inputValues.OutputImageGeometryPath.createChildPath(cellDataName).createChildPath(imageDataName);

  return CreateScatterChart(dataStructure, messageHandler, shouldCancel, &inputValues)();
}

} // namespace nx::core
