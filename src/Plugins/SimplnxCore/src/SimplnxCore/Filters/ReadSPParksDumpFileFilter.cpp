#include "ReadSPParksDumpFileFilter.hpp"

#include "SimplnxCore/Filters/Algorithms/ReadSPParksDumpFile.hpp"

#include "simplnx/Common/AtomicFile.hpp"
#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/VectorParameter.hpp"
#include "simplnx/Utilities/SIMPLConversion.hpp"

#include <filesystem>

namespace fs = std::filesystem;
using namespace nx::core;
namespace
{
std::atomic_int32_t s_InstanceId = 0;
std::map<int32, SPParksDumpFileCache> s_HeaderCache;
} // namespace

namespace nx::core
{
//------------------------------------------------------------------------------
ReadSPParksDumpFileFilter::ReadSPParksDumpFileFilter()
: m_InstanceId(s_InstanceId.fetch_add(1))
{
  s_HeaderCache[m_InstanceId] = {};
}

//------------------------------------------------------------------------------
ReadSPParksDumpFileFilter::~ReadSPParksDumpFileFilter() noexcept
{
  s_HeaderCache.erase(m_InstanceId);
}

//------------------------------------------------------------------------------
std::string ReadSPParksDumpFileFilter::name() const
{
  return FilterTraits<ReadSPParksDumpFileFilter>::name.str();
}

//------------------------------------------------------------------------------
std::string ReadSPParksDumpFileFilter::className() const
{
  return FilterTraits<ReadSPParksDumpFileFilter>::className;
}

//------------------------------------------------------------------------------
Uuid ReadSPParksDumpFileFilter::uuid() const
{
  return FilterTraits<ReadSPParksDumpFileFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string ReadSPParksDumpFileFilter::humanName() const
{
  return "Read SPParks Dump File";
}

//------------------------------------------------------------------------------
std::vector<std::string> ReadSPParksDumpFileFilter::defaultTags() const
{
  return {className(), "IO", "Input", "Read", "Import", "SSParks"};
}

//------------------------------------------------------------------------------
Parameters ReadSPParksDumpFileFilter::parameters() const
{
  Parameters params;

  // Create the parameter descriptors that are needed for this filter
  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});

  params.insert(std::make_unique<ChoicesParameter>(k_LengthUnit_Key, "Length Unit", "The length unit that will be set into the created image geometry",
                                                   to_underlying(IGeometry::LengthUnit::Micrometer), IGeometry::GetAllLengthUnitStrings()));
  params.insert(std::make_unique<VectorFloat32Parameter>(k_Origin_Key, "Origin (Physical Units)", "Specifies the new origin values in physical units.", std::vector<float32>{0.0f, 0.0f, 0.0f},
                                                         std::vector<std::string>{"X", "Y", "Z"}));

  params.insert(std::make_unique<VectorFloat32Parameter>(k_Spacing_Key, "Spacing (Physical Units)", "Specifies the new spacing values in physical units.", std::vector<float32>{1.0f, 1.0f, 1.0f},
                                                         std::vector<std::string>{"X", "Y", "Z"}));

  params.insert(std::make_unique<BoolParameter>(k_OneBasedArrays_Key, "One based arrays", "Specifies if the origin should be changed", false));

  params.insert(std::make_unique<FileSystemPathParameter>(k_InputFilePath_Key, "Input File", "Input image file", fs::path(""), FileSystemPathParameter::ExtensionsType{{".dump"}},
                                                          FileSystemPathParameter::PathType::InputFile, false));

  params.insertSeparator(Parameters::Separator{"Output Data Object(s)"});
  params.insert(std::make_unique<DataGroupCreationParameter>(k_CreatedImageGeometryPath_Key, "Created Image Geometry", "The path to the created Image Geometry", DataPath({"ImageGeometry"})));
  params.insert(
      std::make_unique<DataObjectNameParameter>(k_CellAttributeMatrixName_Key, "Created Cell Attribute Matrix", "The name of the created cell attribute matrix", ImageGeom::k_CellAttributeMatrixName));
  params.insert(std::make_unique<DataObjectNameParameter>(k_FeatureIdsArrayName_Key, "Created Feature Ids", "", "Feature Ids"));

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType ReadSPParksDumpFileFilter::parametersVersion() const
{
  return 1;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer ReadSPParksDumpFileFilter::clone() const
{
  return std::make_unique<ReadSPParksDumpFileFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult ReadSPParksDumpFileFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler,
                                                                  const std::atomic_bool& shouldCancel, const ExecutionContext& executionContext) const
{
  auto pInputFilePathValue = filterArgs.value<FileSystemPathParameter::ValueType>(k_InputFilePath_Key);
  auto pNewImageGeometryPathValue = filterArgs.value<DataPath>(k_CreatedImageGeometryPath_Key);
  auto pCellAttributeMatrixNameValue = filterArgs.value<std::string>(k_CellAttributeMatrixName_Key);
  auto pFeatureIdsArrayNameValue = filterArgs.value<std::string>(k_FeatureIdsArrayName_Key);

  auto origin = filterArgs.value<VectorFloat32Parameter::ValueType>(k_Origin_Key);
  auto spacing = filterArgs.value<VectorFloat32Parameter::ValueType>(k_Spacing_Key);
  auto oneBaseArrays = filterArgs.value<bool>(k_OneBasedArrays_Key);

  nx::core::Result<OutputActions> resultOutputActions;
  std::vector<PreflightValue> preflightUpdatedValues;

  if(pInputFilePathValue != s_HeaderCache[m_InstanceId].inputFile || s_HeaderCache[m_InstanceId].timeStamp < fs::last_write_time(pInputFilePathValue))
  {
    // Read header portion of SPParks dump file and store into the cache
    Result<SPParksDumpFileCache> result = ReadSPParksDumpFileHeader(pInputFilePathValue, oneBaseArrays);
    if(result.invalid())
    {
      return {ConvertInvalidResult<OutputActions>(std::move(result))};
    }

    // Cache the results from algorithm run
    s_HeaderCache[m_InstanceId] = result.value();
  }

  SPParksDumpFileCache& fileCache = s_HeaderCache[m_InstanceId];

  {
    auto createImageGeometryAction = std::make_unique<CreateImageGeometryAction>(pNewImageGeometryPathValue, fileCache.dimensions, origin, spacing, pCellAttributeMatrixNameValue);
    resultOutputActions.value().appendAction(std::move(createImageGeometryAction));
  }

  {
    std::ranges::reverse(fileCache.dimensions);
    auto createAction = std::make_unique<CreateArrayAction>(DataType::int32, fileCache.dimensions, std::vector<usize>{1},
                                                            pNewImageGeometryPathValue.createChildPath(pCellAttributeMatrixNameValue).createChildPath(pFeatureIdsArrayNameValue));
    resultOutputActions.value().appendAction(std::move(createAction));
  }
  // Return both the resultOutputActions and the preflightUpdatedValues via std::move()
  return {std::move(resultOutputActions), std::move(preflightUpdatedValues)};
}

//------------------------------------------------------------------------------
Result<> ReadSPParksDumpFileFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                                                const std::atomic_bool& shouldCancel, const ExecutionContext& executionContext) const
{

  return {};
}

namespace
{
namespace SIMPL
{
constexpr StringLiteral k_OutputFileKey = "OutputFile";
constexpr StringLiteral k_FeatureIdsArrayPathKey = "FeatureIdsArrayPath";

} // namespace SIMPL
} // namespace

Result<Arguments> ReadSPParksDumpFileFilter::FromSIMPLJson(const nlohmann::json& json)
{
  Arguments args = ReadSPParksDumpFileFilter().getDefaultArguments();

  std::vector<Result<>> results;

  Result<> conversionResult = MergeResults(std::move(results));

  return ConvertResultTo<Arguments>(std::move(conversionResult), std::move(args));
}
} // namespace nx::core
