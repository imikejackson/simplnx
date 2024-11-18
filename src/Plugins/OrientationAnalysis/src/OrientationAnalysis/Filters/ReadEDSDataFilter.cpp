#include "ReadEDSDataFilter.hpp"

#include "H5Support/H5Lite.h"
#include "H5Support/H5ScopedSentinel.h"
#include "H5Support/H5Utilities.h"
#include "OrientationAnalysis/Filters/Algorithms/ReadAngData.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateAttributeMatrixAction.hpp"
#include "simplnx/Filter/Actions/CreateDataGroupAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Filter/Actions/CreateStringArrayAction.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"

#include "simplnx/Utilities/SIMPLConversion.hpp"
#include "simplnx/Utilities/StringUtilities.hpp"

#include <filesystem>

namespace fs = std::filesystem;

using namespace H5Support;
// /AMCW_scans/scans/AMCW_0006_2024-08-09_T15-15-14_sample0_centroid_hfw=153.6um_manual-post-run_mainEBSD+EDS/AMCW_0006_2024-08-09_T15-15-14_sample0_centroid_hfw=153.6um_manual-post-run_mainEBSD+EDS/EDS/
namespace
{

nx::core::DataType getNumericType(hid_t dataTypeIdentifier)
{
  H5SUPPORT_MUTEX_LOCK()

  //  if(dataTypeIdentifier == H5T_STRING)
  //  {
  //    return nx::core::NumericType::;
  //  }

  if(H5Tequal(dataTypeIdentifier, H5T_NATIVE_INT8) > 0)
  {
    return nx::core::DataType::int8;
  }
  if(H5Tequal(dataTypeIdentifier, H5T_NATIVE_UINT8) > 0)
  {
    return nx::core::DataType::uint8;
  }

  if(H5Tequal(dataTypeIdentifier, H5T_NATIVE_INT16) > 0)
  {
    return nx::core::DataType::int16;
  }
  if(H5Tequal(dataTypeIdentifier, H5T_NATIVE_UINT16) > 0)
  {
    return nx::core::DataType::uint16;
  }

  if(H5Tequal(dataTypeIdentifier, H5T_NATIVE_INT32) > 0)
  {
    return nx::core::DataType::int32;
  }
  if(H5Tequal(dataTypeIdentifier, H5T_NATIVE_UINT32) > 0)
  {
    return nx::core::DataType::uint32;
  }

  if(H5Tequal(dataTypeIdentifier, H5T_NATIVE_INT64) > 0)
  {
    return nx::core::DataType::int64;
  }
  if(H5Tequal(dataTypeIdentifier, H5T_NATIVE_UINT64) > 0)
  {
    return nx::core::DataType::uint64;
  }

  if(H5Tequal(dataTypeIdentifier, H5T_NATIVE_FLOAT) > 0)
  {
    return nx::core::DataType::float32;
  }
  if(H5Tequal(dataTypeIdentifier, H5T_NATIVE_DOUBLE) > 0)
  {
    return nx::core::DataType::float64;
  }

  return nx::core::DataType::boolean;
}

void CreateSPCDataArrayActions(const std::string& filePath, const std::string& scanGroupPath, const nx::core::DataPath& attrMatPath, nx::core::Result<nx::core::OutputActions>& resultOutputActions)
{
  std::vector<size_t> tupleDims = {1};

  std::vector<std::string> pathElements = nx::core::StringUtilities::split(scanGroupPath, '/');
  std::string projName = pathElements[0];
  std::string sampleName = pathElements[1];
  std::string areaName = pathElements[2];
  std::string scanName = pathElements[3];

  hid_t fileId = H5Utilities::openFile(filePath, true);
  H5ScopedFileSentinel fileSentinel(fileId, false);

  hid_t scanGroupId = H5Utilities::openHDF5Object(fileId, scanGroupPath);
  fileSentinel.addGroupId(scanGroupId);
  if(scanGroupId <= 0)
  {
    return;
  }

  hid_t edsGroupId = H5Utilities::openHDF5Object(scanGroupId, "EDS");
  fileSentinel.addGroupId(edsGroupId);
  if(edsGroupId <= 0)
  {
    return;
  }

  // Open theSPC dataset
  hid_t dataset_id = H5Dopen2(edsGroupId, "SPC", H5P_DEFAULT);
  H5ScopedObjectSentinel dataset_sentinel(dataset_id, false);

  // Get the datatype of the dataset
  hid_t full_datatype = H5Dget_type(dataset_id);

  // Get the datatype of the dataset
  // Iterate over each member and print its name and type
  std::vector<size_t> cDims = {1ULL};
  int n_members = H5Tget_nmembers(full_datatype);
  for(int field_idx = 0; field_idx < n_members; field_idx++)
  {
    // Get member name
    char* namePtr = H5Tget_member_name(full_datatype, field_idx);
    std::string name(namePtr);
    // Get member datatype and class
    hid_t member_type = H5Tget_member_type(full_datatype, field_idx);
    H5T_class_t type_class = H5Tget_class(member_type);

    nx::core::DataPath dataArrayPath = attrMatPath.createChildPath(name);
    if(type_class == H5T_INTEGER || type_class == H5T_FLOAT)
    {
      auto action = std::make_unique<nx::core::CreateArrayAction>(getNumericType(member_type), tupleDims, cDims, dataArrayPath);
      resultOutputActions.value().appendAction(std::move(action));
    }
    else if(type_class == H5T_STRING)
    {
      std::cout << field_idx << "  H5T_STRING with Name: " << namePtr << std::endl;
    }
    else if(type_class == H5T_ARRAY)
    {
      int ndims = H5Tget_array_ndims(member_type);
      std::vector<hsize_t> c_dims(ndims);
      H5Tget_array_dims(member_type, c_dims.data());
      // std::cout << field_idx << "  H5T_ARRAY with Name: " << namePtr <<  "  nDims: " << c_dims[0] << std::endl;
      auto action = std::make_unique<nx::core::CreateArrayAction>(getNumericType(member_type), tupleDims, std::vector<size_t>{static_cast<size_t>(c_dims[0])}, dataArrayPath);
      resultOutputActions.value().appendAction(std::move(action));
    }
    // Free the member name and close member type
    H5Tclose(member_type);
    free(namePtr);
  }
  herr_t returnError;
  CloseH5T(full_datatype, full_datatype_error, returnError);
}

std::vector<size_t> GetEdsScanDimensions(const std::string& filePath, const std::string& edsGroupPath)
{

  return {0, 0, 0};
}

std::vector<float> GetEdsScanSpacing(const std::string& filePath, const std::string& edsGroupPath)
{
  return {0.0f, 0.0f, 0.0f};
}

} // namespace

namespace nx::core
{
//------------------------------------------------------------------------------
std::string ReadEDSDataFilter::name() const
{
  return FilterTraits<ReadEDSDataFilter>::name.str();
}

//------------------------------------------------------------------------------
std::string ReadEDSDataFilter::className() const
{
  return FilterTraits<ReadEDSDataFilter>::className;
}

//------------------------------------------------------------------------------
Uuid ReadEDSDataFilter::uuid() const
{
  return FilterTraits<ReadEDSDataFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string ReadEDSDataFilter::humanName() const
{
  return "Read EDS Data (.edaxh5)";
}

//------------------------------------------------------------------------------
std::vector<std::string> ReadEDSDataFilter::defaultTags() const
{
  return {className(), "IO", "Input", "Read", "Import", "EDAX", "EDAXH5", "EDS"};
}

//------------------------------------------------------------------------------
Parameters ReadEDSDataFilter::parameters() const
{
  Parameters params;
  // Create the parameter descriptors that are needed for this filter
  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});
  params.insert(std::make_unique<FileSystemPathParameter>(k_InputFile_Key, "Input File", "The input .edaxh5 file path", fs::path(""), FileSystemPathParameter::ExtensionsType{".edaxh5"},
                                                          FileSystemPathParameter::PathType::InputFile));
  params.insert(std::make_unique<StringParameter>(k_Hdf5GroupPath_Key, "HDF5 Path to EDS Group", "", ""));

  params.insertSeparator(Parameters::Separator{"Output Image Geometry"});
  params.insert(std::make_unique<DataGroupCreationParameter>(k_CreatedImageGeometryPath_Key, "Image Geometry", "The path to the created Image Geometry", DataPath({ImageGeom::k_TypeName})));
  params.insertSeparator(Parameters::Separator{"Output Cell Attribute Matrix"});
  params.insert(std::make_unique<DataObjectNameParameter>(k_CellAttributeMatrixName_Key, "Cell Attribute Matrix", "The name of the cell data attribute matrix for the created Image Geometry",
                                                          ImageGeom::k_CellDataName));
  params.insertSeparator(Parameters::Separator{"Output Ensemble Attribute Matrix"});
  params.insert(std::make_unique<DataObjectNameParameter>(k_CellEnsembleAttributeMatrixName_Key, "SPC Meta Data", "The Attribute Matrix where the SPC Meta Data information is stored.", "SPC Data"));

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType ReadEDSDataFilter::parametersVersion() const
{
  return 1;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer ReadEDSDataFilter::clone() const
{
  return std::make_unique<ReadEDSDataFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult ReadEDSDataFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler,
                                                          const std::atomic_bool& shouldCancel) const
{
  auto pInputFileValue = filterArgs.value<FileSystemPathParameter::ValueType>(k_InputFile_Key);
  auto pImageGeometryPath = filterArgs.value<DataPath>(k_CreatedImageGeometryPath_Key);
  auto pCellAttributeMatrixNameValue = filterArgs.value<std::string>(k_CellAttributeMatrixName_Key);
  auto pCellEnsembleAttributeMatrixNameValue = filterArgs.value<std::string>(k_CellEnsembleAttributeMatrixName_Key);
  auto pEdsHdfPath = filterArgs.value<std::string>(k_Hdf5GroupPath_Key);

  PreflightResult preflightResult;

  CreateImageGeometryAction::DimensionType imageGeomDims = GetEdsScanDimensions(pInputFileValue, pEdsHdfPath);
  std::vector<size_t> tupleDims = {imageGeomDims[2], imageGeomDims[1], imageGeomDims[0]};

  CreateImageGeometryAction::SpacingType spacing = GetEdsScanSpacing(pInputFileValue, pEdsHdfPath);
  CreateImageGeometryAction::OriginType origin = {0.0F, 0.0F, 0.0F};

  // These variables should be updated with the latest data generated for each variable during preflight.
  // These will be returned through the preflightResult variable to the
  // user interface. You could make these member variables instead if needed.
  std::stringstream ss;
  std::array<float, 3> halfRes = {spacing[0] * 0.5F, spacing[1] * 0.5F, spacing[2] * 0.5F};

  //  ss << "Grid: " << reader.getGrid() << "\n"
  //     << "X Step: " << reader.getXStep() << "    Y Step: " << reader.getYStep() << "\n"
  //     << "Num Odd Cols: " << reader.getNumOddCols() << "    "
  //     << "Num Even Cols: " << reader.getNumEvenCols() << "    "
  //     << "Num Rows: " << reader.getNumRows() << "\n"
  //     << "Sample Physical Dimensions: " << (reader.getXStep() * reader.getNumOddCols()) << " (W) x " << (reader.getYStep() * reader.getNumRows()) << " (H) microns"
  //     << "\n";
  std::string fileInfo = ss.str();

  std::vector<PreflightValue> preflightUpdatedValues = {{"Edaxh5 File Information", fileInfo}};

  // Define a custom class that generates the changes to the DataStructure.
  auto createImageGeometryAction = std::make_unique<CreateImageGeometryAction>(pImageGeometryPath, CreateImageGeometryAction::DimensionType({imageGeomDims[0], imageGeomDims[1], imageGeomDims[2]}),
                                                                               origin, spacing, pCellAttributeMatrixNameValue, IGeometry::LengthUnit::Micrometer);

  // Assign the createImageGeometryAction to the Result<OutputActions>::actions vector via a push_back
  nx::core::Result<OutputActions> resultOutputActions;
  resultOutputActions.value().appendAction(std::move(createImageGeometryAction));

  DataPath cellAttributeMatrixPath = pImageGeometryPath.createChildPath(pCellAttributeMatrixNameValue);

  // Create the Ensemble AttributeMatrix
  tupleDims = {1};
  DataPath ensembleAttributeMatrixPath = pImageGeometryPath.createChildPath(pCellEnsembleAttributeMatrixNameValue);
  {
    auto createAttributeMatrixAction = std::make_unique<CreateAttributeMatrixAction>(ensembleAttributeMatrixPath, tupleDims);
    resultOutputActions.value().appendAction(std::move(createAttributeMatrixAction));
  }

  CreateSPCDataArrayActions(pInputFileValue, pEdsHdfPath, ensembleAttributeMatrixPath, resultOutputActions);

  // Return both the resultOutputActions and the preflightUpdatedValues via std::move()
  return {std::move(resultOutputActions), std::move(preflightUpdatedValues)};
}

//------------------------------------------------------------------------------
Result<> ReadEDSDataFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                                        const std::atomic_bool& shouldCancel) const
{
  ReadAngDataInputValues inputValues;

  inputValues.InputFile = filterArgs.value<FileSystemPathParameter::ValueType>(k_InputFile_Key);
  inputValues.DataContainerName = filterArgs.value<DataPath>(k_CreatedImageGeometryPath_Key);
  inputValues.CellAttributeMatrixName = filterArgs.value<std::string>(k_CellAttributeMatrixName_Key);
  inputValues.CellEnsembleAttributeMatrixName = filterArgs.value<std::string>(k_CellEnsembleAttributeMatrixName_Key);

  ReadAngData readAngData(dataStructure, messageHandler, shouldCancel, &inputValues);
  return readAngData();
}

namespace
{
namespace SIMPL
{
constexpr StringLiteral k_InputFileKey = "InputFile";
constexpr StringLiteral k_DataContainerNameKey = "DataContainerName";
constexpr StringLiteral k_CellAttributeMatrixNameKey = "CellAttributeMatrixName";
constexpr StringLiteral k_CellEnsembleAttributeMatrixNameKey = "CellEnsembleAttributeMatrixName";
} // namespace SIMPL
} // namespace

Result<Arguments> ReadEDSDataFilter::FromSIMPLJson(const nlohmann::json& json)
{
  Arguments args = ReadEDSDataFilter().getDefaultArguments();

  std::vector<Result<>> results;

  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::InputFileFilterParameterConverter>(args, json, SIMPL::k_InputFileKey, k_InputFile_Key));
  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::DataContainerCreationFilterParameterConverter>(args, json, SIMPL::k_DataContainerNameKey, k_CreatedImageGeometryPath_Key));
  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::LinkedPathCreationFilterParameterConverter>(args, json, SIMPL::k_CellAttributeMatrixNameKey, k_CellAttributeMatrixName_Key));
  results.push_back(
      SIMPLConversion::ConvertParameter<SIMPLConversion::LinkedPathCreationFilterParameterConverter>(args, json, SIMPL::k_CellEnsembleAttributeMatrixNameKey, k_CellEnsembleAttributeMatrixName_Key));

  Result<> conversionResult = MergeResults(std::move(results));

  return ConvertResultTo<Arguments>(std::move(conversionResult), std::move(args));
}
} // namespace nx::core
