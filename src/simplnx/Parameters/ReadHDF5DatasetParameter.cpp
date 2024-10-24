/* ============================================================================
 * Copyright (c) 2022-2022 BlueQuartz Software, LLC
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
 * contributors may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "ReadHDF5DatasetParameter.hpp"

#include "simplnx/Utilities/SIMPLConversion.hpp"

#include <nlohmann/json.hpp>

using namespace nx::core;
namespace
{
constexpr StringLiteral k_DataPathsKey = "data_paths";
constexpr StringLiteral k_InputFileKey = "input_file";
constexpr StringLiteral k_ParentGroupKey = "parent_group";
} // namespace

namespace nx::core
{
// -----------------------------------------------------------------------------
Result<ReadHDF5DatasetParameter::DatasetImportInfo> ReadHDF5DatasetParameter::DatasetImportInfo::ReadJson(const nlohmann::json& json)
{
  DatasetImportInfo data;
  if(!json.contains(k_DatasetPath_Key.view()))
  {
    return MakeErrorResult<DatasetImportInfo>(-500, fmt::format("ImportHDF5DatasetParameter ValueType: Cannot find the key \"{}\" in the ValueType json object.", k_DatasetPath_Key));
  }
  else if(!json[k_DatasetPath_Key.str()].is_string())
  {
    return MakeErrorResult<DatasetImportInfo>(-501, fmt::format("ImportHDF5DatasetParameter ValueType: 'Dataset Path' value is of type {} and is not a string.", json[k_DatasetPath_Key].type_name()));
  }
  data.dataSetPath = json[k_DatasetPath_Key.str()];

  if(!json.contains(k_ComponentDimensions_Key.view()))
  {
    return MakeErrorResult<DatasetImportInfo>(-502, fmt::format("ImportHDF5DatasetParameter ValueType: Cannot find the key \"{}\" in the ValueType json object.", k_ComponentDimensions_Key));
  }
  else if(!json[k_ComponentDimensions_Key.str()].is_string())
  {
    return MakeErrorResult<DatasetImportInfo>(
        -503, fmt::format("ImportHDF5DatasetParameter ValueType: 'Component Dimensions' value is of type {} and is not a string.", json[k_ComponentDimensions_Key].type_name()));
  }
  data.componentDimensions = json[k_ComponentDimensions_Key.str()];

  if(!json.contains(k_TupleDimensions_Key.view()))
  {
    return MakeErrorResult<DatasetImportInfo>(-502, fmt::format("ImportHDF5DatasetParameter ValueType: Cannot find the key \"{}\" in the ValueType json object.", k_TupleDimensions_Key));
  }
  else if(!json[k_TupleDimensions_Key.str()].is_string())
  {
    return MakeErrorResult<DatasetImportInfo>(
        -503, fmt::format("ImportHDF5DatasetParameter ValueType: 'Tuple Dimensions' value is of type {} and is not a string.", json[k_TupleDimensions_Key].type_name()));
  }
  data.tupleDimensions = json[k_TupleDimensions_Key.str()];

  return {data};
}

// -----------------------------------------------------------------------------
nlohmann::json ReadHDF5DatasetParameter::DatasetImportInfo::writeJson() const
{
  nlohmann::json json;
  json[k_DatasetPath_Key.str()] = dataSetPath;
  json[k_ComponentDimensions_Key.str()] = componentDimensions;
  json[k_TupleDimensions_Key.str()] = tupleDimensions;
  return json;
}

// -----------------------------------------------------------------------------
ReadHDF5DatasetParameter::ReadHDF5DatasetParameter(const std::string& name, const std::string& humanName, const std::string& helpText, const ValueType& defaultValue)
: ValueParameter(name, humanName, helpText)
, m_DefaultValue(defaultValue)
{
}

// -----------------------------------------------------------------------------
Uuid ReadHDF5DatasetParameter::uuid() const
{
  return ParameterTraits<ReadHDF5DatasetParameter>::uuid;
}

// -----------------------------------------------------------------------------
IParameter::AcceptedTypes ReadHDF5DatasetParameter::acceptedTypes() const
{
  return {typeid(ValueType)};
}

//------------------------------------------------------------------------------
IParameter::VersionType ReadHDF5DatasetParameter::getVersion() const
{
  return 1;
}

// -----------------------------------------------------------------------------
nlohmann::json ReadHDF5DatasetParameter::toJsonImpl(const std::any& value) const
{
  const auto& datasetImportInfo = GetAnyRef<ValueType>(value);
  nlohmann::json json;
  nlohmann::json dataPathsJson = nlohmann::json::array();
  for(const auto& importInfo : datasetImportInfo.datasets)
  {
    dataPathsJson.push_back(importInfo.writeJson());
  }
  json[k_InputFileKey.str()] = datasetImportInfo.inputFile;
  if(datasetImportInfo.parent.has_value())
  {
    json[k_ParentGroupKey.str()] = datasetImportInfo.parent.value().toString();
  }
  json[k_DataPathsKey.str()] = std::move(dataPathsJson);
  return json;
}

// -----------------------------------------------------------------------------
Result<std::any> ReadHDF5DatasetParameter::fromJsonImpl(const nlohmann::json& json, VersionType version) const
{
  static constexpr StringLiteral prefix = "FilterParameter 'ImportHDF5DatasetParameter' JSON Error: ";
  if(!json.is_object())
  {
    return MakeErrorResult<std::any>(-780, fmt::format("{}JSON value for key '{}' is not an object", prefix, name()));
  }

  if(!json.contains(k_DataPathsKey.view()))
  {
    return MakeErrorResult<std::any>(-781, fmt::format("{}JSON does not contain key '{} / {}'", prefix, name(), k_DataPathsKey.view()));
  }

  if(!json.contains(k_InputFileKey.view()))
  {
    return MakeErrorResult<std::any>(-782, fmt::format("{}JSON does not contain key '{} / {}'", prefix, name(), k_InputFileKey.view()));
  }

  ValueType importData;
  importData.inputFile = json[k_InputFileKey.str()];
  if(json.contains(k_ParentGroupKey.view()))
  {
    importData.parent = DataPath::FromString(json[k_ParentGroupKey.str()].get<std::string>());
  }
  const auto& jsonDataPaths = json[k_DataPathsKey.str()];
  if(!jsonDataPaths.is_null())
  {
    if(!jsonDataPaths.is_array())
    {
      return MakeErrorResult<std::any>(-783, fmt::format("{}JSON value for key '{} / {}' is not an array", prefix, name(), k_DataPathsKey));
    }
    std::vector<Error> errors;
    for(const auto& jsonImportInfo : jsonDataPaths)
    {
      auto importInfo = DatasetImportInfo::ReadJson(jsonImportInfo);
      if(importInfo.invalid())
      {
        return {{nonstd::make_unexpected(std::move(importInfo.errors()))}};
      }
      importData.datasets.push_back(std::move(importInfo.value()));
    }
  }

  return {{std::move(importData)}};
}

// -----------------------------------------------------------------------------
IParameter::UniquePointer ReadHDF5DatasetParameter::clone() const
{
  return std::make_unique<ReadHDF5DatasetParameter>(name(), humanName(), helpText(), m_DefaultValue);
}

// -----------------------------------------------------------------------------
std::any ReadHDF5DatasetParameter::defaultValue() const
{
  return m_DefaultValue;
}

// -----------------------------------------------------------------------------
Result<> ReadHDF5DatasetParameter::validate(const std::any& value) const
{
  [[maybe_unused]] auto data = std::any_cast<ValueType>(value);
  return {};
}

namespace SIMPLConversion
{
namespace
{
constexpr StringLiteral k_DatasetPathKey = "dataset_path";
constexpr StringLiteral k_ComponentDimensionsKey = "component_dimensions";
constexpr StringLiteral k_TupleDimensionsKey = "tuple_dimensions";
} // namespace

Result<ImportHDF5DatasetFilterParameterConverter::ValueType> ImportHDF5DatasetFilterParameterConverter::convert(const nlohmann::json& json1, const nlohmann::json& json2, const nlohmann::json& json3)
{
  if(!json1.is_array())
  {
    return MakeErrorResult<ValueType>(-1, fmt::format("ReadHDF5DatasetParameter json '{}' is not an array", json1.dump()));
  }

  if(!json2.is_string())
  {
    return MakeErrorResult<ValueType>(-3, fmt::format("ImportHDF5DatasetFilterParameter json '{}' is not a string", json2.dump()));
  }

  auto dataContainerNameResult = ReadDataContainerName(json3, "ImportHDF5DatasetFilterParameter");
  if(dataContainerNameResult.invalid())
  {
    return ConvertInvalidResult<ValueType>(std::move(dataContainerNameResult));
  }

  auto attributeMatrixNameResult = ReadAttributeMatrixName(json3, "ImportHDF5DatasetFilterParameter");
  if(attributeMatrixNameResult.invalid())
  {
    return ConvertInvalidResult<ValueType>(std::move(attributeMatrixNameResult));
  }

  ValueType value;
  value.inputFile = json2.get<std::string>();
  value.parent = DataPath({std::move(dataContainerNameResult.value()), std::move(attributeMatrixNameResult.value())});

  for(const auto& iter : json1)
  {
    if(!iter.is_object())
    {
      return MakeErrorResult<ValueType>(-2, fmt::format("ImportHDF5DatasetFilterParameter json '{}' is not an object", iter.dump()));
    }

    ParameterType::DatasetImportInfo info;
    info.dataSetPath = iter[k_DatasetPathKey].get<std::string>();
    info.componentDimensions = iter[k_ComponentDimensionsKey].get<std::string>();
    value.datasets.push_back(std::move(info));
  }

  return {std::move(value)};
}
} // namespace SIMPLConversion
} // namespace nx::core
