#include "ChangeAngleRepresentationFilter.hpp"

#include "simplnx/Common/Numbers.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Utilities/ParallelDataAlgorithm.hpp"

#include <fmt/format.h>

#include "simplnx/Utilities/SIMPLConversion.hpp"

#include <cmath>

using namespace nx::core;

namespace
{

namespace EulerAngleConversionType
{
constexpr uint64 DegreesToRadians = 0;
constexpr uint64 RadiansToDegrees = 1;
} // namespace EulerAngleConversionType

class ChangeAngleRepresentationImpl
{
public:
  ChangeAngleRepresentationImpl(Float32AbstractDataStore& angles, float factor)
  : m_Angles(angles)
  , m_ConvFactor(factor)
  {
  }
  ~ChangeAngleRepresentationImpl() noexcept = default;

  void convert(size_t start, size_t end) const
  {
    for(size_t i = start; i < end; i++)
    {
      m_Angles[i] = m_Angles[i] * m_ConvFactor;
    }
  }

  void operator()(const Range& range) const
  {
    convert(range.min(), range.max());
  }

private:
  Float32AbstractDataStore& m_Angles;
  float32 m_ConvFactor = 0.0F;
};
} // namespace

namespace nx::core
{
//------------------------------------------------------------------------------
std::string ChangeAngleRepresentationFilter::name() const
{
  return FilterTraits<ChangeAngleRepresentationFilter>::name.str();
}

//------------------------------------------------------------------------------
std::string ChangeAngleRepresentationFilter::className() const
{
  return FilterTraits<ChangeAngleRepresentationFilter>::className;
}

//------------------------------------------------------------------------------
Uuid ChangeAngleRepresentationFilter::uuid() const
{
  return FilterTraits<ChangeAngleRepresentationFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string ChangeAngleRepresentationFilter::humanName() const
{
  return "Convert Angles to Degrees or Radians";
}

//------------------------------------------------------------------------------
std::vector<std::string> ChangeAngleRepresentationFilter::defaultTags() const
{
  return {className(), "Processing", "Conversion"};
}

//------------------------------------------------------------------------------
Parameters ChangeAngleRepresentationFilter::parameters() const
{
  Parameters params;

  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});
  params.insert(std::make_unique<ChoicesParameter>(k_ConversionType_Key, "Conversion Type", "Tells the Filter which conversion is being made", 0,
                                                   ChoicesParameter::Choices{"Degrees to Radians", "Radians to Degrees"}));
  params.insert(std::make_unique<ArraySelectionParameter>(k_AnglesArrayPath_Key, "Angles", "The DataArray containing the angles to be converted.", DataPath{},
                                                          ArraySelectionParameter::AllowedTypes{DataType::float32}));

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType ChangeAngleRepresentationFilter::parametersVersion() const
{
  return 1;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer ChangeAngleRepresentationFilter::clone() const
{
  return std::make_unique<ChangeAngleRepresentationFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult ChangeAngleRepresentationFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler,
                                                                        const std::atomic_bool& shouldCancel) const
{
  auto pConversionTypeValue = filterArgs.value<ChoicesParameter::ValueType>(k_ConversionType_Key);

  if(pConversionTypeValue > 1)
  {
    return {MakeErrorResult<OutputActions>(-67001, fmt::format("The conversion type must be either [0|1]. Value given is '{}'", pConversionTypeValue))};
  }

  return {};
}

//------------------------------------------------------------------------------
Result<> ChangeAngleRepresentationFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                                                      const std::atomic_bool& shouldCancel) const
{
  /****************************************************************************
   * Extract the actual input values from the 'filterArgs' object
   ***************************************************************************/
  auto pConversionTypeValue = filterArgs.value<ChoicesParameter::ValueType>(k_ConversionType_Key);
  auto pAnglesDataPathValue = filterArgs.value<DataPath>(k_AnglesArrayPath_Key);

  auto& angles = dataStructure.getDataAs<Float32Array>(pAnglesDataPathValue)->getDataStoreRef();

  float conversionFactor = 1.0f;
  if(pConversionTypeValue == EulerAngleConversionType::DegreesToRadians)
  {
    conversionFactor = static_cast<float>(nx::core::numbers::pi / 180.0f);
  }
  else if(pConversionTypeValue == EulerAngleConversionType::RadiansToDegrees)
  {
    conversionFactor = static_cast<float>(180.0f / nx::core::numbers::pi);
  }

  // Parallel algorithm to find duplicate nodes
  ParallelDataAlgorithm dataAlg;
  dataAlg.setRange(0ULL, angles.getSize());
  dataAlg.execute(::ChangeAngleRepresentationImpl(angles, conversionFactor));

  return {};
}
} // namespace nx::core

namespace
{
namespace SIMPL
{
constexpr StringLiteral k_ConversionTypeKey = "ConversionType";
constexpr StringLiteral k_CellEulerAnglesArrayPathKey = "CellEulerAnglesArrayPath";
} // namespace SIMPL
} // namespace

Result<Arguments> ChangeAngleRepresentationFilter::FromSIMPLJson(const nlohmann::json& json)
{
  Arguments args = ChangeAngleRepresentationFilter().getDefaultArguments();

  std::vector<Result<>> results;

  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::ChoiceFilterParameterConverter>(args, json, SIMPL::k_ConversionTypeKey, k_ConversionType_Key));
  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::DataArraySelectionFilterParameterConverter>(args, json, SIMPL::k_CellEulerAnglesArrayPathKey, k_AnglesArrayPath_Key));

  Result<> conversionResult = MergeResults(std::move(results));

  return ConvertResultTo<Arguments>(std::move(conversionResult), std::move(args));
}
