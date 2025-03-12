//
// Created by Michael Jackson on 3/12/25.
//

#include "IStatsData.hpp"

#include "simplnx/DataStructure/DataStructure.hpp"

#include "fmt/format.h"

#include <numeric>
#include <stdexcept>

namespace nx::core
{

StatsDataArray* StatsDataArray::Create(DataStructure& dataStructure, const std::string_view& name, const std::optional<IdType>& parentId)
{
  auto data = std::shared_ptr<StatsDataArray>(new StatsDataArray(dataStructure, name.data()));
  if(!AttemptToAddObject(dataStructure, data, parentId))
  {
    return nullptr;
  }
  return data.get();
}

DataObject* StatsDataArray::shallowCopy()
{
  return nullptr; // new StatsDataArray(*this);
}

std::shared_ptr<DataObject> StatsDataArray::deepCopy(const DataPath& copyPath)
{
  // TODO: THIS IS BUSTED
  auto& dataStruct = getDataStructureRef();
  if(dataStruct.containsData(copyPath))
  {
    return nullptr;
  }
  // Don't construct with identifier since it will get created when inserting into data structure
  const auto copy = std::shared_ptr<StatsDataArray>(new StatsDataArray(dataStruct, copyPath.getTargetName()));
  if(dataStruct.insert(copy, copyPath.getParent()))
  {
    return copy;
  }
  return nullptr;
}

DataObject::Type StatsDataArray::getDataObjectType() const
{
  return DataObject::Type::StatsDataArray;
}

std::string StatsDataArray::getTypeName() const
{
  return k_TypeName;
}

IArray::ArrayType StatsDataArray::getArrayType() const
{
  return ArrayType::StatsDataArray;
}

void StatsDataArray::resizeTuples(const std::vector<usize>& tupleShape)
{
  auto numTuples = std::accumulate(tupleShape.cbegin(), tupleShape.cend(), static_cast<usize>(1), std::multiplies<>());
  if(numTuples != size())
  {
    m_Arrays.resize(numTuples);
  }
}

} // namespace nx::core
