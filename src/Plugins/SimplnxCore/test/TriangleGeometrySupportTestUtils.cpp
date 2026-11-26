#include "TriangleGeometrySupportTestUtils.hpp"

#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <chrono>

namespace fs = std::filesystem;

namespace nx::core::TriangleGeometrySupportTest
{
namespace
{
DataPath SupportDataPath(SupportKind kind)
{
  return k_MeshGeomPath.createChildPath(kind == SupportKind::FaceNeighbors ? TriangleGeom::k_EltNeighbors : TriangleGeom::k_EltsContainingVert);
}
} // namespace

SupportSnapshot CaptureSupportSnapshot(const DataStructure& dataStructure, SupportKind kind)
{
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<TriangleGeom>(k_MeshGeomPath));
  const auto& geometry = dataStructure.getDataRefAs<TriangleGeom>(k_MeshGeomPath);
  SupportSnapshot snapshot;
  snapshot.association = kind == SupportKind::FaceNeighbors ? geometry.getElementNeighborsId() : geometry.getElementContainingVertId();
  const auto object = dataStructure.getSharedDataAs<SupportListType>(SupportDataPath(kind));
  if(object != nullptr)
  {
    snapshot.object = object;
    snapshot.id = object->getId();
    snapshot.paths = object->getDataPaths();
    for(usize i = 0; i < object->size(); ++i)
    {
      std::vector<SupportIndexType> list;
      list.reserve(object->getNumberOfElements(i));
      for(usize j = 0; j < object->getNumberOfElements(i); ++j)
      {
        list.push_back(object->getElementListPointer(i)[j]);
      }
      snapshot.lists.push_back(std::move(list));
    }
  }
  return snapshot;
}

void RequireSupportSnapshot(const DataStructure& dataStructure, SupportKind kind, const SupportSnapshot& expected)
{
  const auto actual = CaptureSupportSnapshot(dataStructure, kind);
  CHECK(actual.association == expected.association);
  REQUIRE(actual.id == expected.id);
  CHECK(actual.object.lock() == expected.object.lock());
  CHECK(actual.paths == expected.paths);
  CHECK(actual.lists == expected.lists);
}

void CreateSupportData(DataStructure& dataStructure, SupportKind kind, const std::vector<std::vector<SupportIndexType>>& lists, bool associate)
{
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<TriangleGeom>(k_MeshGeomPath));
  auto& geometry = dataStructure.getDataRefAs<TriangleGeom>(k_MeshGeomPath);
  auto* supportPtr = SupportListType::Create(dataStructure, kind == SupportKind::FaceNeighbors ? TriangleGeom::k_EltNeighbors : TriangleGeom::k_EltsContainingVert, geometry.getId());
  REQUIRE(supportPtr != nullptr);
  supportPtr->allocateLists(std::vector<uint16>(lists.size(), 0));
  for(usize i = 0; i < lists.size(); ++i)
  {
    auto list = lists[i];
    if(!list.empty())
    {
      REQUIRE(supportPtr->setElementList(i, static_cast<uint16>(list.size()), list.data()));
    }
  }
  if(associate)
  {
    if(kind == SupportKind::FaceNeighbors)
    {
      geometry.setElementNeighborsId(supportPtr->getId());
    }
    else
    {
      geometry.setElementContainingVertId(supportPtr->getId());
    }
  }
  REQUIRE(CaptureSupportSnapshot(dataStructure, kind).lists == lists);
}

void ApplySupportFixtureActions(DataStructure& dataStructure, const IFilter& filter, const Arguments& args)
{
  auto preflight = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflight.output_actions);
  const auto applied = preflight.output_actions.value().applyAll(dataStructure, IDataAction::Mode::Execute);
  SIMPLNX_RESULT_REQUIRE_VALID(applied);
}

void RequireInjectedFailure(const Result<>& result, const std::string& message)
{
  SIMPLNX_RESULT_REQUIRE_INVALID(result);
  REQUIRE(result.errors().size() == 1);
  CHECK(result.errors()[0].code == k_InjectedError);
  CHECK(result.errors()[0].message == message);
}

void RequireSupportFixtureOocStores(const DataStructure& dataStructure, const std::vector<DataPath>& paths)
{
  for(const auto& path : paths)
  {
    CAPTURE(path);
    const auto* arrayPtr = dataStructure.getDataAs<IDataArray>(path);
    REQUIRE(arrayPtr != nullptr);
    REQUIRE(arrayPtr->getDataFormat() == "HDF5-OOC");
    REQUIRE(arrayPtr->getIDataStoreRef().getStoreType() == IDataStore::StoreType::OutOfCore);
  }
}

SupportSaveFile::SupportSaveFile(const std::string& label)
: path(fs::temp_directory_path() / fmt::format("support-ownership-{}-{}-0.dream3d", label, std::chrono::steady_clock::now().time_since_epoch().count()))
{
  REQUIRE_FALSE(fs::exists(path));
}

SupportSaveFile::~SupportSaveFile() noexcept
{
  std::error_code error;
  fs::remove(path, error);
}
} // namespace nx::core::TriangleGeometrySupportTest
