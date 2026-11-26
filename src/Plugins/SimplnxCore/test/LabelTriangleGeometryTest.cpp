#include "TriangleGeometrySupportTestUtils.hpp"

#include "SimplnxCore/Filters/Algorithms/LabelTriangleGeometry.hpp"
#include "SimplnxCore/Filters/LabelTriangleGeometryFilter.hpp"
#include "SimplnxCore/SimplnxCore_test_dirs.hpp"

#include "simplnx/Core/Application.hpp"
#include "simplnx/DataStructure/Geometry/TriangleGeom.hpp"
#include "simplnx/DataStructure/IO/HDF5/DataStructureReader.hpp"
#include "simplnx/DataStructure/IO/HDF5/DataStructureWriter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <catch2/catch.hpp>

#include <array>
#include <filesystem>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

using namespace nx::core;
using namespace nx::core::TriangleGeometrySupportTest;
namespace fs = std::filesystem;

namespace
{
fs::path KExemplarDataFilePath = fs::path(fmt::format("{}/label_triangle_geometry_test_v2/label_triangle_geometry_test.dream3d", nx::core::unit_test::k_TestFilesDir));
fs::path KBaseDataFilePath = fs::path(fmt::format("{}/label_triangle_geometry_test_v2/data_to_generate_test/combined_stls.dream3d", nx::core::unit_test::k_TestFilesDir));

constexpr StringLiteralType k_CreatedAMName = "Cell Feature AM";
constexpr StringLiteralType k_NumTrianglesName = "NumTriangles";
constexpr StringLiteralType k_RegionIdsName = "Region Ids";

const DataPath k_TriangleGeomPath({"TriangleGeometry"});
const DataPath k_CreatedRegionIdsPath = k_TriangleGeomPath.createChildPath(constants::k_FaceDataWithSpaces).createChildPath(k_RegionIdsName);
const DataPath k_CellFeatureAMPath = k_TriangleGeomPath.createChildPath(k_CreatedAMName);
const DataPath k_NumTrianglesPath = k_CellFeatureAMPath.createChildPath(k_NumTrianglesName);
} // namespace

TEST_CASE("SimplnxCore::LabelTriangleGeometryFilter: Valid Filter Execution", "[SimplnxCore][LabelTriangleGeometryFilter]")
{
  unit_test::LoadPlugins();

  const nx::core::unit_test::TestFileSentinel testDataSentinel(nx::core::unit_test::k_TestFilesDir, "label_triangle_geometry_test_v2.tar.gz", "label_triangle_geometry_test_v2", true, true);

  DataStructure dataStructure = unit_test::LoadDataStructure(KBaseDataFilePath);

  {
    // Instantiate the filter, a DataStructure object and an Arguments Object
    LabelTriangleGeometryFilter filter;
    Arguments args;

    // Create default Parameters for the filter.
    args.insertOrAssign(LabelTriangleGeometryFilter::k_TriangleGeomPath_Key, std::make_any<DataPath>(::k_TriangleGeomPath));
    args.insertOrAssign(LabelTriangleGeometryFilter::k_CreatedRegionIdsPath_Key, std::make_any<DataPath>(::k_CreatedRegionIdsPath));
    args.insertOrAssign(LabelTriangleGeometryFilter::k_TriangleAttributeMatrixName_Key, std::make_any<std::string>(::k_CreatedAMName));
    args.insertOrAssign(LabelTriangleGeometryFilter::k_NumTrianglesName_Key, std::make_any<std::string>(::k_NumTrianglesName));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);

    // This is in here because the exemplar face attribute matrix is not sized correctly. This will
    // correct that value allowing the test to proceed normally.
    auto& exemplarContourTriGeom = dataStructure.getDataRefAs<TriangleGeom>(k_TriangleGeomPath);
    auto vertexAttributeMatrixResizeResult = exemplarContourTriGeom.getVertexAttributeMatrix()->resizeTuples({144});
    SIMPLNX_RESULT_REQUIRE_VALID(vertexAttributeMatrixResizeResult);

    // Execute the filter and check the result
    auto executeResult = filter.execute(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);
  }

  DataStructure exemplarDataStructure = unit_test::LoadDataStructure(KExemplarDataFilePath);

  unit_test::CompareExemplarToGeneratedData(dataStructure, exemplarDataStructure, ::k_CellFeatureAMPath, ::k_TriangleGeomPath.toString());
  unit_test::CompareExemplarToGeneratedData(dataStructure, exemplarDataStructure, ::k_TriangleGeomPath.createChildPath(constants::k_FaceDataWithSpaces), ::k_TriangleGeomPath.toString());

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::LabelTriangleGeometryFilter: SIMPL Backwards Compatibility", "[SimplnxCore][LabelTriangleGeometryFilter][BackwardsCompatibility]")
{
  const fs::path fixturePath = fs::path(nx::core::unit_test::k_SourceDir.view()) / "test" / "simpl_conversion" / "6_5" / "LabelTriangleGeometryFilter.json";
  const auto args = unit_test::LoadSIMPLConversionArguments(fixturePath, FilterTraits<LabelTriangleGeometryFilter>::uuid);
  CHECK(args.value<DataPath>(LabelTriangleGeometryFilter::k_TriangleGeomPath_Key) == DataPath({"DataContainer"}));
  CHECK(args.value<DataPath>(LabelTriangleGeometryFilter::k_CreatedRegionIdsPath_Key) == DataPath({"DataContainer", "FaceData", "RegionIds"}));
  CHECK(args.value<std::string>(LabelTriangleGeometryFilter::k_TriangleAttributeMatrixName_Key) == "TriangleData");
  CHECK(args.value<std::string>(LabelTriangleGeometryFilter::k_NumTrianglesName_Key) == "NumTriangles");
}

namespace
{
const DataPath k_RegionIdsPath = k_MeshFacePath.createChildPath("RegionIds");
const DataPath k_MeshFeaturePath = k_MeshGeomPath.createChildPath("Features");
const DataPath k_TriangleCountPath = k_MeshFeaturePath.createChildPath("NumTriangles");

Arguments CreateLabelFixture(DataStructure& dataStructure, bool incidence, bool neighbors)
{
  auto* geometryPtr = TriangleGeom::Create(dataStructure, "Triangles");
  REQUIRE(geometryPtr != nullptr);
  auto* faceDataPtr = AttributeMatrix::Create(dataStructure, "FaceData", {1}, geometryPtr->getId());
  REQUIRE(faceDataPtr != nullptr);
  geometryPtr->setFaceAttributeMatrix(*faceDataPtr);
  auto& vertices = CreateSupportFixtureArray<float32>(dataStructure, k_MeshGeomPath.createChildPath("Vertices"), 3, 3, geometryPtr->getId());
  auto& faces = CreateSupportFixtureArray<SupportIndexType>(dataStructure, k_MeshGeomPath.createChildPath("Faces"), 1, 3, geometryPtr->getId());
  vertices[3] = 1.0F;
  vertices[7] = 1.0F;
  faces[1] = 1;
  faces[2] = 2;
  geometryPtr->setVertices(vertices);
  geometryPtr->setFaceList(faces);
  if(incidence)
  {
    CreateSupportData(dataStructure, SupportKind::VertexIncidence, {{0}, {0}, {0}});
  }
  if(neighbors)
  {
    CreateSupportData(dataStructure, SupportKind::FaceNeighbors, {{}});
  }
  Arguments args;
  args.insertOrAssign(LabelTriangleGeometryFilter::k_TriangleGeomPath_Key, std::make_any<DataPath>(k_MeshGeomPath));
  args.insertOrAssign(LabelTriangleGeometryFilter::k_CreatedRegionIdsPath_Key, std::make_any<DataPath>(k_RegionIdsPath));
  args.insertOrAssign(LabelTriangleGeometryFilter::k_TriangleAttributeMatrixName_Key, std::make_any<std::string>("Features"));
  args.insertOrAssign(LabelTriangleGeometryFilter::k_NumTrianglesName_Key, std::make_any<std::string>("NumTriangles"));
  return args;
}

Result<> RunLabel(DataStructure& dataStructure, bool cancelled)
{
  LabelTriangleGeometryInputValues inputs{.triangle_geom_path = k_MeshGeomPath, .region_ids_path = k_RegionIdsPath, .triangle_am_path = k_MeshFeaturePath, .num_triangles_path = k_TriangleCountPath};
  const std::atomic_bool cancel = cancelled;
  return LabelTriangleGeometry(dataStructure, {}, cancel, &inputs)();
}

void CheckLabelInputs(const DataStructure& dataStructure)
{
  const auto* verticesPtr = dataStructure.getDataAs<Float32ArrayType>(k_MeshGeomPath.createChildPath("Vertices"));
  const auto* facesPtr = dataStructure.getDataAs<DataArray<SupportIndexType>>(k_MeshGeomPath.createChildPath("Faces"));
  REQUIRE(verticesPtr != nullptr);
  REQUIRE(facesPtr != nullptr);
  REQUIRE(verticesPtr->getNumberOfTuples() == 3);
  REQUIRE(facesPtr->getNumberOfTuples() == 1);
  const std::array<float32, 9> expectedVertices = {0, 0, 0, 1, 0, 0, 0, 1, 0};
  for(usize i = 0; i < expectedVertices.size(); ++i)
  {
    CHECK((*verticesPtr)[i] == expectedVertices[i]);
  }
  for(usize i = 0; i < 3; ++i)
  {
    CHECK((*facesPtr)[i] == i);
  }
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<TriangleGeom>(k_MeshGeomPath));
  const auto& geometry = dataStructure.getDataRefAs<TriangleGeom>(k_MeshGeomPath);
  CHECK(geometry.getVertices() == verticesPtr);
  CHECK(geometry.getFaces() == facesPtr);
}

void CheckLabelOutputs(const DataStructure& dataStructure)
{
  const auto* regionsPtr = dataStructure.getDataAs<Int32ArrayType>(k_RegionIdsPath);
  const auto* countsPtr = dataStructure.getDataAs<UInt64ArrayType>(k_TriangleCountPath);
  REQUIRE(regionsPtr != nullptr);
  REQUIRE(countsPtr != nullptr);
  REQUIRE(regionsPtr->getNumberOfTuples() == 1);
  CHECK((*regionsPtr)[0] == 1);
  REQUIRE(countsPtr->getTupleShape() == ShapeType{3});
  const std::array<uint64, 3> expected = {0, 1, 0};
  for(usize i = 0; i < expected.size(); ++i)
  {
    CHECK((*countsPtr)[i] == expected[i]);
  }
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<AttributeMatrix>(k_MeshFeaturePath));
  CHECK(dataStructure.getDataRefAs<AttributeMatrix>(k_MeshFeaturePath).getShape() == ShapeType{3});
}

void CheckLabelSave(DataStructure& dataStructure, const SupportSnapshot& incidence, const SupportSnapshot& neighbors)
{
  SupportSaveFile file("label");
  // The writer omits supplied connectivity and leaves the caller's objects in place.
  const auto save = hdf5::DataStructureWriter::WriteFile(dataStructure, file.path);
  SIMPLNX_RESULT_REQUIRE_VALID(save);
  RequireSupportSnapshot(dataStructure, SupportKind::VertexIncidence, incidence);
  RequireSupportSnapshot(dataStructure, SupportKind::FaceNeighbors, neighbors);
  CheckLabelInputs(dataStructure);
  CheckLabelOutputs(dataStructure);
  const auto read = hdf5::DataStructureReader::ReadFile(file.path);
  SIMPLNX_RESULT_REQUIRE_VALID(read);
  CheckLabelInputs(read.value());
  CheckLabelOutputs(read.value());
  RequireSupportSnapshot(read.value(), SupportKind::VertexIncidence, {});
  RequireSupportSnapshot(read.value(), SupportKind::FaceNeighbors, {});
  unit_test::CheckArraysInheritTupleDims(read.value());
}

void CheckLabelSuccess(bool incidence, bool neighbors, bool cancelled, bool save, bool ooc)
{
  DataStructure dataStructure;
  const auto args = CreateLabelFixture(dataStructure, incidence, neighbors);
  const auto incidenceBefore = CaptureSupportSnapshot(dataStructure, SupportKind::VertexIncidence);
  const auto neighborsBefore = CaptureSupportSnapshot(dataStructure, SupportKind::FaceNeighbors);
  LabelTriangleGeometryFilter filter;
  if(cancelled)
  {
    ApplySupportFixtureActions(dataStructure, filter, args);
    const auto result = RunLabel(dataStructure, true);
    SIMPLNX_RESULT_REQUIRE_VALID(result);
  }
  else
  {
    const auto preflight = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflight.output_actions);
    const auto result = filter.execute(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(result.result);
  }
  CheckLabelInputs(dataStructure);
  CheckLabelOutputs(dataStructure);
  if(ooc)
  {
    RequireSupportFixtureOocStores(dataStructure, {k_MeshGeomPath.createChildPath("Vertices"), k_MeshGeomPath.createChildPath("Faces"), k_RegionIdsPath, k_TriangleCountPath});
  }
  if(save)
  {
    CheckLabelSave(dataStructure, incidenceBefore, neighborsBefore);
  }
  else
  {
    RequireSupportSnapshot(dataStructure, SupportKind::VertexIncidence, incidenceBefore);
    RequireSupportSnapshot(dataStructure, SupportKind::FaceNeighbors, neighborsBefore);
  }
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}
} // namespace

TEST_CASE("SimplnxCore::LabelTriangleGeometryFilter: support ownership", "[SimplnxCore][LabelTriangleGeometryFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  const bool incidence = GENERATE(false, true);
  const bool neighbors = GENERATE(false, true);
  const bool cancelled = GENERATE(false, true);
  CAPTURE(incidence, neighbors, cancelled);
  CheckLabelSuccess(incidence, neighbors, cancelled, false, false);
}

TEST_CASE("SimplnxCore::LabelTriangleGeometryFilter: cleanup on resize failure", "[SimplnxCore][LabelTriangleGeometryFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  const bool incidence = GENERATE(false, true);
  const bool neighbors = GENERATE(false, true);
  CAPTURE(incidence, neighbors);
  DataStructure dataStructure;
  const auto args = CreateLabelFixture(dataStructure, incidence, neighbors);
  const auto incidenceBefore = CaptureSupportSnapshot(dataStructure, SupportKind::VertexIncidence);
  const auto neighborsBefore = CaptureSupportSnapshot(dataStructure, SupportKind::FaceNeighbors);
  ApplySupportFixtureActions(dataStructure, LabelTriangleGeometryFilter{}, args);
  auto store = std::make_shared<InjectedFailureStore<uint64>>(ShapeType{1}, ShapeType{1}, 73);
  store->reject_resize = true;
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<UInt64ArrayType>(k_TriangleCountPath));
  const auto replacement = dataStructure.getDataRefAs<UInt64ArrayType>(k_TriangleCountPath).setDataStore(store);
  SIMPLNX_RESULT_REQUIRE_VALID(replacement);
  const auto result = RunLabel(dataStructure, false);
  RequireInjectedFailure(result, "AttributeMatrix 'Features' failed to resize child array 'Triangles/Features/NumTriangles' to shape [3]: injected resize failure");
  RequireSupportSnapshot(dataStructure, SupportKind::VertexIncidence, incidenceBefore);
  RequireSupportSnapshot(dataStructure, SupportKind::FaceNeighbors, neighborsBefore);
  CheckLabelInputs(dataStructure);
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Int32ArrayType>(k_RegionIdsPath));
  CHECK(dataStructure.getDataRefAs<Int32ArrayType>(k_RegionIdsPath)[0] == 1);
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<AttributeMatrix>(k_MeshFeaturePath));
  CHECK(dataStructure.getDataRefAs<AttributeMatrix>(k_MeshFeaturePath).getShape() == ShapeType{3});
  CHECK(store->getTupleShape() == ShapeType{1});
  CHECK((*store)[0] == 73);
  unit_test::CheckArraysInheritTupleDims(dataStructure, {k_TriangleCountPath});
}

TEST_CASE("SimplnxCore::LabelTriangleGeometryFilter: support creation failure", "[SimplnxCore][LabelTriangleGeometryFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  // 0: incidence count; 1: third read/orphan neighbor; 2: reused incidence/first neighbor read.
  const usize mode = GENERATE(0, 1, 2);
  CAPTURE(mode);
  DataStructure dataStructure;
  const auto args = CreateLabelFixture(dataStructure, mode == 2, false);
  const auto incidenceBefore = CaptureSupportSnapshot(dataStructure, SupportKind::VertexIncidence);
  ApplySupportFixtureActions(dataStructure, LabelTriangleGeometryFilter{}, args);
  auto store = std::make_shared<InjectedFailureStore<SupportIndexType>>(ShapeType{1}, ShapeType{3}, 0);
  (*store)[1] = 1;
  (*store)[2] = 2;
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<DataArray<SupportIndexType>>(k_MeshGeomPath.createChildPath("Faces")));
  const auto replacement = dataStructure.getDataRefAs<DataArray<SupportIndexType>>(k_MeshGeomPath.createChildPath("Faces")).setDataStore(store);
  SIMPLNX_RESULT_REQUIRE_VALID(replacement);
  store->reads.clear();
  store->failing_read = mode == 1 ? 3 : 1;
  const auto result = RunLabel(dataStructure, false);
  RequireInjectedFailure(result, "injected face read failure");
  const std::vector<std::pair<usize, usize>> expectedReads(store->failing_read, {0, 3});
  REQUIRE(store->reads == expectedReads);
  RequireSupportSnapshot(dataStructure, SupportKind::VertexIncidence, incidenceBefore);
  RequireSupportSnapshot(dataStructure, SupportKind::FaceNeighbors, {});
  CheckLabelInputs(dataStructure);
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::LabelTriangleGeometryFilter: unassociated named child", "[SimplnxCore][LabelTriangleGeometryFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  const bool neighbors = GENERATE(false, true);
  CAPTURE(neighbors);
  DataStructure dataStructure;
  const auto args = CreateLabelFixture(dataStructure, false, false);
  CreateSupportData(dataStructure, neighbors ? SupportKind::FaceNeighbors : SupportKind::VertexIncidence, {{42, 7}}, false);
  const auto incidenceBefore = CaptureSupportSnapshot(dataStructure, SupportKind::VertexIncidence);
  const auto neighborsBefore = CaptureSupportSnapshot(dataStructure, SupportKind::FaceNeighbors);
  ApplySupportFixtureActions(dataStructure, LabelTriangleGeometryFilter{}, args);
  const auto result = RunLabel(dataStructure, false);
  SIMPLNX_RESULT_REQUIRE_INVALID(result);
  REQUIRE(result.errors().size() == 1);
  CHECK(result.errors()[0].code == (neighbors ? -2232 : -2231));
  RequireSupportSnapshot(dataStructure, SupportKind::VertexIncidence, incidenceBefore);
  RequireSupportSnapshot(dataStructure, SupportKind::FaceNeighbors, neighborsBefore);
  CheckLabelInputs(dataStructure);
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::LabelTriangleGeometryFilter: retained support save compatibility", "[SimplnxCore][LabelTriangleGeometryFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  const bool incidence = GENERATE(false, true);
  const bool neighbors = GENERATE(false, true);
  CAPTURE(incidence, neighbors);
  CheckLabelSuccess(incidence, neighbors, false, true, false);
}

TEST_CASE("SimplnxCore::LabelTriangleGeometryFilter: real store support ownership", "[SimplnxCore][.OocStoreContract]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  REQUIRE(Application::GetOrCreateInstance()->getIOManager("HDF5-OOC") != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceOutOfCore, 1);
  const bool incidence = GENERATE(false, true);
  const bool neighbors = GENERATE(false, true);
  CAPTURE(incidence, neighbors);
  CheckLabelSuccess(incidence, neighbors, false, true, true);
  std::cout << fmt::format("Label ownership/save/recovery complete incidence={} neighbors={} vertices=HDF5-OOC faces=HDF5-OOC regions=HDF5-OOC counts=HDF5-OOC\n", incidence, neighbors);
}
