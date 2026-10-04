#include "SimplnxCore/Filters/ApplyTransformationToGeometryFilter.hpp"
#include "SimplnxCore/SimplnxCore_test_dirs.hpp"

#include "simplnx/Core/Application.hpp"
#include "simplnx/DataStructure/Geometry/RectGridGeom.hpp"
#include "simplnx/DataStructure/Geometry/TriangleGeom.hpp"
#include "simplnx/DataStructure/Geometry/VertexGeom.hpp"
#include "simplnx/DataStructure/IO/HDF5/DataStructureWriter.hpp"
#include "simplnx/DataStructure/NeighborList.hpp"
#include "simplnx/DataStructure/StringArray.hpp"
#include "simplnx/Parameters/ArrayCreationParameter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Parameters/DynamicTableParameter.hpp"
#include "simplnx/Parameters/VectorParameter.hpp"
#include "simplnx/Pipeline/Pipeline.hpp"
#include "simplnx/Pipeline/PipelineFilter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Utilities/DataStoreUtilities.hpp"
#include "simplnx/Utilities/ImageRotationUtilities.hpp"
#include "simplnx/Utilities/Parsing/HDF5/IO/FileIO.hpp"

#include <catch2/catch.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

using namespace nx::core;
using namespace nx::core::Constants;
using namespace nx::core::UnitTest;

namespace
{

} // namespace

namespace apply_transformation_to_geometry
{
const nx::core::ChoicesParameter::ValueType k_PrecomputedTransformationMatrixIdx = 1ULL;
const nx::core::ChoicesParameter::ValueType k_ManualTransformationMatrixIdx = 2ULL;
const nx::core::ChoicesParameter::ValueType k_RotationIdx = 3ULL;
const nx::core::ChoicesParameter::ValueType k_TranslationIdx = 4ULL;
const nx::core::ChoicesParameter::ValueType k_ScaleIdx = 5ULL;

const nx::core::ChoicesParameter::ValueType k_NearestNeighborInterpolationIdx = 0ULL;
const nx::core::ChoicesParameter::ValueType k_LinearInterpolationIdx = 1ULL;

const std::string k_InputGeometryName("InputData");
const std::string k_InputNodeGeometryName("InputNodeData");
const DataPath k_InputCellAttrMatrixPath(DataPath({k_InputGeometryName, "VertexData"}));
const std::string k_Rotation45XGeometryName("Rotation45X");
const std::string k_Rotation45YGeometryName("Rotation45Y");
const std::string k_Rotation45ZGeometryName("Rotation45Z");
const std::string k_Rotation90XGeometryName("Rotation90X");
const std::string k_Rotation90YGeometryName("Rotation90Y");
const std::string k_Rotation90ZGeometryName("Rotation90Z");
const std::string k_ScaleGeometryName("Scale");
const std::string k_ScaleNodeGeometryName("Scale_Node");
const std::string k_TranslationGeometryName("Translation");
const std::string k_TranslationNodeGeometryName("Translation_Node");
const std::string k_ManualGeometryName("Manual");
const std::string k_ManualNodeGeometryName("Manual_Node");
const std::string k_PrecomputedGeometryName("Precomputed");
const std::string k_PrecomputedNodeGeometryName("Precomputed_Node");
const std::string k_Rotation45XNodeGeometryName("Rotation_45X_Node");
const std::string k_Rotation45YNodeGeometryName("Rotation_45Y_Node");
const std::string k_Rotation45ZNodeGeometryName("Rotation_45Z_Node");
const std::string k_Rotation90XNodeGeometryName("Rotation_90X_Node");
const std::string k_Rotation90YNodeGeometryName("Rotation_90Y_Node");
const std::string k_Rotation90ZNodeGeometryName("Rotation_90Z_Node");
const std::string k_Rotation45XGlobalGeometryName("Rotation45X_Global");
const std::string k_Rotation45YGlobalGeometryName("Rotation45Y_Global");
const std::string k_Rotation45ZGlobalGeometryName("Rotation45Z_Global");
const std::string k_Rotation90XGlobalGeometryName("Rotation90X_Global");
const std::string k_Rotation90YGlobalGeometryName("Rotation90Y_Global");
const std::string k_Rotation90ZGlobalGeometryName("Rotation90Z_Global");
const std::string k_ScaleGlobalGeometryName("Scale_Global");
const std::string k_TranslationGlobalGeometryName("Translation_Global");
const std::string k_ManualGlobalGeometryName("Manual_Global");
const std::string k_PrecomputedGlobalGeometryName("Precomputed_Global");
const DataPath k_PrecomputedTransformationMatrixPath({"Transformation Matrices", "Precomputed"});
const std::string k_ExemplaryNNDataName("Data_NN");
const std::string k_ExemplaryLinearDataName("Data_L");

const std::string k_SharedVertexListName("SharedVertexList");

const int32 k_CellAttrMatrixUnusedWarning = -5555;

void CompareImageGeometries(const DataStructure& dataStructure, const ImageGeom& exemplaryGeom, const ImageGeom& calculatedGeom, const std::string& exemplaryDataName)
{
  UnitTest::CompareImageGeometry(&exemplaryGeom, &calculatedGeom);

  REQUIRE_NOTHROW(exemplaryGeom.getCellDataRef());
  REQUIRE_NOTHROW(calculatedGeom.getCellDataRef());
  auto exemplaryAM = exemplaryGeom.getCellDataRef();
  auto calculatedAM = calculatedGeom.getCellDataRef();
  REQUIRE(exemplaryAM.getShape() == calculatedAM.getShape());

  const DataPath exemplarPath({exemplaryGeom.getName(), k_Cell_Data, exemplaryDataName});
  const DataPath calculatedPath({calculatedGeom.getName(), k_Cell_Data, "Data"});
  const auto& exemplarData = dataStructure.getDataRefAs<IDataArray>(exemplarPath);
  const auto& calculatedData = dataStructure.getDataRefAs<IDataArray>(calculatedPath);
  UnitTest::CompareDataArrays<int32>(exemplarData, calculatedData);
}

// Build deterministic cell data in XYZ geometry order and ZYX tuple order.
ImageGeom* CreateSmallImage(DataStructure& dataStructure, const SizeVec3& dims = {3, 2, 1}, const FloatVec3& origin = {1.0F, 2.0F, 3.0F}, const FloatVec3& spacing = {1.0F, 1.0F, 1.0F},
                            int32 valueScale = 1)
{
  auto* imageGeom = ImageGeom::Create(dataStructure, "Image");
  REQUIRE(imageGeom != nullptr);
  imageGeom->setDimensions(dims);
  imageGeom->setOrigin(origin);
  imageGeom->setSpacing(spacing);
  const ShapeType tupleShape = {dims[2], dims[1], dims[0]};
  auto* cellAM = AttributeMatrix::Create(dataStructure, "Cell Data", tupleShape, imageGeom->getId());
  REQUIRE(cellAM != nullptr);
  imageGeom->setCellData(*cellAM);
  auto store = DataStoreUtilities::CreateDataStore<int32>(dataStructure, DataPath({"Image", "Cell Data", "Data"}), tupleShape, {1});
  REQUIRE(store != nullptr);
  auto* data = Int32Array::Create(dataStructure, "Data", store, cellAM->getId());
  REQUIRE(data != nullptr);
  for(usize tupleIdx = 0; tupleIdx < data->getNumberOfTuples(); tupleIdx++)
  {
    (*data)[tupleIdx] = valueScale * static_cast<int32>(tupleIdx + 1);
  }
  return imageGeom;
}

// M90 sends (x,y,z) to (-y,x,z); rotation is about the global coordinate origin.
Arguments SmallImageArguments(ChoicesParameter::ValueType transform = k_ManualTransformationMatrixIdx, ChoicesParameter::ValueType interpolation = k_NearestNeighborInterpolationIdx)
{
  Arguments args;
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(DataPath({"Image"})));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(DataPath({"Image", "Cell Data"})));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<ChoicesParameter::ValueType>(transform));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_InterpolationType_Key, std::make_any<ChoicesParameter::ValueType>(interpolation));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TranslateGeometryToGlobalOrigin_Key, std::make_any<bool>(false));
  const DynamicTableParameter::ValueType m90 = {{0, -1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_ManualTransformationMatrix_Key, std::make_any<DynamicTableParameter::ValueType>(m90));
  return args;
}

void CreateCellMask(DataStructure& dataStructure)
{
  auto* cellAM = dataStructure.getDataAs<AttributeMatrix>(DataPath({"Image", "Cell Data"}));
  REQUIRE(cellAM != nullptr);
  auto store = DataStoreUtilities::CreateDataStore<bool>(dataStructure, DataPath({"Image", "Cell Data", "Mask"}), cellAM->getShape(), {1});
  REQUIRE(store != nullptr);
  auto* mask = BoolArray::Create(dataStructure, "Mask", store, cellAM->getId());
  REQUIRE(mask != nullptr);
  const std::array<bool, 6> values = {true, false, true, true, false, true};
  REQUIRE(mask->getNumberOfTuples() == values.size());
  for(usize tupleIdx = 0; tupleIdx < values.size(); tupleIdx++)
  {
    (*mask)[tupleIdx] = values[tupleIdx];
  }
}

void CreateCellNames(DataStructure& dataStructure)
{
  auto* cellAM = dataStructure.getDataAs<AttributeMatrix>(DataPath({"Image", "Cell Data"}));
  REQUIRE(cellAM != nullptr);
  std::vector<std::string> names;
  for(usize tupleIdx = 0; tupleIdx < cellAM->getNumberOfTuples(); tupleIdx++)
  {
    names.emplace_back(1, static_cast<char>('a' + tupleIdx));
  }
  auto* namesArray = StringArray::CreateWithValues(dataStructure, "Names", cellAM->getShape(), names, cellAM->getId());
  REQUIRE(namesArray != nullptr);
}

void CreateCellNeighborLists(DataStructure& dataStructure)
{
  auto* cellAM = dataStructure.getDataAs<AttributeMatrix>(DataPath({"Image", "Cell Data"}));
  REQUIRE(cellAM != nullptr);
  auto* lists = NeighborList<int32>::Create(dataStructure, "NL", cellAM->getShape(), cellAM->getId());
  REQUIRE(lists != nullptr);
  for(usize tupleIdx = 0; tupleIdx < cellAM->getNumberOfTuples(); tupleIdx++)
  {
    const auto cellIdx = static_cast<int32>(tupleIdx);
    const std::vector<int32> values = cellIdx == 0 ? std::vector<int32>{} : std::vector<int32>{cellIdx, 10 * cellIdx};
    lists->setList(cellIdx, values);
  }
}

} // namespace apply_transformation_to_geometry

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Translation_Node", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const nx::core::UnitTest::TestFileSentinel testDataSentinel1(nx::core::unit_test::k_TestFilesDir, "apply_transformation_to_geometry_v2.tar.gz", "apply_transformation_to_geometry_v2.dream3d");

  auto baseDataFilePath = fs::path(fmt::format("{}/apply_transformation_to_geometry_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);
  const DataPath inputGeometryPath({apply_transformation_to_geometry::k_InputNodeGeometryName});
  {
    const ApplyTransformationToGeometryFilter filter;
    Arguments args;

    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(inputGeometryPath));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(apply_transformation_to_geometry::k_InputCellAttrMatrixPath));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_TranslationIdx));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Translation_Key, std::make_any<nx::core::VectorFloat32Parameter::ValueType>({100.0F, 50.0F, -100.0F}));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)
    REQUIRE(preflightResult.outputActions.warnings().size() == 1);
    REQUIRE(preflightResult.outputActions.warnings()[0].code == apply_transformation_to_geometry::k_CellAttrMatrixUnusedWarning);

    // Execute the filter and check the result
    auto executeResult = filter.execute(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
  }
#ifdef SIMPLNX_WRITE_TEST_OUTPUT
  WriteTestDataStructure(dataStructure, fmt::format("{}/apply_transformation_to_geometry_translation.dream3d", unit_test::k_BinaryTestOutputDir));
#endif
  {
    const DataPath exemplarPath({apply_transformation_to_geometry::k_TranslationNodeGeometryName, apply_transformation_to_geometry::k_SharedVertexListName});
    const DataPath calculatedPath({apply_transformation_to_geometry::k_InputNodeGeometryName, apply_transformation_to_geometry::k_SharedVertexListName});
    const auto& exemplarData = dataStructure.getDataRefAs<IDataArray>(exemplarPath);
    const auto& calculatedData = dataStructure.getDataRefAs<IDataArray>(calculatedPath);
    UnitTest::CompareDataArrays<float32>(exemplarData, calculatedData);
  }
  nx::core::HDF5::FileIO fileWriter = nx::core::HDF5::FileIO::WriteFile(fmt::format("{}/ApplyTransformationToGeometryFilter_translation.dream3d", unit_test::k_BinaryTestOutputDir));

  auto resultH5 = HDF5::DataStructureWriter::WriteFile(dataStructure, fileWriter);
  SIMPLNX_RESULT_REQUIRE_VALID(resultH5);

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Rotation_Node", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const nx::core::UnitTest::TestFileSentinel testDataSentinel1(nx::core::unit_test::k_TestFilesDir, "apply_transformation_to_geometry_v2.tar.gz", "apply_transformation_to_geometry_v2.dream3d");

  auto baseDataFilePath = fs::path(fmt::format("{}/apply_transformation_to_geometry_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);
  const DataPath inputGeometryPath({apply_transformation_to_geometry::k_InputNodeGeometryName});

  VectorFloat32Parameter::ValueType rotation45X = {1.0F, 0.0F, 0.0F, 45.0F};
  VectorFloat32Parameter::ValueType rotation45Y = {0.0F, 1.0F, 0.0F, 45.0F};
  VectorFloat32Parameter::ValueType rotation45Z = {0.0F, 0.0F, 1.0F, 45.0F};
  VectorFloat32Parameter::ValueType rotation90X = {1.0F, 0.0F, 0.0F, 90.0F};
  VectorFloat32Parameter::ValueType rotation90Y = {0.0F, 1.0F, 0.0F, 90.0F};
  VectorFloat32Parameter::ValueType rotation90Z = {0.0F, 0.0F, 1.0F, 90.0F};
  auto [exemplaryGeomName, rotation] = GENERATE_REF(
      std::make_tuple(apply_transformation_to_geometry::k_Rotation45XNodeGeometryName, rotation45X), std::make_tuple(apply_transformation_to_geometry::k_Rotation45YNodeGeometryName, rotation45Y),
      std::make_tuple(apply_transformation_to_geometry::k_Rotation45ZNodeGeometryName, rotation45Z), std::make_tuple(apply_transformation_to_geometry::k_Rotation90XNodeGeometryName, rotation90X),
      std::make_tuple(apply_transformation_to_geometry::k_Rotation90YNodeGeometryName, rotation90Y), std::make_tuple(apply_transformation_to_geometry::k_Rotation90ZNodeGeometryName, rotation90Z));

  {
    const ApplyTransformationToGeometryFilter filter;
    Arguments args;

    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(inputGeometryPath));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(apply_transformation_to_geometry::k_InputCellAttrMatrixPath));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_RotationIdx));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Rotation_Key, std::make_any<nx::core::VectorFloat32Parameter::ValueType>(rotation));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TranslateGeometryToGlobalOrigin_Key, std::make_any<nx::core::BoolParameter::ValueType>(true));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)
    REQUIRE(preflightResult.outputActions.warnings().size() == 1);
    REQUIRE(preflightResult.outputActions.warnings()[0].code == apply_transformation_to_geometry::k_CellAttrMatrixUnusedWarning);

    // Execute the filter and check the result
    auto executeResult = filter.execute(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
  }
#ifdef SIMPLNX_WRITE_TEST_OUTPUT
  WriteTestDataStructure(dataStructure, fmt::format("{}/apply_transformation_to_geometry_rotation.dream3d", unit_test::k_BinaryTestOutputDir));
#endif
  {
    const DataPath exemplarPath({exemplaryGeomName, apply_transformation_to_geometry::k_SharedVertexListName});
    const DataPath calculatedPath({apply_transformation_to_geometry::k_InputNodeGeometryName, apply_transformation_to_geometry::k_SharedVertexListName});
    const auto& exemplarData = dataStructure.getDataRefAs<IDataArray>(exemplarPath);
    const auto& calculatedData = dataStructure.getDataRefAs<IDataArray>(calculatedPath);
    UnitTest::CompareDataArrays<float32>(exemplarData, calculatedData);
  }
  nx::core::HDF5::FileIO fileWriter = nx::core::HDF5::FileIO::WriteFile(fmt::format("{}/ApplyTransformationToGeometryFilter_rotation.dream3d", unit_test::k_BinaryTestOutputDir));

  auto resultH5 = HDF5::DataStructureWriter::WriteFile(dataStructure, fileWriter);
  SIMPLNX_RESULT_REQUIRE_VALID(resultH5);

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Scale_Node", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const nx::core::UnitTest::TestFileSentinel testDataSentinel1(nx::core::unit_test::k_TestFilesDir, "apply_transformation_to_geometry_v2.tar.gz", "apply_transformation_to_geometry_v2.dream3d");

  auto baseDataFilePath = fs::path(fmt::format("{}/apply_transformation_to_geometry_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);
  const DataPath inputGeometryPath({apply_transformation_to_geometry::k_InputNodeGeometryName});
  {
    const ApplyTransformationToGeometryFilter filter;
    Arguments args;

    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(inputGeometryPath));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(apply_transformation_to_geometry::k_InputCellAttrMatrixPath));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_ScaleIdx));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Scale_Key, std::make_any<nx::core::VectorFloat32Parameter::ValueType>({0.5F, 1.5F, 10.0F}));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)
    REQUIRE(preflightResult.outputActions.warnings().size() == 1);
    REQUIRE(preflightResult.outputActions.warnings()[0].code == apply_transformation_to_geometry::k_CellAttrMatrixUnusedWarning);

    // Execute the filter and check the result
    auto executeResult = filter.execute(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
  }

#ifdef SIMPLNX_WRITE_TEST_OUTPUT
  WriteTestDataStructure(dataStructure, fmt::format("{}/apply_transformation_to_geometry_scale.dream3d", unit_test::k_BinaryTestOutputDir));
#endif

  {
    const DataPath exemplarPath({apply_transformation_to_geometry::k_ScaleNodeGeometryName, apply_transformation_to_geometry::k_SharedVertexListName});
    const DataPath calculatedPath({apply_transformation_to_geometry::k_InputNodeGeometryName, apply_transformation_to_geometry::k_SharedVertexListName});
    const auto& exemplarData = dataStructure.getDataRefAs<IDataArray>(exemplarPath);
    const auto& calculatedData = dataStructure.getDataRefAs<IDataArray>(calculatedPath);
    UnitTest::CompareDataArrays<float32>(exemplarData, calculatedData);
  }
  nx::core::HDF5::FileIO fileWriter = nx::core::HDF5::FileIO::WriteFile(fmt::format("{}/ApplyTransformationToGeometryFilter_scale.dream3d", unit_test::k_BinaryTestOutputDir));

  auto resultH5 = HDF5::DataStructureWriter::WriteFile(dataStructure, fileWriter);
  SIMPLNX_RESULT_REQUIRE_VALID(resultH5);

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Scale: Origin_And_Spacing Check", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const std::string k_SelectedImageGeomName = "Image Geometry";
  const DataPath k_SelectedImageGeomPath = DataPath({k_SelectedImageGeomName});
  const std::string k_SelectedAttrMatrixName = "Cell Data";
  const DataPath k_SelectedAttrMatrixPath = DataPath({k_SelectedImageGeomName, k_SelectedAttrMatrixName});

  DataStructure ds;
  ImageGeom* geom = ImageGeom::Create(ds, "Image Geometry");

  Vec3<float32> origin = {51.0F, 23.0F, -64.0F};
  Vec3<float32> spacing = {3.0F, 5.0F, 8.0F};
  std::vector<float32> scaleFactor;
  SECTION("Scale Increasing")
  {
    scaleFactor = {2.0F, 2.0F, 2.0F};
  }
  SECTION("Scale Decreasing")
  {
    scaleFactor = {0.4F, 0.4F, 0.4F};
  }
  geom->setOrigin(origin);
  geom->setSpacing(spacing);
  AttributeMatrix* cellAM = AttributeMatrix::Create(ds, "Cell Data", {}, geom->getId());

  const ApplyTransformationToGeometryFilter filter;
  Arguments args;

  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(k_SelectedImageGeomPath));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(k_SelectedAttrMatrixPath));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_ScaleIdx));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Scale_Key, std::make_any<nx::core::VectorFloat32Parameter::ValueType>(scaleFactor));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_InterpolationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(0));

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(ds, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)

  // Execute the filter and check the result
  auto executeResult = filter.execute(ds, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)

  Vec3<float32> newOrigin = {1.0f, 1.0f, 1.0f};
  std::transform(origin.begin(), origin.end(), scaleFactor.begin(), newOrigin.begin(), std::multiplies<>());
  Vec3<float32> newSpacing = {1.0f, 1.0f, 1.0f};
  std::transform(spacing.begin(), spacing.end(), scaleFactor.begin(), newSpacing.begin(), std::multiplies<>());
  REQUIRE(geom->getOrigin() == newOrigin);
  REQUIRE(geom->getSpacing() == newSpacing);

  UnitTest::CheckArraysInheritTupleDims(ds);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Manual_Node", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const nx::core::UnitTest::TestFileSentinel testDataSentinel1(nx::core::unit_test::k_TestFilesDir, "apply_transformation_to_geometry_v2.tar.gz", "apply_transformation_to_geometry_v2.dream3d");

  auto baseDataFilePath = fs::path(fmt::format("{}/apply_transformation_to_geometry_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);
  const DataPath inputGeometryPath({apply_transformation_to_geometry::k_InputNodeGeometryName});
  {
    const ApplyTransformationToGeometryFilter filter;
    Arguments args;

    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(inputGeometryPath));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(apply_transformation_to_geometry::k_InputCellAttrMatrixPath));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key,
                        std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_ManualTransformationMatrixIdx));
    // This should reflect the geometry across the x-axis.
    const DynamicTableParameter::ValueType dynamicTable{{{-1.0, 0, 0, 0}, {0, 1.0, 0, 0}, {0, 0, 1.0, 0}, {0, 0, 0, 1.0}}};
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_ManualTransformationMatrix_Key, std::make_any<nx::core::DynamicTableParameter::ValueType>(dynamicTable));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)
    REQUIRE(preflightResult.outputActions.warnings().size() == 1);
    REQUIRE(preflightResult.outputActions.warnings()[0].code == apply_transformation_to_geometry::k_CellAttrMatrixUnusedWarning);

    // Execute the filter and check the result
    auto executeResult = filter.execute(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
  }

#ifdef SIMPLNX_WRITE_TEST_OUTPUT
  WriteTestDataStructure(dataStructure, fmt::format("{}/apply_transformation_to_geometry_manual.dream3d", unit_test::k_BinaryTestOutputDir));
#endif

  {
    const DataPath exemplarPath({apply_transformation_to_geometry::k_ManualNodeGeometryName, apply_transformation_to_geometry::k_SharedVertexListName});
    const DataPath calculatedPath({apply_transformation_to_geometry::k_InputNodeGeometryName, apply_transformation_to_geometry::k_SharedVertexListName});
    const auto& exemplarData = dataStructure.getDataRefAs<IDataArray>(exemplarPath);
    const auto& calculatedData = dataStructure.getDataRefAs<IDataArray>(calculatedPath);
    UnitTest::CompareDataArrays<float32>(exemplarData, calculatedData);
  }
  nx::core::HDF5::FileIO fileWriter = nx::core::HDF5::FileIO::WriteFile(fmt::format("{}/ApplyTransformationToGeometryFilter_manual.dream3d", unit_test::k_BinaryTestOutputDir));

  auto resultH5 = HDF5::DataStructureWriter::WriteFile(dataStructure, fileWriter);
  SIMPLNX_RESULT_REQUIRE_VALID(resultH5);

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Precomputed_Node", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const nx::core::UnitTest::TestFileSentinel testDataSentinel1(nx::core::unit_test::k_TestFilesDir, "apply_transformation_to_geometry_v2.tar.gz", "apply_transformation_to_geometry_v2.dream3d");

  auto baseDataFilePath = fs::path(fmt::format("{}/apply_transformation_to_geometry_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);
  const DataPath inputGeometryPath({apply_transformation_to_geometry::k_InputNodeGeometryName});
  {
    const ApplyTransformationToGeometryFilter filter;
    Arguments args;

    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(inputGeometryPath));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(apply_transformation_to_geometry::k_InputCellAttrMatrixPath));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key,
                        std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_PrecomputedTransformationMatrixIdx));
    const DataPath precomputedPath({apply_transformation_to_geometry::k_InputNodeGeometryName, "Precomputed AM", "TransformationMatrix"});
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_ComputedTransformationMatrix_Key, std::make_any<DataPath>(precomputedPath));

    // Preflight the filter and check result
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);
    REQUIRE(preflightResult.outputActions.warnings().size() == 1);
    REQUIRE(preflightResult.outputActions.warnings()[0].code == apply_transformation_to_geometry::k_CellAttrMatrixUnusedWarning);

    // Execute the filter and check the result
    auto executeResult = filter.execute(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
  }

#ifdef SIMPLNX_WRITE_TEST_OUTPUT
  WriteTestDataStructure(dataStructure, fmt::format("{}/apply_transformation_to_geometry_manual.dream3d", unit_test::k_BinaryTestOutputDir));
#endif

  {
    const DataPath exemplarPath({apply_transformation_to_geometry::k_PrecomputedNodeGeometryName, apply_transformation_to_geometry::k_SharedVertexListName});
    const DataPath calculatedPath({apply_transformation_to_geometry::k_InputNodeGeometryName, apply_transformation_to_geometry::k_SharedVertexListName});
    const auto& exemplarData = dataStructure.getDataRefAs<IDataArray>(exemplarPath);
    const auto& calculatedData = dataStructure.getDataRefAs<IDataArray>(calculatedPath);
    UnitTest::CompareDataArrays<float32>(exemplarData, calculatedData);
  }
  nx::core::HDF5::FileIO fileWriter = nx::core::HDF5::FileIO::WriteFile(fmt::format("{}/ApplyTransformationToGeometryFilter_precomputed.dream3d", unit_test::k_BinaryTestOutputDir));

  auto resultH5 = HDF5::DataStructureWriter::WriteFile(dataStructure, fileWriter);
  SIMPLNX_RESULT_REQUIRE_VALID(resultH5);

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

/*******************************************************************************
 * @brief This section is for Image Geometry with Nearest Neighbor Interpolation
 ******************************************************************************/
TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Translation_Image", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const nx::core::UnitTest::TestFileSentinel testDataSentinel1(nx::core::unit_test::k_TestFilesDir, "apply_transformation_to_geometry_v2.tar.gz", "apply_transformation_to_geometry_v2.dream3d");

  auto baseDataFilePath = fs::path(fmt::format("{}/apply_transformation_to_geometry_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);

  auto [translateGeomToGlobalOrigin, exemplaryGeomName, interpolationIdx] =
      GENERATE_REF(std::make_tuple(false, apply_transformation_to_geometry::k_TranslationGeometryName, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_TranslationGlobalGeometryName, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_TranslationGeometryName, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_TranslationGlobalGeometryName, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx));

  std::string interpolationTypeStr;
  std::string exemplaryDataName;
  if(interpolationIdx == apply_transformation_to_geometry::k_LinearInterpolationIdx)
  {
    interpolationTypeStr = "Linear";
    exemplaryDataName = apply_transformation_to_geometry::k_ExemplaryLinearDataName;
  }
  else
  {
    interpolationTypeStr = "Nearest Neighbor";
    exemplaryDataName = apply_transformation_to_geometry::k_ExemplaryNNDataName;
  }

  DYNAMIC_SECTION(fmt::format("Geometry Name = {}, Interpolation Type = {}, Translate To Global Origin = {}", exemplaryGeomName, interpolationTypeStr, translateGeomToGlobalOrigin))
  {
    const DataPath inputGeometryPath({apply_transformation_to_geometry::k_InputGeometryName});
    const DataPath inputCellAMPath = inputGeometryPath.createChildPath(k_Cell_Data);

    {
      const ApplyTransformationToGeometryFilter filter;
      Arguments args;

      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(inputGeometryPath));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_TranslationIdx));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_InterpolationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(interpolationIdx));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(inputCellAMPath));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Translation_Key, std::make_any<nx::core::VectorFloat32Parameter::ValueType>({-10.0F, 10.0F, 20.0F}));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TranslateGeometryToGlobalOrigin_Key, std::make_any<nx::core::BoolParameter::ValueType>(translateGeomToGlobalOrigin));

      // Preflight the filter and check result
      auto preflightResult = filter.preflight(dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)
      // Execute the filter and check the result
      auto executeResult = filter.execute(dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
    }
#ifdef SIMPLNX_WRITE_TEST_OUTPUT
    WriteTestDataStructure(dataStructure, fmt::format("{}/apply_transformation_to_geometry_translation.dream3d", unit_test::k_BinaryTestOutputDir));
#endif

    REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(DataPath({exemplaryGeomName})));
    REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(inputGeometryPath));
    auto exemplaryGeom = dataStructure.getDataRefAs<ImageGeom>(DataPath({exemplaryGeomName}));
    auto calculatedGeom = dataStructure.getDataRefAs<ImageGeom>(inputGeometryPath);
    apply_transformation_to_geometry::CompareImageGeometries(dataStructure, exemplaryGeom, calculatedGeom, exemplaryDataName);

    UnitTest::CheckArraysInheritTupleDims(dataStructure);
  }
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Rotation_Image", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const nx::core::UnitTest::TestFileSentinel testDataSentinel1(nx::core::unit_test::k_TestFilesDir, "apply_transformation_to_geometry_v2.tar.gz", "apply_transformation_to_geometry_v2.dream3d");

  const auto scenario = GENERATE(from_range(UnitTest::SelectAlgorithmTestScenariosForInMemoryStores()));
  CAPTURE(scenario);

  auto baseDataFilePath = fs::path(fmt::format("{}/apply_transformation_to_geometry_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);
  const DataPath inputGeometryPath({apply_transformation_to_geometry::k_InputGeometryName});
  const DataPath inputCellAMPath = inputGeometryPath.createChildPath(k_Cell_Data);

  VectorFloat32Parameter::ValueType rotation45X = {1.0F, 0.0F, 0.0F, 45.0F};
  VectorFloat32Parameter::ValueType rotation45Y = {0.0F, 1.0F, 0.0F, 45.0F};
  VectorFloat32Parameter::ValueType rotation45Z = {0.0F, 0.0F, 1.0F, 45.0F};
  VectorFloat32Parameter::ValueType rotation90X = {1.0F, 0.0F, 0.0F, 90.0F};
  VectorFloat32Parameter::ValueType rotation90Y = {0.0F, 1.0F, 0.0F, 90.0F};
  VectorFloat32Parameter::ValueType rotation90Z = {0.0F, 0.0F, 1.0F, 90.0F};
  auto [translateGeomToGlobalOrigin, exemplaryGeomName, rotation, interpolationIdx] =
      GENERATE_REF(std::make_tuple(false, apply_transformation_to_geometry::k_Rotation45XGeometryName, rotation45X, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation45XGlobalGeometryName, rotation45X, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation45YGeometryName, rotation45Y, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation45YGlobalGeometryName, rotation45Y, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation45ZGeometryName, rotation45Z, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation45ZGlobalGeometryName, rotation45Z, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation90XGeometryName, rotation90X, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation90XGlobalGeometryName, rotation90X, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation90YGeometryName, rotation90Y, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation90YGlobalGeometryName, rotation90Y, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation90ZGeometryName, rotation90Z, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation90ZGlobalGeometryName, rotation90Z, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation45XGeometryName, rotation45X, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation45XGlobalGeometryName, rotation45X, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation45YGeometryName, rotation45Y, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation45YGlobalGeometryName, rotation45Y, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation45ZGeometryName, rotation45Z, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation45ZGlobalGeometryName, rotation45Z, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation90XGeometryName, rotation90X, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation90XGlobalGeometryName, rotation90X, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation90YGeometryName, rotation90Y, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation90YGlobalGeometryName, rotation90Y, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_Rotation90ZGeometryName, rotation90Z, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_Rotation90ZGlobalGeometryName, rotation90Z, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx));

  std::string interpolationTypeStr;
  std::string exemplaryDataName;
  if(interpolationIdx == apply_transformation_to_geometry::k_LinearInterpolationIdx)
  {
    interpolationTypeStr = "Linear";
    exemplaryDataName = apply_transformation_to_geometry::k_ExemplaryLinearDataName;
  }
  else
  {
    interpolationTypeStr = "Nearest Neighbor";
    exemplaryDataName = apply_transformation_to_geometry::k_ExemplaryNNDataName;
  }

  DYNAMIC_SECTION(fmt::format("Geometry Name = {}, Rotation = [{}, {}, {}, {}], Interpolation Type = {}, Translate To Global Origin = {}", exemplaryGeomName, rotation[0], rotation[1], rotation[2],
                              rotation[3], interpolationTypeStr, translateGeomToGlobalOrigin))
  {
    {
      const ApplyTransformationToGeometryFilter filter;
      Arguments args;

      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(inputGeometryPath));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_RotationIdx));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_InterpolationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(interpolationIdx));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(inputCellAMPath));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Rotation_Key, std::make_any<nx::core::VectorFloat32Parameter::ValueType>(rotation));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TranslateGeometryToGlobalOrigin_Key, std::make_any<nx::core::BoolParameter::ValueType>(translateGeomToGlobalOrigin));

      // Preflight the filter and check result
      auto preflightResult = filter.preflight(dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)
      UnitTest::RequireAutomaticCreateArrayActions(preflightResult.outputActions, dataStructure.getDataRefAs<AttributeMatrix>(inputCellAMPath).getSize());
      std::vector<DataPath> expectedArrayPaths;
      const DataPath outputCellPath({fmt::format(".{}", inputGeometryPath.getTargetName()), inputCellAMPath.getTargetName()});
      for([[maybe_unused]] const auto& [id, child] : dataStructure.getDataRefAs<AttributeMatrix>(inputCellAMPath))
      {
        expectedArrayPaths.push_back(outputCellPath.createChildPath(child->getName()));
      }
      UnitTest::RequireAutomaticCreateArrayActions(preflightResult.outputActions, expectedArrayPaths);

      // Execute the filter and check the result
      UnitTest::AlgorithmTestScope algorithmTestScope(scenario);
      auto executeResult = algorithmTestScope.executeFilter(filter, dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
    }
#ifdef SIMPLNX_WRITE_TEST_OUTPUT
    WriteTestDataStructure(dataStructure, fmt::format("{}/apply_transformation_to_geometry_rotation.dream3d", unit_test::k_BinaryTestOutputDir));
#endif

    REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(DataPath({exemplaryGeomName})));
    REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(inputGeometryPath));
    auto exemplaryGeom = dataStructure.getDataRefAs<ImageGeom>(DataPath({exemplaryGeomName}));
    auto calculatedGeom = dataStructure.getDataRefAs<ImageGeom>(inputGeometryPath);
    apply_transformation_to_geometry::CompareImageGeometries(dataStructure, exemplaryGeom, calculatedGeom, exemplaryDataName);

    UnitTest::CheckArraysInheritTupleDims(dataStructure);
  }
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Scale_Image", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const nx::core::UnitTest::TestFileSentinel testDataSentinel1(nx::core::unit_test::k_TestFilesDir, "apply_transformation_to_geometry_v2.tar.gz", "apply_transformation_to_geometry_v2.dream3d");

  auto baseDataFilePath = fs::path(fmt::format("{}/apply_transformation_to_geometry_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);
  const DataPath inputGeometryPath({apply_transformation_to_geometry::k_InputGeometryName});
  const DataPath inputCellAMPath = inputGeometryPath.createChildPath(k_Cell_Data);

  auto [translateGeomToGlobalOrigin, exemplaryGeomName, interpolationIdx] =
      GENERATE_REF(std::make_tuple(false, apply_transformation_to_geometry::k_ScaleGeometryName, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_ScaleGlobalGeometryName, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_ScaleGeometryName, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_ScaleGlobalGeometryName, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx));

  std::string interpolationTypeStr;
  std::string exemplaryDataName;
  if(interpolationIdx == apply_transformation_to_geometry::k_LinearInterpolationIdx)
  {
    interpolationTypeStr = "Linear";
    exemplaryDataName = apply_transformation_to_geometry::k_ExemplaryLinearDataName;
  }
  else
  {
    interpolationTypeStr = "Nearest Neighbor";
    exemplaryDataName = apply_transformation_to_geometry::k_ExemplaryNNDataName;
  }

  DYNAMIC_SECTION(fmt::format("Geometry Name = {}, Interpolation Type = {}, Translate To Global Origin = {}", exemplaryGeomName, interpolationTypeStr, translateGeomToGlobalOrigin))
  {
    {
      const ApplyTransformationToGeometryFilter filter;
      Arguments args;

      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(inputGeometryPath));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_ScaleIdx));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_InterpolationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(interpolationIdx));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(inputCellAMPath));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Scale_Key, std::make_any<nx::core::VectorFloat32Parameter::ValueType>({0.05F, 0.05F, 0.05F}));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TranslateGeometryToGlobalOrigin_Key, std::make_any<nx::core::BoolParameter::ValueType>(translateGeomToGlobalOrigin));

      // Preflight the filter and check result
      auto preflightResult = filter.preflight(dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)

      // Execute the filter and check the result
      auto executeResult = filter.execute(dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
    }

#ifdef SIMPLNX_WRITE_TEST_OUTPUT
    WriteTestDataStructure(dataStructure, fmt::format("{}/apply_transformation_to_geometry_scale.dream3d", unit_test::k_BinaryTestOutputDir));
#endif

    REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(DataPath({exemplaryGeomName})));
    REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(inputGeometryPath));
    auto exemplaryGeom = dataStructure.getDataRefAs<ImageGeom>(DataPath({exemplaryGeomName}));
    auto calculatedGeom = dataStructure.getDataRefAs<ImageGeom>(inputGeometryPath);
    apply_transformation_to_geometry::CompareImageGeometries(dataStructure, exemplaryGeom, calculatedGeom, exemplaryDataName);

    UnitTest::CheckArraysInheritTupleDims(dataStructure);
  }
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Manual_Image", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const nx::core::UnitTest::TestFileSentinel testDataSentinel1(nx::core::unit_test::k_TestFilesDir, "apply_transformation_to_geometry_v2.tar.gz", "apply_transformation_to_geometry_v2.dream3d");

  const auto scenario = GENERATE(from_range(UnitTest::SelectAlgorithmTestScenariosForInMemoryStores()));
  CAPTURE(scenario);

  auto baseDataFilePath = fs::path(fmt::format("{}/apply_transformation_to_geometry_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);
  const DataPath inputGeometryPath({apply_transformation_to_geometry::k_InputGeometryName});
  const DataPath inputCellAMPath = inputGeometryPath.createChildPath(k_Cell_Data);

  auto [translateGeomToGlobalOrigin, exemplaryGeomName, interpolationIdx] =
      GENERATE_REF(std::make_tuple(false, apply_transformation_to_geometry::k_ManualGeometryName, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_ManualGlobalGeometryName, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_ManualGeometryName, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_ManualGlobalGeometryName, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx));

  std::string interpolationTypeStr;
  std::string exemplaryDataName;
  if(interpolationIdx == apply_transformation_to_geometry::k_LinearInterpolationIdx)
  {
    interpolationTypeStr = "Linear";
    exemplaryDataName = apply_transformation_to_geometry::k_ExemplaryLinearDataName;
  }
  else
  {
    interpolationTypeStr = "Nearest Neighbor";
    exemplaryDataName = apply_transformation_to_geometry::k_ExemplaryNNDataName;
  }

  DYNAMIC_SECTION(fmt::format("Geometry Name = {}, Interpolation Type = {}, Translate To Global Origin = {}", exemplaryGeomName, interpolationTypeStr, translateGeomToGlobalOrigin))
  {
    {
      const ApplyTransformationToGeometryFilter filter;
      Arguments args;

      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(inputGeometryPath));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key,
                          std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_ManualTransformationMatrixIdx));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_InterpolationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(interpolationIdx));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(inputCellAMPath)); // This should reflect the geometry across the x-axis.
      const DynamicTableParameter::ValueType dynamicTable{{{-1.0, 0, 0, 0}, {0, 1.0, 0, 0}, {0, 0, 1.0, 0}, {0, 0, 0, 1.0}}};
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_ManualTransformationMatrix_Key, std::make_any<nx::core::DynamicTableParameter::ValueType>(dynamicTable));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TranslateGeometryToGlobalOrigin_Key, std::make_any<nx::core::BoolParameter::ValueType>(translateGeomToGlobalOrigin));

      // Preflight the filter and check result
      auto preflightResult = filter.preflight(dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)

      // Execute the filter and check the result
      UnitTest::AlgorithmTestScope algorithmTestScope(scenario);
      auto executeResult = algorithmTestScope.executeFilter(filter, dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
    }

#ifdef SIMPLNX_WRITE_TEST_OUTPUT
    WriteTestDataStructure(dataStructure, fmt::format("{}/apply_transformation_to_geometry_manual.dream3d", unit_test::k_BinaryTestOutputDir));
#endif

    REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(DataPath({exemplaryGeomName})));
    REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(inputGeometryPath));
    auto exemplaryGeom = dataStructure.getDataRefAs<ImageGeom>(DataPath({exemplaryGeomName}));
    auto calculatedGeom = dataStructure.getDataRefAs<ImageGeom>(inputGeometryPath);
    apply_transformation_to_geometry::CompareImageGeometries(dataStructure, exemplaryGeom, calculatedGeom, exemplaryDataName);

    UnitTest::CheckArraysInheritTupleDims(dataStructure);
  }
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Precomputed_Image", "[SimplnxCore][ApplyTransformationToGeometryFilter]")
{
  UnitTest::LoadPlugins();

  const nx::core::UnitTest::TestFileSentinel testDataSentinel1(nx::core::unit_test::k_TestFilesDir, "apply_transformation_to_geometry_v2.tar.gz", "apply_transformation_to_geometry_v2.dream3d");

  const auto scenario = GENERATE(from_range(UnitTest::SelectAlgorithmTestScenariosForInMemoryStores()));
  CAPTURE(scenario);

  auto baseDataFilePath = fs::path(fmt::format("{}/apply_transformation_to_geometry_v2.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = UnitTest::LoadDataStructure(baseDataFilePath);
  const DataPath inputGeometryPath({apply_transformation_to_geometry::k_InputGeometryName});
  const DataPath inputCellAMPath = inputGeometryPath.createChildPath(k_Cell_Data);

  auto [translateGeomToGlobalOrigin, exemplaryGeomName, interpolationIdx] =
      GENERATE_REF(std::make_tuple(false, apply_transformation_to_geometry::k_PrecomputedGeometryName, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_PrecomputedGlobalGeometryName, apply_transformation_to_geometry::k_LinearInterpolationIdx),
                   std::make_tuple(false, apply_transformation_to_geometry::k_PrecomputedGeometryName, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx),
                   std::make_tuple(true, apply_transformation_to_geometry::k_PrecomputedGlobalGeometryName, apply_transformation_to_geometry::k_NearestNeighborInterpolationIdx));

  std::string interpolationTypeStr;
  std::string exemplaryDataName;
  if(interpolationIdx == apply_transformation_to_geometry::k_LinearInterpolationIdx)
  {
    interpolationTypeStr = "Linear";
    exemplaryDataName = apply_transformation_to_geometry::k_ExemplaryLinearDataName;
  }
  else
  {
    interpolationTypeStr = "Nearest Neighbor";
    exemplaryDataName = apply_transformation_to_geometry::k_ExemplaryNNDataName;
  }

  DYNAMIC_SECTION(fmt::format("Geometry Name = {}, Precomputed Matrix Path = {}, Interpolation Type = {}, Translate To Global Origin = {}", exemplaryGeomName,
                              apply_transformation_to_geometry::k_PrecomputedTransformationMatrixPath.toString(), interpolationTypeStr, translateGeomToGlobalOrigin))
  {
    {
      const ApplyTransformationToGeometryFilter filter;
      Arguments args;

      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(inputGeometryPath));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key,
                          std::make_any<nx::core::ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_PrecomputedTransformationMatrixIdx));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_InterpolationType_Key, std::make_any<nx::core::ChoicesParameter::ValueType>(interpolationIdx));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key, std::make_any<DataPath>(inputCellAMPath));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_ComputedTransformationMatrix_Key, std::make_any<DataPath>(apply_transformation_to_geometry::k_PrecomputedTransformationMatrixPath));
      args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TranslateGeometryToGlobalOrigin_Key, std::make_any<nx::core::BoolParameter::ValueType>(translateGeomToGlobalOrigin));

      // Preflight the filter and check result
      auto preflightResult = filter.preflight(dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions)

      // Execute the filter and check the result
      UnitTest::AlgorithmTestScope algorithmTestScope(scenario);
      auto executeResult = algorithmTestScope.executeFilter(filter, dataStructure, args);
      SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result)
    }

#ifdef SIMPLNX_WRITE_TEST_OUTPUT
    WriteTestDataStructure(dataStructure, fmt::format("{}/apply_transformation_to_geometry_manual.dream3d", unit_test::k_BinaryTestOutputDir));
#endif

    REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(DataPath({exemplaryGeomName})));
    REQUIRE_NOTHROW(dataStructure.getDataRefAs<ImageGeom>(inputGeometryPath));
    auto exemplaryGeom = dataStructure.getDataRefAs<ImageGeom>(DataPath({exemplaryGeomName}));
    auto calculatedGeom = dataStructure.getDataRefAs<ImageGeom>(inputGeometryPath);
    apply_transformation_to_geometry::CompareImageGeometries(dataStructure, exemplaryGeom, calculatedGeom, exemplaryDataName);

    UnitTest::CheckArraysInheritTupleDims(dataStructure);
  }
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter: SIMPL Backwards Compatibility", "[SimplnxCore][ApplyTransformationToGeometryFilter][BackwardsCompatibility]")
{
  auto app = Application::GetOrCreateInstance();
  UnitTest::LoadPlugins();
  auto filterList = app->getFilterList();

  const fs::path conversionDir = fs::path(nx::core::unit_test::k_SourceDir.view()) / "test" / "simpl_conversion";

  const std::vector<std::pair<std::string, fs::path>> fixtures = {
      {"SIMPL 6.5 (UUID)", conversionDir / "6_5" / "ApplyTransformationToGeometryFilter.json"},
      {"SIMPL 6.4 (Filter_Name)", conversionDir / "6_4" / "ApplyTransformationToGeometryFilter.json"},
  };

  for(const auto& [label, fixturePath] : fixtures)
  {
    DYNAMIC_SECTION(label)
    {
      auto pipelineResult = Pipeline::FromSIMPLFile(fixturePath, filterList);
      REQUIRE(pipelineResult.valid());

      auto& pipeline = pipelineResult.value();
      REQUIRE(pipeline.size() == 1);

      auto* pipelineFilter = dynamic_cast<PipelineFilter*>(pipeline.at(0));
      REQUIRE(pipelineFilter != nullptr);

      const IFilter* filter = pipelineFilter->getFilter();
      REQUIRE(filter != nullptr);
      REQUIRE(filter->uuid() == FilterTraits<ApplyTransformationToGeometryFilter>::uuid);

      CHECK(pipelineFilter->getComments().empty());

      const Arguments args = pipelineFilter->getArguments();
      // Complex type (FloatVec3p1FilterParameterConverter) - verified by successful pipeline loading
      CHECK(args.value<ChoicesParameter::ValueType>(ApplyTransformationToGeometryFilter::k_TransformationType_Key) == 0);
      CHECK(args.value<ChoicesParameter::ValueType>(ApplyTransformationToGeometryFilter::k_InterpolationType_Key) == 0);
      // Complex type (DynamicTableFilterParameterConverter) - verified by successful pipeline loading
      // Complex type (FloatVec3FilterParameterConverter) - verified by successful pipeline loading
      // Complex type (FloatVec3FilterParameterConverter) - verified by successful pipeline loading
      CHECK(args.value<DataPath>(ApplyTransformationToGeometryFilter::k_ComputedTransformationMatrix_Key) == DataPath({"DataContainer", "CellData", "TestArray"}));
      CHECK(args.value<DataPath>(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key) == DataPath({"DataContainer"}));
      CHECK(args.value<DataPath>(ApplyTransformationToGeometryFilter::k_CellAttributeMatrixPath_Key) == DataPath({"DataContainer", "CellData"}));
    }
  }
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:NoTransform_Image", "[SimplnxCore][ApplyTransformationToGeometryFilter][AT-1]")
{
  UnitTest::LoadPlugins();
  DataStructure dataStructure;
  apply_transformation_to_geometry::CreateSmallImage(dataStructure, {3, 2, 1}, {1, 2, 3}, {0.5F, 1.5F, 2.5F}, 10);
  const ApplyTransformationToGeometryFilter filter;
  auto args = apply_transformation_to_geometry::SmallImageArguments(0, 0);
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);
  REQUIRE(preflightResult.outputActions.warnings().size() == 1);
  REQUIRE(preflightResult.outputActions.warnings()[0].code == 82001);
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const auto* imageGeom = dataStructure.getDataAs<ImageGeom>(DataPath({"Image"}));
  REQUIRE(imageGeom != nullptr);
  REQUIRE(imageGeom->getDimensions() == SizeVec3(3, 2, 1));
  REQUIRE(imageGeom->getOrigin() == FloatVec3(1, 2, 3));
  REQUIRE(imageGeom->getSpacing() == FloatVec3(0.5F, 1.5F, 2.5F));
  REQUIRE(dataStructure.getDataAs<ImageGeom>(DataPath({".Image"})) == nullptr);
  const auto* data = dataStructure.getDataAs<Int32Array>(DataPath({"Image", "Cell Data", "Data"}));
  REQUIRE(data != nullptr);
  const std::array<int32, 6> expected = {10, 20, 30, 40, 50, 60};
  REQUIRE(data->getNumberOfTuples() == expected.size());
  for(usize tupleIdx = 0; tupleIdx < expected.size(); tupleIdx++)
  {
    CAPTURE(tupleIdx);
    REQUIRE((*data)[tupleIdx] == expected[tupleIdx]);
  }
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:NoTransform_SaveMatrix_Node", "[SimplnxCore][ApplyTransformationToGeometryFilter][AT-2]")
{
  UnitTest::LoadPlugins();
  DataStructure dataStructure;
  auto* vertexGeom = VertexGeom::Create(dataStructure, "Vertices");
  REQUIRE(vertexGeom != nullptr);
  auto* vertices = Float32Array::CreateWithStore<Float32DataStore>(dataStructure, "SharedVertexList", {2}, {3}, vertexGeom->getId());
  REQUIRE(vertices != nullptr);
  vertexGeom->setVertices(*vertices);
  const std::array<float32, 6> expectedVertices = {1, 2, 3, -4, 5, -6};
  for(usize valueIdx = 0; valueIdx < expectedVertices.size(); valueIdx++)
  {
    (*vertices)[valueIdx] = expectedVertices[valueIdx];
  }
  const ApplyTransformationToGeometryFilter filter;
  Arguments args;
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(DataPath({"Vertices"})));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<ChoicesParameter::ValueType>(0));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SaveTransformMatrix_Key, std::make_any<bool>(true));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformMatrixOutputPath_Key, std::make_any<DataPath>(DataPath({"Transformation Matrix"})));
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const auto* matrix = dataStructure.getDataAs<Float32Array>(DataPath({"Transformation Matrix"}));
  REQUIRE(matrix != nullptr);
  const std::array<float32, 16> expectedMatrix = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  REQUIRE(matrix->getSize() == expectedMatrix.size());
  for(usize valueIdx = 0; valueIdx < expectedMatrix.size(); valueIdx++)
  {
    CAPTURE(valueIdx);
    CHECK((*matrix)[valueIdx] == expectedMatrix[valueIdx]);
  }
  const auto* outputVertices = dataStructure.getDataAs<Float32Array>(DataPath({"Vertices", "SharedVertexList"}));
  REQUIRE(outputVertices != nullptr);
  REQUIRE(outputVertices->getNumberOfTuples() == 2);
  for(usize valueIdx = 0; valueIdx < expectedVertices.size(); valueIdx++)
  {
    REQUIRE((*outputVertices)[valueIdx] == expectedVertices[valueIdx]);
  }
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Translation_Scale_Image_NoInterpolation", "[SimplnxCore][ApplyTransformationToGeometryFilter][AT-3]")
{
  UnitTest::LoadPlugins();
  DataStructure dataStructure;
  apply_transformation_to_geometry::CreateSmallImage(dataStructure, {3, 2, 1}, {1, 2, 3}, {0.5F, 1.5F, 2.5F}, 10);
  auto args = apply_transformation_to_geometry::SmallImageArguments(0, 2);
  FloatVec3 expectedOrigin;
  FloatVec3 expectedSpacing;
  SECTION("Translation")
  {
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_TranslationIdx));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Translation_Key, std::make_any<VectorFloat32Parameter::ValueType>(VectorFloat32Parameter::ValueType{7, -8, 9}));
    expectedOrigin = {8, -6, 12}; // (1,2,3) + (7,-8,9).
    expectedSpacing = {0.5F, 1.5F, 2.5F};
  }
  SECTION("Scale")
  {
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_ScaleIdx));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Scale_Key, std::make_any<VectorFloat32Parameter::ValueType>(VectorFloat32Parameter::ValueType{2, 3, 4}));
    expectedOrigin = {2, 6, 12}; // (1,2,3) * (2,3,4), component by component.
    expectedSpacing = {1, 4.5F, 10};
  }
  const ApplyTransformationToGeometryFilter filter;
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);
  const auto* imageGeom = dataStructure.getDataAs<ImageGeom>(DataPath({"Image"}));
  REQUIRE(imageGeom != nullptr);
  REQUIRE(imageGeom->getDimensions() == SizeVec3(3, 2, 1));
  REQUIRE(imageGeom->getOrigin() == expectedOrigin);
  REQUIRE(imageGeom->getSpacing() == expectedSpacing);
  const auto* data = dataStructure.getDataAs<Int32Array>(DataPath({"Image", "Cell Data", "Data"}));
  REQUIRE(data != nullptr);
  const std::array<int32, 6> expected = {10, 20, 30, 40, 50, 60};
  REQUIRE(data->getNumberOfTuples() == expected.size());
  for(usize tupleIdx = 0; tupleIdx < expected.size(); tupleIdx++)
  {
    REQUIRE((*data)[tupleIdx] == expected[tupleIdx]);
  }
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Linear_Integer_Rounding", "[SimplnxCore][ApplyTransformationToGeometryFilter][AT-4]")
{
  UnitTest::LoadPlugins();
  const auto scenario = GENERATE(from_range(UnitTest::SelectAlgorithmTestScenariosForInMemoryStores()));
  CAPTURE(scenario);
  UnitTest::AlgorithmTestScope scope(scenario);

  SECTION("Unit interpolation rounds only the final value")
  {
    const ImageRotationUtilities::RotateImageGeometryWithTrilinearInterpolation<int32> interpolator(nullptr, nullptr, ImageRotationUtilities::RotateArgs{},
                                                                                                    ImageRotationUtilities::Matrix4fR::Identity(), nullptr);
    const std::vector<ImageRotationUtilities::AccumulationValueType<int32>> firstCorners = {0, 1, 2, 1, 0, 1, 2, 1};
    const std::vector<ImageRotationUtilities::AccumulationValueType<int32>> secondCorners = {2, 3, 3, 2, 2, 3, 3, 2};
    // Both Z planes are equal. First: lerp(0.5, 1.5, 0.5) = 1.
    // Second: lerp(2.6, 2.4, 0.3) = 2.54, which rounds to 3.
    const int32 firstValue = interpolator.calculateInterpolatedValue(firstCorners, Eigen::Vector3f(0.5F, 0.5F, 0.25F), 1, 0);
    const int32 secondValue = interpolator.calculateInterpolatedValue(secondCorners, Eigen::Vector3f(0.6F, 0.3F, 0.7F), 1, 0);
    CHECK(firstValue == 1);
    CHECK(secondValue == 3);
  }

  SECTION("Rotated integer field agrees with the analytic linear field")
  {
    DataStructure dataStructure;
    apply_transformation_to_geometry::CreateSmallImage(dataStructure, {6, 5, 4}, {2, 3, 4});
    auto* input = dataStructure.getDataAs<Int32Array>(DataPath({"Image", "Cell Data", "Data"}));
    REQUIRE(input != nullptr);
    for(usize z = 0; z < 4; z++)
    {
      for(usize y = 0; y < 5; y++)
      {
        for(usize x = 0; x < 6; x++)
        {
          (*input)[x + 6 * (y + 5 * z)] = static_cast<int32>(1 + x + 10 * y + 100 * z);
        }
      }
    }
    auto args = apply_transformation_to_geometry::SmallImageArguments(apply_transformation_to_geometry::k_RotationIdx, apply_transformation_to_geometry::k_LinearInterpolationIdx);
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Rotation_Key, std::make_any<VectorFloat32Parameter::ValueType>(VectorFloat32Parameter::ValueType{0, 0, 1, 45}));
    const ApplyTransformationToGeometryFilter filter;
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);
    auto executeResult = scope.executeFilter(filter, dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

    const auto* outputGeom = dataStructure.getDataAs<ImageGeom>(DataPath({"Image"}));
    const auto* output = dataStructure.getDataAs<Int32Array>(DataPath({"Image", "Cell Data", "Data"}));
    REQUIRE(outputGeom != nullptr);
    REQUIRE(output != nullptr);
    const auto dims = outputGeom->getDimensions();
    const auto origin = outputGeom->getOrigin();
    const auto spacing = outputGeom->getSpacing();
    REQUIRE(output->getNumberOfTuples() == dims[0] * dims[1] * dims[2]);
    const float64 c = std::sqrt(0.5);
    usize checkedCells = 0;
    for(usize z = 0; z < dims[2]; z++)
    {
      for(usize y = 0; y < dims[1]; y++)
      {
        for(usize x = 0; x < dims[0]; x++)
        {
          // Invert Rz(45): src = ((X+Y)/sqrt(2), (Y-X)/sqrt(2), Z).
          // Subtract the source origin and half a cell to get continuous sample indices.
          const float64 destX = origin[0] + (static_cast<float64>(x) + 0.5) * spacing[0];
          const float64 destY = origin[1] + (static_cast<float64>(y) + 0.5) * spacing[1];
          const float64 destZ = origin[2] + (static_cast<float64>(z) + 0.5) * spacing[2];
          const std::array<float64, 3> f = {c * (destX + destY) - 2.5, c * (destY - destX) - 3.5, destZ - 4.5};
          if(f[0] < -1.0e-4 || f[0] > 5.0 + 1.0e-4 || f[1] < -1.0e-4 || f[1] > 4.0 + 1.0e-4 || f[2] < -1.0e-4 || f[2] > 3.0 + 1.0e-4)
          {
            continue;
          }
          const float64 expected = 1 + f[0] + 10 * f[1] + 100 * f[2];
          const usize tupleIdx = x + dims[0] * (y + dims[1] * z);
          CAPTURE(x, y, z, expected, (*output)[tupleIdx]);
          REQUIRE(std::abs(static_cast<float64>((*output)[tupleIdx]) - expected) <= 0.501);
          checkedCells++;
        }
      }
    }
    REQUIRE(checkedCells >= 50);
    UnitTest::CheckArraysInheritTupleDims(dataStructure);
  }
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Linear_Rejects_Bool_String_NeighborList", "[SimplnxCore][ApplyTransformationToGeometryFilter][AT-5-Linear]")
{
  UnitTest::LoadPlugins();
  DataStructure dataStructure;
  apply_transformation_to_geometry::CreateSmallImage(dataStructure);
  int32 expectedError = 0;
  SECTION("Bool")
  {
    apply_transformation_to_geometry::CreateCellMask(dataStructure);
    expectedError = -82023;
  }
  SECTION("String")
  {
    apply_transformation_to_geometry::CreateCellNames(dataStructure);
    expectedError = -82021;
  }
  SECTION("NeighborList")
  {
    apply_transformation_to_geometry::CreateCellNeighborLists(dataStructure);
    expectedError = -82022;
  }
  const ApplyTransformationToGeometryFilter filter;
  auto args = apply_transformation_to_geometry::SmallImageArguments(apply_transformation_to_geometry::k_ManualTransformationMatrixIdx, apply_transformation_to_geometry::k_LinearInterpolationIdx);
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.outputActions);
  REQUIRE(preflightResult.outputActions.errors().size() == 1);
  REQUIRE(preflightResult.outputActions.errors()[0].code == expectedError);
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:NearestNeighbor_Copies_Bool_String_NeighborList", "[SimplnxCore][ApplyTransformationToGeometryFilter][AT-5-NearestNeighbor]")
{
  UnitTest::LoadPlugins();
  const auto scenario = GENERATE(from_range(UnitTest::SelectAlgorithmTestScenariosForInMemoryStores()));
  CAPTURE(scenario);
  UnitTest::AlgorithmTestScope scope(scenario);
  DataStructure dataStructure;
  const ApplyTransformationToGeometryFilter filter;
  auto args = apply_transformation_to_geometry::SmallImageArguments();

  SECTION("M90 copies all cell array types")
  {
    apply_transformation_to_geometry::CreateSmallImage(dataStructure);
    apply_transformation_to_geometry::CreateCellMask(dataStructure);
    apply_transformation_to_geometry::CreateCellNames(dataStructure);
    apply_transformation_to_geometry::CreateCellNeighborLists(dataStructure);
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);
    auto executeResult = scope.executeFilter(filter, dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

    const auto* imageGeom = dataStructure.getDataAs<ImageGeom>(DataPath({"Image"}));
    REQUIRE(imageGeom != nullptr);
    REQUIRE(imageGeom->getDimensions() == SizeVec3(2, 3, 1));
    REQUIRE(imageGeom->getOrigin() == FloatVec3(-4, 1, 3));
    REQUIRE(imageGeom->getSpacing() == FloatVec3(1, 1, 1));
    const auto* data = dataStructure.getDataAs<Int32Array>(DataPath({"Image", "Cell Data", "Data"}));
    const auto* mask = dataStructure.getDataAs<BoolArray>(DataPath({"Image", "Cell Data", "Mask"}));
    const auto* names = dataStructure.getDataAs<StringArray>(DataPath({"Image", "Cell Data", "Names"}));
    const auto* lists = dataStructure.getDataAs<NeighborList<int32>>(DataPath({"Image", "Cell Data", "NL"}));
    REQUIRE(data != nullptr);
    REQUIRE(mask != nullptr);
    REQUIRE(names != nullptr);
    REQUIRE(lists != nullptr);
    REQUIRE(data->getNumberOfTuples() == 6);
    REQUIRE(mask->getNumberOfTuples() == 6);
    REQUIRE(names->getNumberOfTuples() == 6);
    REQUIRE(lists->getNumberOfTuples() == 6);
    // M90 reverses source Y along output X, and source X becomes output Y.
    // Output tuples therefore copy source indices 3,0,4,1,5,2.
    const std::array<int32, 6> expectedData = {4, 1, 5, 2, 6, 3};
    const std::array<bool, 6> expectedMask = {true, true, false, false, true, true};
    const std::array<std::string, 6> expectedNames = {"d", "a", "e", "b", "f", "c"};
    const std::array<std::vector<int32>, 6> expectedLists = {std::vector<int32>{3, 30}, {}, {4, 40}, {1, 10}, {5, 50}, {2, 20}};
    for(usize tupleIdx = 0; tupleIdx < 6; tupleIdx++)
    {
      CAPTURE(tupleIdx);
      REQUIRE((*data)[tupleIdx] == expectedData[tupleIdx]);
      REQUIRE((*mask)[tupleIdx] == expectedMask[tupleIdx]);
      REQUIRE((*names)[tupleIdx] == expectedNames[tupleIdx]);
      REQUIRE(lists->getList(static_cast<int32>(tupleIdx)) == expectedLists[tupleIdx]);
    }
  }

  SECTION("Outside cells have zero data and empty strings and lists")
  {
    const SizeVec3 sourceDims = {4, 3, 1};
    const FloatVec3 sourceOrigin = {1, 2, 3};
    const FloatVec3 sourceSpacing = {1, 1, 1};
    apply_transformation_to_geometry::CreateSmallImage(dataStructure, sourceDims, sourceOrigin, sourceSpacing);
    apply_transformation_to_geometry::CreateCellNames(dataStructure);
    apply_transformation_to_geometry::CreateCellNeighborLists(dataStructure);
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_RotationIdx));
    args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Rotation_Key, std::make_any<VectorFloat32Parameter::ValueType>(VectorFloat32Parameter::ValueType{0, 0, 1, 45}));
    auto preflightResult = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);
    auto executeResult = scope.executeFilter(filter, dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

    const auto* outputGeom = dataStructure.getDataAs<ImageGeom>(DataPath({"Image"}));
    const auto* data = dataStructure.getDataAs<Int32Array>(DataPath({"Image", "Cell Data", "Data"}));
    const auto* names = dataStructure.getDataAs<StringArray>(DataPath({"Image", "Cell Data", "Names"}));
    const auto* lists = dataStructure.getDataAs<NeighborList<int32>>(DataPath({"Image", "Cell Data", "NL"}));
    REQUIRE(outputGeom != nullptr);
    REQUIRE(data != nullptr);
    REQUIRE(names != nullptr);
    REQUIRE(lists != nullptr);
    const auto dims = outputGeom->getDimensions();
    const auto origin = outputGeom->getOrigin();
    const auto spacing = outputGeom->getSpacing();
    const usize numTuples = dims[0] * dims[1] * dims[2];
    REQUIRE(data->getNumberOfTuples() == numTuples);
    REQUIRE(names->getNumberOfTuples() == numTuples);
    REQUIRE(lists->getNumberOfTuples() == numTuples);
    const float64 c = std::sqrt(0.5);
    usize insideCells = 0;
    usize outsideCells = 0;
    for(usize z = 0; z < dims[2]; z++)
    {
      for(usize y = 0; y < dims[1]; y++)
      {
        for(usize x = 0; x < dims[0]; x++)
        {
          // A destination cell center is origin + (index + 1/2) * spacing.
          // Rz(-45) maps it to src = ((X+Y)/sqrt(2), (Y-X)/sqrt(2), Z).
          // Source centers lie at origin + (i+1/2)*spacing, so their nearest
          // index is floor((src-origin)/spacing). Check signed indices before
          // flattening; an index outside [0,dim) has no source cell.
          const float64 destX = origin[0] + (static_cast<float64>(x) + 0.5) * spacing[0];
          const float64 destY = origin[1] + (static_cast<float64>(y) + 0.5) * spacing[1];
          const float64 destZ = origin[2] + (static_cast<float64>(z) + 0.5) * spacing[2];
          const std::array<int64, 3> sourceCell = {static_cast<int64>(std::floor((c * (destX + destY) - sourceOrigin[0]) / sourceSpacing[0])),
                                                   static_cast<int64>(std::floor((c * (destY - destX) - sourceOrigin[1]) / sourceSpacing[1])),
                                                   static_cast<int64>(std::floor((destZ - sourceOrigin[2]) / sourceSpacing[2]))};
          const usize tupleIdx = x + dims[0] * (y + dims[1] * z);
          CAPTURE(x, y, z, sourceCell);
          const bool isOutside = sourceCell[0] < 0 || sourceCell[0] >= static_cast<int64>(sourceDims[0]) || sourceCell[1] < 0 || sourceCell[1] >= static_cast<int64>(sourceDims[1]) ||
                                 sourceCell[2] < 0 || sourceCell[2] >= static_cast<int64>(sourceDims[2]);
          if(isOutside)
          {
            REQUIRE((*data)[tupleIdx] == 0);
            REQUIRE((*names)[tupleIdx].empty());
            REQUIRE(lists->getList(static_cast<int32>(tupleIdx)).empty());
            outsideCells++;
          }
          else
          {
            // The fixture has Data[i]=i+1, Names[i]='a'+i, NL[i]={i,10*i}
            // except NL[0], which is empty. X is the fastest source axis.
            const auto sourceIdx = static_cast<int32>(sourceCell[0] + 4 * (sourceCell[1] + 3 * sourceCell[2]));
            const std::string expectedName(1, static_cast<char>('a' + sourceIdx));
            const std::vector<int32> expectedList = sourceIdx == 0 ? std::vector<int32>{} : std::vector<int32>{sourceIdx, 10 * sourceIdx};
            REQUIRE((*data)[tupleIdx] == sourceIdx + 1);
            REQUIRE((*names)[tupleIdx] == expectedName);
            REQUIRE(lists->getList(static_cast<int32>(tupleIdx)) == expectedList);
            insideCells++;
          }
        }
      }
    }
    REQUIRE(outsideCells >= 1);
    REQUIRE(insideCells >= 6);
  }
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:Image_ChildNamesPreserved", "[SimplnxCore][ApplyTransformationToGeometryFilter][AT-6]")
{
  UnitTest::LoadPlugins();
  const auto scenario = GENERATE(from_range(UnitTest::SelectAlgorithmTestScenariosForInMemoryStores()));
  CAPTURE(scenario);
  UnitTest::AlgorithmTestScope scope(scenario);
  DataStructure dataStructure;
  auto* imageGeom = apply_transformation_to_geometry::CreateSmallImage(dataStructure);
  auto* featureAM = AttributeMatrix::Create(dataStructure, "ImageFeatures", {2}, imageGeom->getId());
  REQUIRE(featureAM != nullptr);
  auto store = DataStoreUtilities::CreateDataStore<float32>(dataStructure, DataPath({"Image", "ImageFeatures", "Avg"}), {2}, {1});
  REQUIRE(store != nullptr);
  auto* averages = Float32Array::Create(dataStructure, "Avg", store, featureAM->getId());
  REQUIRE(averages != nullptr);
  (*averages)[0] = 1.5F;
  (*averages)[1] = 2.5F;
  const ApplyTransformationToGeometryFilter filter;
  auto args = apply_transformation_to_geometry::SmallImageArguments();
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);
  auto executeResult = scope.executeFilter(filter, dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  const auto* outputAverages = dataStructure.getDataAs<Float32Array>(DataPath({"Image", "ImageFeatures", "Avg"}));
  REQUIRE(outputAverages != nullptr);
  REQUIRE(outputAverages->getNumberOfTuples() == 2);
  REQUIRE((*outputAverages)[0] == 1.5F);
  REQUIRE((*outputAverages)[1] == 2.5F);
  REQUIRE(dataStructure.getDataAs<AttributeMatrix>(DataPath({"Image", ".ImageFeatures"})) == nullptr);
  const auto* data = dataStructure.getDataAs<Int32Array>(DataPath({"Image", "Cell Data", "Data"}));
  REQUIRE(data != nullptr);
  // As in the M90 fixture: source indices 3,0,4,1,5,2.
  const std::array<int32, 6> expected = {4, 1, 5, 2, 6, 3};
  REQUIRE(data->getNumberOfTuples() == expected.size());
  for(usize tupleIdx = 0; tupleIdx < expected.size(); tupleIdx++)
  {
    REQUIRE((*data)[tupleIdx] == expected[tupleIdx]);
  }
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ApplyTransformationToGeometryFilter:RectGrid_Rejected", "[SimplnxCore][ApplyTransformationToGeometryFilter][AT-7]")
{
  UnitTest::LoadPlugins();
  const bool globalOrigin = GENERATE(false, true);
  CAPTURE(globalOrigin);
  DataStructure dataStructure;
  auto* gridGeom = RectGridGeom::Create(dataStructure, "Grid");
  REQUIRE(gridGeom != nullptr);
  gridGeom->setDimensions({3, 2, 1});
  auto* xBounds = Float32Array::CreateWithStore<Float32DataStore>(dataStructure, "X Bounds", {4}, {1}, gridGeom->getId());
  auto* yBounds = Float32Array::CreateWithStore<Float32DataStore>(dataStructure, "Y Bounds", {3}, {1}, gridGeom->getId());
  auto* zBounds = Float32Array::CreateWithStore<Float32DataStore>(dataStructure, "Z Bounds", {2}, {1}, gridGeom->getId());
  REQUIRE(xBounds != nullptr);
  REQUIRE(yBounds != nullptr);
  REQUIRE(zBounds != nullptr);
  const std::array<float32, 4> xValues = {0, 1, 3, 6};
  const std::array<float32, 3> yValues = {0, 2, 5};
  const std::array<float32, 2> zValues = {0, 4};
  for(usize boundIdx = 0; boundIdx < xValues.size(); boundIdx++)
  {
    (*xBounds)[boundIdx] = xValues[boundIdx];
  }
  for(usize boundIdx = 0; boundIdx < yValues.size(); boundIdx++)
  {
    (*yBounds)[boundIdx] = yValues[boundIdx];
  }
  for(usize boundIdx = 0; boundIdx < zValues.size(); boundIdx++)
  {
    (*zBounds)[boundIdx] = zValues[boundIdx];
  }
  auto boundsResult = gridGeom->setBounds(xBounds, yBounds, zBounds);
  SIMPLNX_RESULT_REQUIRE_VALID(boundsResult);
  const ApplyTransformationToGeometryFilter filter;
  Arguments args;
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(DataPath({"Grid"})));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TransformationType_Key, std::make_any<ChoicesParameter::ValueType>(apply_transformation_to_geometry::k_TranslationIdx));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_Translation_Key, std::make_any<VectorFloat32Parameter::ValueType>(VectorFloat32Parameter::ValueType{1, 2, 3}));
  args.insertOrAssign(ApplyTransformationToGeometryFilter::k_TranslateGeometryToGlobalOrigin_Key, std::make_any<bool>(globalOrigin));
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.outputActions);
  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}
