#include "TriangleGeometrySupportTestUtils.hpp"

#include "SimplnxCore/Filters/Algorithms/FeatureFaceCurvature.hpp"
#include "SimplnxCore/Filters/FeatureFaceCurvatureFilter.hpp"
#include "SimplnxCore/SimplnxCore_test_dirs.hpp"

#include "simplnx/Core/Application.hpp"
#include "simplnx/DataStructure/DataStore.hpp"
#include "simplnx/DataStructure/Geometry/TriangleGeom.hpp"
#include "simplnx/DataStructure/IO/HDF5/DataStructureReader.hpp"
#include "simplnx/DataStructure/IO/HDF5/DataStructureWriter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <catch2/catch.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

using namespace nx::core;
using namespace nx::core::TriangleGeometrySupportTest;
namespace fs = std::filesystem;

namespace
{
enum class FloatComparison
{
  Standard,
  Absolute
};

template <FloatComparison ComparisonValue = FloatComparison::Standard>
void CompareDataArrays(const DataStructure& dataStructure, const DataPath& arrayPath1, const DataPath& arrayPath2)
{
  const auto* dataArray1Ptr = dataStructure.getDataAs<Float64ArrayType>(arrayPath1);
  const auto* dataArray2Ptr = dataStructure.getDataAs<Float64ArrayType>(arrayPath2);

  REQUIRE(dataArray1Ptr != nullptr);
  REQUIRE(dataArray2Ptr != nullptr);

  const auto length = dataArray1Ptr->getSize();
  REQUIRE(dataArray1Ptr->getSize() == dataArray2Ptr->getSize());

  const Float64AbstractDataStoreType& dataStore1 = dataArray1Ptr->getDataStoreRef();
  const Float64AbstractDataStoreType& dataStore2 = dataArray2Ptr->getDataStoreRef();

  for(usize i = 0; i < length; i++)
  {
    INFO(fmt::format("i = {} | d1 = {} | d2 = {}", i, dataStore1[i], dataStore2[i]));
    if constexpr(ComparisonValue == FloatComparison::Standard)
    {
      REQUIRE(unit_test::CloseEnough<float64>(dataStore1[i], dataStore2[i]));
    }
    else
    {
      REQUIRE(unit_test::CloseEnoughAbs<float64>(dataStore1[i], dataStore2[i]));
    }
  }
}
} // namespace

TEST_CASE("SimplnxCore::FeatureFaceCurvatureFilter: Test Algorithm", "[FeatureFaceCurvatureFilter]")
{
  unit_test::LoadPlugins();

  DataPath triangleGeomPath({"TriangleDataContainer"});
  DataPath faceAttribMatrixPath = triangleGeomPath.createChildPath("FaceData");

  DataPath kPrincipalCurvature1Path = faceAttribMatrixPath.createChildPath("PrincipalCurvature1_D3D-2");
  DataPath kPrincipalCurvature2Path = faceAttribMatrixPath.createChildPath("PrincipalCurvature2_D3D-2");
  DataPath kPrincipalDirection1Path = faceAttribMatrixPath.createChildPath("PrincipalDirection1_D3D-2");
  DataPath kPrincipalDirection2Path = faceAttribMatrixPath.createChildPath("PrincipalDirection2_D3D-2");
  DataPath kGaussianCurvaturePath = faceAttribMatrixPath.createChildPath("GaussianCurvatures-2");
  DataPath kMeanCurvaturePath = faceAttribMatrixPath.createChildPath("MeanCurvatures-2");
  DataPath kWeingartenMatrixPath = faceAttribMatrixPath.createChildPath("WeingartenMatrix-2");

  const nx::core::unit_test::TestFileSentinel testDataSentinel(nx::core::unit_test::k_TestFilesDir, "6_5_Goldfeather.tar.gz", "6_5_Goldfeather");

  FeatureFaceCurvatureFilter filter;

  // Read Exemplar DREAM3D File Filter
  auto exemplarFilePath = fs::path(fmt::format("{}/6_5_Goldfeather/6_5_Goldfeather.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = unit_test::LoadDataStructure(exemplarFilePath);

  Arguments args;
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_TriangleGeomPath_Key, std::make_any<DataPath>(triangleGeomPath));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_NeighborhoodRing_Key, std::make_any<int32>(3));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputePrincipalDirection_Key, std::make_any<bool>(true));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeGaussianCurvature_Key, std::make_any<bool>(true));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeMeanCurvaturePath_Key, std::make_any<bool>(true));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeWeingartenMatrix_Key, std::make_any<bool>(true));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_UseFaceNormals_Key, std::make_any<bool>(true));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceAttribMatrixPath_Key, std::make_any<DataPath>(faceAttribMatrixPath));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceLabelsPath_Key, std::make_any<DataPath>(faceAttribMatrixPath.createChildPath("FaceLabels")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FeatureFaceIdsPath_Key, std::make_any<DataPath>(faceAttribMatrixPath.createChildPath("FeatureFaceId")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceNormalsPath_Key, std::make_any<DataPath>(faceAttribMatrixPath.createChildPath("FaceNormals")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceCentroidsPath_Key, std::make_any<DataPath>(faceAttribMatrixPath.createChildPath("FaceCentroids")));

  args.insertOrAssign(FeatureFaceCurvatureFilter::k_PrincipalCurvature1Path_Key, std::make_any<DataPath>(kPrincipalCurvature1Path));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_PrincipalCurvature2Path_Key, std::make_any<DataPath>(kPrincipalCurvature2Path));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_PrincipalDirection1Path_Key, std::make_any<DataPath>(kPrincipalDirection1Path));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_PrincipalDirection2Path_Key, std::make_any<DataPath>(kPrincipalDirection2Path));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_GaussianCurvaturePath_Key, std::make_any<DataPath>(kGaussianCurvaturePath));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_MeanCurvaturePath_Key, std::make_any<DataPath>(kMeanCurvaturePath));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_WeingartenMatrixPath_Key, std::make_any<DataPath>(kWeingartenMatrixPath));

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);

  // Execute the filter and check the result
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  // Compare arrays
  {
    INFO("Principal Curvature 1");
    DataPath path1 = faceAttribMatrixPath.createChildPath("PrincipalCurvature1_D3D");
    ::CompareDataArrays(dataStructure, path1, kPrincipalCurvature1Path);
  }

  {
    INFO("Principal Curvature 2");
    DataPath path1 = faceAttribMatrixPath.createChildPath("PrincipalCurvature2_D3D");
    ::CompareDataArrays(dataStructure, path1, kPrincipalCurvature2Path);
  }

  /*
   * The principal directions can be sign flipped. This is due to Eigen::SelfAdjointEigenSolver
   * which returns unique eigenvectors up to a sign.
   */
  {
    INFO("Principal Direction 1");
    DataPath path1 = faceAttribMatrixPath.createChildPath("PrincipalDirection1_D3D");
    const DataPath& path2 = kPrincipalDirection1Path;
    const auto* dataArray1Ptr = dataStructure.getDataAs<Float64ArrayType>(path1);
    const auto* dataArray2Ptr = dataStructure.getDataAs<Float64ArrayType>(path2);
    CompareDataArrays<FloatComparison::Absolute>(dataStructure, path1, path2);
  }

  {
    INFO("Principal Direction 2");
    DataPath path1 = faceAttribMatrixPath.createChildPath("PrincipalDirection2_D3D");
    const DataPath& path2 = kPrincipalDirection2Path;
    const auto* dataArray1Ptr = dataStructure.getDataAs<Float64ArrayType>(path1);
    const auto* dataArray2Ptr = dataStructure.getDataAs<Float64ArrayType>(path2);
    CompareDataArrays<FloatComparison::Absolute>(dataStructure, path1, path2);
  }

  {
    INFO("Gaussian Curvatures");
    DataPath path1 = faceAttribMatrixPath.createChildPath("GaussianCurvatures");
    ::CompareDataArrays(dataStructure, path1, kGaussianCurvaturePath);
  }

  {
    INFO("Mean Curvatures");
    DataPath path1 = faceAttribMatrixPath.createChildPath("MeanCurvatures");
    ::CompareDataArrays(dataStructure, path1, kMeanCurvaturePath);
  }

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::FeatureFaceCurvatureFilter: NaN patch data", "[FeatureFaceCurvatureFilter][NaNPatch]")
{
  unit_test::LoadPlugins();
  const bool useNormals = GENERATE(false, true);
  const bool optionalOutputs = GENERATE(false, true);
  const usize badComponent = GENERATE(0, 1, 2, 3, 4, 5);
  // Nine faces leave eight valid samples for the seven-coefficient cubic fit.
  const usize badFace = GENERATE(8, 0, 4);
  const usize validCount = GENERATE(8, 1, 0);
  CAPTURE(useNormals, optionalOutputs, badComponent, badFace, validCount);

  DataStructure dataStructure;
  auto* geomPtr = TriangleGeom::Create(dataStructure, "Triangles");
  REQUIRE(geomPtr != nullptr);
  auto* faceDataPtr = AttributeMatrix::Create(dataStructure, "FaceData", {9}, geomPtr->getId());
  REQUIRE(faceDataPtr != nullptr);
  geomPtr->setFaceAttributeMatrix(*faceDataPtr);
  auto* verticesPtr = Float32ArrayType::Create(dataStructure, "Vertices", std::make_unique<DataStore<float32>>(ShapeType{19}, ShapeType{3}, 0.0F), geomPtr->getId());
  auto* facesPtr = DataArray<IGeometry::MeshIndexType>::Create(dataStructure, "Faces", std::make_unique<DataStore<IGeometry::MeshIndexType>>(ShapeType{9}, ShapeType{3}, 0), geomPtr->getId());
  REQUIRE(verticesPtr != nullptr);
  REQUIRE(facesPtr != nullptr);
  geomPtr->setVertices(*verticesPtr);
  geomPtr->setFaceList(*facesPtr);
  auto* labelsPtr = Int32ArrayType::Create(dataStructure, "Labels", std::make_unique<DataStore<int32>>(ShapeType{9}, ShapeType{2}, 1), faceDataPtr->getId());
  auto* idsPtr = Int32ArrayType::Create(dataStructure, "Ids", std::make_unique<DataStore<int32>>(ShapeType{9}, ShapeType{1}, 1), faceDataPtr->getId());
  auto* normalsPtr = Float64ArrayType::Create(dataStructure, "Normals", std::make_unique<DataStore<float64>>(ShapeType{9}, ShapeType{3}, 0.0), faceDataPtr->getId());
  auto* centroidsPtr = Float64ArrayType::Create(dataStructure, "Centroids", std::make_unique<DataStore<float64>>(ShapeType{9}, ShapeType{3}, 0.0), faceDataPtr->getId());
  REQUIRE(labelsPtr != nullptr);
  REQUIRE(idsPtr != nullptr);
  REQUIRE(normalsPtr != nullptr);
  REQUIRE(centroidsPtr != nullptr);

  for(usize t = 0; t < 9; ++t)
  {
    // All faces share vertex 0, so every one-ring patch is the sorted set {0,...,8}.
    // The other vertices are 1.5*c +/- (-cy,cx), giving centroid c and normal +Z.
    const float64 x = 1.0 + static_cast<float64>(t % 3);
    const float64 y = 1.0 + static_cast<float64>((t / 3) * (t / 3));
    verticesPtr->setComponent(2 * t + 1, 0, static_cast<float32>(1.5 * x + y));
    verticesPtr->setComponent(2 * t + 1, 1, static_cast<float32>(1.5 * y - x));
    verticesPtr->setComponent(2 * t + 2, 0, static_cast<float32>(1.5 * x - y));
    verticesPtr->setComponent(2 * t + 2, 1, static_cast<float32>(1.5 * y + x));
    facesPtr->setComponent(t, 1, 2 * t + 1);
    facesPtr->setComponent(t, 2, 2 * t + 2);
    labelsPtr->setComponent(t, 1, 2);
    normalsPtr->setComponent(t, 2, 1.0);
    centroidsPtr->setComponent(t, 0, x);
    centroidsPtr->setComponent(t, 1, y);
    const bool invalid = validCount == 0 || (validCount == 1 ? t != badFace : t == badFace);
    if(invalid)
    {
      auto* dataPtr = badComponent < 3 ? normalsPtr : centroidsPtr;
      dataPtr->setComponent(t, badComponent % 3, std::numeric_limits<float64>::quiet_NaN());
    }
  }

  const DataPath geomPath({"Triangles"});
  const DataPath facePath = geomPath.createChildPath("FaceData");
  Arguments args;
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_TriangleGeomPath_Key, std::make_any<DataPath>(geomPath));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceAttribMatrixPath_Key, std::make_any<DataPath>(facePath));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_NeighborhoodRing_Key, std::make_any<int32>(1));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_UseFaceNormals_Key, std::make_any<bool>(useNormals));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputePrincipalDirection_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeGaussianCurvature_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeMeanCurvaturePath_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeWeingartenMatrix_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceLabelsPath_Key, std::make_any<DataPath>(facePath.createChildPath("Labels")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FeatureFaceIdsPath_Key, std::make_any<DataPath>(facePath.createChildPath("Ids")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceNormalsPath_Key, std::make_any<DataPath>(facePath.createChildPath("Normals")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceCentroidsPath_Key, std::make_any<DataPath>(facePath.createChildPath("Centroids")));
  const std::array outputKeys = {FeatureFaceCurvatureFilter::k_PrincipalCurvature1Path_Key, FeatureFaceCurvatureFilter::k_PrincipalCurvature2Path_Key,
                                 FeatureFaceCurvatureFilter::k_GaussianCurvaturePath_Key,   FeatureFaceCurvatureFilter::k_MeanCurvaturePath_Key,
                                 FeatureFaceCurvatureFilter::k_WeingartenMatrixPath_Key,    FeatureFaceCurvatureFilter::k_PrincipalDirection1Path_Key,
                                 FeatureFaceCurvatureFilter::k_PrincipalDirection2Path_Key};
  for(usize i = 0; i < outputKeys.size(); ++i)
  {
    args.insertOrAssign(outputKeys[i], std::make_any<DataPath>(facePath.createChildPath("Output" + std::to_string(i))));
  }
  FeatureFaceCurvatureFilter filter;
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  for(usize i = 0; i < (optionalOutputs ? outputKeys.size() : 2); ++i)
  {
    const auto* outputPtr = dataStructure.getDataAs<Float64ArrayType>(facePath.createChildPath("Output" + std::to_string(i)));
    REQUIRE(outputPtr != nullptr);
    for(usize t = 0; t < 9; ++t)
    {
      for(usize c = 0; c < outputPtr->getNumberOfComponents(); ++c)
      {
        CAPTURE(i, t, c);
        const float64 value = outputPtr->getDataStoreRef().getComponentValue(t, c);
        if(validCount < 2 || t == badFace)
        {
          CHECK(std::isnan(value));
        }
        else if(i < 5)
        {
          // Closed-form oracle: z=0 has zero Hessian, principal/mean/Gaussian
          // curvatures and Weingarten matrix, including after removing a sample.
          CHECK(value == Approx(0.0).margin(1.0e-12));
        }
        else
        {
          CHECK(std::isfinite(value));
        }
      }
    }
  }
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::FeatureFaceCurvatureFilter: Single triangle face", "[SimplnxCore][FeatureFaceCurvatureFilter][SingleTriangle]")
{
  const bool useNormals = GENERATE(false, true);
  const bool optionalOutputs = GENERATE(false, true);
  CAPTURE(useNormals, optionalOutputs);

  DataStructure dataStructure;
  auto* geomPtr = TriangleGeom::Create(dataStructure, "Triangles");
  REQUIRE(geomPtr != nullptr);
  auto* faceDataPtr = AttributeMatrix::Create(dataStructure, "FaceData", {1}, geomPtr->getId());
  REQUIRE(faceDataPtr != nullptr);
  geomPtr->setFaceAttributeMatrix(*faceDataPtr);
  auto* verticesPtr = Float32ArrayType::Create(dataStructure, "Vertices", std::make_unique<DataStore<float32>>(ShapeType{3}, ShapeType{3}, 0.0F), geomPtr->getId());
  auto* facesPtr = DataArray<IGeometry::MeshIndexType>::Create(dataStructure, "Faces", std::make_unique<DataStore<IGeometry::MeshIndexType>>(ShapeType{1}, ShapeType{3}, 0), geomPtr->getId());
  REQUIRE(verticesPtr != nullptr);
  REQUIRE(facesPtr != nullptr);
  geomPtr->setVertices(*verticesPtr);
  geomPtr->setFaceList(*facesPtr);
  auto* labelsPtr = Int32ArrayType::Create(dataStructure, "Labels", std::make_unique<DataStore<int32>>(ShapeType{1}, ShapeType{2}, 1), faceDataPtr->getId());
  auto* idsPtr = Int32ArrayType::Create(dataStructure, "Ids", std::make_unique<DataStore<int32>>(ShapeType{1}, ShapeType{1}, 1), faceDataPtr->getId());
  auto* normalsPtr = Float64ArrayType::Create(dataStructure, "Normals", std::make_unique<DataStore<float64>>(ShapeType{1}, ShapeType{3}, 0.0), faceDataPtr->getId());
  auto* centroidsPtr = Float64ArrayType::Create(dataStructure, "Centroids", std::make_unique<DataStore<float64>>(ShapeType{1}, ShapeType{3}, 0.0), faceDataPtr->getId());
  REQUIRE(labelsPtr != nullptr);
  REQUIRE(idsPtr != nullptr);
  REQUIRE(normalsPtr != nullptr);
  REQUIRE(centroidsPtr != nullptr);

  // Triangle (0,0,0), (1,0,0), (0,1,0): centroid (1/3,1/3,0), normal +Z.
  verticesPtr->setComponent(1, 0, 1.0F);
  verticesPtr->setComponent(2, 1, 1.0F);
  facesPtr->setComponent(0, 1, 1);
  facesPtr->setComponent(0, 2, 2);
  labelsPtr->setComponent(0, 1, 2);
  normalsPtr->setComponent(0, 2, 1.0);
  centroidsPtr->setComponent(0, 0, 1.0 / 3.0);
  centroidsPtr->setComponent(0, 1, 1.0 / 3.0);

  const DataPath geomPath({"Triangles"});
  const DataPath facePath = geomPath.createChildPath("FaceData");
  Arguments args;
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_TriangleGeomPath_Key, std::make_any<DataPath>(geomPath));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceAttribMatrixPath_Key, std::make_any<DataPath>(facePath));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_NeighborhoodRing_Key, std::make_any<int32>(1));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_UseFaceNormals_Key, std::make_any<bool>(useNormals));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputePrincipalDirection_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeGaussianCurvature_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeMeanCurvaturePath_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeWeingartenMatrix_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceLabelsPath_Key, std::make_any<DataPath>(facePath.createChildPath("Labels")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FeatureFaceIdsPath_Key, std::make_any<DataPath>(facePath.createChildPath("Ids")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceNormalsPath_Key, std::make_any<DataPath>(facePath.createChildPath("Normals")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceCentroidsPath_Key, std::make_any<DataPath>(facePath.createChildPath("Centroids")));
  const std::array outputKeys = {FeatureFaceCurvatureFilter::k_PrincipalCurvature1Path_Key, FeatureFaceCurvatureFilter::k_PrincipalCurvature2Path_Key,
                                 FeatureFaceCurvatureFilter::k_GaussianCurvaturePath_Key,   FeatureFaceCurvatureFilter::k_MeanCurvaturePath_Key,
                                 FeatureFaceCurvatureFilter::k_WeingartenMatrixPath_Key,    FeatureFaceCurvatureFilter::k_PrincipalDirection1Path_Key,
                                 FeatureFaceCurvatureFilter::k_PrincipalDirection2Path_Key};
  for(usize i = 0; i < outputKeys.size(); ++i)
  {
    args.insertOrAssign(outputKeys[i], std::make_any<DataPath>(facePath.createChildPath("Output" + std::to_string(i))));
  }
  FeatureFaceCurvatureFilter filter;
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.output_actions);
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  // Independent oracle: one centroid and normal cannot determine surface curvature.
  // The insufficient-patch convention requires NaN in every requested component.
  for(usize i = 0; i < (optionalOutputs ? outputKeys.size() : 2); ++i)
  {
    const auto* outputPtr = dataStructure.getDataAs<Float64ArrayType>(facePath.createChildPath("Output" + std::to_string(i)));
    REQUIRE(outputPtr != nullptr);
    REQUIRE(outputPtr->getNumberOfTuples() == 1);
    for(usize c = 0; c < outputPtr->getNumberOfComponents(); ++c)
    {
      CAPTURE(i, c);
      CHECK(std::isnan(outputPtr->getDataStoreRef().getComponentValue(0, c)));
    }
  }
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

namespace
{
const std::array k_CurvatureOutputKeys = {FeatureFaceCurvatureFilter::k_PrincipalCurvature1Path_Key, FeatureFaceCurvatureFilter::k_PrincipalCurvature2Path_Key,
                                          FeatureFaceCurvatureFilter::k_GaussianCurvaturePath_Key,   FeatureFaceCurvatureFilter::k_MeanCurvaturePath_Key,
                                          FeatureFaceCurvatureFilter::k_WeingartenMatrixPath_Key,    FeatureFaceCurvatureFilter::k_PrincipalDirection1Path_Key,
                                          FeatureFaceCurvatureFilter::k_PrincipalDirection2Path_Key};

DataPath CurvatureOutputPath(usize index)
{
  return k_MeshFacePath.createChildPath("Output" + std::to_string(index));
}

Arguments CreateCurvatureFixture(DataStructure& dataStructure, bool incidence, bool optionalOutputs, bool neighbors = false)
{
  auto* geometryPtr = TriangleGeom::Create(dataStructure, "Triangles");
  REQUIRE(geometryPtr != nullptr);
  auto* faceDataPtr = AttributeMatrix::Create(dataStructure, "FaceData", {9}, geometryPtr->getId());
  REQUIRE(faceDataPtr != nullptr);
  geometryPtr->setFaceAttributeMatrix(*faceDataPtr);
  auto& vertices = CreateSupportFixtureArray<float32>(dataStructure, k_MeshGeomPath.createChildPath("Vertices"), 19, 3, geometryPtr->getId());
  auto& faces = CreateSupportFixtureArray<SupportIndexType>(dataStructure, k_MeshGeomPath.createChildPath("Faces"), 9, 3, geometryPtr->getId());
  geometryPtr->setVertices(vertices);
  geometryPtr->setFaceList(faces);
  auto& labels = CreateSupportFixtureArray<int32>(dataStructure, k_MeshFacePath.createChildPath("Labels"), 9, 2, faceDataPtr->getId());
  auto& ids = CreateSupportFixtureArray<int32>(dataStructure, k_MeshFacePath.createChildPath("Ids"), 9, 1, faceDataPtr->getId());
  auto& normals = CreateSupportFixtureArray<float64>(dataStructure, k_MeshFacePath.createChildPath("Normals"), 9, 3, faceDataPtr->getId());
  auto& centroids = CreateSupportFixtureArray<float64>(dataStructure, k_MeshFacePath.createChildPath("Centroids"), 9, 3, faceDataPtr->getId());
  for(usize t = 0; t < 9; ++t)
  {
    const float64 x = 1.0 + static_cast<float64>(t % 3);
    const float64 y = 1.0 + static_cast<float64>((t / 3) * (t / 3));
    vertices.setComponent(2 * t + 1, 0, static_cast<float32>(1.5 * x + y));
    vertices.setComponent(2 * t + 1, 1, static_cast<float32>(1.5 * y - x));
    vertices.setComponent(2 * t + 2, 0, static_cast<float32>(1.5 * x - y));
    vertices.setComponent(2 * t + 2, 1, static_cast<float32>(1.5 * y + x));
    faces.setComponent(t, 1, 2 * t + 1);
    faces.setComponent(t, 2, 2 * t + 2);
    labels.setComponent(t, 0, 11);
    labels.setComponent(t, 1, 13);
    ids[t] = 1;
    normals.setComponent(t, 2, 1.0);
    centroids.setComponent(t, 0, x);
    centroids.setComponent(t, 1, y);
  }
  if(incidence)
  {
    // Literal independent vertex incidence, in vertex order 0..18.
    CreateSupportData(dataStructure, SupportKind::VertexIncidence, {{0, 1, 2, 3, 4, 5, 6, 7, 8}, {0}, {0}, {1}, {1}, {2}, {2}, {3}, {3}, {4}, {4}, {5}, {5}, {6}, {6}, {7}, {7}, {8}, {8}});
  }
  if(neighbors)
  {
    CreateSupportData(dataStructure, SupportKind::FaceNeighbors, {{}, {}, {}, {}, {}, {}, {}, {}, {}});
  }
  Arguments args;
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_TriangleGeomPath_Key, std::make_any<DataPath>(k_MeshGeomPath));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceAttribMatrixPath_Key, std::make_any<DataPath>(k_MeshFacePath));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_NeighborhoodRing_Key, std::make_any<int32>(1));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_UseFaceNormals_Key, std::make_any<bool>(true));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputePrincipalDirection_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeGaussianCurvature_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeMeanCurvaturePath_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_ComputeWeingartenMatrix_Key, std::make_any<bool>(optionalOutputs));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceLabelsPath_Key, std::make_any<DataPath>(k_MeshFacePath.createChildPath("Labels")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FeatureFaceIdsPath_Key, std::make_any<DataPath>(k_MeshFacePath.createChildPath("Ids")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceNormalsPath_Key, std::make_any<DataPath>(k_MeshFacePath.createChildPath("Normals")));
  args.insertOrAssign(FeatureFaceCurvatureFilter::k_FaceCentroidsPath_Key, std::make_any<DataPath>(k_MeshFacePath.createChildPath("Centroids")));
  for(usize i = 0; i < k_CurvatureOutputKeys.size(); ++i)
  {
    args.insertOrAssign(k_CurvatureOutputKeys[i], std::make_any<DataPath>(CurvatureOutputPath(i)));
  }
  return args;
}

Result<> RunCurvature(DataStructure& dataStructure, bool cancelled)
{
  FeatureFaceCurvatureInputValues inputs{};
  inputs.n_ring_count = 1;
  inputs.triangle_geom_path = k_MeshGeomPath;
  inputs.surface_mesh_principal_curvature1_path = CurvatureOutputPath(0);
  inputs.surface_mesh_feature_face_ids_path = k_MeshFacePath.createChildPath("Ids");
  inputs.surface_mesh_face_labels_path = k_MeshFacePath.createChildPath("Labels");
  inputs.surface_mesh_face_normals_path = k_MeshFacePath.createChildPath("Normals");
  inputs.surface_mesh_triangle_centroids_path = k_MeshFacePath.createChildPath("Centroids");
  inputs.use_normals_for_curve_fitting = true;
  inputs.surface_mesh_principal_curvature1s_path = CurvatureOutputPath(0);
  inputs.surface_mesh_principal_curvature2s_path = CurvatureOutputPath(1);
  inputs.surface_mesh_gaussian_curvatures_path = CurvatureOutputPath(2);
  inputs.surface_mesh_mean_curvatures_path = CurvatureOutputPath(3);
  inputs.surface_mesh_weingarten_matrix_path = CurvatureOutputPath(4);
  inputs.surface_mesh_principal_direction1s_path = CurvatureOutputPath(5);
  inputs.surface_mesh_principal_direction2s_path = CurvatureOutputPath(6);
  const std::atomic_bool cancel = cancelled;
  return FeatureFaceCurvature(dataStructure, {}, cancel, &inputs)();
}

void CheckCurvatureInputs(const DataStructure& dataStructure)
{
  const auto* verticesPtr = dataStructure.getDataAs<Float32ArrayType>(k_MeshGeomPath.createChildPath("Vertices"));
  const auto* facesPtr = dataStructure.getDataAs<DataArray<SupportIndexType>>(k_MeshGeomPath.createChildPath("Faces"));
  const auto* labelsPtr = dataStructure.getDataAs<Int32ArrayType>(k_MeshFacePath.createChildPath("Labels"));
  const auto* idsPtr = dataStructure.getDataAs<Int32ArrayType>(k_MeshFacePath.createChildPath("Ids"));
  const auto* normalsPtr = dataStructure.getDataAs<Float64ArrayType>(k_MeshFacePath.createChildPath("Normals"));
  const auto* centroidsPtr = dataStructure.getDataAs<Float64ArrayType>(k_MeshFacePath.createChildPath("Centroids"));
  REQUIRE(verticesPtr != nullptr);
  REQUIRE(facesPtr != nullptr);
  REQUIRE(labelsPtr != nullptr);
  REQUIRE(idsPtr != nullptr);
  REQUIRE(normalsPtr != nullptr);
  REQUIRE(centroidsPtr != nullptr);
  REQUIRE(verticesPtr->getNumberOfTuples() == 19);
  REQUIRE(facesPtr->getNumberOfTuples() == 9);
  for(usize c = 0; c < 3; ++c)
  {
    CHECK(verticesPtr->getDataStoreRef().getComponentValue(0, c) == 0.0F);
  }
  for(usize t = 0; t < 9; ++t)
  {
    const float64 x = 1.0 + static_cast<float64>(t % 3);
    const float64 y = 1.0 + static_cast<float64>((t / 3) * (t / 3));
    const std::array<float32, 6> expectedVertices = {static_cast<float32>(1.5 * x + y), static_cast<float32>(1.5 * y - x), 0, static_cast<float32>(1.5 * x - y), static_cast<float32>(1.5 * y + x), 0};
    for(usize i = 0; i < expectedVertices.size(); ++i)
    {
      CHECK((*verticesPtr)[(2 * t + 1) * 3 + i] == expectedVertices[i]);
    }
    CHECK(facesPtr->getDataStoreRef().getComponentValue(t, 0) == 0);
    CHECK(facesPtr->getDataStoreRef().getComponentValue(t, 1) == 2 * t + 1);
    CHECK(facesPtr->getDataStoreRef().getComponentValue(t, 2) == 2 * t + 2);
    CHECK(labelsPtr->getDataStoreRef().getComponentValue(t, 0) == 11);
    CHECK(labelsPtr->getDataStoreRef().getComponentValue(t, 1) == 13);
    CHECK((*idsPtr)[t] == 1);
    const std::array<float64, 3> expectedNormal = {0, 0, 1};
    const std::array<float64, 3> expectedCentroid = {x, y, 0};
    for(usize c = 0; c < 3; ++c)
    {
      CHECK(normalsPtr->getDataStoreRef().getComponentValue(t, c) == expectedNormal[c]);
      CHECK(centroidsPtr->getDataStoreRef().getComponentValue(t, c) == expectedCentroid[c]);
    }
  }
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<TriangleGeom>(k_MeshGeomPath));
  const auto& geometry = dataStructure.getDataRefAs<TriangleGeom>(k_MeshGeomPath);
  CHECK(geometry.getVertices() == verticesPtr);
  CHECK(geometry.getFaces() == facesPtr);
}

std::vector<float64> CheckCurvatureOutputs(const DataStructure& dataStructure, bool optionalOutputs)
{
  std::vector<float64> values;
  for(usize i = 0; i < (optionalOutputs ? k_CurvatureOutputKeys.size() : 2); ++i)
  {
    const auto* outputPtr = dataStructure.getDataAs<Float64ArrayType>(CurvatureOutputPath(i));
    REQUIRE(outputPtr != nullptr);
    REQUIRE(outputPtr->getNumberOfTuples() == 9);
    for(usize j = 0; j < outputPtr->getSize(); ++j)
    {
      const float64 value = (*outputPtr)[j];
      // A plane has zero Hessian. Direction signs remain unconstrained.
      if(i < 5)
      {
        CHECK(value == Approx(0.0).margin(1.0e-12));
      }
      else
      {
        CHECK(std::isfinite(value));
      }
      values.push_back(value);
    }
  }
  return values;
}

void CheckCurvatureSave(DataStructure& dataStructure, const SupportSnapshot& incidence)
{
  SupportSaveFile file("curvature");
  const auto outputs = CheckCurvatureOutputs(dataStructure, false);
  // The writer omits supplied connectivity and leaves the caller's objects in place.
  const auto save = hdf5::DataStructureWriter::WriteFile(dataStructure, file.path);
  SIMPLNX_RESULT_REQUIRE_VALID(save);
  RequireSupportSnapshot(dataStructure, SupportKind::VertexIncidence, incidence);
  CheckCurvatureInputs(dataStructure);
  CHECK(CheckCurvatureOutputs(dataStructure, false) == outputs);
  const auto read = hdf5::DataStructureReader::ReadFile(file.path);
  SIMPLNX_RESULT_REQUIRE_VALID(read);
  CheckCurvatureInputs(read.value());
  CHECK(CheckCurvatureOutputs(read.value(), false) == outputs);
  RequireSupportSnapshot(read.value(), SupportKind::VertexIncidence, {});
  RequireSupportSnapshot(read.value(), SupportKind::FaceNeighbors, {});
  unit_test::CheckArraysInheritTupleDims(read.value());
}

void CheckCurvatureSuccess(bool incidence, bool cancelled, bool optionalOutputs, bool save, bool ooc)
{
  DataStructure dataStructure;
  const auto args = CreateCurvatureFixture(dataStructure, incidence, optionalOutputs, !save);
  const auto incidenceBefore = CaptureSupportSnapshot(dataStructure, SupportKind::VertexIncidence);
  const auto neighborsBefore = CaptureSupportSnapshot(dataStructure, SupportKind::FaceNeighbors);
  FeatureFaceCurvatureFilter filter;
  if(cancelled)
  {
    ApplySupportFixtureActions(dataStructure, filter, args);
    const auto result = RunCurvature(dataStructure, true);
    SIMPLNX_RESULT_REQUIRE_VALID(result);
  }
  else
  {
    const auto preflight = filter.preflight(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(preflight.output_actions);
    const auto result = filter.execute(dataStructure, args);
    SIMPLNX_RESULT_REQUIRE_VALID(result.result);
    CheckCurvatureOutputs(dataStructure, optionalOutputs);
  }
  CheckCurvatureInputs(dataStructure);
  RequireSupportSnapshot(dataStructure, SupportKind::FaceNeighbors, neighborsBefore);
  if(ooc)
  {
    std::vector<DataPath> paths = {k_MeshGeomPath.createChildPath("Vertices"), k_MeshGeomPath.createChildPath("Faces"),   k_MeshFacePath.createChildPath("Labels"),
                                   k_MeshFacePath.createChildPath("Ids"),      k_MeshFacePath.createChildPath("Normals"), k_MeshFacePath.createChildPath("Centroids")};
    for(usize i = 0; i < (optionalOutputs ? k_CurvatureOutputKeys.size() : 2); ++i)
    {
      paths.push_back(CurvatureOutputPath(i));
    }
    RequireSupportFixtureOocStores(dataStructure, paths);
  }
  if(save)
  {
    CheckCurvatureSave(dataStructure, incidenceBefore);
  }
  else
  {
    RequireSupportSnapshot(dataStructure, SupportKind::VertexIncidence, incidenceBefore);
  }
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}
} // namespace

TEST_CASE("SimplnxCore::FeatureFaceCurvatureFilter: support ownership", "[SimplnxCore][FeatureFaceCurvatureFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  const bool incidence = GENERATE(false, true);
  CAPTURE(incidence);
  CheckCurvatureSuccess(incidence, false, false, false, false);
}

TEST_CASE("SimplnxCore::FeatureFaceCurvatureFilter: finite planar control", "[SimplnxCore][FeatureFaceCurvatureFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  CheckCurvatureSuccess(false, false, true, false, false);
}

TEST_CASE("SimplnxCore::FeatureFaceCurvatureFilter: cancelled cleanup", "[SimplnxCore][FeatureFaceCurvatureFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  const bool incidence = GENERATE(false, true);
  CAPTURE(incidence);
  CheckCurvatureSuccess(incidence, true, false, false, false);
}

TEST_CASE("SimplnxCore::FeatureFaceCurvatureFilter: support creation failure", "[SimplnxCore][FeatureFaceCurvatureFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  DataStructure dataStructure;
  const auto args = CreateCurvatureFixture(dataStructure, false, false);
  ApplySupportFixtureActions(dataStructure, FeatureFaceCurvatureFilter{}, args);
  for(usize i = 0; i < 2; ++i)
  {
    REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64ArrayType>(CurvatureOutputPath(i)));
    dataStructure.getDataRefAs<Float64ArrayType>(CurvatureOutputPath(i)).fill(73.0);
  }
  auto store = std::make_shared<InjectedFailureStore<SupportIndexType>>(ShapeType{9}, ShapeType{3}, 0);
  for(usize t = 0; t < 9; ++t)
  {
    (*store)[3 * t + 1] = 2 * t + 1;
    (*store)[3 * t + 2] = 2 * t + 2;
  }
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<DataArray<SupportIndexType>>(k_MeshGeomPath.createChildPath("Faces")));
  const auto replacement = dataStructure.getDataRefAs<DataArray<SupportIndexType>>(k_MeshGeomPath.createChildPath("Faces")).setDataStore(store);
  SIMPLNX_RESULT_REQUIRE_VALID(replacement);
  store->reads.clear();
  store->failing_read = 1;
  const auto result = RunCurvature(dataStructure, false);
  RequireInjectedFailure(result, "injected face read failure");
  const std::vector<std::pair<usize, usize>> expectedReads = {{0, 27}};
  REQUIRE(store->reads == expectedReads);
  RequireSupportSnapshot(dataStructure, SupportKind::VertexIncidence, {});
  RequireSupportSnapshot(dataStructure, SupportKind::FaceNeighbors, {});
  for(usize i = 0; i < 2; ++i)
  {
    REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64ArrayType>(CurvatureOutputPath(i)));
    const auto& output = dataStructure.getDataRefAs<Float64ArrayType>(CurvatureOutputPath(i));
    for(usize t = 0; t < output.getSize(); ++t)
    {
      CHECK(output[t] == 73.0);
    }
  }
  CheckCurvatureInputs(dataStructure);
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::FeatureFaceCurvatureFilter: unassociated named child", "[SimplnxCore][FeatureFaceCurvatureFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  DataStructure dataStructure;
  const auto args = CreateCurvatureFixture(dataStructure, false, false);
  CreateSupportData(dataStructure, SupportKind::VertexIncidence, {{42, 7}}, false);
  const auto before = CaptureSupportSnapshot(dataStructure, SupportKind::VertexIncidence);
  ApplySupportFixtureActions(dataStructure, FeatureFaceCurvatureFilter{}, args);
  const auto result = RunCurvature(dataStructure, false);
  SIMPLNX_RESULT_REQUIRE_INVALID(result);
  REQUIRE(result.errors().size() == 1);
  CHECK(result.errors()[0].code == -2231);
  RequireSupportSnapshot(dataStructure, SupportKind::VertexIncidence, before);
  RequireSupportSnapshot(dataStructure, SupportKind::FaceNeighbors, {});
  CheckCurvatureInputs(dataStructure);
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::FeatureFaceCurvatureFilter: retained support save compatibility", "[SimplnxCore][FeatureFaceCurvatureFilter]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceInCore, 1);
  const bool incidence = GENERATE(false, true);
  CAPTURE(incidence);
  CheckCurvatureSuccess(incidence, false, false, true, false);
}

TEST_CASE("SimplnxCore::FeatureFaceCurvatureFilter: real store support ownership", "[SimplnxCore][.OocStoreContract]")
{
  REQUIRE(Application::GetOrCreateInstance() != nullptr);
  REQUIRE(Application::GetOrCreateInstance()->getIOManager("HDF5-OOC") != nullptr);
  const unit_test::PreferencesSentinel preferences(DataStorageMode::ForceOutOfCore, 1);
  const bool incidence = GENERATE(false, true);
  CAPTURE(incidence);
  CheckCurvatureSuccess(incidence, false, false, true, true);
  std::cout << fmt::format(
      "Curvature ownership/save/recovery complete incidence={} vertices=HDF5-OOC faces=HDF5-OOC labels=HDF5-OOC ids=HDF5-OOC normals=HDF5-OOC centroids=HDF5-OOC outputs=HDF5-OOC\n", incidence);
}
