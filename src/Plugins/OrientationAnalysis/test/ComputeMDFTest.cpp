#include "OrientationAnalysis/Filters/ComputeMDFFilter.hpp"
#include "OrientationAnalysis/OrientationAnalysis_test_dirs.hpp"

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <catch2/catch.hpp>

using namespace nx::core;

namespace
{
DataStructure CreateDataStructure()
{
  DataStructure dataStructure;

  auto* imageGeom = ImageGeom::Create(dataStructure, "ImageGeom");
  imageGeom->setDimensions({4, 4, 1});

  auto* cellAM = AttributeMatrix::Create(dataStructure, "Cell Data", {16}, imageGeom->getId());
  UnitTest::CreateTestDataArray<float32>(dataStructure, "EulerAngles", {16}, {3}, cellAM->getId());
  UnitTest::CreateTestDataArray<int32>(dataStructure, "Phases", {16}, {1}, cellAM->getId());
  UnitTest::CreateTestDataArray<int32>(dataStructure, "FeatureIds", {16}, {1}, cellAM->getId());

  auto* ensembleAM = AttributeMatrix::Create(dataStructure, "Cell Ensemble Data", {2}, imageGeom->getId());
  auto* crystalStructures = UnitTest::CreateTestDataArray<uint32>(dataStructure, "CrystalStructures", {2}, {1}, ensembleAM->getId());
  (*crystalStructures)[0] = 999;
  (*crystalStructures)[1] = 1;

  return dataStructure;
}
} // namespace

TEST_CASE("OrientationAnalysis::ComputeMDF: Preflight", "[OrientationAnalysis][ComputeMDF]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure = CreateDataStructure();

  DataPath imageGeomPath({"ImageGeom"});
  DataPath eulerAnglesPath = imageGeomPath.createChildPath("Cell Data").createChildPath("EulerAngles");
  DataPath phasesPath = imageGeomPath.createChildPath("Cell Data").createChildPath("Phases");
  DataPath featureIdsPath = imageGeomPath.createChildPath("Cell Data").createChildPath("FeatureIds");
  DataPath crystalStructuresPath = imageGeomPath.createChildPath("Cell Ensemble Data").createChildPath("CrystalStructures");
  DataPath outputGroupPath({"MDF Data"});

  ComputeMDFFilter filter;
  Arguments args;

  args.insertOrAssign(ComputeMDFFilter::k_ImageGeometryPath_Key, std::make_any<GeometrySelectionParameter::ValueType>(imageGeomPath));
  args.insertOrAssign(ComputeMDFFilter::k_CellEulerAnglesArrayPath_Key, std::make_any<ArraySelectionParameter::ValueType>(eulerAnglesPath));
  args.insertOrAssign(ComputeMDFFilter::k_CellPhasesArrayPath_Key, std::make_any<ArraySelectionParameter::ValueType>(phasesPath));
  args.insertOrAssign(ComputeMDFFilter::k_FeatureIdsArrayPath_Key, std::make_any<ArraySelectionParameter::ValueType>(featureIdsPath));
  args.insertOrAssign(ComputeMDFFilter::k_CrystalStructuresArrayPath_Key, std::make_any<ArraySelectionParameter::ValueType>(crystalStructuresPath));
  args.insertOrAssign(ComputeMDFFilter::k_Halfwidth_Key, std::make_any<Float32Parameter::ValueType>(5.0f));
  args.insertOrAssign(ComputeMDFFilter::k_NumCurvePoints_Key, std::make_any<Int32Parameter::ValueType>(13));
  args.insertOrAssign(ComputeMDFFilter::k_UseAreaWeights_Key, std::make_any<BoolParameter::ValueType>(false));
  args.insertOrAssign(ComputeMDFFilter::k_OutputGroupPath_Key, std::make_any<DataGroupCreationParameter::ValueType>(outputGroupPath));

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  REQUIRE_NOTHROW(dataStructure.getDataRefAs<IDataArray>(outputGroupPath.createChildPath("Phase-1").createChildPath("MDF")));
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<IDataArray>(outputGroupPath.createChildPath("Phase-1").createChildPath("Angle Distribution").createChildPath("Angles")));
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<IDataArray>(outputGroupPath.createChildPath("Phase-1").createChildPath("Angle Distribution").createChildPath("MDF Density")));
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<IDataArray>(outputGroupPath.createChildPath("Phase-1").createChildPath("Angle Distribution").createChildPath("Random Density")));

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}
