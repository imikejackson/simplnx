#include "OrientationAnalysis/Filters/ComputeMDFFilter.hpp"
#include "OrientationAnalysis/OrientationAnalysis_test_dirs.hpp"

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/Common/Constants.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"

#include <catch2/catch.hpp>

#include <EbsdLib/LaueOps/LaueOps.h>
#include <EbsdLib/Orientation/Euler.hpp>
#include <EbsdLib/Orientation/Quaternion.hpp>
#include <EbsdLib/Orientation/Rodrigues.hpp>

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

// Builds an 8x8x1 single-slice ImageGeom. When splitFeatures is true the left half (x < 4) is
// feature 1 with Euler {0,0,0} and the right half (x >= 4) is feature 2 with Euler {45deg,0,0}
// (a 45 degree rotation about the sample Z axis). When false every cell is a single feature with
// identical orientation, i.e. no feature boundaries exist.
DataStructure CreateBicrystalDataStructure(bool splitFeatures)
{
  DataStructure dataStructure;

  auto* imageGeom = ImageGeom::Create(dataStructure, "ImageGeom");
  imageGeom->setDimensions({8, 8, 1});

  const usize numCells = 64;
  auto* cellAM = AttributeMatrix::Create(dataStructure, "Cell Data", {numCells}, imageGeom->getId());
  auto* eulerAngles = UnitTest::CreateTestDataArray<float32>(dataStructure, "EulerAngles", {numCells}, {3}, cellAM->getId());
  auto* phases = UnitTest::CreateTestDataArray<int32>(dataStructure, "Phases", {numCells}, {1}, cellAM->getId());
  auto* featureIds = UnitTest::CreateTestDataArray<int32>(dataStructure, "FeatureIds", {numCells}, {1}, cellAM->getId());

  const auto fortyFiveDegrees = static_cast<float32>(45.0 * Constants::k_PiOver180D);
  for(usize y = 0; y < 8; y++)
  {
    for(usize x = 0; x < 8; x++)
    {
      const usize cellIndex = y * 8 + x;
      (*phases)[cellIndex] = 1;
      const bool secondFeature = splitFeatures && (x >= 4);
      (*featureIds)[cellIndex] = secondFeature ? 2 : 1;
      (*eulerAngles)[cellIndex * 3 + 0] = secondFeature ? fortyFiveDegrees : 0.0f;
      (*eulerAngles)[cellIndex * 3 + 1] = 0.0f;
      (*eulerAngles)[cellIndex * 3 + 2] = 0.0f;
    }
  }

  auto* ensembleAM = AttributeMatrix::Create(dataStructure, "Cell Ensemble Data", {2}, imageGeom->getId());
  auto* crystalStructures = UnitTest::CreateTestDataArray<uint32>(dataStructure, "CrystalStructures", {2}, {1}, ensembleAM->getId());
  (*crystalStructures)[0] = 999;
  (*crystalStructures)[1] = 1;

  return dataStructure;
}

// Fills the given Arguments with the standard set of paths for the bicrystal DataStructure.
Arguments CreateBicrystalArgs(float32 halfwidth, int32 numCurvePoints)
{
  DataPath imageGeomPath({"ImageGeom"});
  DataPath cellDataPath = imageGeomPath.createChildPath("Cell Data");
  DataPath outputGroupPath({"MDF Data"});

  Arguments args;
  args.insertOrAssign(ComputeMDFFilter::k_ImageGeometryPath_Key, std::make_any<GeometrySelectionParameter::ValueType>(imageGeomPath));
  args.insertOrAssign(ComputeMDFFilter::k_CellEulerAnglesArrayPath_Key, std::make_any<ArraySelectionParameter::ValueType>(cellDataPath.createChildPath("EulerAngles")));
  args.insertOrAssign(ComputeMDFFilter::k_CellPhasesArrayPath_Key, std::make_any<ArraySelectionParameter::ValueType>(cellDataPath.createChildPath("Phases")));
  args.insertOrAssign(ComputeMDFFilter::k_FeatureIdsArrayPath_Key, std::make_any<ArraySelectionParameter::ValueType>(cellDataPath.createChildPath("FeatureIds")));
  args.insertOrAssign(ComputeMDFFilter::k_CrystalStructuresArrayPath_Key,
                      std::make_any<ArraySelectionParameter::ValueType>(imageGeomPath.createChildPath("Cell Ensemble Data").createChildPath("CrystalStructures")));
  args.insertOrAssign(ComputeMDFFilter::k_Halfwidth_Key, std::make_any<Float32Parameter::ValueType>(halfwidth));
  args.insertOrAssign(ComputeMDFFilter::k_NumCurvePoints_Key, std::make_any<Int32Parameter::ValueType>(numCurvePoints));
  args.insertOrAssign(ComputeMDFFilter::k_UseAreaWeights_Key, std::make_any<BoolParameter::ValueType>(false));
  args.insertOrAssign(ComputeMDFFilter::k_OutputGroupPath_Key, std::make_any<DataGroupCreationParameter::ValueType>(outputGroupPath));
  return args;
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

TEST_CASE("OrientationAnalysis::ComputeMDF: Cubic bicrystal", "[OrientationAnalysis][ComputeMDF]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure = CreateBicrystalDataStructure(true);

  ComputeMDFFilter filter;
  Arguments args = CreateBicrystalArgs(10.0f, 200);

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  DataPath outputGroupPath({"MDF Data"});
  DataPath phaseGroupPath = outputGroupPath.createChildPath("Phase-1");
  DataPath angleDistPath = phaseGroupPath.createChildPath("Angle Distribution");

  // Independently compute the MDF bin that EbsdLib assigns the 45deg@[001] misorientation using
  // the SAME quaternion convention as the algorithm (miso = q1.conjugate() * q2).
  auto orientationOps = ebsdlib::LaueOps::GetAllOrientationOps();
  ebsdlib::LaueOps::Pointer cubicOps = orientationOps[1];
  const double fortyFiveDegrees = 45.0 * Constants::k_PiOver180D;
  ebsdlib::QuatD quat1 = ebsdlib::EulerDType(0.0, 0.0, 0.0).toQuaternion();
  ebsdlib::QuatD quat2 = ebsdlib::EulerDType(fortyFiveDegrees, 0.0, 0.0).toQuaternion();
  ebsdlib::QuatD misoQuat = quat1.conjugate() * quat2;
  const int expectedBin = cubicOps->getMisoBin(cubicOps->getMDFFZRod(misoQuat.toRodrigues()));
  const usize expectedMdfSize = cubicOps->getMDFSize();
  REQUIRE(expectedMdfSize == 5832);

  // 1. MDF array resized to the CubicOps MDF size.
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64Array>(phaseGroupPath.createChildPath("MDF")));
  const auto& mdfArray = dataStructure.getDataRefAs<Float64Array>(phaseGroupPath.createChildPath("MDF"));
  REQUIRE(mdfArray.getNumberOfTuples() == expectedMdfSize);

  // 2. argmax(MDF) equals the independently computed bin.
  usize mdfArgMax = 0;
  float64 mdfMax = mdfArray[0];
  for(usize i = 1; i < mdfArray.getNumberOfTuples(); i++)
  {
    if(mdfArray[i] > mdfMax)
    {
      mdfMax = mdfArray[i];
      mdfArgMax = i;
    }
  }
  // The 18^3 = 5832 cell MDF grid covers the full Rodrigues cube, but many cells lie outside the
  // misorientation fundamental zone and fold (via getMDFFZRod, applied inside binCenter /
  // determineRodriguesVector) onto the SAME physical misorientation as an in-FZ cell. The KDE peak
  // therefore lands on whichever cell center folds closest to the input misorientation, which is
  // not necessarily the raw accumulation index. The physically meaningful check is that the peak
  // bin, folded to the FZ and re-binned, equals the canonical bin EbsdLib assigns the input
  // misorientation.
  double binCenterSeed[3] = {0.5, 0.5, 0.5};
  ebsdlib::RodriguesDType peakRod = cubicOps->determineRodriguesVector(binCenterSeed, static_cast<int>(mdfArgMax));
  const int canonicalPeakBin = cubicOps->getMisoBin(peakRod);
  REQUIRE(canonicalPeakBin == expectedBin);

  // 3. The measured "MDF Density" curve peaks near 45 degrees (Angles are in DEGREES).
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64Array>(angleDistPath.createChildPath("Angles")));
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64Array>(angleDistPath.createChildPath("MDF Density")));
  const auto& anglesArray = dataStructure.getDataRefAs<Float64Array>(angleDistPath.createChildPath("Angles"));
  const auto& mdfDensityArray = dataStructure.getDataRefAs<Float64Array>(angleDistPath.createChildPath("MDF Density"));

  usize densityArgMax = 0;
  float64 densityMax = mdfDensityArray[0];
  for(usize i = 1; i < mdfDensityArray.getNumberOfTuples(); i++)
  {
    if(mdfDensityArray[i] > densityMax)
    {
      densityMax = mdfDensityArray[i];
      densityArgMax = i;
    }
  }
  const float64 curveStep = anglesArray[1] - anglesArray[0];
  REQUIRE(std::abs(anglesArray[densityArgMax] - 45.0) <= 2.0 * curveStep);

  // 4. The random (Mackenzie) reference distribution has unit mean.
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64Array>(angleDistPath.createChildPath("Random Density")));
  const auto& randomDensityArray = dataStructure.getDataRefAs<Float64Array>(angleDistPath.createChildPath("Random Density"));
  float64 randomSum = 0.0;
  for(usize i = 0; i < randomDensityArray.getNumberOfTuples(); i++)
  {
    randomSum += randomDensityArray[i];
  }
  const float64 randomMean = randomSum / static_cast<float64>(randomDensityArray.getNumberOfTuples());
  REQUIRE(std::abs(randomMean - 1.0) < 1.0e-6);

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("OrientationAnalysis::ComputeMDF: no boundaries warns and zero-fills", "[OrientationAnalysis][ComputeMDF]")
{
  UnitTest::LoadPlugins();

  DataStructure dataStructure = CreateBicrystalDataStructure(false);

  ComputeMDFFilter filter;
  Arguments args = CreateBicrystalArgs(10.0f, 200);

  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  // Execution must SUCCEED even though there are no feature boundaries.
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  DataPath outputGroupPath({"MDF Data"});
  DataPath phaseGroupPath = outputGroupPath.createChildPath("Phase-1");
  DataPath angleDistPath = phaseGroupPath.createChildPath("Angle Distribution");

  // MDF array is all zeros.
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64Array>(phaseGroupPath.createChildPath("MDF")));
  const auto& mdfArray = dataStructure.getDataRefAs<Float64Array>(phaseGroupPath.createChildPath("MDF"));
  for(usize i = 0; i < mdfArray.getNumberOfTuples(); i++)
  {
    REQUIRE(mdfArray[i] == 0.0);
  }

  // Measured "MDF Density" is zero everywhere.
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64Array>(angleDistPath.createChildPath("MDF Density")));
  const auto& mdfDensityArray = dataStructure.getDataRefAs<Float64Array>(angleDistPath.createChildPath("MDF Density"));
  for(usize i = 0; i < mdfDensityArray.getNumberOfTuples(); i++)
  {
    REQUIRE(mdfDensityArray[i] == 0.0);
  }

  // The random (Mackenzie) reference distribution is still populated with unit mean.
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Float64Array>(angleDistPath.createChildPath("Random Density")));
  const auto& randomDensityArray = dataStructure.getDataRefAs<Float64Array>(angleDistPath.createChildPath("Random Density"));
  float64 randomSum = 0.0;
  for(usize i = 0; i < randomDensityArray.getNumberOfTuples(); i++)
  {
    randomSum += randomDensityArray[i];
  }
  const float64 randomMean = randomSum / static_cast<float64>(randomDensityArray.getNumberOfTuples());
  REQUIRE(std::abs(randomMean - 1.0) < 1.0e-6);

  UnitTest::CheckArraysInheritTupleDims(dataStructure);
}
