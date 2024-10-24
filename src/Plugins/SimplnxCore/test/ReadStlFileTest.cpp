#include "SimplnxCore/Filters/ReadStlFileFilter.hpp"
#include "SimplnxCore/SimplnxCore_test_dirs.hpp"

#include "simplnx/DataStructure/Geometry/TriangleGeom.hpp"
#include "simplnx/Parameters/ArrayCreationParameter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Utilities/Parsing/HDF5/Writers/FileWriter.hpp"

#include <catch2/catch.hpp>

#include <filesystem>
#include <string>

namespace fs = std::filesystem;
using namespace nx::core;
using namespace nx::core::Constants;

TEST_CASE("SimplnxCore::ReadStlFileFilter:Valid_File", "[SimplnxCore][ReadStlFileFilter]")
{
  const nx::core::UnitTest::TestFileSentinel testDataSentinel(nx::core::unit_test::k_CMakeExecutable, nx::core::unit_test::k_TestFilesDir, "ReadSTLFileTest.tar.gz", "ReadSTLFileTest");

  // Instantiate the filter, a DataStructure object and an Arguments Object
  DataStructure dataStructure;
  Arguments args;
  ReadStlFileFilter filter;

  DataPath triangleGeomDataPath({"[Triangle Geometry]"});

  std::string inputFile = fmt::format("{}/ReadSTLFileTest/ASTMD638_specimen.stl", unit_test::k_TestFilesDir);

  // Create default Parameters for the filter.
  args.insertOrAssign(ReadStlFileFilter::k_StlFilePath_Key, std::make_any<FileSystemPathParameter::ValueType>(fs::path(inputFile)));
  args.insertOrAssign(ReadStlFileFilter::k_CreatedTriangleGeometryPath_Key, std::make_any<DataPath>(triangleGeomDataPath));

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  // Execute the filter and check the result
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(executeResult.result);

  TriangleGeom& triangleGeom = dataStructure.getDataRefAs<TriangleGeom>(triangleGeomDataPath);
  REQUIRE(triangleGeom.getNumberOfFaces() == 92);
  REQUIRE(triangleGeom.getNumberOfVertices() == 48);

  // Write the DataStructure out to the file system
#ifdef SIMPLNX_WRITE_TEST_OUTPUT
  WriteTestDataStructure(dataStructure, fs::path(fmt::format("{}/StlFileReaderTest.dream3d", unit_test::k_BinaryTestOutputDir)));
#endif
}

TEST_CASE("SimplnxCore::ReadStlFileFilter:STLParseError", "[SimplnxCore][ReadStlFileFilter]")
{
  const nx::core::UnitTest::TestFileSentinel testDataSentinel(nx::core::unit_test::k_CMakeExecutable, nx::core::unit_test::k_TestFilesDir, "ReadSTLFileTest.tar.gz", "ReadSTLFileTest");

  // Instantiate the filter, a DataStructure object and an Arguments Object
  DataStructure dataStructure;
  Arguments args;
  ReadStlFileFilter filter;

  DataPath triangleGeomDataPath({"[Triangle Geometry]"});

  std::string inputFile = fmt::format("{}/ReadSTLFileTest/stl_test_wrong_num_triangles.stl", unit_test::k_TestFilesDir);

  // Create default Parameters for the filter.
  args.insertOrAssign(ReadStlFileFilter::k_StlFilePath_Key, std::make_any<FileSystemPathParameter::ValueType>(fs::path(inputFile)));
  args.insertOrAssign(ReadStlFileFilter::k_CreatedTriangleGeometryPath_Key, std::make_any<DataPath>(triangleGeomDataPath));

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  // Execute the filter and check the result
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(executeResult.result);

  REQUIRE(executeResult.result.errors().front().code == -1108);
}

TEST_CASE("SimplnxCore::ReadStlFileFilter:TriangleParseError", "[SimplnxCore][ReadStlFileFilter]")
{
  const nx::core::UnitTest::TestFileSentinel testDataSentinel(nx::core::unit_test::k_CMakeExecutable, nx::core::unit_test::k_TestFilesDir, "ReadSTLFileTest.tar.gz", "ReadSTLFileTest");

  // Instantiate the filter, a DataStructure object and an Arguments Object
  DataStructure dataStructure;
  Arguments args;
  ReadStlFileFilter filter;

  DataPath triangleGeomDataPath({"[Triangle Geometry]"});

  std::string inputFile = fmt::format("{}/ReadSTLFileTest/stl_test_2_TriangleParseError.stl", unit_test::k_TestFilesDir);

  // Create default Parameters for the filter.
  args.insertOrAssign(ReadStlFileFilter::k_StlFilePath_Key, std::make_any<FileSystemPathParameter::ValueType>(fs::path(inputFile)));
  args.insertOrAssign(ReadStlFileFilter::k_CreatedTriangleGeometryPath_Key, std::make_any<DataPath>(triangleGeomDataPath));

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  // Execute the filter and check the result
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(executeResult.result);

  REQUIRE(executeResult.result.errors().front().code == -1106);
}

TEST_CASE("SimplnxCore::ReadStlFileFilter:AttributeParseError", "[SimplnxCore][ReadStlFileFilter]")
{
  const nx::core::UnitTest::TestFileSentinel testDataSentinel(nx::core::unit_test::k_CMakeExecutable, nx::core::unit_test::k_TestFilesDir, "ReadSTLFileTest.tar.gz", "ReadSTLFileTest");

  // Instantiate the filter, a DataStructure object and an Arguments Object
  DataStructure dataStructure;
  Arguments args;
  ReadStlFileFilter filter;

  DataPath triangleGeomDataPath({"[Triangle Geometry]"});

  std::string inputFile = fmt::format("{}/ReadSTLFileTest/stl_test_2_AttributeParseError.stl", unit_test::k_TestFilesDir);

  // Create default Parameters for the filter.
  args.insertOrAssign(ReadStlFileFilter::k_StlFilePath_Key, std::make_any<FileSystemPathParameter::ValueType>(fs::path(inputFile)));
  args.insertOrAssign(ReadStlFileFilter::k_CreatedTriangleGeometryPath_Key, std::make_any<DataPath>(triangleGeomDataPath));

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflightResult.outputActions);

  // Execute the filter and check the result
  auto executeResult = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(executeResult.result);

  REQUIRE(executeResult.result.errors().front().code == -1107);
}
