#include "CellDataAliasTestUtils.hpp"

#include "SimplnxCore/Filters/Algorithms/ErodeDilateBadData.hpp"
#include "SimplnxCore/Filters/ErodeDilateBadDataFilter.hpp"
#include "SimplnxCore/SimplnxCore_test_dirs.hpp"

#include <catch2/catch.hpp>

#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataStore.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Parameters/MultiArraySelectionParameter.hpp"
#include "simplnx/UnitTest/UnitTestCommon.hpp"
#include "simplnx/Utilities/Parsing/HDF5/IO/FileIO.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;
using namespace nx::core;
using namespace nx::core::constants;
using namespace nx::core::unit_test;

namespace
{
constexpr ChoicesParameter::ValueType k_Dilate = 0ULL;
constexpr ChoicesParameter::ValueType k_Erode = 1ULL;

const std::string k_EbsdScanDataName("EBSD Scan Data");

const DataPath k_InputData({"Input Data"});
const DataPath k_EbsdScanDataDataPath = k_InputData.createChildPath(k_EbsdScanDataName);
const DataPath k_FeatureIdsDataPath = k_EbsdScanDataDataPath.createChildPath("FeatureIds");
const StringLiteralType k_MiscData = "Misc";

constexpr usize k_NumTuples = 32;
// ImageGeom dimensions are ordered X, Y, Z. AttributeMatrix and DataArray tuple shapes are ordered
// slowest-to-fastest, i.e. Z, Y, X. Both describe the same 32 cells, so keep them as separate constants
// rather than reusing one that happens to multiply out to the same tuple count.
const SizeVec3Type k_GeometryDimensions{4, 4, 2};
const ShapeType k_TupleShape{k_GeometryDimensions[2], k_GeometryDimensions[1], k_GeometryDimensions[0]};
const DataPath k_DataPath({::k_ImageGeometry, ::k_CellData, k_MiscData});
const DataPath k_ImageFeatureIdsPath({::k_ImageGeometry, ::k_CellData, k_FeatureIds});

using DirectionType = std::array<bool, 3>;
constexpr DirectionType k_XDir{true, false, false};
constexpr DirectionType k_XYDir{true, true, false};
constexpr DirectionType k_XYZDir{true, true, true};
constexpr DirectionType k_XZDir{true, false, true};
constexpr DirectionType k_YDir{false, true, false};
constexpr DirectionType k_YZDir{false, true, true};
constexpr DirectionType k_ZDir{false, false, true};

using ExemplarDataType = std::array<int32, k_NumTuples>;
// Exemplar Dilate data for A/B testing
constexpr ExemplarDataType k_ExemplarFeatureIdsDilateX1{0, 0, 1, 2, 2, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateX1{0, 0, 2, 3, 4, 5, 6, 7, 8, 10, 10, 10, 13, 13, 14, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateX2{0, 0, 0, 2, 2, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 0, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateX2{0, 0, 0, 3, 4, 5, 6, 7, 10, 10, 10, 10, 13, 13, 14, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 31, 31, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateXY1{0, 0, 1, 2, 0, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 0, 5, 6, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateXY1{0, 0, 2, 3, 0, 5, 10, 7, 8, 13, 10, 10, 13, 13, 14, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 31, 28, 29, 31, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateXY2{0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 3, 3, 3, 0, 5, 5, 0, 0, 5, 0, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateXY2{0, 0, 10, 3, 0, 13, 10, 10, 13, 13, 10, 10, 13, 13, 14, 14, 16, 17, 18, 19, 20, 21, 22, 31, 24, 25, 31, 31, 28, 31, 31, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateXYZ1{0, 0, 1, 2, 0, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 3, 3, 3, 3, 5, 5, 0, 0, 5, 0, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateXYZ1{0, 0, 2, 3, 0, 5, 10, 7, 8, 13, 10, 10, 13, 13, 14, 31, 0, 17, 18, 19, 20, 21, 22, 23, 24, 25, 10, 31, 28, 13, 31, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateXYZ2{0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 3, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateXYZ2{0, 0, 10, 3, 0, 13, 10, 10, 13, 13, 10, 10, 13, 13, 14, 31, 0, 0, 18, 19, 0, 21, 10, 31, 24, 13, 10, 31, 13, 13, 31, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateXZ1{0, 0, 1, 2, 2, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 3, 3, 3, 3, 5, 5, 0, 5, 5, 0, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateXZ1{0, 0, 2, 3, 4, 5, 6, 7, 8, 10, 10, 10, 13, 13, 14, 31, 0, 17, 18, 19, 20, 21, 22, 23, 24, 25, 10, 27, 28, 13, 31, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateXZ2{0, 0, 0, 2, 2, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 3, 3, 3, 3, 5, 0, 0, 0, 0, 0, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateXZ2{0, 0, 0, 3, 4, 5, 6, 7, 10, 10, 10, 10, 13, 13, 14, 31, 0, 0, 18, 19, 20, 21, 22, 23, 24, 10, 10, 10, 13, 13, 31, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateY1{0, 1, 1, 2, 0, 1, 0, 2, 1, 0, 0, 2, 2, 0, 0, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 0, 5, 6, 6, 0};
constexpr ExemplarDataType k_ExemplarDataDilateY1{0, 1, 2, 3, 0, 5, 10, 7, 8, 13, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 31, 28, 29, 30, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateY2{0, 1, 0, 2, 0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0, 3, 4, 4, 4, 4, 3, 3, 3, 0, 5, 5, 5, 0, 5, 6, 6, 0};
constexpr ExemplarDataType k_ExemplarDataDilateY2{0, 1, 10, 3, 0, 13, 10, 7, 0, 13, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 31, 24, 25, 26, 31, 28, 29, 30, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateYZ1{0, 1, 1, 2, 0, 1, 0, 2, 1, 0, 0, 2, 2, 0, 0, 0, 0, 4, 4, 4, 3, 3, 3, 3, 5, 5, 0, 0, 5, 0, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateYZ1{0, 1, 2, 3, 0, 5, 10, 7, 8, 13, 10, 11, 12, 13, 14, 31, 0, 17, 18, 19, 20, 21, 22, 23, 24, 25, 10, 31, 28, 13, 14, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateYZ2{0, 1, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 4, 4, 0, 3, 0, 0, 5, 0, 0, 0, 5, 0, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateYZ2{0, 1, 10, 3, 0, 13, 10, 7, 0, 13, 10, 31, 12, 13, 14, 31, 0, 17, 18, 19, 0, 21, 10, 31, 24, 13, 10, 31, 28, 13, 14, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateZ1{0, 1, 1, 2, 2, 1, 2, 2, 1, 1, 0, 2, 2, 0, 0, 0, 0, 4, 4, 4, 3, 3, 3, 3, 5, 5, 0, 5, 5, 0, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateZ1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 31, 0, 17, 18, 19, 20, 21, 22, 23, 24, 25, 10, 27, 28, 13, 14, 31};

constexpr ExemplarDataType k_ExemplarFeatureIdsDilateZ2{0, 1, 1, 2, 2, 1, 2, 2, 1, 1, 0, 2, 2, 0, 0, 0, 0, 4, 4, 4, 3, 3, 3, 3, 5, 5, 0, 5, 5, 0, 0, 0};
constexpr ExemplarDataType k_ExemplarDataDilateZ2{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 31, 0, 17, 18, 19, 20, 21, 22, 23, 24, 25, 10, 27, 28, 13, 14, 31};

// Exemplar Erode data for A/B testing
constexpr ExemplarDataType k_ExemplarFeatureIdsErodeX1{1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 6};
constexpr ExemplarDataType k_ExemplarDataErodeX1{1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 11, 12, 12, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 30};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeX2{1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 6};
constexpr ExemplarDataType k_ExemplarDataErodeX2{1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 11, 12, 12, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 30};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeXY1{1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 5};
constexpr ExemplarDataType k_ExemplarDataErodeXY1{1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 11, 12, 9, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 27};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeXY2{1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 5};
constexpr ExemplarDataType k_ExemplarDataErodeXY2{1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 11, 12, 9, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 27};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeXYZ1{1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 3};
constexpr ExemplarDataType k_ExemplarDataErodeXYZ1{1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 11, 12, 9, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 15};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeXYZ2{1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 3};
constexpr ExemplarDataType k_ExemplarDataErodeXYZ2{1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 11, 12, 9, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 15};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeXZ1{1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 3};
constexpr ExemplarDataType k_ExemplarDataErodeXZ1{1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 11, 12, 12, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 15};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeXZ2{1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 3};
constexpr ExemplarDataType k_ExemplarDataErodeXZ2{1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 11, 12, 12, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 15};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeY1{2, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 0, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 5};
constexpr ExemplarDataType k_ExemplarDataErodeY1{4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 6, 11, 12, 9, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 27};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeY2{2, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 2, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 5};
constexpr ExemplarDataType k_ExemplarDataErodeY2{4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 6, 11, 12, 9, 6, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 27};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeYZ1{2, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 6, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 3};
constexpr ExemplarDataType k_ExemplarDataErodeYZ1{4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 6, 11, 12, 9, 30, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 15};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeYZ2{2, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 6, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 3};
constexpr ExemplarDataType k_ExemplarDataErodeYZ2{4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 6, 11, 12, 9, 30, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 15};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeZ1{4, 1, 1, 2, 2, 1, 2, 2, 1, 1, 5, 2, 2, 6, 6, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 3};
constexpr ExemplarDataType k_ExemplarDataErodeZ1{16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 26, 11, 12, 29, 30, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 15};

constexpr ExemplarDataType k_ExemplarFeatureIdsErodeZ2{4, 1, 1, 2, 2, 1, 2, 2, 1, 1, 5, 2, 2, 6, 6, 3, 4, 4, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 3};
constexpr ExemplarDataType k_ExemplarDataErodeZ2{16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 26, 11, 12, 29, 30, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 15};

/**
 * @brief The expected FeatureIds and Misc output for a single (operation, directions, iterations) combination.
 */
struct ExemplarRecord
{
  ChoicesParameter::ValueType operation;
  DirectionType directions;
  int32 iterations;
  const ExemplarDataType& expected_feature_ids;
  const ExemplarDataType& expected_data;
};

// One row per valid parameter combination: 2 operations x 7 direction combinations x 2 iteration counts.
// The all-directions-off combination is rejected by preflight and therefore has no row.
// clang-format off
const std::array<ExemplarRecord, 28> k_Exemplars{{
    {.operation=k_Dilate, .directions=k_XDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsDilateX1, .expected_data=k_ExemplarDataDilateX1},
    {.operation=k_Dilate, .directions=k_XDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsDilateX2, .expected_data=k_ExemplarDataDilateX2},
    {.operation=k_Dilate, .directions=k_XYDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsDilateXY1, .expected_data=k_ExemplarDataDilateXY1},
    {.operation=k_Dilate, .directions=k_XYDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsDilateXY2, .expected_data=k_ExemplarDataDilateXY2},
    {.operation=k_Dilate, .directions=k_XYZDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsDilateXYZ1, .expected_data=k_ExemplarDataDilateXYZ1},
    {.operation=k_Dilate, .directions=k_XYZDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsDilateXYZ2, .expected_data=k_ExemplarDataDilateXYZ2},
    {.operation=k_Dilate, .directions=k_XZDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsDilateXZ1, .expected_data=k_ExemplarDataDilateXZ1},
    {.operation=k_Dilate, .directions=k_XZDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsDilateXZ2, .expected_data=k_ExemplarDataDilateXZ2},
    {.operation=k_Dilate, .directions=k_YDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsDilateY1, .expected_data=k_ExemplarDataDilateY1},
    {.operation=k_Dilate, .directions=k_YDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsDilateY2, .expected_data=k_ExemplarDataDilateY2},
    {.operation=k_Dilate, .directions=k_YZDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsDilateYZ1, .expected_data=k_ExemplarDataDilateYZ1},
    {.operation=k_Dilate, .directions=k_YZDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsDilateYZ2, .expected_data=k_ExemplarDataDilateYZ2},
    {.operation=k_Dilate, .directions=k_ZDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsDilateZ1, .expected_data=k_ExemplarDataDilateZ1},
    {.operation=k_Dilate, .directions=k_ZDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsDilateZ2, .expected_data=k_ExemplarDataDilateZ2},
    {.operation=k_Erode, .directions=k_XDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsErodeX1, .expected_data=k_ExemplarDataErodeX1},
    {.operation=k_Erode, .directions=k_XDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsErodeX2, .expected_data=k_ExemplarDataErodeX2},
    {.operation=k_Erode, .directions=k_XYDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsErodeXY1, .expected_data=k_ExemplarDataErodeXY1},
    {.operation=k_Erode, .directions=k_XYDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsErodeXY2, .expected_data=k_ExemplarDataErodeXY2},
    {.operation=k_Erode, .directions=k_XYZDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsErodeXYZ1, .expected_data=k_ExemplarDataErodeXYZ1},
    {.operation=k_Erode, .directions=k_XYZDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsErodeXYZ2, .expected_data=k_ExemplarDataErodeXYZ2},
    {.operation=k_Erode, .directions=k_XZDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsErodeXZ1, .expected_data=k_ExemplarDataErodeXZ1},
    {.operation=k_Erode, .directions=k_XZDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsErodeXZ2, .expected_data=k_ExemplarDataErodeXZ2},
    {.operation=k_Erode, .directions=k_YDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsErodeY1, .expected_data=k_ExemplarDataErodeY1},
    {.operation=k_Erode, .directions=k_YDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsErodeY2, .expected_data=k_ExemplarDataErodeY2},
    {.operation=k_Erode, .directions=k_YZDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsErodeYZ1, .expected_data=k_ExemplarDataErodeYZ1},
    {.operation=k_Erode, .directions=k_YZDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsErodeYZ2, .expected_data=k_ExemplarDataErodeYZ2},
    {.operation=k_Erode, .directions=k_ZDir, .iterations=1, .expected_feature_ids=k_ExemplarFeatureIdsErodeZ1, .expected_data=k_ExemplarDataErodeZ1},
    {.operation=k_Erode, .directions=k_ZDir, .iterations=2, .expected_feature_ids=k_ExemplarFeatureIdsErodeZ2, .expected_data=k_ExemplarDataErodeZ2},
}};
// clang-format on

DataStructure CreateTestData()
{
  DataStructure dataStructure;
  auto* geomPtr = ImageGeom::Create(dataStructure, ::k_ImageGeometry);
  geomPtr->setDimensions(k_GeometryDimensions);

  auto* cellDataPtr = AttributeMatrix::Create(dataStructure, ::k_CellData, k_TupleShape, geomPtr->getId());
  geomPtr->setCellData(*cellDataPtr);

  // Feature IDs
  auto featureIdsPtr = std::make_shared<Int32DataStoreType>(k_TupleShape, ShapeType{1}, 0);
  auto* featureIdsArrayPtr = Int32ArrayType::Create(dataStructure, ::k_FeatureIds, featureIdsPtr, cellDataPtr->getId());

  // Index 0, 14, 31
  auto& featureIds = featureIdsArrayPtr->getDataStoreRef();
  featureIds[0] = 0;
  featureIds[1] = 1;
  featureIds[2] = 1;
  featureIds[3] = 2;

  featureIds[4] = 2;
  featureIds[5] = 1;
  featureIds[6] = 2;
  featureIds[7] = 2;

  featureIds[8] = 1;
  featureIds[9] = 1;
  featureIds[10] = 0;
  featureIds[11] = 2;

  featureIds[12] = 2;
  featureIds[13] = 0;
  featureIds[14] = 0;
  featureIds[15] = 3;
  // Z
  featureIds[16] = 4;
  featureIds[17] = 4;
  featureIds[18] = 4;
  featureIds[19] = 4;

  featureIds[20] = 3;
  featureIds[21] = 3;
  featureIds[22] = 3;
  featureIds[23] = 3;

  featureIds[24] = 5;
  featureIds[25] = 5;
  featureIds[26] = 5;
  featureIds[27] = 5;

  featureIds[28] = 5;
  featureIds[29] = 6;
  featureIds[30] = 6;
  featureIds[31] = 0;

  // Misc DataArray
  auto dataStorePtr = std::make_shared<Int32DataStoreType>(k_TupleShape, ShapeType{1}, 0);
  auto* miscArrayPtr = Int32ArrayType::Create(dataStructure, k_MiscData, dataStorePtr, cellDataPtr->getId());

  auto& dataStore = miscArrayPtr->getDataStoreRef();
  for(usize i = 0; i < dataStore.size(); i++)
  {
    dataStore[i] = i;
  }

  return dataStructure;
}

/**
 * @brief Verifies that an array listed in IgnoredDataArrayPaths still holds its original values.
 */
void CheckPathIgnored(const DataStructure& dataStructure)
{
  const auto* ignored = dataStructure.getDataAs<Int32ArrayType>(k_DataPath);
  REQUIRE(ignored != nullptr);
  REQUIRE(ignored->getNumberOfTuples() == k_NumTuples);
  REQUIRE(ignored->getNumberOfComponents() == 1);
  for(usize i = 0; i < k_NumTuples; ++i)
  {
    CHECK((*ignored)[i] == static_cast<int32>(i));
  }
}

/**
 * @brief Compares the generated FeatureIds and Misc arrays against the expected values for this combination
 * of operation, enabled directions, and iteration count.
 */
void CheckOutput(const DataStructure& dataStructure, ChoicesParameter::ValueType operation, const DirectionType& directions, int32 iterations)
{
  const auto exemplarIterPtr = std::find_if(k_Exemplars.cbegin(), k_Exemplars.cend(),
                                            [&](const ExemplarRecord& record) { return record.operation == operation && record.directions == directions && record.iterations == iterations; });
  if(exemplarIterPtr == k_Exemplars.cend())
  {
    FAIL(fmt::format("No expected output is tabulated for operation {} with directions X={} Y={} Z={} and {} iteration(s)", operation, directions[0], directions[1], directions[2], iterations));
    return;
  }

  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Int32ArrayType>(k_ImageFeatureIdsPath));
  const Int32AbstractDataStoreType& featureIds = dataStructure.getDataRefAs<Int32ArrayType>(k_ImageFeatureIdsPath).getDataStoreRef();
  REQUIRE_NOTHROW(dataStructure.getDataRefAs<Int32ArrayType>(k_DataPath));
  const Int32AbstractDataStoreType& dataStore = dataStructure.getDataRefAs<Int32ArrayType>(k_DataPath).getDataStoreRef();

  REQUIRE(featureIds.size() == exemplarIterPtr->expected_feature_ids.size());
  REQUIRE(dataStore.size() == exemplarIterPtr->expected_data.size());

  for(usize i = 0; i < dataStore.size(); i++)
  {
    REQUIRE(featureIds[i] == exemplarIterPtr->expected_feature_ids[i]);
    REQUIRE(dataStore[i] == exemplarIterPtr->expected_data[i]);
  }
}

Arguments MakeErodeDilateArgs(ChoicesParameter::ValueType operation, int32 numIterations, const DirectionType& directions, const DataPath& geometryPath, const DataPath& featureIdsPath,
                              const MultiArraySelectionParameter::ValueType& ignoredPaths = {})
{
  Arguments args;

  args.insertOrAssign(ErodeDilateBadDataFilter::k_Operation_Key, std::make_any<ChoicesParameter::ValueType>(operation));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_NumIterations_Key, std::make_any<int32>(numIterations));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_XDirOn_Key, std::make_any<bool>(directions[0]));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_YDirOn_Key, std::make_any<bool>(directions[1]));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_ZDirOn_Key, std::make_any<bool>(directions[2]));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_CellFeatureIdsArrayPath_Key, std::make_any<DataPath>(featureIdsPath));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_IgnoredDataArrayPaths_Key, std::make_any<MultiArraySelectionParameter::ValueType>(ignoredPaths));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(geometryPath));
  return args;
}

void RunFilter(DataStructure& dataStructure, const Arguments& args)
{
  const ErodeDilateBadDataFilter filter;
  const auto preflight = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(preflight.output_actions);
  const auto execute = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(execute.result);
}

void CheckDirectionalOutputs(ChoicesParameter::ValueType operation)
{
  const auto directions = GENERATE(k_XDir, k_YDir, k_ZDir, k_XYDir, k_XZDir, k_YZDir, k_XYZDir);
  const int32 iterations = GENERATE(1, 2);
  CAPTURE(directions, iterations);
  auto dataStructure = CreateTestData();
  RunFilter(dataStructure, MakeErodeDilateArgs(operation, iterations, directions, DataPath({k_ImageGeometry}), k_ImageFeatureIdsPath));
  CheckOutput(dataStructure, operation, directions, iterations);
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}
} // namespace

TEST_CASE("SimplnxCore::ErodeDilateBadDataFilter(Erode)", "[SimplnxCore][ErodeDilateBadDataFilter]")
{
  unit_test::LoadPlugins();

  const nx::core::unit_test::TestFileSentinel testDataSentinel(nx::core::unit_test::k_TestFilesDir, "6_6_erode_dilate_test.tar.gz", "6_6_erode_dilate_test");

  // Read Exemplar DREAM3D File Filter
  auto exemplarFilePath = fs::path(fmt::format("{}/6_6_erode_dilate_test/6_6_erode_dilate_bad_data.dream3d", unit_test::k_TestFilesDir));
  DataStructure dataStructure = LoadDataStructure(exemplarFilePath);

  RunFilter(dataStructure, MakeErodeDilateArgs(k_Erode, 2, k_XYZDir, k_InputData, k_FeatureIdsDataPath));

// Write the DataStructure out to the file system
#ifdef SIMPLNX_WRITE_TEST_OUTPUT
  WriteTestDataStructure(dataStructure, fs::path(fmt::format("{}/7_0_erode_dilate_bad_data.dream3d", unit_test::k_BinaryTestOutputDir)));
#endif

  const std::string kExemplarDataContainerName("Exemplar Bad Data Erode");
  const DataPath kErodeCellAttributeMatrixDataPath = DataPath({kExemplarDataContainerName, "EBSD Scan Data"});

  unit_test::CompareExemplarToGeneratedData(dataStructure, dataStructure, k_EbsdScanDataDataPath, kExemplarDataContainerName);

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ErodeDilateBadDataFilter(Erode) Expanded", "[SimplnxCore][ErodeDilateBadDataFilter]")
{
  unit_test::LoadPlugins();

  CheckDirectionalOutputs(nx::core::detail::k_ErodeIndex);
}

TEST_CASE("SimplnxCore::ErodeDilateBadDataFilter(Dilate) Expanded", "[SimplnxCore][ErodeDilateBadDataFilter]")
{
  unit_test::LoadPlugins();

  CheckDirectionalOutputs(nx::core::detail::k_DilateIndex);
}

TEST_CASE("SimplnxCore::ErodeDilateBadDataFilter Ignored Path", "[SimplnxCore][ErodeDilateBadDataFilter]")
{
  unit_test::LoadPlugins();

  DataStructure dataStructure = CreateTestData();
  const ChoicesParameter::ValueType operation = GENERATE(k_Dilate, k_Erode);
  RunFilter(dataStructure, MakeErodeDilateArgs(operation, 1, k_XYZDir, DataPath({k_ImageGeometry}), k_ImageFeatureIdsPath, {k_DataPath}));

  CheckPathIgnored(dataStructure);
  const auto& expectedIds = operation == k_Dilate ? k_ExemplarFeatureIdsDilateXYZ1 : k_ExemplarFeatureIdsErodeXYZ1;
  const auto* actualIds = dataStructure.getDataAs<Int32ArrayType>(k_ImageFeatureIdsPath);
  REQUIRE(actualIds != nullptr);
  REQUIRE(actualIds->getSize() == expectedIds.size());
  for(usize i = 0; i < expectedIds.size(); ++i)
  {
    CHECK((*actualIds)[i] == expectedIds[i]);
  }

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ErodeDilateBadDataFilter: Preflight Feature Data Modification Warning", "[SimplnxCore][ErodeDilateBadDataFilter]")
{
  unit_test::LoadPlugins();

  DataStructure dataStructure = CreateTestData();

  const ErodeDilateBadDataFilter filter;
  const auto args = MakeErodeDilateArgs(k_Dilate, 1, k_XYZDir, DataPath({k_ImageGeometry}), k_ImageFeatureIdsPath);

  const auto preflightResult = filter.preflight(dataStructure, args);
  CellDataAliasTest::CheckFeatureDataWarning(preflightResult);
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ErodeDilateBadDataFilter No Direction", "[SimplnxCore][ErodeDilateBadDataFilter]")
{
  unit_test::LoadPlugins();

  DataStructure dataStructure = CreateTestData();
  const DirectionType directions = {false, false, false};
  const ChoicesParameter::ValueType operation = GENERATE(k_Dilate, k_Erode);
  const int32 numIterations = GENERATE(1, 2);

  const ErodeDilateBadDataFilter filter;
  const auto args = MakeErodeDilateArgs(operation, numIterations, directions, DataPath({k_ImageGeometry}), k_ImageFeatureIdsPath);

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.output_actions);

  REQUIRE(preflightResult.output_actions.errors()[0].code == -14601);

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ErodeDilateBadDataFilter No Dimensions", "[SimplnxCore][ErodeDilateBadDataFilter]")
{
  unit_test::LoadPlugins();

  DataStructure dataStructure = CreateTestData();
  const DirectionType directions = {true, true, true};
  const ChoicesParameter::ValueType operation = k_Dilate;
  const int32 numIterations = 1;
  const DataPath geomPath({k_ImageGeometry});

  auto* imageGeomPtr = dataStructure.getDataAs<ImageGeom>(geomPath);
  imageGeomPtr->setDimensions(SizeVec3Type());

  const ErodeDilateBadDataFilter filter;
  const auto args = MakeErodeDilateArgs(operation, numIterations, directions, geomPath, k_ImageFeatureIdsPath);

  // Preflight the filter and check result
  auto preflightResult = filter.preflight(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_INVALID(preflightResult.output_actions);

  REQUIRE(preflightResult.output_actions.errors()[0].code == -14602);

  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ErodeDilateBadDataFilter: Strict Feature Ids", "[SimplnxCore][ErodeDilateBadDataFilter][StrictFeatureIds]")
{
  unit_test::LoadPlugins();
  DataStructure dataStructure;
  CellDataAliasTest::BuildStrictInput(dataStructure);

  ErodeDilateBadDataFilter filter;
  auto args = filter.getDefaultArguments();
  args.insertOrAssign(ErodeDilateBadDataFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(CellDataAliasTest::k_StrictGeometry));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_CellFeatureIdsArrayPath_Key, std::make_any<DataPath>(CellDataAliasTest::k_StrictIds));

  int32 expectedCode = 0;
  SECTION("Correct Cell Attribute Matrix")
  {
    args.insertOrAssign(ErodeDilateBadDataFilter::k_IgnoredDataArrayPaths_Key, std::make_any<std::vector<DataPath>>(std::vector<DataPath>{CellDataAliasTest::k_StrictIds}));
  }
  SECTION("Feature Ids in another Attribute Matrix with the same tuple shape")
  {
    args.insertOrAssign(ErodeDilateBadDataFilter::k_CellFeatureIdsArrayPath_Key, std::make_any<DataPath>(CellDataAliasTest::k_OtherValues));
    expectedCode = -14605;
  }
  SECTION("Ignored array in another Attribute Matrix")
  {
    args.insertOrAssign(ErodeDilateBadDataFilter::k_IgnoredDataArrayPaths_Key, std::make_any<std::vector<DataPath>>(std::vector<DataPath>{CellDataAliasTest::k_OtherValues}));
    expectedCode = -14604;
  }
  const auto result = filter.preflight(dataStructure, args);
  CellDataAliasTest::CheckStrictPreflight(result, expectedCode);
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ErodeDilateBadDataFilter: Ignored Array Cell AM Alias", "[SimplnxCore][ErodeDilateBadDataFilter][StrictFeatureIds]")
{
  unit_test::LoadPlugins();
  DataStructure dataStructure;
  CellDataAliasTest::BuildAliasInput(dataStructure);

  ErodeDilateBadDataFilter filter;
  auto args = filter.getDefaultArguments();
  args.insertOrAssign(ErodeDilateBadDataFilter::k_SelectedImageGeometryPath_Key, std::make_any<DataPath>(CellDataAliasTest::k_AliasGeometry));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_CellFeatureIdsArrayPath_Key, std::make_any<DataPath>(CellDataAliasTest::k_AliasCells.createChildPath("Ids")));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_IgnoredDataArrayPaths_Key, std::make_any<std::vector<DataPath>>(std::vector<DataPath>{CellDataAliasTest::k_AliasCells.createChildPath("Ignored")}));
  args.insertOrAssign(ErodeDilateBadDataFilter::k_NumIterations_Key, std::make_any<int32>(1));
  const auto preflight = filter.preflight(dataStructure, args);
  CellDataAliasTest::CheckAliasPreflight(preflight);
  const auto execute = filter.execute(dataStructure, args);
  SIMPLNX_RESULT_REQUIRE_VALID(execute.result);
  CellDataAliasTest::CheckAliasOutput(dataStructure, {0, 0, 0, 0});
  unit_test::CheckArraysInheritTupleDims(dataStructure);
}

TEST_CASE("SimplnxCore::ErodeDilateBadDataFilter: SIMPL Backwards Compatibility", "[SimplnxCore][ErodeDilateBadDataFilter][BackwardsCompatibility]")
{
  const auto version = GENERATE("6_4", "6_5");
  CAPTURE(version);
  const fs::path fixturePath = fs::path(unit_test::k_SourceDir.view()) / "test" / "simpl_conversion" / version / "ErodeDilateBadDataFilter.json";
  const auto args = unit_test::LoadSIMPLConversionArguments(fixturePath, FilterTraits<ErodeDilateBadDataFilter>::uuid);
  CHECK(args.value<ChoicesParameter::ValueType>(ErodeDilateBadDataFilter::k_Operation_Key) == k_Dilate);
  CHECK(args.value<int32>(ErodeDilateBadDataFilter::k_NumIterations_Key) == 5);
  CHECK(args.value<bool>(ErodeDilateBadDataFilter::k_XDirOn_Key) == true);
  CHECK(args.value<bool>(ErodeDilateBadDataFilter::k_YDirOn_Key) == true);
  CHECK(args.value<bool>(ErodeDilateBadDataFilter::k_ZDirOn_Key) == true);
  CHECK(args.value<DataPath>(ErodeDilateBadDataFilter::k_SelectedImageGeometryPath_Key) == DataPath({"DataContainer"}));
  CHECK(args.value<DataPath>(ErodeDilateBadDataFilter::k_CellFeatureIdsArrayPath_Key) == DataPath({"DataContainer", "CellData", "TestArray"}));
  CHECK(args.value<std::vector<DataPath>>(ErodeDilateBadDataFilter::k_IgnoredDataArrayPaths_Key) ==
        std::vector<DataPath>{DataPath({"DataContainer", "CellData", "Confidence"}), DataPath({"DataContainer", "CellData", "Intensity"})});
}
