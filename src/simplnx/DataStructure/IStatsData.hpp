#pragma once

#include "simplnx/DataStructure/IArray.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace nx::core
{

namespace StatisticsType
{
inline constexpr uint32_t Feature_SizeVBoverA = 0;      //!<
inline constexpr uint32_t Feature_SizeVCoverA = 1;      //!<
inline constexpr uint32_t Feature_SizeVNeighbors = 2;   //!<
inline constexpr uint32_t Feature_SizeVOmega3 = 3;      //!<
inline constexpr uint32_t Feature_SizeVClustering = 4;  //!<
inline constexpr uint32_t UnknownStatisticsGroup = 999; //!<

} // namespace StatisticsType

namespace DistributionType
{
inline constexpr uint32_t Beta = 0;                    //!<
inline constexpr uint32_t LogNormal = 1;               //!<
inline constexpr uint32_t Power = 2;                   //!<
inline constexpr uint32_t RDFFrequency = 3;            //!<
inline constexpr uint32_t RDFMaxMin = 4;               //!<
inline constexpr uint32_t UnknownDistributionType = 5; //!<
inline constexpr uint32_t Count = 6;                   //!<

enum ColumnCount
{
  BetaColumnCount = 2,        //!<
  LogNormalColumnCount = 2,   //!<
  PowerLawColumnCount = 2,    //!<
  RawDistDataColumnCount = 1, //!<
  UnknownColumCount = 0       //!<
};
} // namespace DistributionType

class IStatsData
{
public:
  IStatsData() = default;
  virtual ~IStatsData() = default;

  IStatsData(const IStatsData&) = default;
  IStatsData(IStatsData&&) = default;

  IStatsData& operator=(const IStatsData&) = default;
  IStatsData& operator=(IStatsData&&) noexcept = default;

  [[nodiscard]] virtual std::string getName() const = 0;
};

class PrimaryStatsData : public IStatsData
{
public:
  PrimaryStatsData() = default;
  ~PrimaryStatsData() override = default;

  PrimaryStatsData(const PrimaryStatsData&) = default;
  PrimaryStatsData(PrimaryStatsData&&) = default;

  PrimaryStatsData& operator=(const PrimaryStatsData&) = default;
  PrimaryStatsData& operator=(PrimaryStatsData&&) noexcept = default;

  [[nodiscard]] std::string getName() const override
  {
    return "PrimaryStatsData";
  }

  using Float32Array = std::vector<float>;
  using VectorOfFloatArray = std::vector<Float32Array>;

private:
  std::array<float, 3> m_FeatureDiameterInfo;
  float m_BoundaryArea = {};
  VectorOfFloatArray m_FeatureSizeDistribution = {};
  uint32_t m_FeatureSize_DistType = nx::core::DistributionType::LogNormal;
  Float32Array m_BinNumbers = {};
  VectorOfFloatArray m_FeatureSize_BOverA = {};
  uint32_t m_BOverA_DistType = nx::core::DistributionType::Beta;
  VectorOfFloatArray m_FeatureSize_COverA = {};
  uint32_t m_COverA_DistType = nx::core::DistributionType::Beta;
  VectorOfFloatArray m_FeatureSize_Neighbors = {};
  uint32_t m_Neighbors_DistType = nx::core::DistributionType::LogNormal;
  VectorOfFloatArray m_FeatureSize_Omegas = {};
  uint32_t m_Omegas_DistType = nx::core::DistributionType::Beta;
  Float32Array m_MisorientationBins = {};
  VectorOfFloatArray m_MDF_Weights = {};
  Float32Array m_ODF = {};
  VectorOfFloatArray m_ODF_Weights = {};
  Float32Array m_AxisOrientation = {};
  VectorOfFloatArray m_AxisODF_Weights = {};
};

class PrecipitateStatsData : public IStatsData
{
public:
  PrecipitateStatsData() = default;
  ~PrecipitateStatsData() override = default;

  PrecipitateStatsData(const PrecipitateStatsData&) = default;
  PrecipitateStatsData(PrecipitateStatsData&&) = default;

  PrecipitateStatsData& operator=(const PrecipitateStatsData&) = default;
  PrecipitateStatsData& operator=(PrecipitateStatsData&&) noexcept = default;
  [[nodiscard]] std::string getName() const override
  {
    return "PrecipitateStatsData";
  }
};

class MatrixStatsData : public IStatsData
{
public:
  MatrixStatsData() = default;
  ~MatrixStatsData() override = default;

  MatrixStatsData(const MatrixStatsData&) = default;
  MatrixStatsData(MatrixStatsData&&) = default;

  MatrixStatsData& operator=(const MatrixStatsData&) = default;
  MatrixStatsData& operator=(MatrixStatsData&&) noexcept = default;
  [[nodiscard]] std::string getName() const override
  {
    return "MatrixStatsData";
  }
};

class BoundaryStatsData : public IStatsData
{
public:
  BoundaryStatsData() = default;
  ~BoundaryStatsData() override = default;

  BoundaryStatsData(const BoundaryStatsData&) = default;
  BoundaryStatsData(BoundaryStatsData&&) = default;

  BoundaryStatsData& operator=(const BoundaryStatsData&) = default;
  BoundaryStatsData& operator=(BoundaryStatsData&&) noexcept = default;
  [[nodiscard]] std::string getName() const override
  {
    return "BoundaryStatsData";
  }
};

class TransformationStatsData : public IStatsData
{
public:
  TransformationStatsData() = default;
  ~TransformationStatsData() override = default;

  TransformationStatsData(const TransformationStatsData&) = default;
  TransformationStatsData(TransformationStatsData&&) = default;

  TransformationStatsData& operator=(const TransformationStatsData&) = default;
  TransformationStatsData& operator=(TransformationStatsData&&) noexcept = default;
  [[nodiscard]] std::string getName() const override
  {
    return "TransformationStatsData";
  }
};

class StatsDataArray : public IArray
{

public:
  using value_type = std::unique_ptr<IStatsData>;

  using collection_type = std::vector<value_type>;
  using reference = value_type&;
  using const_reference = const value_type&;
  using iterator = typename collection_type::iterator;
  using const_iterator = typename collection_type::const_iterator;

  StatsDataArray() = delete;

  ~StatsDataArray() override = default;

  StatsDataArray(const StatsDataArray&) = default;
  StatsDataArray(StatsDataArray&&) = default;

  StatsDataArray& operator=(const StatsDataArray&) = delete;
  StatsDataArray& operator=(StatsDataArray&&) noexcept = delete;

  static StatsDataArray* Create(DataStructure& dataStructure, const std::string_view& name, const std::optional<IdType>& parentId = {});

  DataObject* shallowCopy() override;
  std::shared_ptr<DataObject> deepCopy(const DataPath& copyPath) override;
  DataObject::Type getDataObjectType() const override;
  std::string getTypeName() const override;

  /**
   * @brief This method sets the shape of the dimensions to `tupleShape`.
   *
   * There are 3 possibilities when using this function:
   * [1] The number of tuples of the new shape is *LESS* than the original. In this
   * case a memory allocation will take place and the first 'N' elements of data
   * will be copied into the new array. The remaining data is *LOST*
   *
   * [2] The number of tuples of the new shape is *EQUAL* to the original. In this
   * case the shape is set and the function returns.
   *
   * [3] The number of tuples of the new shape is *GREATER* than the original. In
   * this case a new array is allocated and all the data from the original array
   * is copied into the new array and the remaining elements are initialized to
   * the default initialization value.
   *
   * @param tupleShape The new shape of the data where the dimensions are "C" ordered
   * from *slowest* to *fastest*.
   */
  void resizeTuples(const std::vector<usize>& tupleShape) override;

  /**
   * @brief Returns an enumeration of the class or subclass. Used for quick comparison or type deduction
   * @return
   */
  ArrayType getArrayType() const override;

  /**
   * @brief Returns the number of elements.
   * @return usize
   */
  usize getSize() const override
  {
    return m_Arrays.size();
  };

  /**
   * @brief Returns the number of elements.
   * @return usize
   */
  usize size() const override
  {
    return m_Arrays.size();
  };

  /**
   * @brief Returns if there are any elements in the array object
   * @return bool, true if the DataArray has a size() == 0
   */
  bool empty() const override
  {
    return m_Arrays.empty();
  };

  /**
   * @brief Returns the tuple shape.
   * @return
   */
  ShapeType getTupleShape() const override
  {
    return ShapeType{m_Arrays.size()};
  };

  /**
   * @brief Returns the component shape.
   * @return
   */
  ShapeType getComponentShape() const override
  {
    return ShapeType{1};
  };

  /**
   * @brief Returns the number of tuples.
   * @return usize
   */
  usize getNumberOfTuples() const override
  {
    return m_Arrays.size();
  };

  /**
   * @brief Returns the number of components per tuple.
   * @return usize
   */
  usize getNumberOfComponents() const override
  {
    return 1;
  }

  void addShape(std::unique_ptr<IStatsData> stats)
  {
    m_Arrays.push_back(std::move(stats)); // Store in vector
  }

protected:
  StatsDataArray(DataStructure& dataStructure, std::string name)
  : IArray(dataStructure, std::move(name))
  {
  }

private:
  std::vector<std::unique_ptr<IStatsData>> m_Arrays; // Vector of unique_ptr to Shape
};

} // namespace nx::core
