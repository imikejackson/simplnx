#pragma once

#include "SimplnxCore/SimplnxCore_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/Filter/IFilter.hpp"
#include "simplnx/Parameters/FileSystemPathParameter.hpp"

#include <filesystem>

namespace fs = std::filesystem;

namespace nx::core
{

struct SIMPLNXCORE_EXPORT ReadSPParksDumpFileInputValues
{
  fs::path InputFilePath;
  DataPath ImageGeometryPath;
  std::string CellAttributeMatrixName;
  std::string DensityArrayName;
  bool OneBasedArrays;
};

struct SIMPLNXCORE_EXPORT SPParksDumpFileCache
{
  fs::path inputFile;
  std::vector<usize> dimensions = {1, 1, 1};
  fs::file_time_type timeStamp;

  void flush()
  {
    inputFile.clear();
    timeStamp = fs::file_time_type();
    dimensions = {1, 1, 1};
  }
};

Result<SPParksDumpFileCache> ReadSPParksDumpFileHeader(const std::filesystem::path& inputFilePath, bool oneBasedArrays);

class SIMPLNXCORE_EXPORT ReadSPParksDumpFile
{
public:
  ReadSPParksDumpFile(DataStructure& dataStructure, const IFilter::MessageHandler& msgHandler, const std::atomic_bool& shouldCancel, ReadSPParksDumpFileInputValues* inputValues);
  ~ReadSPParksDumpFile() noexcept;

  ReadSPParksDumpFile(const ReadSPParksDumpFile&) = delete;
  ReadSPParksDumpFile(ReadSPParksDumpFile&&) noexcept = delete;
  ReadSPParksDumpFile& operator=(const ReadSPParksDumpFile&) = delete;
  ReadSPParksDumpFile& operator=(ReadSPParksDumpFile&&) noexcept = delete;

  Result<> operator()();

  const std::atomic_bool& getCancel();

private:
  DataStructure& m_DataStructure;
  const ReadSPParksDumpFileInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;
};

} // namespace nx::core
