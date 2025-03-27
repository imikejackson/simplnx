#include "ReadSPParksDumpFile.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataGroup.hpp"
#include "simplnx/Utilities/DataArrayUtilities.hpp"
#include "simplnx/Utilities/Parsing/Text/CsvParser.hpp"
#include "simplnx/Utilities/StringUtilities.hpp"

using namespace nx::core;
namespace fs = std::filesystem;

namespace
{
const int32 k_VolBinaryAllocateMismatch = -91504;
} // namespace

namespace nx::core
{
// -----------------------------------------------------------------------------
Result<SPParksDumpFileCache> ReadSPParksDumpFileHeader(const std::filesystem::path& inputFilePath, bool oneBasedArrays)
{
  /*
  ITEM: TIMESTEP
  210    21000.6
  ITEM: NUMBER OF ATOMS
  106480
  ITEM: BOX BOUNDS
  0.5 44.5
  0.5 44.5
  0.5 55.5
  ITEM: ATOMS id type x y z
  */

  int32_t oneBase = 0;
  if(oneBasedArrays)
  {
    oneBase = 1;
  }

  std::ifstream m_InStream(inputFilePath, std::ios::binary);

  // We are going to reuse the 'buf' variable
  std::array<char, 1024> buffer;
  buffer.fill(0);
  if(CsvParser::ReadLine(m_InStream, buffer.data(), 1024) == 0)
  {
  } // ITEM: TIMESTEP
  buffer.fill(0);
  if(CsvParser::ReadLine(m_InStream, buffer.data(), 1024) == 0)
  {
  } // 210    21000.6
  buffer.fill(0);
  if(CsvParser::ReadLine(m_InStream, buffer.data(), 1024) == 0)
  {
  } // ITEM: NUMBER OF ATOMS
  buffer.fill(0);
  if(CsvParser::ReadLine(m_InStream, buffer.data(), 1024) == 0)
  {
  } // 106480

  std::string lineStr(buffer.data());
  lineStr = StringUtilities::simplified(lineStr);

  // Parse out the number of atoms
  int64_t numAtoms = 0;
  Result<int64> result = ConvertTo<int64>::convert(lineStr);
  if(result.invalid())
  {
    return MakeErrorResult<SPParksDumpFileCache>(-98640, "Could not parse numAtoms from file");
  }
  numAtoms = result.value();

  buffer.fill(0);
  if(CsvParser::ReadLine(m_InStream, buffer.data(), 1024) == 0)
  {
  } // ITEM: BOX BOUNDS

  float low = 0.0f;
  float high = 0.0f;

  m_InStream >> low >> high;
  usize nx = static_cast<usize>(floor(high) - ceil(low)) + oneBase;

  m_InStream >> low >> high;
  usize ny = static_cast<usize>(floor(high) - ceil(low)) + oneBase;

  m_InStream >> low >> high;
  usize nz = static_cast<usize>(floor(high) - ceil(low)) + oneBase;

  if(numAtoms != nx * ny * nz)
  {
    return MakeErrorResult<SPParksDumpFileCache>(-260210, fmt::format("Number of sites does not match the calculated number of sites {} != {} * {} * {}", numAtoms, nx, ny, nz));
  }

  SPParksDumpFileCache fileCache;
  fileCache.dimensions = {nx, ny, nz};
  fileCache.inputFile = inputFilePath;
  fileCache.timeStamp = fs::last_write_time(inputFilePath);

  return {std::move(fileCache)};
}
} // namespace nx::core

// -----------------------------------------------------------------------------
ReadSPParksDumpFile::ReadSPParksDumpFile(DataStructure& dataStructure, const IFilter::MessageHandler& msgHandler, const std::atomic_bool& shouldCancel, ReadSPParksDumpFileInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(msgHandler)
{
}

// -----------------------------------------------------------------------------
ReadSPParksDumpFile::~ReadSPParksDumpFile() noexcept = default;

// -----------------------------------------------------------------------------
const std::atomic_bool& ReadSPParksDumpFile::getCancel()
{
  return m_ShouldCancel;
}

// -----------------------------------------------------------------------------
Result<> ReadSPParksDumpFile::operator()()
{
  // const DataPath densityArrayPath = m_InputValues->ImageGeometryPath.createChildPath(m_InputValues->CellAttributeMatrixName).createChildPath(m_InputValues->DensityArrayName);
  // auto& density = m_DataStructure.getDataAs<Float32Array>(densityArrayPath)->getDataStoreRef();
  //
  // auto filesize = static_cast<usize>(fs::file_size(m_InputValues->InputFilePath));
  // const usize allocatedBytes = density.getSize() * sizeof(float32);
  //
  // if(filesize < allocatedBytes)
  // {
  //   return {MakeErrorResult(k_VolBinaryAllocateMismatch, fmt::format("Binary file size ({}) is smaller than the number of allocated bytes ({}).", filesize, allocatedBytes))};
  // }
  //
  // m_MessageHandler(IFilter::Message::Type::Info, "Reading Data from .vol File.....");
  // return ImportFromBinaryFile(m_InputValues->InputFilePath, density);
  return {};
}
