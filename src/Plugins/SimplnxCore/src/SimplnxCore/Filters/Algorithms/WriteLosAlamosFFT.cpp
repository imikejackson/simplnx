#include "WriteLosAlamosFFT.hpp"

#include "simplnx/Common/Constants.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/Utilities/FilterUtilities.hpp"

#include <fstream>

namespace fs = std::filesystem;
using namespace nx::core;

namespace
{
using ull = unsigned long long int;
}

// -----------------------------------------------------------------------------
WriteLosAlamosFFT::WriteLosAlamosFFT(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, WriteLosAlamosFFTInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(mesgHandler)
{
}

// -----------------------------------------------------------------------------
WriteLosAlamosFFT::~WriteLosAlamosFFT() noexcept = default;

// -----------------------------------------------------------------------------
const std::atomic_bool& WriteLosAlamosFFT::getCancel()
{
  return m_ShouldCancel;
}

// -----------------------------------------------------------------------------
Result<> WriteLosAlamosFFT::operator()()
{
  /**
   * Header print function was unimplemented in original filter so it was omitted here and condensed to just the writeFile() function
   */

  // Make sure any directory path is also available as the user may have just typed
  // in a path without actually creating the full path
  Result<> createDirectoriesResult = nx::core::CreateOutputDirectories(m_InputValues->OutputFile.parent_path());
  if(createDirectoriesResult.invalid())
  {
    return createDirectoriesResult;
  }

  std::ofstream file = std::ofstream(m_InputValues->OutputFile, std::ios_base::out | std::ios_base::binary);
  if(!file.is_open())
  {
    return MakeErrorResult(-73450, fmt::format("Error creating and opening output file at path: {}", m_InputValues->OutputFile.string()));
  }

  SizeVec3 dims = m_DataStructure.getDataAs<ImageGeom>(m_InputValues->ImageGeomPath)->getDimensions();

  auto& cellEulerAngles = m_DataStructure.getDataAs<Float32Array>(m_InputValues->CellEulerAnglesArrayPath)->getDataStoreRef();
  auto& cellPhases = m_DataStructure.getDataAs<Int32Array>(m_InputValues->CellPhasesArrayPath)->getDataStoreRef();
  auto& featureIds = m_DataStructure.getDataAs<Int32Array>(m_InputValues->FeatureIdsArrayPath)->getDataStoreRef();

  float phi1 = 0.0f, phi = 0.0f, phi2 = 0.0f;

  auto start = std::chrono::steady_clock::now();
  usize total = dims[0] * dims[1] * dims[2];
  for(usize z = 0; z < dims[2]; ++z)
  {
    for(usize y = 0; y < dims[1]; ++y)
    {
      for(usize x = 0; x < dims[0]; ++x)
      {
        usize index = (z * dims[0] * dims[1]) + (dims[0] * y) + x;
        if(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() > 1000)
        {
          m_MessageHandler(fmt::format("Writing in progress: {} of {} written. Percentage: {}%", index, total, static_cast<uint32>((static_cast<float64>(index) / total) * 100.0)));
          start = std::chrono::steady_clock::now();
        }

        phi1 = cellEulerAngles[index * 3] * 180.0f * Constants::k_1OverPiF;
        phi = cellEulerAngles[index * 3 + 1] * 180.0f * Constants::k_1OverPiF;
        phi2 = cellEulerAngles[index * 3 + 2] * 180.0f * Constants::k_1OverPiF;
        file << fmt::format("{:.3f} {:.3f} {:.3f} {} {} {} {} {}\n", phi1, phi, phi2, static_cast<::ull>(x + 1), static_cast<::ull>(y + 1), static_cast<::ull>(z + 1), featureIds[index],
                            cellPhases[index]);
      }
    }
  }

  file.flush();
  file.close();

  return {};
}
