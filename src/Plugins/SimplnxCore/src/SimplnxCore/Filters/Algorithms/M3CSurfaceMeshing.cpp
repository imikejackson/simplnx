#include "M3CSurfaceMeshing.hpp"

#include "SimplnxCore/Filters/Algorithms/TupleTransfer.hpp"

#include "simplnx/Common/Array.hpp"
#include "simplnx/Common/Range.hpp"
#include "simplnx/Common/Result.hpp"
#include "simplnx/Common/Types.hpp"
#include "simplnx/DataStructure/AbstractDataStore.hpp"
#include "simplnx/DataStructure/AttributeMatrix.hpp"
#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/DataStructure/Geometry/IGeometry.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"
#include "simplnx/DataStructure/Geometry/TriangleGeom.hpp"
#include "simplnx/DataStructure/IArray.hpp"
#include "simplnx/DataStructure/IDataArray.hpp"
#include "simplnx/DataStructure/IO/Generic/ITemporaryRecordStore.hpp"
#include "simplnx/Filter/IFilter.hpp"
#include "simplnx/Utilities/AlgorithmDispatch.hpp"
#include "simplnx/Utilities/BoundedRecordPageCache.hpp"
#include "simplnx/Utilities/DataStoreUtilities.hpp"
#include "simplnx/Utilities/InMemoryTemporaryRecordStore.hpp"
#include "simplnx/Utilities/Meshing/TriangleUtilities.hpp"
#include "simplnx/Utilities/ParallelDataAlgorithm.hpp"

#include <fmt/format.h>
#include <nonstd/span.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <span>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace nx::core;

// The M3C core closely follows legacy DREAM3D M3CEntireVolume. Its flat,
// 1-based arrays preserve legacy topology and output ordering.
namespace
{
/**
 * @class TemporaryRecordVector
 * @brief Fixed-record scratch vector with a bounded typed page cache.
 * @tparam T Specifies the fixed scratch-record type.
 *
 * Candidate nodes, cube offsets, and triangle metadata use this wrapper. A
 * genuine OOC dispatch cannot fall back to resident scratch.
 */
template <typename T>
class TemporaryRecordVector
{
public:
  /**
   * @brief Creates the record store and its bounded typed cache.
   * @param recordCount Specifies initial fixed-record count.
   * @param requireExternalStore Prevents resident fallback during genuine OOC.
   * @param shouldCancel Stops before storage allocation when true.
   * @param recordsPerPage Specifies records per backing-store request.
   * @param cachePages Specifies maximum resident pages.
   * @return Initialized record vector, or a provider/allocation error.
   */
  static Result<TemporaryRecordVector> create(const uint64 recordCount, const bool requireExternalStore, const std::atomic_bool& shouldCancel, const uint64 recordsPerPage = 4096,
                                              const usize cachePages = 8)
  {
    if(recordCount > 0 && shouldCancel)
    {
      return MakeErrorResult<TemporaryRecordVector>(-90540, "M3C temporary-record vector creation was canceled.");
    }
    if(recordsPerPage == 0 || cachePages == 0 || recordsPerPage > std::numeric_limits<uint64>::max() / sizeof(T))
    {
      return MakeErrorResult<TemporaryRecordVector>(-90541, "M3C temporary-record vector configuration is invalid.");
    }

    TemporaryRecordStoreConfig config;
    config.recordSize = sizeof(T);
    config.maxRecordsPerBatch = recordsPerPage;
    config.initialRecordCount = recordCount;
    auto storeResult = DataStoreUtilities::GetIOCollection().createTemporaryRecordStore(config);

    std::unique_ptr<ITemporaryRecordStore> store;
    if(storeResult.valid())
    {
      store = std::move(storeResult.value());
    }
    else if(!requireExternalStore)
    {
      auto fallbackResult = InMemoryTemporaryRecordStore::Create(config);
      if(fallbackResult.invalid())
      {
        return ConvertInvalidResult<TemporaryRecordVector>(std::move(fallbackResult));
      }
      store = std::move(fallbackResult.value());
    }
    else
    {
      return ConvertInvalidResult<TemporaryRecordVector>(std::move(storeResult));
    }
    if(store == nullptr)
    {
      return MakeErrorResult<TemporaryRecordVector>(-90542, "M3C temporary-record provider returned a null store.");
    }

    TemporaryRecordVector vector;
    vector.m_Store = std::move(store);
    try
    {
      vector.m_Cache = std::make_unique<BoundedRecordPageCache<T>>(*vector.m_Store, recordsPerPage, cachePages);
    } catch(const std::bad_alloc&)
    {
      return MakeErrorResult<TemporaryRecordVector>(-90543, "M3C temporary-record vector could not allocate its bounded page cache.");
    }
    return {std::move(vector)};
  }

  TemporaryRecordVector() = default;
  TemporaryRecordVector(TemporaryRecordVector&&) noexcept = default;
  TemporaryRecordVector& operator=(TemporaryRecordVector&&) noexcept = default;
  TemporaryRecordVector(const TemporaryRecordVector&) = delete;
  TemporaryRecordVector& operator=(const TemporaryRecordVector&) = delete;

  /**
   * @brief Returns the owned byte-record store.
   * @return Store used for bulk operations.
   */
  ITemporaryRecordStore& store() noexcept
  {
    return *m_Store;
  }

  /**
   * @brief Returns the owned typed page cache.
   * @return Cache used for localized random access.
   */
  BoundedRecordPageCache<T>& cache() noexcept
  {
    return *m_Cache;
  }

  /**
   * @brief Writes dirty pages before the next algorithm phase.
   * @param shouldCancel Stops cache flushing when true.
   * @return Error from cache flushing, or success.
   */
  Result<> flush(const std::atomic_bool& shouldCancel)
  {
    return m_Cache->flush(shouldCancel);
  }

private:
  std::unique_ptr<ITemporaryRecordStore> m_Store{};
  std::unique_ptr<BoundedRecordPageCache<T>> m_Cache{};
};

/**
 * @struct M3CCandidateNodeRecord
 * @brief External scratch record for one possible M3C node.
 * @p Type records whether/how the candidate is used; @p PruneReferences records
 * whether dropped or surviving triangles reference it; @p CompactId is assigned
 * after counting all live candidates so output vertices can be written densely.
 */
struct M3CCandidateNodeRecord
{
  int8 Type = 0;
  uint8 PruneReferences = 0;
  std::array<std::byte, 6> Padding{};
  uint64 CompactId = 0;
};
static_assert(std::is_trivially_copyable_v<M3CCandidateNodeRecord>);

// A SiteIdType indexes the padded Feature Id grid. The 64-bit type prevents overflow
// when a large grid derives seven candidate-node IDs from each site.
using SiteIdType = int64;
// This sentinel marks candidate slots that are not real mesh nodes.
constexpr uint32 k_UnusedNodeId = std::numeric_limits<uint32>::max();

constexpr int k_NumNeigh = 26;

// --- M3C working structs (mirror SIMPL/Geometry/MeshStructs.h SurfaceMesh::M3C) ---
struct Node
{
  std::array<float, 3> Coord{};
};
struct VoxelCoord
{
  std::array<float, 3> Coord{};
};
/**
 * @struct Neighbor
 * @brief Stores one-based indexes for 26 neighboring sites.
 */
struct Neighbor
{
  std::array<SiteIdType, 27> NeighId{}; // 1-based; index 0 unused. 64-bit: these index the FeatureId grid.
};
/**
 * @struct Face
 * @brief Stores one marching square's edges and center node.
 */
struct Face
{
  // Recompute corner sites to keep this largest working array compact. Edge IDs
  // are 32-bit mesh indexes. Face-center node IDs retain 64-bit site indexes.
  std::array<uint32, 4> EdgeId{};
  SiteIdType FaceCenterNode = 0; // face-center node id, -1 if none
  int8 NEdge = 0;
  int8 Effect = 0; // 0 = useless square, 1 = straddles >=2 labels
};
/**
 * @struct Segment
 * @brief Stores one oriented face-edge segment and its labels.
 */
struct Segment
{
  std::array<int64, 2> NodeId{};
  std::array<int, 2> NSpin{}; // labels on left/right of the arrow
};
/**
 * @struct Triangle
 * @brief Stores one generated triangle and its adjacent labels.
 */
struct Triangle
{
  std::array<int64, 3> NodeId{};
  std::array<int, 2> NSpin{};
};

// Coordinates are pure functions of the padded site index. Compute them on
// demand to avoid full-volume coordinate arrays.
/**
 * @struct SiteCoords
 * @brief Calculates padded-grid coordinates on demand.
 */
struct SiteCoords
{
  usize FileDim0 = 0;
  usize FileDim1 = 0;
  usize FileNsp = 0; // FileDim0 * FileDim1
  std::array<float, 3> Res{};
  std::array<float, 3> Origin{};

  /**
   * @brief Calculates one site coordinate.
   * @param site Specifies a one-based padded site index.
   * @return Coordinate in image units.
   */
  VoxelCoord operator[](const int64 site) const
  {
    const usize linear = static_cast<usize>(site - 1);
    // Subtract the ghost shell so padded site (1,1,1) maps to the image origin.
    const int64 xIndex = static_cast<int64>(linear % FileDim0) - 1;
    const int64 yIndex = static_cast<int64>((linear / FileDim0) % FileDim1) - 1;
    const int64 zIndex = static_cast<int64>(linear / FileNsp) - 1;
    // A site is a CELL CENTER, not the cell's lower corner. initialize_nodes (legacy) places the 7
    // candidate nodes of a site at +half-spacing offsets. The marching cube spans from one site
    // to its (+1,+1,+1) neighbor - so the interface between two adjacent cells falls on the plane
    // midway between their centers, which is exactly their shared face. Returning the lower corner
    // here instead shifted every vertex by half a cell, placing the mesh partly outside the volume.
    return VoxelCoord{
        {((static_cast<float>(xIndex) + 0.5f) * Res[0]) + Origin[0], ((static_cast<float>(yIndex) + 0.5f) * Res[1]) + Origin[1], ((static_cast<float>(zIndex) + 0.5f) * Res[2]) + Origin[2]}};
  }
};

/**
 * @struct NodeCoords
 * @brief Calculates seven M3C candidate-node coordinates per site.
 */
struct NodeCoords
{
  SiteCoords Sites{};

  /**
   * @brief Calculates one candidate-node coordinate.
   * @param nodeId Specifies a zero-based candidate-node index.
   * @return Edge-midpoint, face-center, or body-center coordinate.
   *
   * Each site has three positive-edge midpoints, three positive-face centers,
   * and one body center. Their order matches the legacy node layout.
   */
  Node operator[](const int64 nodeId) const
  {
    const int64 site = (nodeId / 7) + 1;
    const int kind = static_cast<int>(nodeId % 7);
    const VoxelCoord siteCoord = Sites[site];
    const float halfX = Sites.Res[0] / 2.0f;
    const float halfY = Sites.Res[1] / 2.0f;
    const float halfZ = Sites.Res[2] / 2.0f;
    Node node{{siteCoord.Coord[0], siteCoord.Coord[1], siteCoord.Coord[2]}};
    switch(kind)
    {
    case 0:
      node.Coord[0] += halfX;
      break;
    case 1:
      node.Coord[1] += halfY;
      break;
    case 2:
      node.Coord[2] += halfZ;
      break;
    case 3:
      node.Coord[0] += halfX;
      node.Coord[1] += halfY;
      break;
    case 4:
      node.Coord[0] += halfX;
      node.Coord[2] += halfZ;
      break;
    case 5:
      node.Coord[1] += halfY;
      node.Coord[2] += halfZ;
      break;
    case 6:
      node.Coord[0] += halfX;
      node.Coord[1] += halfY;
      node.Coord[2] += halfZ;
      break;
    default:
      break;
    }
    return node;
  }
};

// These 20 multi-material marching-square cases match the legacy slice and
// whole-volume algorithms. Each case contains at most four edges.
// k_EdgeTable2d maps node-slot pairs to edges. Slots 0 through 3 are edge
// midpoints, and slot 4 is the face center.
// k_NsTable2d maps each edge to the two corner labels on its sides.
// clang-format off
constexpr std::array<std::array<int, 8>, 20> k_EdgeTable2d = {{
    {-1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1}, {0, 1, -1, -1, -1, -1, -1, -1},   {-1, -1, -1, -1, -1, -1, -1, -1},
    {0, 2, -1, -1, -1, -1, -1, -1},   {1, 2, -1, -1, -1, -1, -1, -1},   {0, 4, 2, 4, 1, 4, -1, -1},       {-1, -1, -1, -1, -1, -1, -1, -1}, {3, 0, -1, -1, -1, -1, -1, -1},
    {3, 1, -1, -1, -1, -1, -1, -1},   {3, 4, 0, 4, 1, 4, -1, -1},       {2, 3, -1, -1, -1, -1, -1, -1},   {3, 4, 0, 4, 2, 4, -1, -1},       {3, 4, 1, 4, 2, 4, -1, -1},
    {3, 0, 1, 2, -1, -1, -1, -1},     {0, 1, 2, 3, -1, -1, -1, -1},     {0, 1, 2, 3, -1, -1, -1, -1},     {3, 0, 1, 2, -1, -1, -1, -1},     {3, 4, 1, 4, 0, 4, 2, 4}}};

constexpr std::array<std::array<int, 8>, 20> k_NsTable2d = {{
    {-1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1}, {1, 0, -1, -1, -1, -1, -1, -1},   {-1, -1, -1, -1, -1, -1, -1, -1},
    {1, 0, -1, -1, -1, -1, -1, -1},   {2, 1, -1, -1, -1, -1, -1, -1},   {1, 0, 3, 2, 2, 1, -1, -1},       {-1, -1, -1, -1, -1, -1, -1, -1}, {0, 3, -1, -1, -1, -1, -1, -1},
    {0, 3, -1, -1, -1, -1, -1, -1},   {0, 3, 1, 0, 2, 1, -1, -1},       {3, 2, -1, -1, -1, -1, -1, -1},   {0, 3, 1, 0, 3, 2, -1, -1},       {0, 3, 2, 1, 3, 2, -1, -1},
    {0, 3, 2, 1, -1, -1, -1, -1},     {1, 0, 3, 2, -1, -1, -1, -1},     {1, 0, 3, 2, -1, -1, -1, -1},     {0, 3, 2, 1, -1, -1, -1, -1},     {0, 3, 2, 1, 1, 0, 3, 2}}};
// clang-format on

// -----------------------------------------------------------------------------
// Copy FeatureIds into a 1-based working grid, wrapping it in a ghost shell of
// negative labels when addSurfaceLayer is true, fill voxel coordinates,
// and renumber any FeatureId==0 to maxGrainId. Returns maxGrainId (the value that
// zeros were remapped to; callers revert it on output). Transcribed from
// M3CEntireVolume::initialize_micro_from_grainIds.
// -----------------------------------------------------------------------------
// The single sentinel used for every cell of the ghost shell. Any negative value works; only the
// sign is tested. It must be the SAME for all ghost cells - see initializeMicro.
constexpr int32 k_GhostLabel = -3;

int initializeMicro(const bool addSurfaceLayer, const std::array<usize, 3>& dims, const std::array<usize, 3>& fileDim, const AbstractDataStore<int32>& grainIds, std::vector<int32>& featureIds)
{
  int maxGrainId = 0;

  if(!addSurfaceLayer)
  {
    const usize totalPoints = dims[0] * dims[1] * dims[2];
    for(usize i = 0; i < totalPoints; ++i)
    {
      featureIds[i + 1] = grainIds[i];
      if(featureIds[i + 1] > maxGrainId)
      {
        maxGrainId = featureIds[i + 1];
      }
    }
  }
  else
  {
    // Wrap the volume in a one-cell ghost shell. Ghost cells carry a NEGATIVE sentinel label so the
    // marching-cubes code treats them as "outside the volume"; only the sign is ever tested.
    //
    // Legacy used six DISTINCT sentinels here (-3 bottom z-slice, -4/-7 the y-row pads, -5/-6 the
    // per-row x-end pads, -8 top z-slice) to record which face or edge of the shell a ghost cell
    // belonged to. Nothing reads that back, but the marching cubes compares labels for INEQUALITY,
    // so neighboring ghost cells with different sentinels looked like a material interface and were
    // triangulated - generating surface outside the volume along the shell's own internal seams.
    // A single shared sentinel leaves the shell internally uniform, so the only interfaces it can
    // produce are the real ghost-to-feature ones that form the volume's exterior surface.
    usize index = 0;
    usize gIdx = 0;

    // Bottom wrapping slice
    for(usize i = 0; i < (fileDim[0] * fileDim[1]); ++i)
    {
      featureIds[++index] = k_GhostLabel;
    }
    // Bulk of the volume, wrapped per-plane and per-row
    for(usize zIndex = 0; zIndex < dims[2]; ++zIndex)
    {
      for(usize i = 0; i < fileDim[0]; ++i)
      {
        featureIds[++index] = k_GhostLabel;
      }
      for(usize yIndex = 0; yIndex < dims[1]; ++yIndex)
      {
        featureIds[++index] = k_GhostLabel; // leading surface voxel for this row
        for(usize xIndex = 0; xIndex < dims[0]; ++xIndex)
        {
          featureIds[++index] = grainIds[gIdx++];
          if(featureIds[index] > maxGrainId)
          {
            maxGrainId = featureIds[index];
          }
        }
        featureIds[++index] = k_GhostLabel; // trailing surface voxel for this row
      }
      for(usize i = 0; i < fileDim[0]; ++i)
      {
        featureIds[++index] = k_GhostLabel;
      }
    }
    // Top wrapping slice
    for(usize i = 0; i < (fileDim[0] * fileDim[1]); ++i)
    {
      featureIds[++index] = k_GhostLabel;
    }
  }

  // Reserve one positive label for input Feature Id 0.
  maxGrainId = maxGrainId + 1;

  featureIds[0] = 0; // Point 0 is garbage

  // Renumber zero labels without changing negative ghost cells. Coordinates are
  // computed on demand by SiteCoords and NodeCoords.
  const usize totalPoints = fileDim[0] * fileDim[1] * fileDim[2];
  for(usize id = 1; id <= totalPoints; id++)
  {
    if(featureIds[id] == 0)
    {
      featureIds[id] = maxGrainId;
    }
  }
  return maxGrainId;
}

/**
 * @struct NeighborAccessor
 * @brief Reconstructs 26 legacy neighbors for one padded site.
 *
 * The ghost shell makes toroidal border indexes harmless. Cache the returned
 * Neighbor when a caller needs multiple neighbor indexes.
 */
struct NeighborAccessor
{
  SiteIdType Ns = 0;
  SiteIdType Nsp = 0;
  int XDim = 0;

  /**
   * @brief Calculates the 26 neighbors of one padded site.
   * @param site_id Specifies a one-based padded site index.
   * @return Neighbor indexes in the legacy order.
   */
  Neighbor operator[](const SiteIdType site_id) const
  {
    // Recover the legacy loop coordinates for this one-based site.
    const SiteIdType within = (site_id - 1) % Nsp;          // == rowOffset + (xIndex - 1)
    const int xIndex = static_cast<int>(within % XDim) + 1; // 1..XDim
    const SiteIdType rowOffset = within - (xIndex - 1);     // multiple of XDim, 0..Nsp-XDim
    const SiteIdType planeOffset = ((site_id - 1) / Nsp) * Nsp;

    Neighbor siteNeighbors{};
    siteNeighbors.NeighId[0] = 0; // index 0 unused

    // same plane
    siteNeighbors.NeighId[1] = planeOffset + rowOffset + (xIndex % XDim) + 1;
    siteNeighbors.NeighId[2] = planeOffset + ((rowOffset - XDim + Nsp) % Nsp) + (xIndex % XDim) + 1;
    siteNeighbors.NeighId[3] = planeOffset + ((rowOffset - XDim + Nsp) % Nsp) + xIndex;
    siteNeighbors.NeighId[4] = planeOffset + ((rowOffset - XDim + Nsp) % Nsp) + ((xIndex - 2 + XDim) % XDim) + 1;
    siteNeighbors.NeighId[5] = planeOffset + rowOffset + ((xIndex - 2 + XDim) % XDim) + 1;
    siteNeighbors.NeighId[6] = planeOffset + ((rowOffset + XDim) % Nsp) + ((xIndex - 2 + XDim) % XDim) + 1;
    siteNeighbors.NeighId[7] = planeOffset + ((rowOffset + XDim) % Nsp) + xIndex;
    siteNeighbors.NeighId[8] = planeOffset + ((rowOffset + XDim) % Nsp) + (xIndex % XDim) + 1;

    // upper plane
    siteNeighbors.NeighId[9] = ((planeOffset - Nsp + Ns) % Ns) + rowOffset + xIndex;
    siteNeighbors.NeighId[10] = ((planeOffset - Nsp + Ns) % Ns) + rowOffset + (xIndex % XDim) + 1;
    siteNeighbors.NeighId[11] = ((planeOffset - Nsp + Ns) % Ns) + ((rowOffset - XDim + Nsp) % Nsp) + (xIndex % XDim) + 1;
    siteNeighbors.NeighId[12] = ((planeOffset - Nsp + Ns) % Ns) + ((rowOffset - XDim + Nsp) % Nsp) + xIndex;
    siteNeighbors.NeighId[13] = ((planeOffset - Nsp + Ns) % Ns) + ((rowOffset - XDim + Nsp) % Nsp) + ((xIndex - 2 + XDim) % XDim) + 1;
    siteNeighbors.NeighId[14] = ((planeOffset - Nsp + Ns) % Ns) + rowOffset + ((xIndex - 2 + XDim) % XDim) + 1;
    siteNeighbors.NeighId[15] = ((planeOffset - Nsp + Ns) % Ns) + ((rowOffset + XDim) % Nsp) + ((xIndex - 2 + XDim) % XDim) + 1;
    siteNeighbors.NeighId[16] = ((planeOffset - Nsp + Ns) % Ns) + ((rowOffset + XDim) % Nsp) + xIndex;
    siteNeighbors.NeighId[17] = ((planeOffset - Nsp + Ns) % Ns) + ((rowOffset + XDim) % Nsp) + (xIndex % XDim) + 1;

    // lower plane
    siteNeighbors.NeighId[18] = ((planeOffset + Nsp) % Ns) + rowOffset + xIndex;
    siteNeighbors.NeighId[19] = ((planeOffset + Nsp) % Ns) + rowOffset + (xIndex % XDim) + 1;
    siteNeighbors.NeighId[20] = ((planeOffset + Nsp) % Ns) + ((rowOffset - XDim + Nsp) % Nsp) + (xIndex % XDim) + 1;
    siteNeighbors.NeighId[21] = ((planeOffset + Nsp) % Ns) + ((rowOffset - XDim + Nsp) % Nsp) + xIndex;
    siteNeighbors.NeighId[22] = ((planeOffset + Nsp) % Ns) + ((rowOffset - XDim + Nsp) % Nsp) + ((xIndex - 2 + XDim) % XDim) + 1;
    siteNeighbors.NeighId[23] = ((planeOffset + Nsp) % Ns) + rowOffset + ((xIndex - 2 + XDim) % XDim) + 1;
    siteNeighbors.NeighId[24] = ((planeOffset + Nsp) % Ns) + ((rowOffset + XDim) % Nsp) + ((xIndex - 2 + XDim) % XDim) + 1;
    siteNeighbors.NeighId[25] = ((planeOffset + Nsp) % Ns) + ((rowOffset + XDim) % Nsp) + xIndex;
    siteNeighbors.NeighId[26] = ((planeOffset + Nsp) % Ns) + ((rowOffset + XDim) % Nsp) + (xIndex % XDim) + 1;
    return siteNeighbors;
  }
};

/**
 * @brief Calculates the four corner sites of a marching square.
 * @param squareId Encodes the source site and square orientation.
 * @param neighbors Provides padded-grid neighbor indexes.
 * @return Corner site indexes in marching-square order.
 *
 * On-demand calculation keeps four site indexes out of every Face record.
 */
std::array<SiteIdType, 4> squareCorners(const SiteIdType squareId, const NeighborAccessor& neighbors)
{
  const SiteIdType site = (squareId / 3) + 1;
  const int ord = static_cast<int>(squareId % 3);
  const Neighbor siteNeighbors = neighbors[site];
  switch(ord)
  {
  case 0: // top (same z)
    return {site, siteNeighbors.NeighId[1], siteNeighbors.NeighId[8], siteNeighbors.NeighId[7]};
  case 1: // back (same y)
    return {site, siteNeighbors.NeighId[1], siteNeighbors.NeighId[19], siteNeighbors.NeighId[18]};
  default: // left (same x)
    return {siteNeighbors.NeighId[7], site, siteNeighbors.NeighId[18], siteNeighbors.NeighId[25]};
  }
}

/**
 * @brief Initializes three empty marching squares per padded site.
 * @param squares Receives empty edge and flag fields.
 * @param numSitesDim3 Specifies padded site count.
 *
 * Candidate coordinates are calculated on demand. The node-type vector uses
 * value initialization, so neither data set needs a separate initialization pass.
 */
void initializeSquares(std::vector<Face>& squares, const SiteIdType numSitesDim3)
{
  for(SiteIdType sqId = 0; sqId < 3 * numSitesDim3; sqId++)
  {
    for(int j = 0; j < 4; j++)
    {
      squares[sqId].EdgeId[j] = k_UnusedNodeId;
    }
    squares[sqId].NEdge = 0;
    squares[sqId].FaceCenterNode = -1;
    squares[sqId].Effect = 0;
  }
}

/**
 * @namespace m3c_node_type
 * @brief Defines node categories consumed by mesh-smoothing algorithms.
 *
 * These values match the legacy SurfaceMesh NodeType contract.
 */
namespace m3c_node_type
{
constexpr int8 k_Unused = 0;
constexpr int8 k_Default = 2;
constexpr int8 k_TriplePoint = 3;
constexpr int8 k_QuadPoint = 4;
constexpr int8 k_SurfaceDefault = 12;
constexpr int8 k_SurfaceTriplePoint = 13;
constexpr int8 k_SurfaceQuadPoint = 14;
} // namespace m3c_node_type

/**
 * @brief Classifies four corner labels into a marching-square case.
 * @param tns Provides four corner labels in square order.
 * @return Case index from 0 through 19.
 */
int getSquareIndex(const std::array<int, 4>& tns)
{
  std::array<int, 6> aBit{};
  aBit[0] = (tns[0] == tns[1]) ? 0 : 1;
  aBit[1] = (tns[1] == tns[2]) ? 0 : 1;
  aBit[2] = (tns[2] == tns[3]) ? 0 : 1;
  aBit[3] = (tns[3] == tns[0]) ? 0 : 1;
  aBit[4] = (tns[0] == tns[2]) ? 0 : 1;
  aBit[5] = (tns[1] == tns[3]) ? 0 : 1;

  int tempIndex = (8 * aBit[3]) + (4 * aBit[2]) + (2 * aBit[1]) + (1 * aBit[0]);
  if(tempIndex == 15)
  {
    const int subIndex = (2 * aBit[4]) + (1 * aBit[5]);
    if(subIndex != 0)
    {
      tempIndex = tempIndex + subIndex + 1;
    }
  }
  return tempIndex;
}

/**
 * @brief Resolves the case-15 saddle from eight in-plane neighbors.
 * @param tnst Provides four corner site indexes.
 * @param featureIds Provides padded Feature Id values.
 * @param neighbors Calculates padded-grid neighbors.
 * @param sqid Is unused by the legacy-compatible calculation.
 * @return Zero or one to select the case-15 topology.
 *
 * The algorithm matches M3CSliceBySlice by connecting the corner with the fewest
 * in-plane same-label neighbors. The 26-neighbor volume variant can create ties
 * that are resolved arbitrarily and produce spurious handles.
 */
int treatAnomaly(const std::array<SiteIdType, 4>& tnst, const std::vector<int32>& featureIds, const NeighborAccessor& neighbors, SiteIdType /*sqid*/)
{
  std::array<int, 4> numNeigh = {0, 0, 0, 0};

  for(int i = 0; i < 4; i++)
  {
    const SiteIdType csite = tnst[i];
    const int cspin = featureIds[csite];
    const Neighbor siteNeighbors = neighbors[csite]; // Cache the 8 in-plane neighbors read below.
    for(int j = 1; j <= 8; j++)
    {
      const SiteIdType nsite = siteNeighbors.NeighId[j];
      const int nspin = featureIds[nsite];
      if(cspin == nspin)
      {
        numNeigh[i] = numNeigh[i] + 1;
      }
    }
  }

  int min = 1000;
  int minid = -1;
  for(int ii = 0; ii < 4; ii++)
  {
    if(numNeigh[ii] < min)
    {
      min = numNeigh[ii];
      minid = ii;
    }
  }

  int tempFlag = 0;
  if(minid == -1 || minid == 1 || minid == 3)
  {
    tempFlag = 0;
  }
  else
  {
    tempFlag = 1;
  }
  return tempFlag;
}

/**
 * @brief Maps edge-table slots to candidate-node indexes.
 * @param cst Specifies the square origin site.
 * @param ord Specifies the square orientation.
 * @param nidx Provides two edge-table node slots.
 * @param nid Receives two candidate-node indexes.
 * @param numSitesDim2 Specifies padded sites per Z plane.
 * @param xDim1 Specifies padded X dimension.
 */
void getNodes(const SiteIdType cst, const int ord, const std::array<int, 2>& nidx, std::array<SiteIdType, 2>& nid, const SiteIdType numSitesDim2, const int xDim1)
{
  for(int ii = 0; ii < 2; ii++)
  {
    const int tempIndex = nidx[ii];
    if(ord == 0)
    {
      switch(tempIndex)
      {
      case 0:
        nid[ii] = 7 * (cst - 1);
        break;
      case 1:
        nid[ii] = (7 * cst) + 1;
        break;
      case 2:
        nid[ii] = 7 * (cst + xDim1 - 1);
        break;
      case 3:
        nid[ii] = (7 * (cst - 1)) + 1;
        break;
      case 4:
        nid[ii] = (7 * (cst - 1)) + 3;
        break;
      }
    }
    else if(ord == 1)
    {
      switch(tempIndex)
      {
      case 0:
        nid[ii] = 7 * (cst - 1);
        break;
      case 1:
        nid[ii] = (7 * cst) + 2;
        break;
      case 2:
        nid[ii] = 7 * (cst + numSitesDim2 - 1);
        break;
      case 3:
        nid[ii] = (7 * (cst - 1)) + 2;
        break;
      case 4:
        nid[ii] = (7 * (cst - 1)) + 4;
        break;
      }
    }
    else
    {
      switch(tempIndex)
      {
      case 0:
        nid[ii] = (7 * (cst - 1)) + 1;
        break;
      case 1:
        nid[ii] = (7 * (cst - 1)) + 2;
        break;
      case 2:
        nid[ii] = (7 * (cst + numSitesDim2 - 1)) + 1;
        break;
      case 3:
        nid[ii] = (7 * (cst + xDim1 - 1)) + 2;
        break;
      case 4:
        nid[ii] = (7 * (cst - 1)) + 5;
        break;
      }
    }
  }
}

/**
 * @brief Maps square-corner slots to two edge-side labels.
 * @param featureIds Provides padded Feature Id values.
 * @param cst Specifies the square origin site.
 * @param ord Specifies the square orientation.
 * @param pID Provides two square-corner slots.
 * @param pSpin Receives the two Feature Id values.
 * @param numSitesDim2 Specifies padded sites per Z plane.
 * @param xDim1 Specifies padded X dimension.
 */
void getSpins(const std::vector<int32>& featureIds, const SiteIdType cst, const int ord, const std::array<int, 2>& pID, std::array<int, 2>& pSpin, const SiteIdType numSitesDim2, const int xDim1)
{
  for(int i = 0; i < 2; i++)
  {
    const int pixTemp = pID[i];
    if(ord == 0)
    {
      switch(pixTemp)
      {
      case 0:
        pSpin[i] = featureIds[cst];
        break;
      case 1:
        pSpin[i] = featureIds[cst + 1];
        break;
      case 2:
        pSpin[i] = featureIds[cst + xDim1 + 1];
        break;
      case 3:
        pSpin[i] = featureIds[cst + xDim1];
        break;
      }
    }
    else if(ord == 1)
    {
      switch(pixTemp)
      {
      case 0:
        pSpin[i] = featureIds[cst];
        break;
      case 1:
        pSpin[i] = featureIds[cst + 1];
        break;
      case 2:
        pSpin[i] = featureIds[cst + numSitesDim2 + 1];
        break;
      case 3:
        pSpin[i] = featureIds[cst + numSitesDim2];
        break;
      }
    }
    else if(ord == 2)
    {
      switch(pixTemp)
      {
      case 0:
        pSpin[i] = featureIds[cst + xDim1];
        break;
      case 1:
        pSpin[i] = featureIds[cst];
        break;
      case 2:
        pSpin[i] = featureIds[cst + numSitesDim2];
        break;
      case 3:
        pSpin[i] = featureIds[cst + numSitesDim2 + xDim1];
        break;
      }
    }
  }
}

/**
 * @brief Counts face edges and marks effective squares.
 * @param squares Receives each square's effect flag.
 * @param featureIds Provides padded Feature Id values.
 * @param neighbors Calculates padded-grid neighbors.
 * @param numSitesDim3 Specifies padded site count.
 * @param shouldCancel Stops before later squares when true.
 * @return Count accumulated before completion or cancellation.
 *
 * The count permits one exact allocation before edge generation.
 */
int64 getNumberFEdges(std::vector<Face>& squares, const std::vector<int32>& featureIds, const NeighborAccessor& neighbors, const SiteIdType numSitesDim3, const std::atomic_bool& shouldCancel)
{
  int64 sumEdge = 0;
  for(SiteIdType k = 0; k < (3 * numSitesDim3); k++)
  {
    if(shouldCancel)
    {
      return sumEdge;
    }
    const std::array<SiteIdType, 4> tnsite = squareCorners(k, neighbors);
    std::array<int, 4> tnspin{};
    int numGhostCorners = 0;
    for(int cornerIdx = 0; cornerIdx < 4; cornerIdx++)
    {
      tnspin[cornerIdx] = featureIds[tnsite[cornerIdx]];
      if(tnspin[cornerIdx] < 0)
      {
        numGhostCorners++;
      }
    }
    if(numGhostCorners != 4)
    {
      squares[k].Effect = 1; // mark as effective (can be marching-cubed)
    }

    if(numGhostCorners != 4)
    {
      int sqIndex = getSquareIndex(tnspin);
      if(sqIndex == 15)
      {
        sqIndex = sqIndex + treatAnomaly(tnsite, featureIds, neighbors, k);
      }

      int numCEdge = 0;
      if(sqIndex == 0)
      {
        numCEdge = 0;
      }
      else if(sqIndex == 19)
      {
        numCEdge = 4;
      }
      else if(sqIndex == 15 || sqIndex == 16 || sqIndex == 17 || sqIndex == 18)
      {
        numCEdge = 2;
      }
      else if(sqIndex == 7 || sqIndex == 11 || sqIndex == 13 || sqIndex == 14)
      {
        if(numGhostCorners == 3)
        {
          numCEdge = 2;
        }
        else if(numGhostCorners == 1)
        {
          // A single negative corner is not a valid legacy square case.
          numCEdge = 0;
        }
        else
        {
          numCEdge = 3;
        }
      }
      else
      {
        numCEdge = 1;
      }
      sumEdge = sumEdge + numCEdge;
    }
  }
  return sumEdge;
}

/**
 * @brief Creates face edges and classifies their candidate nodes.
 * @param squares Receives edge indexes and face-center nodes.
 * @param featureIds Provides padded Feature Id values.
 * @param neighbors Calculates padded-grid neighbors.
 * @param nodeType Receives candidate-node categories.
 * @param faceEdges Receives face-edge records.
 * @param numSitesDim3 Specifies padded site count.
 * @param numSitesDim2 Specifies padded sites per Z plane.
 * @param xDim Specifies padded X dimension.
 * @param shouldCancel Stops before later squares when true.
 */
void getNodesFEdges(std::vector<Face>& squares, const std::vector<int32>& featureIds, const NeighborAccessor& neighbors, std::vector<int8>& nodeType, std::vector<Segment>& faceEdges,
                    const SiteIdType numSitesDim3, const SiteIdType numSitesDim2, const int xDim, const std::atomic_bool& shouldCancel)
{
  int64 eid = 0;
  for(SiteIdType k = 0; k < (3 * numSitesDim3); k++)
  {
    if(shouldCancel)
    {
      return;
    }
    const SiteIdType cubeOrigin = (k / 3) + 1;
    const int sqOrder = static_cast<int>(k % 3);

    const std::array<SiteIdType, 4> tnsite = squareCorners(k, neighbors);
    std::array<int, 4> tnspin{};
    int numGhostCorners = 0;
    for(int cornerIdx = 0; cornerIdx < 4; cornerIdx++)
    {
      tnspin[cornerIdx] = featureIds[tnsite[cornerIdx]];
      if(tnspin[cornerIdx] < 0)
      {
        numGhostCorners++;
      }
    }

    int edgeCount = 0;
    if(numGhostCorners != 4)
    {
      int sqIndex = getSquareIndex(tnspin);
      if(sqIndex == 15)
      {
        sqIndex = sqIndex + treatAnomaly(tnsite, featureIds, neighbors, k);
      }
      if(sqIndex != 0)
      {
        for(int j = 0; j < 8; j = j + 2)
        {
          if(k_EdgeTable2d[sqIndex][j] != -1)
          {
            std::array<int, 2> nodeIndex = {k_EdgeTable2d[sqIndex][j], k_EdgeTable2d[sqIndex][j + 1]};
            const std::array<int, 2> pixIndex = {k_NsTable2d[sqIndex][j], k_NsTable2d[sqIndex][j + 1]};
            std::array<SiteIdType, 2> nodeID{};
            std::array<int, 2> pixSpin{};
            getNodes(cubeOrigin, sqOrder, nodeIndex, nodeID, numSitesDim2, xDim);
            getSpins(featureIds, cubeOrigin, sqOrder, pixIndex, pixSpin, numSitesDim2, xDim);

            if(pixSpin[0] > 0 || pixSpin[1] > 0)
            {
              faceEdges[eid].NodeId[0] = nodeID[0];
              faceEdges[eid].NodeId[1] = nodeID[1];
              faceEdges[eid].NSpin[0] = pixSpin[0];
              faceEdges[eid].NSpin[1] = pixSpin[1];
              squares[k].EdgeId[edgeCount] = static_cast<uint32>(eid);
              edgeCount++;
              eid++;
            }
            else
            {
              // Pure exterior edges do not create output mesh nodes.
              nodeType[nodeID[0]] = m3c_node_type::k_Unused;
              nodeType[nodeID[1]] = m3c_node_type::k_Unused;
            }

            // Face centers represent triple or quad points. Other slots represent
            // default interface nodes.
            for(int ii = 0; ii < 2; ii++)
            {
              if(nodeIndex[ii] == 4)
              {
                if(sqIndex == 7 || sqIndex == 11 || sqIndex == 13 || sqIndex == 14)
                {
                  const SiteIdType tnode = nodeID[ii];
                  squares[k].FaceCenterNode = tnode;
                  nodeType[tnode] = m3c_node_type::k_TriplePoint;
                }
                else if(sqIndex == 19)
                {
                  const SiteIdType tnode = nodeID[ii];
                  squares[k].FaceCenterNode = tnode;
                  nodeType[tnode] = m3c_node_type::k_QuadPoint;
                }
              }
              else
              {
                // Every interior edge endpoint is a real mesh node. Without this
                // promotion, compaction can remove a node that stored edges reference.
                const SiteIdType tnode = nodeID[ii];
                nodeType[tnode] = m3c_node_type::k_Default;
              }
            }
          }
        }
      }
    }
    // Each square has at most four edges.
    squares[k].NEdge = static_cast<int8>(edgeCount);
  }
}

/**
 * @struct FaceEdgeLoops
 * @brief Holds the closed face-edge loops of one cube.
 * @details Each case handler walks the face edges of a cube loop by loop to
 * make triangles. A handler needs the edges in loop order and the size of
 * each loop. This struct keeps those results together, so one helper
 * computes them once for every handler.
 */
struct FaceEdgeLoops
{
  /** Lists the face edges in loop order. A handler reads this list to follow one loop. */
  std::vector<SiteIdType> BurntList{};
  /** Gives the edge count of each loop. The index is the loop number, so element 0 stays unused. */
  std::vector<int> Count{};
  /** Gives one more than the highest loop number. A handler counts from 1 up to this value. */
  int LoopId = 0;
};

/**
 * @brief Groups the face edges of one cube into closed loops.
 * @details The helper starts at each unburnt edge. It then adds every edge that
 * has the same two labels and one common node. The helper does not change the
 * edge records. The handlers that call it start each triangle fan at a known
 * node, so they do not need a common edge direction.
 * @param afe Provides cube face-edge indexes.
 * @param faceEdges Provides face-edge records. The helper only reads them.
 * @param nfedge Specifies cube face-edge count.
 * @return Edges in loop order, loop sizes and loop count.
 */
FaceEdgeLoops burnFaceEdgeLoops(const std::span<const SiteIdType> afe, const std::span<const Segment> faceEdges, const int nfedge)
{
  std::vector<int> burnt(nfedge, 0);
  std::vector<SiteIdType> burntList(nfedge, -1);

  int loopID = 1;
  int tail = 0;
  int head = 0;

  for(int i = 0; i < nfedge; i++)
  {
    const SiteIdType cedge = afe[i];
    if(burnt[i] == 0)
    {
      burnt[i] = loopID;
      burntList[tail] = cedge;
      int coin = 0;
      do
      {
        const SiteIdType chaser = burntList[tail];
        const int cspin1 = faceEdges[chaser].NSpin[0];
        const int cspin2 = faceEdges[chaser].NSpin[1];
        const SiteIdType cnode1 = faceEdges[chaser].NodeId[0];
        const SiteIdType cnode2 = faceEdges[chaser].NodeId[1];

        for(int j = 0; j < nfedge; j++)
        {
          const SiteIdType nedge = afe[j];
          if(burnt[j] == 0)
          {
            const int nspin1 = faceEdges[nedge].NSpin[0];
            const int nspin2 = faceEdges[nedge].NSpin[1];
            const SiteIdType nnode1 = faceEdges[nedge].NodeId[0];
            const SiteIdType nnode2 = faceEdges[nedge].NodeId[1];
            const int spinFlag = (((cspin1 == nspin1) && (cspin2 == nspin2)) || ((cspin1 == nspin2) && (cspin2 == nspin1))) ? 1 : 0;
            int nodeFlag = 0;
            if(((cnode1 == nnode1) && (cnode2 != nnode2)) || ((cnode1 == nnode2) && (cnode2 != nnode1)) || ((cnode2 == nnode1) && (cnode1 != nnode2)) || ((cnode2 == nnode2) && (cnode1 != nnode1)))
            {
              nodeFlag = 1;
            }
            if(spinFlag == 1 && nodeFlag == 1)
            {
              head = head + 1;
              burntList[head] = nedge;
              burnt[j] = loopID;
            }
          }
        }

        if(tail == head)
        {
          coin = 0;
          tail = tail + 1;
          head = tail;
          loopID++;
        }
        else
        {
          tail = tail + 1;
          coin = 1;
        }
      } while(coin);
    }
  }

  std::vector<int> count(loopID, 0);
  for(int k = 1; k < loopID; k++)
  {
    for(int kk = 0; kk < nfedge; kk++)
    {
      if(k == burnt[kk])
      {
        count[k] = count[k] + 1;
      }
    }
  }

  return FaceEdgeLoops{std::move(burntList), std::move(count), loopID};
}

/**
 * @brief Groups the face edges of one cube into closed loops and turns them to one direction.
 * @details The helper matches the head node of the current edge against each
 * node of a candidate edge. A match on the tail node means the candidate edge
 * points the wrong way. The helper then writes through @p faceEdges and swaps the two
 * labels and the two nodes of that edge. Thus all edges of a loop point the
 * same way. The case 0 handlers need this direction, because they make a
 * triangle fan directly from the loop order.
 * @param afe Provides cube face-edge indexes.
 * @param faceEdges Provides face-edge records. The helper writes the turned records back.
 * @param nfedge Specifies cube face-edge count.
 * @return Edges in loop order, loop sizes and loop count.
 */
FaceEdgeLoops burnAndOrientFaceEdgeLoops(const std::span<const SiteIdType> afe, const std::span<Segment> faceEdges, const int nfedge)
{
  std::vector<int> burnt(nfedge, 0);
  std::vector<SiteIdType> burntList(nfedge, -1);

  int loopID = 1;
  int tail = 0;
  int head = 0;

  for(int i = 0; i < nfedge; i++)
  {
    const SiteIdType cedge = afe[i];
    if(burnt[i] == 0)
    {
      burnt[i] = loopID;
      burntList[tail] = cedge;
      int coin = 0;
      do
      {
        const SiteIdType chaser = burntList[tail];
        const int cspin1 = faceEdges[chaser].NSpin[0];
        const int cspin2 = faceEdges[chaser].NSpin[1];
        const SiteIdType cnode1 = faceEdges[chaser].NodeId[0];
        const SiteIdType cnode2 = faceEdges[chaser].NodeId[1];

        for(int j = 0; j < nfedge; j++)
        {
          const SiteIdType nedge = afe[j];
          if(burnt[j] == 0)
          {
            const int nspin1 = faceEdges[nedge].NSpin[0];
            const int nspin2 = faceEdges[nedge].NSpin[1];
            const SiteIdType nnode1 = faceEdges[nedge].NodeId[0];
            const SiteIdType nnode2 = faceEdges[nedge].NodeId[1];
            const int spinFlag = (((cspin1 == nspin1) && (cspin2 == nspin2)) || ((cspin1 == nspin2) && (cspin2 == nspin1))) ? 1 : 0;
            int nodeFlag = 0;
            int flip = 0;
            if((cnode2 == nnode1) && (cnode1 != nnode2))
            {
              nodeFlag = 1;
              flip = 0;
            }
            else if((cnode2 == nnode2) && (cnode1 != nnode1))
            {
              nodeFlag = 1;
              flip = 1;
            }
            else
            {
              nodeFlag = 0;
              flip = 0;
            }
            if(spinFlag == 1 && nodeFlag == 1)
            {
              head = head + 1;
              burntList[head] = nedge;
              burnt[j] = loopID;
              if(flip == 1)
              {
                faceEdges[nedge].NSpin[0] = nspin2;
                faceEdges[nedge].NSpin[1] = nspin1;
                faceEdges[nedge].NodeId[0] = nnode2;
                faceEdges[nedge].NodeId[1] = nnode1;
              }
            }
          }
        }

        if(tail == head)
        {
          coin = 0;
          tail = tail + 1;
          head = tail;
          loopID++;
        }
        else
        {
          tail = tail + 1;
          coin = 1;
        }
      } while(coin);
    }
  }

  std::vector<int> count(loopID, 0);
  for(int k = 1; k < loopID; k++)
  {
    for(int kk = 0; kk < nfedge; kk++)
    {
      if(k == burnt[kk])
      {
        count[k] = count[k] + 1;
      }
    }
  }

  return FaceEdgeLoops{std::move(burntList), std::move(count), loopID};
}

/**
 * @brief Counts triangles for a cube without face centers.
 * @param afe Provides cube face-edge indexes.
 * @param faceEdges Provides oriented face-edge records.
 * @param nfedge Specifies cube face-edge count.
 * @return Triangle count after closed-loop fan triangulation.
 */
int getNumberCase0Triangles(const std::span<const SiteIdType> afe, const std::span<Segment> faceEdges, const int nfedge)
{
  const FaceEdgeLoops loops = burnAndOrientFaceEdgeLoops(afe, faceEdges, nfedge);
  const std::vector<int>& count = loops.Count;
  const int loopID = loops.LoopId;

  int numTri = 0;
  for(int jj = 1; jj < loopID; jj++)
  {
    const int numN = count[jj];
    if(numN == 3)
    {
      numTri = numTri + 1;
    }
    else if(numN > 3)
    {
      numTri = numTri + (numN - 2);
    }
  }
  return numTri;
}

/**
 * @brief Counts triangles for a cube with two face centers.
 * @param afe Provides cube face-edge indexes.
 * @param faceEdges Provides mutable oriented face-edge records.
 * @param nfedge Specifies cube face-edge count.
 * @param afc Provides face-center node indexes.
 * @param nfctr Is fixed at two and is unused.
 * @return Triangle count after open- and closed-loop triangulation.
 *
 * Valid label data extends each chase loop by one edge. Loop guards bound
 * malformed or non-manifold input instead of overrunning a buffer.
 */
int getNumberCase2Triangles(const std::span<const SiteIdType> afe, const std::span<Segment> faceEdges, const int nfedge, const std::array<SiteIdType, 6>& afc, int /*nfctr*/)
{
  const FaceEdgeLoops loops = burnFaceEdgeLoops(afe, faceEdges, nfedge);
  const std::vector<SiteIdType>& burntList = loops.BurntList;
  const std::vector<int>& count = loops.Count;
  const int loopID = loops.LoopId;

  int numTri = 0;
  const SiteIdType start = afc[0];
  int toIndex = 0;
  int from = 0;

  for(int j1 = 1; j1 < loopID; j1++)
  {
    int openL = 0;
    int flip = 0;
    SiteIdType startEdge = -1;
    const int numN = count[j1];
    toIndex = toIndex + numN;
    from = toIndex - numN;
    std::vector<SiteIdType> burntLoop(static_cast<usize>(numN) + 2, 0);

    for(int i1 = from; i1 < toIndex; i1++)
    {
      const SiteIdType cedge = burntList[i1];
      const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
      const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
      if(start == cnode1)
      {
        openL = 1;
        startEdge = cedge;
        flip = 0;
      }
      else if(start == cnode2)
      {
        openL = 1;
        startEdge = cedge;
        flip = 1;
      }
    }

    if(openL == 1)
    {
      if(flip == 1)
      {
        const SiteIdType tnode = faceEdges[startEdge].NodeId[0];
        const int tspin = faceEdges[startEdge].NSpin[0];
        faceEdges[startEdge].NodeId[0] = faceEdges[startEdge].NodeId[1];
        faceEdges[startEdge].NodeId[1] = tnode;
        faceEdges[startEdge].NSpin[0] = faceEdges[startEdge].NSpin[1];
        faceEdges[startEdge].NSpin[1] = tspin;
      }

      burntLoop[0] = startEdge;
      int index = 1;
      SiteIdType endNode = faceEdges[startEdge].NodeId[1];
      SiteIdType chaser = startEdge;
      do
      {
        const int passStart = index; // chase-loop guard: detect a pass that fails to extend the chain
        for(int burntEdgeIdx = from; burntEdgeIdx < toIndex; burntEdgeIdx++)
        {
          const SiteIdType cedge = burntList[burntEdgeIdx];
          const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
          const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
          if((cedge != chaser) && (endNode == cnode1))
          {
            burntLoop[index] = cedge;
            index++;
          }
          else if((cedge != chaser) && (endNode == cnode2))
          {
            burntLoop[index] = cedge;
            index++;
            const SiteIdType tnode = faceEdges[cedge].NodeId[0];
            const int tspin = faceEdges[cedge].NSpin[0];
            faceEdges[cedge].NodeId[0] = faceEdges[cedge].NodeId[1];
            faceEdges[cedge].NodeId[1] = tnode;
            faceEdges[cedge].NSpin[0] = faceEdges[cedge].NSpin[1];
            faceEdges[cedge].NSpin[1] = tspin;
          }
          if(index >= numN)
          {
            break; // chain complete; also caps degenerate multi-match passes so burntLoop cannot overrun
          }
        }
        if(index == passStart)
        {
          break; // degenerate input: the pass matched no edge, so the chain can never close
        }
        chaser = burntLoop[index - 1];
        endNode = faceEdges[chaser].NodeId[1];
      } while(index < numN);

      if((numN + 1) == 3)
      {
        numTri = numTri + 1;
      }
      else if((numN + 1) > 3)
      {
        numTri = numTri + ((numN + 1) - 2);
      }
    }
    else
    {
      const SiteIdType startEdge2 = burntList[from];
      burntLoop[0] = startEdge2;
      int index = 1;
      SiteIdType endNode = faceEdges[startEdge2].NodeId[1];
      SiteIdType chaser = startEdge2;
      do
      {
        const int passStart = index; // chase-loop guard: detect a pass that fails to extend the chain
        for(int burntEdgeIdx = from; burntEdgeIdx < toIndex; burntEdgeIdx++)
        {
          const SiteIdType cedge = burntList[burntEdgeIdx];
          const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
          const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
          if((cedge != chaser) && (endNode == cnode1))
          {
            burntLoop[index] = cedge;
            index++;
          }
          else if((cedge != chaser) && (endNode == cnode2))
          {
            burntLoop[index] = cedge;
            index++;
            const SiteIdType tnode = faceEdges[cedge].NodeId[0];
            const int tspin = faceEdges[cedge].NSpin[0];
            faceEdges[cedge].NodeId[0] = faceEdges[cedge].NodeId[1];
            faceEdges[cedge].NodeId[1] = tnode;
            faceEdges[cedge].NSpin[0] = faceEdges[cedge].NSpin[1];
            faceEdges[cedge].NSpin[1] = tspin;
          }
          if(index >= numN)
          {
            break; // chain complete; also caps degenerate multi-match passes so burntLoop cannot overrun
          }
        }
        if(index == passStart)
        {
          break; // degenerate input: the pass matched no edge, so the chain can never close
        }
        chaser = burntLoop[index - 1];
        endNode = faceEdges[chaser].NodeId[1];
      } while(index < numN);

      if(numN == 3)
      {
        numTri = numTri + 1;
      }
      else if(numN > 3)
      {
        numTri = numTri + (numN - 2);
      }
    }
  }
  return numTri;
}

/**
 * @brief Counts triangles for a cube with three or more face centers.
 * @param afe Provides cube face-edge indexes.
 * @param faceEdges Provides mutable oriented face-edge records.
 * @param nfedge Specifies cube face-edge count.
 * @param afc Provides face-center node indexes.
 * @param nfctr Specifies face-center count.
 * @return Triangle count after body-center and closed-loop triangulation.
 */
int getNumberCaseMTriangles(const std::span<const SiteIdType> afe, const std::span<Segment> faceEdges, const int nfedge, const std::array<SiteIdType, 6>& afc, const int nfctr)
{
  const FaceEdgeLoops loops = burnFaceEdgeLoops(afe, faceEdges, nfedge);
  const std::vector<SiteIdType>& burntList = loops.BurntList;
  const std::vector<int>& count = loops.Count;
  const int loopID = loops.LoopId;

  int numTri = 0;
  int toIndex = 0;
  int from = 0;

  for(int j1 = 1; j1 < loopID; j1++)
  {
    int openL = 0;
    int flip = 0;
    SiteIdType startEdge = -1;
    const int numN = count[j1];
    toIndex = toIndex + numN;
    from = toIndex - numN;
    std::vector<SiteIdType> burntLoop(static_cast<usize>(numN) + 2, 0);

    for(int i1 = from; i1 < toIndex; i1++)
    {
      const SiteIdType cedge = burntList[i1];
      const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
      const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
      for(int n1 = 0; n1 < nfctr; n1++)
      {
        const SiteIdType start = afc[n1];
        if(start == cnode1)
        {
          openL = 1;
          startEdge = cedge;
          flip = 0;
        }
        else if(start == cnode2)
        {
          openL = 1;
          startEdge = cedge;
          flip = 1;
        }
      }
    }

    if(openL == 1)
    {
      if(flip == 1)
      {
        const SiteIdType tnode = faceEdges[startEdge].NodeId[0];
        const int tspin = faceEdges[startEdge].NSpin[0];
        faceEdges[startEdge].NodeId[0] = faceEdges[startEdge].NodeId[1];
        faceEdges[startEdge].NodeId[1] = tnode;
        faceEdges[startEdge].NSpin[0] = faceEdges[startEdge].NSpin[1];
        faceEdges[startEdge].NSpin[1] = tspin;
      }

      burntLoop[0] = startEdge;
      int index = 1;
      SiteIdType endNode = faceEdges[startEdge].NodeId[1];
      SiteIdType chaser = startEdge;
      do
      {
        const int passStart = index; // chase-loop guard: detect a pass that fails to extend the chain
        for(int burntEdgeIdx = from; burntEdgeIdx < toIndex; burntEdgeIdx++)
        {
          const SiteIdType cedge = burntList[burntEdgeIdx];
          const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
          const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
          if((cedge != chaser) && (endNode == cnode1))
          {
            burntLoop[index] = cedge;
            index++;
          }
          else if((cedge != chaser) && (endNode == cnode2))
          {
            burntLoop[index] = cedge;
            index++;
            const SiteIdType tnode = faceEdges[cedge].NodeId[0];
            const int tspin = faceEdges[cedge].NSpin[0];
            faceEdges[cedge].NodeId[0] = faceEdges[cedge].NodeId[1];
            faceEdges[cedge].NodeId[1] = tnode;
            faceEdges[cedge].NSpin[0] = faceEdges[cedge].NSpin[1];
            faceEdges[cedge].NSpin[1] = tspin;
          }
          if(index >= numN)
          {
            break; // chain complete; also caps degenerate multi-match passes so burntLoop cannot overrun
          }
        }
        if(index == passStart)
        {
          break; // degenerate input: the pass matched no edge, so the chain can never close
        }
        chaser = burntLoop[index - 1];
        endNode = faceEdges[chaser].NodeId[1];
      } while(index < numN);

      if((numN + 2) == 3)
      {
        numTri = numTri + 1;
      }
      else if((numN + 2) > 3)
      {
        numTri = numTri + ((numN + 2) - 2);
      }
    }
    else
    {
      const SiteIdType startEdge2 = burntList[from];
      burntLoop[0] = startEdge2;
      int index = 1;
      SiteIdType endNode = faceEdges[startEdge2].NodeId[1];
      SiteIdType chaser = startEdge2;
      do
      {
        const int passStart = index; // chase-loop guard: detect a pass that fails to extend the chain
        for(int burntEdgeIdx = from; burntEdgeIdx < toIndex; burntEdgeIdx++)
        {
          const SiteIdType cedge = burntList[burntEdgeIdx];
          const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
          const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
          if((cedge != chaser) && (endNode == cnode1))
          {
            burntLoop[index] = cedge;
            index++;
          }
          else if((cedge != chaser) && (endNode == cnode2))
          {
            burntLoop[index] = cedge;
            index++;
            const SiteIdType tnode = faceEdges[cedge].NodeId[0];
            const int tspin = faceEdges[cedge].NSpin[0];
            faceEdges[cedge].NodeId[0] = faceEdges[cedge].NodeId[1];
            faceEdges[cedge].NodeId[1] = tnode;
            faceEdges[cedge].NSpin[0] = faceEdges[cedge].NSpin[1];
            faceEdges[cedge].NSpin[1] = tspin;
          }
          if(index >= numN)
          {
            break; // chain complete; also caps degenerate multi-match passes so burntLoop cannot overrun
          }
        }
        if(index == passStart)
        {
          break; // degenerate input: the pass matched no edge, so the chain can never close
        }
        chaser = burntLoop[index - 1];
        endNode = faceEdges[chaser].NodeId[1];
      } while(index < numN);

      if(numN == 3)
      {
        numTri = numTri + 1;
      }
      else if(numN > 3)
      {
        numTri = numTri + (numN - 2);
      }
    }
  }
  return numTri;
}

/**
 * @brief Counts all triangles and classifies body-center nodes.
 * @param featureIds Provides padded Feature Id values.
 * @param squares Provides marching-square records.
 * @param neighbors Calculates padded-grid neighbors.
 * @param nodeType Receives body-center node categories.
 * @param faceEdges Provides mutable oriented face-edge records.
 * @param numSitesDim3 Specifies padded site count.
 * @param numSitesDim2 Specifies padded sites per Z plane.
 * @param xDim Specifies padded X dimension.
 * @param shouldCancel Stops before later cubes when true.
 * @return Triangle count accumulated before completion or cancellation.
 */
int64 getNumberTriangles(const std::vector<int32>& featureIds, const std::vector<Face>& squares, const NeighborAccessor& neighbors, std::vector<int8>& nodeType, std::vector<Segment>& faceEdges,
                         const SiteIdType numSitesDim3, const SiteIdType numSitesDim2, const int xDim, const std::atomic_bool& shouldCancel)
{
  int64 nTri0 = 0;
  int64 nTri2 = 0;
  int64 nTriM = 0;

  for(SiteIdType i = 1; i <= (numSitesDim3 - numSitesDim2); i++)
  {
    if(shouldCancel)
    {
      return 0;
    }
    int cubeFlag = 0;
    std::array<SiteIdType, 6> sqID{};
    sqID[0] = 3 * (i - 1);
    sqID[1] = (3 * (i - 1)) + 1;
    sqID[2] = (3 * (i - 1)) + 2;
    sqID[3] = (3 * i) + 2;
    sqID[4] = (3 * (i + xDim - 1)) + 1;
    sqID[5] = 3 * (i + numSitesDim2 - 1);
    const SiteIdType bodyCenterNode = (7 * (i - 1)) + 6;
    int nFC = 0;
    int nFE = 0;
    int eff = 0;
    std::array<SiteIdType, 6> arrayFC{};
    for(int ii = 0; ii < 6; ii++)
    {
      arrayFC[ii] = -1;
    }
    int fcid = 0;
    for(int ii = 0; ii < 6; ii++)
    {
      const SiteIdType tsq = sqID[ii];
      const SiteIdType tFCnode = squares[tsq].FaceCenterNode;
      if(tFCnode != -1)
      {
        arrayFC[fcid] = tFCnode;
        fcid++;
      }
      nFE = nFE + squares[tsq].NEdge;
      eff = eff + squares[tsq].Effect;
    }
    nFC = fcid;
    if(eff > 0)
    {
      cubeFlag = 1;
    }

    if(nFC >= 3)
    {
      const std::array<SiteIdType, 4> corners1 = squareCorners(sqID[0], neighbors);
      const std::array<SiteIdType, 4> corners2 = squareCorners(sqID[5], neighbors);
      std::array<int, 8> arraySpin{};
      for(int j = 0; j < 4; j++)
      {
        arraySpin[j] = featureIds[corners1[j]];
        arraySpin[j + 4] = featureIds[corners2[j]];
      }
      int nds = 0;
      int nburnt = 0;
      for(int k = 0; k < 8; k++)
      {
        const int cspin = arraySpin[k];
        if(cspin != -1)
        {
          nds++;
          arraySpin[k] = -1;
          nburnt++;
          for(int kk = 0; kk < 8; kk++)
          {
            if(cspin == arraySpin[kk])
            {
              arraySpin[kk] = -1;
              nburnt++;
            }
          }
        }
      }
      (void)nburnt;
      // Five or more labels can meet at a body center. NodeType supports only
      // the "four or more" category used by downstream mesh consumers.
      nodeType[bodyCenterNode] = static_cast<int8>(std::min(nds, static_cast<int>(m3c_node_type::k_QuadPoint)));
    }

    if(cubeFlag == 1 && nFE > 2)
    {
      std::vector<SiteIdType> arrayFE(nFE);
      int tindex = 0;
      for(int i1 = 0; i1 < 6; i1++)
      {
        const SiteIdType tsq = sqID[i1];
        const int tnfe = static_cast<int>(static_cast<uint8>(squares[tsq].NEdge));
        for(int i2 = 0; i2 < tnfe; i2++)
        {
          arrayFE[tindex] = squares[tsq].EdgeId[i2];
          tindex++;
        }
      }

      // Square cases determine face-center count. A cube can have zero or two
      // through six centers. One crossing cannot terminate inside one cube.
      if(nFC == 0)
      {
        nTri0 = nTri0 + getNumberCase0Triangles(arrayFE, faceEdges, nFE);
      }
      else if(nFC == 2)
      {
        nTri2 = nTri2 + getNumberCase2Triangles(arrayFE, faceEdges, nFE, arrayFC, nFC);
      }
      else if(nFC > 2 && nFC <= 6)
      {
        nTriM = nTriM + getNumberCaseMTriangles(arrayFE, faceEdges, nFE, arrayFC, nFC);
      }
    }
  }
  return nTri0 + nTri2 + nTriM;
}

/**
 * @brief Generates triangles for a cube without face centers.
 * @param triangles Receives triangle records.
 * @param mCubeID Receives the source cube for each triangle.
 * @param afe Provides cube face-edge indexes.
 * @param nodeCoords Is retained by the legacy call shape and is unused.
 * @param faceEdges Provides mutable oriented face-edge records.
 * @param nfedge Specifies cube face-edge count.
 * @param tin Specifies the first output triangle index.
 * @param tout Receives the next unused output triangle index.
 * @param tcrd1 Is retained by the legacy call shape and is unused.
 * @param tcrd2 Is retained by the legacy call shape and is unused.
 * @param mcid Specifies the source cube index.
 */
void getCase0Triangles(std::vector<Triangle>& triangles, std::vector<SiteIdType>& mCubeID, const std::span<const SiteIdType> afe, const NodeCoords& nodeCoords, const std::span<Segment> faceEdges,
                       const int nfedge, const int64 tin, int64& tout, const std::array<double, 3>& tcrd1, const std::array<double, 3>& tcrd2, const SiteIdType mcid)
{

  const FaceEdgeLoops loops = burnAndOrientFaceEdgeLoops(afe, faceEdges, nfedge);
  const std::vector<SiteIdType>& burntList = loops.BurntList;
  const std::vector<int>& count = loops.Count;
  const int loopID = loops.LoopId;

  int sumN = 0;
  int64 ctid = tin;

  for(int jj = 1; jj < loopID; jj++)
  {
    const int numN = count[jj];
    sumN = sumN + numN;
    const int from = sumN - numN;
    std::vector<SiteIdType> loop(numN);
    for(int mm = 0; mm < numN; mm++)
    {
      loop[mm] = burntList[from + mm];
    }

    if(numN == 3)
    {
      const SiteIdType te0 = loop[0], te1 = loop[1], te2 = loop[2];
      const SiteIdType tv0 = faceEdges[te0].NodeId[0];
      const SiteIdType tv1 = faceEdges[te1].NodeId[0];
      const SiteIdType tv2 = faceEdges[te2].NodeId[0];
      triangles[ctid].NodeId[0] = tv0;
      triangles[ctid].NodeId[1] = tv1;
      triangles[ctid].NodeId[2] = tv2;
      triangles[ctid].NSpin[0] = faceEdges[te0].NSpin[0];
      triangles[ctid].NSpin[1] = faceEdges[te0].NSpin[1];
      mCubeID[ctid] = mcid;
      ctid++;
    }
    else if(numN > 3)
    {
      const int numT = numN - 2;
      int cnumT = 0;
      int front = 0;
      int back = numN - 1;

      const SiteIdType te0 = loop[front];
      const SiteIdType te1 = loop[back];
      SiteIdType tv0 = faceEdges[te0].NodeId[0];
      SiteIdType tv1 = faceEdges[te0].NodeId[1];
      SiteIdType tv2 = faceEdges[te1].NodeId[0];
      triangles[ctid].NodeId[0] = tv0;
      triangles[ctid].NodeId[1] = tv1;
      triangles[ctid].NodeId[2] = tv2;
      triangles[ctid].NSpin[0] = faceEdges[te0].NSpin[0];
      triangles[ctid].NSpin[1] = faceEdges[te0].NSpin[1];
      mCubeID[ctid] = mcid;
      SiteIdType newNode0 = tv2;
      cnumT++;
      ctid++;

      do
      {
        if((cnumT % 2) != 0)
        {
          front = front + 1;
          const SiteIdType currentEdge = loop[front];
          tv0 = faceEdges[currentEdge].NodeId[0];
          tv1 = faceEdges[currentEdge].NodeId[1];
          tv2 = newNode0;
          triangles[ctid].NodeId[0] = tv0;
          triangles[ctid].NodeId[1] = tv1;
          triangles[ctid].NodeId[2] = tv2;
          triangles[ctid].NSpin[0] = faceEdges[currentEdge].NSpin[0];
          triangles[ctid].NSpin[1] = faceEdges[currentEdge].NSpin[1];
          mCubeID[ctid] = mcid;
          newNode0 = tv1;
          cnumT++;
          ctid++;
        }
        else
        {
          back = back - 1;
          const SiteIdType currentEdge = loop[back];
          tv0 = faceEdges[currentEdge].NodeId[0];
          tv1 = faceEdges[currentEdge].NodeId[1];
          tv2 = newNode0;
          triangles[ctid].NodeId[0] = tv0;
          triangles[ctid].NodeId[1] = tv1;
          triangles[ctid].NodeId[2] = tv2;
          triangles[ctid].NSpin[0] = faceEdges[currentEdge].NSpin[0];
          triangles[ctid].NSpin[1] = faceEdges[currentEdge].NSpin[1];
          mCubeID[ctid] = mcid;
          newNode0 = tv0;
          cnumT++;
          ctid++;
        }
      } while(cnumT < numT);
    }
  }
  tout = ctid;
}

/**
 * @brief Generates triangles for a cube with two face centers.
 * @param triangles Receives triangle records.
 * @param mCubeID Receives the source cube for each triangle.
 * @param afe Provides cube face-edge indexes.
 * @param nodeCoords Is retained by the legacy call shape and is unused.
 * @param faceEdges Provides mutable oriented face-edge records.
 * @param nfedge Specifies cube face-edge count.
 * @param afc Provides face-center node indexes.
 * @param nfctr Is fixed at two and is unused.
 * @param tin Specifies the first output triangle index.
 * @param tout Receives the next unused output triangle index.
 * @param tcrd1 Is retained by the legacy call shape and is unused.
 * @param tcrd2 Is retained by the legacy call shape and is unused.
 * @param mcid Specifies the source cube index.
 */
void getCase2Triangles(std::vector<Triangle>& triangles, std::vector<SiteIdType>& mCubeID, const std::span<const SiteIdType> afe, const NodeCoords& nodeCoords, const std::span<Segment> faceEdges,
                       const int nfedge, const std::array<SiteIdType, 6>& afc, int /*nfctr*/, const int64 tin, int64& tout, const std::array<double, 3>& tcrd1, const std::array<double, 3>& tcrd2,
                       const SiteIdType mcid)
{

  const FaceEdgeLoops loops = burnFaceEdgeLoops(afe, faceEdges, nfedge);
  const std::vector<SiteIdType>& burntList = loops.BurntList;
  const std::vector<int>& count = loops.Count;
  const int loopID = loops.LoopId;

  const SiteIdType start = afc[0];
  int toIndex = 0;
  int from = 0;
  int64 ctid = tin;

  for(int j1 = 1; j1 < loopID; j1++)
  {
    int openL = 0;
    int flip = 0;
    SiteIdType startEdge = -1;
    const int numN = count[j1];
    toIndex = toIndex + numN;
    from = toIndex - numN;
    std::vector<SiteIdType> burntLoop(static_cast<usize>(numN) + 2, 0);

    for(int i1 = from; i1 < toIndex; i1++)
    {
      const SiteIdType cedge = burntList[i1];
      const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
      const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
      if(start == cnode1)
      {
        openL = 1;
        startEdge = cedge;
        flip = 0;
      }
      else if(start == cnode2)
      {
        openL = 1;
        startEdge = cedge;
        flip = 1;
      }
    }

    if(openL == 1)
    {
      if(flip == 1)
      {
        const SiteIdType tnode = faceEdges[startEdge].NodeId[0];
        const int tspin = faceEdges[startEdge].NSpin[0];
        faceEdges[startEdge].NodeId[0] = faceEdges[startEdge].NodeId[1];
        faceEdges[startEdge].NodeId[1] = tnode;
        faceEdges[startEdge].NSpin[0] = faceEdges[startEdge].NSpin[1];
        faceEdges[startEdge].NSpin[1] = tspin;
      }
      burntLoop[0] = startEdge;
      int index = 1;
      SiteIdType endNode = faceEdges[startEdge].NodeId[1];
      SiteIdType chaser = startEdge;
      do
      {
        const int passStart = index; // chase-loop guard: detect a pass that fails to extend the chain
        for(int burntEdgeIdx = from; burntEdgeIdx < toIndex; burntEdgeIdx++)
        {
          const SiteIdType cedge = burntList[burntEdgeIdx];
          const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
          const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
          if((cedge != chaser) && (endNode == cnode1))
          {
            burntLoop[index] = cedge;
            index++;
          }
          else if((cedge != chaser) && (endNode == cnode2))
          {
            burntLoop[index] = cedge;
            index++;
            const SiteIdType tnode = faceEdges[cedge].NodeId[0];
            const int tspin = faceEdges[cedge].NSpin[0];
            faceEdges[cedge].NodeId[0] = faceEdges[cedge].NodeId[1];
            faceEdges[cedge].NodeId[1] = tnode;
            faceEdges[cedge].NSpin[0] = faceEdges[cedge].NSpin[1];
            faceEdges[cedge].NSpin[1] = tspin;
          }
          if(index >= numN)
          {
            break; // chain complete; also caps degenerate multi-match passes so burntLoop cannot overrun
          }
        }
        if(index == passStart)
        {
          break; // degenerate input: the pass matched no edge, so the chain can never close
        }
        chaser = burntLoop[index - 1];
        endNode = faceEdges[chaser].NodeId[1];
      } while(index < numN);

      if(numN == 2)
      {
        const SiteIdType te0 = burntLoop[0], te1 = burntLoop[1];
        const SiteIdType tv0 = faceEdges[te0].NodeId[0];
        const SiteIdType tv1 = faceEdges[te1].NodeId[0];
        const SiteIdType tv2 = faceEdges[te1].NodeId[1];
        triangles[ctid].NodeId[0] = tv0;
        triangles[ctid].NodeId[1] = tv1;
        triangles[ctid].NodeId[2] = tv2;
        triangles[ctid].NSpin[0] = faceEdges[te0].NSpin[0];
        triangles[ctid].NSpin[1] = faceEdges[te0].NSpin[1];
        mCubeID[ctid] = mcid;
        ctid++;
      }
      else if(numN > 2)
      {
        const int numT = numN - 1;
        int cnumT = 0;
        int front = 0;
        int back = numN;
        const SiteIdType te0 = burntLoop[front];
        const SiteIdType te1 = burntLoop[back - 1];
        SiteIdType tv0 = faceEdges[te0].NodeId[0];
        SiteIdType tv1 = faceEdges[te0].NodeId[1];
        SiteIdType tv2 = faceEdges[te1].NodeId[1];
        triangles[ctid].NodeId[0] = tv0;
        triangles[ctid].NodeId[1] = tv1;
        triangles[ctid].NodeId[2] = tv2;
        triangles[ctid].NSpin[0] = faceEdges[te0].NSpin[0];
        triangles[ctid].NSpin[1] = faceEdges[te0].NSpin[1];
        mCubeID[ctid] = mcid;
        SiteIdType newNode0 = tv2;
        cnumT++;
        ctid++;
        do
        {
          if((cnumT % 2) != 0)
          {
            front = front + 1;
            const SiteIdType currentEdge = burntLoop[front];
            tv0 = faceEdges[currentEdge].NodeId[0];
            tv1 = faceEdges[currentEdge].NodeId[1];
            tv2 = newNode0;
            triangles[ctid].NodeId[0] = tv0;
            triangles[ctid].NodeId[1] = tv1;
            triangles[ctid].NodeId[2] = tv2;
            triangles[ctid].NSpin[0] = faceEdges[currentEdge].NSpin[0];
            triangles[ctid].NSpin[1] = faceEdges[currentEdge].NSpin[1];
            mCubeID[ctid] = mcid;
            newNode0 = tv1;
            cnumT++;
            ctid++;
          }
          else
          {
            back = back - 1;
            const SiteIdType currentEdge = burntLoop[back];
            tv0 = faceEdges[currentEdge].NodeId[0];
            tv1 = faceEdges[currentEdge].NodeId[1];
            tv2 = newNode0;
            triangles[ctid].NodeId[0] = tv0;
            triangles[ctid].NodeId[1] = tv1;
            triangles[ctid].NodeId[2] = tv2;
            triangles[ctid].NSpin[0] = faceEdges[currentEdge].NSpin[0];
            triangles[ctid].NSpin[1] = faceEdges[currentEdge].NSpin[1];
            mCubeID[ctid] = mcid;
            newNode0 = tv0;
            cnumT++;
            ctid++;
          }
        } while(cnumT < numT);
      }
    }
    else
    {
      const SiteIdType startEdge2 = burntList[from];
      burntLoop[0] = startEdge2;
      int index = 1;
      SiteIdType endNode = faceEdges[startEdge2].NodeId[1];
      SiteIdType chaser = startEdge2;
      do
      {
        const int passStart = index; // chase-loop guard: detect a pass that fails to extend the chain
        for(int burntEdgeIdx = from; burntEdgeIdx < toIndex; burntEdgeIdx++)
        {
          const SiteIdType cedge = burntList[burntEdgeIdx];
          const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
          const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
          if((cedge != chaser) && (endNode == cnode1))
          {
            burntLoop[index] = cedge;
            index++;
          }
          else if((cedge != chaser) && (endNode == cnode2))
          {
            burntLoop[index] = cedge;
            index++;
            const SiteIdType tnode = faceEdges[cedge].NodeId[0];
            const int tspin = faceEdges[cedge].NSpin[0];
            faceEdges[cedge].NodeId[0] = faceEdges[cedge].NodeId[1];
            faceEdges[cedge].NodeId[1] = tnode;
            faceEdges[cedge].NSpin[0] = faceEdges[cedge].NSpin[1];
            faceEdges[cedge].NSpin[1] = tspin;
          }
          if(index >= numN)
          {
            break; // chain complete; also caps degenerate multi-match passes so burntLoop cannot overrun
          }
        }
        if(index == passStart)
        {
          break; // degenerate input: the pass matched no edge, so the chain can never close
        }
        chaser = burntLoop[index - 1];
        endNode = faceEdges[chaser].NodeId[1];
      } while(index < numN);

      if(numN == 3)
      {
        const SiteIdType te0 = burntLoop[0], te1 = burntLoop[1], te2 = burntLoop[2];
        const SiteIdType tv0 = faceEdges[te0].NodeId[0];
        const SiteIdType tv1 = faceEdges[te1].NodeId[0];
        const SiteIdType tv2 = faceEdges[te2].NodeId[0];
        triangles[ctid].NodeId[0] = tv0;
        triangles[ctid].NodeId[1] = tv1;
        triangles[ctid].NodeId[2] = tv2;
        triangles[ctid].NSpin[0] = faceEdges[te0].NSpin[0];
        triangles[ctid].NSpin[1] = faceEdges[te0].NSpin[1];
        mCubeID[ctid] = mcid;
        ctid++;
      }
      else if(numN > 3)
      {
        const int numT = numN - 2;
        int cnumT = 0;
        int front = 0;
        int back = numN - 1;
        const SiteIdType te0 = burntLoop[front];
        const SiteIdType te1 = burntLoop[back];
        SiteIdType tv0 = faceEdges[te0].NodeId[0];
        SiteIdType tv1 = faceEdges[te0].NodeId[1];
        SiteIdType tv2 = faceEdges[te1].NodeId[0];
        triangles[ctid].NodeId[0] = tv0;
        triangles[ctid].NodeId[1] = tv1;
        triangles[ctid].NodeId[2] = tv2;
        triangles[ctid].NSpin[0] = faceEdges[te0].NSpin[0];
        triangles[ctid].NSpin[1] = faceEdges[te0].NSpin[1];
        mCubeID[ctid] = mcid;
        SiteIdType newNode0 = tv2;
        cnumT++;
        ctid++;
        do
        {
          if((cnumT % 2) != 0)
          {
            front = front + 1;
            const SiteIdType currentEdge = burntLoop[front];
            tv0 = faceEdges[currentEdge].NodeId[0];
            tv1 = faceEdges[currentEdge].NodeId[1];
            tv2 = newNode0;
            triangles[ctid].NodeId[0] = tv0;
            triangles[ctid].NodeId[1] = tv1;
            triangles[ctid].NodeId[2] = tv2;
            triangles[ctid].NSpin[0] = faceEdges[currentEdge].NSpin[0];
            triangles[ctid].NSpin[1] = faceEdges[currentEdge].NSpin[1];
            mCubeID[ctid] = mcid;
            newNode0 = tv1;
            cnumT++;
            ctid++;
          }
          else
          {
            back = back - 1;
            const SiteIdType currentEdge = burntLoop[back];
            tv0 = faceEdges[currentEdge].NodeId[0];
            tv1 = faceEdges[currentEdge].NodeId[1];
            tv2 = newNode0;
            triangles[ctid].NodeId[0] = tv0;
            triangles[ctid].NodeId[1] = tv1;
            triangles[ctid].NodeId[2] = tv2;
            triangles[ctid].NSpin[0] = faceEdges[currentEdge].NSpin[0];
            triangles[ctid].NSpin[1] = faceEdges[currentEdge].NSpin[1];
            mCubeID[ctid] = mcid;
            newNode0 = tv0;
            cnumT++;
            ctid++;
          }
        } while(cnumT < numT);
      }
    }
  }
  tout = ctid;
}

/**
 * @brief Generates triangles for a cube with three or more face centers.
 * @param triangles Receives triangle records.
 * @param mCubeID Receives the source cube for each triangle.
 * @param afe Provides cube face-edge indexes.
 * @param nodeCoords Is retained by the legacy call shape and is unused.
 * @param faceEdges Provides mutable oriented face-edge records.
 * @param nfedge Specifies cube face-edge count.
 * @param afc Provides face-center node indexes.
 * @param nfctr Specifies face-center count.
 * @param tin Specifies the first output triangle index.
 * @param tout Receives the next unused output triangle index.
 * @param ccn Specifies the body-center candidate node.
 * @param tcrd1 Is retained by the legacy call shape and is unused.
 * @param tcrd2 Is retained by the legacy call shape and is unused.
 * @param mcid Specifies the source cube index.
 *
 * Open loops use a fan from the body-center node.
 */
void getCaseMTriangles(std::vector<Triangle>& triangles, std::vector<SiteIdType>& mCubeID, const std::span<const SiteIdType> afe, const NodeCoords& nodeCoords, const std::span<Segment> faceEdges,
                       const int nfedge, const std::array<SiteIdType, 6>& afc, const int nfctr, const int64 tin, int64& tout, const SiteIdType ccn, const std::array<double, 3>& tcrd1,
                       const std::array<double, 3>& tcrd2, const SiteIdType mcid)
{

  const FaceEdgeLoops loops = burnFaceEdgeLoops(afe, faceEdges, nfedge);
  const std::vector<SiteIdType>& burntList = loops.BurntList;
  const std::vector<int>& count = loops.Count;
  const int loopID = loops.LoopId;

  int toIndex = 0;
  int from = 0;
  int64 ctid = tin;

  for(int j1 = 1; j1 < loopID; j1++)
  {
    int openL = 0;
    int flip = 0;
    SiteIdType startEdge = -1;
    const int numN = count[j1];
    toIndex = toIndex + numN;
    from = toIndex - numN;
    std::vector<SiteIdType> burntLoop(static_cast<usize>(numN) + 2, 0);

    for(int i1 = from; i1 < toIndex; i1++)
    {
      const SiteIdType cedge = burntList[i1];
      const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
      const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
      for(int n1 = 0; n1 < nfctr; n1++)
      {
        const SiteIdType start = afc[n1];
        if(start == cnode1)
        {
          openL = 1;
          startEdge = cedge;
          flip = 0;
        }
        else if(start == cnode2)
        {
          openL = 1;
          startEdge = cedge;
          flip = 1;
        }
      }
    }

    if(openL == 1)
    {
      if(flip == 1)
      {
        const SiteIdType tnode = faceEdges[startEdge].NodeId[0];
        const int tspin = faceEdges[startEdge].NSpin[0];
        faceEdges[startEdge].NodeId[0] = faceEdges[startEdge].NodeId[1];
        faceEdges[startEdge].NodeId[1] = tnode;
        faceEdges[startEdge].NSpin[0] = faceEdges[startEdge].NSpin[1];
        faceEdges[startEdge].NSpin[1] = tspin;
      }
      burntLoop[0] = startEdge;
      int index = 1;
      SiteIdType endNode = faceEdges[startEdge].NodeId[1];
      SiteIdType chaser = startEdge;
      do
      {
        const int passStart = index; // chase-loop guard: detect a pass that fails to extend the chain
        for(int burntEdgeIdx = from; burntEdgeIdx < toIndex; burntEdgeIdx++)
        {
          const SiteIdType cedge = burntList[burntEdgeIdx];
          const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
          const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
          if((cedge != chaser) && (endNode == cnode1))
          {
            burntLoop[index] = cedge;
            index++;
          }
          else if((cedge != chaser) && (endNode == cnode2))
          {
            burntLoop[index] = cedge;
            index++;
            const SiteIdType tnode = faceEdges[cedge].NodeId[0];
            const int tspin = faceEdges[cedge].NSpin[0];
            faceEdges[cedge].NodeId[0] = faceEdges[cedge].NodeId[1];
            faceEdges[cedge].NodeId[1] = tnode;
            faceEdges[cedge].NSpin[0] = faceEdges[cedge].NSpin[1];
            faceEdges[cedge].NSpin[1] = tspin;
          }
          if(index >= numN)
          {
            break; // chain complete; also caps degenerate multi-match passes so burntLoop cannot overrun
          }
        }
        if(index == passStart)
        {
          break; // degenerate input: the pass matched no edge, so the chain can never close
        }
        chaser = burntLoop[index - 1];
        endNode = faceEdges[chaser].NodeId[1];
      } while(index < numN);

      // Open loops use a fan from the body-center node.
      for(int iii = 0; iii < numN; iii++)
      {
        const SiteIdType currentEdge = burntLoop[iii];
        const SiteIdType tn0 = faceEdges[currentEdge].NodeId[0];
        const SiteIdType tn1 = faceEdges[currentEdge].NodeId[1];
        const int ts0 = faceEdges[currentEdge].NSpin[0];
        const int ts1 = faceEdges[currentEdge].NSpin[1];
        triangles[ctid].NodeId[0] = ccn;
        triangles[ctid].NodeId[1] = tn0;
        triangles[ctid].NodeId[2] = tn1;
        triangles[ctid].NSpin[0] = ts0;
        triangles[ctid].NSpin[1] = ts1;
        mCubeID[ctid] = mcid;
        ctid++;
      }
    }
    else
    {
      const SiteIdType startEdge2 = burntList[from];
      burntLoop[0] = startEdge2;
      int index = 1;
      SiteIdType endNode = faceEdges[startEdge2].NodeId[1];
      SiteIdType chaser = startEdge2;
      do
      {
        const int passStart = index; // chase-loop guard: detect a pass that fails to extend the chain
        for(int burntEdgeIdx = from; burntEdgeIdx < toIndex; burntEdgeIdx++)
        {
          const SiteIdType cedge = burntList[burntEdgeIdx];
          const SiteIdType cnode1 = faceEdges[cedge].NodeId[0];
          const SiteIdType cnode2 = faceEdges[cedge].NodeId[1];
          if((cedge != chaser) && (endNode == cnode1))
          {
            burntLoop[index] = cedge;
            index++;
          }
          else if((cedge != chaser) && (endNode == cnode2))
          {
            burntLoop[index] = cedge;
            index++;
            const SiteIdType tnode = faceEdges[cedge].NodeId[0];
            const int tspin = faceEdges[cedge].NSpin[0];
            faceEdges[cedge].NodeId[0] = faceEdges[cedge].NodeId[1];
            faceEdges[cedge].NodeId[1] = tnode;
            faceEdges[cedge].NSpin[0] = faceEdges[cedge].NSpin[1];
            faceEdges[cedge].NSpin[1] = tspin;
          }
          if(index >= numN)
          {
            break; // chain complete; also caps degenerate multi-match passes so burntLoop cannot overrun
          }
        }
        if(index == passStart)
        {
          break; // degenerate input: the pass matched no edge, so the chain can never close
        }
        chaser = burntLoop[index - 1];
        endNode = faceEdges[chaser].NodeId[1];
      } while(index < numN);

      if(numN == 3)
      {
        const SiteIdType te0 = burntLoop[0], te1 = burntLoop[1], te2 = burntLoop[2];
        const SiteIdType tv0 = faceEdges[te0].NodeId[0];
        const SiteIdType tv1 = faceEdges[te1].NodeId[0];
        const SiteIdType tv2 = faceEdges[te2].NodeId[0];
        triangles[ctid].NodeId[0] = tv0;
        triangles[ctid].NodeId[1] = tv1;
        triangles[ctid].NodeId[2] = tv2;
        triangles[ctid].NSpin[0] = faceEdges[te0].NSpin[0];
        triangles[ctid].NSpin[1] = faceEdges[te0].NSpin[1];
        mCubeID[ctid] = mcid;
        ctid++;
      }
      else if(numN > 3)
      {
        const int numT = numN - 2;
        int cnumT = 0;
        int front = 0;
        int back = numN - 1;
        const SiteIdType te0 = burntLoop[front];
        const SiteIdType te1 = burntLoop[back];
        SiteIdType tv0 = faceEdges[te0].NodeId[0];
        SiteIdType tv1 = faceEdges[te0].NodeId[1];
        SiteIdType tv2 = faceEdges[te1].NodeId[0];
        triangles[ctid].NodeId[0] = tv0;
        triangles[ctid].NodeId[1] = tv1;
        triangles[ctid].NodeId[2] = tv2;
        triangles[ctid].NSpin[0] = faceEdges[te0].NSpin[0];
        triangles[ctid].NSpin[1] = faceEdges[te0].NSpin[1];
        mCubeID[ctid] = mcid;
        SiteIdType newNode0 = tv2;
        cnumT++;
        ctid++;
        do
        {
          if((cnumT % 2) != 0)
          {
            front = front + 1;
            const SiteIdType currentEdge = burntLoop[front];
            tv0 = faceEdges[currentEdge].NodeId[0];
            tv1 = faceEdges[currentEdge].NodeId[1];
            tv2 = newNode0;
            triangles[ctid].NodeId[0] = tv0;
            triangles[ctid].NodeId[1] = tv1;
            triangles[ctid].NodeId[2] = tv2;
            triangles[ctid].NSpin[0] = faceEdges[currentEdge].NSpin[0];
            triangles[ctid].NSpin[1] = faceEdges[currentEdge].NSpin[1];
            mCubeID[ctid] = mcid;
            newNode0 = tv1;
            cnumT++;
            ctid++;
          }
          else
          {
            back = back - 1;
            const SiteIdType currentEdge = burntLoop[back];
            tv0 = faceEdges[currentEdge].NodeId[0];
            tv1 = faceEdges[currentEdge].NodeId[1];
            tv2 = newNode0;
            triangles[ctid].NodeId[0] = tv0;
            triangles[ctid].NodeId[1] = tv1;
            triangles[ctid].NodeId[2] = tv2;
            triangles[ctid].NSpin[0] = faceEdges[currentEdge].NSpin[0];
            triangles[ctid].NSpin[1] = faceEdges[currentEdge].NSpin[1];
            mCubeID[ctid] = mcid;
            newNode0 = tv0;
            cnumT++;
            ctid++;
          }
        } while(cnumT < numT);
      }
    }
  }
  tout = ctid;
}

// -----------------------------------------------------------------------------
// Fill the pre-sized triangle array cube-by-cube. Transcribed from
// M3CEntireVolume::get_triangles.
// -----------------------------------------------------------------------------
// Sharp Bounding Box Edges support.
//
// M3C's candidate nodes sit on a half-cell lattice: an edge-midpoint node has cell-center coordinates on
// two axes and a cell-face coordinate on the third, a face-center node has one cell-center coordinate and
// a body center none. Along a bounding-box edge the marching square straddling it has one real corner and
// three ghost corners, and the case table joins its two edge midpoints with a diagonal: a 45 degree
// chamfer half a cell deep on both walls. The chamfer vertices are exactly the OUTERMOST row of wall
// nodes, because on a wall the only nodes within half a cell of a neighboring wall are those whose
// cell-center coordinate lies in the first or last cell along that axis. Snapping that row onto the
// neighboring wall plane extends both walls to the edge line, where the two rows coincide and are
// merged; the chamfer triangles then reference a repeated node and are dropped.
//
// Everything is decided on the integer half-cell lattice, never on float coordinates, so the pass is
// exact and independent of spacing and origin.
struct HalfCellLattice
{
  const NodeCoords& Coordinates;

  // Position of candidate node `nodeId` in half-cell units from the volume origin, i.e., the node's coordinate
  // is origin + u * spacing / 2. The bounding planes are u == 0 and u == 2 * dims; the outermost rows of
  // cell-center nodes are u == 1 and u == 2 * dims - 1.
  std::array<int64, 3> operator()(const SiteIdType nodeId) const
  {
    const SiteCoords& sites = Coordinates.Sites;
    const usize linear = static_cast<usize>(nodeId / 7);
    const int kind = static_cast<int>(nodeId % 7);
    // Same padded-index decomposition as SiteCoords::operator[] (real cell (0,0,0) is padded (1,1,1)).
    const int64 xIndex = static_cast<int64>(linear % sites.FileDim0) - 1;
    const int64 yIndex = static_cast<int64>((linear / sites.FileDim0) % sites.FileDim1) - 1;
    const int64 zIndex = static_cast<int64>(linear / sites.FileNsp) - 1;
    // Which axes carry the +half-spacing offset for this node kind (see NodeCoords::operator[]).
    const bool offX = (kind == 0 || kind == 3 || kind == 4 || kind == 6);
    const bool offY = (kind == 1 || kind == 3 || kind == 5 || kind == 6);
    const bool offZ = (kind == 2 || kind == 4 || kind == 5 || kind == 6);
    return {(2 * xIndex) + 1 + (offX ? 1 : 0), (2 * yIndex) + 1 + (offY ? 1 : 0), (2 * zIndex) + 1 + (offZ ? 1 : 0)};
  }
};

// Result of the sharp-edge pass: coordinate overrides for the nodes it moved (every other node keeps
// nodeCoords[nodeId]) and the number of chamfer triangles it removed.
struct SharpEdgeResult
{
  std::unordered_map<SiteIdType, Node> SnappedCoords;
  std::unordered_map<SiteIdType, SiteIdType> MergedInto;
  std::unordered_set<SiteIdType> Touched;
  int64 NumFacesRemoved = 0;
};

/**
 * @brief Remaps one triangle and rejects faces collapsed by the sharp-edge pass.
 * @param triangle Receives representative candidate IDs.
 * @param result Supplies node merges and coordinate overrides.
 * @param nodeCoords Supplies coordinates for unchanged candidates.
 * @return True if the remapped triangle survives.
 */
bool remapSharpEdgeTriangle(Triangle& triangle, const SharpEdgeResult& result, const NodeCoords& nodeCoords)
{
  const auto finalCoord = [&result, &nodeCoords](const SiteIdType nodeId) -> Node {
    const auto snappedIter = result.SnappedCoords.find(nodeId);
    return (snappedIter != result.SnappedCoords.end()) ? snappedIter->second : nodeCoords[nodeId];
  };
  int numTouched = 0;
  for(int corner = 0; corner < 3; corner++)
  {
    const auto mergedIter = result.MergedInto.find(triangle.NodeId[corner]);
    if(mergedIter != result.MergedInto.end())
    {
      triangle.NodeId[corner] = mergedIter->second;
    }
    if(result.Touched.count(triangle.NodeId[corner]) != 0)
    {
      numTouched++;
    }
  }
  if(triangle.NodeId[0] == triangle.NodeId[1] || triangle.NodeId[1] == triangle.NodeId[2] || triangle.NodeId[0] == triangle.NodeId[2])
  {
    return false;
  }
  if(numTouched == 3)
  {
    const Node nodeA = finalCoord(triangle.NodeId[0]);
    const Node nodeB = finalCoord(triangle.NodeId[1]);
    const Node nodeC = finalCoord(triangle.NodeId[2]);
    const double abx = static_cast<double>(nodeB.Coord[0]) - nodeA.Coord[0];
    const double aby = static_cast<double>(nodeB.Coord[1]) - nodeA.Coord[1];
    const double abz = static_cast<double>(nodeB.Coord[2]) - nodeA.Coord[2];
    const double acx = static_cast<double>(nodeC.Coord[0]) - nodeA.Coord[0];
    const double acy = static_cast<double>(nodeC.Coord[1]) - nodeA.Coord[1];
    const double acz = static_cast<double>(nodeC.Coord[2]) - nodeA.Coord[2];
    const double crossX = (aby * acz) - (abz * acy);
    const double crossY = (abz * acx) - (abx * acz);
    const double crossZ = (abx * acy) - (aby * acx);
    if(crossX == 0.0 && crossY == 0.0 && crossZ == 0.0)
    {
      return false;
    }
  }
  return true;
}

/**
 * @brief Snaps boundary nodes and removes triangles collapsed along box edges.
 * @tparam NodeTypes Specifies dense or sparse mutable node-type storage.
 * @param triangles Contains surviving triangles and receives the compacted faces.
 * @param mCubeID Contains matching source cubes and receives the compacted cube IDs.
 * @param nodeType Contains exterior-promoted types and receives retired-node markers.
 * @param numCandidateNodes Specifies the dense candidate count when candidateIds is empty.
 * @param nodeCoords Supplies the original node coordinates and half-cell lattice.
 * @param dims Specifies the three image dimensions in cells.
 * @param candidateIds Selects sparse candidates in ascending order, or leaves the dense range selected.
 * @return Coordinate overrides and merge records for output generation.
 * @pre Exterior node promotion is complete, and node compaction has not started.
 */
template <typename NodeTypes>
SharpEdgeResult sharpenBoundingBoxEdges(std::vector<Triangle>& triangles, std::vector<SiteIdType>& mCubeID, NodeTypes& nodeType, const SiteIdType numCandidateNodes, const NodeCoords& nodeCoords,
                                        const std::array<usize, 3>& dims, const nonstd::span<const SiteIdType> candidateIds = {})
{
  SharpEdgeResult result;
  const HalfCellLattice lattice{nodeCoords};
  const SiteCoords& sites = nodeCoords.Sites;
  const std::array<int64, 3> wallHi = {2 * static_cast<int64>(dims[0]), 2 * static_cast<int64>(dims[1]), 2 * static_cast<int64>(dims[2])};
  // Lattice positions packed into one integer for hashing.
  const auto packLattice = [&wallHi](const std::array<int64, 3>& latticePosition) -> uint64 {
    return static_cast<uint64>((((latticePosition[2] * (wallHi[1] + 1)) + latticePosition[1]) * (wallHi[0] + 1)) + latticePosition[0]);
  };

  // Pass 1: for every boundary node decide its snapped lattice position; nodes landing on the same
  // position are merged into the first (lowest id) one to get there, which keeps the pass deterministic.
  std::unordered_map<uint64, SiteIdType> representativeByPosition;
  auto& mergedInto = result.MergedInto;
  const SiteIdType count = candidateIds.empty() ? numCandidateNodes : static_cast<SiteIdType>(candidateIds.size());
  for(SiteIdType index = 0; index < count; index++)
  {
    const SiteIdType nodeId = candidateIds.empty() ? index : candidateIds[static_cast<usize>(index)];
    if(nodeType[static_cast<usize>(nodeId)] < 10)
    {
      continue; // interior node, or unused candidate
    }
    std::array<int64, 3> latticePosition = lattice(nodeId);
    bool onWall = false;
    for(usize ax = 0; ax < 3; ax++)
    {
      onWall = onWall || latticePosition[ax] == 0 || latticePosition[ax] == wallHi[ax];
    }
    if(!onWall)
    {
      continue; // cannot happen for a promoted node; guards the lattice arithmetic
    }
    std::array<bool, 3> snappedAxis = {false, false, false};
    for(usize ax = 0; ax < 3; ax++)
    {
      // A one-cell-thick axis has a single cell-center row that is half a cell from BOTH of its bounding
      // planes; there is no unambiguous edge to snap it to, so that axis is left chamfered.
      if(dims[ax] < 2)
      {
        continue;
      }
      if(latticePosition[ax] == 1)
      {
        latticePosition[ax] = 0;
        snappedAxis[ax] = true;
      }
      else if(latticePosition[ax] == wallHi[ax] - 1)
      {
        latticePosition[ax] = wallHi[ax];
        snappedAxis[ax] = true;
      }
    }
    const auto [representativeIter, inserted] = representativeByPosition.try_emplace(packLattice(latticePosition), nodeId);
    if(inserted)
    {
      if(snappedAxis[0] || snappedAxis[1] || snappedAxis[2])
      {
        // Keep the node's own float coordinates on the axes that did not move, and put it EXACTLY on the
        // plane value the rest of simplnx derives for the volume bounds on the axes that did.
        Node node = nodeCoords[nodeId];
        for(usize ax = 0; ax < 3; ax++)
        {
          if(snappedAxis[ax])
          {
            node.Coord[ax] = (latticePosition[ax] == 0) ? sites.Origin[ax] : sites.Origin[ax] + (static_cast<float>(dims[ax]) * sites.Res[ax]);
          }
        }
        result.SnappedCoords.emplace(nodeId, node);
      }
    }
    else
    {
      const SiteIdType representative = representativeIter->second;
      mergedInto.emplace(nodeId, representative);
      nodeType[static_cast<usize>(representative)] = std::max(nodeType[static_cast<usize>(representative)], nodeType[static_cast<usize>(nodeId)]);
      nodeType[static_cast<usize>(nodeId)] = m3c_node_type::k_Unused;
    }
  }

  if(mergedInto.empty())
  {
    return result;
  }

  // Every node the pass touched (moved or merged into). Used to find the degenerate triangles left on the
  // edge lines, and afterwards to clear any of these nodes no surviving triangle references.
  auto& touched = result.Touched;
  for(const auto& [nodeId, node] : result.SnappedCoords)
  {
    touched.insert(nodeId);
  }
  for(const auto& [nodeId, representative] : mergedInto)
  {
    touched.insert(representative);
  }

  // Pass 2: remap the triangles' node ids and drop the ones the merge collapsed. A chamfer triangle has
  // two vertices on the same cell of the edge line, so after the merge it repeats a node id. The only
  // other way a triangle can lose its area here is for all three vertices to end up on one edge line
  // (exactly collinear, so the cross product is exactly zero); that is checked only for triangles made
  // entirely of touched nodes, which is the sole place it can arise.
  const int64 nTriangle = static_cast<int64>(triangles.size());
  int64 survivingCount = 0;
  for(int64 i = 0; i < nTriangle; i++)
  {
    Triangle triangle = triangles[static_cast<usize>(i)];
    if(!remapSharpEdgeTriangle(triangle, result, nodeCoords))
    {
      continue;
    }
    triangles[static_cast<usize>(survivingCount)] = triangle;
    mCubeID[static_cast<usize>(survivingCount)] = mCubeID[static_cast<usize>(i)];
    survivingCount++;
  }
  triangles.resize(static_cast<usize>(survivingCount));
  mCubeID.resize(static_cast<usize>(survivingCount));
  result.NumFacesRemoved = nTriangle - survivingCount;

  // Pass 3: a touched node is normally still referenced by the wall triangles on either side of the
  // edge, but if every triangle that used it was a chamfer (possible once the Bounding Box Skin prune
  // has removed the walls around it) it is now an orphan and must not be emitted.
  auto orphanCandidates = touched;
  for(const Triangle& triangle : triangles)
  {
    for(const SiteIdType nodeId : triangle.NodeId)
    {
      orphanCandidates.erase(nodeId);
    }
  }
  for(const SiteIdType orphan : orphanCandidates)
  {
    nodeType[static_cast<usize>(orphan)] = m3c_node_type::k_Unused;
    // Coordinate overrides also classify collapsed faces during streamed regeneration.
  }
  return result;
}

// -----------------------------------------------------------------------------
/**
 * @brief Fills pre-sized triangle arrays in cube order.
 * @param siteCoords Calculates padded-site coordinates.
 * @param triangles Receives triangle records.
 * @param mCubeID Receives the source cube for each triangle.
 * @param squares Provides marching-square records.
 * @param nodeCoords Calculates candidate-node coordinates.
 * @param faceEdges Provides mutable oriented face-edge records.
 * @param numSitesDim3 Specifies padded site count.
 * @param numSitesDim2 Specifies padded sites per Z plane.
 * @param xDim Specifies padded X dimension.
 * @param shouldCancel Stops before later cubes when true.
 */
void getTriangles(const SiteCoords& siteCoords, std::vector<Triangle>& triangles, std::vector<SiteIdType>& mCubeID, const std::vector<Face>& squares, const NodeCoords& nodeCoords,
                  std::vector<Segment>& faceEdges, const SiteIdType numSitesDim3, const SiteIdType numSitesDim2, const int xDim, const std::atomic_bool& shouldCancel)
{
  int64 tidIn = 0;
  int64 tidOut = 0;

  for(SiteIdType i = 1; i <= (numSitesDim3 - numSitesDim2); i++)
  {
    if(shouldCancel)
    {
      return;
    }
    int cubeFlag = 0;
    std::array<SiteIdType, 6> sqID{};
    sqID[0] = 3 * (i - 1);
    sqID[1] = (3 * (i - 1)) + 1;
    sqID[2] = (3 * (i - 1)) + 2;
    sqID[3] = (3 * i) + 2;
    sqID[4] = (3 * (i + xDim - 1)) + 1;
    sqID[5] = 3 * (i + numSitesDim2 - 1);
    int nFC = 0;
    int nFE = 0;
    int eff = 0;
    const SiteIdType bodyCtr = (7 * (i - 1)) + 6;
    std::array<SiteIdType, 6> arrayFC{};
    for(int ii = 0; ii < 6; ii++)
    {
      arrayFC[ii] = -1;
    }
    int fcid = 0;
    for(int ii = 0; ii < 6; ii++)
    {
      const SiteIdType tsq = sqID[ii];
      const SiteIdType tFCnode = squares[tsq].FaceCenterNode;
      if(tFCnode != -1)
      {
        arrayFC[fcid] = tFCnode;
        fcid++;
      }
      nFE = nFE + squares[tsq].NEdge;
      eff = eff + squares[tsq].Effect;
    }
    nFC = fcid;
    if(eff > 0)
    {
      cubeFlag = 1;
    }

    if(cubeFlag == 1 && nFE > 2)
    {
      std::array<double, 3> coord1{};
      std::array<double, 3> coord2{};
      for(int k = 0; k < 3; k++)
      {
        coord1[k] = siteCoords[i].Coord[k];
        coord2[k] = siteCoords[i + 1 + xDim + numSitesDim2].Coord[k];
      }
      std::vector<SiteIdType> arrayFE(nFE);
      int tindex = 0;
      for(int i1 = 0; i1 < 6; i1++)
      {
        const SiteIdType tsq = sqID[i1];
        const int tnfe = static_cast<int>(static_cast<uint8>(squares[tsq].NEdge));
        for(int i2 = 0; i2 < tnfe; i2++)
        {
          arrayFE[tindex] = squares[tsq].EdgeId[i2];
          tindex++;
        }
      }

      if(nFC == 0)
      {
        getCase0Triangles(triangles, mCubeID, arrayFE, nodeCoords, faceEdges, nFE, tidIn, tidOut, coord1, coord2, i);
        tidIn = tidOut;
      }
      else if(nFC == 2)
      {
        getCase2Triangles(triangles, mCubeID, arrayFE, nodeCoords, faceEdges, nFE, arrayFC, nFC, tidIn, tidOut, coord1, coord2, i);
        tidIn = tidOut;
      }
      else if(nFC > 2 && nFC <= 6)
      {
        getCaseMTriangles(triangles, mCubeID, arrayFE, nodeCoords, faceEdges, nFE, arrayFC, nFC, tidIn, tidOut, bodyCtr, coord1, coord2, i);
        tidIn = tidOut;
      }
    }
  }
}

/**
 * @brief Converts a padded site to an original cell index.
 * @param site Specifies the one-based padded site.
 * @param fileDim Specifies padded grid dimensions.
 * @param dims Specifies original grid dimensions.
 * @return Original zero-based cell index, or SIZE_MAX for a ghost site.
 */
usize paddedSiteToOriginalCell(const int64 site, const std::array<usize, 3>& fileDim, const std::array<usize, 3>& dims)
{
  const usize linear = static_cast<usize>(site - 1);
  const usize paddedX = linear % fileDim[0];
  const usize paddedY = (linear / fileDim[0]) % fileDim[1];
  const usize paddedZ = linear / (fileDim[0] * fileDim[1]);
  if(paddedX >= 1 && paddedX <= dims[0] && paddedY >= 1 && paddedY <= dims[1] && paddedZ >= 1 && paddedZ <= dims[2])
  {
    return ((paddedZ - 1) * dims[0] * dims[1]) + ((paddedY - 1) * dims[0]) + (paddedX - 1);
  }
  return std::numeric_limits<usize>::max();
}

/**
 * @brief Finds a non-ghost source cell for one working label.
 * @param workLabel Specifies the renumbered Feature Id.
 * @param cubeSite Specifies the cube origin site.
 * @param neighbors Provides cube neighbors.
 * @param featureIds Provides padded Feature Id values.
 * @param fileDim Specifies padded grid dimensions.
 * @param dims Specifies original grid dimensions.
 * @return Original cell index, or SIZE_MAX when the label is exterior.
 */
usize findSourceCell(const int workLabel, const int64 cubeSite, const NeighborAccessor& neighbors, const std::vector<int32>& featureIds, const std::array<usize, 3>& fileDim,
                     const std::array<usize, 3>& dims)
{
  const Neighbor siteNeighbors = neighbors[cubeSite]; // cache: 7 neighbors of the cube site read below
  const std::array<int64, 8> cornerSites = {cubeSite,
                                            siteNeighbors.NeighId[1],
                                            siteNeighbors.NeighId[7],
                                            siteNeighbors.NeighId[8],
                                            siteNeighbors.NeighId[18],
                                            siteNeighbors.NeighId[19],
                                            siteNeighbors.NeighId[25],
                                            siteNeighbors.NeighId[26]};
  for(const int64 site : cornerSites)
  {
    if(featureIds[site] == workLabel)
    {
      const usize original = paddedSiteToOriginalCell(site, fileDim, dims);
      if(original != std::numeric_limits<usize>::max())
      {
        return original;
      }
    }
  }
  return std::numeric_limits<usize>::max();
}

/**
 * @brief Finalizes mesh topology, transfers arrays, and repairs winding.
 * @param dataStructure Provides input and output objects.
 * @param inputValues Specifies output paths and options.
 * @param messageHandler Receives progress messages.
 * @param shouldCancel Stops later finalization stages when true.
 * @param triangles Provides generated triangle records.
 * @param mCubeID Provides triangle cube indexes.
 * @param fedges Provides face-edge scratch records.
 * @param nodeType Provides candidate node types.
 * @param featureIds Provides padded Feature Id values.
 * @param nodeCoords Calculates node coordinates.
 * @param neighbors Provides padded-grid neighbors.
 * @param numSites Specifies padded-grid site count.
 * @param fileDim Specifies padded grid dimensions.
 * @param dims Specifies original grid dimensions.
 * @param maxGrainId Specifies the reserved zero-label value.
 * @return Error during output or transfer, or success after cancellation.
 */
Result<> finalizeMesh(DataStructure& dataStructure, const M3CSurfaceMeshingInputValues* inputValues, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel,
                      std::vector<Triangle>& triangles, std::vector<SiteIdType>& mCubeID, std::vector<Segment>& fedges, std::vector<int8>& nodeType, std::vector<int32>& featureIds,
                      const NodeCoords& nodeCoords, const NeighborAccessor& neighbors, const SiteIdType numSites, const std::array<usize, 3>& fileDim, const std::array<usize, 3>& dims,
                      const int maxGrainId)
{
  const int64 nTriangle = static_cast<int64>(triangles.size());

  // Faces suppressed by the prune below. Stays 0 when the mode is Off, or when the mode is on but
  // nothing matched. Read after the winding-repair pass to decide which (if either) of the two
  // Bounding Box Skin warnings to emit.
  int64 numFacesPruned = 0;

  // Bounding Box Skin option, 'Background-Backed Walls Only' mode: drop faces whose output Face Labels would be {-1, 0}. In the
  // internal representation that is one negative ghost label paired with maxGrainId, which
  // is the renumbered zero-feature (see toFaceLabel below). Pruning the scratch vectors here
  // means the output TriangleGeom is sized from the surviving count and never over-allocated.
  if(inputValues->BoundingBoxSkinMode == BoundingBoxSkinMode::k_BackgroundBackedWallsOnly)
  {
    messageHandler.sendInfoMessage("Omitting bounding box skin faces...");

    // True when this triangle is a bounding-box wall face backed by background (i.e., its output
    // Face Labels would be {-1, 0}). M3C's single sequential pass over `triangles` (below) is the
    // only place this predicate is evaluated, so -- unlike QuickSurfaceMesh's SkipWallFace and
    // SurfaceNets' SkipPaddingQuad -- there is no second pass it must stay in agreement with.
    const auto skipBackgroundSkinFace = [maxGrainId](const Triangle& triangle) -> bool {
      const int spinA = triangle.NSpin[0];
      const int spinB = triangle.NSpin[1];
      return (spinA < 0 && spinB == maxGrainId) || (spinB < 0 && spinA == maxGrainId);
    };

    // Count how many triangles will be dropped before allocating droppedNodeIds: pushing onto an
    // unreserved vector here causes dozens of reallocations (and a transient ~1.5x peak) on a
    // representative dataset. This adds a second pass over `triangles`, but it only evaluates the
    // same boolean predicate the compaction loop below already does -- no allocation -- so it is
    // negligible next to the reallocations it avoids.
    int64 numToDrop = 0;
    for(int64 i = 0; i < nTriangle; i++)
    {
      if(skipBackgroundSkinFace(triangles[static_cast<usize>(i)]))
      {
        numToDrop++;
      }
    }

    int64 survivingCount = 0;
    // node_id values touched by a DROPPED triangle, recorded before the in-place compaction below
    // overwrites them. Used to narrow the nodeType clear (see below) to exactly the nodes the prune
    // itself orphaned, at a cost of O(3 * droppedCount) instead of a second full 7*numSites mask.
    std::vector<SiteIdType> droppedNodeIds;
    droppedNodeIds.reserve(static_cast<usize>(3 * numToDrop));
    for(int64 i = 0; i < nTriangle; i++)
    {
      const Triangle& triangle = triangles[static_cast<usize>(i)];
      if(skipBackgroundSkinFace(triangle))
      {
        droppedNodeIds.push_back(triangle.NodeId[0]);
        droppedNodeIds.push_back(triangle.NodeId[1]);
        droppedNodeIds.push_back(triangle.NodeId[2]);
      }
      else
      {
        triangles[static_cast<usize>(survivingCount)] = triangle;
        mCubeID[static_cast<usize>(survivingCount)] = mCubeID[static_cast<usize>(i)];
        survivingCount++;
      }
    }
    triangles.resize(static_cast<usize>(survivingCount));
    mCubeID.resize(static_cast<usize>(survivingCount));
    numFacesPruned = nTriangle - survivingCount;

    // Clear nodeType only for nodes the prune itself orphaned: referenced by a DROPPED triangle and
    // by no SURVIVING triangle. Nodes referenced only by survivors are left untouched, and
    // "pre-existing" candidates that no triangle -- dropped or surviving -- ever referenced are left
    // exactly as they were; the option no longer sweeps up orphan candidates it had no hand in creating.
    if(!droppedNodeIds.empty())
    {
      std::sort(droppedNodeIds.begin(), droppedNodeIds.end());
      droppedNodeIds.erase(std::unique(droppedNodeIds.begin(), droppedNodeIds.end()), droppedNodeIds.end());

      std::vector<bool> referencedBySurvivor(droppedNodeIds.size(), false);
      for(const auto& triangle : triangles)
      {
        for(const SiteIdType nodeId : triangle.NodeId)
        {
          const auto droppedNodeIter = std::lower_bound(droppedNodeIds.begin(), droppedNodeIds.end(), nodeId);
          if(droppedNodeIter != droppedNodeIds.end() && *droppedNodeIter == nodeId)
          {
            referencedBySurvivor[static_cast<usize>(droppedNodeIter - droppedNodeIds.begin())] = true;
          }
        }
      }
      for(usize i = 0; i < droppedNodeIds.size(); i++)
      {
        if(!referencedBySurvivor[i])
        {
          nodeType[static_cast<usize>(droppedNodeIds[i])] = m3c_node_type::k_Unused;
        }
      }
    }
  }

  // Promote surface nodes to their exterior variant (+10). A triangle that borders the outside of the
  // volume has exactly one negative feature label (nSpin[0]*nSpin[1] < 0), so each of its nodes lies on
  // the volume boundary. This is the only output-relevant effect of the legacy triangle-side/inner-edge
  // connectivity pass: the per-triangle edge ids, edgePlace flags, and unique inner-edge list it also
  // built never appear in the output (Triangle Geometry + Face Labels + Node Types), so that machinery
  // has been removed.
  for(usize j = 0; j < triangles.size(); j++)
  {
    if(triangles[j].NSpin[0] * triangles[j].NSpin[1] < 0)
    {
      for(int i = 0; i < 3; i++)
      {
        const SiteIdType nodeId = triangles[j].NodeId[i];
        if(nodeType[nodeId] < 10)
        {
          nodeType[nodeId] = static_cast<int8>(nodeType[nodeId] + 10);
        }
      }
    }
  }

  // Sharp Bounding Box Edges: snap the outermost wall rows onto the box edges and drop the chamfer
  // triangles (see sharpenBoundingBoxEdges). Runs on the scratch vectors, so the output TriangleGeom is
  // sized from the surviving count exactly as for the skin prune above.
  SharpEdgeResult sharpEdges;
  if(inputValues->SharpBoundingBoxEdges)
  {
    messageHandler.sendInfoMessage("Sharpening bounding box edges...");
    sharpEdges = sharpenBoundingBoxEdges(triangles, mCubeID, nodeType, 7 * numSites, nodeCoords, dims);
    messageHandler.sendInfoMessage(fmt::format("Sharpened bounding box edges: removed {} chamfer triangles", sharpEdges.NumFacesRemoved));
  }

  const int64 nTriangleFinal = static_cast<int64>(triangles.size());

  // The face-edge segments are no longer needed; release before the memory-heavy output + winding stages.
  std::vector<Segment>().swap(fedges);

  if(shouldCancel)
  {
    return {};
  }

  messageHandler.sendInfoMessage("Writing surface mesh...");
  // Node-id compaction without a dense 7*numSites candidate->id map. A candidate's compacted id is simply
  // the number of real nodes (nodeType > 0) that precede it; we answer that from a coarse per-block
  // prefix over nodeType plus a small in-block scan (saves ~3.8 GB at 512^3 vs a uint32 map). This is
  // valid because the prefix is built here, after the skin prune and the sharp-edge pass have cleared
  // the nodes they retire and the surface-node promotion has added +10 to the rest.
  const SiteIdType numCandidateNodes = 7 * numSites;
  constexpr SiteIdType nodeBlock = 128;
  const SiteIdType numNodeBlocks = (numCandidateNodes + nodeBlock - 1) / nodeBlock;
  std::vector<uint32> nodeBlockBase(static_cast<usize>(numNodeBlocks));
  int64 realNodeRunning = 0;
  for(SiteIdType blockIdx = 0; blockIdx < numNodeBlocks; blockIdx++)
  {
    nodeBlockBase[static_cast<usize>(blockIdx)] = static_cast<uint32>(realNodeRunning);
    const SiteIdType lowIndex = blockIdx * nodeBlock;
    const SiteIdType highIndex = std::min<SiteIdType>(lowIndex + nodeBlock, numCandidateNodes);
    for(SiteIdType candidateId = lowIndex; candidateId < highIndex; candidateId++)
    {
      if(nodeType[candidateId] > 0)
      {
        realNodeRunning++;
      }
    }
  }
  const int64 nNodes = realNodeRunning;
  const auto compactedNodeId = [&](const SiteIdType candidateId) -> int64 {
    int64 compactId = nodeBlockBase[static_cast<usize>(candidateId / nodeBlock)];
    for(SiteIdType cc = (candidateId / nodeBlock) * nodeBlock; cc < candidateId; cc++)
    {
      if(nodeType[cc] > 0)
      {
        compactId++;
      }
    }
    return compactId;
  };

  auto& triangleGeom = dataStructure.getDataRefAs<TriangleGeom>(inputValues->TriangleGeometryPath);
  Result<> resizeResult = triangleGeom.resizeVertexList(static_cast<usize>(nNodes));
  if(resizeResult.invalid())
  {
    return resizeResult;
  }
  resizeResult = triangleGeom.resizeFaceList(static_cast<usize>(nTriangleFinal));
  if(resizeResult.invalid())
  {
    return resizeResult;
  }
  resizeResult = triangleGeom.getVertexAttributeMatrix()->resizeTuples({static_cast<usize>(nNodes)});
  if(resizeResult.invalid())
  {
    return resizeResult;
  }
  resizeResult = triangleGeom.getFaceAttributeMatrix()->resizeTuples({static_cast<usize>(nTriangleFinal)});
  if(resizeResult.invalid())
  {
    return resizeResult;
  }

  auto& vertexStore = triangleGeom.getVertices()->getDataStoreRef();
  auto& triStore = triangleGeom.getFaces()->getDataStoreRef();
  auto& faceLabels = dataStructure.getDataRefAs<Int32Array>(inputValues->FaceLabelsDataPath).getDataStoreRef();
  auto& nodeTypesOut = dataStructure.getDataRefAs<Int8Array>(inputValues->NodeTypesDataPath).getDataStoreRef();
  resizeResult = faceLabels.resizeTuples({static_cast<usize>(nTriangleFinal)});
  if(resizeResult.invalid())
  {
    return resizeResult;
  }
  resizeResult = nodeTypesOut.resizeTuples({static_cast<usize>(nNodes)});
  if(resizeResult.invalid())
  {
    return resizeResult;
  }

  // Emit real candidates in ascending order. This order preserves the legacy
  // compact node numbering without a candidate-to-node map.
  int64 vtxRunning = 0;
  for(SiteIdType i = 0; i < numCandidateNodes; i++)
  {
    if(nodeType[i] > 0)
    {
      const auto snappedIt = sharpEdges.SnappedCoords.find(i);
      const Node nodeCoord = (snappedIt != sharpEdges.SnappedCoords.end()) ? snappedIt->second : nodeCoords[i];
      vertexStore[(static_cast<usize>(vtxRunning) * 3) + 0] = nodeCoord.Coord[0];
      vertexStore[(static_cast<usize>(vtxRunning) * 3) + 1] = nodeCoord.Coord[1];
      vertexStore[(static_cast<usize>(vtxRunning) * 3) + 2] = nodeCoord.Coord[2];
      nodeTypesOut[static_cast<usize>(vtxRunning)] = nodeType[i];
      vtxRunning++;
    }
  }

  // FaceLabels matches QuickSurfaceMesh and SurfaceNets. Negative ghost labels
  // become -1, and the reserved zero label becomes 0. The smaller label is first
  // because downstream filters require this order. Winding repair uses the same order.
  const auto toFaceLabel = [maxGrainId](const int nSpin) -> int32 { return (nSpin < 0) ? -1 : ((nSpin == maxGrainId) ? 0 : nSpin); };

  // Triangles: remap to compacted node ids and write the ordered FaceLabels.
  for(int64 i = 0; i < nTriangleFinal; i++)
  {
    triStore[(static_cast<usize>(i) * 3) + 0] = static_cast<IGeometry::MeshIndexType>(compactedNodeId(triangles[i].NodeId[0]));
    triStore[(static_cast<usize>(i) * 3) + 1] = static_cast<IGeometry::MeshIndexType>(compactedNodeId(triangles[i].NodeId[1]));
    triStore[(static_cast<usize>(i) * 3) + 2] = static_cast<IGeometry::MeshIndexType>(compactedNodeId(triangles[i].NodeId[2]));

    const int32 labelA = toFaceLabel(triangles[i].NSpin[0]);
    const int32 labelB = toFaceLabel(triangles[i].NSpin[1]);
    faceLabels[(static_cast<usize>(i) * 2) + 0] = (labelA <= labelB) ? labelA : labelB;
    faceLabels[(static_cast<usize>(i) * 2) + 1] = (labelA <= labelB) ? labelB : labelA;
  }

  // Transfer selected arrays to both face sides. Each side uses a source cell
  // whose working label matches that side. TupleTransfer skips exterior sides.
  if(!inputValues->SelectedCellDataArrayPaths.empty() || !inputValues->SelectedFeatureDataArrayPaths.empty())
  {
    messageHandler.sendInfoMessage("Transferring attribute arrays to the mesh faces...");
    std::vector<std::shared_ptr<AbstractTupleTransfer>> transfers;
    for(usize i = 0; i < inputValues->SelectedCellDataArrayPaths.size(); i++)
    {
      AddTupleTransferInstance(dataStructure, inputValues->SelectedCellDataArrayPaths[i], inputValues->CreatedDataArrayPaths[i], transfers);
    }
    const usize numCellArrays = inputValues->SelectedCellDataArrayPaths.size();
    for(usize i = 0; i < inputValues->SelectedFeatureDataArrayPaths.size(); i++)
    {
      AddFeatureTupleTransferInstance(dataStructure, inputValues->SelectedFeatureDataArrayPaths[i], inputValues->CreatedDataArrayPaths[numCellArrays + i], inputValues->FeatureIdsArrayPath, transfers);
    }

    for(int64 i = 0; i < nTriangleFinal; i++)
    {
      // Use the FaceLabels order so each transferred component aligns with its label.
      const int32 labelA = toFaceLabel(triangles[i].NSpin[0]);
      const int32 labelB = toFaceLabel(triangles[i].NSpin[1]);
      const bool side0IsComp0 = (labelA <= labelB);
      const int nSpinComp0 = side0IsComp0 ? triangles[i].NSpin[0] : triangles[i].NSpin[1];
      const int nSpinComp1 = side0IsComp0 ? triangles[i].NSpin[1] : triangles[i].NSpin[0];
      const usize cell0 = findSourceCell(nSpinComp0, mCubeID[i], neighbors, featureIds, fileDim, dims);
      const usize cell1 = findSourceCell(nSpinComp1, mCubeID[i], neighbors, featureIds, fileDim, dims);
      for(const auto& transfer : transfers)
      {
        transfer->quickSurfaceTransfer(static_cast<usize>(i), cell0, cell1, faceLabels);
      }
    }
  }

  // Winding repair reads only the output geometry and FaceLabels. Release working
  // buffers before adjacency allocation to reduce peak memory.
  std::vector<Triangle>().swap(triangles);
  std::vector<SiteIdType>().swap(mCubeID);
  std::vector<int8>().swap(nodeType);
  std::vector<int32>().swap(featureIds);
  std::vector<uint32>().swap(nodeBlockBase);

  // M3C does not guarantee globally consistent normals. Optional repair uses
  // triangle connectivity to make winding consistent with FaceLabels.
  if(inputValues->RepairTriangleWinding)
  {
    messageHandler.sendInfoMessage("Generating connectivity and triangle neighbors...");
    triangleGeom.findElementNeighbors(true);
    const auto optionalId = triangleGeom.getElementNeighborsId();
    if(optionalId.has_value())
    {
      const auto& connectivity = dataStructure.getDataRefAs<IGeometry::ElementDynamicList>(optionalId.value());
      messageHandler.sendInfoMessage("Repairing windings...");
      Result<> windingResult = MeshingUtilities::RepairTriangleWinding(triangleGeom.getFaces()->getDataStoreRef(), connectivity,
                                                                       dataStructure.getDataAs<Int32Array>(inputValues->FaceLabelsDataPath)->getDataStoreRef(), shouldCancel, messageHandler);
      const auto containingVertId = triangleGeom.getElementContainingVertId();
      if(containingVertId.has_value())
      {
        dataStructure.removeData(containingVertId.value());
      }
      dataStructure.removeData(optionalId.value());
      if(windingResult.invalid())
      {
        return windingResult;
      }
    }
  }

  if(inputValues->BoundingBoxSkinMode == BoundingBoxSkinMode::k_BackgroundBackedWallsOnly)
  {
    // An entirely-background volume has nothing but {-1, 0} faces, so omitting the skin
    // legitimately produces an empty mesh. Report it rather than returning silently. nNodes is
    // expected to be zero here as well (every node was orphaned by the prune and cleared); it is
    // passed through so the warning stays honest if that ever changes.
    if(nTriangleFinal == 0)
    {
      return MeshingUtilities::MakeEmptyMeshWarning(inputValues->TriangleGeometryPath, dataStructure.getDataRefAs<Int32Array>(inputValues->FeatureIdsArrayPath).getNumberOfTuples(),
                                                    static_cast<usize>(nNodes));
    }
    // A fully-indexed volume (no Feature Id 0) makes the option a no-op: nothing was pruned, and
    // the user otherwise gets byte-identical output with no feedback that the option had no effect.
    if(numFacesPruned == 0)
    {
      return MeshingUtilities::MakeNoFacesPrunedWarning(inputValues->TriangleGeometryPath);
    }
  }

  return {};
}
} // namespace

namespace nx::core
{
M3CSurfaceMeshing::M3CSurfaceMeshing(DataStructure& dataStructure, M3CSurfaceMeshingInputValues* inputValues, const std::atomic_bool& shouldCancel, const IFilter::MessageHandler& mesgHandler)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(mesgHandler)
{
}

M3CSurfaceMeshing::~M3CSurfaceMeshing() noexcept = default;

Result<> M3CSurfaceMeshing::operator()()
{
  // Every dynamic cell and mesh array selects the residency path. A selected
  // input or created output can be disk-backed while Feature Ids remain resident.
  std::vector<const IArray*> dispatchTargets;
  const auto appendArray = [this, &dispatchTargets](const DataPath& path) {
    if(const auto* arrayPtr = m_DataStructure.getDataAs<IDataArray>(path); arrayPtr != nullptr)
    {
      dispatchTargets.push_back(arrayPtr);
    }
  };
  appendArray(m_InputValues->FeatureIdsArrayPath);
  appendArray(m_InputValues->NodeTypesDataPath);
  appendArray(m_InputValues->FaceLabelsDataPath);
  for(const auto& path : m_InputValues->SelectedCellDataArrayPaths)
  {
    appendArray(path);
  }
  for(const auto& path : m_InputValues->CreatedDataArrayPaths)
  {
    appendArray(path);
  }
  const auto& triangleGeom = m_DataStructure.getDataRefAs<TriangleGeom>(m_InputValues->TriangleGeometryPath);
  dispatchTargets.push_back(triangleGeom.getVertices());
  dispatchTargets.push_back(triangleGeom.getFaces());

  const bool usesOutOfCoreStore = AnyOutOfCore(AlgorithmArrayTargets(dispatchTargets));
  const bool useOutOfCorePath = !ForceInCoreAlgorithm() && (usesOutOfCoreStore || ForceOocAlgorithm());
  RecordAlgorithmPathExecution(useOutOfCorePath ? AlgorithmPath::OutOfCore : AlgorithmPath::InCore, usesOutOfCoreStore);

  const auto& featureIdsStore = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->FeatureIdsArrayPath).getDataStoreRef();
  Result<> sentinelCheck = MeshingUtilities::ValidateFeatureIdsAgainstSentinels(featureIdsStore, m_InputValues->FeatureIdsArrayPath, true, m_ShouldCancel, m_MessageHandler);
  if(sentinelCheck.invalid())
  {
    return sentinelCheck;
  }

  if(useOutOfCorePath)
  {
    return runOutOfCore(dispatchTargets, usesOutOfCoreStore);
  }

  // Default: the multithreaded sliding-window sweep (runWindowed(parallel=true)). Peak per-site scratch
  // is O(sliceArea) instead of O(volume), and the per-cube work runs across all cores. It is watertight
  // and correct, with byte-identical vertices, FaceLabels, and NodeTypes to the serial path, but a
  // slightly different (still valid) triangulation of the same interfaces -- the legacy per-cube loop
  // triangulation depends on cross-cube edge-flip propagation, which is inherently serial. The parallel
  // output is deterministic (each cube depends only on its own inputs, independent of thread scheduling).
  //
  // Two serial reference paths are kept for validation/debugging, selected via environment variables:
  //   M3C_SERIAL=1        -> runWindowed(false): serial sliding window (same tessellation as legacy)
  //   M3C_WHOLE_VOLUME=1  -> runEntireVolume():  serial whole-volume (O(volume) memory)
  // Both serial paths are byte-identical to each other.

  if(const char* wholeVolPtr = std::getenv("M3C_WHOLE_VOLUME"); wholeVolPtr != nullptr && std::string_view(wholeVolPtr) == "1")
  {
    return runEntireVolume();
  }
  if(const char* serialPtr = std::getenv("M3C_SERIAL"); serialPtr != nullptr && std::string_view(serialPtr) == "1")
  {
    return runWindowed(false);
  }
  return runWindowed(true);
}

Result<> M3CSurfaceMeshing::runOutOfCore(const std::vector<const IArray*>& dispatchTargets, const bool usesOutOfCoreStore)
{
  if(dispatchTargets.empty())
  {
    return MakeErrorResult(-90544, "M3C out-of-core dispatch did not receive any dynamic storage targets.");
  }

  const auto& imageGeom = m_DataStructure.getDataRefAs<ImageGeom>(m_InputValues->GridGeomDataPath);
  const auto& featureIds = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->FeatureIdsArrayPath);
  const auto& featureIdsStore = featureIds.getDataStoreRef();
  const SizeVec3 gridDims = imageGeom.getDimensions();
  const std::array<usize, 3> dims = {gridDims[0], gridDims[1], gridDims[2]};
  if(dims[0] == 0 || dims[1] == 0 || dims[2] == 0 || dims[0] > std::numeric_limits<usize>::max() / dims[1] || dims[0] * dims[1] > std::numeric_limits<usize>::max() / dims[2])
  {
    return MakeErrorResult(-90546, "M3C out-of-core input dimensions are zero or overflow the cell count.");
  }
  const usize cellCount = dims[0] * dims[1] * dims[2];
  if(featureIdsStore.getNumberOfTuples() != cellCount)
  {
    return MakeErrorResult(-90547, "M3C out-of-core FeatureIds tuple count does not match the Image Geometry.");
  }

  // First bounded pass preserves initializeMicro's zero-feature renumbering
  // without retaining a second copy of the cell data. Large bounded batches
  // are heap-backed so this path remains within the default Windows stack.
  constexpr usize kFeatureIdBulkValues = 65536;
  std::vector<int32> maxScanBuffer(kFeatureIdBulkValues);
  int32 maxGrainId = 0;
  for(usize offset = 0; offset < cellCount;)
  {
    if(m_ShouldCancel)
    {
      return {};
    }
    const usize count = std::min(kFeatureIdBulkValues, cellCount - offset);
    auto readResult = featureIdsStore.copyIntoBuffer(offset, nonstd::span<int32>(maxScanBuffer.data(), count));
    if(readResult.invalid())
    {
      return readResult;
    }
    for(usize index = 0; index < count; index++)
    {
      maxGrainId = std::max(maxGrainId, maxScanBuffer[index]);
    }
    offset += count;
  }
  if(maxGrainId == std::numeric_limits<int32>::max())
  {
    return MakeErrorResult(-90548, "M3C out-of-core FeatureIds maximum cannot be incremented to reserve the zero feature.");
  }
  maxGrainId++;

  const std::array<usize, 3> fileDim = {dims[0] + 2, dims[1] + 2, dims[2] + 2};
  if(fileDim[0] < dims[0] || fileDim[1] < dims[1] || fileDim[2] < dims[2] || fileDim[0] > std::numeric_limits<usize>::max() / fileDim[1] ||
     fileDim[0] * fileDim[1] > std::numeric_limits<usize>::max() / fileDim[2])
  {
    return MakeErrorResult(-90549, "M3C out-of-core padded dimensions overflow.");
  }
  if(fileDim[0] > static_cast<usize>(std::numeric_limits<int>::max()) || fileDim[0] * fileDim[1] > static_cast<usize>(std::numeric_limits<SiteIdType>::max()) ||
     fileDim[0] * fileDim[1] * fileDim[2] > static_cast<usize>(std::numeric_limits<SiteIdType>::max()))
  {
    return MakeErrorResult(-90550, "M3C out-of-core padded dimensions cannot be represented by its signed site/index arithmetic.");
  }
  const usize paddedSiteCount = fileDim[0] * fileDim[1] * fileDim[2];
  const SiteIdType numSites = static_cast<SiteIdType>(paddedSiteCount);
  const usize paddedSitesPerPlane = fileDim[0] * fileDim[1];
  const SiteIdType numSitesPerPlane = static_cast<SiteIdType>(paddedSitesPerPlane);
  if(numSites > std::numeric_limits<SiteIdType>::max() / 7)
  {
    return MakeErrorResult(-90551, "M3C out-of-core square or candidate-node count overflows its site index type.");
  }

  // Four LRU Z slices cover local cube and slice-plane anomaly lookups while keeping the
  // padded ghost shell implicit. The cache also handles the NeighborAccessor's
  // toroidal border indices without materializing a padded volume.
  const usize sourceSliceSize = dims[0] * dims[1];
  std::array<std::vector<int32>, 4> sourceSlices{};
  std::array<int64, 4> sourceSliceZ = {-1, -1, -1, -1};
  std::array<uint64, 4> sourceSliceUse{};
  uint64 sourceUseCounter = 0;
  const auto sourceValue = [&](const SiteIdType site) -> Result<int32> {
    const usize linear = static_cast<usize>(site - 1);
    const usize xIndex = linear % fileDim[0];
    const usize yIndex = (linear / fileDim[0]) % fileDim[1];
    const usize zIndex = linear / (fileDim[0] * fileDim[1]);
    if(zIndex == 0 || zIndex + 1 == fileDim[2] || yIndex == 0 || yIndex + 1 == fileDim[1] || xIndex == 0 || xIndex + 1 == fileDim[0])
    {
      return {k_GhostLabel};
    }
    const int64 sourceZ = static_cast<int64>(zIndex - 1);
    usize slot = 0;
    while(slot < sourceSlices.size() && sourceSliceZ[slot] != sourceZ)
    {
      slot++;
    }
    if(slot == sourceSlices.size())
    {
      slot = static_cast<usize>(std::min_element(sourceSliceUse.begin(), sourceSliceUse.end()) - sourceSliceUse.begin());
      try
      {
        sourceSlices[slot].resize(sourceSliceSize);
      } catch(const std::bad_alloc&)
      {
        return MakeErrorResult<int32>(-90560, "M3C out-of-core rolling FeatureIds slice allocation failed.");
      }
      const usize sourceOffset = static_cast<usize>(sourceZ) * sourceSliceSize;
      auto readResult = featureIdsStore.copyIntoBuffer(sourceOffset, nonstd::span<int32>(sourceSlices[slot].data(), sourceSliceSize));
      if(readResult.invalid())
      {
        return ConvertInvalidResult<int32>(std::move(readResult));
      }
      sourceSliceZ[slot] = sourceZ;
    }
    sourceSliceUse[slot] = ++sourceUseCounter;
    const int32 value = sourceSlices[slot][((yIndex - 1) * dims[0]) + (xIndex - 1)];
    return {value == 0 ? maxGrainId : value};
  };

  const NeighborAccessor neighbors{numSites, numSitesPerPlane, static_cast<int>(fileDim[0])};
  const FloatVec3 spacing = imageGeom.getSpacing();
  const FloatVec3 origin = imageGeom.getOrigin();
  const SiteCoords siteCoords{fileDim[0], fileDim[1], fileDim[0] * fileDim[1], {spacing[0], spacing[1], spacing[2]}, {origin[0], origin[1], origin[2]}};
  const NodeCoords nodeCoords{siteCoords};
  const uint64 candidateCount = static_cast<uint64>(7 * numSites);
  const SiteIdType lastCube = numSites - numSitesPerPlane;
  if(lastCube < 0 || static_cast<uint64>(lastCube) == std::numeric_limits<uint64>::max())
  {
    return MakeErrorResult(-90552, "M3C out-of-core cube-count record range overflows.");
  }
  const uint64 cubeRecordCount = static_cast<uint64>(lastCube) + 1;
  auto candidateResult = TemporaryRecordVector<M3CCandidateNodeRecord>::create(candidateCount, usesOutOfCoreStore, m_ShouldCancel);
  if(candidateResult.invalid())
  {
    return ConvertResult(std::move(candidateResult));
  }
  auto triangleCountResult = TemporaryRecordVector<int64>::create(cubeRecordCount, usesOutOfCoreStore, m_ShouldCancel);
  if(triangleCountResult.invalid())
  {
    return ConvertResult(std::move(triangleCountResult));
  }
  auto candidateNodes = std::move(candidateResult.value());
  auto triangleCounts = std::move(triangleCountResult.value());
  const M3CCandidateNodeRecord unusedNode{};
  auto fillNodesResult = candidateNodes.store().fill(0, candidateCount, nonstd::span<const std::byte>(reinterpret_cast<const std::byte*>(&unusedNode), sizeof(unusedNode)), m_ShouldCancel);
  if(fillNodesResult.invalid())
  {
    return fillNodesResult;
  }
  const int64 zeroTriangleCount = 0;
  auto fillCountsResult =
      triangleCounts.store().fill(0, cubeRecordCount, nonstd::span<const std::byte>(reinterpret_cast<const std::byte*>(&zeroTriangleCount), sizeof(zeroTriangleCount)), m_ShouldCancel);
  if(fillCountsResult.invalid())
  {
    return fillCountsResult;
  }

  const auto setNodeType = [&](const SiteIdType nodeId, const int8 type) -> Result<> {
    auto nodeResult = candidateNodes.cache().read(static_cast<uint64>(nodeId), m_ShouldCancel);
    if(nodeResult.invalid())
    {
      return ConvertResult(std::move(nodeResult));
    }
    auto node = nodeResult.value();
    node.Type = type;
    return candidateNodes.cache().write(static_cast<uint64>(nodeId), node, m_ShouldCancel);
  };

  // Reconstruct one marching square in fixed local storage. The edge ids are
  // local to the caller; only the candidate-node classification survives pass
  // one and it lives in the external record vector.
  const auto buildSquare = [&](const SiteIdType squareId, Face& square, std::array<Segment, 64>& segments, int& segmentCount, const bool writeNodeTypes) -> Result<> {
    square = {};
    for(auto& edge : square.EdgeId)
    {
      edge = k_UnusedNodeId;
    }
    square.FaceCenterNode = -1;
    const SiteIdType cubeOrigin = (squareId / 3) + 1;
    const int squareOrder = static_cast<int>(squareId % 3);
    const auto corners = squareCorners(squareId, neighbors);
    std::array<int, 4> spins{};
    int ghostCorners = 0;
    for(int index = 0; index < 4; index++)
    {
      auto spinResult = sourceValue(corners[index]);
      if(spinResult.invalid())
      {
        return ConvertResult(std::move(spinResult));
      }
      spins[index] = spinResult.value();
      ghostCorners += spins[index] < 0 ? 1 : 0;
    }
    if(ghostCorners != 4)
    {
      square.Effect = 1;
    }
    if(ghostCorners == 4)
    {
      return {};
    }
    int squareIndex = getSquareIndex(spins);
    if(squareIndex == 15)
    {
      std::array<int, 4> neighborCounts = {0, 0, 0, 0};
      for(int corner = 0; corner < 4; corner++)
      {
        const Neighbor cornerNeighbors = neighbors[corners[corner]];
        // Match the eight slice-plane neighbors used by treatAnomaly().
        for(int neighborIndex = 1; neighborIndex <= 8; neighborIndex++)
        {
          auto neighborSpin = sourceValue(cornerNeighbors.NeighId[neighborIndex]);
          if(neighborSpin.invalid())
          {
            return ConvertResult(std::move(neighborSpin));
          }
          neighborCounts[corner] += spins[corner] == neighborSpin.value() && neighborSpin.value() > 0 ? 1 : 0;
        }
      }
      int minimum = 1000;
      int minimumIndex = -1;
      for(int corner = 0; corner < 4; corner++)
      {
        if(neighborCounts[corner] < minimum)
        {
          minimum = neighborCounts[corner];
          minimumIndex = corner;
        }
      }
      squareIndex += minimumIndex == 0 || minimumIndex == 2 ? 1 : 0;
    }
    if(squareIndex == 0)
    {
      return {};
    }
    for(int edgeIndex = 0; edgeIndex < 8; edgeIndex += 2)
    {
      if(k_EdgeTable2d[squareIndex][edgeIndex] == -1)
      {
        continue;
      }
      const std::array<int, 2> nodeIndex = {k_EdgeTable2d[squareIndex][edgeIndex], k_EdgeTable2d[squareIndex][edgeIndex + 1]};
      const std::array<int, 2> pixelIndex = {k_NsTable2d[squareIndex][edgeIndex], k_NsTable2d[squareIndex][edgeIndex + 1]};
      std::array<SiteIdType, 2> nodeIds{};
      getNodes(cubeOrigin, squareOrder, nodeIndex, nodeIds, numSitesPerPlane, static_cast<int>(fileDim[0]));
      const std::array<int, 2> pixelSpins = {spins[pixelIndex[0]], spins[pixelIndex[1]]};
      if(pixelSpins[0] > 0 || pixelSpins[1] > 0)
      {
        if(segmentCount >= static_cast<int>(segments.size()))
        {
          return MakeErrorResult(-90561, "M3C out-of-core local square edge buffer overflowed.");
        }
        segments[static_cast<usize>(segmentCount)] = Segment{{nodeIds[0], nodeIds[1]}, {pixelSpins[0], pixelSpins[1]}};
        square.EdgeId[square.NEdge++] = static_cast<uint32>(segmentCount++);
      }
      else if(writeNodeTypes)
      {
        auto firstResult = setNodeType(nodeIds[0], m3c_node_type::k_Unused);
        if(firstResult.invalid())
        {
          return firstResult;
        }
        auto secondResult = setNodeType(nodeIds[1], m3c_node_type::k_Unused);
        if(secondResult.invalid())
        {
          return secondResult;
        }
      }
      for(int node = 0; node < 2; node++)
      {
        if(nodeIndex[node] == 4 && (squareIndex == 7 || squareIndex == 11 || squareIndex == 13 || squareIndex == 14 || squareIndex == 19))
        {
          square.FaceCenterNode = nodeIds[node];
        }
        if(writeNodeTypes)
        {
          int8 type = m3c_node_type::k_Default;
          if(nodeIndex[node] == 4)
          {
            if(squareIndex == 19)
            {
              type = m3c_node_type::k_QuadPoint;
            }
            else if(squareIndex == 7 || squareIndex == 11 || squareIndex == 13 || squareIndex == 14)
            {
              type = m3c_node_type::k_TriplePoint;
            }
            else
            {
              continue;
            }
          }
          auto typeResult = setNodeType(nodeIds[node], type);
          if(typeResult.invalid())
          {
            return typeResult;
          }
        }
      }
    }
    return {};
  };

  // Match the legacy edge-stage visitation order when classifying
  // candidate nodes. No resident square/edge vector survives this pass.
  for(SiteIdType squareId = 0; squareId < 3 * numSites; squareId++)
  {
    if(m_ShouldCancel)
    {
      return {};
    }
    Face square{};
    std::array<Segment, 64> segments{};
    int segmentCount = 0;
    auto squareResult = buildSquare(squareId, square, segments, segmentCount, true);
    if(squareResult.invalid())
    {
      return squareResult;
    }
  }

  // A cube references only six squares and at most 24 face edges, so its
  // count is rebuilt locally. Count records preserve the default parallel
  // path's cube indexing/order for the later prefix + generation pass.
  uint64 numFacesPruned = 0;
  uint64 maximumTrianglesPerCube = 0;
  const auto skipBackgroundSkinFace = [maxGrainId](const Triangle& triangle) {
    const int spinA = triangle.NSpin[0];
    const int spinB = triangle.NSpin[1];
    return (spinA < 0 && spinB == maxGrainId) || (spinB < 0 && spinA == maxGrainId);
  };
  std::vector<Triangle> edgeTriangles;
  std::vector<SiteIdType> edgeCubes;
  const HalfCellLattice lattice{nodeCoords};
  const auto cubeTouchesBoxEdge = [&](const SiteIdType cube) {
    const usize linear = static_cast<usize>(cube - 1);
    const std::array<usize, 3> position = {linear % fileDim[0], (linear / fileDim[0]) % fileDim[1], linear / (fileDim[0] * fileDim[1])};
    int wallAxes = 0;
    int nearWallAxes = 0;
    for(usize axis = 0; axis < 3; axis++)
    {
      wallAxes += (position[axis] == 0 || position[axis] == dims[axis]) ? 1 : 0;
      nearWallAxes += (position[axis] <= 1 || position[axis] >= dims[axis] - 1) ? 1 : 0;
    }
    return wallAxes > 0 && nearWallAxes >= 2;
  };
  const auto touchesBoxEdge = [&](const Triangle& triangle) {
    for(const SiteIdType nodeId : triangle.NodeId)
    {
      const auto position = lattice(nodeId);
      int wallAxes = 0;
      int nearWallAxes = 0;
      for(usize axis = 0; axis < 3; axis++)
      {
        const int64 upper = 2 * static_cast<int64>(dims[axis]);
        const bool onWall = position[axis] == 0 || position[axis] == upper;
        wallAxes += onWall ? 1 : 0;
        nearWallAxes += (onWall || (dims[axis] >= 2 && (position[axis] == 1 || position[axis] == upper - 1))) ? 1 : 0;
      }
      if(wallAxes > 0 && nearWallAxes >= 2)
      {
        return true;
      }
    }
    return false;
  };
  for(SiteIdType cube = 1; cube <= lastCube; cube++)
  {
    if(m_ShouldCancel)
    {
      return {};
    }
    const std::array<SiteIdType, 6> squareIds = {
        3 * (cube - 1), (3 * (cube - 1)) + 1, (3 * (cube - 1)) + 2, (3 * cube) + 2, (3 * (cube + static_cast<SiteIdType>(fileDim[0]) - 1)) + 1, 3 * (cube + numSitesPerPlane - 1)};
    std::array<Face, 6> squares{};
    std::array<Segment, 64> segments{};
    int segmentCount = 0;
    std::array<SiteIdType, 6> faceCenters = {-1, -1, -1, -1, -1, -1};
    int faceCenterCount = 0;
    int edgeCount = 0;
    int effectiveCount = 0;
    for(int square = 0; square < 6; square++)
    {
      auto squareResult = buildSquare(squareIds[square], squares[square], segments, segmentCount, false);
      if(squareResult.invalid())
      {
        return squareResult;
      }
      if(squares[square].FaceCenterNode != -1)
      {
        faceCenters[faceCenterCount++] = squares[square].FaceCenterNode;
      }
      edgeCount += squares[square].NEdge;
      effectiveCount += squares[square].Effect;
    }
    if(faceCenterCount >= 3)
    {
      const auto firstCorners = squareCorners(squareIds[0], neighbors);
      const auto lastCorners = squareCorners(squareIds[5], neighbors);
      int uniqueSpins = 0;
      std::array<int, 8> cubeSpins{};
      for(int index = 0; index < 4; index++)
      {
        auto firstSpin = sourceValue(firstCorners[index]);
        auto lastSpin = sourceValue(lastCorners[index]);
        if(firstSpin.invalid() || lastSpin.invalid())
        {
          return firstSpin.invalid() ? ConvertResult(std::move(firstSpin)) : ConvertResult(std::move(lastSpin));
        }
        cubeSpins[index] = firstSpin.value();
        cubeSpins[index + 4] = lastSpin.value();
      }
      for(int index = 0; index < 8; index++)
      {
        const int spin = cubeSpins[index];
        if(spin != -1)
        {
          uniqueSpins++;
          cubeSpins[index] = -1;
          for(int other = 0; other < 8; other++)
          {
            if(cubeSpins[other] == spin)
            {
              cubeSpins[other] = -1;
            }
          }
        }
      }
      auto bodyResult = setNodeType((7 * (cube - 1)) + 6, static_cast<int8>(std::min(uniqueSpins, static_cast<int>(m3c_node_type::k_QuadPoint))));
      if(bodyResult.invalid())
      {
        return bodyResult;
      }
    }
    int64 count = 0;
    if(effectiveCount > 0 && edgeCount > 2)
    {
      std::array<SiteIdType, 64> edgeIds{};
      int edgeIndex = 0;
      for(const auto& square : squares)
      {
        for(int index = 0; index < square.NEdge; index++)
        {
          edgeIds[edgeIndex++] = square.EdgeId[index];
        }
      }
      if(faceCenterCount == 0)
      {
        count = getNumberCase0Triangles(edgeIds, segments, edgeCount);
      }
      else if(faceCenterCount == 2)
      {
        count = getNumberCase2Triangles(edgeIds, segments, edgeCount, faceCenters, faceCenterCount);
      }
      else if(faceCenterCount > 2 && faceCenterCount <= 6)
      {
        count = getNumberCaseMTriangles(edgeIds, segments, edgeCount, faceCenters, faceCenterCount);
      }
      maximumTrianglesPerCube = std::max(maximumTrianglesPerCube, static_cast<uint64>(count));

      if((m_InputValues->BoundingBoxSkinMode == BoundingBoxSkinMode::k_BackgroundBackedWallsOnly || (m_InputValues->SharpBoundingBoxEdges && cubeTouchesBoxEdge(cube))) && count > 0)
      {
        std::vector<Triangle> countTriangles(static_cast<usize>(count));
        std::vector<SiteIdType> countCubes(static_cast<usize>(count));
        std::array<double, 3> coord1{};
        std::array<double, 3> coord2{};
        for(int component = 0; component < 3; component++)
        {
          coord1[component] = siteCoords[cube].Coord[component];
          coord2[component] = siteCoords[cube + 1 + static_cast<SiteIdType>(fileDim[0]) + numSitesPerPlane].Coord[component];
        }
        int64 generatedCount = 0;
        if(faceCenterCount == 0)
        {
          getCase0Triangles(countTriangles, countCubes, edgeIds, nodeCoords, segments, edgeCount, 0, generatedCount, coord1, coord2, cube);
        }
        else if(faceCenterCount == 2)
        {
          getCase2Triangles(countTriangles, countCubes, edgeIds, nodeCoords, segments, edgeCount, faceCenters, faceCenterCount, 0, generatedCount, coord1, coord2, cube);
        }
        else
        {
          getCaseMTriangles(countTriangles, countCubes, edgeIds, nodeCoords, segments, edgeCount, faceCenters, faceCenterCount, 0, generatedCount, (7 * (cube - 1)) + 6, coord1, coord2, cube);
        }
        if(generatedCount != count)
        {
          return MakeErrorResult(-90559, "M3C out-of-core pruning generation disagrees with its counted triangle range.");
        }

        int64 survivingCount = 0;
        for(const Triangle& triangle : countTriangles)
        {
          const bool dropTriangle = m_InputValues->BoundingBoxSkinMode == BoundingBoxSkinMode::k_BackgroundBackedWallsOnly && skipBackgroundSkinFace(triangle);
          const uint8 referenceFlag = dropTriangle ? uint8{1} : uint8{2};
          for(const SiteIdType nodeId : triangle.NodeId)
          {
            auto nodeResult = candidateNodes.cache().read(static_cast<uint64>(nodeId), m_ShouldCancel);
            if(nodeResult.invalid())
            {
              return ConvertResult(std::move(nodeResult));
            }
            auto node = nodeResult.value();
            node.PruneReferences = static_cast<uint8>(node.PruneReferences | referenceFlag);
            auto writeResult = candidateNodes.cache().write(static_cast<uint64>(nodeId), node, m_ShouldCancel);
            if(writeResult.invalid())
            {
              return writeResult;
            }
          }
          if(dropTriangle)
          {
            numFacesPruned++;
          }
          else
          {
            survivingCount++;
            if(m_InputValues->SharpBoundingBoxEdges && touchesBoxEdge(triangle))
            {
              edgeTriangles.push_back(triangle);
              edgeCubes.push_back(cube);
            }
          }
        }
        count = survivingCount;
      }
    }
    auto writeResult = triangleCounts.cache().write(static_cast<uint64>(cube), count, m_ShouldCancel);
    if(writeResult.invalid())
    {
      return writeResult;
    }
  }
  SharpEdgeResult sharpEdges;
  if(m_InputValues->SharpBoundingBoxEdges && !edgeTriangles.empty())
  {
    // Only edge-adjacent triangles remain resident. Their count grows with the sum of the three dimensions, not the volume.
    std::unordered_map<SiteIdType, int8> edgeNodeTypes;
    std::unordered_map<SiteIdType, int64> cubeCountChanges;
    for(usize face = 0; face < edgeTriangles.size(); face++)
    {
      const auto& triangle = edgeTriangles[face];
      cubeCountChanges[edgeCubes[face]]--;
      for(const SiteIdType nodeId : triangle.NodeId)
      {
        auto [entry, inserted] = edgeNodeTypes.try_emplace(nodeId, 0);
        if(inserted)
        {
          auto node = candidateNodes.cache().read(static_cast<uint64>(nodeId), m_ShouldCancel);
          if(node.invalid())
          {
            return ConvertResult(std::move(node));
          }
          entry->second = node.value().Type;
        }
        if((triangle.NSpin[0] < 0) != (triangle.NSpin[1] < 0) && entry->second < 10)
        {
          entry->second = static_cast<int8>(entry->second + 10);
        }
      }
    }
    std::vector<SiteIdType> edgeNodeIds;
    edgeNodeIds.reserve(edgeNodeTypes.size());
    for(const auto& [nodeId, type] : edgeNodeTypes)
    {
      edgeNodeIds.push_back(nodeId);
    }
    // Ascending candidate order preserves the in-core representative and vertex ordering.
    std::sort(edgeNodeIds.begin(), edgeNodeIds.end());
    sharpEdges = sharpenBoundingBoxEdges(edgeTriangles, edgeCubes, edgeNodeTypes, 0, nodeCoords, dims, nonstd::span<const SiteIdType>(edgeNodeIds));
    for(const SiteIdType cube : edgeCubes)
    {
      cubeCountChanges[cube]++;
    }
    for(const auto& [cube, change] : cubeCountChanges)
    {
      auto count = triangleCounts.cache().read(static_cast<uint64>(cube), m_ShouldCancel);
      if(count.invalid())
      {
        return ConvertResult(std::move(count));
      }
      auto write = triangleCounts.cache().write(static_cast<uint64>(cube), count.value() + change, m_ShouldCancel);
      if(write.invalid())
      {
        return write;
      }
    }
    for(const auto& [nodeId, type] : edgeNodeTypes)
    {
      auto node = candidateNodes.cache().read(static_cast<uint64>(nodeId), m_ShouldCancel);
      if(node.invalid())
      {
        return ConvertResult(std::move(node));
      }
      auto record = node.value();
      record.Type = type;
      auto write = candidateNodes.cache().write(static_cast<uint64>(nodeId), record, m_ShouldCancel);
      if(write.invalid())
      {
        return write;
      }
    }
  }
  std::vector<Triangle>().swap(edgeTriangles);
  std::vector<SiteIdType>().swap(edgeCubes);
  auto flushNodesResult = candidateNodes.flush(m_ShouldCancel);
  if(flushNodesResult.invalid())
  {
    return flushNodesResult;
  }
  auto flushCountsResult = triangleCounts.flush(m_ShouldCancel);
  if(flushCountsResult.invalid())
  {
    return flushCountsResult;
  }

  // Convert counts in place to the 1-based deterministic cube offsets used by
  // the parallel path. The values remain in external storage for pass 2.
  uint64 triangleTotal = 0;
  for(SiteIdType cube = 1; cube <= lastCube; cube++)
  {
    if(m_ShouldCancel)
    {
      return {};
    }
    auto countResult = triangleCounts.cache().read(static_cast<uint64>(cube), m_ShouldCancel);
    if(countResult.invalid())
    {
      return ConvertResult(std::move(countResult));
    }
    const int64 count = countResult.value();
    if(count < 0 || static_cast<uint64>(count) > static_cast<uint64>(std::numeric_limits<int64>::max()) - triangleTotal ||
       triangleTotal > static_cast<uint64>(std::numeric_limits<usize>::max()) - static_cast<uint64>(count))
    {
      return MakeErrorResult(-90553, "M3C out-of-core triangle count or offset overflows its output range.");
    }
    auto offsetResult = triangleCounts.cache().write(static_cast<uint64>(cube), static_cast<int64>(triangleTotal), m_ShouldCancel);
    if(offsetResult.invalid())
    {
      return offsetResult;
    }
    triangleTotal += static_cast<uint64>(count);
  }
  auto offsetFlushResult = triangleCounts.flush(m_ShouldCancel);
  if(offsetFlushResult.invalid())
  {
    return offsetFlushResult;
  }

  // Candidate IDs compact in ascending candidate order, exactly matching the
  // original assign_new_nodeID traversal. Vertex/NodeTypes output waits until
  // generation has marked exterior nodes.
  uint64 nodeTotal = 0;
  for(uint64 candidate = 0; candidate < candidateCount; candidate++)
  {
    if(m_ShouldCancel)
    {
      return {};
    }
    auto nodeResult = candidateNodes.cache().read(candidate, m_ShouldCancel);
    if(nodeResult.invalid())
    {
      return ConvertResult(std::move(nodeResult));
    }
    auto node = nodeResult.value();
    bool nodeChanged = false;
    if(m_InputValues->BoundingBoxSkinMode == BoundingBoxSkinMode::k_BackgroundBackedWallsOnly && node.PruneReferences == uint8{1} && node.Type > 0)
    {
      node.Type = m3c_node_type::k_Unused;
      nodeChanged = true;
    }
    if(node.Type > 0)
    {
      if(nodeTotal >= std::numeric_limits<usize>::max())
      {
        return MakeErrorResult(-90554, "M3C out-of-core compacted node count overflows its output range.");
      }
      node.CompactId = nodeTotal++;
      nodeChanged = true;
    }
    if(nodeChanged)
    {
      auto writeResult = candidateNodes.cache().write(candidate, node, m_ShouldCancel);
      if(writeResult.invalid())
      {
        return writeResult;
      }
    }
  }
  auto compactFlushResult = candidateNodes.flush(m_ShouldCancel);
  if(compactFlushResult.invalid())
  {
    return compactFlushResult;
  }

  auto& triangleGeom = m_DataStructure.getDataRefAs<TriangleGeom>(m_InputValues->TriangleGeometryPath);
  Result<> resizeResult = triangleGeom.resizeVertexList(nodeTotal);
  if(resizeResult.invalid())
  {
    return resizeResult;
  }
  resizeResult = triangleGeom.resizeFaceList(triangleTotal);
  if(resizeResult.invalid())
  {
    return resizeResult;
  }
  resizeResult = triangleGeom.getVertexAttributeMatrix()->resizeTuples({static_cast<usize>(nodeTotal)});
  if(resizeResult.invalid())
  {
    return resizeResult;
  }
  resizeResult = triangleGeom.getFaceAttributeMatrix()->resizeTuples({static_cast<usize>(triangleTotal)});
  if(resizeResult.invalid())
  {
    return resizeResult;
  }
  auto& faceStore = triangleGeom.getFaces()->getDataStoreRef();
  auto& faceLabelsStore = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->FaceLabelsDataPath).getDataStoreRef();
  auto& nodeTypesStore = m_DataStructure.getDataRefAs<Int8Array>(m_InputValues->NodeTypesDataPath).getDataStoreRef();
  resizeResult = faceLabelsStore.resizeTuples({static_cast<usize>(triangleTotal)});
  if(resizeResult.invalid())
  {
    return resizeResult;
  }
  resizeResult = nodeTypesStore.resizeTuples({static_cast<usize>(nodeTotal)});
  if(resizeResult.invalid())
  {
    return resizeResult;
  }

  std::vector<std::shared_ptr<AbstractTupleTransfer>> transfers;
  for(usize index = 0; index < m_InputValues->SelectedCellDataArrayPaths.size(); index++)
  {
    AddTupleTransferInstance(m_DataStructure, m_InputValues->SelectedCellDataArrayPaths[index], m_InputValues->CreatedDataArrayPaths[index], transfers);
  }
  const usize cellArrayCount = m_InputValues->SelectedCellDataArrayPaths.size();
  for(usize index = 0; index < m_InputValues->SelectedFeatureDataArrayPaths.size(); index++)
  {
    AddFeatureTupleTransferInstance(m_DataStructure, m_InputValues->SelectedFeatureDataArrayPaths[index], m_InputValues->CreatedDataArrayPaths[cellArrayCount + index],
                                    m_InputValues->FeatureIdsArrayPath, transfers);
  }
  const auto outputLabel = [maxGrainId](const int spin) { return spin < 0 ? int32{-1} : (spin == maxGrainId ? int32{0} : spin); };
  const auto sourceCell = [&](const int label, const SiteIdType cube) -> Result<usize> {
    const Neighbor siteNeighbors = neighbors[cube];
    const std::array<SiteIdType, 8> corners = {
        cube, siteNeighbors.NeighId[1], siteNeighbors.NeighId[7], siteNeighbors.NeighId[8], siteNeighbors.NeighId[18], siteNeighbors.NeighId[19], siteNeighbors.NeighId[25], siteNeighbors.NeighId[26]};
    for(const SiteIdType site : corners)
    {
      auto value = sourceValue(site);
      if(value.invalid())
      {
        return ConvertInvalidResult<usize>(std::move(value));
      }
      if(value.value() == label)
      {
        const usize linear = static_cast<usize>(site - 1);
        const usize xIndex = linear % fileDim[0];
        const usize yIndex = (linear / fileDim[0]) % fileDim[1];
        const usize zIndex = linear / (fileDim[0] * fileDim[1]);
        if(xIndex >= 1 && xIndex <= dims[0] && yIndex >= 1 && yIndex <= dims[1] && zIndex >= 1 && zIndex <= dims[2])
        {
          return {((zIndex - 1) * sourceSliceSize) + ((yIndex - 1) * dims[0]) + (xIndex - 1)};
        }
      }
    }
    return {std::numeric_limits<usize>::max()};
  };

  constexpr usize kFaceBatch = 16384;
  std::vector<IGeometry::MeshIndexType> faceValues(kFaceBatch * 3);
  std::vector<int32> labelValues(kFaceBatch * 2);
  std::vector<QuickSurfaceTransferData> transferValues(kFaceBatch);
  const usize localCapacity = std::max<usize>(1, static_cast<usize>(maximumTrianglesPerCube));
  std::vector<Triangle> localTriangles(localCapacity);
  std::vector<SiteIdType> localCubes(localCapacity);
  for(SiteIdType cube = 1; cube <= lastCube; cube++)
  {
    if(m_ShouldCancel)
    {
      return {};
    }
    auto offset = triangleCounts.cache().read(static_cast<uint64>(cube), m_ShouldCancel);
    if(offset.invalid())
    {
      return ConvertResult(std::move(offset));
    }
    const usize destination = static_cast<usize>(offset.value());
    const std::array<SiteIdType, 6> squareIds = {
        3 * (cube - 1), (3 * (cube - 1)) + 1, (3 * (cube - 1)) + 2, (3 * cube) + 2, (3 * (cube + static_cast<SiteIdType>(fileDim[0]) - 1)) + 1, 3 * (cube + numSitesPerPlane - 1)};
    std::array<Face, 6> squares{};
    std::array<Segment, 64> segments{};
    std::array<SiteIdType, 64> edgeIds{};
    std::array<SiteIdType, 6> centers = {-1, -1, -1, -1, -1, -1};
    int segmentCount = 0;
    int edgeCount = 0;
    int centerCount = 0;
    int effectiveCount = 0;
    for(int square = 0; square < 6; square++)
    {
      auto result = buildSquare(squareIds[square], squares[square], segments, segmentCount, false);
      if(result.invalid())
      {
        return result;
      }
      if(squares[square].FaceCenterNode != -1)
      {
        centers[centerCount++] = squares[square].FaceCenterNode;
      }
      effectiveCount += squares[square].Effect;
      for(int edge = 0; edge < squares[square].NEdge; edge++)
      {
        edgeIds[edgeCount++] = squares[square].EdgeId[edge];
      }
    }
    usize generated = 0;
    if(effectiveCount > 0 && edgeCount > 2)
    {
      std::array<double, 3> coord1{};
      std::array<double, 3> coord2{};
      for(int component = 0; component < 3; component++)
      {
        coord1[component] = siteCoords[cube].Coord[component];
        coord2[component] = siteCoords[cube + 1 + static_cast<SiteIdType>(fileDim[0]) + numSitesPerPlane].Coord[component];
      }
      int64 end = 0;
      if(centerCount == 0)
      {
        getCase0Triangles(localTriangles, localCubes, edgeIds, nodeCoords, segments, edgeCount, 0, end, coord1, coord2, cube);
      }
      else if(centerCount == 2)
      {
        getCase2Triangles(localTriangles, localCubes, edgeIds, nodeCoords, segments, edgeCount, centers, centerCount, 0, end, coord1, coord2, cube);
      }
      else if(centerCount > 2)
      {
        getCaseMTriangles(localTriangles, localCubes, edgeIds, nodeCoords, segments, edgeCount, centers, centerCount, 0, end, (7 * (cube - 1)) + 6, coord1, coord2, cube);
      }
      generated = static_cast<usize>(end);
    }
    if(m_InputValues->BoundingBoxSkinMode == BoundingBoxSkinMode::k_BackgroundBackedWallsOnly)
    {
      usize survivingCount = 0;
      for(usize index = 0; index < generated; index++)
      {
        if(!skipBackgroundSkinFace(localTriangles[index]))
        {
          localTriangles[survivingCount] = localTriangles[index];
          localCubes[survivingCount] = localCubes[index];
          survivingCount++;
        }
      }
      generated = survivingCount;
    }
    if(m_InputValues->SharpBoundingBoxEdges)
    {
      usize survivingCount = 0;
      for(usize index = 0; index < generated; index++)
      {
        Triangle triangle = localTriangles[index];
        if(remapSharpEdgeTriangle(triangle, sharpEdges, nodeCoords))
        {
          localTriangles[survivingCount] = triangle;
          localCubes[survivingCount] = localCubes[index];
          survivingCount++;
        }
      }
      generated = survivingCount;
    }
    uint64 expectedEnd = triangleTotal;
    if(cube != lastCube)
    {
      auto nextOffset = triangleCounts.cache().read(static_cast<uint64>(cube + 1), m_ShouldCancel);
      if(nextOffset.invalid() || nextOffset.value() < 0)
      {
        return nextOffset.invalid() ? ConvertResult(std::move(nextOffset)) : MakeErrorResult(-90558, "M3C out-of-core triangle offset is negative.");
      }
      expectedEnd = static_cast<uint64>(nextOffset.value());
    }
    if(generated > localCapacity || expectedEnd < static_cast<uint64>(destination) || generated != expectedEnd - static_cast<uint64>(destination))
    {
      return MakeErrorResult(-90555, "M3C out-of-core generation disagrees with its counted triangle range.");
    }
    for(usize start = 0; start < generated; start += kFaceBatch)
    {
      const usize count = std::min(kFaceBatch, generated - start);
      for(usize local = 0; local < count; local++)
      {
        const Triangle& triangle = localTriangles[start + local];
        for(int vertex = 0; vertex < 3; vertex++)
        {
          auto nodeResult = candidateNodes.cache().read(static_cast<uint64>(triangle.NodeId[vertex]), m_ShouldCancel);
          if(nodeResult.invalid() || nodeResult.value().Type <= 0)
          {
            return nodeResult.invalid() ? ConvertResult(std::move(nodeResult)) : MakeErrorResult(-90556, "M3C out-of-core triangle references an unused candidate node.");
          }
          auto node = nodeResult.value();
          if((triangle.NSpin[0] < 0) != (triangle.NSpin[1] < 0) && node.Type < 10)
          {
            node.Type = static_cast<int8>(node.Type + 10);
            auto write = candidateNodes.cache().write(static_cast<uint64>(triangle.NodeId[vertex]), node, m_ShouldCancel);
            if(write.invalid())
            {
              return write;
            }
          }
          faceValues[(local * 3) + vertex] = node.CompactId;
        }
        const int32 labelA = outputLabel(triangle.NSpin[0]);
        const int32 labelB = outputLabel(triangle.NSpin[1]);
        const bool aFirst = labelA <= labelB;
        labelValues[local * 2] = aFirst ? labelA : labelB;
        labelValues[(local * 2) + 1] = aFirst ? labelB : labelA;
        auto first = sourceCell(aFirst ? triangle.NSpin[0] : triangle.NSpin[1], cube);
        auto second = sourceCell(aFirst ? triangle.NSpin[1] : triangle.NSpin[0], cube);
        if(first.invalid() || second.invalid())
        {
          return first.invalid() ? ConvertResult(std::move(first)) : ConvertResult(std::move(second));
        }
        transferValues[local] = {destination + start + local, first.value(), second.value(), labelValues[local * 2], labelValues[(local * 2) + 1]};
      }
      auto faceWrite = faceStore.copyFromBuffer((destination + start) * 3, nonstd::span<const IGeometry::MeshIndexType>(faceValues.data(), count * 3));
      auto labelWrite = faceLabelsStore.copyFromBuffer((destination + start) * 2, nonstd::span<const int32>(labelValues.data(), count * 2));
      if(faceWrite.invalid() || labelWrite.invalid())
      {
        return faceWrite.invalid() ? faceWrite : labelWrite;
      }
      for(const auto& transfer : transfers)
      {
        auto result = transfer->quickSurfaceTransferBatch(nonstd::span<const QuickSurfaceTransferData>(transferValues.data(), count));
        if(result.invalid())
        {
          return result;
        }
      }
    }
  }
  auto promotionFlush = candidateNodes.flush(m_ShouldCancel);
  if(promotionFlush.invalid())
  {
    return promotionFlush;
  }

  auto& vertexStore = triangleGeom.getVertices()->getDataStoreRef();
  constexpr usize kVertexBatch = 16384;
  std::vector<float32> vertexValues(kVertexBatch * 3);
  std::vector<int8> typeValues(kVertexBatch);
  for(uint64 candidate = 0; candidate < candidateCount;)
  {
    if(m_ShouldCancel)
    {
      return {};
    }
    const uint64 end = std::min(candidateCount, candidate + static_cast<uint64>(kVertexBatch));
    usize count = 0;
    uint64 firstCompact = 0;
    for(uint64 current = candidate; current < end; current++)
    {
      auto nodeResult = candidateNodes.cache().read(current, m_ShouldCancel);
      if(nodeResult.invalid())
      {
        return ConvertResult(std::move(nodeResult));
      }
      const auto& node = nodeResult.value();
      if(node.Type > 0)
      {
        if(count == 0)
        {
          firstCompact = node.CompactId;
        }
        const auto snapped = sharpEdges.SnappedCoords.find(static_cast<SiteIdType>(current));
        const Node coordinate = snapped == sharpEdges.SnappedCoords.end() ? nodeCoords[static_cast<SiteIdType>(current)] : snapped->second;
        vertexValues[count * 3] = coordinate.Coord[0];
        vertexValues[(count * 3) + 1] = coordinate.Coord[1];
        vertexValues[(count * 3) + 2] = coordinate.Coord[2];
        typeValues[count++] = node.Type;
      }
    }
    if(count > 0)
    {
      auto vertexWrite = vertexStore.copyFromBuffer(firstCompact * 3, nonstd::span<const float32>(vertexValues.data(), count * 3));
      auto typeWrite = nodeTypesStore.copyFromBuffer(firstCompact, nonstd::span<const int8>(typeValues.data(), count));
      if(vertexWrite.invalid() || typeWrite.invalid())
      {
        return vertexWrite.invalid() ? vertexWrite : typeWrite;
      }
    }
    candidate = end;
  }
  if(m_InputValues->RepairTriangleWinding)
  {
    auto& ioCollection = DataStoreUtilities::GetIOCollection();
    if(ioCollection.hasExternalSortCapability() && ioCollection.hasTemporaryRecordStoreCapability())
    {
      auto result = MeshingUtilities::RepairTriangleWindingExternal(faceStore, faceLabelsStore, m_ShouldCancel, m_MessageHandler);
      if(result.invalid())
      {
        return result;
      }
    }
    else if(usesOutOfCoreStore)
    {
      return MakeErrorResult(-90557, "M3C out-of-core winding repair requires external-sort and temporary-record-store providers.");
    }
    else
    {
      triangleGeom.findElementNeighbors(true);
      const auto optionalId = triangleGeom.getElementNeighborsId();
      if(optionalId.has_value())
      {
        const auto& connectivity = m_DataStructure.getDataRefAs<IGeometry::ElementDynamicList>(optionalId.value());
        auto result = MeshingUtilities::RepairTriangleWinding(faceStore, connectivity, faceLabelsStore, m_ShouldCancel, m_MessageHandler);
        m_DataStructure.removeData(triangleGeom.getElementContainingVertId().value());
        m_DataStructure.removeData(triangleGeom.getElementNeighborsId().value());
        if(result.invalid())
        {
          return result;
        }
      }
    }
  }
  if(m_InputValues->BoundingBoxSkinMode == BoundingBoxSkinMode::k_BackgroundBackedWallsOnly)
  {
    if(triangleTotal == 0)
    {
      return MeshingUtilities::MakeEmptyMeshWarning(m_InputValues->TriangleGeometryPath, featureIds.getNumberOfTuples(), static_cast<usize>(nodeTotal));
    }
    if(numFacesPruned == 0)
    {
      return MeshingUtilities::MakeNoFacesPrunedWarning(m_InputValues->TriangleGeometryPath);
    }
  }
  return {};
}

Result<> M3CSurfaceMeshing::runEntireVolume()
{
  // M3C node coordinates require the uniform spacing of ImageGeom.
  const auto& imageGeom = m_DataStructure.getDataRefAs<ImageGeom>(m_InputValues->GridGeomDataPath);
  const auto& featureIds = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->FeatureIdsArrayPath);
  const auto& featureIdsStore = featureIds.getDataStoreRef();

  SizeVec3 gridDims = imageGeom.getDimensions();
  std::array<usize, 3> dims = {gridDims[0], gridDims[1], gridDims[2]};
  const FloatVec3 spacing = imageGeom.getSpacing();
  const FloatVec3 imgOrigin = imageGeom.getOrigin();
  const std::array<float, 3> res = {spacing[0], spacing[1], spacing[2]};
  const std::array<float, 3> origin = {imgOrigin[0], imgOrigin[1], imgOrigin[2]};

  // NX inputs do not include an exterior layer. Add one for surface closure.
  constexpr bool addSurfaceLayer = true;
  std::array<usize, 3> fileDim = {dims[0] + 2, dims[1] + 2, dims[2] + 2};
  const usize totalPoints = fileDim[0] * fileDim[1] * fileDim[2];
  // SiteIdType supports grids with more than 2^31 voxels. The 32-bit edge and node
  // storage still limits mesh size to approximately 2^32 elements.
  const SiteIdType numSites = static_cast<SiteIdType>(totalPoints);
  const usize paddedSitesPerPlane = fileDim[0] * fileDim[1];
  const SiteIdType numSitesPerPlane = static_cast<SiteIdType>(paddedSitesPerPlane);

  // Read Feature Ids directly into the padded grid. Zero-label renumbering changes
  // only this working copy.
  m_MessageHandler.sendInfoMessage("Initializing working grid and ghost layer...");
  std::vector<int32> point(totalPoints + 1, 0);
  const int maxGrainId = initializeMicro(addSurfaceLayer, dims, fileDim, featureIdsStore, point);

  // Calculate coordinates and neighbors on demand to avoid three full-volume arrays.
  const SiteCoords siteCoords{fileDim[0], fileDim[1], fileDim[0] * fileDim[1], {res[0], res[1], res[2]}, {origin[0], origin[1], origin[2]}};
  const NodeCoords nodeCoords{siteCoords};
  const NeighborAccessor neighbors{numSites, numSitesPerPlane, static_cast<int>(fileDim[0])};

  // Each site owns top, back, and left squares. Seven candidate node types start unused.
  m_MessageHandler.sendInfoMessage("Initializing candidate nodes and squares...");
  std::vector<Face> squares(static_cast<usize>(3) * numSites);
  std::vector<int8> nodeType(static_cast<usize>(7) * numSites, 0);
  initializeSquares(squares, numSites);

  if(m_ShouldCancel)
  {
    return {};
  }

  // Count face edges before their exact allocation.
  m_MessageHandler.sendInfoMessage("Counting face edges...");
  const int64 nFEdge = getNumberFEdges(squares, point, neighbors, numSites, m_ShouldCancel);

  m_MessageHandler.sendInfoMessage("Finding nodes and edges on each square...");
  std::vector<Segment> fedges(static_cast<usize>(nFEdge < 0 ? 0 : nFEdge));
  getNodesFEdges(squares, point, neighbors, nodeType, fedges, numSites, numSitesPerPlane, static_cast<int>(fileDim[0]), m_ShouldCancel);

  if(m_ShouldCancel)
  {
    return {};
  }

  // Count triangles before their exact allocation.
  m_MessageHandler.sendInfoMessage("Counting triangles...");
  const int64 nTriangle = getNumberTriangles(point, squares, neighbors, nodeType, fedges, numSites, numSitesPerPlane, static_cast<int>(fileDim[0]), m_ShouldCancel);

  m_MessageHandler.sendInfoMessage("Generating triangles...");
  std::vector<Triangle> triangles(static_cast<usize>(nTriangle < 0 ? 0 : nTriangle));
  std::vector<SiteIdType> mCubeID(static_cast<usize>(nTriangle < 0 ? 0 : nTriangle), 0);
  getTriangles(siteCoords, triangles, mCubeID, squares, nodeCoords, fedges, numSites, numSitesPerPlane, static_cast<int>(fileDim[0]), m_ShouldCancel);

  if(m_ShouldCancel)
  {
    return {};
  }

  return finalizeMesh(m_DataStructure, m_InputValues, m_MessageHandler, m_ShouldCancel, triangles, mCubeID, fedges, nodeType, point, nodeCoords, neighbors, numSites, fileDim, dims, maxGrainId);
}

Result<> M3CSurfaceMeshing::runWindowed(const bool parallel) const
{
  // The sliding window keeps marching-square scratch to two Z slices. Serial
  // execution matches runEntireVolume. Parallel execution can change triangulation.
  const auto& imageGeom = m_DataStructure.getDataRefAs<ImageGeom>(m_InputValues->GridGeomDataPath);
  const auto& featureIds = m_DataStructure.getDataRefAs<Int32Array>(m_InputValues->FeatureIdsArrayPath);
  const auto& featureIdsStore = featureIds.getDataStoreRef();

  SizeVec3 gridDims = imageGeom.getDimensions();
  std::array<usize, 3> dims = {gridDims[0], gridDims[1], gridDims[2]};
  const FloatVec3 spacing = imageGeom.getSpacing();
  const FloatVec3 imgOrigin = imageGeom.getOrigin();
  const std::array<float, 3> res = {spacing[0], spacing[1], spacing[2]};
  const std::array<float, 3> origin = {imgOrigin[0], imgOrigin[1], imgOrigin[2]};

  constexpr bool addSurfaceLayer = true;
  std::array<usize, 3> fileDim = {dims[0] + 2, dims[1] + 2, dims[2] + 2};
  const usize totalPoints = fileDim[0] * fileDim[1] * fileDim[2];
  const SiteIdType numSites = static_cast<SiteIdType>(totalPoints);
  const usize paddedSitesPerPlane = fileDim[0] * fileDim[1];
  const SiteIdType numSitesPerPlane = static_cast<SiteIdType>(paddedSitesPerPlane);
  const int xDim = static_cast<int>(fileDim[0]);

  m_MessageHandler.sendInfoMessage("Initializing working grid and ghost layer...");
  std::vector<int32> point(totalPoints + 1, 0);
  const int maxGrainId = initializeMicro(addSurfaceLayer, dims, fileDim, featureIdsStore, point);

  const SiteCoords siteCoords{fileDim[0], fileDim[1], fileDim[0] * fileDim[1], {res[0], res[1], res[2]}, {origin[0], origin[1], origin[2]}};
  const NodeCoords nodeCoords{siteCoords};
  const NeighborAccessor neighbors{numSites, numSitesPerPlane, xDim};

  // Mesh-scale vectors grow with output size. Only per-site squares are windowed.
  // Face edges reserve once from an early rate estimate. Triangles and cube IDs
  // resize once after the first pass counts their exact size.
  std::vector<int8> nodeType(static_cast<usize>(7) * numSites, 0);
  std::vector<Segment> fedges;
  std::vector<Triangle> triangles;
  std::vector<SiteIdType> mCubeID;

  // The square window holds two Z slices. It slides one slice as the cube sweep
  // advances and keeps absolute square IDs mapped to local slots.
  const SiteIdType sliceSquares = 3 * numSitesPerPlane; // squares per z-slice
  std::vector<Face> window(static_cast<usize>(2) * sliceSquares);
  SiteIdType winBaseSite = 1;
  auto winIndex = [&winBaseSite](const SiteIdType squareId) -> usize { return static_cast<usize>(squareId - (3 * (winBaseSite - 1))); };

  // Build square edges and node types in ascending square order. This preserves
  // global face-edge IDs and the legacy effect flag.
  auto computeSquares = [&](const SiteIdType kLo, const SiteIdType kHi, const bool appendEdges, int64& eid) {
    for(SiteIdType k = kLo; k < kHi; k++)
    {
      Face& sqk = window[winIndex(k)];
      for(int j = 0; j < 4; j++)
      {
        sqk.EdgeId[j] = k_UnusedNodeId;
      }
      sqk.NEdge = 0;
      sqk.FaceCenterNode = -1;
      sqk.Effect = 0;

      const SiteIdType cubeOrigin = (k / 3) + 1;
      const int sqOrder = static_cast<int>(k % 3);
      const std::array<SiteIdType, 4> tnsite = squareCorners(k, neighbors);
      std::array<int, 4> tnspin{};
      int numGhostCorners = 0;
      for(int cornerIdx = 0; cornerIdx < 4; cornerIdx++)
      {
        tnspin[cornerIdx] = point[tnsite[cornerIdx]];
        if(tnspin[cornerIdx] < 0)
        {
          numGhostCorners++;
        }
      }
      if(numGhostCorners != 4)
      {
        sqk.Effect = 1;
      }

      int edgeCount = 0;
      if(numGhostCorners != 4)
      {
        int sqIndex = getSquareIndex(tnspin);
        if(sqIndex == 15)
        {
          sqIndex = sqIndex + treatAnomaly(tnsite, point, neighbors, k);
        }
        if(sqIndex != 0)
        {
          for(int j = 0; j < 8; j = j + 2)
          {
            if(k_EdgeTable2d[sqIndex][j] != -1)
            {
              std::array<int, 2> nodeIndex = {k_EdgeTable2d[sqIndex][j], k_EdgeTable2d[sqIndex][j + 1]};
              const std::array<int, 2> pixIndex = {k_NsTable2d[sqIndex][j], k_NsTable2d[sqIndex][j + 1]};
              std::array<SiteIdType, 2> nodeID{};
              std::array<int, 2> pixSpin{};
              getNodes(cubeOrigin, sqOrder, nodeIndex, nodeID, numSitesPerPlane, xDim);
              getSpins(point, cubeOrigin, sqOrder, pixIndex, pixSpin, numSitesPerPlane, xDim);

              if(pixSpin[0] > 0 || pixSpin[1] > 0)
              {
                Segment seg{};
                seg.NodeId[0] = nodeID[0];
                seg.NodeId[1] = nodeID[1];
                seg.NSpin[0] = pixSpin[0];
                seg.NSpin[1] = pixSpin[1];
                sqk.EdgeId[edgeCount] = static_cast<uint32>(eid);
                if(appendEdges)
                {
                  fedges.push_back(seg);
                }
                eid++;
                edgeCount++;
              }
              else
              {
                nodeType[nodeID[0]] = m3c_node_type::k_Unused;
                nodeType[nodeID[1]] = m3c_node_type::k_Unused;
              }

              for(int ii = 0; ii < 2; ii++)
              {
                if(nodeIndex[ii] == 4)
                {
                  if(sqIndex == 7 || sqIndex == 11 || sqIndex == 13 || sqIndex == 14)
                  {
                    const SiteIdType tnode = nodeID[ii];
                    sqk.FaceCenterNode = tnode;
                    nodeType[tnode] = m3c_node_type::k_TriplePoint;
                  }
                  else if(sqIndex == 19)
                  {
                    const SiteIdType tnode = nodeID[ii];
                    sqk.FaceCenterNode = tnode;
                    nodeType[tnode] = m3c_node_type::k_QuadPoint;
                  }
                }
                else
                {
                  // Interior edge endpoints must remain real nodes after compaction.
                  const SiteIdType tnode = nodeID[ii];
                  nodeType[tnode] = m3c_node_type::k_Default;
                }
              }
            }
          }
        }
      }
      sqk.NEdge = static_cast<int8>(edgeCount);
    }
  };

  // Case functions flip shared face edges while tracing loops. Count every cube
  // before generating any triangle to preserve the serial ordering.
  int64 nTriangle = 0;

  const int64 totalSlices = (numSitesPerPlane > 0) ? (numSites / numSitesPerPlane) : 1;
  const int64 progressStep = std::max<int64>(1, totalSlices / 20); // ~20 progress updates per sweep

  // Estimate face-edge capacity after early slices. One reservation avoids repeated
  // multi-gigabyte reallocations and transient memory spikes.
  bool fedgesReserved = false;
  const int64 reserveAfterSlices = std::max<int64>(4, totalSlices / 16);
  auto maybeReserveFedges = [&](const int64 eid) {
    const int64 sliceIdx = winBaseSite / numSitesPerPlane;
    if(fedgesReserved || sliceIdx < reserveAfterSlices || sliceIdx >= totalSlices)
    {
      return;
    }
    fedgesReserved = true;
    // The window has computed two slices beyond the slide counter.
    const int64 slicesComputed = sliceIdx + 2;
    const auto projected = static_cast<usize>(static_cast<double>(eid) / static_cast<double>(slicesComputed) * static_cast<double>(totalSlices) * 1.05);
    if(projected > fedges.capacity())
    {
      fedges.reserve(projected);
    }
  };

  auto sweep = [&](const bool appendEdges, const bool generate) {
    winBaseSite = 1;
    int64 eid = 0;
    int64 tidRun = 0;
    computeSquares(0, std::min<SiteIdType>(2 * sliceSquares, 3 * numSites), appendEdges, eid);

    for(SiteIdType i = 1; i <= (numSites - numSitesPerPlane); i++)
    {
      if(m_ShouldCancel)
      {
        return;
      }

      // Slide the window to cover the current and next site planes.
      while(i >= winBaseSite + numSitesPerPlane)
      {
        std::memmove(window.data(), window.data() + sliceSquares, static_cast<usize>(sliceSquares) * sizeof(Face));
        winBaseSite += numSitesPerPlane;
        const int64 sliceIdx = winBaseSite / numSitesPerPlane;
        if(sliceIdx % progressStep == 0)
        {
          m_MessageHandler.sendInfoMessage(fmt::format("Sweeping z-slices ({}): slice {} / {}", generate ? "pass 2, generating triangles" : "pass 1, counting", sliceIdx, totalSlices));
        }
        const SiteIdType newLoSquare = 3 * (winBaseSite + numSitesPerPlane - 1);
        const SiteIdType newHiSquare = std::min<SiteIdType>(3 * (winBaseSite + (2 * numSitesPerPlane) - 1), 3 * numSites);
        if(newLoSquare < newHiSquare)
        {
          computeSquares(newLoSquare, newHiSquare, appendEdges, eid);
        }
        if(appendEdges)
        {
          maybeReserveFedges(eid);
        }
      }

      std::array<SiteIdType, 6> sqID{};
      sqID[0] = 3 * (i - 1);
      sqID[1] = (3 * (i - 1)) + 1;
      sqID[2] = (3 * (i - 1)) + 2;
      sqID[3] = (3 * i) + 2;
      sqID[4] = (3 * (i + xDim - 1)) + 1;
      sqID[5] = 3 * (i + numSitesPerPlane - 1);

      std::array<SiteIdType, 6> arrayFC{};
      for(int ii = 0; ii < 6; ii++)
      {
        arrayFC[ii] = -1;
      }
      int fcid = 0;
      int nFE = 0;
      int eff = 0;
      for(int ii = 0; ii < 6; ii++)
      {
        const Face& sqf = window[winIndex(sqID[ii])];
        if(sqf.FaceCenterNode != -1)
        {
          arrayFC[fcid] = sqf.FaceCenterNode;
          fcid++;
        }
        nFE = nFE + sqf.NEdge;
        eff = eff + sqf.Effect;
      }
      const int nFC = fcid;
      const int cubeFlag = (eff > 0) ? 1 : 0;
      const SiteIdType bodyCenterNode = (7 * (i - 1)) + 6;

      // The count pass assigns the body-center type when three or more face
      // centers meet. This timing matches the whole-volume path.
      if(!generate && nFC >= 3)
      {
        const std::array<SiteIdType, 4> corners1 = squareCorners(sqID[0], neighbors);
        const std::array<SiteIdType, 4> corners2 = squareCorners(sqID[5], neighbors);
        std::array<int, 8> arraySpin{};
        for(int j = 0; j < 4; j++)
        {
          arraySpin[j] = point[corners1[j]];
          arraySpin[j + 4] = point[corners2[j]];
        }
        int nds = 0;
        for(int k = 0; k < 8; k++)
        {
          const int cspin = arraySpin[k];
          if(cspin != -1)
          {
            nds++;
            arraySpin[k] = -1;
            for(int kk = 0; kk < 8; kk++)
            {
              if(cspin == arraySpin[kk])
              {
                arraySpin[kk] = -1;
              }
            }
          }
        }
        // NodeType uses k_QuadPoint for four or more labels.
        nodeType[bodyCenterNode] = static_cast<int8>(std::min(nds, static_cast<int>(m3c_node_type::k_QuadPoint)));
      }

      if(cubeFlag != 1 || nFE <= 2)
      {
        continue;
      }

      std::vector<SiteIdType> arrayFE(nFE);
      int tindex = 0;
      for(int i1 = 0; i1 < 6; i1++)
      {
        const Face& sqf = window[winIndex(sqID[i1])];
        const int tnfe = static_cast<int>(static_cast<uint8>(sqf.NEdge));
        for(int i2 = 0; i2 < tnfe; i2++)
        {
          arrayFE[tindex] = sqf.EdgeId[i2];
          tindex++;
        }
      }

      if(!generate)
      {
        // The first pass counts triangles and applies whole-volume edge flips.
        if(nFC == 0)
        {
          nTriangle += getNumberCase0Triangles(arrayFE, fedges, nFE);
        }
        else if(nFC == 2)
        {
          nTriangle += getNumberCase2Triangles(arrayFE, fedges, nFE, arrayFC, nFC);
        }
        else if(nFC > 2 && nFC <= 6)
        {
          nTriangle += getNumberCaseMTriangles(arrayFE, fedges, nFE, arrayFC, nFC);
        }
        continue;
      }

      // The second pass writes triangles to the pre-sized arrays in cube order.
      std::array<double, 3> coord1{};
      std::array<double, 3> coord2{};
      for(int k = 0; k < 3; k++)
      {
        coord1[k] = siteCoords[i].Coord[k];
        coord2[k] = siteCoords[i + 1 + xDim + numSitesPerPlane].Coord[k];
      }
      const int64 tin = tidRun;
      int64 tout = tin;
      if(nFC == 0)
      {
        getCase0Triangles(triangles, mCubeID, arrayFE, nodeCoords, fedges, nFE, tin, tout, coord1, coord2, i);
      }
      else if(nFC == 2)
      {
        getCase2Triangles(triangles, mCubeID, arrayFE, nodeCoords, fedges, nFE, arrayFC, nFC, tin, tout, coord1, coord2, i);
      }
      else
      {
        getCaseMTriangles(triangles, mCubeID, arrayFE, nodeCoords, fedges, nFE, arrayFC, nFC, tin, tout, bodyCenterNode, coord1, coord2, i);
      }
      tidRun = tout;
    }
  };

  if(!parallel)
  {
    m_MessageHandler.sendInfoMessage("Sweeping z-slices (pass 1: face edges + triangle count)...");
    sweep(true, false);
    if(m_ShouldCancel)
    {
      return {};
    }
    triangles.resize(static_cast<usize>(nTriangle < 0 ? 0 : nTriangle));
    mCubeID.resize(static_cast<usize>(nTriangle < 0 ? 0 : nTriangle), 0);

    m_MessageHandler.sendInfoMessage("Sweeping z-slices (pass 2: generating triangles)...");
    sweep(false, true);
  }
  else
  {
    // The serial edge stage preserves vertices, labels, and node types. Parallel
    // cubes read shared squares and flip private edge copies. This avoids shared
    // mutation. Cross-cube flips are omitted, so triangulation can differ while
    // interfaces remain valid and watertight.
    const SiteIdType lastCube = numSites - numSitesPerPlane;
    const usize numCubes = (lastCube >= 1) ? static_cast<usize>(lastCube) : 0; // cubes are 1..lastCube

    // Each cube uses private edges. Counting sets body-center node types. Generation
    // writes triangles at the precomputed offset.
    auto perCube = [&](const SiteIdType cubeSite, const bool doGenerate, const int64 triOffset) -> int64 {
      std::array<SiteIdType, 6> sqID{};
      sqID[0] = 3 * (cubeSite - 1);
      sqID[1] = (3 * (cubeSite - 1)) + 1;
      sqID[2] = (3 * (cubeSite - 1)) + 2;
      sqID[3] = (3 * cubeSite) + 2;
      sqID[4] = (3 * (cubeSite + xDim - 1)) + 1;
      sqID[5] = 3 * (cubeSite + numSitesPerPlane - 1);
      std::array<SiteIdType, 6> arrayFC{};
      for(int ii = 0; ii < 6; ii++)
      {
        arrayFC[ii] = -1;
      }
      int fcid = 0;
      int nFE = 0;
      int eff = 0;
      for(int ii = 0; ii < 6; ii++)
      {
        const Face& sqf = window[winIndex(sqID[ii])];
        if(sqf.FaceCenterNode != -1)
        {
          arrayFC[fcid] = sqf.FaceCenterNode;
          fcid++;
        }
        nFE = nFE + sqf.NEdge;
        eff = eff + sqf.Effect;
      }
      const int nFC = fcid;
      const SiteIdType bodyCenterNode = (7 * (cubeSite - 1)) + 6;
      if(!doGenerate && nFC >= 3)
      {
        const std::array<SiteIdType, 4> corners1 = squareCorners(sqID[0], neighbors);
        const std::array<SiteIdType, 4> corners2 = squareCorners(sqID[5], neighbors);
        std::array<int, 8> arraySpin{};
        for(int j = 0; j < 4; j++)
        {
          arraySpin[j] = point[corners1[j]];
          arraySpin[j + 4] = point[corners2[j]];
        }
        int nds = 0;
        for(int k = 0; k < 8; k++)
        {
          const int cspin = arraySpin[k];
          if(cspin != -1)
          {
            nds++;
            arraySpin[k] = -1;
            for(int kk = 0; kk < 8; kk++)
            {
              if(cspin == arraySpin[kk])
              {
                arraySpin[kk] = -1;
              }
            }
          }
        }
        // NodeType uses k_QuadPoint for four or more labels.
        nodeType[bodyCenterNode] = static_cast<int8>(std::min(nds, static_cast<int>(m3c_node_type::k_QuadPoint)));
      }
      if(eff <= 0 || nFE <= 2)
      {
        return 0;
      }
      // Private face-edge indexes start at zero. A cube has at most 24 edges,
      // and the fixed local buffer holds 64.
      std::array<SiteIdType, 64> localAFE{};
      std::array<Segment, 64> localEdges{};
      int tindex = 0;
      for(int i1 = 0; i1 < 6; i1++)
      {
        const Face& sqf = window[winIndex(sqID[i1])];
        const int tnfe = static_cast<int>(static_cast<uint8>(sqf.NEdge));
        for(int i2 = 0; i2 < tnfe; i2++)
        {
          localEdges[static_cast<usize>(tindex)] = fedges[sqf.EdgeId[i2]];
          localAFE[static_cast<usize>(tindex)] = tindex;
          tindex++;
        }
      }
      if(!doGenerate)
      {
        if(nFC == 0)
        {
          return getNumberCase0Triangles(localAFE, localEdges, nFE);
        }
        if(nFC == 2)
        {
          return getNumberCase2Triangles(localAFE, localEdges, nFE, arrayFC, nFC);
        }
        if(nFC > 2 && nFC <= 6)
        {
          return getNumberCaseMTriangles(localAFE, localEdges, nFE, arrayFC, nFC);
        }
        return 0;
      }
      std::array<double, 3> coord1{};
      std::array<double, 3> coord2{};
      for(int k = 0; k < 3; k++)
      {
        coord1[k] = siteCoords[cubeSite].Coord[k];
        coord2[k] = siteCoords[cubeSite + 1 + xDim + numSitesPerPlane].Coord[k];
      }
      const int64 tin = triOffset;
      int64 tout = tin;
      if(nFC == 0)
      {
        getCase0Triangles(triangles, mCubeID, localAFE, nodeCoords, localEdges, nFE, tin, tout, coord1, coord2, cubeSite);
      }
      else if(nFC == 2)
      {
        getCase2Triangles(triangles, mCubeID, localAFE, nodeCoords, localEdges, nFE, arrayFC, nFC, tin, tout, coord1, coord2, cubeSite);
      }
      else
      {
        getCaseMTriangles(triangles, mCubeID, localAFE, nodeCoords, localEdges, nFE, arrayFC, nFC, tin, tout, bodyCenterNode, coord1, coord2, cubeSite);
      }
      return tout - triOffset;
    };

    // Slide the serial edge window until it covers the target and next site planes.
    auto advanceWindowTo = [&](const SiteIdType targetBaseSite, const bool appendEdges, int64& eid) {
      while(winBaseSite < targetBaseSite)
      {
        std::memmove(window.data(), window.data() + sliceSquares, static_cast<usize>(sliceSquares) * sizeof(Face));
        winBaseSite += numSitesPerPlane;
        const SiteIdType newLoSquare = 3 * (winBaseSite + numSitesPerPlane - 1);
        const SiteIdType newHiSquare = std::min<SiteIdType>(3 * (winBaseSite + (2 * numSitesPerPlane) - 1), 3 * numSites);
        if(newLoSquare < newHiSquare)
        {
          computeSquares(newLoSquare, newHiSquare, appendEdges, eid);
        }
        if(appendEdges)
        {
          maybeReserveFedges(eid);
        }
      }
    };

    std::vector<int64> triOffset(numCubes + 1, 0); // 1-based per-cube; pass 1 fills counts, then prefix -> offsets

    m_MessageHandler.sendInfoMessage("Sweeping z-slices (parallel pass 1: counting)...");
    {
      winBaseSite = 1;
      int64 eid = 0;
      computeSquares(0, std::min<SiteIdType>(2 * sliceSquares, 3 * numSites), true, eid);
      for(SiteIdType sliceBase = 1; sliceBase <= lastCube; sliceBase += numSitesPerPlane)
      {
        if(m_ShouldCancel)
        {
          return {};
        }
        advanceWindowTo(sliceBase, true, eid);
        const SiteIdType cubeEnd = std::min<SiteIdType>(sliceBase + numSitesPerPlane, lastCube + 1);
        ParallelDataAlgorithm alg;
        alg.setRange(static_cast<usize>(sliceBase), static_cast<usize>(cubeEnd));
        alg.execute([&](const Range& range) {
          for(usize idx = range.min(); idx < range.max(); idx++)
          {
            triOffset[idx] = perCube(static_cast<SiteIdType>(idx), false, 0);
          }
        });
      }
    }

    int64 total = 0;
    for(usize i = 1; i <= numCubes; i++)
    {
      const int64 cnt = triOffset[i];
      triOffset[i] = total;
      total += cnt;
    }
    triangles.resize(static_cast<usize>(total));
    mCubeID.resize(static_cast<usize>(total), 0);

    m_MessageHandler.sendInfoMessage("Sweeping z-slices (parallel pass 2: generating triangles)...");
    {
      winBaseSite = 1;
      int64 eid = 0;
      computeSquares(0, std::min<SiteIdType>(2 * sliceSquares, 3 * numSites), false, eid);
      for(SiteIdType sliceBase = 1; sliceBase <= lastCube; sliceBase += numSitesPerPlane)
      {
        if(m_ShouldCancel)
        {
          return {};
        }
        advanceWindowTo(sliceBase, false, eid);
        const SiteIdType cubeEnd = std::min<SiteIdType>(sliceBase + numSitesPerPlane, lastCube + 1);
        ParallelDataAlgorithm alg;
        alg.setRange(static_cast<usize>(sliceBase), static_cast<usize>(cubeEnd));
        alg.execute([&](const Range& range) {
          for(usize idx = range.min(); idx < range.max(); idx++)
          {
            perCube(static_cast<SiteIdType>(idx), true, triOffset[idx]);
          }
        });
      }
    }
  }

  if(m_ShouldCancel)
  {
    return {};
  }

  return finalizeMesh(m_DataStructure, m_InputValues, m_MessageHandler, m_ShouldCancel, triangles, mCubeID, fedges, nodeType, point, nodeCoords, neighbors, numSites, fileDim, dims, maxGrainId);
}
} // namespace nx::core
