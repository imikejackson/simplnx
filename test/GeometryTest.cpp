#include "GeometryTestUtilities.hpp"

#include "complex/DataStructure/DataStore.hpp"
#include "complex/DataStructure/DataStructure.hpp"
#include "complex/DataStructure/Geometry/EdgeGeom.hpp"
#include "complex/DataStructure/Geometry/HexahedralGeom.hpp"
#include "complex/DataStructure/Geometry/ImageGeom.hpp"
#include "complex/DataStructure/Geometry/QuadGeom.hpp"
#include "complex/DataStructure/Geometry/RectGridGeom.hpp"
#include "complex/DataStructure/Geometry/TetrahedralGeom.hpp"
#include "complex/DataStructure/Geometry/TriangleGeom.hpp"
#include "complex/DataStructure/Geometry/VertexGeom.hpp"
#include "complex/Utilities/Math/GeometryMath.hpp"

#include <catch2/catch.hpp>

using namespace complex;

////////////////////////////////////
// Begin generic geometry testing //
////////////////////////////////////
void testAbstractGeometry(IGeometry* geom)
{
  SECTION("abstract geometry")
  {
    SECTION("units")
    {
      const auto units = IGeometry::LengthUnit::Fathom;
      geom->setUnits(units);
      REQUIRE(geom->getUnits() == units);
    }
    SECTION("dimensionality")
    {
      const uint32 uDims = 20;
      geom->setUnitDimensionality(uDims);
      REQUIRE(geom->getUnitDimensionality() == uDims);

      const uint32 sDims = 14;
      geom->setSpatialDimensionality(sDims);
      REQUIRE(geom->getSpatialDimensionality() == sDims);
    }
  }
}

void testGeom2D(INodeGeometry2D* geom)
{
  SECTION("abstract geometry 2D")
  {
    const usize vertId = 2;
    const Point3D<float32> coord = {0.5f, 0.0f, 2.0f};

    // Vertices
    {
      auto vertices = createVertexList(geom);
      REQUIRE(vertices != nullptr);
      geom->setVertices(*vertices);
      REQUIRE(geom->getVertices() == vertices);
      const usize numVertices = 10;
      geom->resizeVertexList(numVertices);
      REQUIRE(geom->getNumberOfVertices() == numVertices);

      geom->setVertexCoordinate(vertId, coord);
      REQUIRE(geom->getVertexCoordinate(vertId) == coord);
    }

    // edges
    {
      auto edges = createEdgeList(geom);
      geom->setEdgeList(*edges);
      REQUIRE(geom->getEdges() == edges);
      const usize numEdges = 5;
      geom->resizeEdgeList(numEdges);
      REQUIRE(geom->getNumberOfEdges() == numEdges);
      const usize edgeId = 3;
      std::array<usize, 2> verts = {vertId, vertId + 1};
      geom->setEdgePointIds(edgeId, verts);

      std::array<Point3Df, 2> edge_verts;
      std::array<usize, 2> vertsOut = {0, 0};
      geom->getEdgePointIds(edgeId, vertsOut);
      for(usize i = 0; i < 2; i++)
      {
        REQUIRE(verts[i] == vertsOut[i]);
      }
      geom->getEdgeCoordinates(edgeId, edge_verts);
      REQUIRE(edge_verts[0] == coord);
    }
  }
}

void testGeom3D(INodeGeometry3D* geom)
{
  SECTION("abstract geometry 3D")
  {
    const usize vertId = 2;
    const Point3D<float32> coord = {0.5f, 0.0f, 2.0f};

    // vertices
    {
      auto vertices = createVertexList(geom);
      geom->setVertices(*vertices);
      REQUIRE(geom->getVertices() == vertices);
      const usize numVertices = 10;
      geom->resizeVertexList(numVertices);
      REQUIRE(geom->getNumberOfVertices() == numVertices);

      geom->setVertexCoordinate(vertId, coord);
      REQUIRE(geom->getVertexCoordinate(vertId) == coord);
    }
    // edges
    {
      auto edges = createEdgeList(geom);
      geom->setEdgeList(*edges);
      REQUIRE(geom->getEdges() == edges);
      const usize numEdges = 5;
      geom->resizeEdgeList(numEdges);
      REQUIRE(geom->getNumberOfEdges() == numEdges);
      const usize edgeId = 3;
      std::array<usize, 2> verts = {vertId, vertId + 1};
      geom->setEdgePointIds(edgeId, verts);
      usize vertsOut[2];
      geom->getEdgePointIds(edgeId, vertsOut);
      for(usize i = 0; i < 2; i++)
      {
        REQUIRE(verts[i] == vertsOut[i]);
      }
      std::array<Point3Df, 2> edge_verts;

      geom->getEdgeCoordinates(edgeId, edge_verts);
      REQUIRE(edge_verts[0] == coord);
    }
    // faces
    {
    }
  }
}

void testGeomGrid(IGridGeometry* geom)
{
  SECTION("abstract geometry grid")
  {
    const usize xDim = 10;
    const usize yDim = 150;
    const usize zDim = 50;
    const SizeVec3 dims = {xDim, yDim, zDim};
    geom->setDimensions(dims);
    REQUIRE(geom->getDimensions() == dims);

    REQUIRE(geom->getNumXCells() == xDim);
    REQUIRE(geom->getNumYCells() == yDim);
    REQUIRE(geom->getNumZCells() == zDim);
  }
}

/////////////////////////////////////
// Begin geometry-specific testing //
/////////////////////////////////////
TEST_CASE("EdgeGeomTest")
{
  DataStructure dataStructure;
  auto geom = createGeom<EdgeGeom>(dataStructure);

  testAbstractGeometry(geom);

  SECTION("type as string")
  {
    REQUIRE(geom->getTypeName() == "EdgeGeom");
  }
}

TEST_CASE("HexahedralGeomTest")
{
  DataStructure dataStructure;
  auto geom = createGeom<HexahedralGeom>(dataStructure);

  testGeom3D(geom);

  SECTION("type as string")
  {
    REQUIRE(geom->getTypeName() == "HexahedralGeom");
  }
}

TEST_CASE("ImageGeomTest")
{
  DataStructure dataStructure;
  auto geom = createGeom<ImageGeom>(dataStructure);

  testGeomGrid(geom);

  SECTION("type as string")
  {
    REQUIRE(geom->getTypeName() == "ImageGeom");
  }
}

TEST_CASE("QuadGeomTest")
{
  DataStructure dataStructure;
  auto geom = createGeom<QuadGeom>(dataStructure);

  testGeom2D(geom);

  SECTION("type as string")
  {
    REQUIRE(geom->getTypeName() == "QuadGeom");
  }
}

TEST_CASE("RectGridGeomTest")
{
  DataStructure dataStructure;
  auto geom = createGeom<RectGridGeom>(dataStructure);

  testGeomGrid(geom);

  SECTION("type as string")
  {
    REQUIRE(geom->getTypeName() == "RectGridGeom");
  }
}

TEST_CASE("TetrahedralGeomTest")
{
  DataStructure dataStructure;
  auto geom = createGeom<TetrahedralGeom>(dataStructure);

  testGeom3D(geom);

  SECTION("type as string")
  {
    REQUIRE(geom->getTypeName() == "TetrahedralGeom");
  }
}

TEST_CASE("TriangleGeomTest")
{
  DataStructure dataStructure;
  auto geom = createGeom<TriangleGeom>(dataStructure);

  testGeom2D(geom);

  SECTION("type as string")
  {
    REQUIRE(geom->getTypeName() == "TriangleGeom");
  }
}

TEST_CASE("VertexGeomTest")
{
  DataStructure dataStructure;
  auto geom = createGeom<VertexGeom>(dataStructure);

  testAbstractGeometry(geom);

  SECTION("type as string")
  {
    REQUIRE(geom->getTypeName() == "VertexGeom");
  }
}


int RayTriangleIntersect(const FloatVec3& orig, const FloatVec3& dir, const FloatVec3& vert0, const FloatVec3& vert1, const FloatVec3& vert2, FloatVec3& bcoord)
{
  constexpr float k_Epsilon = 1e-8;

  FloatVec3 edge1, edge2, tvec, pvec, qvec;
  double det, inv_det;

  /* find vectors for two edges sharing vert0 */
  edge1 = vert1 - vert0;
  edge2 = vert2 - vert0;

  /* begin calculating determinant - also used to calculate U parameter */
  pvec = dir.cross(edge2);

  /* if determinant is near zero, ray lies in plane of triangle */
  det = edge1.dot(pvec);

  /* calculate distance from vert0 to ray origin */
  tvec = orig - vert0;

  inv_det = 1.0 / det;

  qvec = tvec.cross(edge1);

  if(det > k_Epsilon)
  {
    bcoord[0] = tvec.dot(pvec);
    if(bcoord[0] < 0.0 || bcoord[0] > det)
    {
      return 0;
    }

    /* calculate V parameter and test bounds */
    // bcoord[1] = DOT(dir, qvec);
    bcoord[1] = dir.dot(qvec);
    if(bcoord[1] < 0.0 || bcoord[0] + bcoord[1] > det)
    {
      return 0;
    }
  }
  else if(det < -k_Epsilon)
  {
    /* calculate U parameter and test bounds */
    bcoord[0] = tvec.dot(pvec);
    if(bcoord[0] > 0.0 || bcoord[0] < det)
    {
      return 0;
    }

    /* calculate V parameter and test bounds */
    bcoord[1] = dir.dot(qvec);
    if(bcoord[1] > 0.0 || bcoord[0] + bcoord[1] < det)
    {
      return 0;
    }
  }
  else
  {
    return false;
  } /* ray is 0 to the plane of the triangle */

  // float t = edge2.dotProduct(qvec) * inv_det;
  float u = bcoord[0] * inv_det;
  float v = bcoord[1] * inv_det;

  bcoord[0] = 1 - u - v;
  bcoord[1] = u;
  bcoord[2] = v;

  return 1;
}

#define EPSILON 0.000001
#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2];

/* code rewritten to do tests on the sign of the determinant */
/* the division is at the end in the code                    */
int intersect_triangle1(const FloatVec3& orig, const FloatVec3& dir, const FloatVec3& vert0, const FloatVec3& vert1, const FloatVec3& vert2, FloatVec3& bcoord)
{
  double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  double det,inv_det;

  /* find vectors for two edges sharing vert0 */
  SUB(edge1, vert1, vert0);
  SUB(edge2, vert2, vert0);

  /* begin calculating determinant - also used to calculate U parameter */
  CROSS(pvec, dir, edge2);

  /* if determinant is near zero, ray lies in plane of triangle */
  det = DOT(edge1, pvec);
  std::cout << "  det: " << det << "\n";
  if (det > EPSILON)
  {
    /* calculate distance from vert0 to ray origin */
    SUB(tvec, orig, vert0);

    /* calculate U parameter and test bounds */
    bcoord[1] = DOT(tvec, pvec);
    if (bcoord[1] < 0.0 || bcoord[1] > det)
      return 0;

    /* prepare to test V parameter */
    CROSS(qvec, tvec, edge1);

    /* calculate V parameter and test bounds */
    bcoord[2] = DOT(dir, qvec);
    if (bcoord[2] < 0.0 || bcoord[1] + bcoord[2] > det)
      return 0;

  }
  else if(det < -EPSILON)
  {
    /* calculate distance from vert0 to ray origin */
    SUB(tvec, orig, vert0);

    /* calculate U parameter and test bounds */
    bcoord[1] = DOT(tvec, pvec);
    /*      printf("*u=%f\n",(float)*u); */
    /*      printf("det=%f\n",det); */
    if (bcoord[1] > 0.0 || bcoord[1] < det)
      return 0;

    /* prepare to test V parameter */
    CROSS(qvec, tvec, edge1);

    /* calculate V parameter and test bounds */
    bcoord[2] = DOT(dir, qvec) ;
    if (bcoord[2] > 0.0 || bcoord[1] + bcoord[2] < det)
      return 0;
  }
  else return 0;  /* ray is parallel to the plane of the triangle */


  inv_det = 1.0 / det;

  /* calculate t, ray intersects triangle */
  bcoord[0] = DOT(edge2, qvec) * inv_det;
  bcoord[1] *= inv_det;
  bcoord[2] *= inv_det;

  return 1;
}


TEST_CASE("RayIntersects Triangle")
{
  // Start Slightly below the Z Plane and send the ray TOWARDS the +Z direction
  FloatVec3 rayOrigin(0.0f, 0.0f, -1.0f);
  FloatVec3 rayDirection(0.0f, 0.0f, 1.0f); //Happens to be normalized :-)

  FloatVec3 v0(1.0f, 1.0f, 0.0f);
  FloatVec3 v1(-1.0f, 1.0f, 0.0f);
  FloatVec3 v2(0.0f, -1.0f, 0.0f);

  FloatVec3 bcoords(0.0f, 0.0f, 0.0f);

  int intersects = intersect_triangle1(rayOrigin, rayDirection, v0, v1, v2, bcoords);
  std::cout << intersects << fmt::format("  {}", fmt::join(bcoords, ",")) << "\n";

  rayDirection = {0.0f, 0.0f, -1.0f}; // Flip the ray
  intersects = intersect_triangle1(rayOrigin, rayDirection, v0, v1, v2, bcoords);
  std::cout << intersects << fmt::format("  {}", fmt::join(bcoords, ",")) << "\n";


//  REQUIRE(intersects == 1);
//
//  // Start above the XY plane, this should NOT intersect
//  rayOrigin = {0.0f, 0.0f, 10.0f};
//  intersects = intersect_triangle1(rayOrigin, rayDirection, v0, v1, v2, bcoords);
//  std::cout << fmt::format("{}", fmt::join(bcoords, ",")) << "\n";
//
//  REQUIRE(intersects == 0);

}
