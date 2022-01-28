/**
MIT License

Copyright (c) 2021 Jaswant Sai Panchumarti

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "tscTriSurfaceCutter.h"

#include <algorithm>
#include <array>
#include <set>
#include <unordered_map>
#include <vector>

#include <vtkArrayDispatch.h>
#include <vtkBoundingBox.h>
#include <vtkCellArray.h>
#include <vtkCellArrayIterator.h>
#include <vtkCellData.h>
#include <vtkCellIterator.h>
#include <vtkCellLocator.h>
#include <vtkCellType.h>
#include <vtkDataArray.h>
#include <vtkDataArrayRange.h>
#include <vtkDelaunay2D.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLine.h>
#include <vtkMath.h>
#include <vtkMathUtilities.h>
#include <vtkMergePoints.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkStaticCellLocator.h>
#include <vtkTriangle.h>

#define TSC_MAX_EDGE_FLIPS 32 // Hard limit to prevent recursion caused stack overflow.

vtkStandardNewMacro(tscTriSurfaceCutter);

tscTriSurfaceCutter::tscTriSurfaceCutter()
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
  this->CreateDefaultLocators();
}

tscTriSurfaceCutter::~tscTriSurfaceCutter() {}

void tscTriSurfaceCutter::SetCuttersData(vtkPolyData* cutters)
{
  this->SetInputData(1, cutters);
}

void tscTriSurfaceCutter::SetCuttersConnection(vtkAlgorithmOutput* output)
{
  this->SetInputConnection(1, output);
}

int tscTriSurfaceCutter::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

int tscTriSurfaceCutter::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

namespace
{
  // Convention: in a triangle - v: vertex, e: edge
  // with vertices v0--v1--v2,
  // e0 = v0--v1, e1 = v1--v2, e2 = v2--v0;
  const std::array<std::pair<vtkIdType, vtkIdType>, 3> TRIEDGES = { std::make_pair(0, 1),
    std::make_pair(1, 2), std::make_pair(2, 0) };
  const vtkIdType TRIOPPEDGE[3] = { 1, 2, 0 };  // edge id opposite to a vertex
  const vtkIdType TRIOPPVERTS[3] = { 2, 0, 1 }; // vertex id opposite to an edge

  enum class PointTrianglePosition
  {
    OnVertex,
    OnEdge,
    Inside,
    Outside,
    Degenerate
  };

  enum class PointLinePosition
  {
    OnVertex,
    Inside,
    Outside,
  };

  enum class IntersectType
  {
    NoIntersection = 0,
    Intersect,
    Junction
  };

  // This information pertains to the cutter polygon and is used quite frequently.
  struct CutterInfo
  {
    vtkBoundingBox bbox;
    vtkNew<vtkDoubleArray> points;
    int cell_type = VTK_POLYGON;
    double normal[3] = {};
  };
  using CuttersInfoType = std::vector<CutterInfo>;
  using dispatchRR =
    vtkArrayDispatch::Dispatch2ByValueType<vtkArrayDispatch::Reals, vtkArrayDispatch::Reals>;
  using dispatchR = vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals>;

  using SegmentType = std::pair<vtkIdType, vtkIdType>;
  using SegmentsType = std::vector<SegmentType>;

  /**
   * @brief When a triangle is cut by a cutter polygon, it births
   *         children triangles situated around the polygon's edge that cut said triangle.
   */
  struct Child
  {
    Child()
      : cx(0)
      , cy(0)
      , bbox()
    {
    }
    Child(const double& cx_, const double& cy_, const vtkIdType* pts, const vtkIdType npts_,
      const double bounds[4])
      : cx(cx_)
      , cy(cy_)
    {
      point_ids->SetNumberOfIds(npts_);
      std::copy(pts, pts + npts_, point_ids->GetPointer(0));
      bbox = vtkBoundingBox(bounds[0], bounds[1], bounds[2], bounds[3], 0.0, 0.0);
    }
    Child(const double& cx_, const double& cy_, vtkIdList* pts, const double bounds[4])
      : cx(cx_)
      , cy(cy_)
    {
      point_ids->DeepCopy(pts);
      bbox = vtkBoundingBox(bounds[0], bounds[1], bounds[2], bounds[3], 0.0, 0.0);
    }
    double cx, cy;               // centroid
    vtkNew<vtkIdList> point_ids; // the ids
    vtkBoundingBox bbox;         // bbox of all points
  };

  struct GetChildrenImpl
  {
    /**
     * @brief Implementation to populate children triangles.
     * @warning This centroid code won't fly for polygons. (Since we know each child is a triangle,
     * we're safe)
     * @tparam PointsArrT T
     * @param points_arr  coordinates that make up parent triangle
     * @param polys       point ids that make up children triangles
     * @param children    this will have the child triangles with bounds and centroids.
     */
    template <typename PointsArrT, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT>>
    void operator()(PointsArrT* points_arr, vtkCellArray* polys, std::vector<Child>& children)
    {
      auto points = vtk::DataArrayTupleRange<3>(points_arr);

      const vtkIdType& num_polys = polys->GetNumberOfCells();
      for (vtkIdType i = 0; i < num_polys; ++i)
      {
        const vtkIdType* point_ids = nullptr;
        vtkIdType num_point_ids{ 0 };
        polys->GetCellAtId(i, num_point_ids, point_ids);

        // a child requires it's bounds and centroid
        double bounds[4] = { VTK_DOUBLE_MAX, VTK_DOUBLE_MIN, VTK_DOUBLE_MAX, VTK_DOUBLE_MIN };
        double pc[3] = {};
        for (vtkIdType v = 0; v < num_point_ids; ++v)
        {
          const vtkIdType& pt = point_ids[v];
          const double x = points[pt][0];
          const double y = points[pt][1];

          bounds[0] = (x < bounds[0]) ? x : bounds[0];
          bounds[1] = (x > bounds[1]) ? x : bounds[1];
          bounds[2] = (y < bounds[2]) ? y : bounds[2];
          bounds[3] = (y > bounds[3]) ? y : bounds[3];

          for (unsigned short dim = 0; dim < 2; ++dim)
          {
            pc[dim] += points[pt][dim];
          }
        }

        for (unsigned short dim = 0; dim < 2; ++dim)
        {
          pc[dim] /= double(num_point_ids);
        }

        children.emplace_back(pc[0], pc[1], point_ids, num_point_ids, bounds);
      }
    }
  };

  /**
   * @brief When a root triangle(line) intersects a line segment, it births child
   * triangles(lines).
   */
  class Parent : public vtkGenericCell
  {
  public:
    static Parent* New() { VTK_STANDARD_NEW_BODY(Parent); }
    vtkTypeMacro(Parent, vtkGenericCell);
    void Reset()
    {
      this->children.clear();
      this->PointIds->Reset();
      this->Points->Reset();
    }
    void UpdateChildren(vtkCellArray* polys)
    {
      GetChildrenImpl worker;
      vtkDataArray* points_arr = this->Points->GetData();
      if (!dispatchR::Execute(points_arr, worker, polys, this->children))
      {
        worker(points_arr, polys, this->children);
      }
    }

    std::vector<Child> children;
    vtkIdType cellId;

  protected:
    Parent()
      : children()
      , cellId(-1)
    {
      this->children.reserve(10);
    }
    ~Parent() override = default;

  private:
    Parent(const Parent&) = delete;
    void operator=(const Parent&) = delete;
  };

  // Determine a point's location in a triangle.
  // is it 1. on an edge of a triangle?
  //       2. inside/outside ?
  //       3. exactly coincident with a vertex?
  //       4. impossible to determine (triangle is degenerate)
  PointTrianglePosition inTriangle(const double p[3], const double p0[3], const double p1[3],
    const double p2[3], double bary_coords[3], const double& tol, vtkIdType& v, vtkIdType& e)
  {
    v = -1;
    e = -1;
    if (!vtkTriangle::BarycentricCoords(p, p0, p1, p2, bary_coords))
      return PointTrianglePosition::Degenerate;

    const double& s = bary_coords[0];
    const double& t = bary_coords[1];
    const double& st = bary_coords[2]; // 1 - s - t

    bool s_eq_0 =
      vtkMathUtilities::FuzzyCompare(s, 0., tol) || vtkMathUtilities::NearlyEqual(s, 0., tol);
    bool t_eq_0 =
      vtkMathUtilities::FuzzyCompare(t, 0., tol) || vtkMathUtilities::NearlyEqual(t, 0., tol);
    bool st_eq_0 =
      vtkMathUtilities::FuzzyCompare(st, 0., tol) || vtkMathUtilities::NearlyEqual(st, 0., tol);

    bool s_eq_1 = vtkMathUtilities::NearlyEqual(s, 1., tol);
    bool t_eq_1 = vtkMathUtilities::NearlyEqual(t, 1., tol);
    bool st_eq_1 = vtkMathUtilities::NearlyEqual(st, 1., tol);

    if (s_eq_0 && t_eq_0 && s_eq_1)
    {
      v = 0;
    }
    else if (s_eq_0 && st_eq_0 && t_eq_1)
    {
      v = 1;
    }
    else if (t_eq_0 && st_eq_0 && st_eq_1)
    {
      v = 2;
    }
    if (v >= 0)
    {
      return PointTrianglePosition::OnVertex;
    }

    // if all barycentric coordinates are in open interval (0, 1), then point is surely inside
    bool s_gt0_lt1 = !std::signbit(s) && std::isless(s, 1) && !s_eq_0;
    bool t_gt0_lt1 = !std::signbit(t) && std::isless(t, 1) && !t_eq_0;
    bool st_gt0_lt1 = !std::signbit(st) && std::isless(st, 1) && !st_eq_0;
    if (s_gt0_lt1 && t_gt0_lt1 && st_gt0_lt1)
    {
      return PointTrianglePosition::Inside;
    }

    // Either almost on an edge (within tolerance limits)
    if (s_eq_0 && t_gt0_lt1 && st_gt0_lt1)
    {
      e = TRIOPPEDGE[0];
    }
    else if (t_eq_0 && s_gt0_lt1 && st_gt0_lt1)
    {
      e = TRIOPPEDGE[1];
    }
    else if (st_eq_0 && s_gt0_lt1 && t_gt0_lt1)
    {
      e = TRIOPPEDGE[2];
    }
    if (e >= 0)
    {
      return PointTrianglePosition::OnEdge;
    }

    // or the triangle is degenerate.
    if (s_eq_0 && t_eq_0 && st_eq_0)
    {
      return PointTrianglePosition::Degenerate;
    }

    // or the point is just outside.
    return PointTrianglePosition::Outside;
  }

  // Determine a point's location on a line segment.
  // is it 1. on one of the end vertices?
  //       2. inside/outside ?
  PointLinePosition onLine(
    const double p[3], const double p1[3], const double p2[3], const double& tol)
  {
    double t{ 0 };
    const double dist2 = vtkLine::DistanceToLine(p, p1, p2, t, nullptr);

    if (dist2 > tol * tol)
    {
      return PointLinePosition::Outside;
    }

    const double& s = 1 - t;

    bool s_eq_0 =
      vtkMathUtilities::FuzzyCompare(s, 0., tol) || vtkMathUtilities::NearlyEqual(s, 0., tol);
    bool t_eq_0 =
      vtkMathUtilities::FuzzyCompare(t, 0., tol) || vtkMathUtilities::NearlyEqual(t, 0., tol);

    bool s_eq_1 = vtkMathUtilities::NearlyEqual(s, 1., tol);
    bool t_eq_1 = vtkMathUtilities::NearlyEqual(t, 1., tol);

    if ((s_eq_0 && t_eq_1) || (s_eq_1 && t_eq_0))
    {
      return PointLinePosition::OnVertex;
    }

    bool s_gt0_lt1 = !std::signbit(s) && std::isless(s, 1) && !s_eq_0;
    bool t_gt0_lt1 = !std::signbit(t) && std::isless(t, 1) && !t_eq_0;

    if (s_gt0_lt1 && t_gt0_lt1)
      return PointLinePosition::Inside;

    return PointLinePosition::Outside;
  }

  // Reasonable z-value of a point inside(or on an edge/vertex) a triangle.
  void interpZ(double p[3], const double p0[3], const double p1[3], const double p2[3])
  {
    double bary_coords[3] = {};
    vtkTriangle::BarycentricCoords(p, p0, p1, p2, bary_coords);
    p[2] = bary_coords[0] * p0[2] + bary_coords[1] * p1[2] + bary_coords[2] * p2[2];
  }

  // Same as above, but when you already have bary_coords
  void interpZ(
    double p[3], const double p0[3], const double p1[3], const double p2[3], double bary_coords[3])
  {
    p[2] = bary_coords[0] * p0[2] + bary_coords[1] * p1[2] + bary_coords[2] * p2[2];
  }

  /**
   * @brief Quantify return status of vtkLine::Intersection(...)
   * This test checks for vertex and edge junctions.
   * For example
   *  Vertex-vertex junction
   *    (u=0 v=0), (u=0 v=1), (u=1 v=0), (u=1 v=0)
   *  Edge-(edge or vertex) juntions
   *    (u=0 v!=0 v!=1), (u=1 v!=0 v!=1)
   *    (u!=0 u!=1 v=0), (u!=0 u!=1 v=1)
   */
  IntersectType robustIntersect(
    double* p1, double* p2, double* q1, double* q2, double px[3], const double& tol)
  {
    double u(0.0), v(0.0);
    const auto vtk_intersect = vtkLine::Intersection(p1, p2, q1, q2, u, v, tol, vtkLine::Absolute);

    if (vtk_intersect == vtkLine::Intersect)
    {
      if ((vtkMathUtilities::NearlyEqual(u, 0.0) || vtkMathUtilities::NearlyEqual(u, 1.0)) ||
        (vtkMathUtilities::NearlyEqual(v, 0.0) || vtkMathUtilities::NearlyEqual(v, 1.0)))
      {
        return IntersectType::Junction;
      }
      else
      {
        for (unsigned short dim = 0; dim < 2; ++dim)
        {
          px[dim] = p1[dim] + u * (p2[dim] - p1[dim]);
        }
        return IntersectType::Intersect;
      }
    }
    else if (vtk_intersect == vtkLine::OnLine)
    {
      double p21[3], c[3], c2[3];
      vtkMath::Subtract(p2, p1, p21);
      double tol2 = /*tolerance * tolerance*/ 1.0e-8 * sqrt(vtkMath::Dot(p21, p21));
      if (vtkLine::DistanceBetweenLines(p1, p2, q1, q2, c, c2, u, v) > tol2)
      {
        return IntersectType::NoIntersection;
      }
      vtkLine::DistanceToLine(p1, q1, q2, u, px);
      if (-tol <= u && u <= (1.0 + tol))
      {
        return IntersectType::Junction;
      }
      vtkLine::DistanceToLine(p2, q1, q2, u, px);
      if (-tol <= u && u <= (1.0 + tol))
      {
        return IntersectType::Junction;
      }
      vtkLine::DistanceToLine(q1, p1, p2, u, px);
      if (-tol <= u && u <= (1.0 + tol))
      {
        return IntersectType::Junction;
      }
      vtkLine::DistanceToLine(q2, p1, p2, u, px);
      if (-tol <= u && u <= (1.0 + tol))
      {
        return IntersectType::Junction;
      }
    }

    return IntersectType::NoIntersection;
  }

  //
  struct GetRootImpl
  {
    /**
     * @brief Implementation to copy points of a root triangle.
     * @tparam PointsArrT1    Dest T
     * @tparam PointsArrT2    Src T
     * @param points1_arr     root triangle points (dst)
     * @param points2_arr     dataset points (src)
     * @param root_point_ids  indices into points1_arr that make up root
     */
    template <typename PointsArrT1, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT1>,
      typename PointsArrT2, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT2>>
    void operator()(PointsArrT1* points1_arr, PointsArrT2* points2_arr, vtkIdList* root_point_ids)
    {
      auto points_1 = vtk::DataArrayTupleRange<3>(points1_arr);
      auto points_2 = vtk::DataArrayTupleRange<3>(points2_arr);

      for (vtkIdType i = 0; i < root_point_ids->GetNumberOfIds(); ++i)
      {
        const vtkIdType& pt = root_point_ids->GetId(i);
        std::copy(points_2[pt].begin(), points_2[pt].end(), points_1[i].begin());
      }
    }
  };

  struct TriIntersect2dImpl
  {
    /**
     * @brief   A brute-force implementation that intersects(in an aribtrary 2D plane @z=c) a
     * triangle with a bunch of line segments. The triangle and the line must lie in the same 2D
     * plane (zTri ~ zLine) It uses hints such as a point's location w.r.t a triangle to hold off
     * from actually testing intersection of line segments unless really needed. If this
     * implementation could be modified to work for a concave polygon, we can cut polygons too!
     * @tparam PointsArrT1  triangle's points' data array T
     * @tparam PointsArrT2  cutters' points' data array T
     * @param points1_arr   triangle's points' data array
     * @param points2_arr   cutters' points' data array
     * @param lines         point ids of line segments
     * @param constraints   list of segments that make up constraints
     * @param tol           numeric tolerance for intersection math.
     *
     */
    template <typename PointsArrT1, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT1>,
      typename PointsArrT2, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT1>>
    void operator()(PointsArrT1* points1_arr, PointsArrT2* points2_arr, const SegmentsType& lines,
      SegmentsType& constraints, const double& tol)
    {
      if (!lines.size())
        return;

      auto points_1 = vtk::DataArrayTupleRange<3>(points1_arr); // triangle's points
      auto points_2 = vtk::DataArrayTupleRange<3>(points2_arr); // cutter's points

      // clang-format off
      double tri_coords_2d[3][3] = 
      { { points_1[0][0], points_1[0][1], 0.0 },
        { points_1[1][0], points_1[1][1], 0.0 }, 
        { points_1[2][0], points_1[2][1], 0.0 } };  // use this for point in triangle, intersection tests.

      double tri_coords_3d[3][3] = 
      { { points_1[0][0], points_1[0][1], points_1[0][2] },
        { points_1[1][0], points_1[1][1], points_1[1][2] },
        { points_1[2][0], points_1[2][1], points_1[2][2] } }; // use this for z-interpolation
      // clang-format on

      std::unordered_map<vtkIdType, vtkIdType> processed; // k: line pt, v: insertLoc

      // a1, a2: line point ids.
      // p1, p2: line point coordinates.
      // b1, b2: triangle point ids per edge.
      // q1, q2: triangle coordinates per edge.
      for (const auto& line : lines)
      {
        const vtkIdType& a1 = line.first;
        const vtkIdType& a2 = line.second;

        double p1[3] = { points_2[a1][0], points_2[a1][1], 0.0 };
        double p2[3] = { points_2[a2][0], points_2[a2][1], 0.0 };

        vtkIdType a1v(-1), a2v(-1);
        vtkIdType a1e(-1), a2e(-1); // not used in current impl.

        double p1_bary_coords[3] = {};
        double p2_bary_coords[3] = {};

        const auto a1Pos = inTriangle(
          p1, tri_coords_2d[0], tri_coords_2d[1], tri_coords_2d[2], p1_bary_coords, tol, a1v, a1e);
        const auto a2Pos = inTriangle(
          p2, tri_coords_2d[0], tri_coords_2d[1], tri_coords_2d[2], p2_bary_coords, tol, a2v, a2e);

        // for convenience
        const bool a1_inside = a1Pos == PointTrianglePosition::Inside;
        const bool a2_inside = a2Pos == PointTrianglePosition::Inside;
        const bool a1_on_edge = a1Pos == PointTrianglePosition::OnEdge;
        const bool a2_on_edge = a2Pos == PointTrianglePosition::OnEdge;
        const bool a1_on_vert = a1Pos == PointTrianglePosition::OnVertex;
        const bool a2_on_vert = a2Pos == PointTrianglePosition::OnVertex;
        const bool a1_outside = a1Pos == PointTrianglePosition::Outside;
        const bool a2_outside = a2Pos == PointTrianglePosition::Outside;

        double px[3] = {}; // new inserted point.
        std::set<vtkIdType>
          hits; // store point ids. when elements are paired up, they form constraints.
        vtkIdType inserted(-1); // keeps location of px when inserted into points1_arr

        // here, try to avoid a call to `robustIntersect()` as much as possible.
        if (a1_inside || a1_on_edge)
        {
          if (processed.find(a1) == processed.end())
          {
            std::copy(p1, p1 + 3, px);
            interpZ(px, tri_coords_3d[0], tri_coords_3d[1], tri_coords_3d[2], p1_bary_coords);
            inserted = points1_arr->InsertNextTuple(px);
            processed[a1] = inserted;
          }
          hits.insert(processed[a1]);
        }

        if (a2_inside || a2_on_edge)
        {
          if (processed.find(a2) == processed.end())
          {
            std::copy(p2, p2 + 3, px);
            interpZ(px, tri_coords_3d[0], tri_coords_3d[1], tri_coords_3d[2], p2_bary_coords);
            inserted = points1_arr->InsertNextTuple(px);
            processed[a2] = inserted;
          }
          hits.insert(processed[a2]);
        }

        if (a1_on_vert)
        {
          if (processed.find(a1) == processed.end())
          {
            inserted = a1v;
            processed[a1] = inserted;
          }
          hits.insert(processed[a1]);
        }

        if (a2_on_vert)
        {
          if (processed.find(a2) == processed.end())
          {
            inserted = a2v;
            processed[a2] = inserted;
          }
          hits.insert(processed[a2]);
        }

        // Call robustIntersect()
        if ((a1_outside && a2_outside) || (a1_outside && !a2_outside) ||
          (!a1_outside && a2_outside))
        {
          for (const auto& edge : TRIEDGES)
          {
            const vtkIdType& b1 = edge.first;
            const vtkIdType& b2 = edge.second;

            double* q1 = tri_coords_2d[b1];
            double* q2 = tri_coords_2d[b2];

            const auto intersect_type = robustIntersect(p1, p2, q1, q2, px, tol);
            if (intersect_type == IntersectType::Intersect)
            {
              interpZ(px, tri_coords_3d[0], tri_coords_3d[1], tri_coords_3d[2]);
              inserted = points1_arr->InsertNextTuple(px);
              hits.insert(inserted);
            }
            else
            {
              // possibly on line or a junction.
              // junctions were handled before performing intersections (ax_on_edge or ax_on_vert)
              inserted = -1;
              const auto q1_on_p1_p2 = onLine(q1, p1, p2, tol);
              const auto q2_on_p1_p2 = onLine(q2, p1, p2, tol);
              const auto p1_on_q1_q2 = onLine(p1, q1, q2, tol);
              const auto p2_on_q1_q2 = onLine(p2, q1, q2, tol);
              if (q1_on_p1_p2 != PointLinePosition::Outside)
              {
                inserted = b1;
              }
              else if (q2_on_p1_p2 != PointLinePosition::Outside)
              {
                inserted = b2;
              }
              else if (p1_on_q1_q2 != PointLinePosition::Outside)
              {
                if (processed.find(a1) == processed.end())
                {
                  std::copy(p1, p1 + 3, px);
                  interpZ(px, tri_coords_3d[0], tri_coords_3d[1], tri_coords_3d[2], p1_bary_coords);
                  inserted = points1_arr->InsertNextTuple(px);
                  processed[a1] = inserted;
                }
                else
                {
                  inserted = processed[a1];
                }
              }
              else if (p2_on_q1_q2 != PointLinePosition::Outside)
              {
                if (processed.find(a2) == processed.end())
                {
                  std::copy(p2, p2 + 3, px);
                  interpZ(px, tri_coords_3d[0], tri_coords_3d[1], tri_coords_3d[2], p2_bary_coords);
                  inserted = points1_arr->InsertNextTuple(px);
                  processed[a2] = inserted;
                }
                else
                {
                  inserted = processed[a2];
                }
              }
              if (inserted >= 0)
              {
                hits.insert(inserted);
              }
            }
          }
        }

        if (hits.size() < 2)
        {
          continue;
        }
        for (auto it = hits.begin(); it != std::prev(hits.end()); ++it)
        {
          const auto& a = it;
          const auto& b = std::next(a);
          constraints.emplace_back(*a, *b);
        }
      }
    }
  };

  /**
   * @brief Robust point-in-polygon test based on
   * [this](https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html#The%20C%20Code)
   * @tparam ArrT         T
   * @param cutter_pt_ids   Indices into points
   * @param points        coordinates
   * @param test_y        test y-coordinate
   * @param test_x        test x-coordinate
   * @return
   */
  bool isInside(double pt[3], int npts, double* pts)
  {
    const double& test_x = pt[0];
    const double& test_y = pt[1];
    bool inside{ false };
    for (int i = 0, j = npts - 1; i < npts; j = i++)
    {
      const double& ix = pts[i * 3];
      const double& iy = pts[i * 3 + 1];
      const double& jx = pts[j * 3];
      const double& jy = pts[j * 3 + 1];

      if (((iy > test_y) != (jy > test_y)) && (test_x < (jx - ix) * (test_y - iy) / (jy - iy) + ix))
        inside = !inside;
    }
    return inside;
  }

  struct PopTrisImpl
  {
    vtkNew<vtkIdList> popped_cell_pts;
    PopTrisImpl() { this->popped_cell_pts->Allocate(3); }

    /**
     * @brief Populate output structures and arrays with children.
     *
     * @tparam PointsArrT    triangles' points' data array type
     * @tparam InOutsArrT     integral data array type
     * @param points_arr     triangles' points' data array
     * @param inside_out      inside out nature of all cutter polygons
     * @param cutters           [<bbox, ptIds>] for each cutter polygon
     * @param parent          the og triangle that we're processing rn
     * @param in_pd           input point data
     * @param out_pd          output point data
     * @param out_polys        output tris
     * @param out_lines       output lines
     * @param in_cd           input cell data
     * @param out_polys_cd     output cell data (for triangles)
     * @param out_lines_cd    output cell data (for lines)
     * @param locator         used to insert unique points in output
     * @param fallthrough     avoid in/out tests
     * @param cell_type       parent cell type
     *
     */
    template <typename PointsArrT>
    void operator()(PointsArrT* points_arr, const bool& inside_out,
      const CuttersInfoType& cutters_cache, Parent* parent, vtkPointData* in_pd,
      vtkPointData* out_pd, vtkCellArray* out_verts, vtkCellArray* out_lines,
      vtkCellArray* out_polys, vtkCellArray* out_strips, vtkCellData* in_cd,
      vtkCellData* out_verts_cd, vtkCellData* out_lines_cd, vtkCellData* out_polys_cd,
      vtkCellData* out_strips_cd, vtkIncrementalPointLocator* locator, const bool& fallthrough,
      const int& cell_type)
    {
      auto points = vtk::DataArrayTupleRange<3>(points_arr);

      vtkSmartPointer<vtkIdList> root_pt_ids = parent->PointIds;
      const vtkIdType& root_id = parent->cellId;

      for (auto& child : parent->children)
      {
        if (!fallthrough)
        {
          // initial condition clearly depends on inside_out. & vs |
          bool rejected = inside_out ? true : false;
          double x[3] = { child.cx, child.cy, 0.0 };
          auto& tri_bbox = child.bbox;

          for (const auto& cache : cutters_cache)
          {
            double* pts = cache.points->GetPointer(0);
            double bds[6] = {};
            int npts = cache.points->GetNumberOfTuples();
            cache.bbox.GetBounds(bds);

            if (inside_out)
            {
              // the hard-way; default
#ifdef USE_VTK_POINT_IN_POLYGON
              double normal[3] = { cache.normal[0], cache.normal[1], cache.normal[2] };
              bool is_inside = (vtkPolygon::PointInPolygon(x, npts, pts, bds, normal) == 1);
#else           
              bool is_inside = isInside(x, npts, pts);
#endif
              rejected &= !is_inside;
            }
            else
            {
              // the easy-way; effective when inside out is unset. (non-default)
              if (!(cache.bbox.Contains(tri_bbox) || cache.bbox.Intersects(tri_bbox)))
              {
                rejected |= false;
              }
              else
              {
                // the hard-way;
#ifdef USE_VTK_POINT_IN_POLYGON
                double normal[3] = { cache.normal[0], cache.normal[1], cache.normal[2] };
                bool is_inside = (vtkPolygon::PointInPolygon(x, npts, pts, bds, normal) == 1);
#else
                bool is_inside = isInside(x, npts, pts);
#endif
                rejected |= is_inside;
              }
            }
          }
          if (rejected)
          {
            continue;
          }
        }

        vtkCellArray* out_cells = nullptr;
        vtkCellData* out_cd = nullptr;
        // insert triangles
        switch (cell_type)
        {
          case VTK_LINE:
            out_cells = out_lines;
            out_cd = out_lines_cd;
            break;
          case VTK_POLYGON:
            out_cells = out_polys;
            out_cd = out_polys_cd;
            break;
          case VTK_POLY_LINE:
            out_cells = out_lines;
            out_cd = out_lines_cd;
            break;
          case VTK_POLY_VERTEX:
            out_cells = out_verts;
            out_cd = out_verts_cd;
            break;
          case VTK_QUAD:
            out_cells = out_polys;
            out_cd = out_polys_cd;
            break;
          case VTK_TRIANGLE_STRIP:
            out_cells = out_strips;
            out_cd = out_strips_cd;
            break;
          case VTK_TRIANGLE:
            out_cells = out_polys;
            out_cd = out_polys_cd;
            break;
          case VTK_VERTEX:
            out_cells = out_verts;
            out_cd = out_verts_cd;
            break;
          default:
            break;
        }
        if (out_cells && out_cd)
        {
          double dist2 = 0.;
          double* closest = nullptr;
          std::vector<double> weights(root_pt_ids->GetNumberOfIds(), 0.),
            pCoords(root_pt_ids->GetNumberOfIds(), 0.);
          int sub_id = 0;
          this->popped_cell_pts->Reset();

          for (const auto& pt : *(child.point_ids))
          {
            vtkIdType new_pt_id(-1);
            const double p[3] = { points[pt][0], points[pt][1], points[pt][2] };
            if (locator->InsertUniquePoint(p, new_pt_id))
            {
              parent->EvaluatePosition(p, closest, sub_id, pCoords.data(), dist2, weights.data());
              out_pd->InterpolatePoint(in_pd, new_pt_id, root_pt_ids, weights.data());
            }
            this->popped_cell_pts->InsertNextId(new_pt_id);
          }
          const auto& new_cell_id = out_cells->InsertNextCell(this->popped_cell_pts);
          out_cd->CopyData(in_cd, root_id, new_cell_id);
        }
      }
    }
  };

  struct ApplyConstraintImpl
  {
    /**
     * @brief Construct a new functor object.
     * @param mesh  input mesh
     * @param vl    a vertex id into mesh
     * @param vm    another vertex id into mesh
     */
    ApplyConstraintImpl(
      vtkSmartPointer<vtkPolyData> mesh_, const vtkIdType& vl_, const vtkIdType& vm_)
      : mesh(mesh_)
      , vl(vl_)
      , vm(vm_)
    {
    }

    vtkSmartPointer<vtkPolyData> mesh;
    vtkIdType vl, vm;

    /**
     * @brief Try to enforce edge vl_--vm_ in mesh_
     * vtkDelaunay2D's constraint section runs into infinite recursion for certain corner cases.
     * This implementation is no angel either. But, it does the job for points confined to a single
     * triangle.
     *
     * @tparam PointsArrT   triangles' points' data array type
     * @param points_arr    triangles' points' data array
     * @param tol           numeric tolerance for intersection math
     */
    template <typename PointsArrT, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT>>
    void operator()(PointsArrT* points_arr, const double& tol)
    {
      auto points = vtk::DataArrayTupleRange<3>(points_arr);

      double pl[3] = { points[vl][0], points[vl][1], 0.0 };
      double pm[3] = { points[vm][0], points[vm][1], 0.0 };

      vtkIdType* neighbors_vl = nullptr;
      vtkIdType num_neis_vl = 0;
      mesh->BuildLinks();
      mesh->GetPointCells(vl, num_neis_vl, neighbors_vl); // tris around vm;

      vtkIdType t0(-1), t1(-1), v1(-1), v2(-1), vm_(-1);
      // find a triangle st its edge opp to 'vl' intersects vl--vm
      bool abort_scan(false);
      for (vtkIdType t = 0; t < num_neis_vl && !abort_scan; ++t)
      {
        t0 = neighbors_vl[t];
        const vtkIdType* point_ids = nullptr;
        vtkIdType num_point_ids{ 0 };
        mesh->GetCellPoints(t0, num_point_ids, point_ids);
        // find an opposite edge vi--vj that intersects vl--vm
        for (const auto& edge : TRIEDGES)
        {
          const vtkIdType& id0 = edge.first;
          const vtkIdType& id1 = edge.second;
          const vtkIdType& vi = point_ids[id0];
          const vtkIdType& vj = point_ids[id1];
          if ((vi == vl && vj != vm) || (vi == vm && vj != vl) || (vj == vl && vi != vm) ||
            (vj == vm && vi != vl))
          {
            continue;
          }
          double pi[3] = { points[vi][0], points[vi][1], 0.0 };
          double pj[3] = { points[vj][0], points[vj][1], 0.0 };

          double px[3] = {};
          const auto intersect_type = robustIntersect(pi, pj, pl, pm, px, tol);
          if (intersect_type != IntersectType::Intersect)
          {
            continue;
          }

          v1 = vi;
          v2 = vj;
          vtkNew<vtkIdList> neighbors_v1v2;
          mesh->GetCellEdgeNeighbors(t0, v1, v2, neighbors_v1v2);
          if (!neighbors_v1v2->GetNumberOfIds()) // shouldn't happen, but eh just to be safe
          {
            continue;
          }
          t1 = neighbors_v1v2->GetId(0);
          abort_scan = true;
          break;
        }
      }

      if (abort_scan)
      {
        const vtkIdType* point_ids = nullptr;
        vtkIdType num_point_ids{ 0 };
        mesh->GetCellPoints(t1, num_point_ids, point_ids);
        for (const auto& edge : TRIEDGES)
        {
          const vtkIdType& id0 = edge.first;
          const vtkIdType& id1 = edge.second;
          if ((point_ids[id0] == v1 && point_ids[id1] == v2) ||
            (point_ids[id1] == v1 && point_ids[id0] == v2))
          {
            vm_ = point_ids[TRIOPPVERTS[id0]];
            vtkIdType* neisVm_;
            vtkIdType numNeisVm_{ 0 };
            mesh->GetPointCells(vm_, numNeisVm_, neisVm_);
            mesh->RemoveReferenceToCell(v1, t1);
            mesh->RemoveReferenceToCell(v2, t0);
            mesh->ResizeCellList(vl, 1);
            mesh->AddReferenceToCell(vl, t1);
            mesh->ResizeCellList(vm_, 1);
            mesh->AddReferenceToCell(vm_, t0);

            const vtkIdType t0New[3] = { v1, vl, vm_ };
            mesh->ReplaceCell(t0, 3, t0New);

            const vtkIdType t1New[3] = { vl, v2, vm_ };
            mesh->ReplaceCell(t1, 3, t1New);
            break;
          }
        }
      }
    }
  };

  struct SurfCutHelper
  {
    /**
     * @brief Construct a surface cutter helper.
     *
     * @param in_pd input point data
     * @param out_pd output point data
     * @param out_polys output tris
     * @param out_lines output lines
     * @param in_cd input cell data
     * @param out_tri_cd output cell data (triangles)
     * @param out_line_cd output cell data (line segments)
     * @param locator to insert unique points.
     */
    SurfCutHelper(const CuttersInfoType& cutters_cache_, const bool& inside_out_,
      vtkPointData* in_pd_, vtkPointData* out_pd_, vtkCellArray* out_verts_,
      vtkCellArray* out_lines_, vtkCellArray* out_polys_, vtkCellArray* out_strips_,
      vtkCellData* in_cd_, vtkCellData* out_verts_cd_, vtkCellData* out_lines_cd_,
      vtkCellData* out_polys_cd_, vtkCellData* out_strips_cd_, vtkIncrementalPointLocator* locator_,
      vtkPoints* in_points_, vtkPoints* in_cutter_points_)
      : cutters_cache(cutters_cache_)
      , inside_out(inside_out_)
      , out_verts(out_verts_)
      , out_lines(out_lines_)
      , out_polys(out_polys_)
      , out_strips(out_strips_)
      , in_cd(in_cd_)
      , out_verts_cd(out_verts_cd_)
      , out_lines_cd(out_lines_cd_)
      , out_polys_cd(out_polys_cd_)
      , out_strips_cd(out_strips_cd_)
      , locator(locator_)
      , in_pd(in_pd_)
      , out_pd(out_pd_)
      , in_points(in_points_)
      , in_cutter_points(in_cutter_points_)
    {
      del2d->SetProjectionPlaneMode(VTK_DELAUNAY_XY_PLANE);
      del2d->SetInputData(input);
      del2d->SetOffset(100.0); // bump this if Delaunay output is concave
    }

    /**
     * @brief Pop triangles into output data structures. (that were birthed after
     * `push()` followed by `update()`)
     *
     */
    void pop(const bool& fallthrough, const int& cell_type)
    {
      vtkDataArray* points_arr = parent->Points->GetData();
      if (!dispatchR::Execute(points_arr, this->pop_tris_worker, inside_out, cutters_cache, parent,
            in_pd, out_pd, out_verts, out_lines, out_polys, out_strips, in_cd, out_verts_cd,
            out_lines_cd, out_polys_cd, out_strips_cd, locator, fallthrough, cell_type))
      {
        this->pop_tris_worker(points_arr, inside_out, cutters_cache, parent, in_pd, out_pd,
          out_verts, out_lines, out_polys, out_strips, in_cd, out_verts_cd, out_lines_cd,
          out_polys_cd, out_strips_cd, locator, fallthrough, cell_type);
      }
    }

    /**
     * @brief Push a 'root' triangle for processing.
     *
     * @param point_ids     Point ids
     * @param root_idx      cell id from input
     */
    void push(vtkIdList* point_ids, const vtkIdType& root_idx, const int& cell_type)
    {
      switch (cell_type)
      {
        case VTK_LINE:
          parent->SetCellTypeToLine();
          break;
        case VTK_POLYGON:
          parent->SetCellTypeToPolygon();
          break;
        case VTK_POLY_LINE:
          parent->SetCellTypeToPolyLine();
          break;
        case VTK_POLY_VERTEX:
          parent->SetCellTypeToPolyVertex();
          break;
        case VTK_QUAD:
          parent->SetCellTypeToQuad();
          break;
        case VTK_TRIANGLE_STRIP:
          parent->SetCellTypeToTriangleStrip();
          break;
        case VTK_TRIANGLE:
          parent->SetCellTypeToTriangle();
          break;
        case VTK_VERTEX:
          parent->SetCellTypeToVertex();
          break;
        default:
          break;
      }

      parent->cellId = root_idx;
      parent->children.clear();
      parent->Points->SetNumberOfPoints(point_ids->GetNumberOfIds());
      parent->PointIds->DeepCopy(point_ids);
      polys->InsertNextCell(point_ids->GetNumberOfIds());

      for (vtkIdType i = 0; i < point_ids->GetNumberOfIds(); ++i)
      {
        polys->InsertCellPoint(i);
      }

      vtkDataArray* points1_arr = parent->Points->GetData();
      vtkDataArray* points2_arr = in_points->GetData();
      GetRootImpl worker;
      if (!dispatchRR::Execute(points1_arr, points2_arr, worker, parent->PointIds))
      {
        worker(points1_arr, points2_arr, parent->PointIds);
      }
    }

    /**
     * @brief Reset helper's internal state in anticipation for a new input
     *
     */
    void reset()
    {
      constraints.clear();
      parent->Reset();
      polys->Reset();
    }

    /**
     * @brief Triangulate vertices, lines from current internal state.
     *
     * @param tol numeric tolerance for delaunay math. (ignored)
     */
    void triangulate(const double& tol)
    {
      if (parent->Points->GetNumberOfPoints() > 3)
      {
        input->SetPoints(parent->Points);
        // del2d->SetTolerance(tol);
        del2d->Update();

        if (del2d->GetOutput())
        {
          output->ShallowCopy(del2d->GetOutput());

          if (output->GetNumberOfPolys())
          {
            output->BuildLinks();
            for (const auto& edge : constraints)
            {
              if (edge.first == edge.second)
              {
                continue;
              }

              unsigned short num_iters = 0;
              while (!output->IsEdge(edge.first, edge.second))
              {
                auto constraintWorker = ApplyConstraintImpl(output, edge.first, edge.second);
                vtkDataArray* points = output->GetPoints()->GetData();
                if (!dispatchR::Execute(points, constraintWorker, tol))
                {
                  constraintWorker(points, tol);
                }
                ++num_iters;
                if (num_iters > TSC_MAX_EDGE_FLIPS)
                {
                  //__debugbreak(); // happens only for extremely degenerate inputs.
                  break;
                }
              }
            }
            polys->ShallowCopy(output->GetPolys());
          }
        }
      }
    }

    /**
     * @brief Do triangle-line intersect tests in a 2D plane.
     *
     * @param lines a bunch of line segments
     * @param tol numeric tolerance for intersection math
     */
    void triIntersect(const SegmentsType& lines, const double& tol)
    {
      TriIntersect2dImpl worker;
      vtkDataArray* points1_arr = parent->Points->GetData();
      vtkDataArray* points2_arr = in_cutter_points->GetData();
      if (!dispatchRR::Execute(points1_arr, points2_arr, worker, lines, constraints, tol))
      {
        worker(points1_arr, points2_arr, lines, constraints, tol);
      }
    }

    /**
     * @brief If child triangles exist, bring them into output data structure.
     *        else, copy parent triangle.
     *
     */
    void update() { parent->UpdateChildren(polys); }

    vtkNew<Parent> parent;
    PopTrisImpl pop_tris_worker;
    SegmentsType constraints;
    const CuttersInfoType& cutters_cache;
    bool inside_out;
    vtkNew<vtkCellArray> polys;
    vtkNew<vtkDelaunay2D> del2d;
    vtkNew<vtkPolyData> input, output;
    vtkCellArray *out_verts = nullptr, *out_lines = nullptr, *out_polys = nullptr,
                 *out_strips = nullptr;
    vtkCellData *in_cd = nullptr, *out_verts_cd = nullptr, *out_lines_cd = nullptr,
                *out_polys_cd = nullptr, *out_strips_cd = nullptr;
    vtkIncrementalPointLocator* locator = nullptr;
    vtkPointData *in_pd = nullptr, *out_pd = nullptr;
    vtkPoints *in_points = nullptr, *in_cutter_points = nullptr;
  };

}

int tscTriSurfaceCutter::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkSmartPointer<vtkPolyData> input = vtkPolyData::GetData(inputVector[0]);
  vtkSmartPointer<vtkPolyData> cutters = vtkPolyData::GetData(inputVector[1]);

  vtkSmartPointer<vtkPolyData> output = vtkPolyData::GetData(outputVector);

  if (!this->Imprint && !this->Remove) // do nothing
  {
    output->DeepCopy(input);
    return 1;
  }

  ///> Allocate output data structures.
  vtkSmartPointer<vtkPoints> in_points = input->GetPoints();
  vtkSmartPointer<vtkPoints> in_cutter_points = cutters->GetPoints();

  vtkNew<vtkPoints> out_pts;
  vtkNew<vtkCellArray> out_verts, out_lines, out_polys, out_strips;
  vtkNew<vtkCellData> out_verts_cd, out_lines_cd, out_polys_cd, out_strips_cd;
  vtkSmartPointer<vtkCellIterator> iter, cutters_iter;
  iter.TakeReference(input->NewCellIterator());
  cutters_iter.TakeReference(cutters->NewCellIterator());

  const vtkIdType& num_points = input->GetNumberOfPoints();
  const vtkIdType& num_cutters_points = cutters->GetNumberOfPoints();
  const vtkIdType& num_cells = input->GetNumberOfCells();
  const vtkIdType& num_verts = input->GetNumberOfVerts();
  const vtkIdType& num_strips = input->GetNumberOfStrips();
  const vtkIdType& num_polys = input->GetNumberOfPolys();
  const vtkIdType& num_cutter_polys = cutters->GetNumberOfPolys();
  const vtkIdType& num_cutter_lines = cutters->GetNumberOfLines();
  const vtkIdType& num_cutter_cells = cutters->GetNumberOfCells();

  if (!(num_points && num_polys))
  {
    vtkWarningMacro(<< "Input mesh is empty.");
    return 1;
  }
  if (!(num_cutters_points && (num_cutter_polys || num_cutter_lines)))
  {
    vtkWarningMacro(<< "Input cutters are empty.");
    return 1;
  }

  vtkSmartPointer<vtkPointData> in_pd = input->GetPointData();
  vtkSmartPointer<vtkPointData> out_pd = output->GetPointData();
  vtkSmartPointer<vtkCellData> in_cd = input->GetCellData();

  out_verts->Allocate(num_verts);
  vtkDebugMacro(<< "Alloc'd " << num_verts << " verts");

  out_lines->Allocate(num_cutters_points - 1 + (num_polys >> 4));
  vtkDebugMacro(<< "Alloc'd " << num_cutters_points - 1 + (num_polys >> 4) << " lines");

  out_polys->Allocate(num_polys + num_cutters_points - 1);
  vtkDebugMacro(<< "Alloc'd " << num_polys + num_cutters_points - 1 << " polys");

  out_strips->Allocate(num_strips);
  vtkDebugMacro(<< "Alloc'd " << num_verts << " strips");

  out_pts->SetDataType(in_points->GetDataType());
  out_pts->Allocate(num_points + num_cutters_points);
  vtkDebugMacro(<< "Alloc'd " << num_points + num_cutters_points << " points");

  out_pd->InterpolateAllocate(in_pd, num_points + num_cutters_points * 3);
  vtkDebugMacro(<< "Alloc'd " << num_points + num_cutters_points * 3 << " point data tuples");

  out_verts_cd->CopyAllocate(in_cd);
  vtkDebugMacro(<< "Alloc'd CellData for verts");

  out_lines_cd->CopyAllocate(in_cd);
  vtkDebugMacro(<< "Alloc'd CellData for lines");

  out_polys_cd->CopyAllocate(in_cd);
  vtkDebugMacro(<< "Alloc'd CellData for polys");

  out_strips_cd->CopyAllocate(in_cd);
  vtkDebugMacro(<< "Alloc'd CellData for strips");
  ///> Finished allocation

  ///> Find candidate cutter-mesh crossings
  CuttersInfoType cutters_cache(num_cutter_cells);
  std::vector<vtkBoundingBox> cells_cache;
  double cutters_bds[6] = {}, in_bds[6] = {};
  std::unordered_map<vtkIdType, SegmentsType> possible_crossings;
  auto pt_ids = vtkSmartPointer<vtkIdList>::New();
  auto points = vtkSmartPointer<vtkPoints>::New();

  vtkDebugMacro(<< "Determine possible cell-edge crossings.");
  in_cutter_points->GetBounds(cutters_bds);
  in_points->GetBounds(in_bds);
  possible_crossings.reserve(num_cells);

  // Cache info. This helps when there are lots of cutting edges
  for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextCell())
  {
    const vtkIdType& cell_id = iter->GetCellId();
    double bds[6] = {};
    input->GetCellBounds(cell_id, bds);
    cells_cache.emplace_back(bds);
  }

  for (cutters_iter->InitTraversal(); !cutters_iter->IsDoneWithTraversal();
       cutters_iter->GoToNextCell())
  {
    const int& cutter_type = cutters_iter->GetCellType();
    switch (cutter_type)
    {
      case VTK_QUAD:
      case VTK_TRIANGLE:
      case VTK_POLYGON:
      case VTK_LINE:
      case VTK_POLY_LINE:
        break;
      default: // Ignore all other cell types from cutter input.
        continue;
    }

    const vtkIdType& cutter_id = cutters_iter->GetCellId();
    double bds[6] = {};
    cutters->GetCellBounds(cutter_id, bds);
    bds[4] = bds[5] = 0.0;

    pt_ids = cutters_iter->GetPointIds();
    points = cutters_iter->GetPoints();
    const vtkIdType& npts = pt_ids->GetNumberOfIds();

    cutters_cache[cutter_id].bbox = vtkBoundingBox(bds);
    cutters_cache[cutter_id].points->SetNumberOfComponents(3);
    cutters_cache[cutter_id].points->SetNumberOfTuples(npts);
    for (unsigned short dim = 0; dim < 3; ++dim)
    {
      cutters_cache[cutter_id].points->CopyComponent(dim, points->GetData(), dim);
    }
    cutters_cache[cutter_id].cell_type = cutter_type;
    vtkPolygon::ComputeNormal(points, cutters_cache[cutter_id].normal);

    for (vtkIdType i = 0; i < npts; ++i)
    {
      const vtkIdType& this_i = i;
      const vtkIdType next_i = (i + 1) % npts;
      const vtkIdType& e0 = pt_ids->GetId(this_i);
      const vtkIdType& e1 = pt_ids->GetId(next_i);
      double p0[3], p1[3] = {};

      in_cutter_points->GetPoint(e0, p0);
      in_cutter_points->GetPoint(e1, p1);

      double edge_bounds[6] = {};
      vtkBoundingBox edge_bbox;
      edge_bbox.AddPoint(p0);
      edge_bbox.AddPoint(p1);
      edge_bbox.GetBounds(edge_bounds);

      // extend up to input's z{min, max}
      std::copy(in_bds + 4, in_bds + 6, edge_bounds + 4);
      edge_bbox.SetBounds(edge_bounds);

      // Do a box hit-test for each cell in input mesh.
      for (vtkIdType cell_id = 0; cell_id < num_cells; ++cell_id)
      {
        const int& cell_type = input->GetCellType(cell_id);
        if (cell_type != VTK_TRIANGLE)
        {
          continue;
        }

        if (cells_cache[cell_id].Intersects(edge_bbox))
        {
          possible_crossings[cell_id].emplace_back(e0, e1);
        }
      }
    }
  }
  ///> Finished pre-processing input.

  ///> Imprint cutters in candidate cells and remove if needed.
  // strategy for triangles: 1. push triangle into helper,
  //                         2. try to intersect with cutters,
  //                         3. pop child triangles into output poly-data
  //                         4. reset helper's data structures to prepare for next triangle.
  this->CreateDefaultLocators();
  SurfCutHelper helper(cutters_cache, this->InsideOut, in_pd, out_pd, out_verts, out_lines,
    out_polys, out_strips, in_cd, out_verts_cd, out_lines_cd, out_polys_cd, out_strips_cd,
    this->PointLocator, in_points, in_cutter_points);

  // roughly, a quarter no. of cells
  int report_every =
    (num_cells >= 4) ? num_cells >> 2 : (num_cells >= 2 ? num_cells >> 1 : num_cells);
  vtkDebugMacro(<< "Report progress every " << report_every << " cells");

  this->PointLocator->SetTolerance(this->Tolerance);
  this->PointLocator->InitPointInsertion(out_pts, in_bds);
  vtkDebugMacro(<< "Init'd point insertion");

  iter.TakeReference(input->NewCellIterator());
  for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextCell())
  {
    const vtkIdType& cell_id = iter->GetCellId();
    const int& cell_type = input->GetCellType(cell_id);
    if (!(cell_id % report_every))
    {
      vtkDebugMacro(<< "Processing " << cell_id);
    }

    // provide
    pt_ids = iter->GetPointIds();
    helper.push(pt_ids, cell_id, cell_type);
    if (cell_type == VTK_TRIANGLE)
    {
      // embed
      const auto& trials = possible_crossings.find(cell_id);
      if ((trials != possible_crossings.end()) && this->Imprint)
      {
        const SegmentsType& edges = trials->second;
        helper.triIntersect(edges, this->Tolerance);
        helper.triangulate(this->Tolerance);
      }
    }
    helper.update();

    // accept/reject
    switch (cell_type)
    {
      case VTK_TRIANGLE:
        helper.pop(!this->Remove, cell_type);
        break;
      default:
        helper.pop(true, cell_type); // fall through.
        break;
    }

    helper.reset();

    if (!(cell_id % report_every))
    {
      this->UpdateProgress(static_cast<double>(cell_id) / num_cells);
    }
  }
  ///> Finished embed and removal.

  ///> Finalize output ds.
  out_pts->Squeeze();
  output->SetPoints(out_pts);
  out_verts->Squeeze();
  output->SetVerts(out_verts);
  out_lines->Squeeze();
  output->SetLines(out_lines);
  out_polys->Squeeze();
  output->SetPolys(out_polys);
  out_strips->Squeeze();
  output->SetStrips(out_strips);

  out_pd->Squeeze();
  vtkDebugMacro("Finalized output point attributes. NumPoints: " << out_pts->GetNumberOfPoints());

  const vtkIdType& num_out_lines = output->GetNumberOfLines();
  const vtkIdType& num_out_polys = output->GetNumberOfPolys();
  const vtkIdType& num_out_cells = num_verts + num_out_lines + num_out_polys + num_strips;

  const vtkIdType nvl = num_verts + num_out_lines;
  const vtkIdType nvlp = nvl + num_out_polys;

  vtkSmartPointer<vtkCellData> out_cd = output->GetCellData();
  out_lines_cd->Squeeze();
  out_polys_cd->Squeeze();
  out_cd->CopyAllocate(out_polys_cd, num_out_cells);
  out_cd->CopyData(out_verts_cd, 0, num_verts, 0);
  out_cd->CopyData(out_lines_cd, num_verts, num_out_lines, 0);
  out_cd->CopyData(out_polys_cd, nvl, num_out_polys, 0);
  out_cd->CopyData(out_strips_cd, nvlp, num_strips, 0);
  out_cd->Squeeze();

  vtkDebugMacro("Finalized output cell attributes.");
  vtkDebugMacro("Verts: " << num_verts << " Lines: " << num_out_lines << " Polys: " << num_out_polys
                          << " Strips: " << num_strips << "Cells: " << num_out_cells);
  ///> Finished finalizing output data structures.

  return 1;
}

void tscTriSurfaceCutter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "InsideOut: " << (this->InsideOut ? "True" : "False") << "\n";
  os << indent << "Tolerance: " << this->Tolerance << "\n";
  os << indent << "Tolerance2: " << this->Tolerance * this->Tolerance << "\n";
  os << indent << "Imprint: " << (this->Imprint ? "True" : "False") << "\n";
  os << indent << "Remove: " << (this->Remove ? "True" : "False") << "\n";

  if (this->PointLocator)
  {
    os << indent << "Point Locator: " << this->PointLocator << "\n";
    this->PointLocator->PrintSelf(os, indent);
  }
  else
  {
    os << indent << "Point Locator: (None)\n";
  }
}

void tscTriSurfaceCutter::CreateDefaultLocators()
{
  if (!this->PointLocator)
  {
    this->PointLocator = vtkSmartPointer<vtkMergePoints>::New();
    this->PointLocator->SetTolerance(this->Tolerance);
  }
}
