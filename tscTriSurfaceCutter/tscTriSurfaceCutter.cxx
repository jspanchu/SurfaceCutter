#include "tscTriSurfaceCutter.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iterator>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
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
#include <vtkGenericCell.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLine.h>
#include <vtkLogger.h>
#include <vtkMath.h>
#include <vtkMergePoints.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStaticCellLocator.h>
#include <vtkThreadedCompositeDataPipeline.h>
#include <vtkTriangle.h>
#include <vtkTypeInt32Array.h>
#include <vtkUnsignedCharArray.h>

using dispatchRR =
  vtkArrayDispatch::Dispatch2ByValueType<vtkArrayDispatch::Reals, vtkArrayDispatch::Reals>;
using dispatchR = vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals>;

using SegmentType = std::pair<vtkIdType, vtkIdType>;
using SegmentsType = std::vector<SegmentType>;
using LoopInfoType = std::pair<vtkBoundingBox, vtkNew<vtkIdList>>;
using LoopsInfoType = std::vector<LoopInfoType>;

vtkStandardNewMacro(tscTriSurfaceCutter);

tscTriSurfaceCutter::tscTriSurfaceCutter()
  : AccelerateCellLocator(true)
  , ColorAcquiredPts(true)
  , ColorLoopEdges(true)
  , Embed(true)
  , InsideOut(true) // default: remove portions outside loop polygons.
  , Remove(true)
  , Tolerance(1.0e-6)
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
  this->DebugOn();
  this->CreateDefaultLocators();
}

tscTriSurfaceCutter::~tscTriSurfaceCutter() {}

void tscTriSurfaceCutter::SetLoopsData(vtkPolyData* loops)
{
  this->SetInputData(1, loops);
}

void tscTriSurfaceCutter::SetLoopsConnection(vtkAlgorithmOutput* output)
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

namespace tsc_detail
{
  struct GetChildrenImpl
  {
    /**
     * @brief Implementation to fill a container with children triangles.
     * @tparam PointsArrT T
     * @param points_arr  coordinates that make up parent triangle and acquired points
     * @param polys       point ids that make up children triangles
     * @param children    collect children into this container
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
          pc[dim] /= double(num_point_ids);

        children.emplace_back(pc[0], pc[1], point_ids, num_point_ids, bounds);
      }
    }
  };

  // Determine a point's location in a triangle.
  // is it 1. on an edge of a triangle?
  //       2. inside/outside ?
  //       3. exactly coincident with a vertex?
  //       4. impossible to determine (triangle is degenerate)
  tsc_detail::PointInTriangle inTriangle(const double p[3], const double p0[3], const double p1[3],
    const double p2[3], double bary_coords[3], const double& tol, vtkIdType& v, vtkIdType& e)
  {
    v = -1;
    e = -1;
    if (!vtkTriangle::BarycentricCoords(p, p0, p1, p2, bary_coords))
      return tsc_detail::PointInTriangle::Degenerate;

    const double& s = bary_coords[0];
    const double& t = bary_coords[1];
    const double& st = bary_coords[2]; // 1 - s - t

    bool s_eq_0 = std::fpclassify(s) == FP_ZERO || std::fabs(s) <= tol;
    bool t_eq_0 = std::fpclassify(t) == FP_ZERO || std::fabs(t) <= tol;
    bool st_eq_0 = std::fpclassify(st) == FP_ZERO || std::fabs(st) <= tol;

    bool s_eq_1 = std::fabs(s - 1.0) <= tol;
    bool t_eq_1 = std::fabs(t - 1.0) <= tol;
    bool st_eq_1 = std::fabs(st - 1.0) <= tol;

    if (s_eq_0 && t_eq_0 && s_eq_1)
      v = 0;
    else if (s_eq_0 && st_eq_0 && t_eq_1)
      v = 1;
    else if (t_eq_0 && st_eq_0 && st_eq_1)
      v = 2;
    if (v >= 0)
    {
      return tsc_detail::PointInTriangle::OnVertex;
    }

    bool s_gt0_lt1 = !std::signbit(s) && std::isless(s, 1) && !s_eq_0;
    bool t_gt0_lt1 = !std::signbit(t) && std::isless(t, 1) && !t_eq_0;
    bool st_gt0_lt1 = !std::signbit(st) && std::isless(st, 1) && !st_eq_0;

    if (s_gt0_lt1 && t_gt0_lt1 && st_gt0_lt1)
      return tsc_detail::PointInTriangle::Inside;

    if (s_eq_0 && t_gt0_lt1 && st_gt0_lt1)
      e = tsc_detail::TRIOPPEDGE[0];
    else if (t_eq_0 && s_gt0_lt1 && st_gt0_lt1)
      e = tsc_detail::TRIOPPEDGE[1];
    else if (st_eq_0 && s_gt0_lt1 && t_gt0_lt1)
      e = tsc_detail::TRIOPPEDGE[2];
    if (e >= 0)
      return tsc_detail::PointInTriangle::OnEdge;

    if (s_eq_0 && t_eq_0 && st_eq_0)
      return tsc_detail::PointInTriangle::Degenerate;

    return tsc_detail::PointInTriangle::Outside;
  }

  // Reasonable z-value of a point inside(or on an edge/vertex) a triangle.
  inline void interpZ(double p[3], const double p0[3], const double p1[3], const double p2[3])
  {
    double bary_coords[3] = {};
    vtkTriangle::BarycentricCoords(p, p0, p1, p2, bary_coords);
    p[2] = bary_coords[0] * p0[2] + bary_coords[1] * p1[2] + bary_coords[2] * p2[2];
  }

  // Same as above, but when you already have bary_coords
  inline void interpZ(
    double p[3], const double p0[3], const double p1[3], const double p2[3], double bary_coords[3])
  {
    p[2] = bary_coords[0] * p0[2] + bary_coords[1] * p1[2] + bary_coords[2] * p2[2];
  }

  // Determine a point's location on a line segment.
  // is it 1. on one of the end vertices?
  //       2. inside/outside ?
  tsc_detail::PointOnLine onLine(
    const double p[3], const double p1[3], const double p2[3], const double& tol)
  {
    double t{ 0 };
    const double dist2 = vtkLine::DistanceToLine(p, p1, p2, t, nullptr);

    if (dist2 > tol * tol)
    {
      return tsc_detail::PointOnLine::Outside;
    }

    const double& s = 1 - t;

    bool s_eq_0 = std::fpclassify(s) == FP_ZERO || std::fabs(s) <= tol;
    bool t_eq_0 = std::fpclassify(t) == FP_ZERO || std::fabs(t) <= tol;

    bool s_eq_1 = std::fabs(s - 1.0) <= tol;
    bool t_eq_1 = std::fabs(t - 1.0) <= tol;

    if ((s_eq_0 && t_eq_1) || (s_eq_1 && t_eq_0))
    {
      return tsc_detail::PointOnLine::OnVertex;
    }

    bool s_gt0_lt1 = !std::signbit(s) && std::isless(s, 1) && !s_eq_0;
    bool t_gt0_lt1 = !std::signbit(t) && std::isless(t, 1) && !t_eq_0;

    if (s_gt0_lt1 && t_gt0_lt1)
      return tsc_detail::PointOnLine::Inside;

    return tsc_detail::PointOnLine::Outside;
  }

  // Quantify return status of vtkLine::Intersection(...)
  tsc_detail::IntersectType robustIntersect(
    double* p1, double* p2, double* q1, double* q2, double px[3], const double& tol)
  {
    double u(0.0), v(0.0);
    int vtk_intersect = 0;

    // Note: this section is based on consts def'd in vtkLine.cxx.
    // Keep in mind to update this when that changes.
    // Just so we're clear:
    // 0: vtkLine::NoIntersect
    // 2: vtkLine::Intersect
    // 3: vtkLine::OnLine (Consider this as tsc_detail::IntersectType::NoIntersection)
    if ((vtk_intersect = vtkLine::Intersection(p1, p2, q1, q2, u, v)) == 2)
    {
      // This test checks for vertex and edge intersections
      // For example
      //  Vertex intersection
      //    (u=0 v=0), (u=0 v=1), (u=1 v=0), (u=1 v=0)
      //  Edge intersection
      //    (u=0 v!=0 v!=1), (u=1 v!=0 v!=1)
      //    (u!=0 u!=1 v=0), (u!=0 u!=1 v=1)
      if ((tol < u) && (u < 1.0 - tol) && (tol < v) && (v < 1.0 - tol))
      {
        for (unsigned short dim = 0; dim < 2; ++dim)
        {
          px[dim] = p1[dim] + u * (p2[dim] - p1[dim]);
        }
        return tsc_detail::IntersectType::Intersect;
      }
      else
      {
        return tsc_detail::IntersectType::NoIntersection;
      }
    }
    else if (vtk_intersect == 3)
    {
      double p21[3], c[3], c2[3];
      vtkMath::Subtract(p2, p1, p21);
      double tol2 = /*tolerance * tolerance*/ 1.0e-8 * sqrt(vtkMath::Dot(p21, p21));
      if (vtkLine::DistanceBetweenLines(p1, p2, q1, q2, c, c2, u, v) > tol2)
      {
        return tsc_detail::IntersectType::NoIntersection;
      }
      vtkLine::DistanceToLine(p1, p2, q1, u, px);
      if (-tol <= u && u <= (1.0 + tol))
      {
        return tsc_detail::IntersectType::Junction;
      }
      vtkLine::DistanceToLine(p2, q1, q2, u, px);
      if (-tol <= u && u <= (1.0 + tol))
      {
        return tsc_detail::IntersectType::Junction;
      }
      vtkLine::DistanceToLine(q1, p1, p2, u, px);
      if (-tol <= u && u <= (1.0 + tol))
      {
        return tsc_detail::IntersectType::Junction;
      }
      vtkLine::DistanceToLine(q2, p1, p2, u, px);
      if (-tol <= u && u <= (1.0 + tol))
      {
        return tsc_detail::IntersectType::Junction;
      }
    }

    return tsc_detail::IntersectType::NoIntersection;
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
     * @brief   A brute-force implementation that intersects(in 2D plane @z=0.0) a triangle with a
     *          bunch of line segments.
     *          It uses hints such as a point's location w.r.t a triangle to avoid actually
     *          testing intersection of line segments.
     *          If this implementation could be modified to work for a concave polygon,
     *          we can cut polygons too!
     * @tparam PointsArrT1  triangle's points' data array T
     * @tparam PointsArrT2  loops' points' data array T
     * @param points1_arr   triangle's points' data array
     * @param points2_arr   loops' points' data array
     * @param lines         point ids of line segments
     * @param constraints   list of segments that make up constraints
     * @param is_acquired   contains points that have been acquired from lines' points
     * @param tol           numeric tolerance for intersection math.
     *
     */
    template <typename PointsArrT1, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT1>,
      typename PointsArrT2, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT1>>
    void operator()(PointsArrT1* points1_arr, PointsArrT2* points2_arr, const SegmentsType& lines,
      SegmentsType& constraints, std::unordered_set<vtkIdType>& is_acquired, const double& tol)
    {
      if (!lines.size())
        return;

      auto points_1 = vtk::DataArrayTupleRange<3>(points1_arr); // triangle's points
      auto points_2 = vtk::DataArrayTupleRange<3>(points2_arr); // loop's points

      // clang-format off
      double p2d[3][3] = 
      { { points_1[0][0], points_1[0][1], 0.0 },
        { points_1[1][0], points_1[1][1], 0.0 }, 
        { points_1[2][0], points_1[2][1], 0.0 } };  // use this for point in triangle, intersection tests.

      double p3d[3][3] = 
      { { points_1[0][0], points_1[0][1], points_1[0][2] },
        { points_1[1][0], points_1[1][1], points_1[1][2] },
        { points_1[2][0], points_1[2][1], points_1[2][2] } }; // use this for z-interpolation
      // clang-format on

      std::unordered_map<vtkIdType, vtkIdType> processed; // k: line pt, v: insertLoc

      // a1, a2: line point ids.
      // p1, p2: line point coordinates.
      // b1, b2: triangle edge point ids.
      // q1, q2: triangle edge coordinates.
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

        const auto a1Pos =
          tsc_detail::inTriangle(p1, p2d[0], p2d[1], p2d[2], p1_bary_coords, tol, a1v, a1e);
        const auto a2Pos =
          tsc_detail::inTriangle(p2, p2d[0], p2d[1], p2d[2], p2_bary_coords, tol, a2v, a2e);

        // for convenience
        const bool a1_inside = a1Pos == tsc_detail::PointInTriangle::Inside;
        const bool a2_inside = a2Pos == tsc_detail::PointInTriangle::Inside;
        const bool a1_on_edge = a1Pos == tsc_detail::PointInTriangle::OnEdge;
        const bool a2_on_edge = a2Pos == tsc_detail::PointInTriangle::OnEdge;
        const bool a1_on_vert = a1Pos == tsc_detail::PointInTriangle::OnVertex;
        const bool a2_on_vert = a2Pos == tsc_detail::PointInTriangle::OnVertex;
        const bool a1_outside = a1Pos == tsc_detail::PointInTriangle::Outside;
        const bool a2_outside = a2Pos == tsc_detail::PointInTriangle::Outside;

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
            tsc_detail::interpZ(px, p3d[0], p3d[1], p3d[2], p1_bary_coords);
            inserted = points1_arr->InsertNextTuple(px);
            processed[a1] = inserted;
            is_acquired.insert(inserted);
          }
          hits.insert(processed[a1]);
        }

        if (a2_inside || a2_on_edge)
        {
          if (processed.find(a2) == processed.end())
          {
            std::copy(p2, p2 + 3, px);
            tsc_detail::interpZ(px, p3d[0], p3d[1], p3d[2], p2_bary_coords);
            inserted = points1_arr->InsertNextTuple(px);
            processed[a2] = inserted;
            is_acquired.insert(inserted);
          }
          hits.insert(processed[a2]);
        }

        if (a1_on_vert)
        {
          if (processed.find(a1) == processed.end())
          {
            inserted = a1v;
            processed[a1] = inserted;
            is_acquired.insert(inserted);
          }
          hits.insert(processed[a1]);
        }

        if (a2_on_vert)
        {
          if (processed.find(a2) == processed.end())
          {
            inserted = a2v;
            processed[a2] = inserted;
            is_acquired.insert(inserted);
          }
          hits.insert(processed[a2]);
        }

        // Call robustIntersect()
        if ((a1_outside && a2_outside) || (a1_outside && !a2_outside) ||
          (!a1_outside && a2_outside))
        {
          for (const auto& edge : tsc_detail::TRIEDGES)
          {
            const vtkIdType& b1 = edge.first;
            const vtkIdType& b2 = edge.second;

            double* q1 = p2d[b1];
            double* q2 = p2d[b2];

            const auto intersect_type = tsc_detail::robustIntersect(p1, p2, q1, q2, px, tol);
            if (intersect_type == tsc_detail::IntersectType::Intersect)
            {
              tsc_detail::interpZ(px, p3d[0], p3d[1], p3d[2]);
              inserted = points1_arr->InsertNextTuple(px);
              hits.insert(inserted);
            }
            else
            {
              // possibly on line or a junction.
              // junctions were handled before performing intersections (ax_on_edge or ax_on_vert)
              inserted = -1;
              const auto q1_on_p1_p2 = tsc_detail::onLine(q1, p1, p2, tol);
              const auto q2_on_p1_p2 = tsc_detail::onLine(q2, p1, p2, tol);
              const auto p1_on_q1_q2 = tsc_detail::onLine(p1, q1, q2, tol);
              const auto p2_on_q1_q2 = tsc_detail::onLine(p2, q1, q2, tol);
              if (q1_on_p1_p2 != tsc_detail::PointOnLine::Outside)
              {
                inserted = b1;
              }
              else if (q2_on_p1_p2 != tsc_detail::PointOnLine::Outside)
              {
                inserted = b2;
              }
              else if (p1_on_q1_q2 != tsc_detail::PointOnLine::Outside)
              {
                if (processed.find(a1) == processed.end())
                {
                  std::copy(p1, p1 + 3, px);
                  tsc_detail::interpZ(px, p3d[0], p3d[1], p3d[2], p1_bary_coords);
                  inserted = points1_arr->InsertNextTuple(px);
                  processed[a1] = inserted;
                  is_acquired.insert(inserted);
                }
                else
                {
                  inserted = processed[a1];
                }
              }
              else if (p2_on_q1_q2 != tsc_detail::PointOnLine::Outside)
              {
                if (processed.find(a2) == processed.end())
                {
                  std::copy(p2, p2 + 3, px);
                  tsc_detail::interpZ(px, p3d[0], p3d[1], p3d[2], p2_bary_coords);
                  inserted = points1_arr->InsertNextTuple(px);
                  processed[a2] = inserted;
                  is_acquired.insert(inserted);
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
   * @param loop_pt_ids   Indices into points
   * @param points        coordinates
   * @param test_y        test y-coordinate
   * @param test_x        test x-coordinate
   * @return
   */
  template <typename ArrT>
  bool isInside(vtkIdList* loop_pt_ids, vtk::detail::ConstTupleIterator<ArrT, 3> points,
    const double& test_x, const double& test_y)
  {
    bool inside{ false };
    const vtkIdType num_loop_pts(loop_pt_ids->GetNumberOfIds());
    for (vtkIdType i = 0, j = num_loop_pts - 1; i < num_loop_pts; j = i++)
    {
      const double& ix = points[loop_pt_ids->GetId(i)][0];
      const double& iy = points[loop_pt_ids->GetId(i)][1];
      const double& jx = points[loop_pt_ids->GetId(j)][0];
      const double& jy = points[loop_pt_ids->GetId(j)][1];

      if (((iy > test_y) != (jy > test_y)) && (test_x < (jx - ix) * (test_y - iy) / (jy - iy) + ix))
        inside = !inside;
    }
    return inside;
  }

  struct PopTrisImpl
  {
    /**
     * @brief Populate output structures and arrays with children.
     *
     * @tparam PointsArrT1    triangles' points' data array type
     * @tparam PointsArrT2    loops' points' data array type
     * @tparam InOutsArrT     integral data array type
     * @param points1_arr     triangles' points' data array
     * @param points2_arr     loops' points' data array
     * @param inside_out      inside out nature of all loop polygons
     * @param loops           [<bbox, ptIds>] for each loop polygon
     * @param is_acquired     indicates points that were acquired
     * @param constraints     a list of line segments that were constrained
     * @param parent          the og triangle that we're processing rn
     * @param in_pd           input point data
     * @param out_pd          output point data
     * @param out_polys        output tris
     * @param out_lines       output lines
     * @param in_cd           input cell data
     * @param out_polys_cd     output cell data (for triangles)
     * @param out_lines_cd    output cell data (for lines)
     * @param locator         used to insert unique points in output
     * @param acquisition     colors acquired points in output
     * @param fallthrough     avoid in/out tests
     *
     */
    template <typename PointsArrT1, typename PointsArrT2>
    void operator()(PointsArrT1* points1_arr, PointsArrT2* points2_arr, const bool& inside_out,
      const LoopsInfoType& loops, const std::unordered_set<vtkIdType>& is_acquired,
      const SegmentsType& constraints, Parent* parent, vtkPointData* in_pd, vtkPointData* out_pd,
      vtkCellArray* out_verts, vtkCellArray* out_lines, vtkCellArray* out_polys,
      vtkCellArray* out_strips, vtkCellData* in_cd, vtkCellData* out_verts_cd,
      vtkCellData* out_lines_cd, vtkCellData* out_polys_cd, vtkCellData* out_strips_cd,
      vtkIncrementalPointLocator* locator, vtkUnsignedCharArray* acquisition,
      const bool& fallthrough, const int& cell_type)
    {
      auto points_1 = vtk::DataArrayTupleRange<3>(points1_arr);
      auto points_2 = vtk::DataArrayTupleRange<3>(points2_arr);

      vtkSmartPointer<vtkIdList> root_pt_ids = parent->PointIds;
      const vtkIdType& num_points = parent->Points->GetNumberOfPoints();
      const vtkIdType& root_id = parent->cellId;
      std::vector<vtkIdType> old_to_new_pt_ids(num_points, -1);

      for (auto& child : parent->children)
      {
        if (!fallthrough)
        {
          // initial condition clearly depends on inside_out. & vs |
          bool rejected = inside_out ? true : false;

          const double& test_x = child.cx;
          const double& test_y = child.cy;
          auto& tri_bbox = child.bbox;

          for (auto& loop : loops)
          {
            if (inside_out)
            {
              // the hard-way; default
              rejected &= !isInside(loop.second.Get(), points_2.cbegin(), test_x, test_y);
            }
            else
            {
              // the easy-way; effective when inside out is unset. (non-default)
              if (!(loop.first.Contains(tri_bbox) || loop.first.Intersects(tri_bbox)))
              {
                rejected |= false;
              }
              else
              {
                // the hard-way;
                rejected |= isInside(loop.second.Get(), points_2.cbegin(), test_x, test_y);
              }
            }
          }
          if (rejected)
          {
            continue;
          }
        }

        double dist2(0.);
        double* closest = nullptr;
        double parametric_coords[3] = {};
        double weights[3] = {};
        int sub_id{ 0 };
        vtkIdType new_pt_id(-1);

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
          out_cells->InsertNextCell(child.point_ids->GetNumberOfIds());
          for (const auto& pt : *(child.point_ids))
          {
            const double p[3] = { points_1[pt][0], points_1[pt][1], points_1[pt][2] };
            if (locator->InsertUniquePoint(p, new_pt_id))
            {
              parent->EvaluatePosition(p, closest, sub_id, parametric_coords, dist2, weights);
              out_pd->InterpolatePoint(in_pd, new_pt_id, root_pt_ids, weights);

              if (is_acquired.find(pt) != is_acquired.end())
                acquisition->InsertNextValue('\001');
              else
                acquisition->InsertNextValue('\000');
            }
            old_to_new_pt_ids[pt] = new_pt_id;
            out_cells->InsertCellPoint(new_pt_id);
          }
          out_cd->InsertNextTuple(root_id, in_cd);
        }
      }

      // insert line segments
      for (const SegmentType& edge : constraints)
      {
        const vtkIdType line[2] = { old_to_new_pt_ids[edge.first], old_to_new_pt_ids[edge.second] };
        if (line[0] < 0 || line[1] < 0)
        {
          //__debugbreak(); // triangle with this constraint was rejected?
          continue;
        }
        else if (line[0] == line[1])
        {
          // __debugbreak(); // degenerate
          continue;
        }

        out_lines->InsertNextCell(2, line);
        out_lines_cd->InsertNextTuple(root_id, in_cd);
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
      // throw std::logic_error("The method or operation is not implemented.");

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
        // find an edge vi--vj that intersects vl--vm
        for (const auto& edge : tsc_detail::TRIEDGES)
        {
          const vtkIdType& id0 = edge.first;
          const vtkIdType& id1 = edge.second;
          const vtkIdType& vi = point_ids[id0];
          const vtkIdType& vj = point_ids[id1];
          double pi[3] = { points[vi][0], points[vi][1], 0.0 };
          double pj[3] = { points[vj][0], points[vj][1], 0.0 };

          double px[3] = {};
          const auto intersect_type = tsc_detail::robustIntersect(pi, pj, pl, pm, px, tol);
          if (intersect_type != tsc_detail::IntersectType::Intersect)
            continue;

          v1 = vi;
          v2 = vj;
          vtkNew<vtkIdList> neighbors_v1v2;
          mesh->GetCellEdgeNeighbors(t0, v1, v2, neighbors_v1v2);
          if (!neighbors_v1v2->GetNumberOfIds()) // shouldn't happen, but eh just to be safe
            continue;
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
        for (const auto& edge : tsc_detail::TRIEDGES)
        {
          const vtkIdType& id0 = edge.first;
          const vtkIdType& id1 = edge.second;
          if ((point_ids[id0] == v1 && point_ids[id1] == v2) ||
            (point_ids[id1] == v1 && point_ids[id0] == v2))
          {
            vm_ = point_ids[tsc_detail::TRIOPPVERTS[id0]];
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
    SurfCutHelper(const LoopsInfoType& loops_info_holder_, const bool& inside_out_,
      vtkPointData* in_pd_, vtkPointData* out_pd_, vtkCellArray* out_verts_,
      vtkCellArray* out_lines_, vtkCellArray* out_polys_, vtkCellArray* out_strips_,
      vtkCellData* in_cd_, vtkCellData* out_verts_cd_, vtkCellData* out_lines_cd_,
      vtkCellData* out_polys_cd_, vtkCellData* out_strips_cd_, vtkIncrementalPointLocator* locator_,
      vtkPoints* in_points_, vtkPoints* in_loop_points_, vtkUnsignedCharArray* acquisition_)
      : loops_info_holder(loops_info_holder_)
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
      , in_loop_points(in_loop_points_)
      , acquisition(acquisition_)
    {
      in_acquisition = vtkUnsignedCharArray::SafeDownCast(in_pd->GetArray("Acquired"));
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
      PopTrisImpl worker;
      vtkDataArray* points1_arr = parent->Points->GetData();
      vtkDataArray* points2_arr = in_loop_points->GetData();
      if (!dispatchRR::Execute(points1_arr, points2_arr, worker, inside_out, loops_info_holder,
            is_acquired, constraints, parent, in_pd, out_pd, out_verts, out_lines, out_polys,
            out_strips, in_cd, out_verts_cd, out_lines_cd, out_polys_cd, out_strips_cd, locator,
            acquisition, fallthrough, cell_type))
        worker(points1_arr, points2_arr, inside_out, loops_info_holder, is_acquired, constraints,
          parent, in_pd, out_pd, out_verts, out_lines, out_polys, out_strips, in_cd, out_verts_cd,
          out_lines_cd, out_polys_cd, out_strips_cd, locator, acquisition, fallthrough, cell_type);
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

        // superpose input's acquisition
        if (in_acquisition)
        {
          if (in_acquisition->GetTypedComponent(point_ids->GetId(i), 0))
          {
            is_acquired.insert(i);
          }
        }
      }

      vtkDataArray* points1_arr = parent->Points->GetData();
      vtkDataArray* points2_arr = in_points->GetData();
      GetRootImpl worker;
      if (!dispatchRR::Execute(points1_arr, points2_arr, worker, parent->PointIds))
        worker(points1_arr, points2_arr, parent->PointIds);
    }

    /**
     * @brief Reset helper's internal state in anticipation for a new input
     *
     */
    void reset()
    {
      constraints.clear();
      is_acquired.clear();
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
                continue;

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
                if (num_iters > 32)
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
     * @brief Do triangle-line intersect tests in a plane @z=0.0
     *
     * @param lines a bunch of line segments
     * @param tol numeric tolerance for intersection math
     */
    void triIntersect(const SegmentsType& lines, const double& tol)
    {
      TriIntersect2dImpl worker;
      vtkDataArray* points1_arr = parent->Points->GetData();
      vtkDataArray* points2_arr = in_loop_points->GetData();
      if (!dispatchRR::Execute(
            points1_arr, points2_arr, worker, lines, constraints, is_acquired, tol))
        worker(points1_arr, points2_arr, lines, constraints, is_acquired, tol);
    }

    /**
     * @brief If child triangles exist, bring them into output data structure.
     *        else, copy parent triangle.
     *
     */
    inline void update() { parent->UpdateChildren(polys); }

    vtkNew<Parent> parent;
    SegmentsType constraints;
    std::unordered_set<vtkIdType> is_acquired;
    const LoopsInfoType& loops_info_holder;
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
    vtkPoints *in_points = nullptr, *in_loop_points = nullptr;
    vtkUnsignedCharArray* acquisition = nullptr;
    vtkUnsignedCharArray* in_acquisition = nullptr;
  };

  Child::Child()
    : cx(0)
    , cy(0)
    , bbox()
  {
  }

  Child::Child(const double& cx_, const double& cy_, const vtkIdType* pts, const vtkIdType npts_,
    const double bounds[4])
    : cx(cx_)
    , cy(cy_)
  {
    point_ids->SetNumberOfIds(npts_);
    std::copy(pts, pts + npts_, point_ids->GetPointer(0));
    bbox = vtkBoundingBox(bounds[0], bounds[1], bounds[2], bounds[3], 0.0, 0.0);
  }

  Child::Child(const double& cx_, const double& cy_, vtkIdList* pts, const double bounds[4])
    : cx(cx_)
    , cy(cy_)
  {
    point_ids->DeepCopy(pts);
    bbox = vtkBoundingBox(bounds[0], bounds[1], bounds[2], bounds[3], 0.0, 0.0);
  }

  vtkStandardNewMacro(Parent);
  Parent::Parent()
    : children()
    , cellId(-1)
  {
    this->children.reserve(10);
  }
  Parent::~Parent() { this->children.clear(); }
  void Parent::Reset()
  {
    this->children.clear();
    this->PointIds->Reset();
    this->Points->Reset();
  }
  void Parent::UpdateChildren(vtkCellArray* polys)
  {
    GetChildrenImpl worker;
    vtkDataArray* points_arr = this->Points->GetData();
    if (!dispatchR::Execute(points_arr, worker, polys, this->children))
      worker(points_arr, polys, this->children);
  }
}

int tscTriSurfaceCutter::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkSmartPointer<vtkPolyData> input =
    vtkPolyData::GetData(inputVector[0]->GetInformationObject(0));
  vtkSmartPointer<vtkPolyData> loops =
    vtkPolyData::GetData(inputVector[1]->GetInformationObject(0));

  vtkSmartPointer<vtkPolyData> output = vtkPolyData::GetData(outputVector->GetInformationObject(0));

  if (!this->Embed && !this->Remove) // do nothing
  {
    output->DeepCopy(input);
    return 1;
  }

  ///> Allocate output data structures.
  vtkSmartPointer<vtkPoints> in_points = input->GetPoints();
  vtkSmartPointer<vtkPoints> in_loop_points = loops->GetPoints();
  vtkSmartPointer<vtkCellArray> in_loop_polys = loops->GetPolys();

  vtkNew<vtkPoints> out_pts;
  vtkNew<vtkCellArray> out_verts, out_lines, out_polys, out_strips;
  vtkNew<vtkCellData> out_verts_cd, out_lines_cd, out_polys_cd, out_strips_cd;

  const vtkIdType& num_points = input->GetNumberOfPoints();
  const vtkIdType& num_loops_points = in_loop_points->GetNumberOfPoints();
  const vtkIdType& num_cells = input->GetNumberOfCells();
  const vtkIdType& num_verts = input->GetNumberOfVerts();
  const vtkIdType& num_strips = input->GetNumberOfStrips();
  const vtkIdType& num_polys = input->GetNumberOfPolys();
  const vtkIdType& num_loop_polys = loops->GetNumberOfPolys();

  if (!(num_points && num_polys))
  {
    vtkErrorMacro(<< "Input mesh is empty.");
    return 1;
  }
  if (!(num_loops_points && num_loop_polys))
  {
    vtkErrorMacro(<< "Input loops are empty.");
    return 1;
  }

  vtkSmartPointer<vtkPointData> in_pd = input->GetPointData();
  vtkSmartPointer<vtkPointData> out_pd = output->GetPointData();
  vtkSmartPointer<vtkCellData> in_cd = input->GetCellData();

  out_verts->Allocate(num_verts);
  vtkDebugMacro(<< "Alloc'd " << num_verts << " verts");

  out_lines->Allocate(num_loops_points - 1 + (num_polys >> 4));
  vtkDebugMacro(<< "Alloc'd " << num_loops_points - 1 + (num_polys >> 4) << " lines");

  out_polys->Allocate(num_polys + num_loops_points - 1);
  vtkDebugMacro(<< "Alloc'd " << num_polys + num_loops_points - 1 << " polys");

  out_strips->Allocate(num_strips);
  vtkDebugMacro(<< "Alloc'd " << num_verts << " strips");

  out_pts->SetDataType(in_points->GetDataType());
  out_pts->Allocate(num_points + num_loops_points);
  vtkDebugMacro(<< "Alloc'd " << num_points + num_loops_points << " points");

  if (in_pd->HasArray("Acquired"))
  {
    vtkLog(INFO, "Superposing acquisition.");
  }
  out_pd->InterpolateAllocate(in_pd, num_points + num_loops_points * 3);
  vtkDebugMacro(<< "Alloc'd " << num_points + num_loops_points * 3 << " point data tuples");

  if (in_cd->HasArray("Constrained"))
  {
    in_cd->RemoveArray("Constrained");
  }

  out_verts_cd->CopyAllocate(in_cd);
  vtkDebugMacro(<< "Alloc'd CellData for verts");

  out_lines_cd->CopyAllocate(in_cd);
  vtkDebugMacro(<< "Alloc'd CellData for lines");

  out_polys_cd->CopyAllocate(in_cd);
  vtkDebugMacro(<< "Alloc'd CellData for polys");

  out_strips_cd->CopyAllocate(in_cd);
  vtkDebugMacro(<< "Alloc'd CellData for strips");

  vtkNew<vtkUnsignedCharArray> acquisition;
  acquisition->SetName("Acquired");
  acquisition->SetNumberOfComponents(1);
  acquisition->Allocate(num_points + num_loops_points);
  vtkDebugMacro(<< "Alloc'd " << num_points + num_loops_points << " towards acquistion");
  ///> Finished allocation

  ///> Extend bounds along Z.
  double in_bounds[6] = {};
  double loops_bounds[6] = {};
  double out_bounds[6] = {};
  in_points->GetBounds(in_bounds);
  in_loop_points->GetBounds(loops_bounds);

  // extend in z to the combined z-bounds of surface and 2d loops
  {
    const double& zmin = std::min(in_bounds[4], loops_bounds[4]);
    const double& zmax = std::max(in_bounds[5], loops_bounds[5]);

    std::copy(in_bounds, in_bounds + 4, out_bounds);
    out_bounds[4] = zmin ? zmin : -10; // zmin, zmax might be zero.
    out_bounds[5] = zmax ? zmax : 10;
  }

  vtkDebugMacro(<< "Output bounds");
  vtkDebugMacro(<< "XMin: " << out_bounds[0] << ", XMax: " << out_bounds[1]);
  vtkDebugMacro(<< "YMin: " << out_bounds[2] << ", YMax: " << out_bounds[3]);
  vtkDebugMacro(<< "ZMin: " << out_bounds[4] << ", ZMax: " << out_bounds[5]);
  ///> Finished extending bounds.

  ///> Find candidate cells
  this->CreateDefaultLocators();
  this->CellLocator->CacheCellBoundsOn();
  this->CellLocator->SetDataSet(input);
  this->CellLocator->BuildLocator();
  vtkDebugMacro(<< "Built cell locator");

  // cache loops and cells that might cross.
  std::unordered_map<vtkIdType, SegmentsType> possible_crossings;
  LoopsInfoType loops_info_holder(num_loop_polys);
  possible_crossings.reserve(num_polys);

  auto loops_iter = vtk::TakeSmartPointer(loops->NewCellIterator());
  vtkIdType loop_id{ 0 };
  vtkNew<vtkGenericCell> loop;
  vtkNew<vtkIdList> cells;
  for (loops_iter->InitTraversal(); !loops_iter->IsDoneWithTraversal();
       loops_iter->GoToNextCell(), ++loop_id)
  {
    loops_iter->GetCell(loop);

    vtkSmartPointer<vtkIdList> loop_pts = loops_info_holder[loop_id].second;
    loop_pts->DeepCopy(loop->GetPointIds());

    const vtkIdType& num_edges = loop->GetNumberOfEdges();
    vtkDebugMacro(<< "Processing loop: " << loop_id << " (" << num_edges << " edges)");

    for (vtkIdType edge_id = 0; edge_id < num_edges; ++edge_id)
    {
      vtkCell* edge = loop->GetEdge(edge_id);
      const vtkIdType& e0 = edge->GetPointId(0);
      const vtkIdType& e1 = edge->GetPointId(1);

      double edge_bounds[6] = {};
      edge->GetBounds(edge_bounds);

      // extend up to combined {min, max} since that is what CellLocator sees.
      std::copy(out_bounds + 4, out_bounds + 6, edge_bounds + 4);

      vtkDebugMacro("Loop edge: " << edge_id << '(' << e0 << ',' << e1 << ')');
      vtkDebugMacro(<< " Xmin:" << edge_bounds[0] << " Xmax:" << edge_bounds[1]
                    << " Ymin:" << edge_bounds[2] << " Ymax:" << edge_bounds[3]
                    << " Zmin:" << edge_bounds[4] << " Zmax:" << edge_bounds[5]);
      this->CellLocator->FindCellsWithinBounds(edge_bounds, cells);
      vtkDebugMacro(<< "Possible crossings: " << cells->GetNumberOfIds());

      for (const auto& cell_id : *cells)
      {
        if (input->GetCellType(cell_id) == VTK_TRIANGLE)
        {
          possible_crossings[cell_id].emplace_back(e0, e1);
        }
      }
    }

    double loop_bounds[6] = {};
    loops->GetCellBounds(loop_id, loop_bounds);
    memset(loop_bounds + 4, 0, 2 * sizeof(double));
    loops_info_holder[loop_id].first = vtkBoundingBox(loop_bounds);
  }
  vtkDebugMacro(<< "Obtained candidate cells.");
  ///> Finished obtaining candidate cells.

  ///> Embed loops in candidate cells and remove if needed.
  // strategy for triangles: 1. push triangle into helper,
  //                         2. try to intersect with loops,
  //                         3. pop child triangles into output poly-data
  //                         4. reset helper's data structures to prepare for next triangle.
  tsc_detail::SurfCutHelper helper(loops_info_holder, this->InsideOut, in_pd, out_pd, out_verts,
    out_lines, out_polys, out_strips, in_cd, out_verts_cd, out_lines_cd, out_polys_cd,
    out_strips_cd, this->PointLocator, in_points, in_loop_points, acquisition);

  // roughly, a quarter no. of cells
  int report_every =
    (num_cells >= 4) ? num_cells >> 2 : (num_cells >= 2 ? num_cells >> 1 : num_cells);
  vtkDebugMacro(<< "Report progress every " << report_every << " cells");

  this->PointLocator->SetTolerance(this->Tolerance);
  this->PointLocator->InitPointInsertion(out_pts, out_bounds);
  vtkDebugMacro(<< "Init'd point insertion");

  auto input_iter = vtk::TakeSmartPointer(input->NewCellIterator());
  auto cell_pt_ids = vtkSmartPointer<vtkIdList>::New();
  for (input_iter->InitTraversal(); !input_iter->IsDoneWithTraversal(); input_iter->GoToNextCell())
  {
    const vtkIdType& cell_id = input_iter->GetCellId();
    const int& cell_type = input->GetCellType(cell_id);
    vtkDebugMacro(<< "Processing " << cell_id);

    // provide
    cell_pt_ids = input_iter->GetPointIds();
    helper.push(cell_pt_ids, cell_id, cell_type);

    switch (cell_type)
    {
      case VTK_TRIANGLE:
      {
        // embed
        const auto& trials = possible_crossings.find(cell_id);
        if ((trials != possible_crossings.end()) && this->Embed)
        {
          const SegmentsType& edges = trials->second;
          vtkDebugMacro(<< "Crosses " << edges.size() << " edges");

          vtkDebugMacro(<< "Intersect, tol: " << this->Tolerance);
          helper.triIntersect(edges, this->Tolerance);

          vtkDebugMacro(<< "Triangulate, tol: " << this->Tolerance);
          helper.triangulate(this->Tolerance);
        }
        break;
      }
      default:
        break;
    }
    helper.update();

    // accept/reject
    vtkDebugMacro(<< "Popping .. ");
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
      this->UpdateProgress(static_cast<double>(cell_id) / num_cells);
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
  if (this->ColorAcquiredPts)
  {
    acquisition->Squeeze();
    out_pd->AddArray(acquisition);
  }
  vtkDebugMacro("Finalized output point attributes. NumPoints: " << out_pts->GetNumberOfPoints());

  const vtkIdType& num_out_lines = output->GetNumberOfLines();
  const vtkIdType& num_out_polys = output->GetNumberOfPolys();
  const vtkIdType& num_out_cells = num_verts + num_out_lines + num_out_polys + num_strips;

  const vtkIdType nvl = num_verts + num_out_lines;
  const vtkIdType nvlp = nvl + num_out_polys;
  const vtkIdType nvlps = nvlp + num_strips;

  vtkSmartPointer<vtkCellData> out_cd = output->GetCellData();
  out_lines_cd->Squeeze();
  out_polys_cd->Squeeze();
  out_cd->CopyAllocate(out_polys_cd, num_out_cells);
  out_cd->CopyData(out_verts_cd, 0, num_verts, 0);
  out_cd->CopyData(out_lines_cd, num_verts, num_out_lines, 0);
  out_cd->CopyData(out_polys_cd, nvl, num_out_polys, 0);
  out_cd->CopyData(out_strips_cd, nvlp, 0);
  out_cd->Squeeze();

  if (this->ColorLoopEdges)
  {
    vtkNew<vtkUnsignedCharArray> constrained;
    constrained->SetName("Constrained");
    constrained->SetNumberOfComponents(1);
    constrained->SetNumberOfTuples(num_out_cells);

    for (vtkIdType tupIdx = 0; tupIdx < num_verts; ++tupIdx)
      constrained->SetTypedComponent(tupIdx, 0, '\000');

    for (vtkIdType tupIdx = num_verts; tupIdx < nvl; ++tupIdx)
      constrained->SetTypedComponent(tupIdx, 0, '\001');

    for (vtkIdType tupIdx = nvl; tupIdx < nvlp; ++tupIdx)
      constrained->SetTypedComponent(tupIdx, 0, '\000');

    for (vtkIdType tupIdx = nvlp; tupIdx < nvlps; ++tupIdx)
      constrained->SetTypedComponent(tupIdx, 0, '\000');

    out_cd->SetScalars(constrained);
  }
  vtkDebugMacro("Finalized output cell attributes.");
  vtkDebugMacro("Verts: " << num_verts << " Lines: " << num_out_lines << " Polys: " << num_out_polys
                          << " Strips: " << num_strips << "Cells: " << num_out_cells);
  ///> Finished finalizing output data structures.

  return 1;
}

void tscTriSurfaceCutter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "AccelerateCellLocator: " << (this->AccelerateCellLocator ? "True" : "False")
     << "\n";
  os << indent << "ColorAcquiredPts: " << (this->ColorAcquiredPts ? "True" : "False") << "\n";
  os << indent << "ColorLoopEdges  : " << (this->ColorLoopEdges ? "True" : "False") << "\n";
  os << indent << "InsideOut: " << (this->InsideOut ? "True" : "False") << "\n";
  os << indent << "Tolerance: " << this->Tolerance << "\n";
  os << indent << "Tolerance2: " << this->Tolerance * this->Tolerance << "\n";
  os << indent << "Embed: " << (this->Embed ? "True" : "False") << "\n";
  os << indent << "Remove: " << (this->Remove ? "True" : "False") << "\n";

  if (this->CellLocator)
  {
    os << indent << "Cell Locator: " << this->CellLocator << "\n";
    this->CellLocator->PrintSelf(os, indent);
  }
  else
  {
    os << indent << "Cell Locator: (None)\n";
  }
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
  if (!this->CellLocator)
  {
    if (this->AccelerateCellLocator)
      this->CellLocator = vtkSmartPointer<vtkStaticCellLocator>::New();
    else
      this->CellLocator = vtkSmartPointer<vtkCellLocator>::New();
  }
  if (!this->PointLocator)
  {
    this->PointLocator = vtkSmartPointer<vtkMergePoints>::New();
  }
}
