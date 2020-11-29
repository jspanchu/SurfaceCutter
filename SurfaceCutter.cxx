#include "SurfaceCutter.h"

#include <algorithm>
#include <iterator>
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
#include <vtkCellLocator.h>
#include <vtkDataArrayRange.h>
#include <vtkDelaunay2D.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLine.h>
#include <vtkMergePoints.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>

using dispatchRRI = vtkArrayDispatch::Dispatch3ByValueType<vtkArrayDispatch::Reals,
  vtkArrayDispatch::Reals, vtkArrayDispatch::Integrals>;
using dispatchRR =
vtkArrayDispatch::Dispatch2ByValueType<vtkArrayDispatch::Reals, vtkArrayDispatch::Reals>;
using dispatchR = vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals>;

using SegmentsType = std::vector<std::pair<vtkIdType, vtkIdType>>; // {{p1, p2}, {p2, p3}, ...}

vtkStandardNewMacro(SurfaceCutter);

SurfaceCutter::SurfaceCutter()
  : ColorAcquiredPts(true)
  , Tolerance(1.0e-6)
  , ColorLoopEdges(true)
  , InsideOut(true) // default: remove portions outside loop polygons.
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);

  this->CreateDefaultLocators();

  vtkDebugMacro(<< "Initialized " << this->GetClassNameInternal());
}

SurfaceCutter::~SurfaceCutter()
{
  vtkDebugMacro(<< "Destroyed " << this->GetClassNameInternal());
}

void SurfaceCutter::SetLoops(vtkPointSet* loops)
{
  this->SetInputData(1, loops);
}

void SurfaceCutter::SetLoopsConnection(vtkAlgorithmOutput* output)
{
  this->SetInputConnection(1, output);
}

// anon begin
namespace
{
  enum class BaryCentricType
  {
    OnVertex,
    OnEdge,
    Inside,
    Outside,
    Degenerate
  };

  enum class IntersectType
  {
    PerfectCross,
    TJunction,
    NoIntersection,
    OnLine
  };

  // Convention: in a triangle- v: vertex, e: edge
  // with vertices v0--v1--v2,
  // e0 = v0--v1, e1 = v1--v2, e2 = v2--v0;
  const vtkIdType TRIEDGES[3][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
  const vtkIdType TRIOPPEDGE[3] = { 1, 2, 0 };  // edge id opposite to a vertex
  const vtkIdType TRIOPPVERTS[3] = { 2, 0, 1 }; // vertex id opposite to an edge
  const vtkIdType TRIVERTS[3] = { 0, 1, 2 };

  // Find if a point is on an edge of a triangle? Inside/outside ?
  // Or exactly coincident with a vertex?
  BaryCentricType inTriangle(const double p[3], const double p0[3], const double p1[3],
    const double p2[3], double bCoords[3], const double& tol, vtkIdType& v, vtkIdType& e)
  {
    v = -1;
    e = -1;
    if (!vtkTriangle::BarycentricCoords(p, p0, p1, p2, bCoords))
      return BaryCentricType::Degenerate;

    const double& s = bCoords[0];
    const double& t = bCoords[1];
    const double& st = bCoords[2]; // 1 - s - t

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
      return BaryCentricType::OnVertex;
    }

    bool s_gt0_lt1 = !std::signbit(s) && std::isless(s, 1) && !s_eq_0;
    bool t_gt0_lt1 = !std::signbit(t) && std::isless(t, 1) && !t_eq_0;
    bool st_gt0_lt1 = !std::signbit(st) && std::isless(st, 1) && !st_eq_0;

    if (s_gt0_lt1 && t_gt0_lt1 && st_gt0_lt1)
      return BaryCentricType::Inside;

    if (s_eq_0 && t_gt0_lt1 && st_gt0_lt1)
      e = TRIOPPEDGE[0];
    else if (t_eq_0 && s_gt0_lt1 && st_gt0_lt1)
      e = TRIOPPEDGE[1];
    else if (st_eq_0 && s_gt0_lt1 && t_gt0_lt1)
      e = TRIOPPEDGE[2];
    if (e >= 0)
      return BaryCentricType::OnEdge;

    if (s_eq_0 && t_eq_0 && st_eq_0)
      return BaryCentricType::Degenerate;

    return BaryCentricType::Outside;
  }

  // Reasonable z-value of a point, inside of a triangle / on an edge.
  inline void interpZ(double p[3], const double p0[3], const double p1[3], const double p2[3])
  {
    double bCoords[3] = {};
    vtkTriangle::BarycentricCoords(p, p0, p1, p2, bCoords);
    p[2] = bCoords[0] * p0[2] + bCoords[1] * p1[2] + bCoords[2] * p2[2];
  }

  // Same as above, but when you already have bCoords
  inline void interpZ(
    double p[3], const double p0[3], const double p1[3], const double p2[3], double bCoords[3])
  {
    p[2] = bCoords[0] * p0[2] + bCoords[1] * p1[2] + bCoords[2] * p2[2];
  }

  // Quantify return status of vtkLine::Intersection(...)
  IntersectType robustIntersect(const double* segP1, const double* segP2, const double* a1,
    const double* a2, double px[3], const double& tol)
  {
    double u(0.0), v(0.0);
    int intersect = vtkLine::Intersection(segP1, segP2, a1, a2, u, v);

    double uProj[3] = { 0., 0., 0. };
    double vProj[3] = { 0., 0., 0. };

    for (unsigned short dim = 0; dim < 2; ++dim) // discard z
    {
      uProj[dim] = segP1[dim] + u * (segP2[dim] - segP1[dim]);
      vProj[dim] = a1[dim] + v * (a2[dim] - a1[dim]);
    }

    double closestDist = vtkMath::Distance2BetweenPoints(uProj, vProj);
    double tol2 = tol * tol;

    // Note: this section is dependent on consts def'd in vtkLine.cxx.
    // Keep in mind to update this when that changes.
    // Just so we're clear:
    // 0: VTK_NO_INTERSECTION
    // 2: VTK_YES_INTERSECTION
    // 3: VTK_ON_LINE
    if (intersect == 0)
      return IntersectType::NoIntersection;
    else if (intersect == 2 && closestDist < tol2)
    {
      if (std::fabs(1.0 - u) <= tol && closestDist < tol2)
      {
        std::copy(segP2, segP2 + 3, px);
        return IntersectType::TJunction;
      }
      else if (std::fabs(1.0 - v) <= tol && closestDist < tol2)
      {
        std::copy(a2, a2 + 3, px);
        return IntersectType::TJunction;
      }

      if (std::fpclassify(u) == FP_ZERO && closestDist < tol2)
      {
        std::copy(segP1, segP1 + 3, px);
        return IntersectType::TJunction;
      }
      else if (std::fpclassify(v) == FP_ZERO && closestDist < tol2)
      {
        std::copy(a1, a1 + 3, px);
        return IntersectType::TJunction;
      }

      if (std::fabs(u) <= tol && closestDist < tol2)
      {
        std::copy(segP1, segP1 + 3, px);
        return IntersectType::TJunction;
      }
      else if (std::fabs(v) <= tol && closestDist < tol2)
      {
        std::copy(a1, a1 + 3, px);
        return IntersectType::TJunction;
      }
      else
      {
        std::copy(vProj, vProj + 3, px);
        return IntersectType::PerfectCross;
      }
    }
    else if (intersect == 3 && closestDist < tol2)
    {
      for (unsigned short dim = 0; dim < 3; ++dim)
        px[dim] = a1[dim] + v * (a2[dim] - a1[dim]);

      return IntersectType::OnLine;
    }

    return IntersectType::NoIntersection;
  }

  struct GetEdgeBBoxImpl
  {
    // strict 2d bbox
    template <typename PointsArrT, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT>>
    void operator()(
      PointsArrT* pointsArr, const vtkIdType& v1, const vtkIdType& v2, vtkBoundingBox& bbox)
    {
      auto points = vtk::DataArrayTupleRange<3>(pointsArr);
      using PointsT = vtk::GetAPIType<PointsArrT>;
      PointsT xmin = (points[v1][0] < points[v2][0]) ? points[v1][0] : points[v2][0];
      PointsT xmax = (points[v1][0] > points[v2][0]) ? points[v1][0] : points[v2][0];
      PointsT ymin = (points[v1][1] < points[v2][1]) ? points[v1][1] : points[v2][1];
      PointsT ymax = (points[v1][1] > points[v2][1]) ? points[v1][1] : points[v2][1];
      bbox.SetBounds(xmin, xmax, ymin, ymax, 0.0, 0.0);
    }

    // overload to enforce a certain {zmin, zmax}
    template <typename PointsArrT, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT>>
    void operator()(PointsArrT* pointsArr, const vtkIdType& v1, const vtkIdType& v2,
      const double& zmin, const double& zmax, vtkBoundingBox& bbox)
    {
      auto points = vtk::DataArrayTupleRange<3>(pointsArr);
      using PointsT = vtk::GetAPIType<PointsArrT>;
      PointsT xmin = (points[v1][0] < points[v2][0]) ? points[v1][0] : points[v2][0];
      PointsT xmax = (points[v1][0] > points[v2][0]) ? points[v1][0] : points[v2][0];
      PointsT ymin = (points[v1][1] < points[v2][1]) ? points[v1][1] : points[v2][1];
      PointsT ymax = (points[v1][1] > points[v2][1]) ? points[v1][1] : points[v2][1];
      bbox.SetBounds(xmin, xmax, ymin, ymax, zmin, zmax);
    }
  };

  // A child is born when a parent triangle crosses a loop's edge.
  struct Child
  {
    Child()
      : cx(0)
      , cy(0)
      , npts(0)
      , pts()
      , bbox()
    {
      pts.reserve(3);
    }
    Child(const double& cx_, const double& cy_, const vtkIdType* pts_, const vtkIdType& npts_,
      const double bounds[4])
      : cx(cx_)
      , cy(cy_)
      , npts(npts_)
    {
      pts.assign(pts_, pts_ + npts);
      bbox = vtkBoundingBox(bounds[0], bounds[1], bounds[2], bounds[3], 0.0, 0.0);
    }
    double cx, cy;
    vtkIdType npts;
    std::vector<vtkIdType> pts;
    vtkBoundingBox bbox;
  };

  struct GetChildrenImpl
  {
    template <typename PointsArrT, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT>>
    void operator()(PointsArrT* pointsArr, vtkCellArray* tris, std::vector<Child>& children)
    {
      const vtkIdType& ntris = tris->GetNumberOfCells();
      auto points = vtk::DataArrayTupleRange<3>(pointsArr);

      for (vtkIdType tri = 0; tri < ntris; ++tri)
      {
        const vtkIdType* pts = nullptr;
        vtkIdType npts(0);
        tris->GetCellAtId(tri, npts, pts);

        double bounds[4] = { VTK_DOUBLE_MAX, VTK_DOUBLE_MIN, VTK_DOUBLE_MAX, VTK_DOUBLE_MIN };
        double pc[3] = {};
        for (vtkIdType v = 0; v < npts; ++v)
        {
          const vtkIdType& pt = pts[v];
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
          pc[dim] /= 3.0;

        children.emplace_back(pc[0], pc[1], pts, npts, bounds);
      }
    }
  };

  struct Root
  {
    Root()
      : id(0)
    {
      children.reserve(10);
    }
    Root(const vtkIdType& id_)
      : id(id_)
    {
      children.reserve(10);
    }
    vtkNew<vtkTriangle> repr;
    std::vector<Child> children;
    vtkIdType id;

    void reset()
    {
      repr->Points->Resize(3); // preserve initial 3 points
      children.clear();
      id = -1;
    }

    void update(vtkCellArray* triangles)
    {
      GetChildrenImpl worker;
      vtkDataArray* pointsArr = repr->Points->GetData();
      if (!dispatchR::Execute(pointsArr, worker, triangles, children))
        worker(pointsArr, triangles, children);
    }
  };

  // fill in root triangle
  struct GetRootImpl
  {
    // pointsArr1: dataset points
    // pointsArr2: root triangle points
    template <typename PointsArr1T, typename = vtk::detail::EnableIfVtkDataArray<PointsArr1T>,
      typename PointsArr2T, typename = vtk::detail::EnableIfVtkDataArray<PointsArr2T>>
      void operator()(PointsArr1T* pointsArr1, PointsArr2T* pointsArr2, vtkIdList* rootVerts)
    {
      auto points1 = vtk::DataArrayTupleRange<3>(pointsArr1);
      auto points2 = vtk::DataArrayTupleRange<3>(pointsArr2);

      for (vtkIdType tupIdx = 0; tupIdx < 3; ++tupIdx)
      {
        const vtkIdType& pt = rootVerts->GetId(tupIdx);
        std::copy(points2[pt].begin(), points2[pt].end(), points1[tupIdx].begin());
      }
    }
  };

  struct TriIntersect2dImpl
  {
    // pointsArr1: root triangle points
    // pointsArr2: line points
    template <typename PointsArr1T, typename = vtk::detail::EnableIfVtkDataArray<PointsArr1T>,
      typename PointsArr2T, typename = vtk::detail::EnableIfVtkDataArray<PointsArr1T>>
      void operator()(PointsArr1T* pointsArr1, PointsArr2T* pointsArr2, const SegmentsType& lines,
        SegmentsType& constraints, std::unordered_set<vtkIdType>& isAcquired, const double& tol)
    {
      if (!lines.size())
        return;

      auto points1 = vtk::DataArrayTupleRange<3>(pointsArr1);
      auto points2 = vtk::DataArrayTupleRange<3>(pointsArr2);

      using PointsT1 = vtk::GetAPIType<PointsArr1T>;
      using PointsT2 = vtk::GetAPIType<PointsArr2T>;

      double p2d[3][3] = { { points1[0][0], points1[0][1], 0.0 },
        { points1[1][0], points1[1][1], 0.0 }, { points1[2][0], points1[2][1], 0.0 } };

      double p3d[3][3] = { { points1[0][0], points1[0][1], points1[0][2] },
        { points1[1][0], points1[1][1], points1[1][2] },
        { points1[2][0], points1[2][1], points1[2][2] } };

      std::unordered_map<vtkIdType, vtkIdType> processed; // k: line pt, v: insertLoc
      const double tol2 = tol * tol;

      for (const auto& line : lines)
      {
        const vtkIdType& l1 = line.first;
        const vtkIdType& l2 = line.second;

        double segP1[3] = { points2[l1][0], points2[l1][1], 0.0 };
        double segP2[3] = { points2[l2][0], points2[l2][1], 0.0 };

        double segP1_3d[3] = { points2[l1][0], points2[l1][1], points2[l1][2] };
        double segP2_3d[3] = { points2[l2][0], points2[l2][1], points2[l2][2] };

        vtkIdType l1v(-1), l2v(-1);
        vtkIdType l1e(-1), l2e(-1);

        double bCoords1[3] = {};
        double bCoords2[3] = {};

        const auto l1Pos = inTriangle(segP1, p2d[0], p2d[1], p2d[2], bCoords1, tol, l1v, l1e);
        const auto l2Pos = inTriangle(segP2, p2d[0], p2d[1], p2d[2], bCoords2, tol, l2v, l2e);

        const bool l1Inside = l1Pos == BaryCentricType::Inside;
        const bool l2Inside = l2Pos == BaryCentricType::Inside;
        const bool l1OnEdge = l1Pos == BaryCentricType::OnEdge;
        const bool l2OnEdge = l2Pos == BaryCentricType::OnEdge;
        const bool l1OnVert = l1Pos == BaryCentricType::OnVertex;
        const bool l2OnVert = l2Pos == BaryCentricType::OnVertex;
        const bool l1Outsid = l1Pos == BaryCentricType::Outside;
        const bool l2Outsid = l2Pos == BaryCentricType::Outside;

        double px[3] = {};
        std::set<vtkIdType> hits;
        vtkIdType inserted(-1);

        for (const auto& e : TRIEDGES)
        {
          const vtkIdType& v0 = e[0];
          const vtkIdType& v1 = e[1];

          const double* a1 = p2d[v0];
          const double* a2 = p2d[v1];

          const auto intrsctType = robustIntersect(segP1, segP2, a1, a2, px, tol);
          if (intrsctType == IntersectType::PerfectCross)
          {
            interpZ(px, p3d[0], p3d[1], p3d[2]);
            inserted = pointsArr1->InsertNextTuple(px);
            hits.insert(inserted);
          }
          else if (intrsctType == IntersectType::TJunction)
          {
            double minDist = VTK_DOUBLE_MAX;
            vtkIdType tgt(-1);
            for (const auto& v : TRIVERTS)
            {
              const double* a = p2d[v];
              double dist2 = vtkLine::DistanceToLine(a, segP1, segP2);
              if (dist2 < minDist)
              {
                minDist = dist2;
                tgt = v;
              }
            }
            if (tgt >= 0 && minDist <= tol2)
            {
              hits.insert(tgt);
              double dist2p1 = vtkMath::Distance2BetweenPoints(p2d[tgt], segP1);
              double dist2p2 = vtkMath::Distance2BetweenPoints(p2d[tgt], segP2);
              if (dist2p1 < tol2)
                isAcquired.insert(tgt);
              else if (dist2p2 < tol2)
                isAcquired.insert(tgt);
            }
            else if (l1OnEdge)
            {
              if (processed.find(l1) == processed.end())
              {
                std::copy(segP1, segP1 + 3, px);
                interpZ(px, p3d[0], p3d[1], p3d[2], bCoords1);
                inserted = pointsArr1->InsertNextTuple(px);
                processed[l1] = inserted;
                hits.insert(inserted);
                isAcquired.insert(inserted);
              }
            }
            else if (l2OnEdge)
            {
              if (processed.find(l2) == processed.end())
              {
                std::copy(segP2, segP2 + 3, px);
                interpZ(px, p3d[0], p3d[1], p3d[2], bCoords2);
                inserted = pointsArr1->InsertNextTuple(px);
                processed[l2] = inserted;
                hits.insert(inserted);
                isAcquired.insert(inserted);
              }
            }
          }
        }
        if (hits.size() == 2)
        {
          constraints.emplace_back(*(hits.begin()), *(std::next(hits.begin())));
        }
        else if (hits.size() == 1)
        {
          if (l1Inside || l1OnEdge)
          {
            if (processed.find(l1) == processed.end())
            {
              std::copy(segP1, segP1 + 3, px);
              interpZ(px, p3d[0], p3d[1], p3d[2], bCoords1);
              inserted = pointsArr1->InsertNextTuple(px);
              processed[l1] = inserted;
            }
            constraints.emplace_back(*(hits.begin()), processed[l1]);
            isAcquired.insert(processed[l1]);
          }
          else if (l1OnVert)
          {
            processed[l1] = l1v;
            constraints.emplace_back(*(hits.begin()), l1v);
            isAcquired.insert(l1v);
          }

          if (l2Inside || l2OnEdge)
          {
            if (processed.find(l2) == processed.end())
            {
              std::copy(segP2, segP2 + 3, px);
              interpZ(px, p3d[0], p3d[1], p3d[2], bCoords2);
              inserted = pointsArr1->InsertNextTuple(px);
              processed[l2] = inserted;
            }
            constraints.emplace_back(*(hits.begin()), processed[l2]);
            isAcquired.insert(processed[l2]);
          }
          else if (l2OnVert)
          {
            processed[l2] = l2v;
            constraints.emplace_back(*(hits.begin()), l2v);
            isAcquired.insert(l2v);
          }
        }
        else if (!hits.size())
        {
          // the last hope!
          bool force = (l1Inside && l2Inside) || (l1Inside && l2OnEdge) || (l1OnEdge && l2Inside) ||
            (l1OnEdge && l2OnEdge) || (l1Inside && l2OnVert) || (l1OnVert && l2Inside) ||
            (l1OnEdge && l2OnVert) || (l1OnVert && l2OnEdge) || (l1OnVert && l2OnVert);
          if (force)
          {
            if (processed.find(l1) == processed.end())
            {
              std::copy(segP1, segP1 + 3, px);
              interpZ(px, p3d[0], p3d[1], p3d[2], bCoords1);
              inserted = pointsArr1->InsertNextTuple(px);
              processed[l1] = inserted;
            }
            if (processed.find(l2) == processed.end())
            {
              std::copy(segP2, segP2 + 3, px);
              interpZ(px, p3d[0], p3d[1], p3d[2], bCoords2);
              inserted = pointsArr1->InsertNextTuple(px);
              processed[l2] = inserted;
            }
            constraints.emplace_back(processed[l1], processed[l2]);
            isAcquired.insert(processed[l1]);
            isAcquired.insert(processed[l2]);
          }
        }
      }
    }
  };

  struct PopTrisImpl
  {
    // pointsArr1: root triangle points
    // pointsArr2: loop points
    // inOutsArr : describes inside/out nature for each loop polygon
    // and a few other args ..
    ///> Remove triangles
    /// 1. inOut is set, reject triangles outside a polygon.
    /// 2. inOut is unset, reject triangles inside a polygon.
    template <typename PointsArr1T, typename = vtk::detail::EnableIfVtkDataArray<PointsArr1T>,
      typename PointsArr2T, typename = vtk::detail::EnableIfVtkDataArray<PointsArr2T>,
      typename InOutsArrT, typename = vtk::detail::EnableIfVtkDataArray<InOutsArrT>>
      void operator()(PointsArr1T* pointsArr1, PointsArr2T* pointsArr2, InOutsArrT* inOutsArr,
        std::vector<std::pair<vtkBoundingBox, vtkNew<vtkIdList>>>& loops,
        const std::unordered_set<vtkIdType>& isAcquired, const SegmentsType& constraints, Root& parent,
        vtkPointData* inPd, vtkPointData* outPd, vtkCellArray* outTris, vtkCellArray* outLines,
        vtkCellData* inCd, vtkCellData* outTriCd, vtkCellData* outLineCd,
        vtkIncrementalPointLocator* locator, vtkUnsignedCharArray* acquisition)
    {
      // throw std::logic_error("The method or operation is not implemented.");
      auto points1 = vtk::DataArrayTupleRange<3>(pointsArr1);
      auto points2 = vtk::DataArrayTupleRange<3>(pointsArr2);
      auto inOuts = vtk::DataArrayValueRange<1>(inOutsArr);

      using Points1T = vtk::GetAPIType<PointsArr1T>;
      using Points2T = vtk::GetAPIType<PointsArr2T>;
      using InOutsT = vtk::GetAPIType<InOutsArrT>;

      vtkTriangle* rootTri = parent.repr;
      vtkIdList* rootPts = rootTri->PointIds;
      const vtkIdType& numPoints = rootTri->GetPoints()->GetNumberOfPoints();
      auto& children = parent.children;
      const std::size_t numChildren = children.size();
      const bool childless = (numChildren == 1);
      const vtkIdType& rootId = parent.id;
      std::vector<vtkIdType> toNewIds(numPoints, -1);

      for (auto& child : children)
      {
        const double& testx = child.cx;
        const double& testy = child.cy;
        vtkBoundingBox& triBBox = child.bbox;
        const vtkIdType& npts = child.npts;
        const auto& pts = child.pts;

        bool rejected(false);
        vtkIdType loopId(-1);
        for (auto& loopInfo : loops)
        {
          ++loopId;
          vtkBoundingBox& loopBBox = loopInfo.first;
          vtkIdList* loopPts = loopInfo.second;
          InOutsT inOut = inOuts[loopId];

          rejected = false;

          // the easy-way; effective when inside out is unset. (non-default)
          if (!inOut)
          {
            if (!(loopBBox.Contains(triBBox) || loopBBox.Intersects(triBBox)))
            {
              rejected = false;
              continue;
            }
          }

          // the hard-way; default
          bool inside = false;
          vtkIdType nvert(loopPts->GetNumberOfIds());
          const vtkIdType* loopPtsPtr = loopPts->GetPointer(0);
          for (vtkIdType i = 0, j = nvert - 1; i < nvert; j = i++)
          {
            const double ix = points2[loopPtsPtr[i]][0];
            const double iy = points2[loopPtsPtr[i]][1];
            const double jx = points2[loopPtsPtr[j]][0];
            const double jy = points2[loopPtsPtr[j]][1];

            if (((iy > testy) != (jy > testy)) && (testx < (jx - ix) * (testy - iy) / (jy - iy) + ix))
              inside = !inside;
          }

          bool outside = !inside;
          if (inOut && outside)
          {
            rejected = true;
            break;
          }
          else if (!inOut && inside)
          {
            rejected = true;
            break;
          }
        }

        if (rejected)
          continue; // next triangle

        double dist2(0.);
        double* closest = nullptr;
        double pCoords[3] = {};
        double weights[3] = {};
        int subId(0);
        vtkIdType newPtId(-1);

        outTris->InsertNextCell(npts);
        for (const auto& pt : pts)
        {
          const double p[3] = { points1[pt][0], points1[pt][1], points1[pt][2] };
          if (locator->InsertUniquePoint(p, newPtId))
          {
            if (childless) // avoid computation of weights
              weights[0] = 1.0;
            else
              rootTri->EvaluatePosition(p, closest, subId, pCoords, dist2, weights);

            outPd->InterpolatePoint(inPd, newPtId, rootPts, weights);

            if (isAcquired.find(pt) != isAcquired.end())
              acquisition->InsertNextValue('\001');
            else
              acquisition->InsertNextValue('\000');
          }
          toNewIds[pt] = newPtId;
          outTris->InsertCellPoint(newPtId);
        }
        outTriCd->InsertNextTuple(rootId, inCd);
      }

      for (const auto& edge : constraints)
      {
        const vtkIdType line[2] = { toNewIds[edge.first], toNewIds[edge.second] };
        if (line[0] < 0 || line[1] < 0)
        {
          //__debugbreak(); // triangle with this constraint got rejected ?
          continue;
        }

        outLines->InsertNextCell(2, line);
        outLineCd->InsertNextTuple(rootId, inCd);
      }
    }
  };

  // vtkDelaunay2D's constraint section runs into infinite recursion for certain corner cases.
  // This implementation is no angel either. But, it does the job.
  struct ApplyConstraintImpl
  {
    // Try to enforce edge vl_--vm_ in mesh_
    ApplyConstraintImpl(
      vtkSmartPointer<vtkPolyData> mesh_, const vtkIdType& vl_, const vtkIdType& vm_)
      : mesh(mesh_)
      , vl(vl_)
      , vm(vm_)
    {
    }

    vtkSmartPointer<vtkPolyData> mesh;
    vtkIdType vl, vm;

    template <typename PointsArrT, typename = vtk::detail::EnableIfVtkDataArray<PointsArrT>>
    void operator()(PointsArrT* pointsArr, const double& tol)
    {
      // throw std::logic_error("The method or operation is not implemented.");

      auto points = vtk::DataArrayTupleRange<3>(pointsArr);

      const double pl[3] = { points[vl][0], points[vl][1], 0.0 };
      const double pm[3] = { points[vm][0], points[vm][1], 0.0 };

      vtkIdType* neisVl = nullptr;
      vtkIdType numNeisVl = 0;
      mesh->BuildLinks();
      mesh->GetPointCells(vl, numNeisVl, neisVl); // tris around vm;

      vtkIdType t0(-1), t1(-1), v1(-1), v2(-1), vm_(-1);
      // find a triangle st its edge opp to 'vl' intersects vl--vm
      bool abortScan(false);
      for (vtkIdType t = 0; t < numNeisVl && !abortScan; ++t)
      {
        t0 = neisVl[t];
        const vtkIdType* pts = nullptr;
        vtkIdType npts(0);
        mesh->GetCellPoints(t0, npts, pts);
        // find an edge vi--vj that intersects vl--vm
        for (const auto& edge : TRIEDGES)
        {
          const vtkIdType& id0 = edge[0];
          const vtkIdType& id1 = edge[1];
          const vtkIdType& vi = pts[id0];
          const vtkIdType& vj = pts[id1];
          const double pi[3] = { points[vi][0], points[vi][1], 0.0 };
          const double pj[3] = { points[vj][0], points[vj][1], 0.0 };

          double px[3] = {};
          const auto intrsctType = robustIntersect(pi, pj, pl, pm, px, tol);
          if (intrsctType != IntersectType::PerfectCross)
            continue;

          v1 = vi;
          v2 = vj;
          vtkNew<vtkIdList> v1v2Nei;
          mesh->GetCellEdgeNeighbors(t0, v1, v2, v1v2Nei);
          if (!v1v2Nei->GetNumberOfIds()) // shouldn't happen, but eh just to be safe
            continue;
          t1 = v1v2Nei->GetId(0);
          abortScan = true;
          break;
        }
      }

      if (abortScan)
      {
        const vtkIdType* pts = nullptr;
        vtkIdType npts(0);
        mesh->GetCellPoints(t1, npts, pts);
        for (const auto& edge : TRIEDGES)
        {
          const vtkIdType& id0 = edge[0];
          const vtkIdType& id1 = edge[1];
          if ((pts[id0] == v1 && pts[id1] == v2) || (pts[id1] == v1 && pts[id0] == v2))
          {
            vm_ = pts[TRIOPPVERTS[id0]];
            vtkIdType* neisVm_;
            vtkIdType numNeisVm_(0);
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

  // Combine all the above functors with this helper. This will be used in RequestData(...)
  struct SurfCutHelper
  {
    SurfCutHelper(vtkPointData* inPd_, vtkPointData* outPd_, vtkCellArray* outTris_,
      vtkCellArray* outLines_, vtkCellData* inCd_, vtkCellData* outTriCd_, vtkCellData* outLineCd_,
      vtkIncrementalPointLocator* locator_)
      : inPd(inPd_)
      , outPd(outPd_)
      , outTris(outTris_)
      , outLines(outLines_)
      , inCd(inCd_)
      , outTriCd(outTriCd_)
      , outLineCd(outLineCd_)
      , locator(locator_)
    {
      tris->InsertNextCell(3, TRIVERTS);
      in->SetPoints(parent.repr->Points);
      del2d->SetProjectionPlaneMode(VTK_DELAUNAY_XY_PLANE);
      del2d->SetInputData(in);
      del2d->SetOffset(100.0); // bump this if Delaunay output is concave
    }

    void pop(vtkPoints* points, vtkDataArray* insideOuts,
      std::vector<std::pair<vtkBoundingBox, vtkNew<vtkIdList>>>& loopsInf,
      vtkUnsignedCharArray* acquisition)
    {
      PopTrisImpl worker;
      vtkDataArray* pointsArr1 = parent.repr->Points->GetData();
      vtkDataArray* pointsArr2 = points->GetData();
      if (!dispatchRRI::Execute(pointsArr1, pointsArr2, insideOuts, worker, loopsInf, isAcquired,
        constraints, parent, inPd, outPd, outTris, outLines, inCd, outTriCd, outLineCd, locator,
        acquisition))
        worker(pointsArr1, pointsArr2, insideOuts, loopsInf, isAcquired, constraints, parent, inPd,
          outPd, outTris, outLines, inCd, outTriCd, outLineCd, locator, acquisition);
    }

    void push(vtkPoints* points, const vtkIdType* v0, const vtkIdType* v1, const vtkIdType* v2,
      const vtkIdType& rootIdx)
    {
      parent.id = rootIdx;
      parent.children.clear();

      parent.repr->PointIds->SetId(0, *v0);
      parent.repr->PointIds->SetId(1, *v1);
      parent.repr->PointIds->SetId(2, *v2);

      vtkPoints* rootPoints = parent.repr->Points;
      vtkDataArray* pointsArr1 = rootPoints->GetData();
      vtkDataArray* pointsArr2 = points->GetData();
      GetRootImpl worker;
      if (!dispatchRR::Execute(pointsArr1, pointsArr2, worker, parent.repr->PointIds))
        worker(pointsArr1, pointsArr2, parent.repr->PointIds);
    }

    void reset()
    {
      constraints.clear();
      isAcquired.clear();
      parent.reset();
      tris->Reset();
      tris->InsertNextCell(3, TRIVERTS);
    }

    void triangulate(const double& tol)
    {
      if (parent.repr->Points->GetNumberOfPoints() > 3)
      {
        // del2d->SetTolerance(tol);
        del2d->Update();

        if (del2d->GetOutput())
        {
          out->ShallowCopy(del2d->GetOutput());

          if (out->GetNumberOfPolys())
          {
            out->BuildLinks();
            for (const auto& edge : constraints)
            {
              if (edge.first == edge.second)
                continue;

              unsigned short niters = 0;
              while (!out->IsEdge(edge.first, edge.second))
              {
                auto constraintWorker = ApplyConstraintImpl(out, edge.first, edge.second);
                vtkDataArray* points = out->GetPoints()->GetData();
                if (!dispatchR::Execute(points, constraintWorker, tol))
                {
                  constraintWorker(points, tol);
                }
                ++niters;
                if (niters > 32)
                {
                  //__debugbreak(); // happens only for extremely degenerate inputs.
                  break;
                }
              }
            }
            tris->ShallowCopy(out->GetPolys());
          }
        }
      }
    }

    void triIntersect(vtkPoints* points, const SegmentsType& lines, const double& tol)
    {
      TriIntersect2dImpl worker;
      vtkPoints* rootPoints = parent.repr->Points;
      vtkDataArray* pointsArr1 = rootPoints->GetData();
      vtkDataArray* pointsArr2 = points->GetData();
      if (!dispatchRR::Execute(pointsArr1, pointsArr2, worker, lines, constraints, isAcquired, tol))
        worker(pointsArr1, pointsArr2, lines, constraints, isAcquired, tol);
    }

    inline void update() { parent.update(tris); }

    Root parent;
    SegmentsType constraints;
    std::unordered_set<vtkIdType> isAcquired;
    vtkNew<vtkCellArray> tris;
    vtkNew<vtkDelaunay2D> del2d;
    vtkNew<vtkPolyData> in, out;
    vtkCellArray* outTris, * outLines;
    vtkCellData* inCd, * outTriCd, * outLineCd;
    vtkIncrementalPointLocator* locator;
    vtkPointData* inPd, * outPd;
  };
}
// anon end
int SurfaceCutter::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkSmartPointer<vtkPolyData> input =
    vtkPolyData::GetData(inputVector[0]->GetInformationObject(0));
  vtkSmartPointer<vtkPolyData> loops_ =
    vtkPolyData::GetData(inputVector[1]->GetInformationObject(0));

  vtkSmartPointer<vtkPolyData> output = vtkPolyData::GetData(outputVector->GetInformationObject(0));

  // we're about to do stuff that is not thread-safe (query bounds, add scalars, ..).
  // multiple threads might use the same loop polydata for different surfaces,
  // so copy bare minimum
  vtkNew<vtkPolyData> loops;
  loops->CopyStructure(loops_);

  vtkSmartPointer<vtkAOSDataArrayTemplate<int>> insideOuts;
  if ((insideOuts = vtkAOSDataArrayTemplate<int>::FastDownCast(
    loops_->GetCellData()->GetArray("InsideOuts"))) == nullptr)
  {
    vtkDebugMacro(<< "Loop polygons do not have InsideOuts array. Will resort to "
      << this->GetClassNameInternal() << "::InsideOut = " << this->InsideOut);
    insideOuts = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();
    insideOuts->SetNumberOfComponents(1);
    insideOuts->SetNumberOfTuples(loops->GetNumberOfPolys());
    insideOuts->SetName("InsideOuts");
    insideOuts->FillValue(this->InsideOut);
  }

  vtkSmartPointer<vtkPoints> inPts = input->GetPoints();
  vtkSmartPointer<vtkCellArray> inCells = input->GetPolys();
  vtkSmartPointer<vtkPoints> inLoopPts = loops->GetPoints();
  vtkSmartPointer<vtkCellArray> inLoopPolys = loops->GetPolys();
  vtkNew<vtkPoints> outPts;
  vtkNew<vtkCellArray> outTris, outLines;
  vtkNew<vtkCellData> outTriCd, outLineCd;
  vtkIdType numPts(0), numLoopPts(0), numCells(0), numLoops(0);

  if (!(numPts = inPts->GetNumberOfPoints()) || !(numCells = input->GetNumberOfCells()))
  {
    vtkErrorMacro(<< "Input mesh is empty.");
    return 1;
  }
  if (!(numLoopPts = inLoopPts->GetNumberOfPoints()) || !(numLoops = loops->GetNumberOfPolys()))
  {
    vtkErrorMacro(<< "Input loop is empty.");
    return 1;
  }

  vtkSmartPointer<vtkPointData> inPd = input->GetPointData();
  vtkSmartPointer<vtkPointData> outPd = output->GetPointData();

  vtkSmartPointer<vtkCellData> inCd = input->GetCellData();

  outTris->Allocate(numCells + numLoopPts - 1);
  outLines->Allocate(numLoopPts - 1 + (numCells >> 4));
  outPts->SetDataType(inPts->GetDataType());
  outPts->Allocate(numPts + numLoopPts);
  outPd->InterpolateAllocate(inPd, numPts + numLoopPts * 3);
  outTriCd->CopyAllocate(inCd);
  outLineCd->CopyAllocate(inCd);

  vtkDataArray* inPtsArr = inPts->GetData();
  vtkDataArray* inLoopPtsArr = inLoopPts->GetData();

  double inBounds[6] = {};
  double loopsBnds[6] = {};
  double outBnds[6] = {};
  inPts->GetBounds(inBounds);
  inLoopPts->GetBounds(loopsBnds);

  // extend in z to the combined z-bounds of surface and 2d loops
  {
    const double& zmin = std::min(inBounds[4], loopsBnds[4]);
    const double& zmax = std::max(inBounds[5], loopsBnds[5]);

    std::copy(inBounds, inBounds + 4, outBnds);
    outBnds[4] = zmin ? zmin : -10; // zmin, zmax might be zero.
    outBnds[5] = zmax ? zmax : 10;
  }

  vtkDebugMacro(<< "Output bounds");
  vtkDebugMacro(<< "XMin: " << outBnds[0] << ", XMax: " << outBnds[1]);
  vtkDebugMacro(<< "YMin: " << outBnds[2] << ", YMax: " << outBnds[3]);
  vtkDebugMacro(<< "ZMin: " << outBnds[4] << ", ZMax: " << outBnds[5]);

  this->CreateDefaultLocators();

  this->CellLocator->CacheCellBoundsOn();
  this->CellLocator->SetDataSet(input);
  this->CellLocator->BuildLocator();

  this->PointLocator->SetTolerance(this->Tolerance);
  this->PointLocator->InitPointInsertion(outPts, outBnds);

  // cache loops and cells that might cross.
  std::unordered_map<vtkIdType, SegmentsType> canCross;
  std::vector<std::pair<vtkBoundingBox, vtkNew<vtkIdList>>> loopsInf(numLoops);
  vtkNew<vtkIdList> cells;
  GetEdgeBBoxImpl edgeBBoxWorker;
  canCross.reserve(numCells);
  double edgeBnds[6] = {};

  auto loopsIter = vtk::TakeSmartPointer(inLoopPolys->NewIterator());
  for (loopsIter->GoToFirstCell(); !loopsIter->IsDoneWithTraversal(); loopsIter->GoToNextCell())
  {
    const vtkIdType& loopId = loopsIter->GetCurrentCellId();

    vtkSmartPointer<vtkIdList> loopPtIds = loopsInf[loopId].second;
    loopsIter->GetCurrentCell(loopPtIds);

    const vtkIdType& nedges = loopPtIds->GetNumberOfIds();
    for (vtkIdType lEdge = 0; lEdge < nedges; ++lEdge)
    {
      const vtkIdType& i0 = lEdge;
      const vtkIdType i1 = (lEdge + 1) % nedges ? (lEdge + 1) : 0;
      const vtkIdType& e0 = loopPtIds->GetId(i0);
      const vtkIdType& e1 = loopPtIds->GetId(i1);

      vtkBoundingBox lEdgeBBox;
      if (!dispatchR::Execute(inLoopPtsArr, edgeBBoxWorker, e0, e1, lEdgeBBox))
        edgeBBoxWorker(inLoopPtsArr, e0, e1, lEdgeBBox);

      lEdgeBBox.GetBounds(edgeBnds);

      // extend up to combined {min, max} since that is what CellLocator sees.
      edgeBnds[4] = outBnds[4];
      edgeBnds[5] = outBnds[5];
      this->CellLocator->FindCellsWithinBounds(edgeBnds, cells);

      edgeBnds[4] = 0.0; // back to 2d.
      edgeBnds[5] = 0.0;
      for (const auto& cell : *cells)
      {
        // rule out triangles if the two edge bounding boxes do not intersect
        const vtkIdType* pts = nullptr;
        vtkIdType npts(0);
        inCells->GetCellAtId(cell, npts, pts);
        for (const auto& tEdge : TRIEDGES)
        {
          const vtkIdType vi = pts[tEdge[0]];
          const vtkIdType vj = pts[tEdge[1]];

          vtkBoundingBox tEdgeBBox;
          if (!dispatchR::Execute(inPtsArr, edgeBBoxWorker, vi, vj, tEdgeBBox))
            edgeBBoxWorker(inPtsArr, vi, vj, tEdgeBBox);

          if (tEdgeBBox.IntersectBox(lEdgeBBox))
          {
            if (!canCross[cell].size())
              canCross[cell].reserve(10);

            canCross[cell].emplace_back(e0, e1);
            break;
          }
        }
      }
    }

    double loopBnds[6] = {};
    loops->GetCellBounds(loopId, loopBnds);
    loopBnds[4] = 0.0;
    loopBnds[5] = 0.0;
    loopsInf[loopId].first = vtkBoundingBox(loopBnds);
  }

  vtkNew<vtkUnsignedCharArray> acquisition;
  acquisition->SetName("Acquired");
  acquisition->SetNumberOfComponents(1);
  acquisition->Allocate(numPts + numLoopPts);

  // roughly, a quarter no. of cells
  int reportEvery = (numCells >= 4) ? numCells >> 2 : (numCells >= 2 ? numCells >> 1 : numCells);

  // strategy: 1. push triangle into helper,
  //           2. try to intersect with loops,
  //           3. pop child triangles into output polydata
  //           4. reset helper's data structures to prepare for next triangle.
  SurfCutHelper helper(
    inPd, outPd, outTris, outLines, inCd, outTriCd, outLineCd, this->PointLocator);
  auto cellsIter = vtk::TakeSmartPointer(inCells->NewIterator());
  for (cellsIter->GoToFirstCell(); !cellsIter->IsDoneWithTraversal(); cellsIter->GoToNextCell())
  {
    const vtkIdType& cellId = cellsIter->GetCurrentCellId();

    const vtkIdType* pts = nullptr;
    vtkIdType npts(0);
    cellsIter->GetCurrentCell(npts, pts);
    if (npts != 3)
      vtkWarningMacro(<< "Cannot operate on cell with " << npts << " points"
        << ". Please use vtkTriangleFilter first.");

    // provide
    helper.push(inPts, pts, pts + 1, pts + 2, cellId);

    // embed
    const auto& crossEdges = canCross.find(cellId);
    if (crossEdges != canCross.end())
    {
      const SegmentsType& loopEdges = crossEdges->second;
      helper.triIntersect(inLoopPts, loopEdges, this->Tolerance);
      helper.triangulate(this->Tolerance);
    }
    helper.update();

    // accept/reject
    helper.pop(inLoopPts, insideOuts, loopsInf, acquisition);
    helper.reset();

    if (!(cellId % reportEvery))
      this->UpdateProgress(static_cast<double>(cellId) / numCells);
  }

  // finalize
  outPts->Squeeze();
  output->SetPoints(outPts);
  outTris->Squeeze();
  output->SetPolys(outTris);
  outLines->Squeeze();
  output->SetLines(outLines);
  outPd->Squeeze();
  if (this->ColorAcquiredPts)
  {
    acquisition->Squeeze();
    outPd->AddArray(acquisition);
  }

  const vtkIdType& numLines = output->GetNumberOfLines();
  const vtkIdType& numTris = output->GetNumberOfPolys();
  const vtkIdType& numOutCells = numLines + numTris;

  vtkSmartPointer<vtkCellData> outCd = output->GetCellData();
  outLineCd->Squeeze();
  outTriCd->Squeeze();
  outCd->CopyAllocate(outLineCd, numOutCells);
  outCd->CopyData(outLineCd, 0, numLines, 0);
  outCd->CopyData(outTriCd, numLines, numTris, 0);
  outCd->Squeeze();

  if (this->ColorLoopEdges)
  {
    vtkNew<vtkUnsignedCharArray> constrained;
    constrained->SetName("Constrained");
    constrained->SetNumberOfComponents(1);
    constrained->SetNumberOfTuples(numOutCells);

    for (vtkIdType tupIdx = 0; tupIdx < numLines; ++tupIdx)
      constrained->SetTypedComponent(tupIdx, 0, '\001');

    for (vtkIdType tupIdx = numLines; tupIdx < numOutCells; ++tupIdx)
      constrained->SetTypedComponent(tupIdx, 0, '\000');

    outCd->SetScalars(constrained);
  }

  return 1;
}

void SurfaceCutter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "ColorAcquiredPts: " << (this->ColorAcquiredPts ? "True" : "False") << "\n";
  os << indent << "ColorLoopEdges  : " << (this->ColorLoopEdges ? "True" : "False") << "\n";
  os << indent << "InsideOut: " << (this->InsideOut ? "True" : "False") << "\n";
  os << indent << "Tolerance: " << this->Tolerance << "\n";
  os << indent << "Tolerance2: " << this->Tolerance * this->Tolerance << "\n";

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

void SurfaceCutter::CreateDefaultLocators()
{
  if (!this->CellLocator)
  {
    this->CellLocator = vtkSmartPointer<vtkCellLocator>::New();
  }
  if (!this->PointLocator)
  {
    this->PointLocator = vtkSmartPointer<vtkMergePoints>::New();
  }
}
