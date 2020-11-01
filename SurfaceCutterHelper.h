#ifndef SurfaceCutterHelper_h__
#define SurfaceCutterHelper_h__

#include <algorithm>
#include <array>
#include <vector>
#include <limits>
#include <map>
#include <numeric>
#include <utility>

#include <vtkAOSDataArrayTemplate.h>
#include <vtkArrayDispatch.h>
#include <vtkBoundingBox.h>
#include <vtkCellData.h>
#include <vtkCellIterator.h>
#include <vtkCleanPolyData.h>
#include <vtkDataArrayAccessor.h>
#include <vtkDelaunay2D.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>

namespace compgeom
{
  /*
   * @brief Point-in-polygon. Handle degenerate cases too. https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
   * @param nvert Number of vertices in the polygon.
   * @param vertx Array containing the x-coordinates of the polygon's vertices.
   * @param verty Array containing the y-coordinates of the polygon's vertices.
   * @param testx x-coordinate of the test point.
   * @param testy y-coordinate of the test point.
   * @return true: point in polygon. false: point outside polygon.
   * Note: Conditionals within are quite expensive. Rule out trivial cases with bounding box checks before call.
   */
  template<typename T1, typename T2, typename T3>
  bool pnpoly(const T1& nvert, const T2* vertx, const T2* verty, const T3& testx, const T3& testy)
  {
    T1 i, j;
    bool c = false;
    for (i = 0, j = nvert - 1; i < nvert; j = i++) {
      if (((verty[i] > testy) != (verty[j] > testy)) &&
        (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]))
        c = !c;
    }
    return c;
  }

  template<typename T>
  inline std::array<T, 3> triCentroid(const std::array<T, 3>& p1, const std::array<T, 3>& p2, const std::array<T, 3>& p3) {
    return std::array<T, 3>{ {
        (p1[0] + p2[0] + p3[0]) / T(3),
          (p1[1] + p2[1] + p3[1]) / T(3),
          (p1[2] + p2[2] + p3[2]) / T(3)} };
  }

  /* >0 for P2 left of the line through P0 to P1
   * =0 for P2 on the line
   * <0 for P2 right of the line
   */
  template<typename T>
  inline T isLeft(const T& x0, const T& y0, const T& x1, const T& y1, const T& x2, const T& y2)
  {
    return (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
  }

  template<typename T>
  std::vector<std::size_t> sortIndices(const std::vector<T>& vec)
  {
    // initial indices.
    std::vector<std::size_t> idx(vec.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort with custom comparator.
    std::stable_sort(idx.begin(), idx.end(),
      [&vec](std::size_t idx1, std::size_t idx2)
      {
        return vec[idx1] < vec[idx2];
      }
    );

    return idx;
  }

  template<typename T>
  std::vector<std::array<T, 3>> clockWise(const std::vector<std::array<T, 3>>& points, const std::array<T, 3>& center, std::vector<std::size_t>& sortedIndices)
  {
    std::vector<T> angles;
    std::vector<std::array<T, 3>> cw;

    sortedIndices.clear();
    sortedIndices.reserve(points.size());
    angles.reserve(points.size());
    cw.reserve(points.size());
    for (const auto& pt : points)
      angles.emplace_back(std::atan2(pt[0] - center[0], pt[1] - center[1]));

    for (const auto& idx : sortIndices(angles))
    {
      sortedIndices.emplace_back(idx);
      cw.emplace_back(points[idx]);
    }
    return cw;
  }

  template<typename T>
  std::vector<std::array<T, 3>> clockWise(const std::vector<std::array<T, 3>>& points, const std::array<T, 3>& center)
  {
    std::vector<std::size_t> sortedIndices;
    return clockWise(points, center, sortedIndices);
  }

  template<typename T>
  std::vector<std::array<T, 3>> counterClockWise(const std::vector<std::array<T, 3>>& points, const std::array<T, 3>& center, std::vector<std::size_t>& sortedIndices)
  {
    std::vector<std::array<T, 3>> ccw;
    ccw = clockWise(points, center, sortedIndices);
    std::reverse(ccw.begin(), ccw.end());
    std::reverse(sortedIndices.begin(), sortedIndices.end());
    return ccw;
  }

  template<typename T>
  std::vector<std::array<T, 3>> counterClockWise(const std::vector<std::array<T, 3>>& points, const std::array<T, 3>& center)
  {
    std::vector<std::size_t> sortedIndices;
    return counterClockWise(points, center, sortedIndices);
  }

  enum class JunctionType {
    CROSS,
    COINCIDENT,
    INTERSECT,
    PARALLEL,
    SKEW
  };

  template<typename T>
  const JunctionType intersect(const std::array<T, 3>& p1, const std::array<T, 3>& p2, const std::array<T, 3>& p3, const std::array<T, 3>& p4, std::array<T, 3>& intersectPt)
  {
    //////////////////////////////////////////////////////////////////////////
    // p1       p4
    //    \    /
    //     \  /
    //      \/  intersectPt.
    //      /\
    //     /  \
    //   p3    p2
    //
    // for r, s in [0,1], the lines intersect if:
    // p1 + r(p2 - p1) == p3 + s(p4 - p3).
    //
    // When evaluated for x, y (2-dimensions), the result is two equations in two unknowns (r, s)

    // x1 + r * (x2 - x1) = x3 + s * (x4 - x3)
    // y1 + r * (y2 - y1) = y3 + s * (y4 - y3)
    //
    // r * (x2 - x1) + s * (x3 - x4) = x3 - x1
    // r * (y2 - y1) + s * (y3 - y4) = y3 - y1
    //
    // In matrix form ..
    //  _                 _      _ _       _       _
    // |                   |    |   |     |         |
    // | x2 - x1   x3 - x4 |    | r |     | x3 - x1 |
    // |                   | *  |   |  =  |         |
    // | y2 - y1   y3 - y4 |    | s |     | y3 - y1 |
    // |_                 _|    |_ _|     |_       _|
    //
    // Apply Cramer's rule ..
    //      (x4 - x3) * (y1 - y3) - (x1 - x3) * (y4 - y3) 
    // r = ------------------------------------------------
    //      (x2 - x1) * (y4 - y3) - (x4 - x3) * (y2 - y1)
    //
    //      (x2 - x1) * (y1 - y3) - (x1 - x3) * (y2 - y1)
    // s = ------------------------------------------------
    //      (x2 - x1) * (y4 - y3) - (x4 - x3) * (y2 - y1)
    //
    //////////////////////////////////////////////////////////////////////////
    T r, s;
    T x1 = p1[0], y1 = p1[1];
    T x2 = p2[0], y2 = p2[1];
    T x3 = p3[0], y3 = p3[1];
    T x4 = p4[0], y4 = p4[1];
    constexpr T eps = std::numeric_limits<T>::epsilon();

    // ends are joint.
    if (std::abs(x2 - x4) < eps && std::abs(y2 - y4) < eps)
    {
      intersectPt = p2;
      return JunctionType::INTERSECT;
    }
    else if (std::abs(x1 - x3) < eps && std::abs(y1 - y3) < eps)
    {
      intersectPt = p1;
      return JunctionType::INTERSECT;
    }
    else if (std::abs(x1 - x4) < eps && std::abs(y1 - y4) < eps)
    {
      intersectPt = p1;
      return JunctionType::INTERSECT;
    }
    else if (std::abs(x2 - x3) < eps && std::abs(y2 - y3) < eps)
    {
      intersectPt = p2;
      return JunctionType::INTERSECT;
    }

    T denom = (x2 - x1) * (y4 - y3) - (x4 - x3) * (y2 - y1);

    // parallel
    if (std::abs(denom) < eps)
    {
      intersectPt[0] = 0;
      intersectPt[1] = 0;
      intersectPt[2] = 0;
      return JunctionType::PARALLEL;
    }

    T numer_r = (x4 - x3) * (y1 - y3) - (x1 - x3) * (y4 - y3);
    T numer_s = (x2 - x1) * (y1 - y3) - (x1 - x3) * (y2 - y1);

    // coincident
    if (std::abs(numer_r) < eps && std::abs(numer_s) < eps && std::abs(denom) < eps)
    {
      intersectPt[0] = (x1 + x2) * 0.5;
      intersectPt[1] = (y1 + y2) * 0.5;
      intersectPt[2] = p1[2];
      return JunctionType::COINCIDENT;
    }

    r = numer_r / denom;
    s = numer_s / denom;

    // intersection along line p1---p2
    bool r_gt_0 = !std::signbit(r) && (std::abs(r) - T(0) > eps); // r is positive and r != -0
    bool r_lt_1 = std::signbit(r - T(1)) && (std::abs(r - T(1)) > eps); // r - 1 is negative and r != 1

    // intersection along line p3---p4
    bool s_gt_0 = !std::signbit(s) && (std::abs(s) - T(0) > eps); // s is positive and s != -0
    bool s_lt_1 = std::signbit(s - T(1)) && ((std::abs(s - T(1)) > eps)); // s - 1 is negative and s != 1

    if ((r_gt_0 && r_lt_1) && (s_gt_0 && s_lt_1)) // 0 < r < 1 && 0 < s < 1 ..
    {
      intersectPt[0] = x1 + r * (x2 - x1);
      intersectPt[1] = y1 + r * (y2 - y1);
      intersectPt[2] = p1[2];
      return JunctionType::CROSS;
    }

    bool r_eq_0 = std::abs(r - T(0)) <= eps;
    bool r_eq_1 = std::abs(r - T(1)) <= eps;
    if (r_eq_0 || r_eq_1) // 0 == r || r == 1 .. (exact)
    {
      intersectPt[0] = x1 + r * (x2 - x1);
      intersectPt[1] = y1 + r * (y2 - y1);
      intersectPt[2] = p1[2];
      return JunctionType::INTERSECT;
    }

    bool s_eq_0 = std::abs(s - T(0)) <= eps;
    bool s_eq_1 = std::abs(s - T(1)) <= eps;
    if (s_eq_0 || s_eq_1) // 0 == s || s == 1 .. (exact)
    {
      intersectPt[0] = x3 + s * (x4 - x3);
      intersectPt[1] = y3 + s * (y4 - y3);
      intersectPt[2] = p1[2];
      return JunctionType::INTERSECT;
    }
    return JunctionType::SKEW;
  }
}

// anon begin
namespace
{
  const std::vector<std::vector<vtkIdType>> TRIEDGESEGS = { {0, 1}, {1, 2}, {2, 0} };
  const std::vector<std::pair<vtkIdType, vtkIdType>> TRIEDGES = { {0, 1}, {1, 2}, {2, 0} };

  template<typename LoopsPointsT>
  struct LoopsMeta {
    std::vector<std::array<LoopsPointsT, 3>> coords;
    std::vector<std::pair<vtkIdType, vtkIdType>> edges;
    std::vector<vtkBoundingBox> edgeBboxes, polyBboxes;
    std::vector<int> insideOut;
    std::vector<std::vector<LoopsPointsT>> xs, ys;
    std::vector<vtkIdType> numLoopPts;

    inline void reserve(const vtkIdType& numCells, const vtkIdType& numPoints)
    {
      coords.reserve(numPoints);
      xs.reserve(numCells);
      ys.reserve(numCells);
      edges.reserve(numPoints + 1);
      edgeBboxes.reserve(numPoints + 1);
      polyBboxes.reserve(numCells);
      insideOut.reserve(numCells);
    }

    inline void clear()
    {
      coords.clear();
      xs.clear();
      ys.clear();
      edges.clear();
      edgeBboxes.clear();
      polyBboxes.clear();
      insideOut.clear();
    }
  };

  template<typename PointsT, typename ScalarsT>
  struct MeshPointMeta {
    MeshPointMeta() {}
    MeshPointMeta(const std::array<PointsT, 3>& coord_, const ScalarsT& scalar_, const int& isAcquired_, const int& isIntersect_)
      : coord(coord_), scalar(scalar_), isAcquired(isAcquired_), isIntersect(isIntersect_) {}
    std::array<PointsT, 3> coord;
    ScalarsT scalar;
    int isAcquired;
    int isIntersect;
  };

  template<typename PointsT, typename ScalarsT>
  struct TriMeta {
    TriMeta()
    {
      discard = false;
    }
    TriMeta(const std::array<std::array<PointsT, 3>, 3>& coords_, const std::array<ScalarsT, 3>& scalars_, const std::array<int, 3>& isAcquired_, const std::array<int, 3>& isIntersect_, const std::array<vtkIdType, 3>& verts_)
      : coords(coords_), scalars(scalars_), isAcquired(isAcquired_), isIntersect(isIntersect_), verts(verts_)
    {
      for (std::size_t iPt = 0; iPt < 3; ++iPt)
      {
        xs[iPt] = coords[iPt][0];
        ys[iPt] = coords[iPt][1];
        double xmin = *std::min_element(xs.begin(), xs.end());
        double xmax = *std::max_element(xs.begin(), xs.end());
        double ymin = *std::min_element(ys.begin(), ys.end());
        double ymax = *std::max_element(ys.begin(), ys.end());
        bbox = vtkBoundingBox(xmin, xmax, ymin, ymax, 0.0, 0.0);
        discard = false;
      }
    }
    std::array<PointsT, 3> xs, ys;
    std::array<std::array<PointsT, 3>, 3> coords;
    std::array<ScalarsT, 3> scalars;
    std::array<int, 3> isAcquired;
    std::array<int, 3> isIntersect;
    std::array<vtkIdType, 3> verts;
    vtkBoundingBox bbox;
    bool discard;
  };

  template<typename PointsT, typename ScalarsT>
  ScalarsT triInterp(std::array<PointsT, 3>& coord, const std::array<std::array<PointsT, 3>, 3>& coords, const std::array<ScalarsT, 3>& scalars)
  {
    std::array<double, 3> baryCoords;
    std::vector<std::array<double, 2>> triCoords2d(3);
    std::array<double, 2> interpAt = { coord[0], coord[1] };
    for (std::size_t triVert = 0; triVert < 3; ++triVert) // get 2d only.
    {
      for (std::size_t iDim = 0; iDim < 2; ++iDim)
        triCoords2d[triVert][iDim] = coords[triVert][iDim];
    }

    vtkTriangle::BarycentricCoords(interpAt.data(), triCoords2d[0].data(), triCoords2d[1].data(), triCoords2d[2].data(), baryCoords.data());

    ScalarsT interpScalar = 0;
    ScalarsT barySum = 0;
    for (std::size_t triVert = 0; triVert < 3; ++triVert)
    {
      interpScalar += baryCoords[triVert] * scalars[triVert];
    }
    coord[2] = PointsT(0);
    for (std::size_t triVert = 0; triVert < 3; ++triVert)
    {
      coord[2] += baryCoords[triVert] * coords[triVert][2];
    }
    return interpScalar;
  }

  template<typename PointsT, typename ScalarsT>
  void extractTriInfo(vtkSmartPointer<vtkDataSet> mesh,
    const std::vector<vtkIdType>& smIds,
    const std::vector<MeshPointMeta<PointsT, ScalarsT>>& meshInfo,
    std::vector<TriMeta<PointsT, ScalarsT>>& trisInfo)
  {
    auto smIter = mesh->NewCellIterator();
    std::array<std::array<PointsT, 3>, 3> coords;
    std::array<ScalarsT, 3> scalars;
    std::array<int, 3> isAcquired, isIntersect;
    std::array<vtkIdType, 3> verts;

    for (smIter->InitTraversal(); !smIter->IsDoneWithTraversal(); smIter->GoToNextCell())
    {
      if (smIter->GetCellType() != VTK_TRIANGLE)
        continue;
      vtkSmartPointer<vtkIdList> ptIds = smIter->GetPointIds();

      for (vtkIdType iPt = 0; iPt < 3; ++iPt)
      {
        verts[iPt] = smIds[ptIds->GetId(iPt)];
        coords[iPt] = meshInfo[verts[iPt]].coord;
        scalars[iPt] = meshInfo[verts[iPt]].scalar;
        isAcquired[iPt] = meshInfo[verts[iPt]].isAcquired;
        isIntersect[iPt] = meshInfo[verts[iPt]].isIntersect;
      }
      trisInfo.emplace_back(coords, scalars, isAcquired, isIntersect, verts);
    }
  }

  template<typename PointsT, typename ScalarsT, typename LoopsPointsT>
  void discard(TriMeta<PointsT, ScalarsT>& triInfo, const LoopsMeta<LoopsPointsT>& loopsInfo)
  {
    ///> Remove triangles
    /// 1. insideOut is set, accept triangles inside a polygon.
    /// 2. insideOut is unset, accept triangles outside a polygon.
    auto triCroid = compgeom::triCentroid(triInfo.coords[0], triInfo.coords[1], triInfo.coords[2]);
    std::size_t iPoly(0);
    for (const auto& polyBbox : loopsInfo.polyBboxes)
    {
      if (!loopsInfo.insideOut[iPoly])
      {
        if (!(polyBbox.Contains(triInfo.bbox) || polyBbox.Intersects(triInfo.bbox)))
        {
          ++iPoly;
          continue;
        }
      }

      bool triInsidePoly = compgeom::pnpoly(loopsInfo.xs[iPoly].size(), loopsInfo.xs[iPoly].data(), loopsInfo.ys[iPoly].data(), triCroid[0], triCroid[1]);

      if (loopsInfo.insideOut[iPoly] && !triInsidePoly)
        triInfo.discard = true;
      else if (!loopsInfo.insideOut[iPoly] && triInsidePoly)
        triInfo.discard = true;

      ++iPoly;
    }
  }

  struct ApplyConstraintImpl
  {
    ApplyConstraintImpl(vtkSmartPointer<vtkPolyData> mesh_, const vtkIdType& p1_, const vtkIdType& p2_)
      : mesh(mesh_), p1(p1_), p2(p2_) {}

    vtkSmartPointer<vtkPolyData> mesh;
    vtkIdType p1, p2;

    template<typename PtsArrayT>
    void operator()(PtsArrayT* points)
    {
      //throw std::logic_error("The method or operation is not implemented.");

      using PtsAccess = vtkDataArrayAccessor<PtsArrayT>;

      using PtsType = PtsAccess::APIType;

      VTK_ASSUME(points->GetNumberOfComponents() == 3);

      PtsAccess pointsAccess(points);

      const vtkIdType& pC = p1;
      const vtkIdType& pD_ = p2;

      std::array<PtsType, 3> pCcoord, pD_coord;
      pointsAccess.Get(pC, pCcoord.data());
      pointsAccess.Get(pD_, pD_coord.data());

      // 1 Get triangles around p1.
      vtkIdType* nbrs = nullptr;
      vtkIdType nNbrs = 0;
      mesh->BuildLinks();
      mesh->GetPointCells(pC, nNbrs, nbrs);
      // 2 Find a triangle st it's edge opposite to pC intersects  pC -- pD
      vtkIdType tri0(-1), tri1(-1), pA(-1), pB(-1), pD(-1);
      bool flipTriFound(false);
      for (vtkIdType iNbr = 0; iNbr < nNbrs; ++iNbr)
      {
        const vtkIdType* points = nullptr;
        vtkIdType npts = 0;
        mesh->GetCellPoints(nbrs[iNbr], npts, points);
        for (const auto& edge : TRIEDGES)
        {
          std::array<PtsType, 3> pA_, pB_;
          pointsAccess.Get(points[edge.first], pA_.data());
          pointsAccess.Get(points[edge.second], pB_.data());
          PtsType pASide = compgeom::isLeft(pCcoord[0], pCcoord[1], pD_coord[0], pD_coord[1], pA_[0], pA_[1]);
          PtsType pBSide = compgeom::isLeft(pCcoord[0], pCcoord[1], pD_coord[0], pD_coord[1], pB_[0], pB_[1]);
          bool sameSide = (pASide >= 0 && pBSide >= 0) || (pASide <= 0 && pBSide <= 0);
          if (!sameSide)
          {
            tri0 = nbrs[iNbr];
            pA = points[edge.first];
            pB = points[edge.second];
            auto cellIds = vtkSmartPointer<vtkIdList>::New();
            mesh->GetCellEdgeNeighbors(tri0, pA, pB, cellIds);
            if (!cellIds->GetNumberOfIds()) // shouldn't happen, but eh just to be safe.
              continue;
            tri1 = cellIds->GetId(0);
            flipTriFound = true;
            break;
          }
        }
        if (flipTriFound)
          break;
      }

      if (flipTriFound)
      {
        const vtkIdType* points = nullptr;
        vtkIdType npts(0);
        mesh->GetCellPoints(tri1, npts, points);
        for (const auto& edge : TRIEDGES)
        {
          if ((points[edge.first] == pA && points[edge.second] == pB) ||
            (points[edge.second] == pA && points[edge.first] == pB))
          {
            vtkIdType oppositeVtx = (edge.second + 1) % 3 ? edge.second + 1 : 0;
            pD = points[oppositeVtx];
            break;
          }
        }
        EdgeFlip(pA, pB, pC, pD, tri0, tri1);
      }
    }

    void EdgeFlip(const vtkIdType& pA, const vtkIdType& pB, const vtkIdType& pC, const vtkIdType& pD, const vtkIdType& triABC, const vtkIdType triABD)
    {
      mesh->RemoveReferenceToCell(pA, triABC);
      mesh->RemoveReferenceToCell(pB, triABD);
      mesh->ResizeCellList(pD, 1);
      mesh->AddReferenceToCell(pD, triABD);
      mesh->ResizeCellList(pC, 1);
      mesh->AddReferenceToCell(pC, triABC);

      const vtkIdType triDCA[3] = { pD, pC, pA };
      mesh->ReplaceCell(triABC, 3, triDCA);

      const vtkIdType triDCB[3] = { pD, pC, pB };
      mesh->ReplaceCell(triABD, 3, triDCB);
    }
  };

  struct SurfCutterImpl
  {
  private:
    std::vector<std::array<vtkIdType, 2>> m_Constraints;
  
  public:
    SurfCutterImpl(vtkSmartPointer<vtkPointSet> inMesh_, vtkSmartPointer<vtkPolyData> inLoops_)
      : inMesh(inMesh_), inLoops(inLoops_)
    {
      scalarsName = "Scalars";
    }

    vtkSmartPointer<vtkPointSet> inMesh;
    vtkSmartPointer<vtkPolyData> outMesh;
    vtkSmartPointer<vtkPolyData> inLoops;
    vtkSmartPointer<vtkPolyData> outLoops;
    std::string scalarsName;

    template<typename PointsT, typename ScalarsT>
    vtkSmartPointer<vtkPolyData> GetSubMesh(const std::vector<MeshPointMeta<PointsT, ScalarsT>>& meshInfo, const std::vector<vtkIdType>& subMeshIds, const std::vector<std::vector<vtkIdType>>& cells)
    {
      auto mesh = vtkSmartPointer<vtkPolyData>::New();

      vtkIdType numPoints = subMeshIds.size();
      if (!numPoints)
        return mesh;

      auto coords = vtkSmartPointer<vtkAOSDataArrayTemplate<PointsT>>::New();
      auto scalars = vtkSmartPointer<vtkAOSDataArrayTemplate<ScalarsT>>::New();
      auto acquisitions = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();
      auto intesections = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();

      coords->SetNumberOfComponents(3);
      scalars->SetNumberOfComponents(1);
      acquisitions->SetNumberOfComponents(1);
      intesections->SetNumberOfComponents(1);

      coords->SetNumberOfTuples(numPoints);
      scalars->SetNumberOfTuples(numPoints);
      acquisitions->SetNumberOfTuples(numPoints);
      intesections->SetNumberOfTuples(numPoints);

      scalars->SetName(scalarsName.c_str());
      acquisitions->SetName("Acquired");
      intesections->SetName("Intersected");

      for (vtkIdType tupleIdx = 0; tupleIdx < numPoints; ++tupleIdx)
      {
        const auto& ptId = subMeshIds[tupleIdx];

        for (int iDim = 0; iDim < 3; ++iDim)
          coords->SetTypedComponent(tupleIdx, iDim, meshInfo[ptId].coord[iDim]);

        scalars->SetValue(tupleIdx, meshInfo[ptId].scalar);
        acquisitions->SetValue(tupleIdx, meshInfo[ptId].isAcquired);
        intesections->SetValue(tupleIdx, meshInfo[ptId].isIntersect);
      }

      auto points = vtkSmartPointer<vtkPoints>::New();
      points->SetData(coords);
      mesh->SetPoints(points);
      mesh->GetPointData()->AddArray(scalars);
      mesh->GetPointData()->AddArray(acquisitions);
      mesh->GetPointData()->AddArray(intesections);
      mesh->Allocate(cells.size());

      for (const auto& cell : cells)
      {
        auto numIds = static_cast<int>(cell.size());
        VTKCellType cellType;
        switch (numIds)
        {
        case 2:
          cellType = VTKCellType::VTK_LINE;
          break;
        case 3:
          cellType = VTKCellType::VTK_TRIANGLE;
          break;
        case 4:
          cellType = VTKCellType::VTK_QUAD;
          break;
        default:
          cellType = VTKCellType::VTK_POLYGON;
          break;
        }
        mesh->InsertNextCell(static_cast<int>(cellType), numIds, cell.data());
      }
      return mesh;
    }

    template<typename PointsT, typename ScalarsT>
    void CreateTriMesh(const std::vector<MeshPointMeta<PointsT, ScalarsT>>& meshInfos, const std::vector<TriMeta<PointsT, ScalarsT>>& tris)
    {
      outMesh = vtkSmartPointer<vtkPolyData>::New();

      vtkIdType numPoints = meshInfos.size();

      auto coords = vtkSmartPointer<vtkAOSDataArrayTemplate<PointsT>>::New();
      auto scalars = vtkSmartPointer<vtkAOSDataArrayTemplate<ScalarsT>>::New();
      auto acquisitions = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();
      auto intesections = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();

      coords->SetNumberOfComponents(3);
      scalars->SetNumberOfComponents(1);
      acquisitions->SetNumberOfComponents(1);
      intesections->SetNumberOfComponents(1);

      coords->SetNumberOfTuples(numPoints);
      scalars->SetNumberOfTuples(numPoints);
      acquisitions->SetNumberOfTuples(numPoints);
      intesections->SetNumberOfTuples(numPoints);

      scalars->SetName(scalarsName.c_str());
      acquisitions->SetName("Acquired");
      intesections->SetName("Intersected");

      for (vtkIdType tupleIdx = 0; tupleIdx < numPoints; ++tupleIdx)
      {
        for (int iDim = 0; iDim < 3; ++iDim)
          coords->SetTypedComponent(tupleIdx, iDim, meshInfos[tupleIdx].coord[iDim]);

        scalars->SetValue(tupleIdx, meshInfos[tupleIdx].scalar);
        acquisitions->SetValue(tupleIdx, meshInfos[tupleIdx].isAcquired);
        intesections->SetValue(tupleIdx, meshInfos[tupleIdx].isIntersect);
      }

      auto points = vtkSmartPointer<vtkPoints>::New();
      points->SetData(coords);
      outMesh->SetPoints(points);
      outMesh->GetPointData()->SetScalars(scalars);
      outMesh->GetPointData()->AddArray(acquisitions);
      outMesh->GetPointData()->AddArray(intesections);
      outMesh->Allocate(tris.size());

      for (const auto& tri : tris)
      {
        if (tri.discard)
          continue;
        outMesh->InsertNextCell(VTK_TRIANGLE, 3, tri.verts.data());
      }
      auto cleanUp = vtkSmartPointer<vtkCleanPolyData>::New();
      cleanUp->SetInputData(outMesh);
      cleanUp->PointMergingOn();
      cleanUp->Update();
      outMesh->ShallowCopy(cleanUp->GetOutput());
    }

    template<typename PointsT, typename ScalarsT>
    void CreateProjectedLoop(const std::vector<MeshPointMeta<PointsT, ScalarsT>>& meshInfos)
    {
      outLoops = vtkSmartPointer<vtkPolyData>::New();

      vtkIdType numPoints = meshInfos.size();

      auto coords = vtkSmartPointer<vtkAOSDataArrayTemplate<PointsT>>::New();
      auto scalars = vtkSmartPointer<vtkAOSDataArrayTemplate<ScalarsT>>::New();
      auto acquisitions = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();
      auto intesections = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();

      coords->SetNumberOfComponents(3);
      scalars->SetNumberOfComponents(1);
      acquisitions->SetNumberOfComponents(1);
      intesections->SetNumberOfComponents(1);

      coords->SetNumberOfTuples(numPoints);
      scalars->SetNumberOfTuples(numPoints);
      acquisitions->SetNumberOfTuples(numPoints);
      intesections->SetNumberOfTuples(numPoints);

      scalars->SetName(scalarsName.c_str());
      acquisitions->SetName("Acquired");
      intesections->SetName("Intersected");

      for (vtkIdType tupleIdx = 0; tupleIdx < numPoints; ++tupleIdx)
      {
        for (int iDim = 0; iDim < 3; ++iDim)
          coords->SetTypedComponent(tupleIdx, iDim, meshInfos[tupleIdx].coord[iDim]);

        scalars->SetValue(tupleIdx, meshInfos[tupleIdx].scalar);
        acquisitions->SetValue(tupleIdx, meshInfos[tupleIdx].isAcquired);
        intesections->SetValue(tupleIdx, meshInfos[tupleIdx].isIntersect);
      }

      auto points = vtkSmartPointer<vtkPoints>::New();
      points->SetData(coords);
      outLoops->SetPoints(points);
      outLoops->GetPointData()->SetScalars(scalars);
      outLoops->GetPointData()->AddArray(acquisitions);
      outLoops->GetPointData()->AddArray(intesections);
      outLoops->Allocate(m_Constraints.size());

      for (const auto& edge : m_Constraints)
      {
        outLoops->InsertNextCell(VTK_LINE, 2, edge.data());
      }
      auto cleanUp = vtkSmartPointer<vtkCleanPolyData>::New();
      cleanUp->SetInputData(outLoops);
      cleanUp->PointMergingOn();
      cleanUp->Update();
      outLoops->ShallowCopy(cleanUp->GetOutput());
    }

    template<typename PointsArrT1, typename ScalarsArrT, typename PointsArrT2>
    void operator()(PointsArrT1* meshPts, ScalarsArrT* scalars, PointsArrT2* loopsPts)
    {
      const vtkIdType& numLoops = inLoops->GetNumberOfCells();
      const vtkIdType& numCells = inMesh->GetNumberOfCells();
      if (!(numLoops && numCells))
        return;

      const vtkIdType& numPoints1 = meshPts->GetNumberOfTuples();
      const vtkIdType& numPoints2 = loopsPts->GetNumberOfTuples();
      const vtkIdType& numScalars = scalars->GetNumberOfTuples();

      if (!numPoints1)
        return;

      if (numPoints1 != numScalars)
        return;

      if (!numPoints2)
        return;

      auto insideOuts = vtkAOSDataArrayTemplate<int>::FastDownCast(inLoops->GetCellData()->GetArray("InsideOuts"));
      if (!insideOuts)
        return;

      scalarsName = scalars->GetName();
      m_Constraints.reserve(numLoops * (numPoints2 + 1)); // no. of intersections, a rough estimate 

      auto acquisitions = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();
      acquisitions->SetNumberOfComponents(1);
      acquisitions->SetNumberOfTuples(numPoints1);
      acquisitions->SetName("Acquired");
      acquisitions->FillValue(0);

      auto intersections = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();
      intersections->SetNumberOfComponents(1);
      intersections->SetNumberOfTuples(numPoints1);
      intersections->SetName("Intersected");
      intersections->FillValue(0);

      using Points1Access = vtkDataArrayAccessor<PointsArrT1>;
      using ScalarsAccess = vtkDataArrayAccessor<ScalarsArrT>;
      using Points2Access = vtkDataArrayAccessor<PointsArrT2>;

      using PointsT1 = Points1Access::APIType;
      using ScalarsT = ScalarsAccess::APIType;
      using PointsT2 = Points2Access::APIType;

      using MMeta = MeshPointMeta<PointsT1, ScalarsT>;
      using TMeta = TriMeta<PointsT1, ScalarsT>;
      using LMeta = LoopsMeta<PointsT2>;

      VTK_ASSUME(meshPts->GetNumberOfComponents() == 3);
      VTK_ASSUME(scalars->GetNumberOfComponents() == 1);
      VTK_ASSUME(loopsPts->GetNumberOfComponents() == 3);

      Points1Access meshAccess(meshPts);
      Points2Access loopsAccess(loopsPts);
      ScalarsAccess scalarsAccess(scalars);

      ///> loops info.
      LMeta loopsInfo;
      loopsInfo.reserve(numLoops, numPoints2);
      for (vtkIdType iPt = 0; iPt < numPoints2; ++iPt)
      {
        std::array<PointsT2, 3> coord;
        loopsAccess.Get(iPt, coord.data());
        loopsInfo.coords.emplace_back(coord);
      }

      vtkSmartPointer<vtkCellIterator> polysIter = inLoops->NewCellIterator();
      vtkIdType iPoly(0);
      for (polysIter->InitTraversal(); !polysIter->IsDoneWithTraversal(); polysIter->GoToNextCell(), ++iPoly)
      {
        vtkSmartPointer<vtkIdList> ptIds = polysIter->GetPointIds();
        const vtkIdType& numIds = ptIds->GetNumberOfIds();
        std::vector<PointsT2> _xs(numIds), _ys(numIds);
        std::vector<std::array<PointsT2, 3>> poly(numIds);

        for (vtkIdType idx = 0; idx < numIds; ++idx)
        {
          // information pertaining to an edge of a loop
          const vtkIdType& p0Id = idx;
          const vtkIdType p1Id = (idx + 1) % numIds ? idx + 1 : 0;
          const vtkIdType& p0 = ptIds->GetId(p0Id);
          const vtkIdType& p1 = ptIds->GetId(p1Id);
          loopsInfo.edges.emplace_back(p0, p1);

          double xmin = std::min(loopsInfo.coords[p0][0], loopsInfo.coords[p1][0]);
          double xmax = std::max(loopsInfo.coords[p0][0], loopsInfo.coords[p1][0]);
          double ymin = std::min(loopsInfo.coords[p0][1], loopsInfo.coords[p1][1]);
          double ymax = std::max(loopsInfo.coords[p0][1], loopsInfo.coords[p1][1]);
          loopsInfo.edgeBboxes.emplace_back(xmin, xmax, ymin, ymax, 0.0, 0.0);

          _xs[idx] = loopsInfo.coords[p0][0];
          _ys[idx] = loopsInfo.coords[p0][1];
        }

        // pertaining to loop polygon
        double xmin = *std::min_element(_xs.begin(), _xs.end());
        double xmax = *std::max_element(_xs.begin(), _xs.end());
        double ymin = *std::min_element(_ys.begin(), _ys.end());
        double ymax = *std::max_element(_ys.begin(), _ys.end());
        loopsInfo.xs.emplace_back(_xs);
        loopsInfo.ys.emplace_back(_ys);
        loopsInfo.insideOut.push_back(insideOuts->GetValue(iPoly));
        loopsInfo.numLoopPts.push_back(numIds);
        loopsInfo.polyBboxes.emplace_back(xmin, xmax, ymin, ymax, 0.0, 0.0);
      }

      // now store loops coordinates as PointsT1.
      std::vector<std::array<PointsT1, 3>> loopsCoords;
      loopsCoords.reserve(numPoints2);
      for (const auto& coord : loopsInfo.coords)
      {
        std::array<PointsT1, 3> newCoord;
        std::transform(coord.begin(), coord.end(), newCoord.begin(), [&](const auto& val) {return static_cast<PointsT1>(val); });
        loopsCoords.emplace_back(newCoord);
      }

      // cache point data in meshInfo, tri metadata in trisInfo
      std::vector<MMeta> meshInfo(numPoints1);
      std::vector<TMeta> trisInfo;
      std::vector<vtkIdType> triPtIds(numPoints1);

      for (vtkIdType tupleIdx = 0; tupleIdx < numPoints1; ++tupleIdx)
      {
        meshAccess.Get(tupleIdx, meshInfo[tupleIdx].coord.data());
        meshInfo[tupleIdx].scalar = scalarsAccess.Get(tupleIdx, 0);
        meshInfo[tupleIdx].isAcquired = acquisitions->GetValue(tupleIdx);
        meshInfo[tupleIdx].isIntersect = intersections->GetValue(tupleIdx);
        triPtIds[tupleIdx] = tupleIdx;
      }
      trisInfo.reserve(numCells);
      extractTriInfo(inMesh, triPtIds, meshInfo, trisInfo);
      // 3 new triangles per loop's point, 3 new intersection points per triangle, 3 new intersection points for edge.
      meshInfo.reserve(meshInfo.size() + numPoints2 * 3 * 3 * 3);
      // same as meshInfo, further 3 sub triangles per each new triangle.
      trisInfo.reserve(trisInfo.size() + numPoints2 * 3 * 3 * 3 * 3);

      // for each tri.
      // interp.
      // build submesh
      // add new tris.
      // for each tri.
      // cut
      // build submesh
      // add new tris.
      // remove inside/outside.

      ///> A The sub-mesh after interpolation might consist of 
      /// 1. the actual triangle.
      /// 2. all loops' points that are enclosed by that triangle.
      const std::size_t& nTris1 = trisInfo.size();
      std::vector<vtkIdType> smIds;
      smIds.reserve(numPoints2);

      for (std::size_t iTri = 0; iTri < nTris1; ++iTri)
      {
        TMeta& triInfo = trisInfo[iTri];
        smIds.clear();
        for (vtkIdType iPt = 0; iPt < 3; ++iPt)
        {
          smIds.push_back(triInfo.verts[iPt]);
        }

        for (auto& coord : loopsCoords)
        {
          if (!triInfo.bbox.ContainsPoint(coord[0], coord[1], 0.0))
            continue;

          if (!compgeom::pnpoly(3, triInfo.xs.data(), triInfo.ys.data(), coord[0], coord[1]))
          {
            // check corners.
            for (std::size_t triVert = 0; triVert < 3; ++triVert)
            {
              bool xSimilar = std::abs(triInfo.xs[triVert] - coord[0]) <= std::numeric_limits<float>::epsilon();
              bool ySimilar = std::abs(triInfo.ys[triVert] - coord[1]) <= std::numeric_limits<float>::epsilon();
              if (xSimilar && ySimilar)
              {
                acquisitions->SetValue(triInfo.verts[triVert], 1);
                intersections->SetValue(triInfo.verts[triVert], 1);
              }
            }
            continue; // move on to next point.
          }

          ///> Fastest and sensible interpolation using bary centric coords.
          ScalarsT newScalar = triInterp<PointsT1, ScalarsT>(coord, triInfo.coords, triInfo.scalars);

          smIds.emplace_back(meshInfo.size());
          meshInfo.emplace_back(coord, newScalar, 1, 1);
        } // end interpolation loop

        if (smIds.size() > 3)
        {
          triInfo.discard = true;
          // setup constraints.
          vtkSmartPointer<vtkPolyData> subMesh = GetSubMesh(meshInfo, smIds, TRIEDGESEGS);
          // triangulate with constraints.
          auto del2d = vtkSmartPointer<vtkDelaunay2D>::New();
          del2d->SetInputData(subMesh);
          del2d->SetSourceData(subMesh);
          del2d->Update();
          subMesh->ShallowCopy(del2d->GetOutput());
          extractTriInfo(subMesh, smIds, meshInfo, trisInfo);
        }
        else
        {
          continue;
        }
      }

      ///> B Now insert loops's edges into triangles. 
      /// 1. Throw loops edges onto triangles and see if any of either intersect.
      /// 2. Triangulate the actual triangle + intersection points.
      /// 3. Identify pairs of intersections as constraints.
      /// 4. Apply constraints on the sub mesh.
      const std::size_t& nTris2 = trisInfo.size();
      //smIds.swap(std::vector<std::size_t>());
      smIds.reserve(numLoops * (numPoints2 + 1)); // no. of intersections, a rough estimate 
      for (std::size_t iTri = 0; iTri < nTris2; ++iTri)
      {
        TMeta& triInfo = trisInfo[iTri];
        if (triInfo.discard)
          continue;

        smIds.clear();

        for (vtkIdType iPt = 0; iPt < 3; ++iPt)
          smIds.push_back(triInfo.verts[iPt]);

        vtkIdType iLoopEdge(-1);
        std::map<vtkIdType, vtkIdType> constraints;
        std::vector<vtkIdType> loopHits;
        std::vector<vtkIdType> triHitAt;
        bool triCrossesLoop(false);
        for (const auto& loopsEdge : loopsInfo.edges) // for each loops' edge
        {
          loopHits.clear();
          triHitAt.clear();
          ++iLoopEdge;
          const auto& pa_0 = loopsCoords[loopsEdge.first];
          const auto& pb_0 = loopsCoords[loopsEdge.second];
          const auto& bbox_0 = loopsInfo.edgeBboxes[iLoopEdge];

          if (!bbox_0.Intersects(triInfo.bbox))
            continue;

          loopHits.reserve(3);
          triHitAt.reserve(3);
          for (vtkIdType iTriEdge = 0; iTriEdge < 3; ++iTriEdge) // for each triangle's edge
          {
            const std::pair<vtkIdType, vtkIdType>& triEdge = TRIEDGES[iTriEdge];
            const auto& pa_1 = triInfo.coords[triEdge.first];
            const auto& pb_1 = triInfo.coords[triEdge.second];
            std::array<PointsT1, 3> pc;

            auto junctionType = compgeom::intersect(pa_0, pb_0, pa_1, pb_1, pc);
            if (junctionType == compgeom::JunctionType::CROSS)
            {
              if (pc[0] == pa_0[0] && pc[1] == pa_0[1])
                continue;

              if (pc[0] == pb_0[0] && pc[1] == pb_0[1])
                continue;

              auto newScalar = triInterp<PointsT1, ScalarsT>(pc, triInfo.coords, triInfo.scalars);

              smIds.emplace_back(meshInfo.size());
              meshInfo.emplace_back(pc, newScalar, 0, 1);

              loopHits.emplace_back(static_cast<vtkIdType>(smIds.size() - 1));
              triHitAt.emplace_back(iTriEdge);
              triCrossesLoop |= true;
            }
            else
            {
              triCrossesLoop |= false;
            }
          } // end for triangle's edge

          if (loopHits.size() == 2)
          {
            constraints[loopHits[0]] = loopHits[1];
          }
          else if (loopHits.size() == 1)
          {
            const vtkIdType crossedTriEdge = triHitAt[0];
            const vtkIdType oppositeVtx = (TRIEDGES[crossedTriEdge].second + 1) % 3 ? TRIEDGES[crossedTriEdge].second + 1 : 0;
            constraints[loopHits[0]] = oppositeVtx;
          }

        } // end for loops' edge
        if (!triCrossesLoop)
        {
          discard(triInfo, loopsInfo);
          continue;
        }
        triInfo.discard = true;
        std::vector<std::array<PointsT1, 3>> smCoords;
        for (const auto& smId : smIds)
        {
          smCoords.emplace_back(meshInfo[smId].coord);
        }
        auto centroid = compgeom::triCentroid(triInfo.coords[0], triInfo.coords[1], triInfo.coords[2]);
        std::vector<std::size_t> ccwIds;
        compgeom::counterClockWise(smCoords, centroid, ccwIds);
        std::vector<std::vector<vtkIdType>> cells;
        cells.emplace_back(std::vector<vtkIdType>(ccwIds.begin(), ccwIds.end()));

        // now triangulate (ear-cut) and apply loop edge constraints
        vtkSmartPointer<vtkPolyData> subMesh = GetSubMesh(meshInfo, smIds, cells);
        auto newPolys = vtkSmartPointer<vtkCellArray>::New();
        auto poly = vtkPolygon::New();
        auto ptIds = vtkIdList::New();
        // initialize polygon
        const auto& pts = cells[0];
        const auto& numIds = pts.size();
        poly->PointIds->SetNumberOfIds(numIds);
        poly->Points->SetNumberOfPoints(numIds);
        for (int i = 0; i < numIds; i++)
        {
          poly->PointIds->SetId(i, pts[i]);
          poly->Points->SetPoint(i, subMesh->GetPoint(pts[i]));
        }
        poly->Triangulate(ptIds);
        const auto& numPts = ptIds->GetNumberOfIds();
        const auto numSimplices = numPts / 3;
        std::array<vtkIdType, 3> triPts;
        for (vtkIdType i = 0; i < numSimplices; i++)
        {
          for (vtkIdType j = 0; j < 3; j++)
          {
            triPts[j] = poly->PointIds->GetId(ptIds->GetId(3 * i + j));
          }
          newPolys->InsertNextCell(3, triPts.data());
        } // for each simplex
        poly->Delete();
        ptIds->Delete();
        subMesh->SetPolys(newPolys);
        subMesh->BuildLinks();
        for (const auto& constraint : constraints)
        {
          unsigned short niters = 0;
          while (!subMesh->IsEdge(constraint.first, constraint.second))
          {
            using dispatcher = vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals>;
            auto constraintWorker = ApplyConstraintImpl(subMesh, constraint.first, constraint.second);
            vtkDataArray* points = subMesh->GetPoints()->GetData();
            if (!dispatcher::Execute(points, constraintWorker))
            {
              constraintWorker(points);
            }
            ++niters;
            if (niters > 32)
            {
              //__debugbreak();
              break;
            }
          }
          m_Constraints.emplace_back(std::array<vtkIdType, 2>({ smIds[constraint.first], smIds[constraint.second] }));
        }

        const auto& oldSz = trisInfo.size();
        extractTriInfo(subMesh, smIds, meshInfo, trisInfo);
        const auto& newSz = trisInfo.size();

        for (std::size_t iInfo = oldSz; iInfo < newSz; ++iInfo)
        {
          discard(trisInfo[iInfo], loopsInfo);
        }
      }

      // collect all triangles and finish
      CreateTriMesh(meshInfo, trisInfo);
      CreateProjectedLoop(meshInfo);
    }
  };
}
// anon end

#endif // SurfaceCutterHelper_h__