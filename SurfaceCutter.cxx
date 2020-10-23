#include "SurfaceCutter.h"

#include <vtkAOSDataArrayTemplate.h>
#include <vtkArrayDispatch.h>
#include <vtkCellData.h>
#include <vtkCellIterator.h>
#include <vtkDataArrayAccessor.h>
#include <vtkDataSet.h>
#include <vtkDelaunay2D.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <numeric>
#include <set>
#include <utility>

vtkStandardNewMacro(SurfaceCutter);

SurfaceCutter::SurfaceCutter()
{
  this->ComputeBoolean2D = true;
  this->InsideOut = false; // remove portion inside polygons.
  this->TagInsertedPoints = true;

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);

  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::SCALARS);

  vtkDebugMacro(<< "Initialized " << this->GetClassNameInternal());
}

SurfaceCutter::~SurfaceCutter()
{
  vtkDebugMacro(<< "Destroyed " << this->GetClassNameInternal());
}

int SurfaceCutter::FillInputPortInformation(int port, vtkInformation* info)
{
  switch (port)
  {
  case 0:
    break;
  case 1:
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    break;
  default:
    break;
  }
  return 1;
}

void SurfaceCutter::SetLoops(vtkDataSet* loops)
{
  this->SetInputData(1, loops);
}

// anon begin
namespace
{ 
  const std::vector<std::vector<vtkIdType>> TRIEDGESEGS = { {0, 1}, {1, 2}, {2, 0} };
  const std::vector<std::vector<vtkIdType>> TRIPTIDS = { {0, 1, 2} };
  const std::vector<std::pair<vtkIdType, vtkIdType>> TRIEDGES = { {0, 1}, {1, 2}, {2, 0} };

  template<typename T>
  struct BBox {
    BBox() {}
    BBox(const std::array<T, 3>& p0, const std::array<T, 3>& p1) {
      left = std::min(p0[0], p1[0]);
      right = std::max(p0[0], p1[0]);
      bottom = std::min(p0[1], p1[1]);
      top = std::max(p0[1], p1[1]);
    }
    BBox(const std::vector<T>& xs, const std::vector<T>& ys)
    {
      left = *std::min_element(xs.begin(), xs.end());
      right = *std::max_element(xs.begin(), xs.end());
      bottom = *std::min_element(ys.begin(), ys.end());
      top = *std::max_element(ys.begin(), ys.end());
    }
    BBox(const std::vector<std::array<T, 3>>& points)
    {
      std::vector<T> xs, ys;
      xs.reserve(points.size()); ys.reserve(points.size());
      for (const auto& point : points)
      {
        xs.emplace_back(point[0]); ys.emplace_back(point[1]);
      }
      left = *std::min_element(xs.begin(), xs.end());
      right = *std::max_element(xs.begin(), xs.end());
      bottom = *std::min_element(ys.begin(), ys.end());
      top = *std::max_element(ys.begin(), ys.end());
    }
    BBox(const std::array<T, 4>& bbox)
    {
      left = bbox[0]; right = bbox[1]; bottom = bbox[2]; top = bbox[3];
    }
    BBox(const std::array<T, 6>& bbox)
    {
      left = bbox[0]; right = bbox[1]; bottom = bbox[2]; top = bbox[3];
    }
    T left;
    T right;
    T bottom;
    T top;
    inline T width() const noexcept { return right - left; }
    inline T height() const noexcept { return top - bottom; }
    inline T diagonal() const { return std::hypot(right - left, top - bottom); }
  };

  template<typename T>
  inline static bool fullyInsideBbox(const BBox<T>& bbox1, const BBox<T>& bbox2) noexcept {
    ///> bbox1 inside bbox2
    if (bbox1.left < bbox2.left && bbox2.right < bbox1.right && bbox1.bottom < bbox2.bottom && bbox2.top < bbox1.top)
    {
      return true;
    }
    ///> bbox2 inside bbox1
    if (bbox2.left < bbox1.left && bbox1.right < bbox2.right && bbox2.bottom < bbox1.bottom && bbox1.top < bbox2.top)
    {
      return true;
    }
    return false;
  }

  template<typename T1, typename T2>
  inline static bool fullyOutsideBbox(const BBox<T1>& bbox1, const BBox<T2>& bbox2, bool yDownward = false) noexcept {
    if (!yDownward)
    {
      return (bbox2.left > bbox1.right || bbox2.right < bbox1.left || bbox2.top < bbox1.bottom || bbox2.bottom > bbox1.top);
    }
    else
    {
      return (bbox2.left > bbox1.right || bbox2.right < bbox1.left || bbox2.top > bbox1.bottom || bbox2.bottom < bbox1.top);
    }
  }

  template<typename T1, typename T2>
  inline static bool intersectBbox(const BBox<T1>& bbox1, const BBox<T2>& bbox2, bool yDownward = false) noexcept {
    if (!yDownward)
    {
      return (bbox1.left <= bbox2.right && bbox1.right >= bbox2.left && bbox1.bottom <= bbox2.top && bbox1.top > bbox2.bottom);
    }
    else
    {
      return (bbox1.left <= bbox2.right && bbox1.right >= bbox2.left && bbox1.top <= bbox2.bottom && bbox1.bottom > bbox2.top);
    }
  }
  
  template<typename T>
  inline static bool pointInBbox(const T& x, const T& y, const BBox<T>& bbox, bool yDownward = false) noexcept {
    if (!yDownward)
      return (bbox.left < x) && (x < bbox.right) && (bbox.bottom < y) && (y < bbox.top);
    else
      return (bbox.left < x) && (x < bbox.right) && (bbox.top < y) && (y < bbox.bottom);
  }

  /**
     * @brief Point-in-polygon. Handle degenerate cases too. https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
     * @param nvert Number of vertices in the polygon.
     * @param vertx Array containing the x-coordinates of the polygon's vertices.
     * @param verty Array containing the y-coordinates of the polygon's vertices.
     * @param testx x-coordinate of the test point.
     * @param testy y-coordinate of the test point.
     * @return true: point in polygon. false: point outside polygon.
     */
  template<typename T1, typename T2, typename T3>
  static bool pnpoly(const T1& nvert, const T2* vertx, const T2* verty, const T3& testx, const T3& testy)
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

  /** >0 for P2 left of the line through P0 to P1
   *  =0 for P2 on the line
   *  <0 for P2 right of the line
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
  static std::vector<std::array<T, 3>> clockWise(const std::vector<std::array<T, 3>>& points, const std::array<T, 3>& center, std::vector<std::size_t>& sortedIndices)
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
  static std::vector<std::array<T, 3>> clockWise(const std::vector<std::array<T, 3>>& points, const std::array<T, 3>& center)
  {
    std::vector<std::size_t> sortedIndices;
    return clockWise(points, center, sortedIndices);
  }

  template<typename T>
  static std::vector<std::array<T, 3>> counterClockWise(const std::vector<std::array<T, 3>>& points, const std::array<T, 3>& center, std::vector<std::size_t>& sortedIndices)
  {
    std::vector<std::array<T, 3>> ccw;
    ccw = clockWise(points, center, sortedIndices);
    std::reverse(ccw.begin(), ccw.end());
    std::reverse(sortedIndices.begin(), sortedIndices.end());
    return ccw;
  }

  template<typename T>
  static std::vector<std::array<T, 3>> counterClockWise(const std::vector<std::array<T, 3>>& points, const std::array<T, 3>& center)
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

  /**
    * @brief 2d intersection of line segments p1---p2, p3---p4. z-coordinate of p1 is copied over.
    * @tparam T data type of points.
    * @param p1 x,y,z coordinates.
    * @param p2 x,y,z coordinates.
    * @param p3 x,y,z coordinates.
    * @param p4 x,y,z coordinates.
    * @param intersectPt x,y,z coordinates.
    * @return true: intersection, false: no intersection.
  */
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

  template<typename ClipPtsType>
  struct ClipperMeta {
    std::vector<std::array<ClipPtsType, 3>> coords;
    std::vector<std::pair<vtkIdType, vtkIdType>> edges;
    std::vector<BBox<ClipPtsType>> edgeBboxes, polyBboxes;
    std::vector<int> invertedOp;
    std::vector<std::vector<ClipPtsType>> xs, ys;

    inline void reserve(const vtkIdType& numCells, const vtkIdType& numPoints)
    {
      coords.reserve(numPoints);
      xs.reserve(numCells);
      ys.reserve(numCells);
      edges.reserve(numPoints + 1);
      edgeBboxes.reserve(numPoints + 1);
      polyBboxes.reserve(numCells);
      invertedOp.reserve(numCells);
    }

    inline void clear()
    {
      coords.clear();
      xs.clear();
      ys.clear();
      edges.clear();
      edgeBboxes.clear();
      polyBboxes.clear();
      invertedOp.clear();
    }
  };

  template<typename MeshPtsType, typename ScalarsType>
  struct MeshMeta {
    MeshMeta() {}
    MeshMeta(const std::array<MeshPtsType, 3>& coord, const ScalarsType& scalar, const int& pCrit, const vtkIdType& originalIdx)
      : coord(coord), scalar(scalar), pCrit(pCrit), originalIdx(originalIdx) {}
    std::array<MeshPtsType, 3> coord;
    ScalarsType scalar;
    int pCrit;
    vtkIdType originalIdx;
  };

  template<typename MeshPtsType, typename ScalarsType>
  struct TriMeta {
    std::array<MeshPtsType, 3> xs, ys;
    std::array<std::array<MeshPtsType, 3>, 3> coords;
    std::array<ScalarsType, 3> scalars;
    std::array<int, 3> pCrits;
    std::array<vtkIdType, 3> verts;
    BBox<MeshPtsType> bbox;
    bool discard = false;
  };

  template<typename MeshPtsType, typename ScalarsType>
  ScalarsType TriInterp(std::array<MeshPtsType, 3>& coord, const std::array<std::array<MeshPtsType, 3>, 3>& coords, const std::array<ScalarsType, 3>& scalars)
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

    ScalarsType interpScalar = 0;
    ScalarsType barySum = 0;
    for (std::size_t triVert = 0; triVert < 3; ++triVert)
    {
      interpScalar += baryCoords[triVert] * scalars[triVert];
    }
    coord[2] = MeshPtsType(0);
    for (std::size_t triVert = 0; triVert < 3; ++triVert)
    {
      coord[2] += baryCoords[triVert] * coords[triVert][2];
    }
    return interpScalar;
  }
  
  template<typename MeshPtsType, typename ScalarsType>
  void ExtractTris(vtkSmartPointer<vtkDataSet> mesh,
    std::vector<std::array<vtkIdType, 3>>& tris,
    std::vector<MeshMeta<MeshPtsType, ScalarsType>>& smInfo,
    std::vector<MeshMeta<MeshPtsType, ScalarsType>>& meshInfo,
    std::vector<TriMeta<MeshPtsType, ScalarsType>>& trisInfo)
  {
    auto smIter = mesh->NewCellIterator();
    for (smIter->InitTraversal(); !smIter->IsDoneWithTraversal(); smIter->GoToNextCell())
    {
      if (smIter->GetCellType() != VTK_TRIANGLE)
        continue;
      vtkSmartPointer<vtkIdList> ptIds = smIter->GetPointIds();

      const vtkIdType& p0 = ptIds->GetId(0);
      const vtkIdType& p1 = ptIds->GetId(1);
      const vtkIdType& p2 = ptIds->GetId(2);
      tris.emplace_back(std::array<vtkIdType, 3>({ smInfo[p0].originalIdx, smInfo[p1].originalIdx, smInfo[p2].originalIdx }));

      std::array<vtkIdType, 3>& tri = tris.back();

      TriMeta<MeshPtsType, ScalarsType> newTriInfo;
      for (vtkIdType iPt = 0; iPt < 3; ++iPt)
      {
        newTriInfo.verts[iPt] = tri[iPt];
        newTriInfo.coords[iPt] = meshInfo[newTriInfo.verts[iPt]].coord;
        newTriInfo.xs[iPt] = newTriInfo.coords[iPt][0];
        newTriInfo.ys[iPt] = newTriInfo.coords[iPt][1];

        newTriInfo.scalars[iPt] = meshInfo[newTriInfo.verts[iPt]].scalar;
        newTriInfo.pCrits[iPt] = meshInfo[newTriInfo.verts[iPt]].pCrit;
      }
      newTriInfo.bbox.left = *std::min_element(newTriInfo.xs.begin(), newTriInfo.xs.end());
      newTriInfo.bbox.right = *std::max_element(newTriInfo.xs.begin(), newTriInfo.xs.end());
      newTriInfo.bbox.bottom = *std::min_element(newTriInfo.ys.begin(), newTriInfo.ys.end());
      newTriInfo.bbox.top = *std::max_element(newTriInfo.ys.begin(), newTriInfo.ys.end());
      trisInfo.emplace_back(newTriInfo);
    }
  }

  template<typename MeshPtsType, typename ScalarsType, typename ClipPtsType>
  static void Discard(TriMeta<MeshPtsType, ScalarsType>& triInfo, const ClipperMeta<ClipPtsType>& clipperInfo)
  {
    ///> Remove/Accept triangles
    /// 1. inverted is set, accept triangles inside a polygon.
    /// 2. inverted is unset, accept triangles outside a polygon.
    auto triCroid = triCentroid(triInfo.coords[0], triInfo.coords[1], triInfo.coords[2]);
    std::size_t iPoly(0);
    for (const auto& polyBbox : clipperInfo.polyBboxes)
    {
      if (!clipperInfo.invertedOp[iPoly])
      {
        if (fullyOutsideBbox(triInfo.bbox, polyBbox))
        {
          ++iPoly;
          continue;
        }
      }

      bool triInsidePoly = pnpoly(clipperInfo.xs[iPoly].size(), clipperInfo.xs[iPoly].data(), clipperInfo.ys[iPoly].data(), triCroid[0], triCroid[1]);

      if (clipperInfo.invertedOp[iPoly] && !triInsidePoly)
        triInfo.discard = true;
      else if (!clipperInfo.invertedOp[iPoly] && triInsidePoly)
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
          PtsType pASide = isLeft(pCcoord[0], pCcoord[1], pD_coord[0], pD_coord[1], pA_[0], pA_[1]);
          PtsType pBSide = isLeft(pCcoord[0], pCcoord[1], pD_coord[0], pD_coord[1], pB_[0], pB_[1]);
          bool sameSide = (pASide >= 0 && pBSide >= 0) || (pASide <= 0 && pBSide <= 0);
          if (!sameSide)
          {
            tri0 = nbrs[iNbr];
            pA = points[edge.first];
            pB = points[edge.second];
            auto cellIds = vtkSmartPointer<vtkIdList>::New();
            mesh->GetCellEdgeNeighbors(tri0, pA, pB, cellIds);
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
    SurfCutterImpl(vtkSmartPointer<vtkDataSet> mesh_, vtkSmartPointer<vtkPolyData> clipper_, const bool& insideOut_, const std::string& scalarsName_)
      : mesh(mesh_), clipper(clipper_), insideOut(insideOut_), scalarsName(scalarsName_)
    {}

    vtkSmartPointer<vtkDataSet> mesh;
    vtkSmartPointer<vtkPolyData> clipper;
    bool insideOut;
    std::string scalarsName;



    template<typename MeshPtsType, typename ScalarsType>
    vtkSmartPointer<vtkDataSet> CreateMesh(const std::vector<MeshMeta<MeshPtsType, ScalarsType>>& meshInfos, const std::vector<std::vector<vtkIdType>>& cells)
    {
      auto mesh = vtkSmartPointer<vtkPolyData>::New();

      vtkIdType numPoints = meshInfos.size();
      if (!numPoints)
        return mesh;

      auto pointsData = vtkSmartPointer<vtkAOSDataArrayTemplate<MeshPtsType>>::New();
      auto scalars = vtkSmartPointer<vtkAOSDataArrayTemplate<ScalarsType>>::New();
      auto acquisition = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();

      pointsData->SetNumberOfComponents(3);
      scalars->SetNumberOfComponents(1);
      acquisition->SetNumberOfComponents(1);

      pointsData->SetNumberOfTuples(numPoints);
      scalars->SetNumberOfTuples(numPoints);
      acquisition->SetNumberOfTuples(numPoints);

      scalars->SetName(scalarsName.c_str());
      acquisition->SetName("Acquired");

      for (vtkIdType tupleIdx = 0; tupleIdx < numPoints; ++tupleIdx)
      {
        for (int iDim = 0; iDim < 3; ++iDim)
          pointsData->SetTypedComponent(tupleIdx, iDim, meshInfos[tupleIdx].coord[iDim]);

        scalars->SetValue(tupleIdx, meshInfos[tupleIdx].scalar);
        acquisition->SetValue(tupleIdx, meshInfos[tupleIdx].pCrit);
      }

      auto pts = vtkSmartPointer<vtkPoints>::New();
      pts->SetData(pointsData);
      mesh->SetPoints(pts);
      mesh->GetPointData()->AddArray(scalars);
      mesh->GetPointData()->AddArray(acquisition);
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

    template<typename MeshPtsType, typename ScalarsType>
    vtkSmartPointer<vtkDataSet> CreateTriMesh(const std::vector<MeshMeta<MeshPtsType, ScalarsType>>& meshInfos, const std::vector<std::array<vtkIdType, 3>>& cells)
    {
      auto mesh = vtkSmartPointer<vtkPolyData>::New();

      vtkIdType numPoints = meshInfos.size();
      if (!numPoints)
        return mesh;

      auto pointsData = vtkSmartPointer<vtkAOSDataArrayTemplate<MeshPtsType>>::New();
      auto scalars = vtkSmartPointer<vtkAOSDataArrayTemplate<ScalarsType>>::New();
      auto acquisition = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();

      pointsData->SetNumberOfComponents(3);
      scalars->SetNumberOfComponents(1);
      acquisition->SetNumberOfComponents(1);

      pointsData->SetNumberOfTuples(numPoints);
      scalars->SetNumberOfTuples(numPoints);
      acquisition->SetNumberOfTuples(numPoints);

      scalars->SetName(scalarsName.c_str());
      acquisition->SetName("Acquired");

      for (vtkIdType tupleIdx = 0; tupleIdx < numPoints; ++tupleIdx)
      {
        for (int iDim = 0; iDim < 3; ++iDim)
          pointsData->SetTypedComponent(tupleIdx, iDim, meshInfos[tupleIdx].coord[iDim]);

        scalars->SetValue(tupleIdx, meshInfos[tupleIdx].scalar);
        acquisition->SetValue(tupleIdx, meshInfos[tupleIdx].pCrit);
      }

      auto pts = vtkSmartPointer<vtkPoints>::New();
      pts->SetData(pointsData);
      mesh->SetPoints(pts);
      mesh->GetPointData()->AddArray(scalars);
      mesh->GetPointData()->AddArray(acquisition);
      mesh->Allocate(cells.size());

      for (const auto& cell : cells)
      {
        mesh->InsertNextCell(VTK_TRIANGLE, 3, cell.data());
      }
      return mesh;
    }

    template<typename MeshPtsArrayT, typename ScalarsArrayT, typename ClipPointsArrayT>
    void operator()(MeshPtsArrayT* meshPts, ScalarsArrayT* scalars, ClipPointsArrayT* clipperPts)
    {
      const vtkIdType& numPolys = clipper->GetNumberOfCells();
      const vtkIdType& numCells = mesh->GetNumberOfCells();
      if (!(numPolys && numCells))
        return;

      const vtkIdType& numClipPoints = clipperPts->GetNumberOfTuples();
      const vtkIdType& numPoints = meshPts->GetNumberOfTuples();
      const vtkIdType& numScalars = scalars->GetNumberOfTuples();

      if (!numClipPoints)
        return;

      if (numPoints != numScalars)
        return;

      if (!numClipPoints)
        return;

      auto invertOps = vtkAOSDataArrayTemplate<int>::FastDownCast(clipper->GetCellData()->GetArray("Invert"));
#if 0
      if (!invertOps)
        return;
#endif

      auto acquisition = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();
      acquisition->SetNumberOfComponents(1);
      acquisition->SetNumberOfTuples(numPoints);
      acquisition->SetName("Acquired");
      acquisition->FillValue(0);
      mesh->GetPointData()->AddArray(acquisition);


      using MeshPtsAccess = vtkDataArrayAccessor<MeshPtsArrayT>;
      using ScalarsAccess = vtkDataArrayAccessor<ScalarsArrayT>;
      using ClipPointsAccess = vtkDataArrayAccessor<ClipPointsArrayT>;
      using IntArrAccess = vtkDataArrayAccessor<vtkAOSDataArrayTemplate<int>>;

      using MeshPtsType = MeshPtsAccess::APIType;
      using ScalarsType = ScalarsAccess::APIType;
      using ClipPointsType = ClipPointsAccess::APIType;

      using MMeta = MeshMeta<MeshPtsType, ScalarsType>;
      using TMeta = TriMeta<MeshPtsType, ScalarsType>;
      using CMeta = ClipperMeta<ClipPointsType>;

      VTK_ASSUME(meshPts->GetNumberOfComponents() == 3);
      VTK_ASSUME(scalars->GetNumberOfComponents() == 1);
      VTK_ASSUME(clipperPts->GetNumberOfComponents() == 3);

      MeshPtsAccess meshPtsAccess(meshPts);
      ClipPointsAccess clipPointsAccess(clipperPts);
      ScalarsAccess scalarsAccess(scalars);
      IntArrAccess invertOpsAccess(invertOps);

      CMeta clipperInfo;
      TMeta triInfo;

      ///> clipper info.
      for (vtkIdType iPt = 0; iPt < numClipPoints; ++iPt)
      {
        std::array<ClipPointsType, 3> coord;
        clipPointsAccess.Get(iPt, coord.data());
        clipperInfo.coords.emplace_back(coord);
      }

      vtkSmartPointer<vtkCellIterator> polysIter = clipper->NewCellIterator();
      vtkIdType iPoly(0);
      for (polysIter->InitTraversal(); !polysIter->IsDoneWithTraversal(); polysIter->GoToNextCell(), ++iPoly)
      {
        vtkSmartPointer<vtkIdList> ptIds = polysIter->GetPointIds();
        const vtkIdType& numIds = ptIds->GetNumberOfIds();
        std::vector<ClipPointsType> _xs(numIds), _ys(numIds);
        std::vector<std::array<ClipPointsType, 3>> poly(numIds);

        for (vtkIdType idx = 0; idx < numIds; ++idx)
        {
          const vtkIdType& p0Id = idx;
          const vtkIdType p1Id = (idx + 1) % numIds ? idx + 1 : 0;
          const vtkIdType& p0 = ptIds->GetId(p0Id);
          const vtkIdType& p1 = ptIds->GetId(p1Id);
          clipperInfo.edges.emplace_back(p0, p1);
          clipperInfo.edgeBboxes.emplace_back(BBox<ClipPointsType>(clipperInfo.coords[p0], clipperInfo.coords[p1]));
          _xs[idx] = clipperInfo.coords[p0][0];
          _ys[idx] = clipperInfo.coords[p0][1];
        }
        clipperInfo.xs.emplace_back(_xs);
        clipperInfo.ys.emplace_back(_ys);
        clipperInfo.polyBboxes.emplace_back(BBox<ClipPointsType>(_xs, _ys));
#if 0
        clipperInfo.invertedOp.emplace_back(invertOpsAccess.Get(iPoly, 0));
#endif
      }

      std::vector<std::array<MeshPtsType, 3>> clipperCoords;
      clipperCoords.reserve(numClipPoints + 100); // extra space for self-intersections.
      for (const auto& clipperCoord : clipperInfo.coords)
      {
        std::array<MeshPtsType, 3> newCoord;
        std::transform(clipperCoord.begin(), clipperCoord.end(), newCoord.begin(), [&](const auto& val) {return static_cast<MeshPtsType>(val); });
        clipperCoords.emplace_back(newCoord);
      }

      std::vector<MMeta> meshInfo(numPoints);
      std::vector<std::array<vtkIdType, 3>> tris;
      std::vector<TMeta> trisInfo;

      for (vtkIdType tupleIdx = 0; tupleIdx < numPoints; ++tupleIdx)
      {
        meshPtsAccess.Get(tupleIdx, meshInfo[tupleIdx].coord.data());
        meshInfo[tupleIdx].scalar = scalarsAccess.Get(tupleIdx, 0);
        meshInfo[tupleIdx].pCrit = acquisition->GetValue(tupleIdx);
        meshInfo[tupleIdx].originalIdx = tupleIdx;
      }
      trisInfo.reserve(trisInfo.size() + numCells);
      tris.reserve(tris.size() + numCells);

      ExtractTris(mesh, tris, meshInfo, meshInfo, trisInfo);

      meshInfo.reserve(meshInfo.size() + numClipPoints * 3 * 3); // 3 new triangles per clip point, 3 new intersection points per triangle.
      trisInfo.reserve(trisInfo.size() + numClipPoints * 3 * 3 * 3); // same as meshInfo, further 3 sub triangles per each new triangle.
      tris.reserve(tris.size() + numClipPoints * 3 * 3 * 3); // same as trisInfo

      // for each tri.
      // interp.
      // build submesh
      // add new tris.
      // for each tri.
      // clip.
      // build submesh
      // add new tris.
      // remove inside/outside.
      const std::size_t numTris1 = trisInfo.size();
      for (std::size_t iTri = 0; iTri < numTris1; ++iTri)
      {
        TMeta& triInfo = trisInfo[iTri];
        std::vector<MMeta> smInfo;
        smInfo.reserve(clipperCoords.size());

        for (vtkIdType iPt = 0; iPt < 3; ++iPt)
        {
          smInfo.emplace_back(MMeta(triInfo.coords[iPt], triInfo.scalars[iPt], triInfo.pCrits[iPt], triInfo.verts[iPt]));
        }

        ///> A The sub-mesh after interpolation might consist of 
        /// 1. the actual triangle.
        /// 2. all polyPoints that are enclosed by that triangle.
        for (auto& clipperCoord : clipperCoords)
        {
          if (!pnpoly(3, triInfo.xs.data(), triInfo.ys.data(), clipperCoord[0], clipperCoord[1]))
          {
            // check corners.
            for (std::size_t triVert = 0; triVert < 3; ++triVert)
            {
              bool xSimilar = std::abs(triInfo.xs[triVert] - clipperCoord[0]) <= std::numeric_limits<float>::epsilon();
              bool ySimilar = std::abs(triInfo.ys[triVert] - clipperCoord[1]) <= std::numeric_limits<float>::epsilon();
              if (xSimilar && ySimilar)
                acquisition->SetValue(triInfo.verts[triVert], 1);
            }
            continue; // move on to next point.
          }

          ///> Fastest and sensible interpolation using bary centric coords.
          ScalarsType newScalar = TriInterp<MeshPtsType, ScalarsType>(clipperCoord, triInfo.coords, triInfo.scalars);
          auto smPt = MMeta(clipperCoord, newScalar, int(1), meshInfo.size());
          smInfo.emplace_back(smPt);
          meshInfo.emplace_back(smPt);
        } // end interpolation loop

        if (smInfo.size() > 3)
        {
          triInfo.discard = true;
          // setup constraints.
          vtkSmartPointer<vtkPolyData> subMesh = vtkPolyData::SafeDownCast(CreateMesh(smInfo, TRIEDGESEGS));
          // triangulate with constraints.
          auto del2d = vtkSmartPointer<vtkDelaunay2D>::New();
          del2d->SetInputData(subMesh);
          del2d->SetSourceData(subMesh);
          del2d->Update();
          subMesh->ShallowCopy(del2d->GetOutput());
          ExtractTris(subMesh, tris, smInfo, meshInfo, trisInfo);
        }
        else
        {
          continue;
        }
      }

      ///> B Now insert clipper's edges into triangles. 
      /// 1. An intersection of any triEdge with clipperEdge. (points)
      /// 2. A triEdge that does not intersect with any clipperEdge. (lines)
      /// 3. All permutations of intersection with triPts. (lines)
      const std::size_t numTris2 = trisInfo.size();
      for (std::size_t iTri = 0; iTri < numTris2; ++iTri)
      {
        TMeta& triInfo = trisInfo[iTri];
        if (triInfo.discard)
          continue;

        std::vector<MMeta> smInfo;
        smInfo.reserve(clipperCoords.size());

        for (vtkIdType iPt = 0; iPt < 3; ++iPt)
        {
          smInfo.emplace_back(MMeta(triInfo.coords[iPt], triInfo.scalars[iPt], triInfo.pCrits[iPt], triInfo.verts[iPt]));
        }

        std::vector<std::array<MeshPtsType, 3>> triEdgeCrossings;
        bool triCrossesClipper(false);
        triEdgeCrossings.reserve(numClipPoints);

        vtkIdType iClipEdge(-1);
        std::map<vtkIdType, vtkIdType> constraints;
        for (const auto& clipperEdge : clipperInfo.edges) // for each clipper's edge
        {
          std::vector<vtkIdType> clipperEdgeCrossings;
          std::vector<vtkIdType> triHitAt;

          ++iClipEdge;
          const auto& pa_0 = clipperCoords[clipperEdge.first];
          const auto& pb_0 = clipperCoords[clipperEdge.second];
          const auto& bbox_0 = clipperInfo.edgeBboxes[iClipEdge];

          if (!intersectBbox(bbox_0, triInfo.bbox))
            continue;

          clipperEdgeCrossings.reserve(3);
          triHitAt.reserve(3);
          for (vtkIdType iTriEdge = 0; iTriEdge < 3; ++iTriEdge) // for each triangle's edge
          {
            const std::pair<vtkIdType, vtkIdType>& triEdge = TRIEDGES[iTriEdge];
            const auto& pa_1 = triInfo.coords[triEdge.first];
            const auto& pb_1 = triInfo.coords[triEdge.second];
            std::array<MeshPtsType, 3> pc;

            auto junctionType = intersect(pa_0, pb_0, pa_1, pb_1, pc);
            if (junctionType == JunctionType::CROSS)
            {
              if (pc[0] == pa_0[0] && pc[1] == pa_0[1])
                continue;

              if (pc[0] == pb_0[0] && pc[1] == pb_0[1])
                continue;

              auto newScalar = TriInterp<MeshPtsType, ScalarsType>(pc, triInfo.coords, triInfo.scalars);

              triCrossesClipper |= true;

              auto smPt = MMeta(pc, newScalar, int(0), meshInfo.size()); // flip 0 to 1, does it make a diff?
              smInfo.emplace_back(smPt);
              meshInfo.emplace_back(smPt);

              triEdgeCrossings.emplace_back(pc);
              clipperEdgeCrossings.emplace_back(static_cast<vtkIdType>(smInfo.size() - 1));
              triHitAt.emplace_back(iTriEdge);
            }
            else
            {
              triCrossesClipper |= false;
            }
          } // end for triangle's edge

          if (clipperEdgeCrossings.size() == 2)
          {
            constraints[clipperEdgeCrossings[0]] = clipperEdgeCrossings[1];
          }
          else if (clipperEdgeCrossings.size() == 1)
          {
            const vtkIdType crossedTriEdge = triHitAt[0];
            const vtkIdType oppositeVtx = (TRIEDGES[crossedTriEdge].second + 1) % 3 ? TRIEDGES[crossedTriEdge].second + 1 : 0;
            constraints[clipperEdgeCrossings[0]] = oppositeVtx;
          }

        } // end for clipper's edge
        if (!triCrossesClipper)
        {
#if 0
          Discard(triInfo, clipperInfo);
#endif
          continue;
        }
        triInfo.discard = true;
        std::vector<std::array<MeshPtsType, 3>> smCoords;
        for (const auto& info : smInfo)
        {
          smCoords.emplace_back(info.coord);
        }
        auto centroid = triCentroid(triInfo.coords[0], triInfo.coords[1], triInfo.coords[2]);
        std::vector<std::size_t> ccwIds_;
        counterClockWise(smCoords, centroid, ccwIds_);
        std::vector<std::vector<vtkIdType>> cells;
        cells.emplace_back(std::vector<vtkIdType>(ccwIds_.begin(), ccwIds_.end()));

        vtkSmartPointer<vtkDataSet> subMesh = CreateMesh<MeshPtsType, ScalarsType>(smInfo, cells);
        auto subMeshTris = vtkSmartPointer<vtkPolyData>::New();
        auto triangulate = vtkSmartPointer<vtkTriangleFilter>::New();
        triangulate->SetInputData(subMesh);
        triangulate->Update();
        subMeshTris->ShallowCopy(triangulate->GetOutput());
        subMeshTris->BuildLinks();
        
        for (const auto& constraint : constraints)
        {
          unsigned short niters = 0;
          while (!subMeshTris->IsEdge(constraint.first, constraint.second))
          {
            using dispatcher = vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals>;
            auto constraintWorker = ApplyConstraintImpl(subMeshTris, constraint.first, constraint.second);
            vtkDataArray* points = subMeshTris->GetPoints()->GetData();
            if (!dispatcher::Execute(points, constraintWorker))
            {
              constraintWorker(points);
            }
            if (niters > 128)
            {
              break;
            }
          }
        }
        const std::size_t oldSz = trisInfo.size();
        ExtractTris(subMeshTris, tris, smInfo, meshInfo, trisInfo);
        const std::size_t newSz = trisInfo.size();

#if 0
        for (std::size_t iInfo = oldSz; iInfo < newSz; ++iInfo)
        {
          Discard(trisInfo[iInfo], clipperInfo);
        }
#endif
      }

      std::vector<std::array<vtkIdType, 3>> newTris;
      newTris.reserve(trisInfo.size());
      for (const auto& triInfo : trisInfo)
      {
        if (triInfo.discard)
          continue;
        newTris.emplace_back(triInfo.verts);
      }
      vtkSmartPointer<vtkDataSet> meshDset = CreateTriMesh(meshInfo, newTris);
      vtkSmartPointer<vtkCellArray> cells = vtkPolyData::SafeDownCast(meshDset)->GetPolys();
      mesh->ShallowCopy(meshDset);
      if (mesh->IsA("vtkPolyData"))
      {
        auto meshPd = vtkPolyData::SafeDownCast(mesh);
        meshPd->SetPolys(cells);
      }
      else if (mesh->IsA("vtkUnstructuredGrid"))
      {
        auto meshUgrid = vtkUnstructuredGrid::SafeDownCast(mesh);
        meshUgrid->SetCells(VTK_TRIANGLE, cells);
      }
    }
  };
}
// anon end

int SurfaceCutter::RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkSmartPointer<vtkDataSet> mesh = vtkDataSet::GetData(inputVector[0]->GetInformationObject(0));
  vtkSmartPointer<vtkPolyData> clipper = vtkPolyData::GetData(inputVector[1]->GetInformationObject(0));
  vtkSmartPointer<vtkDataSet> outMesh = vtkDataSet::GetData(outputVector->GetInformationObject(0));


  if (!clipper->GetNumberOfPoints())
  {
    outMesh->ShallowCopy(mesh);
    return 1;
  }

  vtkDataArray* clipperPts = clipper->GetPoints()->GetData();

  using dispatcher = vtkArrayDispatch::Dispatch3ByValueType<vtkArrayDispatch::Reals, vtkArrayDispatch::Reals, vtkArrayDispatch::Reals>;
  if (mesh->IsA("vtkPolyData"))
  {
    vtkSmartPointer<vtkPolyData> meshPd = vtkPolyData::SafeDownCast(mesh);
    vtkDataArray* points = meshPd->GetPoints()->GetData();
    vtkDataArray* scalars = this->GetInputArrayToProcess(0, meshPd);
    auto worker = SurfCutterImpl(meshPd, clipper, this->InsideOut, scalars->GetName());
    if (!dispatcher::Execute(points, scalars, clipperPts, worker))
    {
      worker(points, scalars, clipperPts);
    }
    vtkSmartPointer<vtkPolyData> outMeshPd = vtkPolyData::SafeDownCast(outMesh);
    outMeshPd->ShallowCopy(meshPd);
  }  
  else if (mesh->IsA("vtkUnstructuredGrid"))
  {
    vtkSmartPointer<vtkUnstructuredGrid> meshUgrid = vtkUnstructuredGrid::SafeDownCast(mesh);
    vtkDataArray* points = meshUgrid->GetPoints()->GetData();
    vtkDataArray* scalars = this->GetInputArrayToProcess(0, meshUgrid);
    auto worker = SurfCutterImpl(meshUgrid, clipper, this->InsideOut, scalars->GetName());
    if (!dispatcher::Execute(points, scalars, clipperPts, worker))
    {
      worker(points, scalars, clipperPts);
    }
    vtkSmartPointer<vtkUnstructuredGrid> outMeshUgrid = vtkUnstructuredGrid::SafeDownCast(outMesh);
    outMeshUgrid->ShallowCopy(meshUgrid);
  }

  outMesh->ShallowCopy(mesh);
  return 1;
}
