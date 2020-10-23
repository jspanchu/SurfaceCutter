#ifndef SurfaceCutterImpl_h__
#define SurfaceCutterImpl_h__


#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <utility>

#include <vtkArrayDispatch.h>
#include <vtkAppendFilter.h>
#include <vtkAOSDataArrayTemplate.h>
#include <vtkCellData.h>
#include <vtkCellIterator.h>
#include <vtkDelaunay2D.h>
#include <vtkDataArray.h>
#include <vtkDataArrayAccessor.h>
#include <vtkExtractCells.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkXMLPolyDataWriter.h>

#include "compgeom/bbox.h"
#include "compgeom/orientation.h"
#include "compgeom/polygon.h"
#include "compgeom/intersect2d.h"

namespace tiara
{
  namespace detail
  {
    static const std::vector<std::vector<vtkIdType>> TRIEDGESEGS = { {0, 1}, {1, 2}, {2, 0} };

    static const std::vector<std::vector<vtkIdType>> TRIPTIDS = { {0, 1, 2} };

    static const std::vector<std::pair<vtkIdType, vtkIdType>> TRIEDGES = { {0, 1}, {1, 2}, {2, 0} };

    static std::string activeScalars;

    template<typename ClipPtsType>
    struct ClipperMeta {
      std::vector<std::array<ClipPtsType, 3>> coords;
      std::vector<std::pair<vtkIdType, vtkIdType>> edges;
      std::vector<compgeom::bbox::BBox<ClipPtsType>> edgeBboxes, polyBboxes;
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
      compgeom::bbox::BBox<MeshPtsType> bbox;
      bool discard = false;
    };

    template<typename MeshPtsType, typename ScalarsType>
    static vtkSmartPointer<vtkDataSet> CreateMesh(const std::vector<tiara::detail::MeshMeta<MeshPtsType, ScalarsType>>& meshInfos, const std::vector<std::array<vtkIdType, 3>>& tris)
    {
      auto mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

      vtkIdType numPoints = meshInfos.size();
      if (!numPoints)
        return mesh;

      auto pointsData = vtkSmartPointer<vtkAOSDataArrayTemplate<MeshPtsType>>::New();
      auto scalars = vtkSmartPointer<vtkAOSDataArrayTemplate<ScalarsType>>::New();
      auto pseudoCrits = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();

      pointsData->SetNumberOfComponents(3);
      scalars->SetNumberOfComponents(1);
      pseudoCrits->SetNumberOfComponents(1);

      pointsData->SetNumberOfTuples(numPoints);
      scalars->SetNumberOfTuples(numPoints);
      pseudoCrits->SetNumberOfTuples(numPoints);

      scalars->SetName(activeScalars.c_str());
      pseudoCrits->SetName("PseudoCritical");

      for (vtkIdType tupleIdx = 0; tupleIdx < numPoints; ++tupleIdx)
      {
        for (int iDim = 0; iDim < 3; ++iDim)
          pointsData->SetTypedComponent(tupleIdx, iDim, meshInfos[tupleIdx].coord[iDim]);

        scalars->SetValue(tupleIdx, meshInfos[tupleIdx].scalar);
        pseudoCrits->SetValue(tupleIdx, meshInfos[tupleIdx].pCrit);
      }

      auto pts = vtkSmartPointer<vtkPoints>::New();
      pts->SetData(pointsData);
      mesh->SetPoints(pts);
      mesh->GetPointData()->AddArray(scalars);
      mesh->GetPointData()->AddArray(pseudoCrits);
      mesh->Allocate(tris.size());

      for (const auto& tri : tris)
      {
        mesh->InsertNextCell(VTK_TRIANGLE, 3, tri.data());
      }
      return mesh;
    }

    template<typename MeshPtsType, typename ScalarsType>
    static vtkSmartPointer<vtkDataSet> CreateMesh(const std::vector<tiara::detail::MeshMeta<MeshPtsType, ScalarsType>>& meshInfos, const std::vector<std::vector<vtkIdType>>& cells)
    {
      auto mesh = vtkSmartPointer<vtkPolyData>::New();

      vtkIdType numPoints = meshInfos.size();
      if (!numPoints)
        return mesh;

      auto pointsData = vtkSmartPointer<vtkAOSDataArrayTemplate<MeshPtsType>>::New();
      auto scalars = vtkSmartPointer<vtkAOSDataArrayTemplate<ScalarsType>>::New();
      auto pseudoCrits = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();

      pointsData->SetNumberOfComponents(3);
      scalars->SetNumberOfComponents(1);
      pseudoCrits->SetNumberOfComponents(1);

      pointsData->SetNumberOfTuples(numPoints);
      scalars->SetNumberOfTuples(numPoints);
      pseudoCrits->SetNumberOfTuples(numPoints);

      scalars->SetName(activeScalars.c_str());
      pseudoCrits->SetName("PseudoCritical");

      for (vtkIdType tupleIdx = 0; tupleIdx < numPoints; ++tupleIdx)
      {
        for (int iDim = 0; iDim < 3; ++iDim)
          pointsData->SetTypedComponent(tupleIdx, iDim, meshInfos[tupleIdx].coord[iDim]);

        scalars->SetValue(tupleIdx, meshInfos[tupleIdx].scalar);
        pseudoCrits->SetValue(tupleIdx, meshInfos[tupleIdx].pCrit);
      }

      auto pts = vtkSmartPointer<vtkPoints>::New();
      pts->SetData(pointsData);
      mesh->SetPoints(pts);
      mesh->GetPointData()->AddArray(scalars);
      mesh->GetPointData()->AddArray(pseudoCrits);
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
    static ScalarsType TriInterp(std::array<MeshPtsType, 3>& coord, const std::array<std::array<MeshPtsType, 3>, 3>& coords, const std::array<ScalarsType, 3>& scalars)
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
      std::vector<tiara::detail::MeshMeta<MeshPtsType, ScalarsType>>& smInfo,
      std::vector<tiara::detail::MeshMeta<MeshPtsType, ScalarsType>>& meshInfo,
      std::vector<tiara::detail::TriMeta<MeshPtsType, ScalarsType>>& trisInfo)
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

        tiara::detail::TriMeta<MeshPtsType, ScalarsType> newTriInfo;
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
    static void Discard(tiara::detail::TriMeta<MeshPtsType, ScalarsType>& triInfo, const tiara::detail::ClipperMeta<ClipPtsType>& clipperInfo)
    {
      ///> Remove/Accept triangles
      /// 1. inverted is set, accept triangles inside a polygon.
      /// 2. inverted is unset, accept triangles outside a polygon.
      auto triCroid = compgeom::polygon::triCentroid(triInfo.coords[0], triInfo.coords[1], triInfo.coords[2]);
      std::size_t iPoly(0);
      for (const auto& polyBbox : clipperInfo.polyBboxes)
      {
        if (!clipperInfo.invertedOp[iPoly])
        {
          if (compgeom::bbox::fullyOutsideBbox(triInfo.bbox, polyBbox))
          {
            ++iPoly;
            continue;
          }
        }

        bool triInsidePoly = compgeom::polygon::pnpoly(clipperInfo.xs[iPoly].size(), clipperInfo.xs[iPoly].data(), clipperInfo.ys[iPoly].data(), triCroid[0], triCroid[1]);

        if (clipperInfo.invertedOp[iPoly] && !triInsidePoly)
          triInfo.discard = true;
        else if (!clipperInfo.invertedOp[iPoly] && triInsidePoly)
          triInfo.discard = true;

        ++iPoly;
      }
    }

    struct DatasetCutter
    {
      DatasetCutter(vtkSmartPointer<vtkDataSet> mesh, vtkSmartPointer<vtkPolyData> clipper)
        : mesh(mesh), clipper(clipper)
      {}

      vtkSmartPointer<vtkDataSet> mesh;
      vtkSmartPointer<vtkPolyData> clipper;

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

        auto pseudoCrits = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();
        pseudoCrits->SetNumberOfComponents(1);
        pseudoCrits->SetNumberOfTuples(numPoints);
        pseudoCrits->SetName("PseudoCritical");
        pseudoCrits->FillValue(0);
        mesh->GetPointData()->AddArray(pseudoCrits);


        using MeshPtsAccess = vtkDataArrayAccessor<MeshPtsArrayT>;
        using ScalarsAccess = vtkDataArrayAccessor<ScalarsArrayT>;
        using ClipPointsAccess = vtkDataArrayAccessor<ClipPointsArrayT>;
        using IntArrAccess = vtkDataArrayAccessor<vtkAOSDataArrayTemplate<int>>;

        using MeshPtsType = MeshPtsAccess::APIType;
        using ScalarsType = ScalarsAccess::APIType;
        using ClipPointsType = ClipPointsAccess::APIType;

        using MMeta = tiara::detail::MeshMeta<MeshPtsType, ScalarsType>;
        using TMeta = tiara::detail::TriMeta<MeshPtsType, ScalarsType>;
        using CMeta = tiara::detail::ClipperMeta<ClipPointsType>;

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
            clipperInfo.edgeBboxes.emplace_back(compgeom::bbox::BBox<ClipPointsType>(clipperInfo.coords[p0], clipperInfo.coords[p1]));
            _xs[idx] = clipperInfo.coords[p0][0];
            _ys[idx] = clipperInfo.coords[p0][1];
          }
          clipperInfo.xs.emplace_back(_xs);
          clipperInfo.ys.emplace_back(_ys);
          clipperInfo.polyBboxes.emplace_back(compgeom::bbox::BBox<ClipPointsType>(_xs, _ys));
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
          meshInfo[tupleIdx].pCrit = pseudoCrits->GetValue(tupleIdx);
          meshInfo[tupleIdx].originalIdx = tupleIdx;
        }
        trisInfo.reserve(trisInfo.size() + numCells);
        tris.reserve(tris.size() + numCells);

        tiara::detail::ExtractTris(mesh, tris, meshInfo, meshInfo, trisInfo);

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
            if (!compgeom::polygon::pnpoly(3, triInfo.xs.data(), triInfo.ys.data(), clipperCoord[0], clipperCoord[1]))
            {
              // check corners.
              for (std::size_t triVert = 0; triVert < 3; ++triVert)
              {
                bool xSimilar = std::abs(triInfo.xs[triVert] - clipperCoord[0]) <= std::numeric_limits<float>::epsilon();
                bool ySimilar = std::abs(triInfo.ys[triVert] - clipperCoord[1]) <= std::numeric_limits<float>::epsilon();
                if (xSimilar && ySimilar)
                  pseudoCrits->SetValue(triInfo.verts[triVert], 1);
              }
              continue; // move on to next point.
            }

            ///> Fastest and sensible interpolation using bary centric coords.
            ScalarsType newScalar = tiara::detail::TriInterp<MeshPtsType, ScalarsType>(clipperCoord, triInfo.coords, triInfo.scalars);
            auto smPt = MMeta(clipperCoord, newScalar, int(1), meshInfo.size());
            smInfo.emplace_back(smPt);
            meshInfo.emplace_back(smPt);
          } // end interpolation loop

          if (smInfo.size() > 3)
          {
            triInfo.discard = true;
            vtkSmartPointer<vtkPolyData> subMesh = vtkPolyData::SafeDownCast(tiara::detail::CreateMesh(smInfo, TRIEDGESEGS));
            auto del2d = vtkSmartPointer<vtkDelaunay2D>::New();
            del2d->SetInputData(subMesh);
            del2d->SetSourceData(subMesh);
            del2d->Update();
            subMesh->ShallowCopy(del2d->GetOutput());
            tiara::detail::ExtractTris(subMesh, tris, smInfo, meshInfo, trisInfo);
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

            if (!compgeom::bbox::intersectBbox(bbox_0, triInfo.bbox))
              continue;

            clipperEdgeCrossings.reserve(3);
            triHitAt.reserve(3);
            for (vtkIdType iTriEdge = 0; iTriEdge < 3; ++iTriEdge) // for each triangle's edge
            {
              const std::pair<vtkIdType, vtkIdType>& triEdge = TRIEDGES[iTriEdge];
              const auto& pa_1 = triInfo.coords[triEdge.first];
              const auto& pb_1 = triInfo.coords[triEdge.second];
              std::array<MeshPtsType, 3> pc;

              auto junctionType = compgeom::intersect2d::intersect(pa_0, pb_0, pa_1, pb_1, pc);
              if (junctionType == compgeom::intersect2d::JunctionType::CROSS)
              {
                if (pc[0] == pa_0[0] && pc[1] == pa_0[1])
                  continue;

                if (pc[0] == pb_0[0] && pc[1] == pb_0[1])
                  continue;

                auto newScalar = tiara::detail::TriInterp<MeshPtsType, ScalarsType>(pc, triInfo.coords, triInfo.scalars);

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
            tiara::detail::Discard(triInfo, clipperInfo);
#endif
            continue;
          }
          triInfo.discard = true;
          std::vector<std::array<MeshPtsType, 3>> smCoords;
          for (const auto& info : smInfo)
          {
            smCoords.emplace_back(info.coord);
          }
          auto centroid = compgeom::polygon::triCentroid(triInfo.coords[0], triInfo.coords[1], triInfo.coords[2]);
          std::vector<std::size_t> ccwIds_;
          compgeom::orientation::counterClockWise(smCoords, centroid, ccwIds_);
          std::vector<std::vector<vtkIdType>> cells;
          cells.emplace_back(std::vector<vtkIdType>(ccwIds_.begin(), ccwIds_.end()));

          vtkSmartPointer<vtkDataSet> subMesh = tiara::detail::CreateMesh<MeshPtsType, ScalarsType>(smInfo, cells);
          auto subMeshTris = vtkSmartPointer<vtkPolyData>::New();
          auto triangulate = vtkSmartPointer<vtkTriangleFilter>::New();
          triangulate->SetInputData(subMesh);
          triangulate->Update();
          subMeshTris->ShallowCopy(triangulate->GetOutput());
          subMeshTris->BuildLinks();
          auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
          writer->SetFileName("SubMeshTri.vtp");
          writer->SetInputData(subMeshTris);
          writer->Write();
          for (const auto& constraint : constraints)
          {
            unsigned short niters = 0;
            while (!subMeshTris->IsEdge(constraint.first, constraint.second))
            {
              InsertConstraint(subMeshTris, constraint.first, constraint.second);
              if (niters > 128)
              {
                __debugbreak();
                break;
              }
            }
          }
          const std::size_t oldSz = trisInfo.size();
          tiara::detail::ExtractTris(subMeshTris, tris, smInfo, meshInfo, trisInfo);
          const std::size_t newSz = trisInfo.size();

#if 0
          for (std::size_t iInfo = oldSz; iInfo < newSz; ++iInfo)
          {
            tiara::detail::Discard(trisInfo[iInfo], clipperInfo);
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
        vtkSmartPointer<vtkDataSet> meshDset = tiara::detail::CreateMesh(meshInfo, newTris);
        vtkSmartPointer<vtkCellArray> cells = vtkUnstructuredGrid::SafeDownCast(meshDset)->GetCells();
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

      static void InsertConstraint(vtkSmartPointer<vtkPolyData> mesh, const vtkIdType& p1, const vtkIdType& p2)
      {
        //throw std::logic_error("The method or operation is not implemented.");

        const vtkIdType& pC = p1;
        const vtkIdType& pD_ = p2;

        std::array<double, 3> pCcoord, pD_coord;
        mesh->GetPoint(pC, pCcoord.data());
        mesh->GetPoint(pD_, pD_coord.data());

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
            std::array<double, 3> pA_, pB_;
            mesh->GetPoint(points[edge.first], pA_.data());
            mesh->GetPoint(points[edge.second], pB_.data());
            double pASide = compgeom::orientation::isLeft(pCcoord[0], pCcoord[1], pD_coord[0], pD_coord[1], pA_[0], pA_[1]);
            double pBSide = compgeom::orientation::isLeft(pCcoord[0], pCcoord[1], pD_coord[0], pD_coord[1], pB_[0], pB_[1]);
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
          EdgeFlip(mesh, pA, pB, pC, pD, tri0, tri1);
        }
      }

      static void EdgeFlip(vtkSmartPointer<vtkPolyData> mesh, const vtkIdType& pA, const vtkIdType& pB, const vtkIdType& pC, const vtkIdType& pD, const vtkIdType& triABC, const vtkIdType triABD)
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
  }
}
#endif // SurfaceCutterImpl_h__