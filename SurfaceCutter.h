/**
 * @class   SurfaceCutter
 * @brief   Cut triangulated surfaces with polygons
 *
 * SurfaceCutter geometrically cuts a triangulated surface.
 * It differs from vtkClipDataSet which is a Scalar-based clip operation.
 *
 * Loop polygons from the second input port will be embedded onto the surface.
 * After a re-triangulation of triangles that contain an embedded edge,
 * the result triangles will be rejected/accepted if necessary. See SetInsideOut()
 * This is decided with a point-in-polygon test. Even handles situation where 
 * a polygon's point might coincide with a triangle's edge or a vertex.
 *
 * @note PointData is interpolated to output.
 * CellData is copied over to both constraint lines, new triangles
 *
 * @sa
 * vtkClipDataSet vtkClipPolyData
 *
 */

#ifndef SurfaceCutter_h__
#define SurfaceCutter_h__

#include "vtkPolyDataAlgorithm.h"

#include "vtkAbstractCellLocator.h"
#include "vtkIncrementalPointLocator.h"
#include "vtkSmartPointer.h"

class SurfaceCutter : public vtkPolyDataAlgorithm
{
public:
  /**
   * Construct object with tolerance 1.0e-6, inside out set to true,
   * color acquired points, color loop edges
   */
  static SurfaceCutter* New();
  vtkTypeMacro(SurfaceCutter, vtkPolyDataAlgorithm);

  //@{
  /**
   * Append an array to output point data that colors acquired points. Default: On
   */
  vtkBooleanMacro(ColorAcquiredPts, bool);
  vtkSetMacro(ColorAcquiredPts, bool);
  vtkGetMacro(ColorAcquiredPts, bool);
  //@}

  //@{
  /**
   * Append an array to output cell data which colors constrained lines. Default: On
   */
  vtkBooleanMacro(ColorLoopEdges, bool);
  vtkSetMacro(ColorLoopEdges, bool);
  vtkGetMacro(ColorLoopEdges, bool);
  //@}

  //@{
  /**
   * After the loop's edges are embedded onto the surface,
   * On: remove stuff outside loop
   * Off: remove stuff inside loop
   */
  vtkBooleanMacro(InsideOut, bool);
  vtkSetMacro(InsideOut, bool);
  vtkGetMacro(InsideOut, bool);
  //@}

  //@{
  /**
   * Tolerance for point merging.
   */
  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);
  //@}

  //@{
  /**
   * Any subclass of vtkAbstractCellLocator that implements the method 'FindCellsWithinBounds()'.
   * Ex: vtkStaticCellLocator, vtkCellLocator. Not vtkOBBTree
   */
  vtkSetSmartPointerMacro(CellLocator, vtkAbstractCellLocator);
  vtkGetSmartPointerMacro(CellLocator, vtkAbstractCellLocator);
  //@}

  //@{
  /**
   * Any subclass of vtkIncrementalPointLocator. Prefer vtkMergePoints (the default) for faster
   * execution.
   */
  vtkSetSmartPointerMacro(PointLocator, vtkIncrementalPointLocator);
  vtkGetSmartPointerMacro(PointLocator, vtkIncrementalPointLocator);
  //@}

  /**
   * Provide loop polydata (collection of polygons)
   */
  void SetLoops(vtkPointSet* loops);

  /**
   * Provide loop polydata (collection of polygons)
   */
  void SetLoopsConnection(vtkAlgorithmOutput* output);

protected:
  SurfaceCutter();
  ~SurfaceCutter();

  bool ColorAcquiredPts;
  bool ColorLoopEdges;
  bool InsideOut;
  double Tolerance;
  vtkSmartPointer<vtkAbstractCellLocator> CellLocator;
  vtkSmartPointer<vtkIncrementalPointLocator> PointLocator;

  int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

private:
  SurfaceCutter(const SurfaceCutter&) = delete;
  void operator=(const SurfaceCutter&) = delete;
};

#endif // SurfaceCutter_h__
