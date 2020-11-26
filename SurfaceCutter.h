/**
 * @class   SurfaceCutter
 * @brief   Cut triangulated surfaces with polygons
 *
 * Cut a triangulated surface with one or more polygons.
 * It differs from vtkClipDataSet which is a Scalar-based clip operation.
 *
 * This filter crops an input vtkPolyData consisting of triangles
 * with loops specified by a second input containing polygons. 
 * Note that this filter can handle concave polygons. It only produces triangles
 * and line segments (which are inherited from given loop's edges)
 * 
 * The result triangles will be rejected/accepted if necessary. See SetInsideOut()
 * This is decided with a point-in-polygon test. It also handles situation where 
 * a polygon's point might coincide with a triangle's edge or a vertex.
 *
 * @note PointData is interpolated to output.
 * CellData is copied over to both constraint lines, new triangles
 *
 * @warning
 * The z-values of the input vtkPolyData and the points defining the loops are
 * assumed to lie at z=constant. In other words, this filter assumes that the data lies
 * in a plane orthogonal to the z axis.
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
  void PrintSelf(ostream& os, vtkIndent indent) override;

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
   * Specify a subclass of vtkAbstractCellLocator which implements the method 'FindCellsWithinBounds()'.
   * Ex: vtkStaticCellLocator, vtkCellLocator. Not vtkOBBTree
   */
  vtkSetSmartPointerMacro(CellLocator, vtkAbstractCellLocator);
  vtkGetSmartPointerMacro(CellLocator, vtkAbstractCellLocator);
  //@}

  //@{
  /**
   * Specify a spatial point locator for merging points. By default, an
   * instance of vtkMergePoints is used.
   */
  vtkSetSmartPointerMacro(PointLocator, vtkIncrementalPointLocator);
  vtkGetSmartPointerMacro(PointLocator, vtkIncrementalPointLocator);
  //@}

  /**
   * Specify the a second vtkPolyData input which defines loops used to cut
   * the input polygonal data. These loops must be manifold, i.e., do not
   * self intersect. The loops are defined from the polygons defined in
   * this second input.
   */
  void SetLoops(vtkPointSet* loops);

  /**
   * Specify the a second vtkPolyData input which defines loops used to cut
   * the input polygonal data. These loops must be manifold, i.e., do not
   * self intersect. The loops are defined from the polygons defined in
   * this second input.
   */
  void SetLoopsConnection(vtkAlgorithmOutput* output);
  
  /**
   * Create default locators. Used to create one when none are specified. 
   * The point locator is used to merge coincident points.
   * The cell locator is used to accelerate cell searches. 
   */
  void CreateDefaultLocators();

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
