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

#pragma once
/**
 *
 * @file  tscTriSurfaceCutter.h
 * @class tscTriSurfaceCutter
 * @brief Cut a triangulated surface with one or more polygons.
 *
 * This filter is geometrically based, unlike vtkClipDataSet and vtkClipPolyData
 * (both of which are scalar-based).
 *
 * It crops an input vtkPolyData consisting of triangles
 * with loops specified by a second input containing polygons.
 *
 * The loop polygons can be concave, can have vertices exactly
 * coincident with a mesh point/edge.
 *
 * It computes an **embedding** of the loop polygons' edges upon the mesh
 * followed by **removal** of triangles *in(out)side the polygons. See SetInsideOut().
 *
 * Linear cells other than triangles will be passed through.
 * Line segments and polylines from input will be marked as constraints.
 *
 * It is possible to output a pure embedding or a pure removal.
 *
 * @note:
 * Input point-data is interpolated to output.
 * Input cell-data is copied to output.
 *
 * @sa
 * vtkClipDataSet vtkClipPolyData
 *
 */

#include <tscTriSurfaceCutterModule.h>

#include "vtkAbstractCellLocator.h"
#include "vtkGenericCell.h"
#include "vtkIncrementalPointLocator.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkPolyData;

class TSCTRISURFACECUTTER_EXPORT tscTriSurfaceCutter : public vtkPolyDataAlgorithm
{
public:
  /**
   * Construct object with tolerance 1.0e-6, inside out set to true,
   * color acquired points, color loop edges
   */
  static tscTriSurfaceCutter* New();
  vtkTypeMacro(tscTriSurfaceCutter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //@{
  /**
   * Accelerate cell searches with a multi-threaded cell locator. Default: On
   */
  vtkBooleanMacro(AccelerateCellLocator, bool);
  vtkSetMacro(AccelerateCellLocator, bool);
  vtkGetMacro(AccelerateCellLocator, bool);
  //@}

  //@{
  /**
   * After the loop's edges are embedded onto the surface,
   * On: remove stuff outside all loop polygons
   * Off: remove stuff inside atleast one loop polygon
   */
  vtkBooleanMacro(InsideOut, bool);
  vtkSetMacro(InsideOut, bool);
  vtkGetMacro(InsideOut, bool);
  //@}

  //@{
  /**
   * Numeric tolerance for point merging, intersection math.
   */
  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);
  //@}

  //@{
  /**
   * Specify a subclass of vtkAbstractCellLocator which implements the method
   * 'FindCellsWithinBounds()'. Ex: vtkStaticCellLocator, vtkCellLocator. Not
   * vtkOBBTree
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

  //@{
  /**
   * Do not respect the very functionality of this filter. Only embed loop
   * polygons onto the mesh
   * @note InsideOut option does not apply here.
   */
  vtkBooleanMacro(Embed, bool);
  vtkSetMacro(Embed, bool);
  vtkGetMacro(Embed, bool);
  //@}

  //@{
  /**
   * Partially respect functionality of this filter. Only remove cells
   * in(out)side loop polygons.
   */
  vtkBooleanMacro(Remove, bool);
  vtkSetMacro(Remove, bool);
  vtkGetMacro(Remove, bool);
  //@}

  /**
   * Specify the a second vtkPolyData input which defines loops used to cut
   * the input polygonal data. These loops must be manifold, i.e., do not
   * self intersect. The loops are defined from the polygons defined in
   * this second input.
   */
  void SetLoopsData(vtkPolyData* loops);

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
  tscTriSurfaceCutter();
  ~tscTriSurfaceCutter() override;

  bool AccelerateCellLocator = true;
  bool Embed = true;
  bool InsideOut = true; // default: remove portions outside loop polygons.
  bool Remove = true;
  double Tolerance = 1.0e-6;

  vtkSmartPointer<vtkAbstractCellLocator> CellLocator;
  vtkSmartPointer<vtkIncrementalPointLocator> PointLocator;

  int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;
  int FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info) override;
  int FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info) override;

private:
  tscTriSurfaceCutter(const tscTriSurfaceCutter&) = delete;
  void operator=(const tscTriSurfaceCutter&) = delete;
};
