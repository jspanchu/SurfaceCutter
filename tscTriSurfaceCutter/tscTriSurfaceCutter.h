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
 * with cutters specified by a second input containing polygons.
 *
 * The cutter polygons can be concave, can have vertices exactly
 * coincident with a mesh point/edge.
 *
 * It computes an imprint of the cutter polygons' edges upon the mesh
 * followed by **removal** of triangles *in(out)side the polygons. See SetInsideOut().
 *
 * Linear cells other than triangles will be passed through.
 * Line segments and polylines from input will be marked as constraints.
 *
 * It is also possible to output a pure imprint or a pure removal.
 *
 * @note:
 * Input point-data is interpolated to output.
 * Input cell-data is copied to output.
 * The cutter is projected onto the surface prior to any removal of cells.
 *
 * @sa
 * vtkClipDataSet vtkClipPolyData vtkImprintFilter
 *
 */

#include <tscTriSurfaceCutterModule.h>

#include "vtkIncrementalPointLocator.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkPolyData;

class TSCTRISURFACECUTTER_EXPORT tscTriSurfaceCutter : public vtkPolyDataAlgorithm
{
public:
  /**
   * Construct object with tolerance 1.0e-6, inside out set to true,
   * color acquired points, color cutter edges
   */
  static tscTriSurfaceCutter* New();
  vtkTypeMacro(tscTriSurfaceCutter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //@{
  /**
   * After the cutter's edges are embedded onto the surface,
   * On: remove stuff outside all cutter polygons
   * Off: remove stuff inside atleast one cutter polygon
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
   * Specify a spatial point locator for merging points. By default, an
   * instance of vtkMergePoints is used.
   */
  vtkSetSmartPointerMacro(PointLocator, vtkIncrementalPointLocator);
  vtkGetSmartPointerMacro(PointLocator, vtkIncrementalPointLocator);
  //@}

  //@{
  /**
   * Do not respect the very functionality of this filter. If enabled, it will only imprint cutter
   *    polygons onto the mesh
   * @note InsideOut option does not apply here.
   */
  vtkBooleanMacro(Imprint, bool);
  vtkSetMacro(Imprint, bool);
  vtkGetMacro(Imprint, bool);
  //@}

  //@{
  /**
   * Partially respect functionality of this filter. Only remove cells
   * in(out)side cutter polygons.
   */
  vtkBooleanMacro(Remove, bool);
  vtkSetMacro(Remove, bool);
  vtkGetMacro(Remove, bool);
  //@}

  /**
   * Specify the a second vtkPolyData input which defines cutters used to cut
   * the input polygonal data. These cutters must be manifold, i.e., do not
   * self intersect. The cutters are defined from the polygons defined in
   * this second input.
   */
  void SetCuttersData(vtkPolyData* cutters);

  /**
   * Specify the a second vtkPolyData input which defines cutters used to cut
   * the input polygonal data. These cutters must be manifold, i.e., do not
   * self intersect. The cutters are defined from the polygons defined in
   * this second input.
   */
  void SetCuttersConnection(vtkAlgorithmOutput* output);

  /**
   * Create default locators. Used to create one when none are specified.
   * The point locator is used to merge coincident points.
   * The cell locator is used to accelerate cell searches.
   */
  void CreateDefaultLocators();

protected:
  tscTriSurfaceCutter();
  ~tscTriSurfaceCutter() override;

  bool Imprint = true;
  bool InsideOut = true; // default: remove portions outside cutter polygons.
  bool Remove = true;
  double Tolerance = 1.0e-6;

  vtkSmartPointer<vtkIncrementalPointLocator> PointLocator;

  int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;
  int FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info) override;
  int FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info) override;

private:
  tscTriSurfaceCutter(const tscTriSurfaceCutter&) = delete;
  void operator=(const tscTriSurfaceCutter&) = delete;
};
