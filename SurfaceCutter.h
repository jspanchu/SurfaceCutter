#ifndef SurfaceCutter_h__
#define SurfaceCutter_h__

#include <vtkCellArray.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSetGet.h>

class SurfaceCutter : public vtkDataSetAlgorithm {
public:
  static SurfaceCutter* New();
  vtkTypeMacro(SurfaceCutter, vtkDataSetAlgorithm);

  vtkBooleanMacro(ComputeBoolean2D, bool);
  vtkSetMacro(ComputeBoolean2D, bool);
  vtkGetMacro(ComputeBoolean2D, bool);
  
  vtkBooleanMacro(InsideOut, bool);
  vtkSetMacro(InsideOut, bool);
  vtkGetMacro(InsideOut, bool);

  vtkBooleanMacro(TagAcquiredEdges, bool);
  vtkSetMacro(TagAcquiredEdges, bool);
  vtkGetMacro(TagAcquiredEdges, bool);

  vtkBooleanMacro(TagAcquiredPoints, bool);
  vtkSetMacro(TagAcquiredPoints, bool);
  vtkGetMacro(TagAcquiredPoints, bool);

  void SetLoops(vtkDataSet* loops);

protected:
  SurfaceCutter();
  ~SurfaceCutter();

  bool ComputeBoolean2D;
  bool InsideOut;
  bool TagAcquiredEdges;
  bool TagAcquiredPoints;

  int FillInputPortInformation(int port, vtkInformation* info) override;
  int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

  void AcquirePoints(vtkDataSet* mesh, vtkPolyData* loops);
  void AcquireEdges(vtkDataSet* mesh, vtkPoints* points, vtkCellArray* edges);
  void ApplyBoolean(vtkDataSet* mesh, vtkPolyData* loops);

private:
  SurfaceCutter(const SurfaceCutter&) = delete;
  void operator=(const SurfaceCutter&) = delete;
};

#endif // SurfaceCutter_h__
