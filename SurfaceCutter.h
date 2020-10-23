#ifndef SurfaceCutter_h__
#define SurfaceCutter_h__

#include <vtkDataSetAlgorithm.h>
#include <vtkSetGet.h>

class SurfaceCutter : public vtkDataSetAlgorithm {
public:
  static SurfaceCutter* New();
  vtkTypeMacro(SurfaceCutter, vtkDataSetAlgorithm);

  vtkBooleanMacro(ComputeBoolean2D, bool);
  vtkSetMacro(ComputeBoolean2D, bool);
  vtkGetMacro(ComputeBoolean2D, bool);

  vtkBooleanMacro(TagAcquiredPoints, bool);
  vtkSetMacro(TagAcquiredPoints, bool);
  vtkGetMacro(TagAcquiredPoints, bool);
  
  vtkBooleanMacro(InsideOut, bool);
  vtkSetMacro(InsideOut, bool);
  vtkGetMacro(InsideOut, bool);

  void SetLoops(vtkDataSet* loops);

protected:
  SurfaceCutter();
  ~SurfaceCutter();

  bool ComputeBoolean2D;
  bool TagAcquiredPoints;
  bool InsideOut;

  int FillInputPortInformation(int port, vtkInformation* info) override;
  int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

private:
  SurfaceCutter(const SurfaceCutter&) = delete;
  void operator=(const SurfaceCutter&) = delete;
};

#endif // SurfaceCutter_h__