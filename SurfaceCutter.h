#ifndef SurfaceCutter_h__
#define SurfaceCutter_h__

#include <vtkDataSetAlgorithm.h>

class SurfaceCutter : public vtkDataSetAlgorithm {
public:
  static SurfaceCutter* New();
  vtkTypeMacro(SurfaceCutter, vtkDataSetAlgorithm);

  vtkBooleanMacro(ComputeBoolean2D, bool);
  vtkSetMacro(ComputeBoolean2D, bool);
  vtkGetMacro(ComputeBoolean2D, bool);
  
  vtkBooleanMacro(InsideOut, bool); // default: remove portions outside loop polygons.
  vtkSetMacro(InsideOut, bool);
  vtkGetMacro(InsideOut, bool);

  vtkBooleanMacro(TagAcquiredPoints, bool);
  vtkSetMacro(TagAcquiredPoints, bool);
  vtkGetMacro(TagAcquiredPoints, bool);

  void SetLoops(vtkDataSet* loops);

protected:
  SurfaceCutter();
  ~SurfaceCutter();

  bool ComputeBoolean2D;
  bool InsideOut;
  bool TagAcquiredPoints;

  int FillInputPortInformation(int port, vtkInformation* info) override;
  int FillOutputPortInformation(int port, vtkInformation* info) override;
  int RequestDataObject(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;
  int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

private:
  SurfaceCutter(const SurfaceCutter&) = delete;
  void operator=(const SurfaceCutter&) = delete;
};

#endif // SurfaceCutter_h__
