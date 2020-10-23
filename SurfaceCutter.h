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

  vtkBooleanMacro(TagInsertedPoints, bool);
  vtkSetMacro(TagInsertedPoints, bool);
  vtkGetMacro(TagInsertedPoints, bool);
  
  vtkBooleanMacro(InsideOut, bool);
  vtkSetMacro(InsideOut, bool);
  vtkGetMacro(InsideOut, bool);

protected:
  SurfaceCutter();
  ~SurfaceCutter();

  bool ComputeBoolean2D;
  bool TagInsertedPoints;
  bool InsideOut;

  int RequestData(vtkInformation* request, vtkInformationVector** inInfo, vtkInformationVector* ouInfo) override;

private:
  SurfaceCutter(const SurfaceCutter&) = delete;
  void operator=(const SurfaceCutter&) = delete;
};

#endif // SurfaceCutter_h__