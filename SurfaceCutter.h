#ifndef SurfaceCutter_h__
#define SurfaceCutter_h__

#include <vtkPointSetAlgorithm.h>

class SurfaceCutter : public vtkPointSetAlgorithm
{
public:
  static SurfaceCutter* New();
  vtkTypeMacro(SurfaceCutter, vtkPointSetAlgorithm);

  vtkBooleanMacro(ColorAcquiredPts, bool);
  vtkSetMacro(ColorAcquiredPts, bool);
  vtkGetMacro(ColorAcquiredPts, bool);

  vtkBooleanMacro(ColorLoopEdges, bool);
  vtkSetMacro(ColorLoopEdges, bool);
  vtkGetMacro(ColorLoopEdges, bool);

  vtkBooleanMacro(InsideOut, bool); // default: remove portions outside loop polygons.
  vtkSetMacro(InsideOut, bool);
  vtkGetMacro(InsideOut, bool);

  void SetLoops(vtkPointSet* loops);
  void SetLoopsConnection(vtkAlgorithmOutput* output);

protected:
  SurfaceCutter();
  ~SurfaceCutter();

  bool ColorAcquiredPts;
  bool ColorLoopEdges;
  bool InsideOut;

  int FillInputPortInformation(int port, vtkInformation* info) override;
  int FillOutputPortInformation(int port, vtkInformation* info) override;
  int RequestDataObject(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;
  int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

private:
  SurfaceCutter(const SurfaceCutter&) = delete;
  void operator=(const SurfaceCutter&) = delete;
};

#endif // SurfaceCutter_h__
