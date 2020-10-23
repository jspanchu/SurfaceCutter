#include "SurfaceCutter.h"
#include <vtkObjectFactory.h>

vtkStandardNewMacro(SurfaceCutter);

SurfaceCutter::SurfaceCutter()
{
  this->ComputeBoolean2D = true;
  this->InsideOut = false; // remove portion inside polygons.
  this->TagInsertedPoints = true;
  vtkDebugMacro(<< "Initialized " << this->GetClassNameInternal());
}

SurfaceCutter::~SurfaceCutter()
{
  vtkDebugMacro(<< "Destroyed " << this->GetClassNameInternal());
}

int SurfaceCutter::RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*)
{

  return 1;
}
