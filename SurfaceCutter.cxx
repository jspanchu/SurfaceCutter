#include "SurfaceCutter.h"

#include <vtkAOSDataArrayTemplate.h>
#include <vtkArrayDispatch.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "SurfaceCutterHelper.h"

vtkStandardNewMacro(SurfaceCutter);

SurfaceCutter::SurfaceCutter()
{
  this->ComputeBoolean2D = true;
  this->InsideOut = true; // remove portion outside polygons.
  this->TagAcquiredPoints = true;

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);

  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::SCALARS);

  vtkDebugMacro(<< "Initialized " << this->GetClassNameInternal());
}

SurfaceCutter::~SurfaceCutter()
{
  vtkDebugMacro(<< "Destroyed " << this->GetClassNameInternal());
}

int SurfaceCutter::FillInputPortInformation(int port, vtkInformation* info)
{
  switch (port)
  {
  case 0:
    break;
  case 1:
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    break;
  default:
    break;
  }
  return 1;
}

void SurfaceCutter::SetLoops(vtkDataSet* loops)
{
  this->SetInputData(1, loops);
}

int SurfaceCutter::RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkSmartPointer<vtkDataSet> input = vtkDataSet::GetData(inputVector[0]->GetInformationObject(0));
  vtkSmartPointer<vtkPolyData> loopsIn = vtkPolyData::GetData(inputVector[1]->GetInformationObject(0));
  vtkSmartPointer<vtkDataSet> output = vtkDataSet::GetData(outputVector->GetInformationObject(0));

  if (!loopsIn->GetNumberOfPoints())
  {
    return 1;
  }

  output->DeepCopy(input);
  auto loops = vtkSmartPointer<vtkPolyData>::New();
  loops->DeepCopy(loopsIn);

  vtkDataArray* insideOuts;
  if ((insideOuts = loops->GetCellData()->GetArray("InsideOuts")) == nullptr)
  {
    vtkDebugMacro(<< "Loop polygons do not have InsideOuts array. Will resort to " << this->GetClassNameInternal() << "::InsideOut = " << this->InsideOut);
    auto _insideOuts = vtkSmartPointer<vtkAOSDataArrayTemplate<int>>::New();
    _insideOuts->SetNumberOfComponents(1);
    _insideOuts->SetNumberOfTuples(loops->GetNumberOfPolys());
    _insideOuts->SetName("InsideOuts");
    _insideOuts->FillValue(this->InsideOut);
    loops->GetCellData()->AddArray(_insideOuts);
  }

  vtkDataArray* loopsPts = loops->GetPoints()->GetData();
  vtkDataArray* scalars;
  bool dummyAdded = false;
  if ((scalars = this->GetInputArrayToProcess(0, input)) == nullptr)
  {
    vtkDebugMacro(<< "Input surface is missing scalars. Will add dummy scalars");
    auto _dummy = vtkSmartPointer<vtkAOSDataArrayTemplate<float>>::New();
    _dummy->SetNumberOfComponents(1);
    _dummy->SetNumberOfTuples(output->GetNumberOfPoints());
    _dummy->SetName("DummyScalars");
    _dummy->FillValue(0.0);
    output->GetPointData()->AddArray(_dummy);
    scalars->ShallowCopy(_dummy);
    dummyAdded = true;
  }

  using dispatcher = vtkArrayDispatch::Dispatch3ByValueType<vtkArrayDispatch::Reals, vtkArrayDispatch::Reals, vtkArrayDispatch::Reals>;
  if (output->IsA("vtkPolyData"))
  {
    vtkSmartPointer<vtkPolyData> meshPd = vtkPolyData::SafeDownCast(output);
    vtkDataArray* points = meshPd->GetPoints()->GetData();

    auto worker = SurfCutterImpl(meshPd, loops, scalars->GetName(), this->InsideOut, this->ComputeBoolean2D);
    if (!dispatcher::Execute(points, scalars, loopsPts, worker))
    {
      worker(points, scalars, loopsPts);
    }
  }  
  else if (output->IsA("vtkUnstructuredGrid"))
  {
    vtkSmartPointer<vtkUnstructuredGrid> meshUgrid = vtkUnstructuredGrid::SafeDownCast(output);
    vtkDataArray* points = meshUgrid->GetPoints()->GetData();
    auto worker = SurfCutterImpl(meshUgrid, loops, scalars->GetName(), this->InsideOut, this->ComputeBoolean2D);
    if (!dispatcher::Execute(points, scalars, loopsPts, worker))
    {
      worker(points, scalars, loopsPts);
    }
  }

  if (dummyAdded)
    output->GetPointData()->RemoveArray("DummyScalars");

  if (!this->TagAcquiredPoints)
    output->GetPointData()->RemoveArray("Acquired");

  return 1;
}