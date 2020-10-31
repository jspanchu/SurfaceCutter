#include "SurfaceCutter.h"

#include <vtkAOSDataArrayTemplate.h>
#include <vtkArrayDispatch.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkXMLPolyDataWriter.h>
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
  this->ComputeProjectedLoop = true;
  this->InsideOut = true; // default: remove portions outside loop polygons.
  this->TagAcquiredPoints = true;
  this->TagIntersections = true;

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);

  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::AttributeTypes::SCALARS);

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

int SurfaceCutter::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int SurfaceCutter::RequestDataObject(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  for (int i = 0; i < this->GetNumberOfOutputPorts(); ++i)
  {
    vtkInformation* outInfo = outputVector->GetInformationObject(i);
    vtkPolyData* output = dynamic_cast<vtkPolyData*>(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));
    if (!output)
    {
      output = vtkPolyData::New();
      outInfo->Set(vtkDataObject::DATA_OBJECT(), output);
      output->FastDelete();
      this->GetOutputPortInformation(i)->Set(
        vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType());
    }
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
  vtkSmartPointer<vtkPolyData> output = vtkPolyData::GetData(outputVector->GetInformationObject(0));

  if (!loopsIn->GetNumberOfPoints())
  {
    return 1;
  }

  auto loops = vtkSmartPointer<vtkPolyData>::New();
  loops->SetPoints(loopsIn->GetPoints());
  loops->SetPolys(loopsIn->GetPolys());
  loops->GetCellData()->PassData(loopsIn->GetCellData());

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
  int assoc = vtkDataObject::FIELD_ASSOCIATION_POINTS;
  if ((scalars = this->GetInputArrayToProcess(0, 0, inputVector, assoc)) == nullptr)
  {
    vtkDebugMacro(<< "Input surface is missing scalars. Will add dummy scalars");
    auto _dummy = vtkSmartPointer<vtkAOSDataArrayTemplate<float>>::New();
    _dummy->SetNumberOfComponents(1);
    _dummy->SetNumberOfTuples(input->GetNumberOfPoints());
    _dummy->SetName("DummyScalars");
    _dummy->FillValue(0.0);
    input->GetPointData()->AddArray(_dummy);
    scalars = _dummy;
    dummyAdded = true;
  }

  using dispatcher = vtkArrayDispatch::Dispatch3ByValueType<vtkArrayDispatch::Reals, vtkArrayDispatch::Reals, vtkArrayDispatch::Reals>;

  if (input->IsA("vtkPolyData"))
  {
    vtkSmartPointer<vtkPolyData> inMesh = vtkPolyData::SafeDownCast(input);
    vtkDataArray* points = inMesh->GetPoints()->GetData();
    auto worker = SurfCutterImpl(inMesh, loops, this->ComputeBoolean2D, this->ComputeProjectedLoop);
    if (!dispatcher::Execute(points, scalars, loopsPts, worker))
    {
      worker(points, scalars, loopsPts);
    }
    output->ShallowCopy(worker.outMesh);
  }  
  else if (input->IsA("vtkUnstructuredGrid"))
  {
    vtkSmartPointer<vtkUnstructuredGrid> inMesh = vtkUnstructuredGrid::SafeDownCast(input);
    vtkDataArray* points = inMesh->GetPoints()->GetData();
    auto worker = SurfCutterImpl(inMesh, loops, this->ComputeBoolean2D, this->ComputeProjectedLoop);
    if (!dispatcher::Execute(points, scalars, loopsPts, worker))
    {
      worker(points, scalars, loopsPts);
    }
    output->ShallowCopy(worker.outMesh);
  }

  if (dummyAdded)
  {
    input->GetPointData()->RemoveArray("DummyScalars");
    output->GetPointData()->RemoveArray("DummyScalars");
  }

  if (!this->TagAcquiredPoints)
    output->GetPointData()->RemoveArray("Acquired");

  if (!this->TagIntersections)
    output->GetPointData()->RemoveArray("Intersected");

  return 1;
}