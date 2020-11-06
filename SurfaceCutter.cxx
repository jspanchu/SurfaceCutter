#include "SurfaceCutter.h"

#include <vtkAOSDataArrayTemplate.h>
#include <vtkArrayDispatch.h>
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "SurfaceCutterHelper.h"

using dispatchRRR = vtkArrayDispatch::Dispatch3ByValueType<vtkArrayDispatch::Reals,
  vtkArrayDispatch::Reals, vtkArrayDispatch::Reals>;

vtkStandardNewMacro(SurfaceCutter);

SurfaceCutter::SurfaceCutter()
{
  this->ColorAcquiredPts = true;
  this->ColorLoopEdges = true;
  this->InsideOut = true; // default: remove portions outside loop polygons.

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);

  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
    vtkDataSetAttributes::AttributeTypes::SCALARS);

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

int SurfaceCutter::RequestDataObject(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  for (int i = 0; i < this->GetNumberOfOutputPorts(); ++i)
  {
    vtkInformation* outInfo = outputVector->GetInformationObject(i);
    vtkPolyData* output = dynamic_cast<vtkPolyData*>(outInfo->Get(vtkDataObject::DATA_OBJECT()));
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

void SurfaceCutter::SetLoops(vtkPointSet* loops)
{
  this->SetInputData(1, loops);
}

void SurfaceCutter::SetLoopsConnection(vtkAlgorithmOutput* output)
{
  this->SetInputConnection(1, output);
}

int SurfaceCutter::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkSmartPointer<vtkPointSet> inMesh =
    vtkPointSet::GetData(inputVector[0]->GetInformationObject(0));
  vtkSmartPointer<vtkPolyData> inLoops =
    vtkPolyData::GetData(inputVector[1]->GetInformationObject(0));

  vtkSmartPointer<vtkPolyData> outMesh =
    vtkPolyData::GetData(outputVector->GetInformationObject(0));
  vtkSmartPointer<vtkPolyData> outLoops =
    vtkPolyData::GetData(outputVector->GetInformationObject(1));

  if (!inLoops->GetNumberOfPoints())
  {
    return 1;
  }

  auto loops = vtkSmartPointer<vtkPolyData>::New();
  loops->SetPoints(inLoops->GetPoints());
  loops->SetPolys(inLoops->GetPolys());
  loops->GetCellData()->PassData(inLoops->GetCellData());

  vtkDataArray* insideOuts;
  if ((insideOuts = loops->GetCellData()->GetArray("InsideOuts")) == nullptr)
  {
    vtkDebugMacro(<< "Loop polygons do not have InsideOuts array. Will resort to "
                  << this->GetClassNameInternal() << "::InsideOut = " << this->InsideOut);
    vtkNew<vtkAOSDataArrayTemplate<int>> _insideOuts;
    _insideOuts->SetNumberOfComponents(1);
    _insideOuts->SetNumberOfTuples(loops->GetNumberOfPolys());
    _insideOuts->SetName("InsideOuts");
    _insideOuts->FillValue(this->InsideOut);
    loops->GetCellData()->AddArray(_insideOuts);
  }

  vtkDataArray* scalars;
  bool dummyAdded = false; // remove dummy scalars from `inMesh` later.
  int assoc = vtkDataObject::FIELD_ASSOCIATION_POINTS;
  if ((scalars = this->GetInputArrayToProcess(0, 0, inputVector, assoc)) == nullptr)
  {
    vtkDebugMacro(<< "Input surface is missing scalars. Will add dummy scalars");
    vtkNew<vtkAOSDataArrayTemplate<int>> _dummy;
    _dummy->SetNumberOfComponents(1);
    _dummy->SetNumberOfTuples(inMesh->GetNumberOfPoints());
    _dummy->SetName("DummyScalars");
    _dummy->FillValue(0.0);
    inMesh->GetPointData()->AddArray(_dummy);
    scalars = _dummy;
    dummyAdded = true;
  }

  vtkDataArray* points = inMesh->GetPoints()->GetData();
  vtkDataArray* loopsPts = loops->GetPoints()->GetData();
  auto worker = SurfCutterImpl(inMesh, loops);
  if (!dispatchRRR::Execute(points, scalars, loopsPts, worker))
  {
    worker(points, scalars, loopsPts);
  }
  outMesh->ShallowCopy(worker.outMesh);
  outLoops->ShallowCopy(worker.outLoops);
  outLoops->GetPointData()->SetActiveScalars("Intersected");

  if (dummyAdded)
  {
    inMesh->GetPointData()->RemoveArray("DummyScalars");
    outMesh->GetPointData()->RemoveArray("DummyScalars");
  }

  if (!this->ColorAcquiredPts)
  {
    outMesh->GetPointData()->RemoveArray("Acquired");
    outLoops->GetPointData()->RemoveArray("Acquired");
  }

  if (!this->ColorLoopEdges)
  {
    outMesh->GetPointData()->RemoveArray("Intersected");
    outLoops->GetPointData()->RemoveArray("Intersected");
  }

  return 1;
}