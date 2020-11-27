#include <SurfaceCutter.h>
#include <vtkActor.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkClipDataSet.h>
#include <vtkCommand.h>
#include <vtkCookieCutter.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkImplicitSelectionLoop.h>
#include <vtkInteractorStyleImage.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkThreshold.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <vtkStaticCellLocator.h>

#include <array>
#include <chrono>
#include <iostream>

static double translateSpeed = 0.1;
static double rotateSpeed = 0.5;

static void ShowUsage(const char* appname)
{
  std::cerr << "Usage: " << appname << " <option(s)>"
            << "Options:\n"
            << "\t-h,--help\tShow this help message\n"
            << "\t-m,--mesh \tSpecify mesh file (*.vtp, *.vtu)\n"
            << "\t-l,--loops \tSpecify loops file (*.vtp)\n"
            << "\t-i,--invert \tInvert 2d boolean. Portions inside loops will be "
               "removed.\n"
            << "\t-t,--translationspeed \tSpeed multiplier for mesh translations along "
               "x, y, z\n"
            << "\t-r,--rotationspeed \tSpeed multiplier for mesh rotation along z\n"
            << "\t   --movable \tMake the mesh movable.\n"
            << "\t   --vtkcookiecutter \tUse vtkCookieCutter instead\n"
            << "\t   --vtkclipdataset \tUse vtkClipDataset with "
               "vtkImplicitSelectionLoop instead\n"
            << "\n"
            << std::endl;
}

static void ShowControls()
{
  std::cout << "Controls:\n"
            << "W:    Z+ | S:     Z-\n"
            << "Up:   Y+ | Down:  Y-\n"
            << "Left: X+ | Right: X-\n"
            << "Z:   CCW | C:     CW (Looking down Z-)\n"
            << std::endl;
}

static bool has_suffix(const std::string& str, const std::string& suffix)
{
  return str.size() >= suffix.size() &&
    str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static void KeypressCallbackFunction(
  vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
  auto iren = static_cast<vtkRenderWindowInteractor*>(caller);
  auto meshTransform = static_cast<vtkTransform*>(clientData);

  std::string key = iren->GetKeySym();

  if (key == "Up")
    meshTransform->Translate(0., translateSpeed, 0.);
  else if (key == "Down")
    meshTransform->Translate(0., -translateSpeed, 0.);
  else if (key == "Left")
    meshTransform->Translate(-translateSpeed, 0., 0.);
  else if (key == "Right")
    meshTransform->Translate(translateSpeed, 0., 0.);
  else if (key == "z")
    meshTransform->RotateZ(rotateSpeed);
  else if (key == "c")
    meshTransform->RotateZ(-rotateSpeed);
  else if (key == "k")
    meshTransform->Scale(1.5, 1.5, 1.0);
  else if (key == "l")
    meshTransform->Scale(0.5, 0.5, 1.0);
  else if (key == "h")
  {
    ShowControls();
    return;
  }

  std::cout << "Recomputing .. ";
  auto start_time = std::chrono::high_resolution_clock::now();
  iren->GetRenderWindow()->Render();
  auto stop_time = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds =
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);

  std::cout << "Elapsed : " << elapsed_seconds.count() << "ms\n";
}

int Parse(char** argv, std::string& meshFile, std::string& loopsFile, bool& insideOut,
  bool& useCookieCutter, bool& useClipDataSet, bool& movable, int argc)
{
  int arg = 0;
  do
  {
    if ((std::string(argv[arg]) == "-h") || (std::string(argv[arg]) == "--help"))
    {
      ShowUsage(argv[0]);
      ShowControls();
      return 1;
    }
    else if ((std::string(argv[arg]) == "-m") || (std::string(argv[arg]) == "--mesh"))
      meshFile = std::string(argv[++arg]);
    else if ((std::string(argv[arg]) == "-l") || (std::string(argv[arg]) == "--loops"))
      loopsFile = std::string(argv[++arg]);
    else if ((std::string(argv[arg]) == "-i") || (std::string(argv[arg]) == "--invert"))
      insideOut = false;
    else if ((std::string(argv[arg]) == "-t") || (std::string(argv[arg]) == "--translationspeed"))
      translateSpeed = std::stod(argv[++arg]);
    else if ((std::string(argv[arg]) == "-r") || (std::string(argv[arg]) == "--rotationspeed"))
      rotateSpeed = std::stod(argv[++arg]);
    else if (std::string(argv[arg]) == "--vtkcookiecutter")
      useCookieCutter = true;
    else if (std::string(argv[arg]) == "--vtkclipdataset")
      useClipDataSet = true;
    else if (std::string(argv[arg]) == "--movable")
      movable = false;
    ++arg;
  } while (arg < argc);

  return 0;
}

int main(int argc, char** argv)
{
  bool useCookieCutter(false), movable(true), useClipDataSet(false), insideOut(true);
  std::string meshFile = "data/Surface.vtp";
  std::string loopsFile = "data/Case5.vtp";

  const int invalidArgs =
    Parse(argv, meshFile, loopsFile, insideOut, useCookieCutter, useClipDataSet, movable, argc);

  if (invalidArgs)
    return EXIT_FAILURE;

  vtkNew<vtkNamedColors> colors;

  vtkSmartPointer<vtkXMLUnstructuredDataReader> reader;
  vtkNew<vtkThreshold> constraints; // highlight them.

  if (has_suffix(meshFile, ".vtp"))
  {
    reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(meshFile.c_str());
  }
  else if (has_suffix(meshFile, ".vtu"))
  {
    reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(meshFile.c_str());
  }
  else
  {
    std::cerr << "Unsupported mesh file extension " << meshFile << "\n";
    return EXIT_FAILURE;
  }

  if (!has_suffix(loopsFile, ".vtp"))
  {
    std::cerr << "Unsupported loops file extension " << loopsFile << "\n";
    return EXIT_FAILURE;
  }

  vtkNew<vtkXMLPolyDataReader> loopsReader;
  loopsReader->SetFileName(loopsFile.c_str());

  vtkNew<vtkTransform> meshTransform;
  vtkNew<vtkTransformFilter> meshTransformer;
  meshTransform->PostMultiply();
  meshTransformer->SetTransform(meshTransform);
  if (has_suffix(meshFile, ".vtu"))
  {
    vtkNew<vtkDataSetSurfaceFilter> getSurf;
    getSurf->SetInputConnection(reader->GetOutputPort());
    meshTransformer->SetInputConnection(0, getSurf->GetOutputPort());
  }
  else
    meshTransformer->SetInputConnection(reader->GetOutputPort());

  vtkSmartPointer<vtkAlgorithm> surfCutter;
  if (!(useCookieCutter || useClipDataSet))
  {
    vtkNew<SurfaceCutter> surfCutter_;
    vtkNew<vtkStaticCellLocator> cellLoc;
    surfCutter_->SetCellLocator(cellLoc);
    surfCutter_->SetInputConnection(0, meshTransformer->GetOutputPort());
    surfCutter_->SetInputConnection(1, loopsReader->GetOutputPort());
    surfCutter_->SetInsideOut(insideOut);
    surfCutter_->Print(std::cout);
    surfCutter = vtkAlgorithm::SafeDownCast(surfCutter_);
    constraints->SetInputConnection(surfCutter->GetOutputPort());
    constraints->SetInputArrayToProcess(
      0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Constrained");
  }
  else if (useCookieCutter)
  {
    surfCutter = vtkSmartPointer<vtkCookieCutter>::New();
    surfCutter->SetInputConnection(0, meshTransformer->GetOutputPort());
    surfCutter->SetInputConnection(1, loopsReader->GetOutputPort());
  }
  else if (useClipDataSet)
  {
    vtkNew<vtkClipDataSet> surfCutter_;
    surfCutter_->SetInputConnection(0, meshTransformer->GetOutputPort());

    vtkNew<vtkImplicitSelectionLoop> clipFunc;
    loopsReader->Update();
    clipFunc->SetLoop(loopsReader->GetOutput()->GetPoints());
    surfCutter_->SetInsideOut(insideOut);
    surfCutter_->SetClipFunction(clipFunc);
    surfCutter = vtkAlgorithm::SafeDownCast(surfCutter_);
  }

  vtkNew<vtkDataSetMapper> meshMapper;
  meshMapper->SetInputConnection(meshTransformer->GetOutputPort(0));
  meshMapper->ScalarVisibilityOff();

  vtkNew<vtkDataSetMapper> cutMeshMapper;
  cutMeshMapper->SetInputConnection(surfCutter->GetOutputPort());
  cutMeshMapper->ScalarVisibilityOff();

  vtkNew<vtkDataSetMapper> constraintsMapper;
  constraintsMapper->SetInputConnection(constraints->GetOutputPort());
  constraintsMapper->ScalarVisibilityOff();

  vtkNew<vtkDataSetMapper> polysMapper;
  polysMapper->SetInputConnection(loopsReader->GetOutputPort());

  vtkNew<vtkActor> meshActor;
  meshActor->SetMapper(meshMapper);
  meshActor->GetProperty()->SetRepresentationToSurface();
  meshActor->GetProperty()->SetOpacity(0.1);
  meshActor->GetProperty()->SetColor(colors->GetColor3d("ivory_black").GetData());
  
  vtkNew<vtkActor> cutMeshActor;
  cutMeshActor->SetMapper(cutMeshMapper);
  cutMeshActor->GetProperty()->SetRepresentationToSurface();
  cutMeshActor->GetProperty()->SetEdgeVisibility(true);
  cutMeshActor->GetProperty()->SetEdgeColor(colors->GetColor3d("indigo").GetData());
  cutMeshActor->GetProperty()->SetColor(colors->GetColor3d("white").GetData());

  vtkNew<vtkActor> constraintsActor;
  constraintsActor->SetMapper(constraintsMapper);
  constraintsActor->GetProperty()->SetRepresentationToSurface();
  constraintsActor->GetProperty()->SetLineWidth(4);
  constraintsActor->GetProperty()->SetColor(colors->GetColor3d("cadmium_red_deep").GetData());

  vtkNew<vtkActor> polysActor;
  polysActor->SetMapper(polysMapper);
  polysActor->GetProperty()->SetRepresentationToWireframe();
  polysActor->GetProperty()->SetLineWidth(1);
  polysActor->GetProperty()->SetColor(colors->GetColor3d("white_smoke").GetData());

  vtkNew<vtkRenderer> renderer;
  renderer->AddActor(cutMeshActor);
  renderer->AddActor(meshActor);
  renderer->AddActor(polysActor);
  renderer->AddActor(constraintsActor);
  renderer->SetBackground(colors->GetColor3d("slate_grey").GetData());

  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->SetSize(640, 480);
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("Example: SurfaceCutter");

  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);

  auto start_time = std::chrono::high_resolution_clock::now();
  renderer->GetActiveCamera()->SetParallelProjection(true);
  renderer->ResetCameraClippingRange();
  renderer->ResetCamera();
  renderWindow->Render();
  auto stop_time = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds =
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);

  vtkSmartPointer<vtkPointSet> mesh, loops;
  mesh = vtkPointSet::SafeDownCast(meshTransformer->GetInputDataObject(0, 0));
  loops = vtkPointSet::SafeDownCast(loopsReader->GetOutput());
  std::cout << "Mesh:- "
            << "\n";
  std::cout << " Cells : " << mesh->GetNumberOfCells() << "\n";
  std::cout << " Points: " << mesh->GetNumberOfPoints() << "\n";
  std::cout << "Loops:- "
            << "\n";
  std::cout << " Cells : " << loops->GetNumberOfCells() << "\n";
  std::cout << " Points: " << loops->GetNumberOfPoints() << "\n";
  std::cout << "Elapsed : " << elapsed_seconds.count() << "ms\n";

  
  vtkNew<vtkXMLPolyDataWriter> writer;
  writer->SetFileName("CutMesh.vtp");
  writer->SetInputConnection(surfCutter->GetOutputPort(0));
  writer->Write();

  vtkNew<vtkInteractorStyleImage> istyle;
  renderWindowInteractor->SetInteractorStyle(istyle);
  renderWindowInteractor->Initialize();
  renderWindowInteractor->CreateRepeatingTimer(1);

  if (movable)
  {
    vtkNew<vtkCallbackCommand> keypressCallback;
    keypressCallback->SetClientData(meshTransform);
    keypressCallback->SetCallback(KeypressCallbackFunction);
    renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);
  }
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}