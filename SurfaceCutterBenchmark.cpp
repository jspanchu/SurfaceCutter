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
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <array>
#include <chrono>
#include <iostream>

static double translateSpeed = 0.05;
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
  else if (key == "w")
    meshTransform->Translate(0., 0., translateSpeed);
  else if (key == "s")
    meshTransform->Translate(0., 0., -translateSpeed);
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
  std::string meshFile = "data/Triangle.vtp";
  std::string loopsFile = "data/Case5.vtp";

  const int invalidArgs =
    Parse(argv, meshFile, loopsFile, insideOut, useCookieCutter, useClipDataSet, movable, argc);

  if (invalidArgs)
    return EXIT_FAILURE;

  vtkSmartPointer<vtkXMLUnstructuredDataReader> reader;

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
  meshTransformer->SetInputConnection(reader->GetOutputPort());

  vtkSmartPointer<vtkAlgorithm> surfCutter;
  if (!(useCookieCutter || useClipDataSet))
  {
    vtkNew<SurfaceCutter> surfCutter_;
    surfCutter_->SetInputConnection(0, meshTransformer->GetOutputPort());
    surfCutter_->SetInputConnection(1, loopsReader->GetOutputPort());
    surfCutter_->SetInsideOut(insideOut);
    surfCutter = vtkAlgorithm::SafeDownCast(surfCutter_);
  }
  else if (useCookieCutter)
  {
    surfCutter = vtkSmartPointer<vtkCookieCutter>::New();
    if (has_suffix(meshFile, ".vtu"))
    {
      vtkNew<vtkDataSetSurfaceFilter> getSurf;
      getSurf->SetInputConnection(reader->GetOutputPort());
      surfCutter->SetInputConnection(0, getSurf->GetOutputPort());
    }
    else
    {
      surfCutter->SetInputConnection(0, meshTransformer->GetOutputPort());
    }
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
  meshMapper->SetInputConnection(surfCutter->GetOutputPort(0));
  meshMapper->SetScalarModeToDefault();
  meshMapper->ScalarVisibilityOn();

  vtkNew<vtkDataSetMapper> projLoopsMapper;
  if (!(useCookieCutter || useClipDataSet))
  {
    projLoopsMapper->SetInputConnection(surfCutter->GetOutputPort(1));
  }
  else
  {
    projLoopsMapper->SetInputConnection(loopsReader->GetOutputPort(0));
  }
  projLoopsMapper->SetScalarModeToDefault();
  projLoopsMapper->ScalarVisibilityOn();

  vtkNew<vtkDataSetMapper> polysMapper;
  polysMapper->SetInputConnection(loopsReader->GetOutputPort());

  vtkNew<vtkActor> meshActor;
  meshActor->GetProperty()->SetRepresentationToSurface();
  meshActor->GetProperty()->EdgeVisibilityOn();
  meshActor->SetMapper(meshMapper);

  vtkNew<vtkActor> projLoopsActor;
  projLoopsActor->GetProperty()->SetRepresentationToWireframe();
  projLoopsActor->SetMapper(projLoopsMapper);
  projLoopsActor->GetProperty()->SetLineWidth(4);

  vtkNew<vtkActor> polysActor;
  polysActor->SetMapper(polysMapper);
  polysActor->GetProperty()->SetRepresentationToWireframe();
  polysActor->GetProperty()->SetLineWidth(1);

  vtkNew<vtkRenderer> renderer;
  renderer->AddActor(meshActor);
  renderer->AddActor(projLoopsActor);
  renderer->AddActor(polysActor);
  renderer->SetBackground(0.39, 0.39, 0.39);

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

  vtkNew<vtkInteractorStyleTrackballCamera> istyle;
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