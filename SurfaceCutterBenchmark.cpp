#include <vtkActor.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkClipDataSet.h>
#include <vtkCommand.h>
#include <vtkCylinderSource.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkImplicitSelectionLoop.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkCookieCutter.h>
#include <SurfaceCutter.h>

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
    << "\t-i,--invert \tInvert 2d boolean. Portions inside loops will be removed.\n"
    << "\t-t,--translationspeed \tSpeed multiplier for mesh translations along x, y, z\n"
    << "\t-r,--rotationspeed \tSpeed multiplier for mesh rotation along z\n"
    << "\t   --movable \tMake the mesh movable.\n"
    << "\t   --vtkcookiecutter \tUse vtkCookieCutter instead\n"
    << "\t   --vtkclipdataset \tUse vtkClipDataset with vtkImplicitSelectionLoop instead\n"
    << "\n"
    << std::endl;
}

static void ShowControls() {
  std::cout <<
    "Controls:\n"
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

static void KeypressCallbackFunction(vtkObject* caller, long unsigned int eventId,
  void* clientData, void* callData)
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
  auto elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);

  std::cout << "Elapsed : " << elapsed_seconds.count() << "ms\n";
}

int main(int argc, char** argv) {

  bool useCookieCutter(false), movable(true), useClipDataSet(false), insideOut(true);
  std::string meshFile = "data/BigSurface.vtp";
  std::string loopsFile = "data/TestPolys2.vtp";

  int arg = 0;
  do
  {
    if ((std::string(argv[arg]) == "-h") || (std::string(argv[arg]) == "--help"))
    {
      ShowUsage(argv[0]);
      ShowControls();
      return EXIT_FAILURE;
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

  vtkSmartPointer<vtkNamedColors> colors =
    vtkSmartPointer<vtkNamedColors>::New();

  // Set the background color.
  std::array<unsigned char, 4> bkg{ {26, 51, 102, 255} };
  colors->SetColor("BkgColor", bkg.data());

  vtkSmartPointer<vtkPolyData> meshPd;
  vtkSmartPointer<vtkUnstructuredGrid> meshUgrid;

  if (has_suffix(meshFile, ".vtp"))
  {
    auto meshReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    meshReader->SetFileName(meshFile.c_str());
    meshReader->Update();
    meshPd = meshReader->GetOutput();
  }
  else if (has_suffix(meshFile, ".vtu"))
  {
    auto meshReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    meshReader->SetFileName(meshFile.c_str());
    meshReader->Update();
    meshUgrid = meshReader->GetOutput();
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

  auto meshTransform = vtkSmartPointer<vtkTransform>::New();
  auto meshTransformFilter = vtkSmartPointer<vtkTransformFilter>::New();
  meshTransformFilter->SetTransform(meshTransform);
  //meshTransform->RotateZ(233.5);
  if (meshPd)
    meshTransformFilter->SetInputData(meshPd);
  else if (meshUgrid)
  {
    auto surf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surf->SetInputData(meshUgrid);
    surf->Update();
    meshTransformFilter->SetInputData(surf->GetOutput());
  }

  auto loopsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  loopsReader->SetFileName(loopsFile.c_str());
  
  auto surfCutter = vtkSmartPointer<vtkAlgorithm>::New();
  if (!(useCookieCutter || useClipDataSet))
  {
    auto surfCutter_ = vtkSmartPointer<SurfaceCutter>::New();
    surfCutter_->SetInputConnection(0, meshTransformFilter->GetOutputPort());
    surfCutter_->SetInputConnection(1, loopsReader->GetOutputPort());
    surfCutter_->SetInsideOut(insideOut);
    surfCutter = vtkAlgorithm::SafeDownCast(surfCutter_);
  }
  else if (useCookieCutter)
  {
    // inside out unavailable yet.
    surfCutter = vtkSmartPointer<vtkCookieCutter>::New();
    surfCutter->SetInputConnection(0, meshTransformFilter->GetOutputPort());
    surfCutter->SetInputConnection(1, loopsReader->GetOutputPort());
    surfCutter->GetInputDataObject(1, 0);
  }
  else if (useClipDataSet)
  {
    auto surfCutter_ = vtkSmartPointer<vtkClipDataSet>::New();
    surfCutter_->SetInputConnection(0, meshTransformFilter->GetOutputPort());
    loopsReader->Update();
    auto clipFunc = vtkSmartPointer<vtkImplicitSelectionLoop>::New();
    clipFunc->SetLoop(loopsReader->GetOutput()->GetPoints());
    surfCutter_->SetInsideOut(insideOut);
    surfCutter_->SetClipFunction(clipFunc);
    surfCutter = vtkAlgorithm::SafeDownCast(surfCutter_);
  }

  auto meshMapper = vtkSmartPointer<vtkDataSetMapper>::New();
  meshMapper->SetInputConnection(surfCutter->GetOutputPort(0));
  meshMapper->SetScalarModeToDefault();
  meshMapper->ScalarVisibilityOn();

  auto projLoopsMapper = vtkSmartPointer<vtkDataSetMapper>::New();
  projLoopsMapper->SetInputConnection(surfCutter->GetOutputPort(1));
  projLoopsMapper->SetScalarModeToDefault();
  projLoopsMapper->ScalarVisibilityOn();

  auto polysMapper = vtkSmartPointer<vtkDataSetMapper>::New();
  polysMapper->SetInputConnection(loopsReader->GetOutputPort());

  auto meshActor = vtkSmartPointer<vtkActor>::New();
  meshActor->GetProperty()->SetRepresentationToSurface();
  meshActor->GetProperty()->EdgeVisibilityOn();
  meshActor->SetMapper(meshMapper);

  auto projLoopsActor = vtkSmartPointer<vtkActor>::New();
  projLoopsActor->GetProperty()->SetRepresentationToWireframe();
  projLoopsActor->SetMapper(projLoopsMapper);
  projLoopsActor->GetProperty()->SetLineWidth(4);

  auto polysActor = vtkSmartPointer<vtkActor>::New();
  polysActor->SetMapper(polysMapper);
  polysActor->GetProperty()->SetRepresentationToWireframe();
  polysActor->GetProperty()->SetLineWidth(1);

  auto renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(meshActor);
  renderer->AddActor(projLoopsActor);
  renderer->AddActor(polysActor);
  //renderer->SetBackground(colors->GetColor3d("BkgColor").GetData());
  renderer->SetBackground(1.0, 1.0, 1.0);

  auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetSize(640, 480);
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("Example: SurfaceCutter");

  auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  auto start_time = std::chrono::high_resolution_clock::now();
  renderer->GetActiveCamera()->SetParallelProjection(true);
  renderer->ResetCameraClippingRange();
  renderer->ResetCamera();
  renderWindow->Render();
  auto stop_time = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);

  vtkSmartPointer<vtkPolyData> mesh, loops;
  mesh = vtkPolyData::SafeDownCast(meshTransformFilter->GetInputDataObject(0, 0));
  std::cout << "Mesh:- " << "\n";
  std::cout << " Cells : " << mesh->GetNumberOfPolys() << "\n";
  std::cout << " Points: " << mesh->GetNumberOfPoints() << "\n";
  std::cout << "Loops:- " << "\n";
  std::cout << " Cells : " << loopsReader->GetOutput()->GetNumberOfPolys() << "\n";
  std::cout << " Points: " << loopsReader->GetOutput()->GetNumberOfPoints() << "\n";
  std::cout << "Elapsed : " << elapsed_seconds.count() << "ms\n";

  auto istyle = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(istyle);
  renderWindowInteractor->Initialize();
  renderWindowInteractor->CreateRepeatingTimer(1);

  if (movable)
  {
    auto keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    keypressCallback->SetClientData(meshTransform);
    keypressCallback->SetCallback(KeypressCallbackFunction);
    renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);
  }
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}