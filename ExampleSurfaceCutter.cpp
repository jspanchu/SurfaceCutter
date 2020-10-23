#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkCylinderSource.h>
#include <vtkDataSetMapper.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkXMLPolyDataReader.h>

#include <SurfaceCutter.h>

#include <array>
#include <iostream>

vtkSmartPointer<vtkTransform> meshTransform = vtkSmartPointer<vtkTransform>::New();
vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

static double rotation = 0;

class RotateCommand : public vtkCommand
{
public:
  vtkTypeMacro(RotateCommand, vtkCommand);

  static RotateCommand* New()
  {
    return new RotateCommand;
  }

  void Execute(vtkObject* caller,
    unsigned long vtkNotUsed(eventId),
    void* callData)
  {
    //rotation += 0.1;
    //std::cout << rotation << std::endl;
    //meshTransform->RotateZ(0.1);
    renderWindow->Render();
  }
};

int main(int, char* []) {
  vtkSmartPointer<vtkNamedColors> colors =
    vtkSmartPointer<vtkNamedColors>::New();

  // Set the background color.
  std::array<unsigned char, 4> bkg{ {26, 51, 102, 255} };
  colors->SetColor("BkgColor", bkg.data());

  //auto mesh = vtkSmartPointer<vtkPolyData>::New();
  auto meshReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  meshReader->SetFileName("Testing/triangle.vtp");

  auto meshTransformFilter = vtkSmartPointer<vtkTransformFilter>::New();
  meshTransformFilter->SetTransform(meshTransform);
  meshTransform->RotateZ(233.5);
  meshTransformFilter->SetInputConnection(meshReader->GetOutputPort());

  auto polysReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  polysReader->SetFileName("Testing/Case5.vtp");
  
  auto surfCutter = vtkSmartPointer<SurfaceCutter>::New();
  surfCutter->SetInputConnection(0, meshTransformFilter->GetOutputPort());
  surfCutter->SetInputConnection(1, polysReader->GetOutputPort());

  auto meshMapper = vtkSmartPointer<vtkDataSetMapper>::New();
  meshMapper->SetInputConnection(surfCutter->GetOutputPort());
  meshMapper->SetScalarModeToDefault();
  meshMapper->ScalarVisibilityOn();

  auto polysMapper = vtkSmartPointer<vtkDataSetMapper>::New();
  polysMapper->SetInputConnection(polysReader->GetOutputPort());

  auto meshActor = vtkSmartPointer<vtkActor>::New();
  meshActor->GetProperty()->SetRepresentationToSurface();
  meshActor->GetProperty()->EdgeVisibilityOn();
  meshActor->SetMapper(meshMapper);

  auto polysActor = vtkSmartPointer<vtkActor>::New();
  polysActor->SetMapper(polysMapper);
  polysActor->GetProperty()->SetRepresentationToWireframe();
  polysActor->GetProperty()->SetLineWidth(6);
  polysActor->GetProperty()->SetOpacity(0.5);

  auto renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(meshActor);
  renderer->AddActor(polysActor);
  renderer->SetBackground(colors->GetColor3d("BkgColor").GetData());

  renderWindow->SetSize(640, 480);
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("Example: SurfaceCutter");

  auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderWindow->Render();
  auto istyle = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(istyle);
  renderWindowInteractor->Initialize();
  renderWindowInteractor->CreateRepeatingTimer(1);

  RotateCommand* rotateCallback = RotateCommand::New();
  renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, rotateCallback);
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}