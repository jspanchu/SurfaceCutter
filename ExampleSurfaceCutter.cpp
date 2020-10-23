#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCylinderSource.h>
#include <vtkNamedColors.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>

#include <SurfaceCutter.h>

#include <array>

int main(int, char* []) {
  vtkSmartPointer<vtkNamedColors> colors =
    vtkSmartPointer<vtkNamedColors>::New();

  // Set the background color.
  std::array<unsigned char, 4> bkg{ {26, 51, 102, 255} };
  colors->SetColor("BkgColor", bkg.data());

  auto mesh = vtkSmartPointer<vtkPolyData>::New();
  auto meshReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  meshReader->SetFileName("Testing/triangle.vtp");

  auto polysReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  polysReader->SetFileName("Testing/Case5.vtp");
  
  auto surfCutter = vtkSmartPointer<SurfaceCutter>::New();
  surfCutter->SetInputConnection(0, meshReader->GetOutputPort());
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

  auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetSize(1920, 1080);
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("Example: SurfaceCutter");

  auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderWindow->Render();
  auto istyle = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(istyle);
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}