# SurfaceCutter
Cut triangulated surfaces with 2D loop polygons.

### Usage
```c++
#include <SurfaceCutter.h>

vtkNew<SurfaceCutter> cutter;

// xxx: filter that outputs a vtkPolyData
cutter->SetInputConnection(0, xxx->GetOutputPort());
// yyy: filter that outputs a vtkPolyData
cutter->SetInputConnection(1, yyy->GetOutputPort());
// remove stuff outside loops.
cutter->InsideOutOn();
// Loop points will be tagged '1'. 
cutter->ColorAcquiredPtsOn();
// Loop edges projected upon surface will be tagged '1'.
cutter->ColorLoopEdgesOn();

cutter->Update();
cutMesh = cutter->GetOutput(0); 

```

## InsideOut = true
![InsideOut = true](illustrations/testInsideOutTrue.png)
![true](illustrations/inOutTrue.gif)

## InsideOut = false
![InsideOut = true](illustrations/testInsideOutFalse.png)
![false](illustrations/inOutFalse.gif)

## Translation/Rotation
With [SurfaceCutter](https://github.com/jaswantp/SurfaceCutter), it should be possible to
compute a cut mesh in realtime.
![Translationx](illustrations/translationx.gif)
![Translationy](illustrations/translationy.gif)
![rotation](illustrations/rotation.gif)


## Usage (benchmark)
```
./benchmark.exe <option(s)>Options:
        -h,--help       Show this help message
        -m,--mesh       Specify mesh file (*.vtp, *.vtu)
        -l,--loops      Specify loops file (*.vtp)
        -i,--invert     Invert 2d boolean. Portions inside loops will be removed.
        -t,--translationspeed   Speed multiplier for mesh translations along x, y, z
        -r,--rotationspeed      Speed multiplier for mesh rotation along z
           --movable    Make the mesh movable.
           --vtkcookiecutter    Use vtkCookieCutter instead
           --vtkclipdataset     Use vtkClipDataset with vtkImplicitSelectionLoop instead


Controls:
W:    Z+ | S:     Z-
Up:   Y+ | Down:  Y-
Left: X+ | Right: X-
Z:   CCW | C:     CW (Looking down Z-)
```
