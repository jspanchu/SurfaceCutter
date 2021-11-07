# SurfaceCutter
Cut triangulated surfaces with 2D loop polygons.

### Build
```shell
mkdir build
cd build
cmake .. -GNinja -DTSC_BUILD_PARAVIEW_PLUGIN=ON -DTSC_WRAP_PYTHON=ON
ninja install
```

### Usage
#### ParaView
Filters -> SCUT -> Triangulation Cutter
#### C++
```c++
#include <tscTriSurfaceCutter.h>

vtkNew<tscTriSurfaceCutter> cutter;

// xxx: filter that outputs a vtkPolyData
cutter->SetInputConnection(0, xxx->GetOutputPort());
// yyy: filter that outputs a vtkPolyData
cutter->SetInputConnection(1, yyy->GetOutputPort());

cutter->Update();
cutMesh = cutter->GetOutput(0); 
```
#### Python
```Python
import tsc
...
cutter = tsc.tscTriSurfaceCutter()
cutter.SetInputConnection(triangulation.GetOutputPort())
cutter.SetLoopsConnection(1, loops.GetOutputPort())
...
```

## InsideOut = true
<a href="https://imgur.com/q1xuGLU"><img src="https://i.imgur.com/q1xuGLU.png" title="source: imgur.com" width="400"/></a>
<a href="https://imgur.com/fs6q4st"><img src="https://i.imgur.com/fs6q4st.gif" title="source: imgur.com" width="400"/></a>

## InsideOut = false
<a href="https://imgur.com/kP0lQTa"><img src="https://i.imgur.com/kP0lQTa.png" title="source: imgur.com" width="400"/></a>
<a href="https://imgur.com/uXsGpo2"><img src="https://i.imgur.com/uXsGpo2.gif" title="source: imgur.com" width="400"/></a>

## **Imprint**: *On* | **Remove**: *Off* | **Invert**: NA
<a href="https://imgur.com/6sFCEFX"><img src="https://i.imgur.com/6sFCEFX.png" title="source: imgur.com" width="400"/></a>

## **Imprint**: *On* | **Remove**: *On* | **Invert**: *On*
<a href="https://imgur.com/L5MPDBc"><img src="https://i.imgur.com/L5MPDBc.jpg" title="source: imgur.com" width="400"/></a>

## **Imprint**: *On* | **Remove**: *On* | **Invert**: *Off*
<a href="https://imgur.com/6NJYc3W"><img src="https://i.imgur.com/6NJYc3W.jpg" title="source: imgur.com" width="400"/></a>

## **Imprint**: *Off* | Remove: *On* | **Invert**: *On*
<a href="https://imgur.com/llD8D7y"><img src="https://i.imgur.com/llD8D7y.jpg" title="source: imgur.com" width="400"/></a>

## **Imprint**: *Off* | Remove: *On* | **Invert**: *Off*
<a href="https://imgur.com/TjEd6Rz"><img src="https://i.imgur.com/TjEd6Rz.jpg" title="source: imgur.com" width="400"/></a>

## Real-time
<a href="https://imgur.com/3ewDGsx"><img src="https://i.imgur.com/3ewDGsx.gif" title="source: imgur.com" width="400"/></a>
<a href="https://imgur.com/ctTqun2"><img src="https://i.imgur.com/ctTqun2.gif" title="source: imgur.com" width="400"/></a>
<a href="https://imgur.com/8DeC1yF"><img src="https://i.imgur.com/8DeC1yF.gif" title="source: imgur.com" width="400"/></a>
<a href="https://imgur.com/A0dHcnR"><img src="https://i.imgur.com/A0dHcnR.gif" title="source: imgur.com" width="400"/></a>



## Benchmark
```
./benchmark <option(s)>Options:
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
K:ScaleUp| H:ScaleDown
```
