cmake_minimum_required(VERSION 3.12)

project (SurfaceCutter)

find_package(VTK COMPONENTS 
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkFiltersCore
  vtkFiltersGeneral
  QUIET
  )
 
if (NOT VTK_FOUND)
  message("Skipping SurfaceCutter: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()

message (STATUS "VTK_VERSION: ${VTK_VERSION}")

add_library (SurfaceCutter 
  "SurfaceCutter.cxx"
  "SurfaceCutter.h"
)
target_include_directories (SurfaceCutter 
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  target_link_libraries(SurfaceCutter PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  target_link_libraries(SurfaceCutter PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS SurfaceCutter
    MODULES ${VTK_LIBRARIES}
    )
endif () 