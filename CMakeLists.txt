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
 "SurfaceCutterImpl.h")
target_include_directories (SurfaceCutter 
  PUBLIC
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

option (BUILD_TESTING "Build tests" ON)

if (BUILD_TESTING)
  find_package(VTK COMPONENTS 
      vtkCommonCore
      vtkCommonDataModel
      vtkCommonExecutionModel
      vtkCommonTransforms
      vtkFiltersCore
      vtkFiltersGeneral
      vtkFiltersModeling
      vtkIOXML
      vtkIOXMLParser
      vtkkwiml
      vtkInteractionStyle
      vtkRenderingAnnotation
      vtkRenderingCore
      vtkRenderingFreeType
      vtkRenderingGL2PSOpenGL2
      vtkRenderingOpenGL2
    QUIET
    )

  add_executable (TestSurfaceCutter 
    "TestSurfaceCutter.cpp"
   )
   
  if (VTK_VERSION VERSION_LESS "8.90.0")
    # old system
    include(${VTK_USE_FILE})
    target_link_libraries(TestSurfaceCutter PRIVATE ${VTK_LIBRARIES} SurfaceCutter)
  else ()
    # include all components
    target_link_libraries(TestSurfaceCutter PRIVATE ${VTK_LIBRARIES} SurfaceCutter)
    # vtk_module_autoinit is needed
    vtk_module_autoinit(
      TARGETS SurfaceCutter
      MODULES ${VTK_LIBRARIES}
      )
  endif () 

endif ()