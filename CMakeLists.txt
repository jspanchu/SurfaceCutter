cmake_minimum_required(VERSION 3.12)

project (SurfaceCutter)

option (BUILD_EXAMPLE "Build example" ON)
if (BUILD_EXAMPLE)
  find_package(VTK COMPONENTS 
      vtkCommonCore
      vtkCommonDataModel
      vtkCommonExecutionModel
      vtkCommonTransforms
      vtkFiltersCore
      vtkFiltersExtraction
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
else()
  find_package(VTK COMPONENTS 
    vtkCommonCore
    vtkCommonDataModel
    vtkCommonExecutionModel
    vtkFiltersCore
    vtkFiltersGeneral
    QUIET
    )
endif ()
 
if (NOT VTK_FOUND)
  message("Skipping ${PROJECT_NAME}: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()

message (STATUS "VTK_VERSION: ${VTK_VERSION}")

add_library (${PROJECT_NAME} 
  "${PROJECT_NAME}.cxx"
  "${PROJECT_NAME}.h"
)
target_include_directories (${PROJECT_NAME} 
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
    )
endif () 

if (BUILD_EXAMPLE)
  add_executable (Example${PROJECT_NAME} 
    "Example${PROJECT_NAME}.cpp"
   )
   
  if (VTK_VERSION VERSION_LESS "8.90.0")
    # old system
    include(${VTK_USE_FILE})
    target_link_libraries(Example${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES} ${PROJECT_NAME})
  else ()
    # include all components
    target_link_libraries(Example${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES} ${PROJECT_NAME})
    # vtk_module_autoinit is needed
    vtk_module_autoinit(
      TARGETS Example${PROJECT_NAME}
      MODULES ${VTK_LIBRARIES}
      )
  endif () 
  file (COPY "${CMAKE_CURRENT_SOURCE_DIR}/Testing" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

endif ()