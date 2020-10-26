cmake_minimum_required(VERSION 3.12)

project (SurfaceCutter)

option (BUILD_EXAMPLE "Build example" ON)
option (BUILD_TEST "Build test" ON)

if (BUILD_EXAMPLE)
  find_package(VTK COMPONENTS 
      vtkCommonCore
      vtkCommonDataModel
      vtkCommonExecutionModel
      vtkCommonTransforms
      vtkFiltersCore
      vtkFiltersExtraction
      vtkFiltersFlowPaths 
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
 "SurfaceCutterHelper.h")
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
  target_compile_options(${PROJECT_NAME} PRIVATE "$<$<CONFIG:Release>:/Z7>")
  target_link_options(${PROJECT_NAME} PRIVATE "$<$<CONFIG:Release>:/DEBUG>")
  target_link_options(${PROJECT_NAME} PRIVATE "$<$<CONFIG:Release>:/OPT:REF>")
  target_link_options(${PROJECT_NAME} PRIVATE "$<$<CONFIG:Release>:/OPT:ICF>")
endif () 

file (COPY "${CMAKE_CURRENT_SOURCE_DIR}/data" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

if (BUILD_EXAMPLE)
  add_executable (benchmark 
    "${PROJECT_NAME}Benchmark.cpp"
   )
   
  if (VTK_VERSION VERSION_LESS "8.90.0")
    # old system
    include(${VTK_USE_FILE})
    target_link_libraries(benchmark PRIVATE ${VTK_LIBRARIES} ${PROJECT_NAME})
  else ()
    # include all components
    target_link_libraries(benchmark PRIVATE ${VTK_LIBRARIES} ${PROJECT_NAME})
    # vtk_module_autoinit is needed
    vtk_module_autoinit(
      TARGETS benchmark
      MODULES ${VTK_LIBRARIES}
      )
  endif ()
  target_compile_options(benchmark PRIVATE "$<$<CONFIG:Release>:/Zi>")
  target_link_options(benchmark PRIVATE "$<$<CONFIG:Release>:/DEBUG>")
  target_link_options(benchmark PRIVATE "$<$<CONFIG:Release>:/OPT:REF>")
  target_link_options(benchmark PRIVATE "$<$<CONFIG:Release>:/OPT:ICF>")
endif ()

if (BUILD_TEST)
  enable_testing()
  add_executable (${PROJECT_NAME}Test
    "${PROJECT_NAME}Test.cpp"
   )
  add_test(NAME TestSurfaceCutter COMMAND TestSurfaceCutter )
   
  if (VTK_VERSION VERSION_LESS "8.90.0")
    # old system
    include(${VTK_USE_FILE})
    target_link_libraries(${PROJECT_NAME}Test PRIVATE ${VTK_LIBRARIES} ${PROJECT_NAME})
  else ()
    # include all components
    target_link_libraries(${PROJECT_NAME}Test PRIVATE ${VTK_LIBRARIES} ${PROJECT_NAME})
    # vtk_module_autoinit is needed
    vtk_module_autoinit(
      TARGETS ${PROJECT_NAME}Test
      MODULES ${VTK_LIBRARIES}
      )
  endif () 
endif ()