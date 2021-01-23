cmake_minimum_required(VERSION 3.12)

project (SurfaceCutter)

option (BUILD_EXAMPLE "Build example" ON)
option (BUILD_TEST "Build test" ON)

find_package (VTK)

if (NOT VTK_FOUND)
  find_package (ParaView)
  if (NOT VTK_FOUND)
    message("Skipping ${PROJECT_NAME}: ${VTK_NOT_FOUND_MESSAGE}")
    return ()
  endif ()
endif ()

message (STATUS "VTK_VERSION: ${VTK_VERSION}")

add_library (${PROJECT_NAME} STATIC
  "${PROJECT_NAME}.cxx"
  "${PROJECT_NAME}.h"
)
target_include_directories (${PROJECT_NAME} 
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
endif () 

if (BUILD_EXAMPLE OR BUILD_TEST)
  file (COPY "${CMAKE_CURRENT_SOURCE_DIR}/data" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
endif ()

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
    vtk_module_autoinit(
      TARGETS benchmark
      MODULES ${VTK_LIBRARIES}
      )
  endif ()
endif ()

if (BUILD_TEST)
  enable_testing()
  add_executable (Test${PROJECT_NAME}
    "Test${PROJECT_NAME}.cpp"
   )
  add_test(NAME TestSurfaceCutter COMMAND TestSurfaceCutter )
   
  target_link_libraries(Test${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES} ${PROJECT_NAME})
  if (VTK_VERSION VERSION_LESS "8.90.0")
    # old system
    include(${VTK_USE_FILE})
  endif () 
endif ()