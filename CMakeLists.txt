cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

## top content
project(TSC 
  VERSION 1.0
  DESCRIPTION "Cut triangulated surfaces with polygons"
  HOMEPAGE_URL https://github.com/jaswantp/SurfaceCutter
  LANGUAGES C CXX # C is needed if MPI is enabled in ParaView.
)
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release'.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif ()

include(CTest)
include(GNUInstallDirs)

option(TSC_BUILD_PARAVIEW_PLUGIN "Build a plugin to use SCut in ParaView." ON)
option(TSC_BUILD_SHARED_LIBS     "Generate shared libraries instead of static ones" ON)
option(TSC_BUILD_TEST            "Build test" ON)
option(TSC_BUILD_BENCH           "Build benchmark app" ON)
option(TSC_WRAP_PYTHON           "Generate and build python wrapper" ON)
option(TSC_VERSIONED_INSTALL     "Install artifacts with versioned names." ON)
set(TSC_CUSTOM_LIBRARY_SUFFIX    "<DEFAULT>"
  CACHE STRING "Custom library file name suffix (defaults to the version number)"
)

mark_as_advanced(TSC_VERSIONED_INSTALL)
mark_as_advanced(TSC_CUSTOM_LIBRARY_SUFFIX)

set(BUILD_SHARED_LIBS ${TSC_BUILD_SHARED_LIBS})
set(tsc_version_suffix "")

# Set up our directory structure for output libraries and binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

if (TSC_BUILD_SHARED_LIBS)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

if(TSC_BUILD_TEST)
  set(BUILD_TESTING ON)
endif()

if (TSC_CUSTOM_LIBRARY_SUFFIX STREQUAL "<DEFAULT>")
  if (TSC_VERSIONED_INSTALL)
    set(TSC_CUSTOM_LIBRARY_SUFFIX "${TSC_VERSION_MAJOR}.${TSC_VERSION_MINOR}")
  else ()
    set(TSC_CUSTOM_LIBRARY_SUFFIX "")
  endif ()
endif ()

if (TSC_VERSIONED_INSTALL)
  set(tsc_version_suffix "-${TSC_VERSION_MAJOR}.${TSC_VERSION_MINOR}")
endif ()

## Find ParaView, otherwise VTK
find_package(ParaView 5.10.0 QUIET)
if(NOT ParaView_FOUND)
  find_package(VTK 9.1.0 REQUIRED)
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
message (STATUS "ParaView_VERSION: ${ParaView_VERSION}")

# Install rpath for C++ and python wrapper.
if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH TRUE)
endif()
if(NOT DEFINED CMAKE_BUILD_WITH_INSTALL_RPATH)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()
if(NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

set(_vtk_module_log ALL)
vtk_module_find_modules(tsc_module_files "${PROJECT_SOURCE_DIR}/tscTriSurfaceCutter")
vtk_module_scan(
  MODULE_FILES          ${tsc_module_files}
  PROVIDES_MODULES      tsc_modules
  REQUIRES_MODULES      tsc_required_modules
  UNRECOGNIZED_MODULES  tsc_unrecognized_modules
  WANT_BY_DEFAULT       ON
)

set(tsc_install_export TSC)
if (NOT DEFINED tsc_cmake_destination)
  set(tsc_cmake_destination
    "${CMAKE_INSTALL_LIBDIR}/cmake/tsc${tsc_version_suffix}")
endif ()

vtk_module_build(
  MODULES               ${tsc_modules}
  INSTALL_EXPORT        ${tsc_install_export}
  INSTALL_HEADERS       ON
  HEADERS_DESTINATION   "${CMAKE_INSTALL_INCLUDEDIR}/tsc${tsc_version_suffix}"
  CMAKE_DESTINATION     "${tsc_cmake_destination}"
  LIBRARY_NAME_SUFFIX   "${TSC_CUSTOM_LIBRARY_SUFFIX}"
  TEST_DIRECTORY_NAME   NONE                                                 # Do not use vtk's testing.
  VERSION               "${TSC_VERSION}"
  SOVERSION             "${TSC_VERSION}"
)

if (TSC_WRAP_PYTHON)
  if(NOT DEFINED TSC_PYTHON_MODULE_DIR)
    vtk_module_python_default_destination(PYTHON_SITE_PACKAGES_SUFFIX)
    set(TSC_PYTHON_MODULE_DIR
      ${PYTHON_SITE_PACKAGES_SUFFIX}
      CACHE
      PATH
      "Where the TSC python module in installed"
    )
  endif()
  vtk_module_wrap_python(
    MODULES             ${tsc_modules}
    TARGET              TSC::Python
    WRAPPED_MODULES     tsc_wrapped_modules
    INSTALL_HEADERS     ON
    MODULE_DESTINATION  "${TSC_PYTHON_MODULE_DIR}"
    CMAKE_DESTINATION   "${tsc_cmake_destination}"
    LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"   # slaps $ORIGIN or @loader_path to CMAKE_INSTALL_RPATH
    PYTHON_PACKAGE      "tsc"
    SOABI               ${Python3_SOABI}
    INSTALL_EXPORT      TSCPython
  )

  # Generate __init__.py
  set(InitContent "from __future__ import absolute_import\n")
  foreach(module ${tsc_modules})
    set (target "")
    string(REPLACE "TSC::" "" target ${module})
    string(APPEND InitContent "from .${target} import *\n")
  endforeach()
  file(GENERATE
    OUTPUT  "${CMAKE_BINARY_DIR}/${TSC_PYTHON_MODULE_DIR}/tsc/__init__.py"
    CONTENT "${InitContent}"
  )

  # Install TSC Python
  export(
    EXPORT      TSCPython
    NAMESPACE   TSC::
    FILE        "${CMAKE_BINARY_DIR}/${tsc_cmake_destination}/TSCPython-targets.cmake"
  )
  install(
    EXPORT       TSCPython
    NAMESPACE    TSC::
    FILE         TSCPython-targets.cmake
    DESTINATION  "${tsc_cmake_destination}"
  )
  install(
    FILES       "${CMAKE_BINARY_DIR}/${TSC_PYTHON_MODULE_DIR}/tsc/__init__.py"
    DESTINATION "${TSC_PYTHON_MODULE_DIR}/tsc/"
  )
endif()

if(TSC_BUILD_PARAVIEW_PLUGIN)
  find_package(ParaView 5.10.0 REQUIRED)
  if (NOT ParaView_FOUND)
    message(FATAL_ERROR "Cannot build plugin. ParaView installation does not exist.")
    return()
  endif()
  # Install location
  if(NOT "$ENV{PV_PLUGIN_PATH}" STREQUAL "")
    set(TSC_INSTALL_PLUGIN_DIR
        $ENV{PV_PLUGIN_PATH}
        CACHE
        PATH
        "Directory where the ParaView plugin will be installed"
        )
  else()
    set(TSC_INSTALL_PLUGIN_DIR
        "${CMAKE_INSTALL_BINDIR}/plugins"
        CACHE
        PATH
        "Directory where the ParaView plugin will be installed"
        )
  endif()
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${TSC_INSTALL_PLUGIN_DIR}/TriSurfaceCutter")
  paraview_plugin_scan(
    PLUGIN_FILES      "${PROJECT_SOURCE_DIR}/paraview/tsc.plugin"
    PROVIDES_PLUGINS  tsc_plugins
    ENABLE_BY_DEFAULT ON
  )

  paraview_plugin_build(
    TARGET    TSCPVTargets
    PLUGINS   ${tsc_plugins}
    AUTOLOAD  ${tsc_plugins}
    INSTALL_HEADERS OFF
    RUNTIME_DESTINATION "${CMAKE_INSTALL_BINDIR}/plugins"
    LIBRARY_DESTINATION "${CMAKE_INSTALL_BINDIR}/plugins"
  )
  
  install(
    TARGETS     TSCPVTargets
    DESTINATION "${TSC_INSTALL_PLUGIN_DIR}"
  )
endif()

if(TSC_BUILD_BENCH)
  add_executable(tscBench tscBench.cpp)
  target_link_libraries(tscBench PUBLIC ${VTK_LIBRARIES} tscTriSurfaceCutter)
  VTK_MODULE_AUTOINIT(TARGETS tscBench MODULES ${VTK_LIBRARIES} tscTriSurfaceCutter)
endif()

# A stripped down simulation of VTK Testing, offline.
if(TSC_BUILD_TEST)
  enable_testing()
  set(TEST_SRC_DIR  ${PROJECT_SOURCE_DIR}/tscTriSurfaceCutter/Testing)
  set(TEST_DATA_FILES_FETCH_DIR ${PROJECT_SOURCE_DIR}/Testing/Data)
  set(TEST_DATA_FILES_INPUT_DIR ${PROJECT_BINARY_DIR}/Testing/Data)
  set(TEST_RUNTIME_DIR ${PROJECT_BINARY_DIR}/Testing)

  add_executable(TesttscTriSurfaceCutter ${TEST_SRC_DIR}/TesttscTriSurfaceCutter.cxx)
  target_link_libraries(TesttscTriSurfaceCutter PUBLIC ${VTK_LIBRARIES} tscTriSurfaceCutter)
  add_test(NAME TesttscTriSurfaceCutter COMMAND TesttscTriSurfaceCutter WORKING_DIRECTORY ${TEST_RUNTIME_DIR})
  
  foreach (i RANGE 5) # hard-coded
    file(COPY ${TEST_DATA_FILES_FETCH_DIR}/InOutFalse/InOutFalse_${i}.vtp DESTINATION ${TEST_DATA_FILES_INPUT_DIR}/InOutFalse)
    file(COPY ${TEST_DATA_FILES_FETCH_DIR}/InOutTrue/InOutTrue_${i}.vtp DESTINATION ${TEST_DATA_FILES_INPUT_DIR}/InOutTrue)
  endforeach()
  
  file(COPY ${TEST_DATA_FILES_FETCH_DIR}/InOutFalse.vtm  DESTINATION ${TEST_DATA_FILES_INPUT_DIR})
  file(COPY ${TEST_DATA_FILES_FETCH_DIR}/InOutTrue.vtm  DESTINATION ${TEST_DATA_FILES_INPUT_DIR})
  file(COPY ${TEST_DATA_FILES_FETCH_DIR}/Triangle.vtp  DESTINATION ${TEST_DATA_FILES_INPUT_DIR})
endif()