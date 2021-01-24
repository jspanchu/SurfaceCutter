if (NOT MSVC) # GCC and Clang
  # warning flags
  list(APPEND TSC_COMPILER_FLAGS -Wall -Wshadow)

  # performance and debug flags
  if(CMAKE_BUILD_TYPE MATCHES Release)
    list(APPEND TSC_COMPILER_FLAGS -march=native -Wfatal-errors)
  endif()

  # hardened linker flags for Clang and GCC
  if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    list(APPEND TSC_LINKER_FLAGS
      -Wl,--as-needed
      -Wl,--no-undefined
      )
  endif()
endif()

set (mod_name tscTriSurfaceCutter)
set (classes tscTriSurfaceCutter)
vtk_module_add_module(${mod_name}
  CLASSES ${classes}
)
vtk_module_compile_options(${mod_name} PRIVATE ${TSC_COMPILER_FLAGS})

if(TSC_LINKER_FLAGS)
  vtk_module_link_options(${mod_name} PRIVATE ${TSC_LINKER_FLAGS})
endif()