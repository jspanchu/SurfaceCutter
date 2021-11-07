paraview_add_plugin(TriSurfaceCutter
    VERSION             ${TSC_VERSION}
    REQUIRED_ON_CLIENT
    REQUIRED_ON_SERVER
    MODULES             ${TSC_MODULES}
    SERVER_MANAGER_XML  TriSurfaceCutter.xml
)

if(NOT "${TTK_INSTALL_PLUGIN_DIR}" STREQUAL "") 
  if(APPLE)
    # On macOS,
    # look into the subdirectory "TriSurfaceCutter"
    # to find the actual plugins.
    get_target_property(TEMP
        TriSurfaceCutter INSTALL_RPATH
        )
    set_target_properties(${TARGET_NAME}
      PROPERTIES
        INSTALL_RPATH "@loader_path/${TTK_PLUGIN_SUBDIR};${TEMP}"
    )
  endif(APPLE)

  install(
    TARGETS
      TriSurfaceCutter
    DESTINATION
      "${TSC_INSTALL_PLUGIN_DIR}/${TSC_PLUGIN_SUBDIR}"
    )
endif()