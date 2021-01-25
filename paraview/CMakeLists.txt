paraview_add_plugin(TriSurfaceCutter
    VERSION             ${TSC_VERSION}
    REQUIRED_ON_CLIENT
    REQUIRED_ON_SERVER
    MODULES             ${tsc_modules}
    # MODULE_FILES        ${tsc_module_files}
    SERVER_MANAGER_XML  TriSurfaceCutter.xml
)

set_target_properties(TriSurfaceCutter
PROPERTIES
  INSTALL_RPATH "$ORIGIN/../../../lib"
)