#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "metis::metis" for configuration "Release"
set_property(TARGET metis::metis APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(metis::metis PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/metis.lib"
  )

list(APPEND _cmake_import_check_targets metis::metis )
list(APPEND _cmake_import_check_files_for_metis::metis "${_IMPORT_PREFIX}/lib/metis.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)