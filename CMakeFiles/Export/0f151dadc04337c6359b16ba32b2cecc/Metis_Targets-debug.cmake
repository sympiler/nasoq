#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "metis::metis" for configuration "Debug"
set_property(TARGET metis::metis APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(metis::metis PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/metis.lib"
  )

list(APPEND _cmake_import_check_targets metis::metis )
list(APPEND _cmake_import_check_files_for_metis::metis "${_IMPORT_PREFIX}/lib/metis.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
