cmake_minimum_required(VERSION 3.1)
set(CMAKE_BUILD_TYPE_INIT Release)
project(nasoq_mex)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Move output to this dir (gptoolbox/mex)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
#SET(CMAKE_C_COMPILER /usr/local/Cellar/gcc/9.1.0/bin/gcc-9)
#SET(CMAKE_CXX_COMPILER /usr/local/Cellar/gcc/9.1.0/bin/g++-9)

# From libigl/cmake/libigl.cmake libigl/cmake/LibiglWindows.cmake
if(MSVC)
  if("${MSVC_RUNTIME}" STREQUAL "")
    set(MSVC_RUNTIME "static")
  endif()
  if(${MSVC_RUNTIME} STREQUAL "static")
    message(STATUS "MSVC -> forcing use of statically-linked runtime.")
    foreach(config ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER ${config} config)
      string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
      string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
    endforeach()
    string(REPLACE "/MDd" "/MTd" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
  else()
    message(STATUS "MSVC -> forcing use of dynamically-linked runtime.")
    foreach(config ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER ${config} config)
      string(REPLACE /MT /MD CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
      string(REPLACE /MT /MD CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
    endforeach()
    string(REPLACE "/MTd" "/MDd" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
  endif()
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(GPTOOLBOX_MEX_ROOT "${PROJECT_SOURCE_DIR}")
set(GPTOOLBOX_MEX_EXTERNAL "${GPTOOLBOX_MEX_ROOT}/external")

# libigl
option(LIBIGL_RESTRICTED_MATLAB "Use Matlab" ON)
include(libigl)


set( CMAKE_SHARED_LIBRARY_SUFFIX ".${Matlab_MEX_EXTENSION}" )  # set suffix to .mexa64
set( CMAKE_SHARED_LIBRARY_PREFIX )      # remove the "lib" prefix
set( CMAKE_SHARED_LIBRARY_CXX_FLAGS )       # remove the -fPIC option. mex does not accept the "-fPIC" option
set( CMAKE_POSITION_INDEPENDENT_CODE ON)
set( CMAKE_INSTALL_RPATH "\$ORIGIN" )

# TODO: remove MKL dependency using OpenBLAS submodule
find_package(MKL REQUIRED)

if(UNIX )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 ")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DMKL  ")

# Adding hints for suitesparse
set(SUITESPARSE_INCLUDE_DIR_HINTS ${SUITESPARSE_INCLUDE_DIR_HINTS} ${SUITE_ROOT_PATH}/include)
set(SUITESPARSE_LIBRARY_DIR_HINTS ${SUITESPARSE_LIBRARY_DIR_HINTS} ${SUITE_ROOT_PATH}/lib)
set(BLA_STATIC TRUE)
find_package(SuiteSparse OPTIONAL_COMPONENTS)

#TODO needs a fix
set( INCLUDE_INSTALL_DIR ${METIS_ROOT_PATH}/include)
set( LIB_INSTALL_DIR ${METIS_ROOT_PATH}/libmetis)
find_package(METIS REQUIRED)
set(METIS_LIBRARY "")
set(METIS_INCLUDE_DIR "")


function(compile_each I_MEX_SOURCES I_LIBS )
  foreach(MEX_SOURCE ${I_MEX_SOURCES})
    message(STATUS "MEX_SOURCE: ${MEX_SOURCE}")
    get_filename_component(MEX_OUTPUT ${MEX_SOURCE} NAME_WE)
    set(MEX_TARGET "${MEX_OUTPUT}-target")
    add_library(${MEX_TARGET} SHARED ${MEX_SOURCE})
    set_target_properties( ${MEX_TARGET} PROPERTIES OUTPUT_NAME ${MEX_OUTPUT})
    # Not sure if this should be (WIN32) or (MSVC OR MSYS OR MINGW)
    # https://stackoverflow.com/a/40152725/148668
    if(WIN32)
      # https://stackoverflow.com/a/11095194/148668
      set_target_properties(${MEX_TARGET} PROPERTIES LINK_FLAGS "/export:mexFunction")
    endif()
    target_link_libraries(${MEX_TARGET} ${I_LIBS} ${METIS_LIBRARIES}
      ${SUITESPARSE_LIBRARIES} ${MKL_LIBRARIES})
    target_include_directories(${MEX_TARGET} PUBLIC ${Matlab_INCLUDE_DIRS}
      ${METIS_INCLUDE_DIR} ${METIS_INCLUDES} ${MKL_INCLUDE_DIR}
      ${SUITESPARSE_INCLUDE_DIRS}) 
    list(APPEND MEX_TARGETS ${MEX_TARGET})
    if(APPLE)
      include(OSXFixDylibReferences)
      OSX_FIX_DYLIB_REFERENCES(${MEX_TARGET} ${MKL_LIBRARIES})
    endif()
  endforeach()
endfunction()

#file( GLOB MEX_SOURCES *.cpp )
list(APPEND CORE_LIBS ${Matlab_LIBRARIES})

# LIBIGL
add_definitions("-DMEX")
compile_each("\
nasoq.cpp;\
"
  "nasoq::nasoq;igl::core;igl_restricted::matlab")

