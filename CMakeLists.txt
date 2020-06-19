cmake_minimum_required(VERSION 3.9)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

project(NASOQ-BIN)
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_C_COMPILER /usr/local/Cellar/gcc/9.1.0/bin/gcc-9)
SET(CMAKE_CXX_COMPILER /usr/local/Cellar/gcc/9.1.0/bin/g++-9)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
set(hack /opt/intel/compilers_and_libraries_2018.3.185/mac/compiler/lib/)
set(hack2 /Users/kazem/programs/metis-5.1.0/build/Darwin-x86_64/libmetis/)
link_directories(${hack} ${hack2})
# TODO: remove MKL dependency using OpenBLAS submodule
find_package(MKL REQUIRED)

# Adding hints for suitesparse
set(SUITESPARSE_INCLUDE_DIR_HINTS ${SUITESPARSE_INCLUDE_DIR_HINTS} ${SUITE_ROOT_PATH}/include)
set(SUITESPARSE_LIBRARY_DIR_HINTS ${SUITESPARSE_LIBRARY_DIR_HINTS} ${SUITE_ROOT_PATH}/lib)
set(BLA_STATIC TRUE)
find_package(SuiteSparse OPTIONAL_COMPONENTS)

#TODO adding METIS submodule

#IF(NOT SUITESPARSE_FOUND)
    set( INCLUDE_INSTALL_DIR ${METIS_ROOT_PATH}/include)
    set( LIB_INSTALL_DIR ${METIS_ROOT_PATH}/libmetis)
    find_package(METIS REQUIRED)
    set(METIS_LIBRARY "")
    set(METIS_INCLUDE_DIR "")
#ELSE()
#    set(ignoreMe "${METIS_ROOT_PATH}")
#ENDIF()
message("====${METIS_LIBRARIES}")
find_package(OpenMP)

find_package (Eigen3 OPTIONAL_COMPONENTS NO_MODULE)

if(UNIX )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 ")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DMKL  -Wno-dev ")

set(NASOQ_INCLUDE_DIR ${MKL_INCLUDE_DIR} ${SUITESPARSE_INCLUDE_DIRS}
        ${METIS_INCLUDE_DIR} ${METIS_INCLUDES}
        "${CMAKE_CURRENT_LIST_DIR}/symbolic/" "${CMAKE_CURRENT_LIST_DIR}/common/"
        "${CMAKE_CURRENT_LIST_DIR}/ldl/" "${CMAKE_CURRENT_LIST_DIR}/matrixMatrix/"
        "${CMAKE_CURRENT_LIST_DIR}/matrixVector/"  "${CMAKE_CURRENT_LIST_DIR}/linear_solver/"
        "${CMAKE_CURRENT_LIST_DIR}/gmres/" "${CMAKE_CURRENT_LIST_DIR}/QP/"
        "${CMAKE_CURRENT_LIST_DIR}/triangularSolve/"
        "${CMAKE_CURRENT_LIST_DIR}/smp-format/")

set(NASOQ_LIB_DIR ${MKL_LIBRARIES} ${SUITESPARSE_LIBRARIES}
        ${METIS_LIBRARY} ${METIS_LIBRARIES})

set(SOURCE_FILES nasoq_driver.cpp)
include_directories(${MKL_INCLUDE_DIR} ${SUITESPARSE_INCLUDE_DIRS}
        ${METIS_INCLUDE_DIR} ${METIS_INCLUDES}
        ${EIGEN3_INCLUDE_DIR}
        "symbolic/" "common/" "ldl/" "matrixMatrix/"
        "matrixVector/"  "linear_solver/" "gmres/" "QP/"
        "triangularSolve/")
add_executable(NASOQ ${SOURCE_FILES}  nasoq_driver.cpp)
target_link_libraries(NASOQ PUBLIC
        ${MKL_LIBRARIES} ${SUITESPARSE_LIBRARIES}
        ${METIS_LIBRARY} ${METIS_LIBRARIES}
        )
if(OpenMP_CXX_FOUND)
    target_link_libraries(NASOQ-BIN PRIVATE OpenMP::OpenMP_CXX)
endif()

### Building the interface library
project(NASOQ)
#file(GLOB_RECURSE SRC_FILES ${NASOQ_INCLUDE_DIR} "*.h" )
add_library (NASOQ  INTERFACE )
target_include_directories(NASOQ INTERFACE  ${NASOQ_INCLUDE_DIR} )
#target_sources(NASOQ ${SRC_FILES})
target_link_libraries(NASOQ INTERFACE
        ${NASOQ_LIB_DIR}
        )
## TODO static and shared

if(TARGET Eigen3::Eigen OR EXISTS "${EIGEN3_INCLUDE_DIR}")
    add_subdirectory(eigen_interface)
endif()
