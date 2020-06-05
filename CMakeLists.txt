cmake_minimum_required(VERSION 3.9)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

project(NASOQ)
set(CMAKE_CXX_STANDARD 11)

# TODO: remove MKL dependency using OpenBLAS submodule
find_package(MKL REQUIRED)

# Adding hints for suitesparse
set(SUITESPARSE_INCLUDE_DIR_HINTS ${SUITESPARSE_INCLUDE_DIR_HINTS} ${SUITE_ROOT_PATH}/include)
set(SUITESPARSE_LIBRARY_DIR_HINTS ${SUITESPARSE_LIBRARY_DIR_HINTS} ${SUITE_ROOT_PATH}/lib)
set(BLA_STATIC TRUE)
find_package(SuiteSparse OPTIONAL_COMPONENTS)

#TODO adding METIS submodule
IF(NOT SUITESPARSE_FOUND)
set( INCLUDE_INSTALL_DIR ${METIS_ROOT_PATH}/include)
set( LIB_INSTALL_DIR ${METIS_ROOT_PATH}/libmetis)
find_package(METIS REQUIRED)
set(METIS_LIBRARY "")
set(METIS_INCLUDE_DIR "")
ELSE()
    set(ignoreMe "${METIS_ROOT_PATH}")
ENDIF()

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
        ${METIS_LIBRARY} ${METIS_LIBRARIES} )

set(SOURCE_FILES nasoq_driver.cpp)
include_directories(${NASOQ_INCLUDE_DIR} )
add_executable(NASOQ ${SOURCE_FILES}  nasoq_driver.cpp)
target_link_libraries(NASOQ PUBLIC
        ${NASOQ_LIB_DIR})
if(OpenMP_CXX_FOUND)
    target_link_libraries(NASOQ PUBLIC OpenMP::OpenMP_CXX)
endif()

IF(TARGET Eigen3::Eigen)
add_subdirectory(eigen_interface)
ENDIF(TARGET Eigen3::Eigen)