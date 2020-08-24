cmake_minimum_required(VERSION 3.14)

#add_subdirectory(nasoq)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)


### NASOQ example
project(NASOQ_Test)
set(CMAKE_CXX_STANDARD 11)

# TODO: remove MKL dependency using OpenBLAS submodule
find_package(MKL REQUIRED)

set(SUITESPARSE_INCLUDE_DIR_HINTS ${SUITESPARSE_INCLUDE_DIR_HINTS} ${SUITE_ROOT_PATH}/include)
set(SUITESPARSE_LIBRARY_DIR_HINTS ${SUITESPARSE_LIBRARY_DIR_HINTS} ${SUITE_ROOT_PATH}/lib)
set(BLA_STATIC TRUE)
find_package(SuiteSparse OPTIONAL_COMPONENTS)

set( INCLUDE_INSTALL_DIR ${METIS_ROOT_PATH}/include)
set( LIB_INSTALL_DIR ${METIS_ROOT_PATH}/libmetis)
find_package(METIS REQUIRED)
set(METIS_LIBRARY "")
set(METIS_INCLUDE_DIR "")

include_directories(${METIS_INCLUDES} ${MKL_INCLUDE_DIR}
        "nasoq/QP/" "nasoq/symbolic/" "nasoq/common/" "nasoq/ldl/"
        "nasoq/triangularSolve/" "nasoq/linear_solver/" "nasoq/gmres/"
        "nasoq/matrixVector/" "nasoq/matrixMatrix/")

set(SOURCE_FILES_NAS nasoq_main.cpp)
add_executable(NASOQ_Test ${SOURCE_FILES_NAS} nasoq_main.cpp)
target_link_libraries(NASOQ_Test PRIVATE
        ${MKL_LIBRARIES} ${SUITESPARSE_LIBRARIES}
        ${METIS_LIBRARY} ${METIS_LIBRARIES}        )
if(OpenMP_CXX_FOUND)
    target_link_libraries(NASOQ_Test PRIVATE OpenMP::OpenMP_CXX)
endif()





### LBL example
project(LBL_Test)
set(CMAKE_CXX_STANDARD 11)

# TODO: remove MKL dependency using OpenBLAS submodule
find_package(MKL REQUIRED)

set(SUITESPARSE_INCLUDE_DIR_HINTS ${SUITESPARSE_INCLUDE_DIR_HINTS} ${SUITE_ROOT_PATH}/include)
set(SUITESPARSE_LIBRARY_DIR_HINTS ${SUITESPARSE_LIBRARY_DIR_HINTS} ${SUITE_ROOT_PATH}/lib)
set(BLA_STATIC TRUE)
find_package(SuiteSparse OPTIONAL_COMPONENTS)

set( INCLUDE_INSTALL_DIR ${METIS_ROOT_PATH}/include)
set( LIB_INSTALL_DIR ${METIS_ROOT_PATH}/libmetis)
find_package(METIS REQUIRED)
set(METIS_LIBRARY "")
set(METIS_INCLUDE_DIR "")

include_directories(${METIS_INCLUDES} ${MKL_INCLUDE_DIR}
        "nasoq/QP/" "nasoq/symbolic/" "nasoq/common/" "nasoq/ldl/"
        "nasoq/triangularSolve/" "nasoq/linear_solver/" "nasoq/gmres/"
        "nasoq/matrixVector/" "nasoq/matrixMatrix/")

set(SOURCE_FILES_NAS lbl_main.cpp)
add_executable(LBL_Test ${SOURCE_FILES_NAS} lbl_main.cpp)
target_link_libraries(LBL_Test PRIVATE
        ${MKL_LIBRARIES} ${SUITESPARSE_LIBRARIES}
        ${METIS_LIBRARY} ${METIS_LIBRARIES}        )
if(OpenMP_CXX_FOUND)
    target_link_libraries(LBL_Test PRIVATE OpenMP::OpenMP_CXX)
endif()
