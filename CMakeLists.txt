cmake_minimum_required(VERSION 3.9)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(NASOQ)

set(CMAKE_CXX_STANDARD 11)

find_package(MKL REQUIRED)
set(BLA_STATIC TRUE)
find_package(SuiteSparse REQUIRED)

find_package(OpenMP)


if(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 ")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL  ")

set(SOURCE_FILES nasoq_driver.cpp)
include_directories(${MKL_INCLUDE_DIR} ${SUITESPARSE_INCLUDE_DIRS} ${METIS_INCLUDE_DIR}
        "symbolic/" "common/" "ldl/" "matrixMatrix/"
        "matrixVector/"  "linear_solver/" "gmres/" "QP/"
        "triangularSolve/")
add_executable(NASOQ ${SOURCE_FILES}  nasoq_driver.cpp)
target_link_libraries(NASOQ PUBLIC
        ${MKL_LIBRARIES} ${SUITESPARSE_LIBRARIES} ${METIS_LIBRARY}
        )
if(OpenMP_CXX_FOUND)
    target_link_libraries(NASOQ PUBLIC OpenMP::OpenMP_CXX)
endif()
