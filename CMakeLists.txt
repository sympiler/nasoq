cmake_minimum_required(VERSION 3.9)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)


project(libnasoq_static)

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

find_package(OpenMP)

find_package (Eigen3 OPTIONAL_COMPONENTS NO_MODULE)

if(UNIX )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 ")
endif()

#set(NASOQ_LIB_INCLUDE_DIR
#        ${MKL_INCLUDE_DIR}
#        ${SUITESPARSE_INCLUDE_DIRS}
#        ${METIS_INCLUDE_DIR}
#        ${METIS_INCLUDES}
#        "${CMAKE_CURRENT_LIST_DIR}/symbolic/"
#        "${CMAKE_CURRENT_LIST_DIR}/common/"
#        "${CMAKE_CURRENT_LIST_DIR}/ldl/"
#        "${CMAKE_CURRENT_LIST_DIR}/matrixMatrix/"
#        "${CMAKE_CURRENT_LIST_DIR}/matrixVector/"
#        "${CMAKE_CURRENT_LIST_DIR}/linear_solver/"
#        "${CMAKE_CURRENT_LIST_DIR}/gmres/"
#        "${CMAKE_CURRENT_LIST_DIR}/QP/"
#        "${CMAKE_CURRENT_LIST_DIR}/triangularSolve/"
#        "${CMAKE_CURRENT_LIST_DIR}/smp-format/")

set(NASOQ_LIB_INCLUDE_DIR
        ${MKL_INCLUDE_DIR}
        ${SUITESPARSE_INCLUDE_DIRS}
        ${METIS_INCLUDE_DIR}
        ${METIS_INCLUDES}
        "${CMAKE_CURRENT_LIST_DIR}/smp-format/"
        "${CMAKE_CURRENT_LIST_DIR}/include/")

#set(NASOQ_INCLUDE_DIR
#        ${MKL_INCLUDE_DIR}
#        "${CMAKE_CURRENT_LIST_DIR}/symbolic/"
#        "${CMAKE_CURRENT_LIST_DIR}/common/"
#        "${CMAKE_CURRENT_LIST_DIR}/ldl/"
#        "${CMAKE_CURRENT_LIST_DIR}/matrixMatrix/"
#        "${CMAKE_CURRENT_LIST_DIR}/matrixVector/"
#        "${CMAKE_CURRENT_LIST_DIR}/linear_solver/"
#        "${CMAKE_CURRENT_LIST_DIR}/gmres/"
#        "${CMAKE_CURRENT_LIST_DIR}/QP/"
#        "${CMAKE_CURRENT_LIST_DIR}/triangularSolve/"
#        "${CMAKE_CURRENT_LIST_DIR}/smp-format/")

set(NASOQ_INCLUDE_DIR
        ${MKL_INCLUDE_DIR}
        include/
        "${CMAKE_CURRENT_LIST_DIR}/smp-format/"
        "${CMAKE_CURRENT_LIST_DIR}/include/")

set(NASOQ_LIB_DIR
        ${MKL_LIBRARIES}
        ${SUITESPARSE_LIBRARIES}
        ${METIS_LIBRARY}
        ${METIS_LIBRARIES})

set(SOURCE_FILES
        src/nasoq.cpp
        src/nasoq_step.cpp
        src/common/def.cpp
        src/common/transpose_unsym.cpp
        src/common/Util.cpp
        src/matrixVector/spmv_CSC.cpp
        src/common/Norm.cpp
        src/QP/qp_utils.cpp
        src/common/SparseUtils.cpp
        src/common/DFS.cpp
        src/common/PostOrder.cpp
        src/common/Etree.cpp
        src/common/Transpose.cpp
        src/symbolic/ColumnCount.cpp
        src/symbolic/supernode_detection.cpp
        src/symbolic/Partitioning.cpp
        src/common/TreeUtils.cpp
        src/symbolic/PostOrderSpliting.cpp
        src/symbolic/InspectionLevel_06.cpp
        src/common/Reach.cpp
        src/symbolic/performanceModel.cpp
        src/symbolic/symbolic_phase.cpp
        src/common/Sym_BLAS.cpp
        src/ldl/parallel_blocked_ldlt.cpp
        src/ldl/parallel_blocked_ldlt_02.cpp
        src/ldl/parallel_blocked_ldlt_03.cpp
        src/ldl/Serial_blocked_ldl.cpp
        src/ldl/Serial_update_ldl.cpp
        src/ldl/Serial_blocked_ldl_02_2.cpp
        src/ldl/Serial_update_ldl_static.cpp
        src/ldl/Parallel_update_ldl_02_2.cpp
        src/ldl/serial_simplicial_ldl.cpp
        src/ldl/Serial_update_simplicial_ldl.cpp
        src/ldl/Parallel_simplicial_ldl.cpp
        src/ldl/Parallel_update_simplicial.cpp
        src/triangularSolve/BLAS.cpp
        src/triangularSolve/Triangular_BCSC.cpp
        src/triangularSolve/Triangular_CSC.cpp
        src/linear_solver/solve_phase.cpp
        src/gmres/mgmres.cpp
        src/matrixMatrix/spmm.cpp
        src/ldl/ldlt_check.cpp
        src/QP/osqp_utils.cpp
        src/QP/updown_test.cpp
        src/QP/nasoq_utils.cpp
        src/QP/linear_solver_wrapper.cpp
        )

message("OpenMP_CXX_LIBRARIES = " ${OpenMP_CXX_LIBRARIES})
message("OpenMP_CXX_INCLUDE_DIRS = " ${OpenMP_CXX_INCLUDE_DIRS})

add_library(nasoq_static STATIC ${SOURCE_FILES})
target_include_directories(nasoq_static PRIVATE ${NASOQ_LIB_INCLUDE_DIR} ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(nasoq_static ${NASOQ_LIB_DIR} ${OpenMP_CXX_LIBRARIES})
#set_target_properties(nasoq_static PROPERTIES
#        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
#set_target_properties(nasoq_static PROPERTIES
#        CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)

project(libnasoq_shared)
add_library(nasoq SHARED ${SOURCE_FILES})
target_include_directories(nasoq PRIVATE ${NASOQ_LIB_INCLUDE_DIR} ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(nasoq ${NASOQ_LIB_DIR} ${OpenMP_CXX_LIBRARIES})
set_target_properties(nasoq PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)


project(NASOQ-BIN)
add_executable(NASOQ-BIN nasoq_driver.cpp)
target_include_directories(NASOQ-BIN PRIVATE ${NASOQ_INCLUDE_DIR} ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(NASOQ-BIN nasoq_static ${OpenMP_CXX_LIBRARIES})

#target_link_libraries(NASOQ-BIN PRIVATE
#        ${NASOQ_LIB_DIR})
#target_include_directories(NASOQ-BIN PRIVATE ${NASOQ_INCLUDE_DIR} )
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(NASOQ-BIN PRIVATE OpenMP::OpenMP_CXX)
#endif()
#
#### Building the interface library
#project(NASOQ)
##file(GLOB_RECURSE SRC_FILES ${NASOQ_INCLUDE_DIR} "*.h" )
#add_library (NASOQ  INTERFACE )
#target_include_directories(NASOQ INTERFACE  ${NASOQ_INCLUDE_DIR} )
##target_sources(NASOQ ${SRC_FILES})
#target_link_libraries(NASOQ INTERFACE
#        ${NASOQ_LIB_DIR}
#        )
#
#if(TARGET Eigen3::Eigen)
#    add_subdirectory(eigen_interface)
#endif()

