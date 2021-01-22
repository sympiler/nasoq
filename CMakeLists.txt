cmake_minimum_required(VERSION 3.9)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

#include(ExternalProject)
include(DownloadProject)

set(CMAKE_CXX_STANDARD 11)
set(NASOQ_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(NASOQ_EXTERNAL ${CMAKE_CURRENT_LIST_DIR}/external)


project(libnasoq_static)



set(OPENBLAS_STAT_LIB "")
set(OPENBLAS_DYN_LIB  "")
set(MKL_FOUND FALSE)
find_package(MKL REQUIRED)
if(NOT MKL_FOUND)
    #set(BLA_VENDOR OpenBLAS)
    set(BLA_STATIC TRUE)
    find_package(BLAS)
        message("++++> ${BLAS_LIBRARIES}")
        if(NOT BLAS_FOUND)
            download_project(PROJ openblas
                    GIT_REPOSITORY https://github.com/xianyi/OpenBLAS.git
                    PREFIX  ${NASOQ_EXTERNAL}/openblas/
                    #BINARY_DIR ${NASOQ_ROOT_DIR}/lib/
## TODO find a way to set the install directory here. Header files will be there
                    UPDATE_DISCONNECTED 1
                    )
            add_subdirectory(${openblas_SOURCE_DIR})
            set(OPENBLAS_INCLUDE_DIR ${openblas_SOURCE_DIR}/lapack-netlib/CBLAS/ )# ${openblas_SOURCE_DIR}/lapack-netlib/CBLAS/include/) # ${openblas_SOURCE_DIR}/lapack-netlib/LAPACKE/include/)
            set(OPENBLAS_LIBRARY_DIR ${openblas_SOURCE_DIR}/lapack-netlib/)
            set(OPENBLAS_STAT_LIB openblas)
            set(OPENBLAS_DYN_LIB openblas_d)
            message("======> ${OPENBLAS_INCLUDE_DIR}")
        else() #
            set(OPENBLAS_STAT_LIB ${BLAS_LIBRARIES})
            message(STATUS "OpenBLAS found")
        endif()
            add_definitions(-DOPENBLAS)
    else()
    add_definitions(-DMKL_BLAS)
endif()


# Adding hints for suitesparse
set(SUITESPARSE_INCLUDE_DIR_HINTS ${SUITESPARSE_INCLUDE_DIR_HINTS} ${SUITE_ROOT_PATH}/include)
set(SUITESPARSE_LIBRARY_DIR_HINTS ${SUITESPARSE_LIBRARY_DIR_HINTS} ${SUITE_ROOT_PATH}/lib)
set(BLA_STATIC TRUE)
find_package(SuiteSparse OPTIONAL_COMPONENTS)


#IF(NOT SUITESPARSE_FOUND)
   # set( INCLUDE_INSTALL_DIR ${METIS_ROOT_PATH}/include)
   # set( LIB_INSTALL_DIR ${METIS_ROOT_PATH}/libmetis)
    find_package(METIS OPTIONAL_COMPONENTS)
   # set(METIS_LIBRARY "")
   # set(METIS_INCLUDE_DIR "")
#ELSE()
#    set(ignoreMe "${METIS_ROOT_PATH}")
#ENDIF()
if(METIS_VERSION_OK)
#    find_program(MAKE_EXE NAMES gmake nmake make)
#    set(metis_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/metis")
#    set(metis_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${metis_INSTALL_DIR})
#set(EXTERNAL_INSTALL_LOCATION ${CMAKE_CURRENT_LIST_DIR}/external)
   # set(GKLIB_PATH ${NASOQ_EXTERNAL}/metis/metis-src/GKLIB/ CACHE PATH "d")
   download_project(PROJ metis
      GIT_REPOSITORY https://github.com/cheshmi/METIS.git  # https://github.com/KarypisLab/METIS.git
          # BUILD_IN_SOURCE 1
          #PREFIX  ${metis_BINARY_DIR}
      #     URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
           PREFIX  ${NASOQ_EXTERNAL}/metis/
           #SOURCE_DIR ${NASOQ_EXTERNAL}/metis/metis-src/
           BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/metis/
          UPDATE_DISCONNECTED 1
          #CMAKE_ARGS  -DGKLIB_PATH=${metis_SOURCE_DIR}/GKlib/
          # BUILD_COMMAND ${MAKE_EXE}
      )
    add_subdirectory(${metis_SOURCE_DIR})
    set(METIS_INCLUDE_DIR ${metis_SOURCE_DIR}/include)
    set(METIS_INCLUDES ${metis_SOURCE_DIR}/include/)
    set(METIS_LIBRARY ${NASOQ_ROOT_DIR}/lib/)
    set(METIS_LIBRARIES ${NASOQ_ROOT_DIR}/lib/)
endif()


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
        ${OPENBLAS_INCLUDE_DIR}
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
        ${NASOQ_ROOT_DIR}/lib/
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
add_dependencies(nasoq_static metis)
if (NOT MKL_FOUND AND NOT BLAS_FOUND)
add_dependencies(nasoq_static openblas_static)
endif ()
#target_include_directories(nasoq_static PRIVATE ${NASOQ_LIB_INCLUDE_DIR})
#target_link_libraries(nasoq_static metis ${NASOQ_LIB_DIR})
target_include_directories(nasoq_static PRIVATE ${NASOQ_LIB_INCLUDE_DIR} ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(nasoq_static metis ${NASOQ_LIB_DIR}
        ${OpenMP_CXX_LIBRARIES}
        ${OPENBLAS_STAT_LIB})
#set_target_properties(nasoq_static PROPERTIES
#        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
#set_target_properties(nasoq_static PROPERTIES
#        CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)

project(libnasoq_shared)
add_library(nasoq SHARED ${SOURCE_FILES})
add_dependencies(nasoq metis)
#target_include_directories(nasoq PRIVATE ${NASOQ_LIB_INCLUDE_DIR})
#target_link_libraries(nasoq ${NASOQ_LIB_DIR})
if (NOT MKL_FOUND AND NOT BLAS_FOUND)
    add_dependencies(nasoq openblas_static)
endif ()
target_include_directories(nasoq PRIVATE ${NASOQ_LIB_INCLUDE_DIR} ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(nasoq metis ${NASOQ_LIB_DIR} ${OpenMP_CXX_LIBRARIES} ${OPENBLAS_STAT_LIB})
set_target_properties(nasoq PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/   )


project(NASOQ-BIN)
add_executable(NASOQ-BIN nasoq_driver.cpp)
target_include_directories(NASOQ-BIN PRIVATE ${NASOQ_INCLUDE_DIR})
link_directories(NASOQ-BIN PRIVATE
        ${NASOQ_LIB_DIR}
        #${NASOQ_ROOT_DIR}/lib
        )
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

