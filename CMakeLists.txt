cmake_minimum_required(VERSION 3.16)

project(nasoq CXX C)

#----------------------------[ global configuration ]--------------------------#

set(CMAKE_CXX_STANDARD 11)

get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(NASOQ_IS_TOPLEVEL OFF)
else()
    set(NASOQ_IS_TOPLEVEL ON)
endif()

if(UNIX )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(NASOQ_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

# TODO: maybe install rules could be a more idiomatic way of doing this?
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/lib/")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/lib/")

#----------------------------------[ options ]---------------------------------#

option(NASOQ_BUILD_CLI "Build NASOQ Command-Line Interface EXEs" ${NASOQ_IS_TOPLEVEL})
option(NASOQ_BUILD_EXAMPLES "Build NASOQ Examples" ${NASOQ_IS_TOPLEVEL})
option(NASOQ_BUILD_TESTS "Build NASOQ Tests" ${NASOQ_IS_TOPLEVEL})
option(NASOQ_BUILD_DOCS "Build NASOQ Documentation" ${NASOQ_IS_TOPLEVEL})
option(NASOQ_WITH_EIGEN "Build NASOQ Eigen interface" ON)
option(NASOQ_USE_CLAPACK "Use CLAPACK as the LAPACK implementaiton" OFF)

set(NASOQ_BLAS_BACKEND "MKL" CACHE STRING "BLAS implementation for NASOQ to use")
set_property(CACHE NASOQ_BLAS_BACKEND PROPERTY STRINGS MKL OpenBLAS)

if(${NASOQ_BLAS_BACKEND} STREQUAL "MKL")
    set(NASOQ_USE_BLAS_MKL      ON)
    set(NASOQ_USE_BLAS_OpenBLAS OFF)
elseif(${NASOQ_BLAS_BACKEND} STREQUAL "OpenBLAS")
    set(NASOQ_USE_BLAS_MKL      OFF)
    set(NASOQ_USE_BLAS_OpenBLAS ON)
else()
    message(FATAL_ERROR "unrecognized value for `NASOQ_BLAS_BACKEND` option: '${NASOQ_BLAS_BACKEND}'")
endif()

if(NASOQ_USE_CLAPACK)
    if((NOT NASOQ_USE_BLAS_OpenBLAS) AND (NOT NASOQ_IS_TOPLEVEL))
        message(FATAL_ERROR "You can only set `NASOQ_USE_CLAPACK` to `ON` when `NASOQ_BLAS_BACKEND` is `OpenBLAS`")
    endif()
endif()

#--------------------------------[ third party ]-------------------------------#

list(PREPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_LIST_DIR}/cmake"
    "${CMAKE_CURRENT_LIST_DIR}/cmake/third_party"
)


if(NASOQ_WITH_EIGEN)
    include(eigen)
endif()

if(NASOQ_USE_BLAS_MKL)
    include(mkl)
elseif(NASOQ_USE_BLAS_OpenBLAS)
    if(NASOQ_USE_CLAPACK)
        set(openblas_WITHOUT_LAPACK ON)
    endif()

    include(openblas)

    if(openblas_WITHOUT_LAPACK AND NOT NASOQ_USE_CLAPACK)
        message(FATAL_ERROR "cannot build LAPACKE for use with OpenBLAS (maybe you don't have a Fortran compiler?) try setting `NASOQ_USE_CLAPACK` to `ON` to use a C-language alternative.")
    endif()
endif()

if(NASOQ_USE_CLAPACK)
    # This is outside of the `NASOQ_USE_BLAS_OpenBLAS` block since there is
    # unit test that compares MLK-LAPACKE to NASOQ's LAPACKE-style wrappers
    # around a couple of LAPACKE functions.
    include(clapack)
endif()

include(metis)

# TODO: SuiteSparse

#-----------------------------------[ nasoq ]----------------------------------#

find_package(OpenMP)


set(NASOQ_INTERNAL_INCLUDE_DIRS
    "${CMAKE_CURRENT_LIST_DIR}/smp-format/"
)

add_library(nasoq "")

target_sources(nasoq PRIVATE
    "src/nasoq.cpp"
    "src/nasoq_step.cpp"
    "src/common/def.cpp"
    "src/common/transpose_unsym.cpp"
    "src/common/Util.cpp"
    "src/matrixVector/spmv_CSC.cpp"
    "src/common/Norm.cpp"
    "src/QP/qp_utils.cpp"
    "src/common/SparseUtils.cpp"
    "src/common/DFS.cpp"
    "src/common/PostOrder.cpp"
    "src/common/Etree.cpp"
    "src/common/Transpose.cpp"
    "src/symbolic/ColumnCount.cpp"
    "src/symbolic/supernode_detection.cpp"
    "src/symbolic/Partitioning.cpp"
    "src/common/TreeUtils.cpp"
    "src/symbolic/PostOrderSpliting.cpp"
    "src/symbolic/InspectionLevel_06.cpp"
    "src/common/Reach.cpp"
    "src/symbolic/performanceModel.cpp"
    "src/symbolic/symbolic_phase.cpp"
    "src/common/Sym_BLAS.cpp"

    "src/ldl/Serial_blocked_ldl.cpp"
    "src/ldl/Serial_update_ldl.cpp"
    "src/ldl/Serial_blocked_ldl_02_2.cpp"
    "src/ldl/Serial_update_ldl_static.cpp"
    "src/ldl/serial_simplicial_ldl.cpp"
    "src/ldl/Serial_update_simplicial_ldl.cpp"
    "src/triangularSolve/BLAS.cpp"
    "src/triangularSolve/Triangular_BCSC.cpp"
    "src/triangularSolve/Triangular_CSC.cpp"
    "src/linear_solver/solve_phase.cpp"
    "src/gmres/mgmres.cpp"
    "src/matrixMatrix/spmm.cpp"
    "src/ldl/ldlt_check.cpp"
    "src/QP/osqp_utils.cpp"
    "src/QP/updown_test.cpp"
    "src/QP/nasoq_utils.cpp"
    "src/QP/linear_solver_wrapper.cpp"
)

if(OpenMP_FOUND)
    target_sources(nasoq PRIVATE
        "src/ldl/parallel_blocked_ldlt.cpp"
        "src/ldl/parallel_blocked_ldlt_02.cpp"
        "src/ldl/parallel_blocked_ldlt_03.cpp"
        "src/ldl/Parallel_update_ldl_02_2.cpp"
        "src/ldl/Parallel_simplicial_ldl.cpp"
        "src/ldl/Parallel_update_simplicial.cpp"
    )
    target_compile_definitions(nasoq PRIVATE "OPENMP")
    target_link_libraries(nasoq PRIVATE OpenMP::OpenMP_CXX)
endif()

if(NASOQ_USE_CLAPACK)
    target_sources(nasoq PRIVATE
        "src/clapacke/clapacke_dlapmt.cpp"
        "src/clapacke/clapacke_dsytrf.cpp"
    )
    target_compile_definitions(nasoq PUBLIC "NASOQ_USE_CLAPACK")
    target_link_libraries(nasoq PRIVATE clapack::clapack)
endif()


target_include_directories(nasoq PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_include_directories(nasoq PRIVATE ${NASOQ_INTERNAL_INCLUDE_DIRS})
target_compile_definitions(nasoq PRIVATE "SYM_REMOV")
target_link_libraries(nasoq PRIVATE metis::metis)

if(NASOQ_USE_BLAS_MKL)
    target_link_libraries(nasoq PRIVATE mkl::mkl)
    target_compile_definitions(nasoq PRIVATE "MKL_BLAS")
elseif(NASOQ_USE_BLAS_OpenBLAS)
    target_link_libraries(nasoq PRIVATE OpenBLAS::OpenBLAS)
    target_compile_definitions(nasoq PRIVATE "OPENBLAS")
endif()

if(NASOQ_WITH_EIGEN)
    add_subdirectory(eigen_interface)
    target_link_libraries(nasoq PUBLIC Eigen3::Eigen)
endif()

add_library(nasoq::nasoq ALIAS nasoq)

#------------------------------------[ exes ]----------------------------------#

if(NASOQ_BUILD_CLI)
    add_executable(NASOQ-BIN "nasoq_driver.cpp")
    target_include_directories(NASOQ-BIN PRIVATE ${NASOQ_INTERNAL_INCLUDE_DIRS})
    target_link_libraries(NASOQ-BIN PRIVATE nasoq::nasoq)
    if(NASOQ_USE_BLAS_MKL)
        target_link_libraries(NASOQ-BIN PRIVATE mkl::mkl)
        target_compile_definitions(NASOQ-BIN PRIVATE "MKL_BLAS")
    elseif(NASOQ_USE_BLAS_OpenBLAS)
	    #        target_link_libraries(NASOQ-BIN PRIVATE OpenBLAS::OpenBLAS)
        target_compile_definitions(NASOQ-BIN PRIVATE "OPENBLAS")
    endif()
    if(OpenMP_FOUND)
        target_compile_definitions(NASOQ-BIN PRIVATE "OPENMP")
        target_link_libraries(NASOQ-BIN  PRIVATE OpenMP::OpenMP_CXX)
    endif()
    set_target_properties(NASOQ-BIN PROPERTIES FOLDER exes)
endif()


if(NASOQ_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(NASOQ_BUILD_TESTS)
    add_subdirectory(tests)
endif()

# =============================================================================
#   Add custom target dox_doxygen to generate doxygen documentation for the 
#   nasoq lib.
#
#   Dependencies: doxygen, graphviz.
# =============================================================================
if(NASOQ_BUILD_DOCS)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message("Doxygen build started")

        add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
endif(NASOQ_BUILD_DOCS)
