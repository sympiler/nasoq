cmake_minimum_required(VERSION 2.7)

message("Set to: ${MKLROOT}")
set(MKL_INC "$ENV{MKLROOT}/include")

if(APPLE)
    SET(CMAKE_C_COMPILER /usr/local/Cellar/gcc/9.1.0/bin/gcc-9)
    SET(CMAKE_CXX_COMPILER /usr/local/Cellar/gcc/9.1.0/bin/g++-9)
    set(MKL_LIB "$ENV{MKLROOT}/lib")
    set(MKL_LIB2 "$ENV{MKLROOT}/../compiler/lib")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
endif()
if(UNIX AND NOT APPLE)
    set(MKL_LIB "$ENV{MKLROOT}/lib/intel64")
    # for iomp lib
    set(MKL_LIB2 "$ENV{MKLROOT}/../compiler/lib/intel64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()


#For installing Suitesparse using MKL BLAS. First -liomp tu Suitesparse_config.mk
# Then make sure source ${MKLROOT}/bin/mklvars.sh intel64 is ran
# then make
set(SUITE_INC "$ENV{SUITEROOT}/include")
set(SUITE_LIB "$ENV{SUITEROOT}/lib")

set(METIS_INC "$ENV{METISROOT}/include")
set(METIS_LIB "$ENV{METISROOT}/libmetis")



set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}    -std=c++11 -m64 ")

project(NASOQ)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL  ")

set(SOURCE_FILES nasoq_driver.cpp)
include_directories(${MKL_INC} ${SUITE_INC} ${METIS_INC}
        "symbolic/" "common/" "ldl/" "matrixMatrix/"
        "matrixVector/"  "linear_solver/" "gmres/" "QP/"
        "triangularSolve/")
link_directories(${MKL_LIB} ${MKL_LIB2} ${SUITE_LIB} ${METIS_LIB} )
add_executable(NASOQ ${SOURCE_FILES}  nasoq_driver.cpp)
target_link_libraries(NASOQ
        #gomp
        iomp5
        #mkl_sequential
        mkl_intel_lp64
        mkl_core
        #mkl_gnu_thread
        mkl_intel_thread
        pthread
        m
        dl
        amd
        cholmod
        suitesparseconfig
        metis
        #pastix_murge
        #pastix
        #matrix_driver
        #scotch
        #scotcherrexit
        )

#add_subdirectory(eigen_interface)