if(TARGET BLAS::BLAS)
    return()
endif()

if(NOT APPLE)
    message(FATAL_ERROR "Accelerate is only available on Apple platforms")
endif()

set(BLA_VENDOR Apple)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
