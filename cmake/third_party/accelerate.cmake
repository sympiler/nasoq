if(TARGET BLAS::BLAS)
    return()
endif()

if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "arm64" OR "${CMAKE_OSX_ARCHITECTURES}" MATCHES "arm64")
    # Use Accelerate on macOS M1
    set(BLA_VENDOR Apple)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
else()
    message(FATAL_ERROR "Accelerate is only support on apple M1")
endif()
