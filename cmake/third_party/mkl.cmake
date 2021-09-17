if(TARGET mkl::mkl)
    return()
endif()

find_package(MKL REQUIRED)

# add a target collecting all of the MKL libraries
add_library(mkl INTERFACE)
target_link_libraries(mkl INTERFACE ${MKL_LIBRARIES})
target_include_directories(mkl INTERFACE "${MKL_INCLUDE_DIR}")

add_library(mkl::mkl ALIAS mkl)
