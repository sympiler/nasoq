project(NASOQ_Eigen)

set(SOURCE_FILES_EIG eigen_nasoq_driver.cpp)
add_executable(NASOQ_Eigen ${SOURCE_FILES_EIG}  eigen_nasoq_driver.cpp)
target_link_libraries(NASOQ_Eigen
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
        )