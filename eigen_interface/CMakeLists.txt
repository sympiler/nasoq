project(NASOQ_Eigen)
set(SOURCE_FILES_EIG eigen_nasoq_driver.cpp)
add_executable(NASOQ_Eigen ${SOURCE_FILES_EIG}  eigen_nasoq_driver.cpp)
target_link_libraries(NASOQ_Eigen PRIVATE
        NASOQ
        )
if(OpenMP_CXX_FOUND)
    target_link_libraries(NASOQ_Eigen PRIVATE OpenMP::OpenMP_CXX)
endif()


project(LBL_Eigen)
set(SOURCE_FILES_LBL_EIG eigen_lbl_driver.cpp)
add_executable(LBL_Eigen ${SOURCE_FILES_LBL_EIG}  eigen_lbl_driver.cpp)
target_link_libraries(LBL_Eigen PRIVATE
        NASOQ
        )
if(OpenMP_CXX_FOUND)
    target_link_libraries(LBL_Eigen PRIVATE OpenMP::OpenMP_CXX)
endif()