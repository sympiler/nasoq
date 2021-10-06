
add_library(nasoq_eigen INTERFACE)
target_include_directories(nasoq_eigen INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(nasoq_eigen INTERFACE
	nasoq::nasoq
	Eigen3::Eigen
)

add_library(nasoq::eigen_interface ALIAS nasoq_eigen)

#------------------------------------------------------------------------------#

if(NASOQ_BUILD_CLI)
	add_executable(NASOQ_Eigen "eigen_nasoq_driver.cpp")
	target_include_directories(NASOQ_Eigen PRIVATE "${PROJECT_SOURCE_DIR}/include")
	target_include_directories(NASOQ_Eigen PRIVATE
			"${PROJECT_SOURCE_DIR}/include/nasoq/")
	target_link_libraries(NASOQ_Eigen PRIVATE nasoq::eigen_interface)
	set_target_properties(NASOQ_Eigen PROPERTIES FOLDER exes)

	add_executable(LBL_Eigen "eigen_lbl_driver.cpp")
	target_include_directories(LBL_Eigen PRIVATE "${PROJECT_SOURCE_DIR}/include")
	target_include_directories(LBL_Eigen PRIVATE
			"${PROJECT_SOURCE_DIR}/include/nasoq")
	target_link_libraries(LBL_Eigen PRIVATE nasoq::eigen_interface)
	set_target_properties(LBL_Eigen PROPERTIES FOLDER exes)
endif()
