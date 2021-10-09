include(catch2)
include(nlohmann_json)

add_executable(NASOQ_Catch2_Tests
    "main.cpp"
	"test_1.cpp"
)
target_include_directories(NASOQ_Catch2_Tests PRIVATE
    "${PROJECT_SOURCE_DIR}/include/nasoq/")
target_link_libraries(NASOQ_Catch2_Tests PRIVATE
    nasoq::eigen_interface
    Catch2::Catch2
    nlohmann_json::nlohmann_json
)

if(NASOQ_USE_BLAS_MKL AND NASOQ_USE_CLAPACK)
    # compares the CLAPACK implementation of the needed LAPACKE functions to the MKL ones
    target_sources(NASOQ_Catch2_Tests PRIVATE "test_clapacke.cpp")
    target_link_libraries(NASOQ_Catch2_Tests PRIVATE mkl::mkl)
endif()

target_compile_definitions(NASOQ_Catch2_Tests PRIVATE -DTEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_target_properties(NASOQ_Catch2_Tests PROPERTIES FOLDER tests)
