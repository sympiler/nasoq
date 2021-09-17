include(Catch2)
include(nlohmann_json)

file(GLOB SRC_FILES "*.cpp" "*.h")
add_executable(NASOQ_Catch2_Tests ${SRC_FILES})
target_link_libraries(NASOQ_Catch2_Tests PRIVATE
    nasoq::eigen_interface
    Catch2::Catch2
    nlohmann_json::nlohmann_json
)

target_compile_definitions(NASOQ_Catch2_Tests PRIVATE -DTEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_target_properties(NASOQ_Catch2_Tests PROPERTIES FOLDER tests)
