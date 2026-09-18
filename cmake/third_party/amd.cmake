#
# Fill-reducing ordering: AMD (SuiteSparse), used when NASOQ_ORDERING is set
# to "AMD" instead of the default "METIS".
#
if(TARGET amd::amd)
    return()
endif()

# A parent project that already builds SuiteSparse's own `amd` CMake target
# (e.g. because it also links CHOLMOD, which depends on AMD) can just reuse
# that instead of finding/fetching a second copy.
if(TARGET amd)
    message(STATUS "Third-party (external): reusing existing target 'amd' as 'amd::amd'")
    add_library(amd::amd ALIAS amd)
    return()
endif()

find_package(AMD CONFIG QUIET)
if(TARGET SuiteSparse::AMD)
    message(STATUS "Third-party (external): creating target 'amd::amd' from system AMD")
    add_library(amd::amd ALIAS SuiteSparse::AMD)
else()
    message(STATUS "Third-party (internal): creating target 'amd::amd'")

    include(FetchContent)
    FetchContent_Declare(
        suitesparse_amd
        GIT_REPOSITORY https://github.com/DrTimothyAldenDavis/SuiteSparse
        GIT_TAG        v7.8.3
        GIT_SHALLOW    TRUE
    )
    FetchContent_GetProperties(suitesparse_amd)
    if(NOT suitesparse_amd_POPULATED)
        FetchContent_Populate(suitesparse_amd)
    endif()

    # AMD only needs a handful of SuiteSparse_config hooks (malloc/printf
    # wrappers), not the rest of SuiteSparse.
    file(GLOB AMD_SOURCES "${suitesparse_amd_SOURCE_DIR}/AMD/Source/*.c")
    set(SUITESPARSE_CONFIG_SOURCE "${suitesparse_amd_SOURCE_DIR}/SuiteSparse_config/SuiteSparse_config.c")

    add_library(amd STATIC ${AMD_SOURCES} ${SUITESPARSE_CONFIG_SOURCE})
    target_include_directories(amd PRIVATE
        "${suitesparse_amd_SOURCE_DIR}/AMD/Include"
        "${suitesparse_amd_SOURCE_DIR}/SuiteSparse_config")
    target_include_directories(amd SYSTEM PUBLIC
        "$<BUILD_INTERFACE:${suitesparse_amd_SOURCE_DIR}/AMD/Include>"
        "$<BUILD_INTERFACE:${suitesparse_amd_SOURCE_DIR}/SuiteSparse_config>")
    # NTIMER: skip SuiteSparse_config's platform timer code, which otherwise
    # needs extra platform-specific libraries (e.g. -lrt on Linux) we'd have
    # to detect/link ourselves; NASOQ doesn't use AMD's timing output.
    target_compile_definitions(amd PRIVATE NTIMER)
    set_target_properties(amd PROPERTIES FOLDER third_party)

    add_library(amd::amd ALIAS amd)
endif()
