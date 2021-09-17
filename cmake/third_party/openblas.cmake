if(TARGET OpenBLAS::OpenBLAS)
    return()
endif()


message(STATUS "Third-party (external): creating target 'OpenBLAS::OpenBLAS'")


if(NOT DEFINED openblas_WITHOUT_LAPACK)
    include(CheckLanguage)
    check_language(Fortran)
    if(NOT CMAKE_Fortran_COMPILER)
        set(openblas_WITHOUT_LAPACK ON)
    else()
        set(openblas_WITHOUT_LAPACK OFF)
    endif()
endif()


include(FetchContent)
FetchContent_Declare(
    openblas
    GIT_REPOSITORY https://github.com/xianyi/OpenBLAS
    GIT_TAG v0.3.13
    GIT_SHALLOW TRUE
)

FetchContent_GetProperties(openblas)
if(NOT openblas_POPULATED)
    FetchContent_Populate(openblas)
endif()


if(NOT EXISTS "${openblas_BINARY_DIR}/CMakeCache.txt")
    # run cmake to create the project files
    execute_process(
        COMMAND
            ${CMAKE_COMMAND} "${openblas_SOURCE_DIR}" "-G" "${CMAKE_GENERATOR}"
            "-DCMAKE_INSTALL_PREFIX=inst"
            "-DCMAKE_CXX_FLAGS_DEBUG=\"${CMAKE_CXX_FLAGS_DEBUG}\""
            "-DCMAKE_CXX_FLAGS_RELEASE=\"${CMAKE_CXX_FLAGS_RELEASE}\""
            "-DUSE_THREAD=OFF"
            "-DUSE_LOCKING=ON"
            "-DBUILD_WITHOUT_LAPACK=${openblas_WITHOUT_LAPACK}"
        WORKING_DIRECTORY "${openblas_BINARY_DIR}"
    )
endif()


add_library(OpenBLAS IMPORTED INTERFACE GLOBAL)

foreach(CONFIG_NAME IN LISTS CMAKE_CONFIGURATION_TYPES)
    string(TOUPPER "${CONFIG_NAME}" CONFIG_NAME_TOUPPER)
    set(OPENBLAS_OUTPUT_FILE "${openblas_BINARY_DIR}/lib/${CONFIG_NAME_TOUPPER}")
    set(OPENBLAS_INSTALL_DIR "${openblas_BINARY_DIR}/inst/${CONFIG_NAME}")

    if(NOT EXISTS "${OPENBLAS_OUTPUT_FILE}")
        execute_process(
            COMMAND ${CMAKE_COMMAND} "--build" "." "--config" "${CONFIG_NAME}"
            WORKING_DIRECTORY "${openblas_BINARY_DIR}"
        )
    endif()

    if(NOT EXISTS "${OPENBLAS_INSTALL_DIR}")
        file(MAKE_DIRECTORY "${openblas_BINARY_DIR}/inst")
    endif()

    if(NOT EXISTS "${OPENBLAS_INSTALL_DIR}/include/openblas/openblas_config.h")
        execute_process(
            COMMAND ${CMAKE_COMMAND} "--install" "." "--config" "${CONFIG_NAME}"
                "--prefix" "${OPENBLAS_INSTALL_DIR}"
                WORKING_DIRECTORY "${openblas_BINARY_DIR}"
        )
    endif()

    add_library(OpenBLAS_${CONFIG_NAME} STATIC IMPORTED GLOBAL)
    set_property(
        TARGET OpenBLAS_${CONFIG_NAME}
        PROPERTY IMPORTED_LOCATION
        "${openblas_BINARY_DIR}/inst/${CONFIG_NAME}/lib/openblas.lib"
    )
    target_include_directories(OpenBLAS_${CONFIG_NAME} INTERFACE "${openblas_BINARY_DIR}/inst/${CONFIG_NAME}/include")
    target_link_libraries(OpenBLAS INTERFACE $<$<CONFIG:${CONFIG_NAME}>:OpenBLAS_${CONFIG_NAME}>)
endforeach()

add_library(OpenBLAS::OpenBLAS ALIAS OpenBLAS)
