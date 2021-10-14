if(TARGET OpenBLAS::OpenBLAS)
    return()
endif()


set(BLA_VENDOR OpenBLAS)
set(WINBLAS "")
#    set(BLA_STATIC TRUE)
find_package(BLAS QUIET)

message(STATUS "Third-party (external): creating target 'OpenBLAS::OpenBLAS'")
if(BLAS_FOUND)
    #TODO: this is a hack to include openblas headers
    include(FetchContent)
    FetchContent_Declare(
            kazlibbin
            GIT_REPOSITORY https://github.com/cheshmi/lib_binary.git
            GIT_SHALLOW TRUE
    )
    FetchContent_GetProperties(kazlibbin)
    if(NOT kazlibbin_POPULATED)
        FetchContent_Populate(kazlibbin)
    endif()
    set(BIN_LIBS "${kazlibbin_SOURCE_DIR}")
    if(UNIX AND NOT APPLE)
        set(OPENBLAS_DIR ${BIN_LIBS}/linux/)
    endif()
    if(UNIX AND APPLE)
        set(OPENBLAS_DIR ${BIN_LIBS}/mac/)
    endif()
    add_library(OpenBLAS INTERFACE)
    target_include_directories(OpenBLAS INTERFACE "${OPENBLAS_DIR}")
    #include_directories( "${OPENBLAS_DIR}")
    if(UNIX AND NOT APPLE)
	    target_link_libraries(OpenBLAS INTERFACE ${BLAS_LIBRARIES} lapacke)
    endif()
    if(UNIX AND APPLE)
	    target_link_libraries(OpenBLAS INTERFACE ${BLAS_LIBRARIES})
    endif()
    add_library(OpenBLAS::OpenBLAS ALIAS OpenBLAS)
else()
    if(NOT DEFINED openblas_WITHOUT_LAPACK)
        include(CheckLanguage)
        check_language(Fortran)
        if(NOT CMAKE_Fortran_COMPILER)
            set(openblas_WITHOUT_LAPACK ON)
        else()
            set(openblas_WITHOUT_LAPACK OFF)
        endif()
    endif()
    set(OPENBLAS_MOCK_SOURCES "") # overriden later in certain cases

    include(FetchContent)
    FetchContent_Declare(
        openblas
        GIT_REPOSITORY https://github.com/xianyi/OpenBLAS
        GIT_TAG v0.3.18
        GIT_SHALLOW TRUE
    )

    FetchContent_GetProperties(openblas)
    if(NOT openblas_POPULATED)
        FetchContent_Populate(openblas)

        if(openblas_WITHOUT_LAPACK)
            if(CMAKE_GENERATOR STREQUAL "Xcode")
                # This hack deserves an explanation.  When compiling OpenBLAS
                # without LAPACK, the `add_library(openblas ...)` call will
                # list only target object files and no sources, which causes
                # Xcode not to generate `ibopenblas.a`.  To fix, we create
                # an empty source file to include in the `add_library(openblas ...)`
                # call, which is enough to satisfy Xcode.  To actually *get* this
                # source file name to the `add_library` call, we use another hack
                # and set it (otherwise empty) `LAPACKE_SOURCES` variable in
                # OpenBLAS, which happens to be included where we need it.
                file(TOUCH "${openblas_SOURCE_DIR}/openblas_mock_source_file.c")
                set(OPENBLAS_MOCK_SOURCES "openblas_mock_source_file.c")
            endif()
        endif()
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
                "-DBUILD_DOUBLE=ON"
                "-DBUILD_WITHOUT_LAPACK=${openblas_WITHOUT_LAPACK}"
                "-DNOFORTRAN=${openblas_WITHOUT_LAPACK}"
                "-DNO_LAPACK=${openblas_WITHOUT_LAPACK}"
                "-DLAPACKE_SOURCES=${OPENBLAS_MOCK_SOURCES}" # hack needed to make Xcode happy
            WORKING_DIRECTORY "${openblas_BINARY_DIR}"
        )
    endif()

    # We need to handle this whether or not the Cmake generator supports
    # multiple configurations.
    if(DEFINED CMAKE_CONFIGURATION_TYPES)
        set(NASOQ_OPENBLAS_CONFIG_TYPES ${CMAKE_CONFIGURATION_TYPES})
        set(NASOQ_OPENBLAS_USE_CONFIG_TYPES ON)
    else()
        set(NASOQ_OPENBLAS_CONFIG_TYPES "none")
        set(NASOQ_OPENBLAS_USE_CONFIG_TYPES OFF)
    endif()

    if(NASOQ_OPENBLAS_USE_CONFIG_TYPES)
        add_library(OpenBLAS IMPORTED INTERFACE GLOBAL)
    endif()

    foreach(CONFIG_NAME IN LISTS NASOQ_OPENBLAS_CONFIG_TYPES)
        if(NASOQ_OPENBLAS_USE_CONFIG_TYPES)
            string(TOUPPER "${CONFIG_NAME}" CONFIG_NAME_TOUPPER)
            set(OPENBLAS_INSTALL_DIR "${openblas_BINARY_DIR}/inst/${CONFIG_NAME}")
            set(OPENBLAS_LIB_DIR "${openblas_BINARY_DIR}/inst/${CONFIG_NAME}/lib")
            set(OPENBLAS_INCLUDE_DIR "${openblas_BINARY_DIR}/inst/${CONFIG_NAME}/include")
            set(OPENBLAS_TARGET_NAME OpenBLAS_${CONFIG_NAME})
            set(OPENBLAS_CMAKE_CONFIG_CMD "--config" "${CONFIG_NAME}")
        else()
            set(OPENBLAS_INSTALL_DIR "${openblas_BINARY_DIR}/inst")
            set(OPENBLAS_LIB_DIR "${openblas_BINARY_DIR}/inst/lib")
            set(OPENBLAS_INCLUDE_DIR "${openblas_BINARY_DIR}/inst/include")
            set(OPENBLAS_TARGET_NAME OpenBLAS)
            set(OPENBLAS_CMAKE_CONFIG_CMD "")
        endif()

        if(MSVC)
            set(OPENBLAS_LIB_PATH "${OPENBLAS_LIB_DIR}/openblas.lib")
        else()
            set(OPENBLAS_LIB_PATH "${OPENBLAS_LIB_DIR}/libopenblas.a")
        endif()

        if(NOT EXISTS "${OPENBLAS_LIB_PATH}")
            execute_process(
                COMMAND ${CMAKE_COMMAND} "--build" "." ${OPENBLAS_CMAKE_CONFIG_CMD}
                WORKING_DIRECTORY "${openblas_BINARY_DIR}"
            )
        endif()

        if(NOT EXISTS "${OPENBLAS_INSTALL_DIR}")
            file(MAKE_DIRECTORY "${openblas_BINARY_DIR}/inst")
        endif()

        if(NOT EXISTS "${OPENBLAS_INSTALL_DIR}/include/openblas/openblas_config.h")
            execute_process(
                COMMAND ${CMAKE_COMMAND} "--install" "." ${OPENBLAS_CMAKE_CONFIG_CMD}
                    "--prefix" "${OPENBLAS_INSTALL_DIR}"
                    WORKING_DIRECTORY "${openblas_BINARY_DIR}"
            )
        endif()

        add_library(${OPENBLAS_TARGET_NAME} STATIC IMPORTED GLOBAL)
        set_property(
            TARGET ${OPENBLAS_TARGET_NAME}
            PROPERTY IMPORTED_LOCATION
            ${OPENBLAS_LIB_PATH}
        )
        target_include_directories(${OPENBLAS_TARGET_NAME} INTERFACE "${OPENBLAS_INCLUDE_DIR}")
        
        if(NASOQ_OPENBLAS_USE_CONFIG_TYPES)
            target_link_libraries(OpenBLAS INTERFACE $<$<CONFIG:${CONFIG_NAME}>:OpenBLAS_${CONFIG_NAME}>)
        endif()
    endforeach()

    add_library(OpenBLAS::OpenBLAS ALIAS OpenBLAS)

endif()
