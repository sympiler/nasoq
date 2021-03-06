cmake_minimum_required(VERSION 3.1)

#add_subdirectory(nasoq)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

    if(NOT DEFINED NASOQ_ROOT)
        set(NASOQ_ROOT_DIR "../") #assumed the example folder is within nasoq
    else()
        set(NASOQ_ROOT_DIR "${NASOQ_ROOT}")
    endif()
set(NASOQ_INC ${NASOQ_ROOT_DIR}/include)
set(NASOQ_LIB ${NASOQ_ROOT_DIR}/lib/)

include_directories(${NASOQ_INC} ${METIS_INCLUDES} ${MKL_INCLUDE_DIR})
link_directories(${NASOQ_LIB})

### NASOQ example
project(NASOQ_Test)
set(CMAKE_CXX_STANDARD 11)

# TODO: remove MKL dependency using OpenBLAS submodule
find_package(MKL REQUIRED)



set(SOURCE_FILES_NAS nasoq_main.cpp)
add_executable(NASOQ_Test ${SOURCE_FILES_NAS} nasoq_main.cpp)
target_link_libraries(NASOQ_Test PRIVATE
        nasoq
        ${MKL_LIBRARIES} ${SUITESPARSE_LIBRARIES}
        ${METIS_LIBRARY} ${METIS_LIBRARIES}        )
if(OpenMP_CXX_FOUND)
    target_link_libraries(NASOQ_Test PRIVATE OpenMP::OpenMP_CXX)
endif()


set(SOURCE_FILES_NAS_ST nasoq_step_main.cpp)
add_executable(NASOQ_Step_Test ${SOURCE_FILES_NAS_ST} )
target_link_libraries(NASOQ_Step_Test PRIVATE
        nasoq
        ${MKL_LIBRARIES} ${SUITESPARSE_LIBRARIES}
        ${METIS_LIBRARY} ${METIS_LIBRARIES}        )
if(OpenMP_CXX_FOUND)
    target_link_libraries(NASOQ_Step_Test PRIVATE OpenMP::OpenMP_CXX)
endif()



### LBL example
project(LBL_Test)
set(CMAKE_CXX_STANDARD 11)

# TODO: remove MKL dependency using OpenBLAS submodule
find_package(MKL REQUIRED)



set(SOURCE_FILES_NAS lbl_main.cpp)
add_executable(LBL_Test ${SOURCE_FILES_NAS} lbl_main.cpp)
target_link_libraries(LBL_Test PRIVATE
        nasoq
        ${MKL_LIBRARIES} ${SUITESPARSE_LIBRARIES}
        ${METIS_LIBRARY} ${METIS_LIBRARIES}        )
if(OpenMP_CXX_FOUND)
    target_link_libraries(LBL_Test PRIVATE OpenMP::OpenMP_CXX)
endif()



### LBL SOMOD
project(LBLSOMOD_Test)
set(CMAKE_CXX_STANDARD 11)

# TODO: remove MKL dependency using OpenBLAS submodule
find_package(MKL REQUIRED)



set(SOURCE_FILES_NASOM lbl_somod_main.cpp)
add_executable(LBLSOMOD_Test ${SOURCE_FILES_NASOM} lbl_somod_main.cpp)
target_link_libraries(LBLSOMOD_Test PRIVATE
        nasoq
        ${MKL_LIBRARIES} ${SUITESPARSE_LIBRARIES}
        ${METIS_LIBRARY} ${METIS_LIBRARIES}        )
if(OpenMP_CXX_FOUND)
    target_link_libraries(LBLSOMOD_Test PRIVATE OpenMP::OpenMP_CXX)
endif()
