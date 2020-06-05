cmake_minimum_required(VERSION 3.9)

project(smp-convertor)
set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES smp_convertor.cpp)
add_executable(smp-convertor ${SOURCE_FILES} smp_convertor.cpp)


add_subdirectory(old)