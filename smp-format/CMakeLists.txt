cmake_minimum_required(VERSION 3.5)

project(smp-convertor)
set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES smp_convertor.cpp)
add_executable(smp-convertor ${SOURCE_FILES} smp_convertor.cpp)

project(bounded2smp)
add_executable(bounded2smp bounded2smp.cpp)


project(ie2smp)
add_executable(ie2smp ie2smp.cpp)

