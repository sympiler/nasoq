#
# Copyright 2020 Adobe. All rights reserved.
# This file is licensed to you under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may obtain a copy
# of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under
# the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR REPRESENTATIONS
# OF ANY KIND, either express or implied. See the License for the specific language
# governing permissions and limitations under the License.
#
if(TARGET clapack::clapack)
    return()
endif()

message(STATUS "Third-party (external): creating target 'clapack::clapack'")


include(FetchContent)
FetchContent_Declare(
    clapack
    URL "https://www.netlib.org/clapack/clapack-3.2.1-CMAKE.tgz"
)
FetchContent_GetProperties(clapack)
if(NOT clapack_POPULATED)
    FetchContent_Populate(clapack)

    file(GLOB clapack_SOURCE_FILES "${clapack_SOURCE_DIR}/SRC/*.c")
    add_library(clapack ${clapack_SOURCE_FILES})
    target_compile_definitions(clapack PRIVATE "NO_BLAS_WRAP")
    target_include_directories(clapack PUBLIC "${clapack_SOURCE_DIR}/INCLUDE")

    # just add in the f2c files needed for NASOQ's use of CLAPACK
    target_sources(clapack PRIVATE
        "${clapack_SOURCE_DIR}/F2CLIBS/libf2c/i_nint.c"
        "${clapack_SOURCE_DIR}/F2CLIBS/libf2c/s_cmp.c"
        "${clapack_SOURCE_DIR}/F2CLIBS/libf2c/s_copy.c"
    )
endif()

add_library(clapack::clapack ALIAS clapack)
