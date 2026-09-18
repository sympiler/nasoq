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

    # CLAPACK's vendored f2c.h defines `integer` as `long int`, which is
    # 8 bytes on 64-bit macOS/Linux. Every LAPACK routine in SRC/*.c calls
    # out to external BLAS (dsyrk_, dgemv_, dscal_, idamax_, etc. are not
    # part of this glob and resolve to OpenBLAS::OpenBLAS at link time; see
    # target_link_libraries(nasoq PRIVATE clapack::clapack OpenBLAS::OpenBLAS)
    # in the top-level CMakeLists.txt). OpenBLAS (like virtually every
    # LP64 gfortran-built BLAS) uses the standard 4-byte Fortran INTEGER
    # for its interface. Passing 8-byte `long*` dimension/increment
    # arguments (n, lda, incx, ...) across that ABI boundary to a callee
    # expecting 4-byte `int*` silently corrupts the calls and makes every
    # NASOQ_USE_CLAPACK + OpenBLAS factorization return NaN. Patch the
    # fetched header so `integer`/`uinteger` are 4 bytes, matching
    # OpenBLAS's interface, before compiling any of CLAPACK's sources.
    file(READ "${clapack_SOURCE_DIR}/INCLUDE/f2c.h" _clapack_f2c_h)
    string(REPLACE
        "typedef long int integer;\ntypedef unsigned long int uinteger;"
        "typedef int integer;\ntypedef unsigned int uinteger;"
        _clapack_f2c_h_patched "${_clapack_f2c_h}")
    # `logical` (used e.g. as dlapmt_'s FORWRD argument) is the same
    # Fortran INTEGER-width type and must be patched for the same reason.
    string(REPLACE
        "typedef long int logical;"
        "typedef int logical;"
        _clapack_f2c_h_patched "${_clapack_f2c_h_patched}")
    if(_clapack_f2c_h_patched STREQUAL _clapack_f2c_h)
        message(FATAL_ERROR "clapack.cmake: failed to patch f2c.h's `integer` "
            "typedef to 4 bytes (upstream f2c.h layout changed?); NASOQ's "
            "CLAPACK-to-OpenBLAS calls would silently corrupt with a 64-bit "
            "`integer` and cause NaN results.")
    endif()
    file(WRITE "${clapack_SOURCE_DIR}/INCLUDE/f2c.h" "${_clapack_f2c_h_patched}")

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
