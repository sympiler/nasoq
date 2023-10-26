# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/catch2-src"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/catch2-build"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/catch2-subbuild/catch2-populate-prefix"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/catch2-subbuild/catch2-populate-prefix/tmp"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/catch2-subbuild/catch2-populate-prefix/src/catch2-populate-stamp"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/catch2-subbuild/catch2-populate-prefix/src"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/catch2-subbuild/catch2-populate-prefix/src/catch2-populate-stamp"
)

set(configSubDirs Debug)
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/catch2-subbuild/catch2-populate-prefix/src/catch2-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/catch2-subbuild/catch2-populate-prefix/src/catch2-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
