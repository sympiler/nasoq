# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/metis-src"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/metis-build"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/metis-subbuild/metis-populate-prefix"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/metis-subbuild/metis-populate-prefix/tmp"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/metis-subbuild/metis-populate-prefix/src/metis-populate-stamp"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/metis-subbuild/metis-populate-prefix/src"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/metis-subbuild/metis-populate-prefix/src/metis-populate-stamp"
)

set(configSubDirs Debug)
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/metis-subbuild/metis-populate-prefix/src/metis-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/metis-subbuild/metis-populate-prefix/src/metis-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
