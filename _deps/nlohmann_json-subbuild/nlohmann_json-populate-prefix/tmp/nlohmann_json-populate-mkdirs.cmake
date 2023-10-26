# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/nlohmann_json-src"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/nlohmann_json-build"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/nlohmann_json-subbuild/nlohmann_json-populate-prefix"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/nlohmann_json-subbuild/nlohmann_json-populate-prefix/tmp"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/nlohmann_json-subbuild/nlohmann_json-populate-prefix/src/nlohmann_json-populate-stamp"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/nlohmann_json-subbuild/nlohmann_json-populate-prefix/src"
  "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/nlohmann_json-subbuild/nlohmann_json-populate-prefix/src/nlohmann_json-populate-stamp"
)

set(configSubDirs Debug)
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/nlohmann_json-subbuild/nlohmann_json-populate-prefix/src/nlohmann_json-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "E:/McMaster University/NASOQ-Coding-Task/nasoq/_deps/nlohmann_json-subbuild/nlohmann_json-populate-prefix/src/nlohmann_json-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
