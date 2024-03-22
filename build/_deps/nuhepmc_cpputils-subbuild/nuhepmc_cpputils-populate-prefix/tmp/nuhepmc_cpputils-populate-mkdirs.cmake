# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/juliatenavidal/Downloads/RadiativeTest/emMCRadCorr/build/_deps/nuhepmc_cpputils-src"
  "/Users/juliatenavidal/Downloads/RadiativeTest/emMCRadCorr/build/_deps/nuhepmc_cpputils-build"
  "/Users/juliatenavidal/Downloads/RadiativeTest/emMCRadCorr/build/_deps/nuhepmc_cpputils-subbuild/nuhepmc_cpputils-populate-prefix"
  "/Users/juliatenavidal/Downloads/RadiativeTest/emMCRadCorr/build/_deps/nuhepmc_cpputils-subbuild/nuhepmc_cpputils-populate-prefix/tmp"
  "/Users/juliatenavidal/Downloads/RadiativeTest/emMCRadCorr/build/_deps/nuhepmc_cpputils-subbuild/nuhepmc_cpputils-populate-prefix/src/nuhepmc_cpputils-populate-stamp"
  "/Users/juliatenavidal/Downloads/RadiativeTest/emMCRadCorr/build/_deps/nuhepmc_cpputils-subbuild/nuhepmc_cpputils-populate-prefix/src"
  "/Users/juliatenavidal/Downloads/RadiativeTest/emMCRadCorr/build/_deps/nuhepmc_cpputils-subbuild/nuhepmc_cpputils-populate-prefix/src/nuhepmc_cpputils-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/juliatenavidal/Downloads/RadiativeTest/emMCRadCorr/build/_deps/nuhepmc_cpputils-subbuild/nuhepmc_cpputils-populate-prefix/src/nuhepmc_cpputils-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/juliatenavidal/Downloads/RadiativeTest/emMCRadCorr/build/_deps/nuhepmc_cpputils-subbuild/nuhepmc_cpputils-populate-prefix/src/nuhepmc_cpputils-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
