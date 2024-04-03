cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(emMCRadCorr)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add your source files
set(SOURCES
    src/Utils.h
    src/RadiativeCorrUtils.cxx
    src/radiate_flux.cxx
)

# Add executable target
add_executable(radiate_flux ${SOURCES})
add_executable(processor_skeleton src/processor_skeleton.cxx)

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
target_link_libraries(radiate_flux PUBLIC ROOT::Physics)
target_link_libraries(radiate_flux PUBLIC ROOT::Tree)
target_link_libraries(radiate_flux PUBLIC ROOT::Hist)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake/Modules)

include(CPM)
CPMAddPackage(
  NAME NuHepMC_CPPUtils
  GIT_TAG main
  GIT_REPOSITORY "https://github.com/NuHepMC/cpputils.git"
  OPTIONS "NuHepMC_CPPUtils_BUILTIN_HEPMC3 ON"
)
target_link_libraries(processor_skeleton PUBLIC NuHepMC::CPPUtils)
target_link_libraries(radiate_flux PUBLIC NuHepMC::CPPUtils)
