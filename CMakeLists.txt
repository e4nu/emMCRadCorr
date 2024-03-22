cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(emMCRadCorr)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#include(cmake/CPM.cmake)

# Add your source files
set(SOURCES
    src/PlottingUtils.cxx
    src/RadiativeCorrUtils.cxx
    src/processor_skeleton.cxx
    src/radiate_flux.cxx
)

# Add executable target
add_executable(radiate_flux ${SOURCES})
add_executable(processor_skeleton ${SOURCES})

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
target_link_libraries(radiate_flux PUBLIC ROOT::Physics)
target_link_libraries(radiate_flux PUBLIC ROOT::Tree)
target_link_libraries(radiate_flux PUBLIC ROOT::Hist)

include(CPM.cmake/cmake/CPM.cmake)
CPMAddPackage(
  NAME NuHepMC_CPPUtils
  GIT_TAG main
  GIT_REPOSITORY "https://github.com/NuHepMC/cpputils.git"
  OPTIONS "NuHepMC_CPPUtils_BUILTIN_HEPMC3 ON"
)
target_link_libraries(processor_skeleton PUBLIC NuHepMC::CPPUtils)

#target_include_directories(radiate_flux PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/utils)

