cmake_minimum_required(VERSION 0.10)
project(emMCRadCorr)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#include(cmake/CPM.cmake)

# Add your source files
set(SOURCES
    src/utils/RadiativeCorrUtils.cxx
    src/config/ConstantsI.cxx
    src/config/ParticleI.cxx
    src/config/RadConstants.cxx
    src/apps/radiate_flux.cxx
)

# Add executable target
add_executable(radiate_flux ${SOURCES})

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

# Include directories
#target_include_directories(radiate_flux PUBLIC src/apps)

#CPMAddPackage(
#  NAME NuHepMC_CPPUtils
#  GIT_TAG main
#  GIT_REPOSITORY "https://github.com/NuHepMC/cpputils.git"
#  OPTIONS "NuHepMC_CPPUtils_BUILTIN_HEPMC3 ON"
#)

# Include directories if needed
#target_link_libraries(<mytarget> PUBLIC NuHepMC::CPPUtils)

#target_include_directories(MyExecutable PUBLIC include)

# Add libraries if needed
# target_link_libraries(MyExecutable <library_name>)

# Optionally, you can add compiler flags
# target_compile_options(MyExecutable PRIVATE -Wall -Wextra -pedantic)

