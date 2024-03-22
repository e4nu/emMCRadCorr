include("/Users/juliatenavidal/Downloads/RadiativeTest/emMCRadCorr/CPM.cmake/cmake/CPM.cmake")
CPMAddPackage("NAME;NuHepMC_CPPUtils;GIT_TAG;main;GIT_REPOSITORY;https://github.com/NuHepMC/cpputils.git;OPTIONS;NuHepMC_CPPUtils_BUILTIN_HEPMC3 ON")
set(NuHepMC_CPPUtils_FOUND TRUE)