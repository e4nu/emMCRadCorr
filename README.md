# emMCRadCorr
This is a package that can be used to add radiative corrections to any MC event generator. It is developed oriented to neutrino event generators running on electron-scattering mode. The package modifies the output of neutrino event generators in NuHEPMC format. It can also be used to modify GENIE gst files.

Radiative corrections are included in a number of steps: 
1. Pre-compute incoming electron radiation
2.  Use radiated incoming electron energy spectra as input flux for your preferred event generator
3.  Generate electron-scattering events on the target of interest
4.  Modify the generator output to account for incoming and outgoing electron radiation
5.  Use the event kinematics to compute the corresponding radiative cross-section weights

<img width="726" alt="image" src="https://github.com/e4nu/emMCRadCorr/assets/36236227/8c2acb2e-68c2-465c-8613-b743efd6d3f5">

