# emMCRadCorr
This is a package that can be used to add radiative corrections to any MC event generator. It is developed oriented to neutrino event generators running on electron-scattering mode. The package modifies the output of neutrino event generators in NuHEPMC format. It can also be used to modify GENIE gst files.

Radiative corrections are included in a number of steps: 
1. Given an incident electron flux (typically mono-chromatic), compute the incoming electron radiation and the resulting modified incident electron flux,
2. Generate electron-scattering events on the target and processes of interest using your preferred event generator,
3. Store the events in NuHEPMC data format,
4. Use the event kinematics to calculate the appropriate weight to account for vertex and vacuum-polarization corrections,
5. Generate the real photon (if any) radiated by the scattered electron,
6. Modify the NuHEPMC event record to include the incoming and outgoing electron radiation and the weight.
   
## Build software
In the FNAL farm, 
```
source emMCRadCorr_gpvm_env.sh;
mkdir build ;
cd build ;
cmake .. ;
make ;
```
This gives you access to the main apps.

## Pre-compute electron flux for event generators
The first step to account for radiative effects is to compute the energy spectra of the real photons emitted by the incoming electron due to internal and external bremsstrahlung emission. 
Internal bremsstrahlung describes photon emission due to the Coulomb field of the target nuclei, whilst external bremsstrahlung describes photon emission in the field of nuclei other than the one participating in the scattering. We only account for emission due to electron radiation, as radiation due to other particles is negligible.

For instance, one can radiate a 4.32 GeV electron beam on a H CLAS6 target using the following command:
```./radiate_flux --output-file radflux_H_4325MeV_simc.root --target 1000010010 --Emin 3.4 --Emax 4.35 --ebeam 4.325 --rad-model "simc" --resolution 0.001```
Where the model "simc" includes external radiation via the method explained in PhysRevC.64.054610. The expected output is:
<p align="center">
  <img width="531" alt="image" src="https://github.com/e4nu/emMCRadCorr/assets/36236227/3f2484e1-bcde-40a6-88c0-3393c7938e4e">
</p>
It is stored in the radflux_H_4325MeV.root as hradflux. For bookkeeping, the standard fluxes with high precision are computed and stored here: ```/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation```:
- radflux_C_1161MeV.root
- radflux_C_2261MeV.root
- radflux_C_4461MeV.root

## Generate electron-scattering events with GENIE in HepMC3 format
Generate GENIE events with "radiated"-flux. In the example above, we use a Q2min of 0.8 GeV2.
```
python eAScatteringGridSubmitter.py --jobs-topdir <Your_directory_scratch> --total-xsec <path_to_xsec> --config-dir /pnfs/genie/persistent/users/jtenavid/e4nu_files/e4nu-GENIE-config --probe-list 11 --e-ntotevents <N_events_to_run> --flux <path>/radflux_H_4325MeV.root,hradflux --e-minenergy-fluxrange 3.4 --e-maxenergy-fluxrange 4.35 --starting-point 4 --no-ghep-output --gst-output --tune <TUNE> --e-tgt-list 1000010010 --store-comitinfo --subjob-disk 3GB --subjob-memory 3GB --mainjob-memory 10GB --mainjob-disk 10GB --git-location https://github.com/e4nu/Generator_Rad_Q2_08.git --submit-jobs
```
to generate the equivalent plot without radiation effects simply run:
```
python eAScatteringGridSubmitter.py --jobs-topdir <Your_directory_scratch> --total-xsec <path_to_xsec> --config-dir /pnfs/genie/persistent/users/jtenavid/e4nu_files/e4nu-GENIE-config --probe-list 11 --e-ntotevents <N_events_to_run> --ebeam-energy 4.325 --starting-point 4 --no-ghep-output --gst-output --tune <TUNE> --e-tgt-list 1000010010 --store-comitinfo --subjob-disk 3GB --subjob-memory 3GB --mainjob-memory 10GB --mainjob-disk 10GB --git-location https://github.com/e4nu/Generator_Rad_Q2_08.git --submit-jobs
```

## Apply full radiative corrections to your generated events
Up till now, we computed the interaction kinematics at the interaction level, effectively ignoring radiative effects. In order to account for full event kinematics, we must re-store the true monochromatic beam energy, account for bremstrahalung radiation due to the outgoing electron, and weight the cross-section accordingly to account for vertex, vacumm and bremstrahalung corrections. This is accomplished with the ```process_radcorr.cxx``` app. The app loops over the event record and accounts for the effects described above. 

Example:
```
./process_radcorr --input-hepmc3-file <name_file.hepmc3> --output-file <final_name.hepmc3> --target <target_pdg> --rad-model "simc" --thickness <thickness> --max-egamma <maxEgamma> 
```

## Validation against H(e,e'p) data
A number of scripts are available in the plotting folder to compute H(e,e'p) radiative corrections and compare it to data from JLab. In order to reproduce the validation plot for rad corrections:
- Run "root script_simulation.C"
- Plot with "root ProduceExtDataFig2.cpp"
Generations must have a Q2 minimum cut of 1.4 GeV^2.

## Automated scripts
The steps above require the user to run a number of steps. These have been implemented in a script system optimized for the FNAL farm. The script launches GENIE jobs in parallel using the correct "radiated" electron flux for a given target thickness. The final generated events contain fully radiated information. An example of the command used to generate the events for the H(e,e'p) validation is:
```
python submit_rad_GENIE_gpvm.py --directory /pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/Test_newHepMC3_14Q2/ --model vanderhaghen --ebeam-energy 4.325 --xsec /pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/Test_newHepMC3_14Q2/total_xsec.xml --nevents 100000 --event-gen-list EMQE --tune G18_10a_00_000 --genie-git-branch hepmc3-Q214 --genie-git-location https://github.com/e4nu/Generator-NuHepMC.git --input-radflux /pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/radflux_H_4325MeV_simc.root
```











