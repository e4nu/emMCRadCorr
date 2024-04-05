# e4nu Scripts

This directory contains python scripts used to run jobs at the FNAL grid

- submit_rad_GENIE_gpvm.py can be used to run full GENIE + emMCRadCorr simulations

Example:
```
python submit_rad_GENIE_gpvm.py --directory /pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/TestSubmission --xsec /pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/TestSubmission/total_xsec.xml --nevents 1000000
```