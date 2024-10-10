#!/bin/bash

# Set up the UPS products needed to build and use GENIE
source /cvmfs/fermilab.opensciencegrid.org/products/genie/bootstrap_genie_ups.sh 
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
setup root v6_26_06b -q e26:p3913:prof
setup cmake v3_23_1 -f Linux64bit+3.10-2.17 
setup ifdhc v2_6_6
export IFDH_CP_MAXRETRIES=0

setup pdfsets v5_9_1b 
setup gdb v8_1 

echo "Setting emMCRadCorr environment variables..."

# Finds the directory where this script is located. This method isn't
# foolproof. See https://stackoverflow.com/a/246128/4081973 if you need
# something more robust for edge cases (e.g., you're calling the script using
# symlinks).
export EMMCRADCORR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
