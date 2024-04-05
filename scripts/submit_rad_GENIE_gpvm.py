#! /usr/bin/env python

"""\
This script is responsible to run all radiative corrections for a given 
Target and energy

Author: 
      Julia Tena Vidal <jtenavidal \st tauex.tau.ac.il>
      Tel Aviv University
"""
import os,sys,optparse, glob, tarfile, re
GENIE = os.environ.get('GENIE')
if GENIE:
    sys.path.append(GENIE+"/src/scripts/production/python/xsec_splines")
    sys.path.append(GENIE+"/src/scripts/production/python/event_generation/")
import eFluxScatteringGenCommands as eAFlux

op = optparse.OptionParser(usage=__doc__)
op.add_option("--git-location", dest="GIT_LOCATION", default="https://github.com/e4nu/emMCRadCorr.git", help="emMCRadCorr code location in github. Defaulted to %default")
op.add_option("--git-branch", dest="BRANCH", default="master", help="Branch name. Default: %default")
op.add_option("--gpvm-group", dest="GROUP", default="genie", help="Group in the gpvm assigned to your user")
op.add_option("--directory", dest="JOBSTD", default=os.getenv('PWD'), help="Output directory (default: %default)")
op.add_option("--ebeam-energy", dest="EnergyBeam", default=2, help="Beam energy, Default: %default")
op.add_option("--model", dest="MODEL", default="simc", help="Rad corr model to use. Default %default")
op.add_option("--target", dest="TARGET", default=1000010010, help="Target used for calculation. Default %default")
op.add_option("--thickness", dest="THICKNESS", default="0", help="Thickness. Specify for your experiment")
op.add_option("--MaxEGamma", dest="MaxEGamma", default=0.2,help="Maximum energy for the emited photons. Default %default*EnergyBeam")
op.add_option("--eResolution", dest="ERES", default=0.0001,help="Experimental electron energy resolution. Default %default")
op.add_option("--output-radflux", dest="OUTFLUX", default="rad_flux.root",help="Name of output ROOT file containing decayed electron flux. Default %default")
op.add_option("--input-radflux", dest="INFLUX", default="",help="Name of INPUT ROOT file containing decayed electron flux. OPTIONAL argument")
## GENIE SPECIFIC OPTIONS
op.add_option("--genie-topdir", dest="GENIE", default=os.getenv('GENIE'), help = "GENIE topdir: %default")
op.add_option("--genie-version", dest="VERSION", default="master", help="Genie version. Default: %default")
op.add_option("--genie-git-location", dest="GENIE_GIT_LOCATION", default="https://github.com/GENIE-MC/Generator", help="Github location from where to get the GENIE Generator code. Defaulted to %default")
op.add_option("--genie-git-branch", dest="BRANCH", default="master", help="Genie version branch name. Default: %default")
op.add_option("--cycle", dest="CYCLE", default="01", help="Cycle (default: %default)")
op.add_option("--arch", dest="ARCH", default='SL6.x86_64', help="arch number, default: %default")
op.add_option("--production", dest="PROD", default="routine_validation", help="Production (default: %default)")
op.add_option("--config-dir", dest="CONF", default='', help="Path to GENIE config dir")
op.add_option("--tune", dest="TUNE",default="G18_10a_00_000",help="Genie model configuration tag. Default %default")
op.add_option("--xsec", dest="XSEC",default="total_xsec.xml",help="Specify name and location of input xsec file for GENIE event generation. Default %default")
op.add_option("--nevents",dest="NEVNT",default=100000,type="int", help="Number of events to generate, default %default")
op.add_option("--nmaxevents",dest="NMax",default=400000,type="int", help="Max number of events per generation, default %NMax")
op.add_option("--event-gen-list",dest="EvGenList",default="EM",help="Event generator list: EM, EMQE, EMMEC, EMRES, EMDIS. Default %default")
op.add_option("--seed",dest="Seed", default=210921029, help="Set stargint point seed. Default %default")
op.add_option("--gen-runid", dest="RunID", default=0, help="Set Starting run id. Default %default")
op.add_option("--gst-output", dest="GSTOutput", default=False, action="store_true",help="Store gst root file.")
op.add_option("--no-ghep-output", dest="NoGHEPOutput", default=False, action="store_true",help="GHEP GENIE files is removed to reduce memory.")

opts, args = op.parse_args()

# Check jobstopdir exists
if not os.path.exists(opts.JOBSTD) :
    os.mkdir(opts.JOBSTD)

#JobSub is made available through the UPS package jobsub_client
os.system("source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup" ) 

if opts.BRANCH: 
    print( ' Cloning emMCRadCorr ' + opts.BRANCH ) 

# Configure grid
EMRADCORRCODE=os.getenv('EMMCRADCORR')
emMCRadCorr_setup_file = EMRADCORRCODE+'/emMCRadCorr_gpvm_env.sh'
emMCRadCorr_pnfs_setup = opts.JOBSTD+"/emMCRadCorr_gpvm_env.sh"
if os.path.exists(emMCRadCorr_pnfs_setup) : 
    os.remove(emMCRadCorr_pnfs_setup)
os.system('cp '+emMCRadCorr_setup_file+' '+opts.JOBSTD )

counter = 0 
name_out_file = "rad_corr"

# xml script
if os.path.exists(opts.JOBSTD+"/grid_submission.xml") : 
    os.remove(opts.JOBSTD+"/grid_submission.xml")
grid = open( opts.JOBSTD+"/grid_submission.xml", 'w' ) 

# 1- write job responsible to run radiated flux
if opts.INFLUX=="" :
    if os.path.exists(opts.JOBSTD+"/rad_flux.sh") :
        os.remove(opts.JOBSTD+"/rad_flux.sh")

    script = open( opts.JOBSTD+"/rad_flux.sh", 'w' ) 
    script.write("#!/bin/bash \n")
    script.write("source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups \n")
    script.write("setup ifdhc v2_6_6 \n")
    script.write("export IFDH_CP_MAXRETRIES=0 \n")
    script.write("setup pdfsets v5_9_1b \n")
    script.write("setup gdb v8_1 \n")
    script.write("cd $CONDOR_DIR_INPUT \n")
    script.write("git clone "+opts.GIT_LOCATION+" ;\n")
    script.write("cd emMCRadCorr ; source emMCRadCorr_gpvm_env.sh ; make ;\n")
    #write main command
    script.write("./radiate_flux --output-file "+opts.OUTFLUX+" --target "+str(opts.TARGET)+" --Emin "+str(opts.EnergyBeam-opts.MaxEGamma*opts.EnergyBeam)+" --Emax "+str(opts.EnergyBeam+0.02)+" --ebeam "+str(opts.EnergyBeam)+" --rad-model "+opts.MODEL+" --resolution "+str(opts.ERES)+" \n")
    script.write("ifdh cp -D $CONDOR_DIR_INPUT/"+opts.OUTFLUX+" "+opts.JOBSTD+" \n")
    grid.write("<serial>\n")
    grid.write("jobsub_submit  -n --memory=1GB --disk=1GB --expected-lifetime=1h  --OS=SL7 --mail_on_error file://"+opts.JOBSTD+"/rad_flux.sh \n")
    grid.write("<serial>\n")

# 2 - Run GENIE jobs on grid
if not os.path.exists(opts.XSEC):
    print("Cross section file does not exist. Abort.")
    exit()

message_thresholds = ""
if os.path.exists(opts.CONF+"/Messenger.xml"):
        message_thresholds = "$CONDOR_DIR_INPUT/conf/Messenger.xml"

# Check whether INCL/G4 have to be configured:
configure_INCL = False
configure_G4 = False 
if "c" in opts.TUNE:
    configure_INCL = True
elif "d" in opts.TUNE:
    configure_G4 = True 

command_dict = {}
grid_setup = os.getenv('GENIE')+'src/scripts/production/python/setup_FNAL.sh' 
genie_setup= os.getenv('GENIE')+'src/scripts/production/python/setup_GENIE.sh'
command_dict.update( eAFlux.eFluxScatteringGenCommands("11",str(opts.TARGET),opts.JOBSTD+opts.OUTFLUX+",hradflux",
                                                       str(opts.EnergyBeam-opts.MaxEGamma*opts.EnergyBeam),
                                                       str(opts.EnergyBeam+0.02),opts.XSEC,opts.NEVNT,opts.TUNE, opts.EvGenList, opts.NMax, 
                                                       opts.Seed, opts.RunID, opts.GSTOutput, opts.NoGHEPOutput,opts.VERSION,
                                                       opts.CONF, opts.ARCH, opts.PROD, opts.CYCLE,"FNAL", opts.GROUP,os.getenv('GENIE_MASTER_DIR'),
                                                       opts.GENIE, opts.JOBSTD,grid_setup,genie_setup,message_thresholds,"4","4","4",opts.BRANCH,
                                                       opts.GENIE_GIT_LOCATION,configure_INCL,configure_G4,True))
command_list = command_dict[4]
command_list_next = command_list
in_serial = False 
if len(command_list) == 1 : # serial
    if in_serial == False: 
        grid.write("<serial>\n")
        in_serial = True
        
        grid.write(command_list[0]+"\n")

        if ( in_serial == True ) : 
            grid.write("</serial>\n")
            in_serial = False
    else : 
        grid.write("<parallel>\n")
        for i in range(len(command_list)) : 
            grid.write(command_list[i]+"\n")
            grid.write("</parallel>\n")

# 3 - Process files
#e_on_1000010010_0.hepmc3
true_nsubruns = opts.NEVNT*1.0/opts.NMax
nsubruns = int(round(opts.NEVNT*1.0/opts.NMax))
if( nsubruns < true_nsubruns ) : nsubruns += 1 
if opts.NEVNT <= opts.NMax : nsubruns = 1
number_files = nsubruns

rad_dir = opts.JOBSTD+"/radcorr/"
if not os.path.exists(rad_dir) : 
    os.mkdir(rad_dir)
name_out_file = "rad_corr"

gst_file_names = []
for i in range(0,number_files):
    gst_file_names.append("e_on_"+str(opts.TARGET)+"_"+str(i)+".hepmc3")

name_out_file = "rad_corr"

if number_files == 1 :
    grid.write("<serial>\n")
else :
    grid.write("<parallel>\n")

for x in range(0,len(gst_file_names)):
    
    if os.path.exists(rad_dir+name_out_file+"_e_on_"+str(opts.TARGET)+"_"+str(x)+".sh"):
        os.remove(rad_dir+name_out_file+"_e_on_"+str(opts.TARGET)+"_"+str(x)+".sh")
    script = open( rad_dir+name_out_file+"_e_on_"+str(opts.TARGET)+"_"+str(x)+".sh", 'w' ) 

    script.write("#!/bin/bash \n")
    script.write("source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups \n")
    script.write("setup ifdhc v2_6_6 \n")
    script.write("export IFDH_CP_MAXRETRIES=0 \n")
    script.write("setup pdfsets v5_9_1b \n")
    script.write("setup gdb v8_1 \n")
    script.write("cd $CONDOR_DIR_INPUT \n")
    script.write("ifdh cp -D "+opts.JOBSTD+"/master-routine_validation_01-eScattering/"+gst_file_names[x]+" $CONDOR_DIR_INPUT/ ;\n \n")
    script.write("git clone "+opts.GIT_LOCATION+" ;\n")
    script.write("cd emMCRadCorr ; source emMCRadCorr_gpvm_env.sh ; make ;\n")
    #write main command
    script.write("./process_radcorr --input-hepmc3-file $CONDOR_DIR_INPUT/"+gst_file_names[x]+" --output-file $CONDOR_DIR_INPUT/rad_corr_"+"e_on_"+str(opts.TARGET)+"_"+str(x)+" --true-EBeam "+str(opts.EnergyBeam)+" --rad-model "+opts.MODEL+" --thickness "+str(opts.THICKNESS)+" --max-egamma "+str(opts.MaxEGamma)+"; \n\n")
    script.write("ifdh cp -D $CONDOR_DIR_INPUT/rad_corr_e_on_"+str(opts.TARGET)+"_"+str(x)+".gst.root "+rad_dir+" \n")

    grid.write("jobsub_submit  -n --memory=4GB --disk=4GB --expected-lifetime=4h  --OS=SL7 --mail_on_error file://"+rad_dir+name_out_file+"_"+str(x)+".sh \n")

    counter += 1

if number_files == 1 :
    grid.write("</serial>\n")
else :
    grid.write("</parallel>\n")

if os.path.exists(opts.JOBSTD+"/fnal_dag_submit.fnal"):
    os.remove(opts.JOBSTD+"/fnal_dag_submit.fnal")

fnal_script = open( opts.JOBSTD+"/fnal_dag_submit.fnal", 'w' ) 
fnal_script.write("#!/bin/bash \n")
fnal_script.write("source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups ;\n")
fnal_script.write("setup fife_utils ;\n")
fnal_script.write("jobsub_submit -G "+opts.GROUP+" --OS=SL7 --memory=10GB --disk=10GB --expected-lifetime=25h -N 1 --role=Analysis --dag file://"+opts.JOBSTD+"/grid_submission.xml;\n")
