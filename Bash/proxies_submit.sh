#!/bin/bash
#!
#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! Name of the job:  
#SBATCH -J MTproxies

#! Which project should be charged:
#! MRC private nodes
#SBATCH -p epid
#SBATCH -A MRC-EPID-SL0-CPU

#! How many whole nodes should be allocated?
#! In your case you should leave this at 1
#SBATCH --nodes=1

#! Specify required run time
#SBATCH --time=02:00:00

#SBATCH --ntasks=20
#! one CPU takes one job

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! ############################################################
#! Modify the settings below to specify the application's environment, location
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
./etc/profile.d/modules.sh                
module purge                               

#! Insert additional module load commands after this line if needed:
module load r-4.0.2-gcc-5.4.0-xyx46xb

#! Work directory (i.e. where the job will run):
workdir="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/"

R < /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/proxies.R --no-save > /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/M_T_BMI_proxies.log
