#!/bin/bash
#!
#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! Name of the job:
#SBATCH -J Proxies

#! Which project should be charged:
#SBATCH -A MRC-EPID-SL0-CPU

#! How many whole nodes should be allocated?
# In your case you should leave this at 1

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task 4

#! Specify required run time
#SBATCH --time=10:00:00

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! Do not change:
#SBATCH -p epid

#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load default-impi                   # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:
export CHROM=${1} 

/rfs/project/rfs-mpB3sSsgAn4/Programs/GCTA/gcta64 \
--bfile /rfs/project/rfs-mpB3sSsgAn4/Studies/UKBB/Imputed/bgen_files/Revised_Imputation_2018/Random25K/UKBB_v3Imp_Chr${CHROM}_WhiteUnrel_Random25K_NoDups \
--remove /rfs/project/rfs-mpB3sSsgAn4/Studies/UKBB/IncExc_Lists/EXCLUDEFOR_White_Euro_Relateds_v9.ImpRev.samples \
--chr ${CHROM} \
--maf 0.001 \
--ld /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/F_Testosterone_SNPs.txt \
--ld-wind 1000 \
--ld-sig 0.000000001 \
--out /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/M_Estradiol/Output_Chr${CHROM}.proxies

