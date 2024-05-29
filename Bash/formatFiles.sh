#!/bin/bash

#!
#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! Name of the job:
#SBATCH -J GCTA_Cojo

#! Which project should be charged:
#SBATCH -A MRC-EPID-SL0-CPU

#! How many whole nodes should be allocated?
# In your case you should leave this at 1
#SBATCH --nodes=1

#! Specify required run time
#SBATCH --time=00:10:00

#SBATCH --exclusive

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
cd /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/F_Testosterone/

perl processProxyOutput.pl Output_Chr1.proxies.rsq.ld Output_Chr1.proxies.snp.ld Output_Chr1.proxies.r.ld Chr1_Formatted_r207.txt 0.6 1
perl processProxyOutput.pl Output_Chr2.proxies.rsq.ld Output_Chr2.proxies.snp.ld Output_Chr2.proxies.r.ld Chr2_Formatted_r207.txt 0.6 2
perl processProxyOutput.pl Output_Chr3.proxies.rsq.ld Output_Chr3.proxies.snp.ld Output_Chr3.proxies.r.ld Chr3_Formatted_r207.txt 0.6 3
perl processProxyOutput.pl Output_Chr4.proxies.rsq.ld Output_Chr4.proxies.snp.ld Output_Chr4.proxies.r.ld Chr4_Formatted_r207.txt 0.6 4
perl processProxyOutput.pl Output_Chr5.proxies.rsq.ld Output_Chr5.proxies.snp.ld Output_Chr5.proxies.r.ld Chr5_Formatted_r207.txt 0.6 5
perl processProxyOutput.pl Output_Chr6.proxies.rsq.ld Output_Chr6.proxies.snp.ld Output_Chr6.proxies.r.ld Chr6_Formatted_r207.txt 0.6 6
perl processProxyOutput.pl Output_Chr7.proxies.rsq.ld Output_Chr7.proxies.snp.ld Output_Chr7.proxies.r.ld Chr7_Formatted_r207.txt 0.6 7
perl processProxyOutput.pl Output_Chr8.proxies.rsq.ld Output_Chr8.proxies.snp.ld Output_Chr8.proxies.r.ld Chr8_Formatted_r207.txt 0.6 8
perl processProxyOutput.pl Output_Chr9.proxies.rsq.ld Output_Chr9.proxies.snp.ld Output_Chr9.proxies.r.ld Chr9_Formatted_r207.txt 0.6 9
perl processProxyOutput.pl Output_Chr10.proxies.rsq.ld Output_Chr10.proxies.snp.ld Output_Chr10.proxies.r.ld Chr10_Formatted_r207.txt 0.6 10
perl processProxyOutput.pl Output_Chr11.proxies.rsq.ld Output_Chr11.proxies.snp.ld Output_Chr11.proxies.r.ld Chr11_Formatted_r207.txt 0.6 11
perl processProxyOutput.pl Output_Chr12.proxies.rsq.ld Output_Chr12.proxies.snp.ld Output_Chr12.proxies.r.ld Chr12_Formatted_r207.txt 0.6 12
perl processProxyOutput.pl Output_Chr13.proxies.rsq.ld Output_Chr13.proxies.snp.ld Output_Chr13.proxies.r.ld Chr13_Formatted_r207.txt 0.6 13
perl processProxyOutput.pl Output_Chr14.proxies.rsq.ld Output_Chr14.proxies.snp.ld Output_Chr14.proxies.r.ld Chr14_Formatted_r207.txt 0.6 14
perl processProxyOutput.pl Output_Chr15.proxies.rsq.ld Output_Chr15.proxies.snp.ld Output_Chr15.proxies.r.ld Chr15_Formatted_r207.txt 0.6 15
perl processProxyOutput.pl Output_Chr16.proxies.rsq.ld Output_Chr16.proxies.snp.ld Output_Chr16.proxies.r.ld Chr16_Formatted_r207.txt 0.6 16
perl processProxyOutput.pl Output_Chr17.proxies.rsq.ld Output_Chr17.proxies.snp.ld Output_Chr17.proxies.r.ld Chr17_Formatted_r207.txt 0.6 17
perl processProxyOutput.pl Output_Chr18.proxies.rsq.ld Output_Chr18.proxies.snp.ld Output_Chr18.proxies.r.ld Chr18_Formatted_r207.txt 0.6 18
perl processProxyOutput.pl Output_Chr19.proxies.rsq.ld Output_Chr19.proxies.snp.ld Output_Chr19.proxies.r.ld Chr19_Formatted_r207.txt 0.6 19
perl processProxyOutput.pl Output_Chr20.proxies.rsq.ld Output_Chr20.proxies.snp.ld Output_Chr20.proxies.r.ld Chr20_Formatted_r207.txt 0.6 20
perl processProxyOutput.pl Output_Chr21.proxies.rsq.ld Output_Chr21.proxies.snp.ld Output_Chr21.proxies.r.ld Chr21_Formatted_r207.txt 0.6 21
perl processProxyOutput.pl Output_Chr22.proxies.rsq.ld Output_Chr22.proxies.snp.ld Output_Chr22.proxies.r.ld Chr22_Formatted_r207.txt 0.6 22

cat *_Formatted_r207.txt > M_Estradiol_r06.txt
