#!/bin/bash
#SBATCH -N 1 
#SBATCH -n 16
#SBATCH --mem 32G
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user drnevich@illinois.edu
#SBATCH -J NBID_cal07-%A
#SBATCH -D /home/groups/hpcbio_shared/cbrooke_lab/single_cell/src
#SBATCH -p hpcbio

# NOTES on options above:
# 1. If you are re-using this script, please change the --mail-user above to your e-mail address!!!
# 2. The -D location was set primarily to control where the slurm-*.out files end up.
#    FULL PATH NAME MUST BE USED HERE!!


#change the directory to the main one containing data/, results/ and src/. All relative file paths in the
#codes below will be to this location

cd /home/groups/hpcbio_shared/cbrooke_lab/single_cell/


#load the program needed; 

module load R/3.6.0-IGB-gcc-8.2.0


Rscript Run_simpleSingleCell_NBID_NCore_VirGenes.R \
        results/test_filter_and_normalize/2020-04-18-cal07_InfectedOnly.rds \
        results/NBID_output/2020-04-18-cal07 $SLURM_CPUS_ON_NODE



echo "done"

