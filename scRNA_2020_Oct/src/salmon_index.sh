#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=36G
#SBATCH --mail-user=drnevich@illinois.edu
#SBATCH --mail-type END,FAIL
#SBATCH -J salmonindex
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/ningwang/2020Sep-RNASeq/src/slurm-out
#SBATCH -p hpcbio

# NOTES on options above:
# 1. If you are re-using this script, please change the --mail-user above to your e-mail address!!!
# 2. The -D location was set primarily to control where the slurm-*.out files end up.
#    FULL PATH NAME MUST BE USED HERE!!

#change the directory to the main one containing data/, results/ and src/. All relative file paths in the
#codes below will be to this location

cd ../../

#load the program needed; name may change in future so check each time with 'module avail Salmon' on command line

module load Salmon/1.2.1-IGB-gcc-8.2.0

#Make the index using the new decoy-aware method.
#See README_ref_modifications.txt for how to make the gentrome.fa.gz

salmon index \
-t data/genome/gentrome.fa.gz \
-i data/genome/salmon_1.2.1_CriGri_1.0_annot104 \
--decoys data/genome/decoys.txt \
-p $SLURM_NTASKS
