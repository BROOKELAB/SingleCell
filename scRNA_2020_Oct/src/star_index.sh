#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH --mail-user=drnevich@illinois.edu
#SBATCH --mail-type END,FAIL
#SBATCH -J starindex
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/ningwang/2020Sep-RNASeq/src/slurm-out
#SBATCH -p hpcbio

# NOTES on options above:
# 1. If you are re-using this script, please change the --mail-user above to your e-mail address!!!
# 2. The -D location was set primarily to control where the slurm-*.out files end up.
#    FULL PATH NAME MUST BE USED HERE!!

#change the directory to the main one containing data/, results/ and src/. All relative file paths in the
#codes below will be to this location

cd ../../

#load the program needed; name may change in future so check each time with 'module avail STAR' on command line

module load STAR/2.7.4a-IGB-gcc-8.2.0

# Create output directory for the STAR index and a scratch directory for STAR to write temp files
# NOTE: the -p option will not only create parental directories if needed, but also will not
# throw an error or overwrite the directory if it already exists

GenomeDir=data/genome/STAR_2.7.4a_CriGri_1.0
mkdir -p $GenomeDir

STAR --runThreadN $SLURM_NTASKS \
     --runMode genomeGenerate \
     --genomeDir $GenomeDir \
     --genomeFastaFiles data/genome/GCF_000223135.1_CriGri_1.0_genomic.fna \
     --limitGenomeGenerateRAM 95000000000 \
     --outTmpDir /scratch/$SLURM_JOB_ID
