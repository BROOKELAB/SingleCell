#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --mail-user=jholmes5@illinois.edu
#SBATCH -J fastqc
#SBATCH -a 1-12

module load FastQC/0.11.5-IGB-gcc-4.9.4-Java-1.8.0_121

cd $SLURM_SUBMIT_DIR

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" PI-filenames.txt)

# Single end
fastqc -outdir=../results/fastqc_raw/ --noextract ../data/raw-seq/${line}

# Paired end
# fastqc -outdir=../results/fastqc_raw/ --noextract ../data/raw-seq/${line}R1_001.fastq
# fastqc -outdir=../results/fastqc_raw/ --noextract ../data/raw-seq/${line}R2_001.fastq
