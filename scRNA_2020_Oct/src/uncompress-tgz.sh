#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --error=uncompress.err
#SBATCH --output=uncompress.out
#SBATCH --job-name uncompress-and-gz-fastq
#SBATCH -p hpcbio

set -euo pipefail

module load pigz/2.3.4-IGB-gcc-4.9.4

pigz -p$SLURM_CPUS_PER_TASK -dc $TGZ | tar -xvf -

for i in *.fastq; do
    pigz -p $SLURM_CPUS_PER_TASK $i
done
