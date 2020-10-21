#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --error=uncompress.err
#SBATCH --output=uncompress.out
#SBATCH --job-name uncompress-and-gz-fastq
#SBATCH --partition=hpcbio

set -euo pipefail

module load pbzip2/1.1.13-IGB-gcc-4.9.4
module load pigz/2.3.4-IGB-gcc-4.9.4

pbzip2 -p$SLURM_CPUS_PER_TASK -dc $TBZ | tar -xvf -

for i in *.fastq; do
    pigz -p $SLURM_CPUS_PER_TASK $i
done
