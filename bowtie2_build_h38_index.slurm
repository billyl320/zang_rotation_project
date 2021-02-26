#!/bin/bash
#SBATCH --job-name=bowtie2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=05:00:00
#SBATCH -o /scratch/tzp6pz/chip/slurm_files/bowtie2_build_h38_index.txt
#SBATCH --partition=standard
#SBATCH --account=zanglab

#Load the Bowtie Module
module load gcc/9.2.0
module load bowtie2/2.2.9

# Indexing a reference genome
bowtie2-build /nv/vol190/zanglab/zw5j/data/index/hg38.ucsc.fa /scratch/tzp6pz/chip/bowtie_index/
