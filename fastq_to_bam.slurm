#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH -t 5:00:00
#SBATCH --mem 128GB
#SBATCH -o /scratch/tzp6pz/chip/slurm_files/fastq_to_bam.txt
#SBATCH -p standard
#SBATCH -A zanglab

# this `for` loop, will take chip-seq fastq files as input and output filtered BAM files ready for peak calling.

for fq in /scratch/tzp6pz/chip/rawdata_LNCaP_AR/*.fastq
do
  echo "running analysis on $fq"
  sh chipseq_analysis_on_input_file.sh $fq $SLURM_CPUS_PER_TASK
done
