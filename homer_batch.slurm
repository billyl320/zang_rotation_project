#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH -t 5:00:00
#SBATCH --mem 128GB
#SBATCH -o /scratch/tzp6pz/chip/slurm_files/ouput_homer_batch.txt
#SBATCH -p standard
#SBATCH -A zanglab

for fq in /scratch/tzp6pz/chip/macs/AR_*_peaks.bed
do
  echo "running Homer analyses on: $fq"
  sh homer.sh $fq $SLURM_CPUS_PER_TASK
done


for fq in /scratch/tzp6pz/chip/macs/AR_*_peaks.bed
do
  echo "Performing R analysis on: $fq"
  sh r_genomic_dist_batch.sh $fq $SLURM_CPUS_PER_TASK
done


echo "SLURM FINISHED"
#
