#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH -t 5:00:00
#SBATCH --mem 128GB
#SBATCH -o /scratch/tzp6pz/rna2/ouput_rna_ballgown_slurm.txt
#SBATCH -p standard
#SBATCH -A zanglab

#this `for` loop, will take rna-seq fastq files as input, output filtered BAM
#files, and use stringtie partially so that we can explort to R for analysis.

module load goolf/7.1.0_3.1.4  R/4.0.0

echo "Running code for Ballgown"

echo "Making TXT file for Ballgown"
sh /scratch/tzp6pz/rna2/ballgown_csv_maker.sh > /scratch/tzp6pz/rna2/my_index/ballgown/ballgown.txt

#making directories for results from R
mkdir -p /scratch/tzp6pz/rna2/results/r_ballgown

#cd /scratch/tzp6pz/rna2/

R CMD BATCH /scratch/tzp6pz/rna2/rna_ballgown.r /scratch/tzp6pz/rna2/output_rna_ballgown.txt

echo "SLURM FINISHED"
#
