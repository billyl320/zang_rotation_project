#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH -t 5:00:00
#SBATCH --mem 128GB
#SBATCH -o /scratch/tzp6pz/chip/slurm_files/ouput_homer.txt
#SBATCH -p standard
#SBATCH -A zanglab

#this `for` loop, will take rna-seq fastq files as input, output filtered BAM
#files, and use stringtie partially so that we can explort to R for analysis.

module load gcc/7.1.0 homer/4.10

echo "Running code for HOMER for "

#making directory for homer results
mkdir -p /scratch/tzp6pz/chip/results/homer

#annotate peaks
#take 8th column and report the percentages of each type
annotatePeaks.pl /scratch/tzp6pz/chip/macs/AR_1_peaks.bed hg38   > /scratch/tzp6pz/chip/results/homer/outputfile.txt

#r script to perform the analysis
module load goolf/7.1.0_3.1.4  R/4.0.0

cd /scratch/tzp6pz/chip/results/homer/

R CMD BATCH r_genomic_distribution_homer.r output_r_genomic_distribution_homer.txt

#motif analysis
#findMotifsGenome.pl /scratch/tzp6pz/chip/macs/AR_1_peaks.bed hg38 /scratch/tzp6pz/chip/results/homer/ -size 200 -mask -preparsedDir /scratch/tzp6pz/chip/results/homer/

echo "SLURM FINISHED"
#
