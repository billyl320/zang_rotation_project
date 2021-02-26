#!/bin/bash/

#augmented from https://hbctraining.github.io/Intro-to-ChIPseq/lessons/04_automation.html
# This script takes a fastq file of RNA-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Hisat2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates using stringtie.
# USAGE: sh rna_hisat2.sh <path to the fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1
cpus=$2

# grab base of filename for naming outputs
base=`basename ${fq} .bed`
echo "Sample name is $base"
echo "Number of CPUS: $cpus"

#r script to perform the analysis
module load goolf/7.1.0_3.1.4  R/4.0.0

cd /scratch/tzp6pz/chip/results/homer/${base}

R CMD BATCH /scratch/tzp6pz/chip/results/homer/r_genomic_distribution_homer.r output_r_genomic_distribution_homer.txt

#
