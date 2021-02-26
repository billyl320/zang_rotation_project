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

module load gcc/7.1.0 homer/4.10

echo "Running code for HOMER for $base"

#making directory for homer results
mkdir -p /scratch/tzp6pz/chip/results/homer
mkdir -p /scratch/tzp6pz/chip/results/homer/${base}

#annotate peaks
#take 8th column and report the percentages of each type
echo "Annotating peaks"
annotatePeaks.pl /scratch/tzp6pz/chip/macs/AR_1_peaks.bed hg38   > /scratch/tzp6pz/chip/results/homer/${base}/outputfile.txt

#motif analysis
echo "Motif Analysis"
findMotifsGenome.pl /scratch/tzp6pz/chip/macs/AR_1_peaks.bed hg38 /scratch/tzp6pz/chip/results/homer/${base} -size 200 -mask -preparsedDir /scratch/tzp6pz/chip/results/homer/${base}


#
