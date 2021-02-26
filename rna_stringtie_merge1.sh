#!/bin/bash/

#augmented from https://hbctraining.github.io/Intro-to-ChIPseq/lessons/04_automation.html
# This script takes a fastq file of RNA-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Hisat2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates using stringtie.
# USAGE: sh rna_hisat2.sh <path to the fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1
cpus=$2

# grab base of filename for naming outputs
base=`basename ${fq} .bam`
echo "Sample name is $base"
echo "Number of CPUS: $cpus"

# directory with bowtie genome files
genome=~/data/indexes/genome_tran
genome_file=~/data/genes/genome.gtf
# make all of the output directories
# The -p option means mkdir will create the whole path if it
# does not exist and refrain from complaining if it does exist
mkdir -p ~/my_index/assembly
mkdir -p ~/my_index/map
#mkdir -p ~/results/samtools/intermediate_bams

## set up file names
align_sorted=~/my_index/map/${base}_sorted.bam
assembly_out=~/my_index/assembly/${base}.gtf

## set up more variables for 2 additional directoties to help clean up the results folder
#intermediate_bams=/scratch/tzp6pz/rna/results/bowtie2/intermediate_bams

# set up the software environment
module load gcc/9.2.0
module load stringtie/2.1.0

mkdir -p ~/my_index
cd ~/my_index/


echo "Assembling via Stringtie"
stringtie -e -l $base -p $cpus -G $genome_file -o $assembly_out $align_sorted

echo "#####################"
echo "#####################"
echo "#####################"
