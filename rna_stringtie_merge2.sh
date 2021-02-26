#!/bin/bash/

#augmented from https://hbctraining.github.io/Intro-to-ChIPseq/lessons/04_automation.html
# This script takes a fastq file of RNA-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Hisat2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates using stringtie.
# USAGE: sh rna_hisat2.sh <path to the fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file
cpus=$2

# directory with bowtie genome files
genome=~/data/indexes/genome_tran
genome_file=~/data/genes/genome.gtf
#assembly file
#string_assembly=/scratch/tzp6pz/rna/results/stringtie/samfiles/assembly/
ass2=~/my_index/assembly
mkdir -p $ass2
#merge files location
merge_out=~/my_index/assembly

# set up the software environment
module load gcc/9.2.0
module load stringtie/2.1.0

###
# https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/

mkdir -p ~/my_index
cd ~/my_index/assembly/

echo "Create TXT file for Stringtie Merge"
sh ~/mergelist_maker.sh > ${merge_out}/mergelist.txt

echo "Starting Stringtie Merge"
cd ${merge_out}
stringtie --merge -p $cpus -G $genome_file -o stringtie_merged.gtf ~/data/mergelist.txt
echo "DID THIS MANUALLY - NOT SURE WHY?"
#stringtie --merge -G ~/data/genes/genome.gtf -o stringtie_merged.gtf mergelist.txt

echo "Comparing using gffcompare"
module load gcc/7.1.0
module load gffcompare/0.11.6
gffcompare -r $genome_file -G -o merged ${merge_out}/stringtie_merged.gtf


#
