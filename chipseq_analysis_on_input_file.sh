#!/bin/bash/

#augmented from https://hbctraining.github.io/Intro-to-ChIPseq/lessons/04_automation.html
# This script takes a fastq file of ChIP-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Bowtie2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates and has multi-mappers and duplicate reads removed using sambamba.
# USAGE: sh chipseq_analysis_on_input_file.sh <path to the fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1
cpus=$2

# grab base of filename for naming outputs
base=`basename $fq .fastq`
echo "Sample name is $base"
echo "Number of CPUS: $cpus"

# directory with bowtie genome index
genome=/scratch/tzp6pz/chip/bowtie_index/genome

# make all of the output directories
# The -p option means mkdir will create the whole path if it
# does not exist and refrain from complaining if it does exist
mkdir -p /scratch/tzp6pz/chip/results/fastqc
mkdir -p /scratch/tzp6pz/chip/results/bowtie2/intermediate_bams

# set up output filenames and locations
fastqc_out=/scratch/tzp6pz/chip/results/fastqc/

## set up file names
align_out=/scratch/tzp6pz/chip/results/bowtie2/${base}_unsorted.sam
align_bam=/scratch/tzp6pz/chip/results/bowtie2/${base}_unsorted.bam
align_sorted=/scratch/tzp6pz/chip/results/bowtie2/${base}_sorted.bam
align_filtered=/scratch/tzp6pz/chip/results/bowtie2/${base}_aln.bam

## set up more variables for 2 additional directoties to help clean up the results folder
bowtie_results=/scratch/tzp6pz/chip/results/bowtie2
intermediate_bams=/scratch/tzp6pz/chip/results/bowtie2/intermediate_bams

# set up the software environment
module load fastqc/0.11.5
module load gcc/9.2.0
module load bowtie2/2.2.9
module load samtools/1.9
module load sambamba

echo "Processing file $fq"

# Run FastQC and move output to the appropriate folder
fastqc $fq -f fastq
mv /scratch/tzp6pz/chip/rawdata_LNCaP_AR/${base}_fastqc.html /scratch/tzp6pz/chip/results/fastqc/

#skipping since did using bowtie2_build_h38_index.slurm
# Indexing a reference genome
#bowtie2-build /nv/vol190/zanglab/zw5j/data/index/hg38.ucsc.fa /scratch/tzp6pz/chip/bowtie_index

# Run bowtie2
echo "Starting Bowtie"
#bowtie2 -p $cpus -q --local -x $genome -U $fq -S $align_out
bowtie2 -p $cpus -x $genome -U $fq -S $align_out

# Create BAM from SAM
echo "Converting SAM to BAM"
samtools view -h -S -b -@ 6 -o $align_bam $align_out

# Sort BAM file by genomic coordinates
echo "Sorting BAM file"
sambamba sort -t 6 -o $align_sorted $align_bam

# Filter out duplicates
echo "Filtering out duplicates"
sambamba view -h -t 6 -f bam -F "[XS] == null and not unmapped " $align_sorted > $align_filtered

# Create indices for all the bam files for visualization and QC
echo "Create indices"
samtools index $align_filtered

# Move intermediate files out of the bowtie2 directory
echo "Cleanup"
mv $bowtie_results/${base}*sorted* $intermediate_bams
