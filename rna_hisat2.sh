#!/bin/bash/

#augmented from https://hbctraining.github.io/Intro-to-ChIPseq/lessons/04_automation.html
# This script takes a fastq file of RNA-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Hisat2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates using stringtie.
# USAGE: sh rna_hisat2.sh <path to the fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1
cpus=$2

# grab base of filename for naming outputs
base=`basename ${fq} _1.fastq`
echo "Sample name is $base"
echo "Number of CPUS: $cpus"

# directory with bowtie genome files
genome=~/data/indexes/genome_tran

# make all of the output directories
# The -p option means mkdir will create the whole path if it
# does not exist and refrain from complaining if it does exist
mkdir -p ~/results/fastqc
mkdir -p ~/my_index/assembly
mkdir -p ~/my_index/map

# set up output filenames and locations
fastqc_out=~/results/fastqc/

## set up file names
align_out=~/my_index/map/${base}.sam
align_bam=~/my_index/map/${base}.bam
align_sorted=~/my_index/map/${base}_sorted.bam

# set up the software environment
module load fastqc/0.11.5
module load gcc/9.2.0
module load hisat2/2.1.0
module load stringtie/2.1.0
module load samtools/1.10

echo "Processing file $fq via fastqc"

# Run FastQC and move output to the appropriate folder
fastqc ~/data/samples/${base}_1.fastq -f fastq
mv ~/data/samples/${base}_1_fastqc.html ~/results/fastqc/

fastqc ~/data/samples/${base}_2.fastq -f fastq
mv ~/data/samples/${base}_2_fastqc.html ~/results/fastqc/


#skipping building hisat2 index since did that in hisat_setup.slurm
mkdir -p ~/my_index
cd ~/my_index

# map each sample using 8 threads
echo "Mapping via hisat"
hisat2 -p $cpus --dta -x $genome -1 ~/data/samples/${base}_1.fastq -2 ~/data/samples/${base}_2.fastq -S $align_out

# Create BAM from SAM
echo "Converting SAM to BAM"
samtools view -h -S -b -@ $cpus -o $align_bam $align_out
# Sort BAM file by genomic coordinates
echo "Sorting BAM file"
samtools sort -@ $cpus -o $align_sorted $align_bam

echo "   "
echo "   "
echo "   "
