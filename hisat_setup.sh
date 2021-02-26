#!/bin/bash/

#augmented from https://hbctraining.github.io/Intro-to-ChIPseq/lessons/04_automation.html
# This script takes a fastq file of RNA-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Hisat2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates using stringtie.
# USAGE: sh rna_hisat2.sh <path to the fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file
cpus=$2

module load gcc/7.1.0
module load python/3.6.8

mkdir -p ~/data/samples/

for fq in 1 2
do
  echo "unzipping on rep1: $fq"
  gzip -d -c /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/LNCaP_DHT_RNA_rep1_${fq}.fastq.gz \
  > ~/data/samples/LNCaP_DHT_RNA_rep1_${fq}.fastq
  gzip -d -c /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/LNCaP_Veh_RNA_rep1_${fq}.fastq.gz \
  > ~/data/samples/LNCaP_Veh_RNA_rep1_${fq}.fastq
done

for fq in 1 2
do
  echo "unzipping on rep2: $fq"
  gzip -d -c /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/LNCaP_DHT_RNA_rep2_${fq}.fastq.gz \
  > ~/data/samples/LNCaP_DHT_RNA_rep2_${fq}.fastq
  gzip -d -c /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/LNCaP_Veh_RNA_rep2_${fq}.fastq.gz \
  > ~/data/samples/LNCaP_Veh_RNA_rep2_${fq}.fastq
done

#setup directory
mkdir -p ~/data/genes/
mkdir -p ~/data/genome/
#Copying reference genome
echo "Copying reference genome"
cp  /home/tzp6pz/genome_files/h38/uscs/genes/genome.gtf ~/data/genes/genome.gtf

echo "Creating .ss and .exon files for Hisat2"
mkdir -p ~/my_index
cd ~/my_index
#create ss file
python /home/tzp6pz/rna_seq/hisat2/python/hisat2_extract_splice_sites.py ~/data/genes/genome.gtf > ~/data/genes/genome.ss

#create exon file
python /home/tzp6pz/rna_seq/hisat2/pythonhisat2_extract_exons.py ~/data/genes/genome.gtf > ~/data/genes/genome.exon

module load gcc/9.2.0
module load hisat2/2.1.0

echo "Building hisat2 index using h38 "
#Build HGFM index with transcripts
hisat2-build -p $SLURM_CPUS_PER_TASK --ss  ~/data/genes/genome.ss --exon  ~/data/genes/genome.exon ~/data/genome/genome.mod.fa ~/data/indexes/genome_tran
