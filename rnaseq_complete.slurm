#!/bin/bash
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --cpus-per-task=28
#SBATCH --mem 256GB
#SBATCH -o ~/output_complete_analysis.txt
#SBATCH -p standard
#SBATCH -A zanglab

#must be done in asscoiated path

#hisat2 setup and building index
echo "Hisat2 Setup and Building Index using h39 USCS format"
sh hisat_setup.sh

echo " "

echo "Fasqc, mapping, and sorting each sample"
for fq in ~/data/samples/*_1.fastq
do
  echo "running analysis on pair: $fq"
  sh rna_hisat2.sh $fq $SLURM_CPUS_PER_TASK
done

echo " "

echo "Running code to get files needed for Differential Expression"

echo "  "

echo "Running code for Stringtie Merge (Part 1)"
for fq in ~/my_index/map/*.bam
do
  echo "running analysis on pair: $fq"
  sh rna_stringtie_merge1.sh $fq $SLURM_CPUS_PER_TASK
done

echo "Running code for Stringtie Merge (Part 2)"
sh rna_stringtie_merge2.sh $SLURM_CPUS_PER_TASK

echo "Running code for Stringtie to get ready to export to R in Ballgown (Part 3)"

for fq in ~/my_index/map/*_sorted.bam
do
  echo "running analysis on pair: $fq"
  sh rna_stringtie_part3.sh $fq $SLURM_CPUS_PER_TASK
done

echo "  "
echo "  "
echo "  "
echo "Need to perform analysis in R via Ballgown using rna_ballgown.slurm"


echo "SLURM FINISHED"

#
