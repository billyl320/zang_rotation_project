#!/bin/bash
#SBATCH --job-name=macs2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH -o /scratch/tzp6pz/chip/slurm_files/macs.txt
#SBATCH --partition=standard
#SBATCH --account=zanglab

mkdir -p /scratch/tzp6pz/chip/qc/macs/

#Load macs2
module load macs2/2.2.7.1
module load intel/18.0
module load goolf/7.1.0_3.1.4  R/3.6.3

for fq in 1 2
do
  echo "running analysis on AR $fq"
  macs2 callpeak -t /scratch/tzp6pz/chip/results/bowtie2/LNCaP_DHT_AR_${fq}_aln.bam \
                 -c /scratch/tzp6pz/chip/results/bowtie2/LNCaP_Veh_AR_${fq}_aln.bam \
                 -f BAM \
                 -g hs \
                 -q 0.01 \
                 -n AR_$fq \
                 --outdir /scratch/tzp6pz/chip/macs
done
# see here for more potential guidance; didn't follow this too much for this step
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html

#create csv file
echo "Created CSV file for R script"
mkdir -p /scratch/tzp6pz/chip/qc/
sh /scratch/tzp6pz/chip/slurm_files/chipqc_csv_maker.sh > /scratch/tzp6pz/chip/qc/output.csv

#run R script for QC measures
echo "Running R script for QC measures"
R CMD BATCH /scratch/tzp6pz/chip/qc/r_chipqc.r /scratch/tzp6pz/chip/qc/r_chipqc_output.txt

#
