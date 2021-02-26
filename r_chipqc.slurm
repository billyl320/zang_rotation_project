#!/bin/bash
#SBATCH --job-name=macs2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH -o /scratch/tzp6pz/chip/slurm_files/macs.txt
#SBATCH --partition=standard
#SBATCH --account=zanglab

#Load macs2
module load macs2/2.2.7.1

for fq in 1 2
do
  echo "running analysis on AR $fq"
  macs2 callpeak -t /scratch/tzp6pz/chip/results/bowtie2/align_DHT_AR_$fq.bam \
                 -c /scratch/tzp6pz/chip/results/bowtie2/align_Veh_AR_$fq.bam \
                 -f BAM \
                 -g hs \
                 -q 0.01 \
                 -n AR_$fq \
                 --outdir /scratch/tzp6pz/chip/macs
done
# see here for more potential guidance; didn't follow this too much for this step
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html
