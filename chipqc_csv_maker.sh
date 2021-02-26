#!/bin/sh
echo SampleID, Factor, bamReads, CountrolID, bamControl, Peaks, PeakCaller, Tissue, Condition

for f in 1 2
do
   echo DHT.Rep$f, DHT, $f, /scratch/tzp6pz/chip/results/bowtie2/LNCaP_DHT_AR_${fq}_aln.bam, \
   DHT-Input$f, /scratch/tzp6pz/chip/results/bowtie2/LNCaP_Veh_AR_${fq}_aln.bam, \
    /scratch/tzp6pz/chip/macs/AR_${f}_peaks.xls, macs, NA, NA
done
