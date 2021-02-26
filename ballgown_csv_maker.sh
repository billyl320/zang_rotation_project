#!/bin/sh
loc=/scratch/tzp6pz/rna2/my_index/ballgown/
#loc=L
for fq in ${loc}*_sorted.gtf
do
  base=`basename ${fq} .gtf`
  echo ${loc}$base, $base
done
