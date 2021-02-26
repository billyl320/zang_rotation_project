#!/bin/sh
#loc=/scratch/tzp6pz/rna2/my_index/assembly/L
loc=L
for fq in ${loc}*.gtf
do
  echo $fq
done
