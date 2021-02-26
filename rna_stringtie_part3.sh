#!/bin/bash/

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1
cpus=$2

# grab base of filename for naming outputs
base=`basename ${fq} .bam`
echo "Sample name is $base"
echo "Number of CPUS: $cpus"

# directory with bowtie genome files
genome=~/data/indexes/genome_tran
# make all of the output directories
# The -p option means mkdir will create the whole path if it
# does not exist and refrain from complaining if it does exist
mkdir -p ~/my_index/assembly
mkdir -p ~/my_index/map

## set up file names
align_bam=~/my_index/map/${base}.bam
merge_out=~/my_index/assembly

#ballgown files location
ball_out=~/my_index/ballgown/${base}
mkdir -p $ball_out

# set up the software environment
module load gcc/9.2.0
module load stringtie/2.1.0

###
###NEED TO START AT ASSEMBLY SECTION OF this
# https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/
echo "Doing final Stringtie operation"
stringtie -B -e -p $cpus -G ${merge_out}/stringtie_merged.gtf -o ${ball_out}/${base}.gtf $align_bam

#
