#r script file to get chip-seq quality control measures

#reference for using the ChIPQC R package
#https://hbctraining.github.io/Intro-to-ChIPseq/lessons/06_combine_chipQC_and_metrics.html

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("DiffBind")
#BiocManager::install("ChIPQC")


#for reference
#/scratch/tzp6pz/proj1/work
# /scratch/tzp6pz/proj1/rawdata_LNCaP_AR/LNCaP_DHT_AR_1.fastq
# align_DHT_AR_1.sam
# /scratch/tzp6pz/proj1/chip_macs/AR_1_peaks.narrowPeak

library(ChIPQC)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("BiocParallel")
register(SerialParam())


samples = read.csv('/home/billy/Documents/start_proj/2_12_2021/chip_v2/qc/qc_file.csv')

exampleExp <-ChIPQC(samples, annotation="hg38")#, consensus=TRUE)
QCmetrics(exampleExp)
