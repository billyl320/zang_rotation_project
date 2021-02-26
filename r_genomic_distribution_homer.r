#r code to calculate each type of binding sites

library(xtable)#for making latex tables

#reading in outputfile from homer's annotatePeaks.pl function (analyzing CHiP-seq data from macs)
out<-read.table(file='outputfile.txt', header=TRUE, sep='\t', fill=TRUE)

#answer to question 4 subquestion 2
table(substr(out[,8], 1, 5))/dim(out)[1]

xtable(table(substr(out[,8], 1, 5))/dim(out)[1])

#answer to question 4 subquestion 3
#true = upstream (promoter)
#false = downstream (enhancer)
table(as.numeric(as.character(out[,10]) )>0 )
table(as.numeric(as.character(out[,10]) )>0 )/sum(table(as.numeric(as.character(out[,10]) )>0 ))

xtable(table(as.numeric(as.character(out[,10]) )>0 )/sum(table(as.numeric(as.character(out[,10]) )>0 )))

#getting results from rna_ballgown.r results
gene1<-which(as.character(out[,16])=='TMPRSS2')
#[1]  604  986 1131 3919 > for AR_1
gene2<-which(as.character(out[,16])=='NKX3-1')
#[1]  157 1918 2168 > for AR_1

#answering questions in question 4 subpart 4
#gene1
table(substr(out[gene1,8], 1, 5))/dim(out[gene1,])[1]
table(as.numeric(as.character(out[gene1,10]) )>0 )
table(as.numeric(as.character(out[gene1,10]) )>0 )/sum(table(as.numeric(as.character(out[gene1,10]) )>0 ))

#gene2
table(substr(out[gene2,8], 1, 5))/dim(out[gene2,])[1]
table(as.numeric(as.character(out[gene2,10]) )>0 )
table(as.numeric(as.character(out[gene2,10]) )>0 )/sum(table(as.numeric(as.character(out[gene2,10]) )>0 ))



#
