#R code for ballgown analysis on RNA-seq data
#adopted from
#https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5032908/
library(ballgown)#for interpreting resutls from stringtie etc.
library(genefilter)#additional gene tools
library(dplyr)#data cleaning functions
library(RColorBrewer)#for nice color options
library(gplots)#for nicer heatmap
library(xtable)#for making tables in latex

# change this to the directory that contains all the StringTie results
setwd("/scratch/tzp6pz/rna2/my_index/")

# create a ballgown object
#bg = ballgown(samples='ballgown', meas='all')

# load the sample information
txt_data <- read.table("ballgown/ballgown.txt", sep=',')
#add phenotypes
txt_data[,3]<-c(1, 2, 1, 2)
txt_data[,4]<-c(1, 1, 0, 0)
colnames(txt_data)<-c('sample', 'ids', 'rep', 'treatment')
#remove whitespace
txt_data[,2]<-gsub(" ","", txt_data$ids, fixed=TRUE)
#txt_data<-txt_data[,-2]
# create a ballgown object
bg <- ballgown(samples=txt_data[,1], pData=txt_data[,-1])
#
#expression matrix
gene_expression = as.data.frame(gexpr(bg))
colnames(gene_expression)<-c("DHT1","DHT2","Veh1","Veh2")
#creating heatmap
#reference https://rpubs.com/tgjohnst/heatmaps_testing_1
jpeg(file="/scratch/tzp6pz/rna2/results/r_ballgown/gene_expression_heatmap.jpg",
    width=2000,
    height=2000
    )
heatmap.2(as.matrix(gene_expression),
          col=rev(brewer.pal(9,"RdBu")),#change color
          scale="row",
          Colv = "NA",#do not cluster by column
          cexRow=2,#bigger font
          cexCol=2,#bigger font
          trace="none" #remove histogram
          )
dev.off()


#generating plots
#code inspired by https://rnabio.org/module-03-expression/0003/04/01/DE_Visualization/
min_nonzero=1

data_columns=c(1:4)
short_names=c("DHT1","DHT2","Veh1","Veh2")
colours()
data_colors=c("tomato1","tomato2","royalblue1","royalblue2")

#
jpeg(file="/scratch/tzp6pz/rna2/results/r_ballgown/gene_epxression_boxplot.jpg",
    width=1000,
    height=1000
    )
boxplot(log2(gene_expression[,data_columns]+min_nonzero),
        col=data_colors,
        names=short_names,
        las=2,
        ylab="log2(FPKM)",
        main="Distribution of FPKMs for all 4 Samples")

dev.off()

#differential expression
#reference https://rstudio-pubs-static.s3.amazonaws.com/289617_cb95459057764fdfb4c42b53c69c6d3f.html
stat_results = stattest(bg, feature='gene', meas='FPKM',getFC=TRUE, covariate='treatment')
head(stat_results)
bg_filt = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)
bg_table = texpr(bg_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])
stat_results = merge(stat_results,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

#up/down regulated genes
#reference https://support.bioconductor.org/p/77144/
sig=which(stat_results$pval<0.05)
stat_results[,"de"] = log2(stat_results[,"fc"])

jpeg(file="/scratch/tzp6pz/rna2/results/r_ballgown/up_down_hist.jpg",
    width=500,
    height=500
    )
hist(stat_results[sig,"de"],
      breaks=50,
      col="seagreen",
      xlab="log2(Fold change) DHT vs Veh",
      main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)
dev.off()

#table with output
sigpi = which(stat_results[,"pval"]<0.05)
sigp = stat_results[sigpi,]
sigde = which(abs(sigp[,"de"]) >= 2)
sig_tn_de = sigp[sigde,]

o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE)
output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
write.table(output, file="/scratch/tzp6pz/rna2/results/r_ballgown/SigDE.txt", sep="\t", row.names=FALSE, quote=FALSE)
#View selected columns of the first 25 lines of output
#output[1:25,c(1,4,5)]
xtable(output, digits=4)

#% latex table generated in R 4.0.0 by xtable 1.8-4 package
#% Wed Feb 24 11:42:08 2021
#\begin{table}[ht]
#\centering
#\begin{tabular}{rllrrrr}
#  \hline
# & gene\_name & id & fc & pval & qval & de \\
#  \hline
#3053 & TMPRSS2 & MSTRG.8523 & 5.4402 & 0.0029 & 0.3725 & 2.4437 \\
#  714 & NKX3-1 & MSTRG.12312 & 4.8245 & 0.0004 & 0.3725 & 2.2704 \\
#   \hline
#\end{tabular}
#\end{table}

#
