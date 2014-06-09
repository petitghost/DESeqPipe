####################################################################
#Date: 2014-05-21
#Author: ghost
#Usage: R -f ~/bin/DEseq.R
#Description: Run DEseq to find differential expression genes
#Note: The data no replicate on each condition 
####################################################################

library("DEseq")
x=read.table("s1_s2.readCount",h=T,row.names=1)
#Add the treatment information that DESeq knows what data each column represents
condition=factor(c("s1","s2"))
# Central data structure in the DESeq
cds=newCountDataSet(x,condition) #Create a CountDataSet object
# Normalisation
cds=estimateSizeFactors(cds) #Estimate the size factors for a CountDataSet
# Variance estimation
cds=estimateDispersions(cds,method="blind",sharingMode="fit-only") #Estimate and fit dispersions for a CountDataSet
# Calling differential expression
res=nbinomTest(cds,"s1","s2") #Test for differences between the base means for two conditions
png("MA.png")
plotMA(res) #Makes a MA-plot (default: Adjustment p-val < 0.1)
#調整：plotMA(res,col = ifelse(res$padj>=0.1, "gray32", "red3"))
write.table(res,file="DEseq_all_gene.txt") #Output a table for all genes
resSig=res[res$padj<0.1, ] # Filter according to some chosen threshold for the false dicovery rate (FDR)
write.table(resSig[order(resSig$pval),],file="DEseq_significant_gene.txt") #Output a table for significant genes
png("p-value_his.png")
hist(res$pval,col=rainbow(20),breaks=100)
dev.off()
