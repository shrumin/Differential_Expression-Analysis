setwd("~/Downloads")

##Design experiment table
colTable <- read.csv("Design_Exp.csv",row.names=1)

##count table
countTable <- read.csv("gene_count_matrix.csv",row.names=1)

###################################################################################################
##Deseq analysis

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=countTable,colData=colTable,design= ~ group)
notAllZero <- (rowSums(counts(dds))>0) 
dds <- dds[notAllZero,]

##differential expression analysis
dds <- DESeq(dds)

###################################################################################################
#visualise data distribution

##dispersion plot
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()
plotDispEsts(dds, main="Dispersion plot")

#regularised log transformation
rld <- rlogTransformation(dds)
lgc.raw <- log2(counts(dds, normalized=FALSE))
lgc.norm <- log2(counts(dds, normalized=TRUE))
par(mfrow = c(1,3))
boxplot(lgc.raw,main="Raw-counts") 
boxplot(lgc.norm,main="Normalised-counts")
boxplot(assay(rld),main="rlog-transformed")


#density plot
library(limma)
plotDensities(lgc.raw,group=dds$group,legend="topright",main="RawDensity")
plotDensities(lgc.norm,group=dds$group,legend="topright",main="NormalisedDensity")
plotDensities(assay(rld),group=dds$group,legend="topright",main="rldDensity")

#pca plot

plotPCA(rld, intgroup=c("group"))

###################################################################################################
#DE result extraction

#AvsB
res_A_B <- results(dds,contrast=c("group","A","B"))
res_A_B_2 <- results(dds,contrast=c("group","A","B"),lfcThreshold=2)
summary(res_A_B,alpha=0.05)
summary(res_A_B_2,alpha=0.05)
#AvsC  
res_A_C <- results(dds,contrast=c("group","A","C"))
res_A_C_2 <- results(dds,contrast=c("group","A","C"),lfcThreshold=2)
summary(res_A_C,alpha=0.05)
summary(res_A_C_2,alpha=0.05)
#BvsC
res_B_C <- results(dds,contrast=c("group","B","C"))
res_B_C_2 <- results(dds,contrast=c("group","B","C"),lfcThreshold=2)
summary(res_B_C,alpha=0.05)
summary(res_B_C_2,alpha=0.05)

##arranging by P adj
res_A_B_sort=res_A_B[order(res_A_B$padj),]
res_A_C_sort=res_A_C[order(res_A_C$padj),]
res_B_C_sort=res_B_C[order(res_B_C$padj),]

##SIGNIFICANT ONES
res_A_B_sig<- subset(res_A_B_sort, 0.001)
write.csv(res_A_B_sig, "/Users/shrumin/Downloads/sig_A_B.csv", row.names=TRUE)
res_A_C_sig<- subset(res_A_C_sort, 0.001)
write.csv(res_A_C_sig, "/Users/shrumin/Downloads/sig_A_C.csv", row.names=TRUE)
res_B_C_sig<- subset(res_B_C_sort, 0.001)
write.csv(res_B_C_sig, "/Users/shrumin/Downloads/sig_B_C.csv", row.names=TRUE)
###################################################################################################
##MA PLOT
resNorm <- lfcShrink(dds, coef="group_B_vs_A", type="normal")
plotMA(resNorm, ylim=c(-2,2), main='M', alpha=0.1)
resNorm1 <- lfcShrink(dds, coef="group_C_vs_A", type="normal")
plotMA(resNorm1, ylim=c(-2,2), main='M', alpha=0.1)

###################################################################################################

