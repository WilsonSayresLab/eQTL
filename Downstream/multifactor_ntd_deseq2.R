# DESeq PCA and differential expression analysis of gene count data

setwd("/Users/hnatri/Dropbox (Personal)/TCGA_LIHC/testing/")

library(tidyverse)

mycounts <- read_csv("TCGA_LIHC_count_summary.csv")
metadata <-  read_csv("TCGA_LIHC_metadata.csv")

mycounts <- as.data.frame(mycounts)
metadata <- as.data.frame(metadata)

names(mycounts)[-1]==metadata$id

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData=mycounts, 
                              colData=metadata, 
                              design=~sex+tissue+sex:tissue, 
                              tidy=TRUE)

# Prefiltering

dds <- dds[ rowSums(counts(dds)) > 1, ]

# Setting reference levels
dds$sex <- factor(dds$sex, levels = c("female","male"))
dds$tissue <- factor(dds$tissue, levels = c("normal","tumor"))

# Run the DESeq pipeline
dds <- DESeq(dds)

res <- results(dds, contrast=c("sex", "female", "male"), tidy=TRUE)
res <- tbl_df(res)
res


# provide the dds object and the number of the coefficient we want to moderate. If a results object is provided, the log2FoldChange column will be swapped out, otherwise lfcShrink returns a vector of shrunken log2 fold changes.
#resLFC <- lfcShrink(dds, coef=2, res=res)
#resLFC

# reorder results based on P value
resOrdered <- res[order(res$pvalue),]
summary(res)

# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

#Note that the results function automatically performs independent filtering based on the mean of normalized counts for each gene, optimizing the number of genes which will have an adjusted p value below a given FDR cutoff, alpha. Independent filtering is further discussed below. By default the argument alpha is set to 0.1. If the adjusted p value cutoff will be a value other than 0.1, alpha should be set to that value:
res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)

#library("IHW")
#resIHW <- results(dds, filterFun=ihw)
#summary(resIHW)

#sum(resIHW$padj < 0.1, na.rm=TRUE)




# Visualization

### In DESeq2, the function plotMA shows the log2 fold changes 
### attributable to a given variable over the mean of normalized counts.


#plotMA(res, main="DESeq2", ylim=c(-2,2))


### A simple function for making this plot is plotCounts, which normalizes counts by sequencing depth and adds a pseudocount
### of 1 2 to allow for log scale plotting. The counts are grouped by the variables in intgroup, where more than
### one variable can be specified. Here we specify the gene which had the smallest p value from the results table
### created above. You can select the gene to plot by rowname or by numeric index.

plotCounts(dds, gene=which.min(res$padj), intgroup=c("sex", "tissue"))


### For customized plotting, an argument returnData specifies that the function 
### should only return a data.frame for plotting with ggplot.

d <- plotCounts(dds, gene=which.min(res$padj), intgroup=c("sex", "tissue"),
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=sex, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
scale_y_log10(breaks=c(25,100,400))

#plotMA(res, ylim=c(-2,2),main="DESeq2")

#plotMA(resLFC, ylim=c(-2,2))

#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]

plotCounts(dds, gene=which.min(res$padj), intgroup=c("sex", "tissue"))

write.csv(as.data.frame(resOrdered), 
          file="tumor_sex_results.csv")



# this gives log2(n + 1)
ntd <- normTransform(dds)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("sex")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



#library("vsn")
#meanSdPlot(assay(ntd))




# heatmap of sample to sample distances
sampleDists <- dist(t(assay(ntd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ntd$sex, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(ntd, intgroup=c("sex"))


# SIMPLE

## Principal component analysis

plotPCA(ntd, intgroup=c("sex", "tissue"))

## Heatmap of sample distances
library("gplots")   # If this fails, run: install.packages("gplots")
library("RColorBrewer")
sampleDists <- dist(t(assay(ntd)))
sampleDistMatrix <- as.matrix( sampleDists )
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
heatmap.2(sampleDistMatrix, trace="none", col=colours)

## Heatmap of 100 most variable genes
#library("genefilter")
topVarGenes <- head(order(rowVars(assay(ntd)), decreasing=TRUE), 35)
#heatmap.2(assay(ntd)[topVarGenes, ], scale="column",
#trace="none", dendrogram="column", margins=c(5, 10),
#col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

# with pheatmap
pheatmap(assay(ntd)[topVarGenes, ], scale="column",
trace="none", dendrogram="column", margins=c(5, 10),
col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))
