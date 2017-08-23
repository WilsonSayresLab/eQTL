# Single factor DE analysis using limma/voom

# /Users/hnatri/Dropbox (Personal)/TCGA_LIHC/testing/
# /Users/heini/Dropbox/TCGA_LIHC/testing/
setwd("/Users/hnatri/Dropbox (Personal)/TCGA_LIHC/testing/")

library(limma)
library(edgeR)

counts <- read.csv("dgelist_female_counts.csv", header=FALSE, sep=",")
design <- read.csv("dgelist_female_design.csv", header=TRUE, sep=",")

design

groups <- c("tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","normal","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","normal","tumor","tumor","normal","normal","tumor","tumor","tumor","normal","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","normal","tumor","tumor","tumor","normal","tumor","normal","tumor","normal","tumor","tumor","tumor","tumor","normal","normal","tumor","tumor","tumor","tumor","normal","normal","tumor","tumor","tumor","normal","tumor","normal","normal")
genes <- read.csv("female_genes.csv", header=TRUE)
genes <- data.frame(genes)

pheno <- read.csv("female_metadata.csv", header=TRUE, sep=",")

#pheno

tissue <- factor(pheno$tissue, levels=c("normal", "tumor"))

#tissue

design_matrix <- model.matrix(~0 + tissue)
colnames(design_matrix) <- c("normal", "tumor")
design_matrix

dge <- DGEList(counts=counts, group=factor(groups), genes=genes)
v <- voom(dge, design_matrix, normalize="quantile")
fit <- lmFit(v, design_matrix)
fit <- eBayes(fit)
topTreat(fit, coef=ncol(design_matrix))

contrast_design <- makeContrasts(tumor-normal, levels=design_matrix)

contrast_design

fit2 <- contrasts.fit(fit, contrast_design)
fit2 <- eBayes(fit2)
topTable(fit2)

############
# above, case - control means to contrast (compare) of case from control and this will instruct the limma to perform the differential analysis. After that, you will ask limma to perform the analysis as shown below.

# and then, you can get the result by using

fit2 <- treat(fit2, lfc=log2(1.2))
topTreat(fit2, coef=ncol(contrast_design))

dt <- decideTests(fit2)
summary(dt)

write.fit(fit2, dt, file="results.txt")

# this "makeContrasts" function can also be used if you would like to compare a few conditions together with respect to a few other conditions.
# Ex. contrast_design <- makeContrats( (s1+s2) - (s3+s4), level= design) when the "design" contains s1, s2, s3, and s4 conditions
