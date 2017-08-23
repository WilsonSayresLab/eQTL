# Single factor DE analysis using limma/voom

# /Users/hnatri/Dropbox (Personal)/TCGA_LIHC/testing/
# /Users/heini/Dropbox/TCGA_LIHC/testing/
setwd("/Users/heini/Dropbox/TCGA_LIHC/testing/")

library(limma)
library(edgeR)

counts <- read.csv("dgelist_TCGA_LIHC_count_summary.csv", header=FALSE, sep=",")
#design <- read.csv("dgelist_female_design.csv", header=TRUE, sep=",")

#design

#groups <- c("tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","normal","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","normal","tumor","tumor","normal","normal","tumor","tumor","tumor","normal","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","normal","tumor","tumor","tumor","normal","tumor","normal","tumor","normal","tumor","tumor","tumor","tumor","normal","normal","tumor","tumor","tumor","tumor","normal","normal","tumor","tumor","tumor","normal","tumor","normal","normal")
genes <- read.csv("female_genes.csv", header=TRUE)
genes <- data.frame(genes)

samplenames <- pheno$id
pheno <- read.csv("TCGA_LIHC_metadata.csv", header=TRUE, sep=",")

#pheno

tissue <- factor(pheno$tissue, levels=c("normal", "tumor"))
sex <- factor(pheno$sex, levels=c("female", "male"))

#tissue

multifactor <-factor(paste(pheno$sex,pheno$tissue,sep="."))
design <- model.matrix(~0+multifactor)
colnames(design) <- levels(multifactor)
design

#design_matrix <- model.matrix(~1 + tissue + sex)
#colnames(design_matrix) <- c("normal", "tumor", "female", "male")
#design_matrix

dge <- DGEList(counts=counts, genes=genes)
colnames(dge) <- samplenames
dge$samples$tissue <- tissue
dge$samples$sex <- sex
v <- voom(dge, design, normalize="quantile")
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTreat(fit, coef=ncol(design), n=Inf)

contrast_design <- makeContrasts(female.tumor-male.tumor, levels=design)

contrast_design

fit2 <- contrasts.fit(fit, contrast_design)
fit2 <- eBayes(fit2)

############
# above, case - control means to contrast (compare) of case from control and this will instruct the limma to perform the differential analysis. After that, you will ask limma to perform the analysis as shown below.

# and then, you can get the result by using

fit3 <- treat(fit2, lfc=log2(1.2))
topTreat(fit3, coef=ncol(contrast_design), n=Inf)

dt <- decideTests(fit3)
summary(dt)

write.fit(fit3, dt, file="TCGA_LIHC_femaletumor_maletumor.txt")

tops <- topTreat(fit3, number=Inf) # get all genes 
tops[which(tops$logFC > 0), ] [1:50,] # up reg top 50 
tops[which(tops$logFC < 0), ] [1:50,] # down reg top 50

plotMD(fit3, column=1, status=dt[,1], main=colnames(fit3)[1], 
       xlim=c(-8,13))

library(gplots)
toptable <- topTreat(fit3, coef=ncol(contrast_design), n=Inf)
toptable
topgenes <- toptable$gene[1:100]
i <- which(v$genes$gene %in% topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale="row",
          labRow=v$genes$gene[i], labCol=tissue, 
          col=mycol, trace="none", density.info="none", 
          margin=c(5,5), dendrogram="column",
          cexRow=0.5, keysize=1) # lhei=c(2,10)


# following code limits the lowest and highest color to 5%, and 95% of your range, respectively
quantile.range <- quantile(v$E[i,], probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)

# use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
color.palette  <- colorRampPalette(c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0"))(length(palette.breaks) - 1)

heatmap.2(
  v$E[i,],
  col    = color.palette,
  breaks = palette.breaks,
  labRow=v$genes$gene[i], labCol=tissue, trace="none", density.info="none", 
  margin=c(5,5), dendrogram="column",
  cexRow=0.5, key=FALSE, lwid=c(0.1,4), lhei=c(0.4,4)
)

# Gradient bar
palette.breaks <- seq(-6, 6, 0.001)
color.palette  <- colorRampPalette(c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0"))(length(palette.breaks) - 1)

plot(NULL, xlim = c(-6, 6), ylim = c(0, 1), bty = "n", yaxt = "none", ylab = "", xlab = "")
segments(palette.breaks[-length(palette.breaks)], rep(0, length(palette.breaks) - 1), palette.breaks[-length(palette.breaks)], rep(1, length(palette.breaks) - 1), lwd = 2, col = color.palette)

# volcano plot

volcanoplot(fit2, coef = ncol(contrast_design), style = "p-value", highlight = 20, names = fit$genes$gene,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

plotMDS(v, top = 100, plot=TRUE, pch = 19, cex = 1,
        ndim = 2) #, gene.selection = "pairwise", labels=tissue

top1000 <- toptable$gene[1:1000]
top1000_file <- "femaletumor_maletumor_top1000.csv"
write.table(top1000, top1000_file, append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = FALSE)
