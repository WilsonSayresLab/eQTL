# DESeq PCA and differential expression analysis of gene count data obtained with HTSeq-counts

library(DEseq2)
library(EDASeq)
library(ggplot2)
library(topGO)

run_deseq <- function(HTSEQ_COUNT_TABLE) {

  # Converting a HTSeq count table to a DESeqDataSet
  design <- "~ condition" # condition tumor or normal
  DESeq2Table <- DESeqDataSetFromHTSeqCount(HTSEQ_COUNT_TABLE, directory = ".", design,
                                            ignoreRank = FALSE, ...)
  # Observing the data
  head(assay(DESeq2Table))
  colSums(assay(DESeq2Table))
  colData(DESeq2Table)
  rowData(DESeq2Table)
  
  # Note that the rowData slot is a GRangesList, which contains all the information 
  # about the exons for each gene, i.e., for each row of the count matrix. It also 
  # contains additional annotation for each gene, i.e. a gene description and a gene symbol.
  
  mcols(rowData(DESeq2Table))
  
  # The colData slot contains all the sample metadata.
  
  colData(DESeq2Table)
  
  # Quality control and normalization of data
  # Let’s check how many genes we capture by counting the number of genes that have non–zero counts in all samples.
  # In a typical RNA-Seq experiment, there will be at least several thousand genes that are expressed in all samples. 
  # If the number of non–zero genes is very low, there is usually something wrong with at least some of the samples 
  # (e.g. the sample preparation was not done properly, or there was some problem with the sequencing).
  
  GeneCounts <- counts(DESeq2Table)
  idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
  sum(idx.nz)
  
  # Normalization
  # Relevel the condition factor to make sure that we have the normal samples as a base level.
  # This will allow us to obtain the fold changes always as normal–tumor.
  colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("Normal", "Tumor"))
  colData(DESeq2Table)$condition <- factor(ifelse(is.na(colData(DESeq2Table)$condition),  "Tumor", "Normal"), levels = c("Normal", "Tumor"))
  
  # We now estimate estimate the size factors an print them
  
  DESeq2Table <- estimateSizeFactors(DESeq2Table)
  sizeFactors(DESeq2Table)
  
  # To assess whether the normalization has worked, we plot the densities of counts for 
  # the different samples. Since most of the genes are (heavily) affected by the experimental 
  # conditions, a succesful normalization will lead to overlapping densities.
  
  multidensity( counts(DESeq2Table, normalized = T)[idx.nz ,],
                xlab="mean counts", xlim=c(0, 1000))
  
  # If the normalization has worked, the ECDFs of the different samples should be overlapping.
  
  multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],
             xlab="mean counts", xlim=c(0, 1000))
  
  # To further assess systematic differences between the samples, we can also plot pairwise 
  # mean–average plots: We plot the average of the log–transformed counts vs the fold change 
  # per gene for each of the sample pairs.
  
  # The function combn helps us to automatically generate all sample pairs and the function 
  # MDPlot from the EDASeq package then generates the the pairwise MA plots.
  
  pdf("pairwiseMAs.pdf")
  MA.idx = t(combn(1:8, 2))
  for( i in  seq_along( MA.idx[,1])){ 
    MDPlot(counts(DESeq2Table, normalized = T)[idx.nz ,], 
           c(MA.idx[i,1],MA.idx[i,2]), 
           main = paste( colnames(DESeq2Table)[MA.idx[i,1]], " vs ",
                         colnames(DESeq2Table)[MA.idx[i,2]] ), ylim = c(-3,3))
  }
  dev.off()
  
  # PCA and sample heatmaps
  # log-transformation withDESeq2 regularized–logarithm transformation rlog.
  # For genes with high counts, the rlog transformation differs not much from an ordinary log2 transformation.
  # For genes with lower counts, however, the values are shrunken towards the genes’ averages 
  # across all samples. Using an empirical Bayesian prior in the form of a ridge penalty, this 
  # is done such that the rlog–transformed data are approximately homoskedastic. Note that the rlog 
  # transformation is provided for applications other than differential testing. For differential 
  # testing it is always recommended to apply the DESeq function to raw counts.
  
  # Note the use of the function to transpose the data matrix. We need this because dist 
  # calculates distances between data rows and our samples constitute the columns. We visualize 
  # the distances in a heatmap, using the function heatmap.2 from the r CRANpkg("gplots") package.
  
  # Note that we have changed the row names of the distance matrix to contain the genotype and the 
  # column to contain the sample ID, so that we have all this information in view when looking at the heatmap.
  
  rld <- rlogTransformation(DESeq2Table, blind=TRUE)
  
  
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  rownames(mat) <-  colData(rld)$condition
  colnames(mat) <-  colData(rld)$sampleNO
  
  hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
  heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
  
  # PCA for 500 genes showing the highest variance
  
  ntop = 500
  
  Pvars <- rowVars(assay(rld))
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                        length(Pvars)))]
  
  PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  # Plotting principal components with qplot
  
  dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                      PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                      sampleNO = colData(rld)$sampleNO,
                      condition = colData(rld)$condition)
  
  (qplot(PC1, PC2, data = dataGG, color =  condition, 
         main = "PC1 vs PC2, top variable genes", size = I(6))
    + labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
           y = paste0("PC2, VarExp:", round(percentVar[2],4)))
    + scale_colour_brewer(type="qual", palette=2)
  )
  
  # Remocving outliers
  
  outliers <- as.character(subset(colnames(DESeq2Table), dataGG$PC1 > 0))
  outliers
  
  DESeq2Table <- DESeq2Table[, !(colnames(DESeq2Table) %in% outliers)] 
  
  # Cluster heatmap and the PCA plot without the outliers.
  
  rld <- rlogTransformation(DESeq2Table, blind=TRUE)
  
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  
  rownames(mat) <-  colData(rld)$condition
  colnames(mat) <-  colData(rld)$sampleNO
  
  hmcol <- colorRampPalette(brewer.pal(9, "OrRd"))(255)
  heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
  
  Pvars <- rowVars(assay(rld))
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                        length(Pvars)))]
  
  PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                      PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                      sampleNO = colData(rld)$sampleNO,
                      condition = colData(rld)$condition)
  
  (qplot(PC1, PC2, data = dataGG, color =  condition, 
         main = "PC1 vs PC2, top variable genes, no outliers", size = I(6))
    + labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
           y = paste0("PC2, VarExp:", round(percentVar[2],4)))
    + scale_colour_brewer(type="qual", palette=4)
  )
  
  # Differential expression analysis with Negative–Binomial-distribution (NB) model
  
  # The first step in the analysis of differential expression, is to obtain an estimate 
  # of the dispersion parameter for each gene
  
  # For weak genes, the Poisson noise is an additional source of noise, which is added 
  # to the dispersion. The function visualizes DESeq2’s dispersion estimates:
    
  DESeq2Table <- estimateDispersions(DESeq2Table)
  plotDispEsts(DESeq2Table)
  
  # Statistical testing of Differential expression
  # The test–performed is Wald test, which is a test for coefficients in a regression model.
  
  DESeq2Table <-  nbinomWaldTest(DESeq2Table)
  DESeq2Res <- results(DESeq2Table, pAdjustMethod = "BH")
  
  table(DESeq2Res$padj < 0.1)
  
  # Independent filtering
  # DESeq2 performs automated independent filtering: For weakly expressed genes, we have 
  # no chance of seeing differential expression, because the low read counts suffer from 
  # so high Poisson noise that any biological effect is drowned in the uncertainties from 
  # the read counting.
  
  # The DESeq2 software automatically performs independent filtering which maximizes the 
  # number of genes which will have a BH–adjusted \(p\)–value less than a critical value 
  # (by default, alpha is set to 0.1). This automatic independent filtering is performed by, 
  # and can be controlled by, the results function. We can observe how the number of rejections 
  # changes for various cutoffs based on mean normalized count. The following optimal threshold 
  # and table of possible values is stored as an attribute of the results object.
  
  plot(metadata(DESeq2Res)$filterNumRej, type="b", xlab="quantiles of 'baseMean'",
       ylab="number of rejections")
  
  # Inspection and correction of \(p\)–values
  # A histogram of \(p\)–values should always be plotted in order to check whether they have 
  # been computed correctly. We also do this here:
  
  hist(DESeq2Res$pvalue, col = "lavender", main = "WT vs Deletion", xlab = "p-values")
  
  # If the histogram is hill shaped:

  # We first remove genes filtered out by independent filtering and the dispersion outliers, 
  # they have NA adj. pvals and NA \(p\)–values respectively.
  
  DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]
  
  DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]
  
  # We now remove the original adjusted \(p\)–values, since we will add the corrected ones later on.
  
  DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
  
  # We can now use z–scores returned by DESeq2as input to fdrtool to re–estimate the \(p\)–values
  
  FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)
  
  FDR.DESeq2Res$param[1, "sd"]
  
  # We now add the new BH–adjusted \(p\)–values add values to the results data frame.
  
  DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
  
  # We can now plot the histogram of the “correct” \(p\)–values.
  
  hist(FDR.DESeq2Res$pval, col = "royalblue4", 
       main = "WT vs Deletion, correct null model", xlab = "CORRECTED p-values")
  
  # Extracting differentially expressed genes
  
  table(DESeq2Res[,"padj"] < 0.1)
  
  plotMA(DESeq2Res)
  
  # Check a validated gene
  # A quick way to visualize the counts for a particular gene is to use the plotCounts 
  # function, which takes as arguments the DESeqDataSet, a gene name, and the group 
  # over which to plot the counts
  
  
  # genexyz <- subset(anSig, SYMBOL == "genexyz")$ENSEMBL
  # plotCounts(DESeq2Table, gene=genexyz, intgroup=c("condition"))
  
  # Gene ontology enrichment analysis
  
  # We first get average gene expressions for each of the genes and then find non–DE 
  # genes that show a similar expression as the DE–genes. These genes are then our background.
  
  overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])
  
  sig_idx <- match(anSig$ENSEMBL, rownames(overallBaseMean))
  
  backG <- c()
  
  for(i in sig_idx){
    ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
    backG <- c(backG, ind)
    
  }
  
  backG <- unique(backG)
  backG <- rownames(overallBaseMean)[backG]
  
  # We now remove DE genes from background and the get the total number of genes in the background.
  
  backG <- setdiff(backG,  anSig$ENSEMBL)
  length(backG)
  
  # Plotting the density of the average expressions, shows that the background matching has worked reasonably well.
  
  multidensity( list( 
    all= log2(DESeq2Res[,"baseMean"]) ,
    foreground =log2(DESeq2Res[anSig$ENSEMBL, "baseMean"]), 
    background =log2(DESeq2Res[backG, "baseMean"])), 
    xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")
  
  # Running topGO
  
  # We first create a factor alg which indicates for every gene in our universe (union of background and DE-genes), 
  # whether it is differentially expressed or not. It has the ENSEMBL ID’s of the genes in our universe as names 
  # and contains 1 if the gene is DE and 0 otherwise.
  
  onts = c( "MF", "BP", "CC" )
  
  geneIDs = rownames(overallBaseMean)
  inUniverse = geneIDs %in% c(anSig$ENSEMBL,  backG) 
  inSelection =  geneIDs %in% anSig$ENSEMBL 
  alg <- factor( as.integer( inSelection[inUniverse] ) )
  names(alg) <- geneIDs[inUniverse]
  
  # The tests use a 0.01 \(p\)–value cutoff by default. We order by the classic algorithm.
  
  tab = as.list(onts)
  names(tab) = onts
  for(i in 1:3){
    
    ## prepare data
    tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=5,
                annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl" )
    
    ## run tests
    resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
    resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
    
    ## look at results
    tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                          Fisher.classic = resultTopGO.classic,
                          orderBy = "Fisher.classic" , topNodes = 200)
    
  }
  
  # We can now look at the results, we look at the top 200 GO categories according to the 
  # “Fisher classic” algorithm. The function “GenTable” produces a table of significant GO 
  # categories. Finally, we bind all resulting tables together and write them to an csv–file
  # that can be view with spreadsheet programs.
  
  topGOResults <- rbind.fill(tab)
  write.csv(topGOResults, file = "topGOResults.csv")
  
}

run_deseq(
  snakemake@input[["HTSEQ_COUNT_TABLE"]],
  snakemake@output[["DESEQ_RESULT"]]
)