library(edgeR)

run_edger <- function(COUNT_MATRIX) {

  # Filtering out non-expressed genes: consider only the genes with an average read count of 10 or more.
  
  means <- rowMeans(geneLevelCounts)
  filter <- means >= 10
  table(filter)
  geneLevelCounts <- geneLevelCounts[filter,]
  dim(geneLevelCounts)
  
  # Visualization of mapped reads/sequencing depth
  
  colors <- c(rep("yellow3",8),rep("blue",3),rep("green",3))
  totCounts <- colSums(geneLevelCounts)
  barplot(totCounts, las=2, col=colors)
  barplot(totCounts, las=2, col=c("red","blue")[laneInfo[,4]])
  
  # The boxplot function provides an easy way to visualize the difference in distribution between each experiment
  
  boxplot(geneLevelCounts, las=2, col=colors)
  boxplot(log2(geneLevelCounts+1), las=2, col=colors)
  
  # NORMALIZATION
  
  # Building the edgeR object
  # DGEList is the function that coverts the count matrix into an edgeR object.
  # In addition to the counts, we need to group the samples according to the variable of interest in our experiment.
  # We can then see the elements that the object contains by using the names function
  
  group <- laneInfo[,2]
  group <- droplevels(group)
  counts <- geneLevelCounts
  cds <- DGEList( counts , group = group )
  names(cds)
  
  # These elements can be accessed using the $ symbol. We give some examples here of looking at what is saved in this object.
  
  head(cds$counts) # original count matrix
  cds$samples # contains a summary of your samples
  sum(cds$all.zeros) # How many genes have 0 counts across all samples
  
  # We can calculate the normalization factors which correct for the different sequencing depth of each library (or library size).
  # By default, the function calcNormFactors normalize the data using the "weighted trimmed mean of M-values" (TMM) method. Other options are RLE and upper-quartile.
  
  cds <- calcNormFactors(cds)
  cds$samples
  
  # If we want the normalized pseudo-counts, useful for instance for cluster analysis, we can get them with the following commands.
  
  scale <- cds$samples$lib.size*cds$samples$norm.factors
  normCounts <- round(t(t(counts)/scale)*mean(scale))
  boxplot(log2(normCounts+1), las=2, col=colors)
  
  # Multi-Dimensional Scaling Plot
  # An MDS plot measures the similarity of the samples and projects this measure into 2-dimensions. This can be useful for quality control and visualization of your samples.
  
  plotMDS(cds, main="MDS Plot for Count Data", labels=colnames(cds$counts))
  
  # DIFFERENTIAL EXPRESSION TESTING
  # A two-class comparison
  
  group <- laneInfo[9:14,2]
  group <- droplevels(group)
  counts <- geneLevelCounts[, 9:14]
  cds <- DGEList( counts , group = group )
  cds <- calcNormFactors(cds, method="upperquartile")
  
  # Estimating Dispersion
  # First, we need to estimate the common dispersion: this assumes that all the genes have the same dispersion.
  
  cds <- estimateCommonDisp(cds, verbose=TRUE)
  
  # To understand what this value means, recall the parameterization for the variance of the negative binomial is ν(μ) = μ+μ2 ·φ. For poisson it’s ν(μ) = μ. The implied standard deviations are the square-roots of the variances. Now, suppose a gene had an average count of 200. Then the sd’s under the two models would be
  
  sqrt(200) # poisson sd
  sqrt(200 + 200^2 * cds$common.dispersion) # negative binomial
  sqrt(200 + 200^2 * cds$common.dispersion) / sqrt(200) # ratio: NB sd is almost 5 times larger
  
  # In real experiments, the assumption of common dispersion is almost never met.
  # Often, we observe a relation between mean counts and dispersion, i.e., the more
  # expressed genes have less dispersion.
  
  # The way edgeR estimates a tagwise (i.e. gene-wise) dispersion parameter is by "shrinking"
  # the gene-wise dispersions toward a common value (the common dispersion estimated in the previous step).
  # Alternatively, one can shrink the gene-wise estimates to a common trend, by estimating a smooth
  # function prior to the shrinkage (using the estimateTrendedDisp function). Here we keep things simple
  # and shrink the estimates to the common value.
  
  cds <- estimateTagwiseDisp(cds)
  
  # We can plot the biological coefficient of variation (i.e. the square root of the dispersions) with the following
  
  plotBCV(cds)
  
  # Another useful plot is the mean-variance relation.
  
  meanVarPlot <- plotMeanVar(cds, show.raw.vars=TRUE, show.tagwise.vars=TRUE, show.binned.common.disp.vars=FALSE, show.ave.raw.vars=FALSE, NBline = TRUE , nbins = 100 , #these are arguments about what is plotted
                             pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = "Mean-Variance Plot" ) #these arguments are to make it look prettier
  
  # Testing for DE
  # The function exactTest performs pair-wise tests for differential expression between two groups.
  # The important parameter is pair which indicates which two groups should be compared. The output
  # of exactTest is a list of elements: we can get the table of the results with the topTags function.
  
  et <- exactTest(cds, pair=levels(group))
  topTags(et)
  top <- topTags(et, n=nrow(cds$counts))$table
  head(top)
  
  # We can store the ID of the DE genes and look at the distribution of the p-values
  
  de <- rownames(top[top$FDR<0.05,])
  length(de)
  head(de)
  hist(top$PValue, breaks=20)
  
  # Visualizing the results
  # We can use the function plotSmear to produce a mean-difference plot
  
  plotSmear(cds , de.tags=de)
  abline(h=c(-2, 2), col="green")
  
  # Another useful plot is the "volcano plot" that relates log-fold-changes and p-values
  
  plot(top$logFC, -log10(top$PValue), pch=20, cex=.5, ylab="-log10(p-value)", xlab="logFC", col=as.numeric(rownames(top) %in% de)+1)
  abline(v=c(-2, 2), col="green")
  
  # Finally, we can plot a heatmap of the top genes and extract lists of up- and down-regulated genes.
  
  heatmap(log(normCounts[de[1:500],9:14]+1), ColSideColor=colors[9:14], col=rev(heat.colors(20)))
  upreg <- rownames(top[rownames(top) %in% de & top$logFC > 0,])
  downreg <- rownames(top[rownames(top) %in% de & top$logFC < 0,])
  
  # Outputting the results
  # We can write the results to a file that we can later open in excel or in other programs
  
  write.table(top, file="two-class-results.txt", sep='\t', quote=FALSE)
  
}

run_edger(
  snakemake@input[["COUNT_MATRIX"]],
  snakemake@output[["EDGER_RESULT"]]
)