# DESeq PCA and differential expression analysis of gene count data obtained with HTSeq-counts

# Converting a HTSeq count table to a DESeqDataSet
design <- "~ condition" # condition tumor or normal
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable, directory = ".", design,
                                          ignoreRank = FALSE, ...)

# Quality control and normalization of data
GeneCounts <- counts(DESeq2Table)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)