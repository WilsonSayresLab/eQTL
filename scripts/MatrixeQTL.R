# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

run_matrixeqtl <- function(GT_MATRIX, EXP_MATRIX, COV_MATRIX, SNP_LOC, GENE_LOC, CIS_OUTPUT_NAME, TRANS_OUTPUT_NAME) {

  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Genotype file name
  SNP_file_name = GT_MATRIX
  snps_location_file_name = SNP_LOC
  
  # Gene expression file name
  expression_file_name = EXP_MATRIX
  gene_location_file_name = GENE_LOC
  
  # Covariates file name
  # Set to character() for no covariates
  covariates_file_name = COV_MATRIX
  
  # Output file names
  output_file_name_cis = CIS_OUTPUT_NAME
  output_file_name_tra = TRANS_OUTPUT_NAME
  
  # The p-value threshold determines which gene-SNP associations are saved in 
  # the output file output_file_name. Note that for larger datasets the threshold 
  # should be lower. Setting the threshold to a high value for a large dataset may 
  # cause excessively large output files.
  pvOutputThreshold_cis = 2e-2
  pvOutputThreshold_tra = 1e-2

  # Define the covariance matrix for the error term. This parameter is rarely used. 
  # If the covariance matrix is a multiple of identity, set it to numeric().
  errorCovariance = numeric()
  
  # Distance for local gene-SNP pairs
  cisDist = 1e6
  
  # The next section of the sample code contains three very similar parts loading 
  # the files with genotype, gene expression, and covariates. In each part one can 
  # set the file delimiter (i.e. tabulation "\t", comma ",", or space " "), the string 
  # representation for missing values, the number of rows with column labels, and the 
  # number of columns with row labels. Finally, one can change the number of the 
  # variables in a slice for the file reading procedure (do not change if not sure).
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  unlink(output_file_name_tra);
  unlink(output_file_name_cis);
  
  ## Results:
  
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
  cat('Detected local eQTLs:', '\n');
  show(me$cis$eqtls)
  cat('Detected distant eQTLs:', '\n');
  show(me$trans$eqtls)
  
  ## Plot the Q-Q plot of local and distant p-values
  
  plot(me)
  
}

run_matrixeqtl(
  snakemake@input[["GT_MATRIX"],["EXP_MATRIX"],["COV_MATRIX"],["SNP_LOC"],["GENE_LOC"],["CIS_OUTPUT_NAME"],["TRANS_OUTPUT_NAME"]],
  snakemake@output[["CIS_RESULT"],["TRANS_RESULT"]]
)