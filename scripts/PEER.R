library(peer)
library(qtl)
library(impute)

run_peer <- function(GT_MATRIX, EXP_MATRIX, COV_MATRIX) {

  # Reding data
  cross < - read.cross(format="csvs",genfile=GT_MATRIX, phefile=EXP_MATRIX, genotypes=c(0,1))
  
  # Building the model
  model=PEER()
  
  # Set the maximum number of unobserved factors to model.
  PEER_setNk(model, n_unobserved_factors)
  
}

run_peer(
  snakemake@input[["GT_MATRIX"],["EXP_MATRIX"],["COV_MATRIX"]],
  snakemake@output[["xyz"]]
  )