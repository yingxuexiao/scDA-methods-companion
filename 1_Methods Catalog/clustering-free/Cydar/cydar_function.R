##Cydar：

  library(SingleCellExperiment)
  library(DAseq)
  library(miloR)
  library(tibble)
  library(tidyr)
  library(dplyr)
  library(igraph)
  library(cydar)
  library(pdist)
  library(reshape2)



run_da_cydar <- function(sce, condition_col,
                      sample_col,
                      reduced.dim,
                      d=20,
                      batch_col=NULL,
                      alpha=0.1,
                      tol=1.0,
                      downsample=10,
                      returnCd=TRUE){

  ## Make design matrix
  design_df <- as_tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  

  ## Make list for each sample

  sample_ls <- split(1:ncol(sce), sce[[sample_col]])

  processed.exprs <- lapply(sample_ls, function(s) reducedDim(sce[,s], reduced.dim)[,1:d])

  cd <- prepareCellData(processed.exprs)

  ## Count cells in hyperspheres
  cd <- cydar::countCells(cd, tol=tol, filter=1, downsample=downsample)

  # do DA testing with edgeR
  cd.dge <- DGEList(assay(cd), lib.size=cd$totals)
  
  sim.design <- model.matrix(design, data=design_df)[colnames(cd.dge),]

  sim.dge <- estimateDisp(cd.dge, sim.design)

  sim.fit <- glmQLFit(sim.dge, sim.design)

  sim.res <- glmQLFTest(sim.fit, coef=2)
  
  # control the spatial FDR
  cydar.res <- sim.res$table

  cydar.res$SpatialFDR <- spatialFDR(intensities(cd), sim.res$table$PValue)

  is.sig <- cydar.res$SpatialFDR <= alpha
  if (returnCd) {
    list(Cd=cd, DAres=cydar.res)
  } else {
    cydar.res
  }
}




