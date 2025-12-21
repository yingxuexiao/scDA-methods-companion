################## run milo   ----------------------------------------------------------------------------
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(scuttle)
library(tibble)
library(pdist)
library(MouseGastrulationData)


##input: sample/condition/celltype/ batch(if have)
##
run_da_milo <- function(
  sce,
  condition ,
  sample ,
  reduced.dim = NULL,
  k = NULL,
  d = NULL,
  prop = NULL,
  returnMilo = TRUE,
  batch = NULL
) {
  
  ## Build the Milo object
  milo <- Milo(sce)

  ## Build an adjacency graph
  milo <- buildGraph(milo, k = k, d = d, reduced.dim = reduced.dim)

  ## Build a neighborhood
  milo <- makeNhoods(milo, prop = prop, k = k, d = d, refined = TRUE, reduced_dim = reduced.dim)

  ## Count
  milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), samples = "sample")

  ## Build the design matrix (automatically adjusted based on whether batch is passed in)
  col_data <- as.data.frame(colData(milo))

  # # Determine whether batch needs to be added
  if (!is.null(batch) && batch %in% colnames(col_data)) {
    milo_design <- col_data[, c("sample", "condition", batch)]
    design_formula <- as.formula(paste("~", batch, "+", condition))
  } else {
    milo_design <- col_data[, c("sample", "condition")]
    design_formula <- as.formula(paste("~", condition))
  }

  milo_design <- distinct(milo_design)
  rownames(milo_design) <- milo_design[[sample]]

  ## Distance calculation
  milo <- calcNhoodDistance(milo, d = d, reduced.dim = reduced.dim)

  ## DA test
  da_results <- testNhoods(
    milo,
    design = design_formula,
    design.df = milo_design,
    reduced.dim = reduced.dim
  )

  ##Comment on the cell types of the neighborhood
  if ("celltype" %in% colnames(colData(milo))) {
    da_results <- annotateNhoods(milo, da_results, coldata_col = "celltype")
  }

  ## Return result
  if (returnMilo) {
    return(da_results)
  } else {
    return(da_results %>% arrange(-SpatialFDR) %>% head())
  }
}








