#sccomp_function.R

library(Seurat)
library(sccomp)
#library(loo)
library(dplyr)
library(forcats)
library(tidyr)


##sample/condition/celltype
run_da_sccomp <- function(Seurat_obj, sample, celltype, contrasts) {
  
  condition <- unique(Seurat_obj@meta.data[["condition"]])

  contrasts <- sapply(condition, function(condition1) {
     sapply(condition, function(condition2) {
         if (condition1 != condition2) {
             paste("condition", condition1, " - condition", condition2, sep = "")
         } else {
             NULL  # Skip same condition comparisons
         }
     })
 })
 
contrasts <- unlist(contrasts)

  # Step 1: Estimate first model
  sccomp_result <- Seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ condition,
      .sample = sample,
      .cell_group = celltype,
      variational_inference = FALSE ,
      bimodal_mean_variability_association = TRUE,
      cores = 1
    ) |>
    sccomp_remove_outliers(cores = 1) |>
    sccomp_test()
  
  # Step 2: Estimate second model with contrasts
  sccomp_comparison <- Seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ 0 + condition,
      .sample = sample,
      .cell_group = celltype,
      variational_inference = FALSE ,
      bimodal_mean_variability_association = TRUE,
      cores = 1,
      verbose = FALSE
    ) |>
    sccomp_test(contrasts = c(contrasts))
  
  # Step 3: Estimate third model
  res <- Seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ condition,
      formula_variability = ~ condition,
      .sample = sample,
      .cell_group = celltype,
      variational_inference = FALSE ,
      bimodal_mean_variability_association = TRUE,
      cores = 1,
      verbose = FALSE
    )
  
  # Step 4: Compare models using leave-one-out cross-validation
  model_with_factor_association <- Seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ condition,
      .sample = sample,
      .cell_group = celltype,
      variational_inference = FALSE ,
      bimodal_mean_variability_association = TRUE,
      cores = 1,
      enable_loo = TRUE
    )
  
  model_without_association <- Seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ 1,
      .sample = sample,
      .cell_group = celltype,
      variational_inference = FALSE ,
      bimodal_mean_variability_association = TRUE,
      cores = 1,
      enable_loo = TRUE
    )
  
  #loo_comparison <- loo_compare(
 #   model_with_factor_association |> attr("fit") |> loo(),
  #  model_without_association |> attr("fit") |> loo()
 #)
  
  # Return all results in a list or appropriate structure
  result <- list(
    sccomp_result = sccomp_result,
    sccomp_comparison = sccomp_comparison,
    res = res
    #,
    #loo_comparison = loo_comparison
  )
  
  return(result)
}




