# CTDS_function.R

library(tidyverse)
library(hablar)
library(Seurat)
library(SingleCellExperiment)
library(knitr)
library(kableExtra)
library(patchwork)

source(file.path(methods_dir, "CTDS", "R", "CTDS_function.R"))

run_da_CTDS <- function(
  obj,
  sample_col, 
  celltype_col ,
  condition_col ,
  batch_col = NULL  # Default does not use batch, will automatically check if it exists
) {
  # 1. Create cell proportion matrix
  seurat_prop <- prop.table(
    table(
      sample = obj@meta.data[["sample"]],
      cell.type = obj@meta.data[["celltype"]]
    ),
    margin = 1
  )

  # 2. Construct metadata (automatically handles if batch is included)
  metadata <- obj@meta.data %>%
    as_tibble() %>%
    dplyr::select(all_of(c(sample_col, condition_col, batch_col))) %>%
    dplyr::distinct()

  # 3. Process proportion matrix
  normprop <- seurat_prop %>%
    tibble::as_tibble() %>%
    dplyr::rename("proportion" = n) %>%
    dplyr::full_join(metadata, by = sample_col) %>%
    hablar::convert(fct(cell.type))

  # 4. CTDS analysis main step
  div.res <- CTDS.score(
    obj,
    sample = sample_col,
    cell.type = celltype_col,
    metadata = metadata
  )

  # 5. ANOVA analysis: Detect significant differences between groups
  # Construct dynamic formula
  formula_str <- paste0("statistic ~ ", condition_col)
  lm.model <- lm(as.formula(formula_str), data = div.res)
  anova_results <- anova(lm.model) %>% broom::tidy()

  # 6. Pairwise comparison (t-test)
  pairwise_results <- pairwise.t.test(
    div.res$statistic,
    div.res[[condition_col]],
    p.adjust.method = "BH"
  ) %>%
    broom::tidy()

  # 7. Return all results
  return(list(
    div_res = div.res,
    anova = anova_results,
    pairwise_results = pairwise_results
  ))
}
