library(Seurat)
library(rcna)
library(ggthemes)
library(ggplot2)
library(dplyr)
library(glue)
library(harmony)
library(patchwork)
library(purrr)
library(Matrix)
library(data.table)
library(tidyr)
library(rlang)

run_da_cna1 <- function(seurat_obj, condition_col = "condition", celltype_col = "celltype") {

  # ---------------------------
  # ---------------------------
  nn_graphs <- grep("_nn$", names(seurat_obj@graphs), value = TRUE)
  if (length(nn_graphs) == 0) stop("No *_nn graph found in seurat_obj@graphs")
  
  if ("integrated_nn" %in% nn_graphs) {
    chosen_graph <- "integrated_nn"
  } else if ("RNA_nn" %in% nn_graphs) {
    chosen_graph <- "RNA_nn"
  } else {
    chosen_graph <- nn_graphs[1]
    warning(paste("Multiple *_nn graphs found, using", chosen_graph))
  }
  
  message("Using graph: ", chosen_graph)
  
  # ---------------------------
  # Automatically determine the number of condition groups
  # ---------------------------
  cond_values <- as.character(unique(seurat_obj@meta.data[[condition_col]]))
  if (length(cond_values) < 2) stop("Need at least two groups in `condition_col` for comparison.")
  
  # ---------------------------
  # Analysis of multiple pairs of combinations
  # ---------------------------
  result_list <- list()
  condition_combinations <- combn(cond_values, 2, simplify = FALSE)
  
  for (combo in condition_combinations) {
    message("Running CNA for: ", combo[1], " vs ", combo[2])
    
    cells_to_keep <- rownames(seurat_obj@meta.data)[
      as.character(seurat_obj@meta.data[[condition_col]]) %in% as.character(combo)
    ]
    sub_seurat <- subset(seurat_obj, cells = cells_to_keep)
    
    sub_seurat@meta.data$condition_numeric <- as.numeric(factor(as.character(sub_seurat@meta.data[[condition_col]])))
    
    sub_seurat <- association.Seurat(
      seurat_obj = sub_seurat, 
      test_var = 'condition_numeric', 
      samplem_key = 'orig.ident', 
      graph_use = chosen_graph, 
      verbose = TRUE,
      batches = NULL,
      covs = NULL
    )
    
    result_list[[paste(combo, collapse = "_vs_")]] <- sub_seurat@meta.data
  }
  
  # ---------------------------
  # Summary
  # ---------------------------
  summary_list <- list()
  
  for (comparison_name in names(result_list)) {
    message("Summarizing: ", comparison_name)
    df <- result_list[[comparison_name]]
    
    conds <- as.character(unique(df[[condition_col]]))
    if (length(conds) < 2) {
      warning("Skipping ", comparison_name, " due to insufficient condition count.")
      next
    }
    
    cond_pairs <- combn(conds, 2, simplify = FALSE)
    
    for (pair in cond_pairs) {
      cond1 <- pair[1]
      cond2 <- pair[2]
      
      df_subset <- df %>% filter(as.character(.data[[condition_col]]) %in% c(cond1, cond2))
      
      summary_tmp <- df_subset %>%
        group_by(!!sym(celltype_col), !!sym(condition_col)) %>%
        summarise(
          n_cells = n(),
          mean_cna = mean(cna_ncorrs, na.rm = TRUE),
          median_cna = median(cna_ncorrs, na.rm = TRUE),
          prop_signif_fdr05 = mean(cna_ncorrs_fdr05 > 0, na.rm = TRUE),
          prop_signif_fdr10 = mean(cna_ncorrs_fdr10 > 0, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        pivot_wider(
          names_from = !!sym(condition_col),
          values_from = c(n_cells, mean_cna, median_cna, prop_signif_fdr05, prop_signif_fdr10),
          names_sep = "."
        ) %>%
        mutate(
          mean_cna_diff = .data[[paste0("mean_cna.", cond2)]] - .data[[paste0("mean_cna.", cond1)]],
          median_cna_diff = .data[[paste0("median_cna.", cond2)]] - .data[[paste0("median_cna.", cond1)]],
          prop_signif_fdr05_diff = .data[[paste0("prop_signif_fdr05.", cond2)]] - .data[[paste0("prop_signif_fdr05.", cond1)]],
          prop_signif_fdr10_diff = .data[[paste0("prop_signif_fdr10.", cond2)]] - .data[[paste0("prop_signif_fdr10.", cond1)]],
          comparison = paste(comparison_name, cond1, "vs", cond2, sep = "_")
        ) %>%
        select(
          !!sym(celltype_col),
          all_of(paste0("n_cells.", cond1)),
          all_of(paste0("n_cells.", cond2)),
          mean_cna_diff,
          median_cna_diff,
          prop_signif_fdr05_diff,
          prop_signif_fdr10_diff,
          comparison
        )
      
      summary_list[[paste(comparison_name, cond1, "vs", cond2, sep = "_")]] <- summary_tmp
    }
  }
  
  final_summary_df <- bind_rows(summary_list)
  return(final_summary_df)
}
