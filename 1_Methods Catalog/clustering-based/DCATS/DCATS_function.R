library(DCATS)
library(tidyverse)
library(Seurat)
library(patchwork)
library(SingleCellExperiment)
library(scater)
library(scran) 
library(dplyr)
library(patchwork)
library(scuttle)
library(Matrix)

# Main function
run_da_dcats <- function(seurat, condition_col = "condition", sample_col = "sample", celltype_col = "celltype") {
   # Internal function da_: Automatically get the appropriate SNN graph structure
  get_snn_graph <- function(seurat_obj) {
    graph_names <- names(seurat_obj@graphs)
    snn_graphs <- graph_names[grepl("_snn$", graph_names)]

    if (length(snn_graphs) == 0) {
      stop("❌ No `_snn` graph structure found, please ensure the Seurat object contains clustering graphs.")
    }

    priority_graph <- c("integrated_snn", "RNA_snn")
    selected_graph <- intersect(priority_graph, snn_graphs)

    if (length(selected_graph) == 0) {
      selected_graph <- snn_graphs[1]
      warning(paste("⚠️ No integrated_snn / RNA_snn found, using default graph: ", selected_graph))
    } else {
      selected_graph <- selected_graph[1]
    }

    message(paste("✅ Using graph structure: ", selected_graph))
    return(seurat_obj@graphs[[selected_graph]])
  }

  # Step 1: Extract variables
  condition_vec <- seurat[[condition_col]][, 1]
  celltype_vec  <- seurat[[celltype_col]][, 1]

  # Step 2: Generate all condition pairs
  generate_condition_pairs <- function(conditions) {
    conditions <- unique(na.omit(as.character(conditions)))
    if (length(conditions) < 2) stop("❌ Not enough conditions, at least two different conditions are required.")
    combn(conditions, 2, simplify = FALSE)
  }
  condition_pairs <- generate_condition_pairs(condition_vec)

  # Construct similarity matrix (specific to DCATS)

  graph_mat <- get_snn_graph(seurat)

  knn_mat <- knn_simMat(graph_mat, celltype_vec)

  # Construct count matrix
  count_mat <- table(condition_vec, celltype_vec)

  results <- list()

  for (pair in condition_pairs) {
    cond1 <- pair[1]
    cond2 <- pair[2]
    
    subset_count_mat <- count_mat[c(cond1, cond2), , drop = FALSE]
    design_mat <- data.frame(condition = factor(c(cond1, cond2)))
    rownames(design_mat) <- c(cond1, cond2)

    da_result <- dcats_GLM(subset_count_mat, design_mat, similarity_mat = knn_mat)
    results[[paste(cond1, cond2, sep = "_vs_")]] <- da_result
  }

  return(results)
}
