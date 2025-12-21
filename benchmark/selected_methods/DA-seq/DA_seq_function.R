#DA_seq_function:

library(stringr)
library(DAseq)
library(Seurat)
library(SingleCellExperiment)

run_daseq_analysis <- function(
  sce,
  k.vec = c(10, 20, 30, 40),
  condition_col = "condition",
  sample_col = "samples",
  celltype_col = "celltype",
  reduced_dim = "PCA",
  d = 30

) {
  # Get the combination of all condition groups
  condition_vec <- colData(sce)[["condition"]]
  all_conditions <- unique(condition_vec)
  combs <- combn(all_conditions, 2, simplify = FALSE)

  result.list <- list()  # Used to store the results of all combinations

  for (pair in combs) {
    cond1 <- pair[1]
    cond2 <- pair[2]
    message(paste0("🔍 Running DA-seq for ", cond1, " vs ", cond2))

    # Screen cells belonging to cond1 and cond2
    selected_cells <- which(condition_vec %in% c(cond1, cond2))
    sce_sub <- sce[, selected_cells]

    # Extract metadata
    condition_sub <- colData(sce_sub)[["condition"]]
    sample_sub <- colData(sce_sub)[["sample"]]
    
    # run DA-seq
    daseq_res <- getDAcells(
      X = reducedDim(sce_sub, reduced_dim)[, 1:d],
      cell.labels = sample_sub,
      labels.1 = unique(sample_sub[condition_sub == cond1]),
      labels.2 = unique(sample_sub[condition_sub == cond2]),
      k.vector = k.vec,
      size = 1,
      do.plot = FALSE
    )

    # add DA tags
    da_status <- rep("non-DA", ncol(sce_sub))
    da_status[daseq_res$da.up] <- "DA-up"
    da_status[daseq_res$da.down] <- "DA-down"
    colData(sce_sub)$da_status <- da_status

    # Table statistics
    tab1 <- table(colData(sce_sub)[[celltype_col]], colData(sce_sub)$da_status)
    tab2 <- prop.table(tab1, margin = 1)
    tab3 <- table(colData(sce_sub)[[condition_col]], colData(sce_sub)$da_status)

    # Visualization
    umap_df <- data.frame(
      UMAP1 = reducedDim(sce_sub, "UMAP")[, 1],
      UMAP2 = reducedDim(sce_sub, "UMAP")[, 2],
      da_status = colData(sce_sub)$da_status,
      cell_type = colData(sce_sub)[[celltype_col]]
    )

    library(ggplot2)
    umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = da_status)) +
      geom_point(size = 0.5) +
      theme_classic() +
      ggtitle(paste0("DA-seq: ", cond1, " vs ", cond2))

    # Save to the result list
    key <- paste0(cond1, "_vs_", cond2)
    result.list[[key]] <- list(
      sce = sce_sub,
      daseq_res = daseq_res,
      table_celltype_count = tab1,
      table_celltype_proportion = tab2,
      table_condition_count = tab3,
      umap_plot = umap_plot,
      group1 = cond1,
      group2 = cond2
    )
  }

 
  message("✅ All the DA-seq results of pair-to-two combinations have been consolidated into result.list!")
  
  return(result.list)
}






























