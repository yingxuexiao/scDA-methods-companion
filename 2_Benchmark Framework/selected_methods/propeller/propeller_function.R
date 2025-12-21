library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(limma)
library(ggplot2)
library(scater)
library(patchwork)
library(edgeR)
library(statmod)


run_da_propeller1 <- function(obj) {
  library(dplyr)
  
  # 1. Run propeller method, assuming the input is a Seurat object
  propeller_result <- propeller(obj)
  
  # 2. Extract the cluster and celltype mapping table
  # If it's a Seurat object, use obj@meta.data
  cluster_celltype_map <- obj@meta.data %>%
    dplyr::select(seurat_clusters, celltype) %>%
    dplyr::distinct()
  
  cluster_celltype_map$seurat_clusters <- as.character(cluster_celltype_map$seurat_clusters)
  
  # 3. Add the cluster column (row names) to the propeller result
  propeller_result$cluster <- rownames(propeller_result)
  
  # 4. Merge DA results with celltype mapping
  merged_result <- dplyr::left_join(propeller_result, cluster_celltype_map,
                                    by = c("cluster" = "seurat_clusters"))
  
  # 5. Reorder columns to put celltype first
  merged_result <- merged_result %>%
    dplyr::select(celltype, cluster, everything())
  
  # 6. Return two results: the raw result and the merged result
  return(list(
    raw_result = propeller_result,
    merged_result = merged_result
  ))
}
