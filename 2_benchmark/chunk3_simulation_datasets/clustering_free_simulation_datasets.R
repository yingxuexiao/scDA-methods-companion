### clustering_free_simulation_datasets_processed.R
library(SingleCellExperiment)
library(scran)
library(tidyverse)
library(DAseq)
library(miloR)
library(tibble)
library(dplyr)
library(igraph)
library(cydar)
library(pdist)
library(reshape2)

# Load datasets
sce_branch <- readRDS('path_to_your_data/branch_data_bm.RDS') 
sce_cluster <- readRDS('path_to_your_data/cluster_data_bm.RDS') 
sce_balanced <- readRDS('path_to_your_data/cluster_balanced_data_bm.RDS') 
sce_linear <- readRDS('path_to_your_data/linear_data_bm.RDS') 

# Set parameters
pop <- "M1"  # Based on Yihai Dong's design
pop_enr <- 0.75  # Set the target cell type enrichment level in DA population
pop_col <- "celltype"  # Column name for cell type in the data
reduced.dim <- "PCA"  # Use PCA for dimensionality reduction
seed <- 123  # Set random seed

# Generate synthetic labels
sce_branch <- add_synthetic_labels_pop(branch_data_bm, pop = pop, pop_column = pop_col, pop_enr = pop_enr, redDim = reduced.dim, seed = seed)
sce_cluster <- add_synthetic_labels_pop(cluster_data_bm, pop = pop, pop_column = pop_col, pop_enr = pop_enr, redDim = reduced.dim, seed = seed)
sce_balanced <- add_synthetic_labels_pop(cluster_balanced_data_bm, pop = pop, pop_column = pop_col, pop_enr = pop_enr, redDim = reduced.dim, seed = seed)
sce_linear <- add_synthetic_labels_pop(linear_data_bm, pop = pop, pop_column = pop_col, pop_enr = pop_enr, redDim = reduced.dim, seed = seed)

# Set DA thresholds
if (pop_enr < 0.5) {
  da_lower <- pop_enr + (pop_enr / 100) * 10
  da_upper <- 1 - da_lower
} else {
  da_upper <- pop_enr - (pop_enr / 100) * 10
  da_lower <- 1 - da_upper
}

# Ensure da_upper is greater than da_lower
if (da_upper <= da_lower) {
  stop("Error: da_upper should be greater than da_lower. Please check the 'pop_enr' value.")
}

# Assign true labels based on Condition2_prob
true_labels <- ifelse(branch_data_bm$Condition2_prob < da_lower, "NegLFC", 
                        ifelse(branch_data_bm$Condition2_prob > da_upper, "PosLFC", "NotDA"))
colData(branch_data_bm)[["true_labels"]] <- true_labels

true_labels_2 <- ifelse(sce_cluster$Condition2_prob < da_lower, "NegLFC", 
                        ifelse(sce_cluster$Condition2_prob > da_upper, "PosLFC", "NotDA"))
colData(sce_cluster)[["true_labels"]] <- true_labels_2

true_labels_3 <- ifelse(sce_balanced$Condition2_prob < da_lower, "NegLFC", 
                        ifelse(sce_balanced$Condition2_prob > da_upper, "PosLFC", "NotDA"))
colData(sce_balanced)[["true_labels"]] <- true_labels_3

true_labels_4 <- ifelse(sce_linear$Condition2_prob < da_lower, "NegLFC", 
                        ifelse(sce_linear$Condition2_prob > da_upper, "PosLFC", "NotDA"))
colData(sce_linear)[["true_labels"]] <- true_labels_4

# Save processed datasets
saveRDS(sce_branch, 'sce_branch_with_DA.rds')
saveRDS(sce_cluster, 'sce_cluster_with_DA.rds')
saveRDS(sce_balanced, 'sce_balanced_with_DA.rds')
saveRDS(sce_linear, 'sce_linear_with_DA.rds')

# Convert expression matrix to dgCMatrix (sparse matrix)
library(Matrix)

# Convert SCE object expression matrix
expr_matrix <- as(as.matrix(counts(sce_cluster)), "dgCMatrix")

# Create Seurat object
seurat_cluster <- CreateSeuratObject(counts = expr_matrix)

# Extract PCA and UMAP data from SCE object
pca_data <- reducedDim(sce_cluster, "PCA")
umap_data <- reducedDim(sce_cluster, "UMAP")

# Add PCA and UMAP to Seurat object
seurat_cluster[["pca"]] <- CreateDimReducObject(embeddings = pca_data, key = "PC_", assay = "RNA")
seurat_cluster[["umap"]] <- CreateDimReducObject(embeddings = umap_data, key = "UMAP_", assay = "RNA")

# View Seurat object
str(seurat_cluster)

# Extract meta.data from SCE object
meta_data <- as.data.frame(colData(sce_cluster))

# Add meta.data to Seurat object
seurat_cluster@meta.data <- meta_data

# View Seurat object
str(seurat_cluster)

#### Change format for meld;cna_python:
# Extract expression matrix
coldata <- data.frame(colData(sce_linear)) %>% rownames_to_column()
write_csv(coldata,  "/Users/xiaoying/Desktop/DA结果分析/simulation_datasets/linear.coldata.csv")

# Extract PCA or other dimensionality reduction information
X_pca <- reducedDim(sce_linear, reduced.dim)

## Save reduced dims
write_csv(as.data.frame(X_pca) %>% rownames_to_column(), "/Users/xiaoying/Desktop/DA结果分析/simulation_datasets/linear.pca.csv")

## After this, use Python to read data and generate corresponding h5ad dataset:
import pandas as pd
import anndata
import numpy as np

X_pca = pd.read_csv('../simulation_datasets/linear.pca.csv', index_col=0)
obs = pd.read_csv('../simulation_datasets/linear.coldata.csv', index_col=0)
assert (X_pca.index == obs.index).all(), \
        "the index of data and meta data is different from each other"
# create the anndata
adata = anndata.AnnData(X_pca, obs=obs, dtype=np.float64)
adata.obs.index.name = "cell"

### All datasets have been saved to the simulation_datasets folder
