# Load necessary libraries
library(Seurat)
library(SeuratData)
library(dplyr)

# 1. Load the ifnb dataset
ifnb <- LoadData("ifnb")
cat("Original dataset dimensions:", dim(ifnb), "\n")

# 2. Define the gene count gradients to generate
gene_gradients <- c(500, 1000, 2000, 4000, 6000, 8000, 10000, 12000, 13000, 14000)

# 3. Create an empty list to store each gene gradient dataset
gene_gradient_datasets <- list()

# 4. Randomly select genes based on the gene gradients and create datasets
for (n_genes in gene_gradients) {
  set.seed(123)  # Set seed for reproducibility
  
  # Ensure not to exceed total gene count
  n_genes_actual <- min(n_genes, nrow(ifnb))
  cat("Processing", n_genes_actual, "genes...\n")
  
  # Randomly select n_genes genes
  selected_genes <- sample(rownames(ifnb), n_genes_actual, replace = FALSE)
  
  # Create the subset dataset
  gene_subset <- ifnb[selected_genes, ]
  
  # 5. Add sample information
  if ("stim" %in% colnames(gene_subset@meta.data)) {
    cell_metadata <- gene_subset@meta.data
    
    # Get cells for each condition
    stim_cells <- which(cell_metadata$stim == "STIM")
    ctrl_cells <- which(cell_metadata$stim == "CTRL")
    
    # Create sample names for each condition
    stim_samples <- paste0("STIM_rep", 1:3)  # 3 STIM samples
    ctrl_samples <- paste0("CTRL_rep", 1:3)  # 3 CTRL samples
    
    # Assign sample information
    cell_metadata$sample <- NA
    cell_metadata$sample[stim_cells] <- sample(stim_samples, length(stim_cells), replace = TRUE)
    cell_metadata$sample[ctrl_cells] <- sample(ctrl_samples, length(ctrl_cells), replace = TRUE)
    
    # Add sample information to gene_subset
    gene_subset$sample <- cell_metadata$sample
    
    cat("Sample distribution for", n_genes_actual, "genes:\n")
    print(table(gene_subset$sample))
    
  } else {
    stop("The 'stim' column is missing in the metadata.")
  }
  
  # 6. Save each dataset to the list
  dataset_name <- paste0("gene_", n_genes_actual)
  gene_gradient_datasets[[dataset_name]] <- gene_subset
}

# 7. Check if datasets were successfully generated
cat("Number of datasets generated:", length(gene_gradient_datasets), "\n")

# 8. Perform Seurat standard preprocessing for each dataset
clustering_results <- list()

for (dataset_name in names(gene_gradient_datasets)) {
  cat("\nProcessing", dataset_name, "\n")
  
  dataset <- gene_gradient_datasets[[dataset_name]]
  
  # Generate Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = dataset@assays$RNA@layers$counts, 
    meta.data = dataset@meta.data
  )
  
  cat("  Cells:", ncol(seurat_obj), "Genes:", nrow(seurat_obj), "\n")
  
  # Standard preprocessing pipeline
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  
  # Dynamically set nfeatures parameter
  n_features <- min(2000, nrow(seurat_obj))
  cat("  Using", n_features, "variable features\n")
  
  seurat_obj <- FindVariableFeatures(seurat_obj, 
                                     selection.method = "vst", 
                                     nfeatures = n_features,
                                     verbose = FALSE)
  
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  
  # Run PCA
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  
  # Dynamically set PCA dimensions
  pca_dims <- 1:min(20, length(seurat_obj[["pca"]]))
  cat("  Using PCA dimensions:", min(pca_dims), "to", max(pca_dims), "\n")
  
  # Clustering analysis
  seurat_obj <- FindNeighbors(seurat_obj, dims = pca_dims, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = pca_dims, verbose = FALSE)
  
  # Save clustering results
  clustering_results[[dataset_name]] <- list(
    seurat_obj = seurat_obj,
    n_genes = nrow(seurat_obj),
    n_cells = ncol(seurat_obj),
    samples = unique(seurat_obj$sample)
  )
}

# 9. Extract 10 Seurat objects
seurat_objects <- list()

for (dataset_name in names(clustering_results)) {
  seurat_obj <- clustering_results[[dataset_name]]$seurat_obj
  seurat_objects[[dataset_name]] <- seurat_obj
}

# Check the number and information of extracted objects
cat("\n=== Final Results ===\n")
cat("Number of Seurat objects:", length(seurat_objects), "\n\n")

# Print basic information for each object
for (name in names(seurat_objects)) {
  obj <- seurat_objects[[name]]
  cat("Dataset:", name, "\n")
  cat("  Cells:", ncol(obj), "\n")
  cat("  Genes:", nrow(obj), "\n")
  cat("  Samples:", length(unique(obj$sample)), "\n")
  cat("  Clusters:", length(unique(Idents(obj))), "\n\n")
}

gene_500 <- seurat_objects$gene_500
gene_1000 <- seurat_objects$gene_1000
gene_2000 <- seurat_objects$gene_2000
gene_4000 <- seurat_objects$gene_4000
gene_6000 <- seurat_objects$gene_6000
gene_8000 <- seurat_objects$gene_8000
gene_10000 <- seurat_objects$gene_10000
gene_12000 <- seurat_objects$gene_12000
gene_13000 <- seurat_objects$gene_13000
gene_14000 <- seurat_objects$gene_14000

library(Seurat)
library(SingleCellExperiment)
library(sceasy)
library(reticulate)

# Set conda environment
use_condaenv("base", required = TRUE)

# Define save path
save_path <- "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/gene_number_simulated_datasets/"

# Process each Seurat object and save
for (dataset_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[dataset_name]]  # Get Seurat object
  
  # Add condition column (copy stim column)
  seurat_obj$condition <- seurat_obj$stim
  
  # Add celltype column (copy seurat_annotations column)
  seurat_obj$celltype <- seurat_obj$seurat_annotations
  
  # Save as .rds file
  rds_file <- paste0(save_path, dataset_name, "_seurat.rds")
  saveRDS(seurat_obj, file = rds_file)
  
  seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")
  
  # Convert to .h5ad format
  h5ad_file <- paste0(save_path, dataset_name, ".h5ad")
  sceasy::convertFormat(
    seurat_obj, from = "seurat", to = "anndata", outFile = h5ad_file
  )
  
  # Convert to SingleCellExperiment object
  seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")  # Ensure Assay type is correct
  sce_obj <- as.SingleCellExperiment(seurat_obj)
  
  # Convert PCA and UMAP data
  if ("pca" %in% names(seurat_obj@reductions)) {
    reducedDims(sce_obj)[["PCA"]] <- Embeddings(seurat_obj, "pca")
  }
  if ("umap" %in% names(seurat_obj@reductions)) {
    reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seurat_obj, "umap")
  }
  
  # Save as .rds file (SingleCellExperiment object)
  sce_rds_file <- paste0(save_path, dataset_name, "_sce.rds")
  saveRDS(sce_obj, file = sce_rds_file)
}

cat("All datasets processed and saved successfully.\n")

##########################################################################################
##########################################################################################
# Extract different numbers of cells from the ifnb dataset to simulate cell number variations more efficiently

library(Seurat)
library(SeuratData)

# 1. Load the ifnb dataset
ifnb <- LoadData("ifnb")

# 2. Define the cell number gradients
cell_gradients <- c(500, 1000, 2000, 4000, 6000, 8000, 10000, 12000, 13000, 13999)

# 3. Create an empty list to store each cell gradient dataset
cell_gradient_datasets <- list()

# 4. Extract cells from the ifnb dataset according to cell number gradients
for (n_cells in cell_gradients) {
  set.seed(123)  # Set seed for reproducibility
  
  # Randomly select n_cells cells
  selected_cells <- sample(Cells(ifnb), n_cells, replace = FALSE)
  n_cell_actual <- min(n_cells, nrow(ifnb))

  # Create the subset dataset
  cell_subset <- subset(ifnb, cells = selected_cells)
  
  # 5. Add sample information
  if ("stim" %in% colnames(cell_subset@meta.data)) {
    cell_metadata <- cell_subset@meta.data
    
    # Get cells for each condition
    stim_cells <- which(cell_metadata$stim == "STIM")
    ctrl_cells <- which(cell_metadata$stim == "CTRL")
    
    # Create sample names for each condition
    stim_samples <- paste0("STIM_rep", 1:3)  # 3 STIM samples
    ctrl_samples <- paste0("CTRL_rep", 1:3)  # 3 CTRL samples
    
    # Assign sample information
    cell_metadata$sample <- NA
    cell_metadata$sample[stim_cells] <- sample(stim_samples, length(stim_cells), replace = TRUE)
    cell_metadata$sample[ctrl_cells] <- sample(ctrl_samples, length(ctrl_cells), replace = TRUE)
    
    # Add sample information to cell_subset
    cell_subset$sample <- cell_metadata$sample
    
    cat("Sample distribution for", n_cell_actual, "genes:\n")
    print(table(cell_subset$sample))
    
  } else {
    stop("The 'stim' column is missing in the metadata.")
  }
  
  # 6. Save each dataset to the list
  dataset_name <- paste0("cell_", n_cell_actual)
  cell_gradient_datasets[[dataset_name]] <- cell_subset
}

# 7. Check if datasets were successfully generated
cat("Number of datasets generated:", length(cell_gradient_datasets), "\n")

# 6. Perform Seurat standard preprocessing for each dataset
clustering_results <- list()

for (dataset_name in names(cell_gradient_datasets)) {
  cat("Processing", dataset_name, "\n")
  
  dataset <- cell_gradient_datasets[[dataset_name]]
  
  # Generate Seurat object
  seurat_obj <- CreateSeuratObject(counts = dataset@assays$RNA@layers$counts, meta.data = dataset@meta.data)
  
  # Standard preprocessing pipeline
  seurat_obj <- NormalizeData(seurat_obj)  # Data normalization
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)  # Find variable genes
  seurat_obj <- ScaleData(seurat_obj)  # Data scaling
  seurat_obj <- RunPCA(seurat_obj)  # Principal component analysis
  
  # Clustering analysis
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)  # Calculate nearest neighbors based on first 20 PCs
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)  # Perform clustering based on neighbors
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)  # UMAP visualization

  # Save clustering results
  clustering_results[[dataset_name]] <- list(
    seurat_obj = seurat_obj,
    n_genes = nrow(seurat_obj),
    n_cells = ncol(seurat_obj),
    samples = unique(seurat_obj$sample)
  )
}

# 9. Extract Seurat objects
seurat_objects <- list()

for (dataset_name in names(clustering_results)) {
  seurat_obj <- clustering_results[[dataset_name]]$seurat_obj
  seurat_objects[[dataset_name]] <- seurat_obj
}

# Check the number of extracted objects
cat("\n=== Final Results ===\n")
cat("Number of Seurat objects:", length(seurat_objects), "\n\n")

# Print each object's basic information
for (name in names(seurat_objects)) {
  obj <- seurat_objects[[name]]
  cat("Dataset:", name, "\n")
  cat("  Cells:", ncol(obj), "\n")
  cat("  Genes:", nrow(obj), "\n")
  cat("  Samples:", length(unique(obj$sample)), "\n")
  cat("  Clusters:", length(unique(Idents(obj))), "\n\n")
}

library(Seurat)
library(SingleCellExperiment)
library(sceasy)
library(reticulate)

# Set conda environment
use_condaenv("base", required = TRUE)

# Define save path
save_path <- "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/cell_number_simulated_datasets/"

# Process each Seurat object and save
for (dataset_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[dataset_name]]  # Get Seurat object
  
  # Add condition column (copy stim column)
  seurat_obj$condition <- seurat_obj$stim
  
  # Add celltype column (copy seurat_annotations column)
  seurat_obj$celltype <- seurat_obj$seurat_annotations
  
  # Save as .rds file
  rds_file <- paste0(save_path, dataset_name, "_seurat.rds")
  saveRDS(seurat_obj, file = rds_file)
  
  seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")
  
  # Convert to .h5ad format
  h5ad_file <- paste0(save_path, dataset_name, ".h5ad")
  sceasy::convertFormat(
    seurat_obj, from = "seurat", to = "anndata", outFile = h5ad_file
  )
  
  # Convert to SingleCellExperiment object
  seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")  # Ensure Assay type is correct
  sce_obj <- as.SingleCellExperiment(seurat_obj)
  
  # Convert PCA and UMAP data
  if ("pca" %in% names(seurat_obj@reductions)) {
    reducedDims(sce_obj)[["PCA"]] <- Embeddings(seurat_obj, "pca")
  }
  if ("umap" %in% names(seurat_obj@reductions)) {
    reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seurat_obj, "umap")
  }
  
  # Save as .rds file (SingleCellExperiment object)
  sce_rds_file <- paste0(save_path, dataset_name, "_sce.rds")
  saveRDS(sce_obj, file = sce_rds_file)
}

cat("All datasets processed and saved successfully.\n")
