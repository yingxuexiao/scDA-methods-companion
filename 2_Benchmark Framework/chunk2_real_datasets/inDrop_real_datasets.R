##inDrop_dataset:

###############################################################################
############################-------     GSE102827   ---------##################

library(data.table)
library(Seurat)
library(dplyr)
mem.maxVSize(10*1024*3)

expr <- fread("/Users/xiaoying/Downloads/GSE102827_MATRIX.csv.gz") |> as.data.frame()
rownames(expr) <- expr[[1]]
expr[[1]] <- NULL

expr <- log1p(expr)

meta <- fread("/Users/xiaoying/Downloads/GSE102827_cell_type_assignments.csv.gz") |> as.data.frame()
rownames(meta) <- meta$V1   
colnames(meta) <- as.character(unlist(meta[1, ]))
meta <- meta[-1, ]
meta <- meta[,-1 ]

all(colnames(expr) == rownames(meta))

library(Seurat)

seurat_obj <- CreateSeuratObject(counts = expr, meta.data = meta)

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1)

seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 4000)

seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)

seurat_obj <- RunTSNE(seurat_obj, dims = 1:30, perplexity = 30)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
DimPlot(seurat_obj, reduction = "tsne", group.by = "celltype")

table(seurat_obj$sample,seurat_obj$maintype)

counts_mat <- GetAssayData(seurat_obj, slot = "counts")

seurat_obj$UMI_total <- colSums(counts_mat)

seurat_obj <- subset(seurat_obj, subset = UMI_total >= 700 & UMI_total <= 15000)

mito_genes <- grep("^mt-", rownames(seurat_obj), value = TRUE)
seurat_obj <- seurat_obj[!rownames(seurat_obj) %in% mito_genes, ]

cluster_sizes <- table(seurat_obj$seurat_clusters)
keep_clusters <- names(cluster_sizes[cluster_sizes >= 100])
seurat_obj <- subset(seurat_obj, idents = keep_clusters)

cells_to_keep <- WhichCells(seurat_obj, expression = !is.na(maintype))
seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

saveRDS(seurat_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/Indrop_datasets/GSE102827_seurat.rds")

seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")
library(sceasy)
library(reticulate)
sceasy::convertFormat(
  seurat_obj,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/Indrop_datasets/GSE102827.h5ad"
)

library(SingleCellExperiment)
seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seurat_obj)

if ("pca" %in% names(seurat_obj@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(seurat_obj, "pca")
}
if ("umap" %in% names(seurat_obj@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seurat_obj, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/Indrop_datasets/GSE102827_sce.rds")









###############################################################################
############################-------     GSE270745   ---------##################
library(Seurat)
library(Matrix)

data_dir <- "/Users/xiaoying/Downloads/GSE270745_RAW/"

samples <- list(
  M577 = "GSM8351344_8594-NM-1_S1_L005",
  M576 = "GSM8351345_8594-NM-2_S1_L005",
  M580 = "GSM8351346_8594-NM-4_S1_L005",
  M578 = "GSM8351347_8594-NM-7_S1_L005",
  M579 = "GSM8351348_8594-NM-8_S1_L005",
  M582 = "GSM8351349_8594-NM-10_S1_L005"
)

for(sample_name in names(samples)){
  prefix <- samples[[sample_name]]
  
  sample_dir <- file.path(data_dir, sample_name)
  if(!dir.exists(sample_dir)) dir.create(sample_dir)
  
  file.copy(
    from = list.files(data_dir, pattern = paste0("^", prefix, "_.*\\.mtx\\.gz$"), full.names = TRUE),
    to = file.path(sample_dir, "matrix.mtx.gz")
  )
  file.copy(
    from = list.files(data_dir, pattern = paste0("^", prefix, "_.*barcodes\\.tsv\\.gz$"), full.names = TRUE),
    to = file.path(sample_dir, "barcodes.tsv.gz")
  )
  file.copy(
    from = list.files(data_dir, pattern = paste0("^", prefix, "_.*features\\.tsv\\.gz$"), full.names = TRUE),
    to = file.path(sample_dir, "features.tsv.gz")
  )
}

library(Seurat)

samples <- c("M577","M576","M580","M578","M579","M582")

seurat_list <- lapply(samples, function(sample_name){
  sample_dir <- file.path(data_dir, sample_name)
  
  counts <- Read10X(data.dir = sample_dir, gene.column = 1)
  
  seurat_obj <- CreateSeuratObject(counts = counts, project = "GSE270745", min.cells = 3, min.features = 200)
  seurat_obj$sample <- sample_name
  
  return(seurat_obj)
})

seurat_combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = samples)
seurat_combined

sample_condition <- data.frame(
  sample = c("M577","M576","M580","M578","M579","M582"),
  condition = c("control","control","Cdiff","control","Cdiff","Cdiff"),
  stringsAsFactors = FALSE
)

head(seurat_combined$sample)

seurat_combined$condition <- sample_condition$condition[match(seurat_combined$sample, sample_condition$sample)]

table(seurat_combined$sample, seurat_combined$condition)

library(Seurat)
library(dplyr)
library(ggplot2)

seurat_combined[["percent.mt"]] <- PercentageFeatureSet(seurat_combined, pattern = "^MT-")
seurat_combined <- subset(seurat_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

seurat_combined <- NormalizeData(seurat_combined, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)

seurat_combined <- ScaleData(seurat_combined, vars.to.regress = c("nFeature_RNA", "percent.mt"))

seurat_combined <- RunPCA(seurat_combined, features = VariableFeatures(object = seurat_combined))

seurat_combined <- FindNeighbors(seurat_combined, dims = 1:20)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)

seurat_combined <- RunUMAP(seurat_combined, dims = 1:20)

DimPlot(seurat_combined, group.by = "sample")
DimPlot(seurat_combined, group.by = "seurat_clusters")

markers_list <- list(
  Goblet = c("Clca1", "Tff3"),
  Absorptive = c("Slc26a3", "Car1"),
  EEC = c("Chga", "Chgb"),
  Tuft = c("Dclk1"),
  DCS = c("Reg4"),
  Stem = c("Lgr5", "Hopx"),
  Progenitors = c("Lrig1", "Birc5")
)

DotPlot(seurat_combined, features = markers_list) +
  RotatedAxis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Idents(seurat_combined) <- "seurat_clusters"

map <- c(
  `0`  = "Goblet",
  `1`  = "Goblet",
  `2`  = "Progenitors",
  `3`  = "Absorptive",
  `4`  = "DCS",
  `5`  = "DCS",
  `6`  = "Absorptive",
  `7`  = "Stem",
  `8`  = "Absorptive",
  `9`  = "Goblet",
  `10` = "Goblet",
  `11` = "Tuft",
  `12` = "Goblet",
  `13` = "EEC"
)

all_clusters <- levels(Idents(seurat_combined))
missing <- setdiff(all_clusters, names(map))
if (length(missing) > 0) {
  add <- setNames(rep("Unknown", length(missing)), missing)
  map <- c(map, add)
}

seurat_combined <- RenameIdents(seurat_combined, map)

seurat_combined$celltype_manual <- Idents(seurat_combined)

table(seurat_combined$celltype_manual)
table(seurat_combined$condition, seurat_combined$celltype_manual)
round(prop.table(table(seurat_combined$condition, seurat_combined$celltype_manual), margin = 1), 3)

DimPlot(seurat_combined, group.by = "celltype_manual", label = TRUE, repel = TRUE)

seurat_combined$celltype <- seurat_combined$celltype_manual
seurat_combined$celltype_manual <- NULL

saveRDS(seurat_combined, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/Indrop_datasets/GSE270745_seurat.rds")

seurat_combined[["RNA"]] <- as(seurat_combined[["RNA"]], "Assay")

library(sceasy)
library(reticulate)
sceasy::convertFormat(
  seurat_combined,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/Indrop_datasets/GSE270745.h5ad"
)

library(SingleCellExperiment)
seurat_combined[["RNA"]] <- as(seurat_combined[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seurat_combined)

if ("pca" %in% names(seurat_combined@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(seurat_combined, "pca")
}
if ("umap" %in% names(seurat_combined@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seurat_combined, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/Indrop_datasets/GSE270745_sce.rds")






