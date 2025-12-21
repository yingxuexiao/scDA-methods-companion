###Smart_seq2_real_datasets:

################################################################################
#####################----------      GSE140228      --------####################
library(Seurat)
library(dplyr)
library(ggplot2)

cell_info_path <- "Downloads/GSE140228_cell_info_Smartseq2.tsv.gz"
cell_info <- read.table(gzfile(cell_info_path), header = TRUE, sep = "\t")

head(cell_info)

gene_info_path <- "Downloads/GSE140228_gene_info_Smartseq2.tsv.gz"
gene_info <- read.table(gzfile(gene_info_path), header = TRUE, sep = "\t")

head(gene_info)

read_counts_path <- "Downloads/GSE140228_read_counts_Smartseq2.csv.gz"
read_counts <- read.csv(gzfile(read_counts_path), row.names = 1)
head(read_counts)

library(Seurat)
cell_info$Barcode <- gsub("-", ".", cell_info$Barcode)
seurat_obj <- CreateSeuratObject(counts = read_counts)

head(cell_info)

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, cell_info[match(colnames(seurat_obj), cell_info$Barcode), -1])

head(seurat_obj@meta.data)

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500)

seurat_obj

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(seurat_obj, reduction = "pca")

ElbowPlot(seurat_obj)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

head(seurat_obj$seurat_clusters)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

DimPlot(seurat_obj, reduction = "umap")
DimPlot(seurat_obj, reduction = "umap",group.by = "celltype_sub")

seurat_obj$sample <- seurat_obj$Sample
seurat_obj$Sample <- NULL

seurat_obj$condition <- seurat_obj$Tissue

seurat_obj$celltype <- seurat_obj$celltype_sub

seurat_obj$NewSample <- apply(seurat_obj@meta.data, 1, function(row) {
  if (!is.na(row["Donor"]) & !is.na(row["Tissue"])) {
    paste(row["Donor"], row["Tissue"], sep = "_")
  } else {
    NA
  }
})

table(seurat_obj$NewSample)

seurat_obj$Sample <- seurat_obj$NewSample

saveRDS(seurat_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/smart-seq2/GSE140228_seurat.rds")

seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  seurat_obj,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/smart-seq2/GSE140228.h5ad"
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
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/smart-seq2/GSE140228_sce.rds")









################################################################################
#####################----------      E-MTAB-8560     --------###################


library(MouseThymusAgeing)
smart.all <- MouseSMARTseqData(samples = c("day1","day2","day3","day4","day5"))
smart.all
library(scater)
smart.all <- logNormCounts(smart.all)
seurat_obj <- as.Seurat(smart.all, counts = "counts", data = "logcounts")
counts_mat <- assay(smart.all, "counts")
rownames(counts_mat)[rownames(counts_mat) == "" | is.na(rownames(counts_mat))]
counts_mat <- counts_mat[!(rownames(counts_mat) == "" | is.na(rownames(counts_mat))), ]
rownames(counts_mat) <- make.unique(rownames(counts_mat))
sce_new <- SingleCellExperiment(
  assays = list(counts = counts_mat),
  colData = colData(smart.all)
)
sce_new
smart.seurat <- as.Seurat(
  sce_new,
  counts = "counts",
  data = NULL
)
smart.seurat[["RNA"]] <- smart.seurat[["originalexp"]]
DefaultAssay(smart.seurat) <- "RNA"
smart.seurat[["originalexp"]] <- NULL
smart.seurat <- NormalizeData(smart.seurat)
smart.seurat <- FindVariableFeatures(smart.seurat, selection.method = "vst", nfeatures = 2000)
smart.seurat <- ScaleData(smart.seurat, features = VariableFeatures(smart.seurat))
smart.seurat <- RunPCA(smart.seurat, features = VariableFeatures(smart.seurat))
smart.seurat <- RunUMAP(smart.seurat, dims = 1:30)
smart.seurat <- FindNeighbors(smart.seurat, dims = 1:30)
smart.seurat <- FindClusters(smart.seurat, resolution = 0.5)
smart.seurat$condition <- smart.seurat$Age
smart.seurat$celltype <- smart.seurat$SubType
smart.seurat$Age <- NULL
smart.seurat$SubType <- NULL
smart.seurat$sample <- paste(smart.seurat$condition,
                              smart.seurat$SortDay,
                              sep = "_")
DimPlot(smart.seurat, reduction = "umap", group.by = "seurat_clusters")
DimPlot(smart.seurat, reduction = "umap", group.by = "condition")
saveRDS(smart.seurat, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/smart-seq2/E_MTAB_8560_seurat.rds")
smart.seurat[["RNA"]] <- as(smart.seurat[["RNA"]], "Assay")
library(sceasy)
library(reticulate)
sceasy::convertFormat(
  smart.seurat,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/smart-seq2/E_MTAB_8560.h5ad"
)
library(SingleCellExperiment)
smart.seurat[["RNA"]] <- as(smart.seurat[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(smart.seurat)
if ("pca" %in% names(smart.seurat@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(smart.seurat, "pca")
}
if ("umap" %in% names(smart.seurat@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(smart.seurat, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/smart-seq2/E_MTAB_8560_sce.rds")











################################################################################
#####################----------      GSE146771     --------####################



metadata <- read.table(gzfile("Downloads/GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz"), 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

tpm_data <- read.table(gzfile("Downloads/GSE146771_CRC.Leukocyte.Smart-seq2.TPM.txt.gz"), 
                       header = TRUE,  stringsAsFactors = FALSE)

head(metadata)

head(tpm_data)

library(Seurat)

seurat_obj <- CreateSeuratObject(counts = tpm_data, meta.data = metadata)

seurat_obj

seurat_obj@meta.data <- seurat_obj@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "Sample", "Tissue","Barcode","CellName" ,"Tissue","Global_Cluster","Sub_Cluster","Sub_ClusterID")]

colnames(seurat_obj@meta.data)

seurat_obj <- NormalizeData(seurat_obj)

seurat_obj

seurat_obj <- FindVariableFeatures(seurat_obj)

head(VariableFeatures(seurat_obj))

seurat_obj <- ScaleData(seurat_obj)

seurat_obj <- RunPCA(seurat_obj)

DimPlot(seurat_obj, reduction = "pca")

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)

seurat_obj <- FindClusters(seurat_obj)

DimPlot(seurat_obj, reduction = "pca", group.by = "seurat_clusters")
DimPlot(seurat_obj, reduction = "pca", group.by = "Sample")
DimPlot(seurat_obj, reduction = "pca", group.by = "Tissue")

seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)

DimPlot(seurat_obj, reduction = "tsne")

DimPlot(seurat_obj, reduction = "pca", group.by = "Global_Cluster",label = TRUE)

DimPlot(seurat_obj, reduction = "tsne", group.by = "Global_Cluster",label = TRUE)

table(seurat_obj$Tissue,seurat_obj$Global_Cluster)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

DimPlot(seurat_obj, reduction = "umap")
DimPlot(seurat_obj, reduction = "umap",group.by = "celltype")

seurat_obj$sample <- paste(seurat_obj$Sample, seurat_obj$Tissue, sep = "_")

head(seurat_obj$sample)

table(seurat_obj$sample)

table(seurat_obj$sample,seurat_obj$Global_Cluster)

seurat_obj$celltype <- seurat_obj$Global_Cluster

seurat_obj$condition <- seurat_obj$Tissue

seurat_obj$Global_Cluster <- NULL
seurat_obj$Tissue <- NULL

saveRDS(seurat_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/smart-seq2/GSE146771_seurat.rds")

seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  seurat_obj,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/smart-seq2/GSE146771.h5ad"
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
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/smart-seq2/GSE146771_sce.rds")
