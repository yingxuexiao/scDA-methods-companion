#10X_real_datasets_processed.R


################################################################################ 
###########################     covid_chua   ###################################

chua_obj <- readRDS("./covid_nbt_main.rds")

obj <- UpdateSeuratObject(obj)

immune_types <- c("B cell", "CTL","MC", "moDC", "MoD-Ma", "Neu", "NK", "NKT",
                  "NKT-p", "nrMa", "pDC", "rMa", "Treg")
immune_obj <- subset(chua_obj, subset = celltype %in% immune_types)
immune_obj@meta.data$condition <- immune_obj@meta.data$severity
immune_obj1 <- NormalizeData(immune_obj)
immune_obj1 <- FindVariableFeatures(immune_obj1)
immune_obj1 <- ScaleData(immune_obj1)
immune_obj1 <- RunPCA(immune_obj1)
immune_obj1 <- FindNeighbors(immune_obj1, dims = 1:30)
immune_obj1 <- FindClusters(immune_obj1, resolution = 0.5)


immune_obj1 <- RunUMAP(immune_obj1, dims = 1:30)
saveRDS(COVID_19_chua,file="/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/COVID_19_chua.rds")
saveRDS(COVID_19_chua,file="/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/COVID_19_chua_seurat.rds")
COVID_19_chua[["RNA"]] <- as(COVID_19_chua[["RNA"]], Class = "Assay")
COVID_19_chua_sce <- as.SingleCellExperiment(COVID_19_chua)
saveRDS(COVID_19_chua_sce,file="/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/COVID_19_chua_sce.rds")

library(Seurat)
saveRDS(obj, file = "./chua_covid_immune_obj_annotated.rds")


################################################################################ 
######################     Liao_covid_data   ###################################

#From https://github.com/zhangzlab/covid_balf, found the preprocessed Seurat object at 
#http://cells.ucsc.edu/covid19-balf/nCoV.rds. 
#Download this dataset, and find the cell annotation file meta.tsv for this dataset on the website #https://covid19-balf.cells.ucsc.edu.

nCoV <- readRDS("~/Downloads/nCoV.rds")

nCoV <-  UpdateSeuratObject(nCoV)

DimPlot(nCoV, reduction = "umap", group.by = "sample", label = FALSE) +
   ggtitle("UMAP by Sample")

meta <- read.delim("/Users/xiaoying/Downloads/meta.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(meta$ID)
head(colnames(nCoV))

sum(meta$ID %in% colnames(nCoV))

rownames(meta) <- meta$ID
meta <- meta[colnames(nCoV), ]
nCoV@meta.data$celltype <- meta$celltype　

df <- nCoV@meta.data[, c("sample_new", "celltype")]

library(dplyr)

condition_map <- c(
  HC1 = "healthy", HC2 = "healthy", HC3 = "healthy", HC4 = "healthy",
  O1 = "moderate", O2 = "moderate", O3 = "moderate",
  S1 = "severe/critical", C1 = "severe/critical", C2 = "severe/critical",
  C3 = "severe/critical", C4 = "severe/critical", C5 = "severe/critical"
)

sample_vec <- nCoV@meta.data$sample_new


nCoV@meta.data$condition <- condition_map[sample_vec]

table(nCoV@meta.data$condition)

saveRDS(ncov, file = "ncov_seurat_object.rds")


################################################################################ 
######################     GSE118697   #########################################

library(Seurat)
library(scran)
library(Matrix)

# Load the UMI count matrices
count_1_month <- read.csv("/Users/xiaoying/Downloads/GSE118697/GSE118697_UMI_count_1_month_organoids.csv.gz", row.names = 1)
count_5_month <- read.csv("/Users/xiaoying/Downloads/GSE118697/GSE118697_UMI_count_5_month_organoids.csv.gz", row.names = 1)

# Load sample metadata
metadata_1_month <- read.csv("/Users/xiaoying/Downloads/GSE118697/GSE118697_sample_barcode_metadata_1_month_organoids.csv.gz")
metadata_5_month <- read.csv("/Users/xiaoying/Downloads/GSE118697/GSE118697_sample_barcode_metadata_5_month_organoids.csv.gz")

seurat_obj_1_month <- CreateSeuratObject(counts = count_1_month, meta.data = metadata_1_month)

seurat_obj_5_month <- CreateSeuratObject(counts = count_5_month, meta.data = metadata_5_month)

seurat_obj_1_month <- NormalizeData(seurat_obj_1_month, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj_1_month <- FindVariableFeatures(seurat_obj_1_month)

seurat_obj_1_month <- ScaleData(seurat_obj_1_month)

seurat_obj_1_month <- RunPCA(seurat_obj_1_month)

seurat_obj_1_month <- RunUMAP(seurat_obj_1_month,dims = 1:30)

seurat_obj_1_month <- FindNeighbors(seurat_obj_1_month, dims = 1:21)

seurat_obj_1_month <- FindClusters(seurat_obj_1_month, resolution = 1.2)

DimPlot(seurat_obj_1_month, reduction = "pca", group.by = "seurat_clusters")

seurat_obj_1_month$celltype <- seurat_obj_1_month$cluster_names
seurat_obj_1_month$sample <- seurat_obj_1_month$sample_name
seurat_obj_1_month$condition <- seurat_obj_1_month@meta.data[["strain"]]

DimPlot(seurat_obj_1_month, reduction = "umap", group.by = "celltype",label = TRUE)

DimPlot(seurat_obj_1_month, reduction = "umap", group.by = "condition",label = TRUE)

DimPlot(seurat_obj_1_month, reduction = "umap", group.by = "sample",label = TRUE)


library(ggplot2)
library(dplyr)


saveRDS(seurat_obj_1_month, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE118697_seurat.rds")

seurat_obj_1_month[["RNA"]] <- as(seurat_obj_1_month[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

# 转换为 .h5ad 文件
sceasy::convertFormat(
  seurat_obj_1_month,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets//GSE118697.h5ad"
)

library(SingleCellExperiment)

seurat_obj_1_month[["RNA"]] <- as(seurat_obj_1_month[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seurat_obj_1_month)

if ("pca" %in% names(seurat_obj_1_month@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(seurat_obj_1_month, "pca")
}
if ("umap" %in% names(seurat_obj_1_month@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seurat_obj_1_month, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE118697_sce.rds")


seurat_obj_5_month <- NormalizeData(seurat_obj_5_month, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj_5_month <- FindVariableFeatures(seurat_obj_5_month)

seurat_obj_5_month <- ScaleData(seurat_obj_5_month)

seurat_obj_5_month <- RunPCA(seurat_obj_5_month)

seurat_obj_5_month <- RunUMAP(seurat_obj_5_month,dims = 1:30)

seurat_obj_5_month <- FindNeighbors(seurat_obj_5_month, dims = 1:21)

seurat_obj_5_month <- FindClusters(seurat_obj_5_month, resolution = 1.2)

DimPlot(seurat_obj_5_month, reduction = "umap", group.by = "seurat_clusters")

seurat_obj_5_month$celltype <- seurat_obj_5_month$cluster_name
table(seurat_obj_5_month$sample_id)

seurat_obj_5_month$sample <- seurat_obj_5_month$sample_id
seurat_obj_5_month$condition <- seurat_obj_5_month$strain_id

DimPlot(seurat_obj_5_month, reduction = "umap", group.by = "celltype",label = TRUE)

DimPlot(seurat_obj_5_month, reduction = "umap", group.by = "condition",label = TRUE)

DimPlot(seurat_obj_5_month, reduction = "umap", group.by = "sample",label = TRUE)

saveRDS(seurat_obj_5_month, file = "/Users/xiaoying/Desktop/GSE118697_5_seurat.rds")

seurat_obj_5_month[["RNA"]] <- as(seurat_obj_5_month[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

# 转换为 .h5ad 文件
sceasy::convertFormat(
  seurat_obj_5_month,
  from = "seurat",
  to = "anndata",
  outFile = "/Users/xiaoying/Desktop/GSE118697_5.h5ad"
)

library(SingleCellExperiment)

seurat_obj_5_month[["RNA"]] <- as(seurat_obj_5_month[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seurat_obj_5_month)

if ("pca" %in% names(seurat_obj_5_month@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(seurat_obj_5_month, "pca")
}
if ("umap" %in% names(seurat_obj_5_month@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seurat_obj_5_month, "umap")
}
saveRDS(sce_obj, file = "/Users/xiaoying/Desktop/GSE118697_5_sce.rds")



################################################################################ 
######################   GSE122043     #########################################

library(Seurat)
library(Matrix)
library(dplyr)

base_path <- "./GSE122043"

sample_dirs <- c("35_e13control", "36_e13control", "37_e14control", "38_e14control", "39_e13lof", "40_e14lof")
sample_names <- c("E13.5_WT_1", "E13.5_WT_2", "E14.5_WT_1", "E14.5_WT_2", "E13.5_LOF", "E14.5_LOF")

seurat_list <- list()

for (i in seq_along(sample_dirs)) {
  sample_path <- file.path(base_path, sample_dirs[i])
  sample_name <- sample_names[i]
  
  data <- Read10X(data.dir = sample_path)
  
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
  
  seurat_obj$sample <- sample_name
  
  seurat_list[[sample_name]] <- seurat_obj
  
  message("✅ ：", sample_name)
  
  seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 1000 & nCount_RNA > 2500 & nCount_RNA < 50000)
  
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA"))  
  
  hvf <- VariableFeatures(seurat_obj)
  gene_var <- HVFInfo(seurat_obj)
  hvf_sel <- rownames(gene_var[gene_var$variance.standardized > 0.8, ])
  seurat_obj <- RunPCA(seurat_obj, features = hvf_sel)

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)
 
  seurat_list[[sample_name]] <- seurat_obj

  message("✅ ：", sample_name)
}


obj13_1 <- seurat_list[["E13.5_WT_1"]]

DotPlot(obj13_1, features = c("Col1a1", "Krt14", "Krt10"), group.by = "seurat_clusters") +
  RotatedAxis() +
  ggtitle("E13.5_WT_1 marker expression per cluster")


obj13_1 <- seurat_list[["E13.5_WT_1"]]
obj13_1$celltype <- "Unknown"

obj13_1$celltype[obj13_1$seurat_clusters %in% c(0, 1, 3)] <- "Dermal"

obj13_1$celltype[obj13_1$seurat_clusters %in% c(2)] <- "Keratinocyte"

obj13_1$celltype[obj13_1$seurat_clusters %in% c(4, 5)] <- "Doublet"



obj13_2 <- seurat_list[["E13.5_WT_2"]]

DotPlot(obj13_2, features = c("Col1a1", "Krt14", "Krt10"), group.by = "seurat_clusters") +
  RotatedAxis() +
  ggtitle("E13.5_WT_1 marker expression per cluster")

obj13_2$celltype <- "Unknown"

obj13_2$celltype[obj13_2$seurat_clusters %in% c(0, 1)] <- "Dermal"

obj13_2$celltype[obj13_2$seurat_clusters %in% c(2)] <- "Keratinocyte"

obj13_2$celltype[obj13_2$seurat_clusters %in% c(3,4, 5,6)] <- "Doublet"



obj14_1 <- seurat_list[["E14.5_WT_1"]]

DotPlot(obj14_1, features = c("Col1a1", "Krt14", "Krt10"), group.by = "seurat_clusters") +
  RotatedAxis() +
  ggtitle("E13.5_WT_1 marker expression per cluster")

FeaturePlot(obj14_1, features = c("Col1a1", "Krt14", "Krt10"))
DimPlot(obj14_1, reduction = "tsne", label = TRUE)

obj14_1$celltype <- "Unknown"

obj14_1$celltype[obj14_1$seurat_clusters %in% c(0, 1 )] <- "Dermal"

obj14_1$celltype[obj14_1$seurat_clusters %in% c(2,3)] <- "Keratinocyte"

obj14_1$celltype[obj14_1$seurat_clusters %in% c(4, 5,6,7)] <- "Doublet"




obj14_2 <- seurat_list[["E14.5_WT_2"]]

DotPlot(obj14_2, features = c("Col1a1", "Krt14", "Krt10"), group.by = "seurat_clusters") +
  RotatedAxis() +
  ggtitle("E14.5_WT_2 marker expression per cluster")

FeaturePlot(obj14_2, features = c("Col1a1", "Krt14", "Krt10"))
DimPlot(obj14_2, reduction = "tsne", label = TRUE)

obj14_2$celltype <- "Unknown"

obj14_2$celltype[obj14_2$seurat_clusters %in% c(0, 2 )] <- "Dermal"

obj14_2$celltype[obj14_2$seurat_clusters %in% c(1,3)] <- "Keratinocyte"

obj14_2$celltype[obj14_2$seurat_clusters %in% c(4, 5,6,7,8)] <- "Doublet"



obj13_1_dermal <- subset(obj13_1, subset = celltype == "Dermal")
obj13_2_dermal <- subset(obj13_2, subset = celltype == "Dermal")
obj14_1_dermal <- subset(obj14_1, subset = celltype == "Dermal")
obj14_2_dermal <- subset(obj14_2, subset = celltype == "Dermal")


wt_dermal_merged <- merge(obj13_1_dermal, y = list(obj13_2_dermal, obj14_1_dermal, obj14_2_dermal),
                          add.cell.ids = c("E13.5_WT_1", "E13.5_WT_2", "E14.5_WT_1", "E14.5_WT_2"),
                          project = "WT_Dermal")


# Normalize
wt_dermal_merged <- NormalizeData(wt_dermal_merged)

wt_dermal_merged <- FindVariableFeatures(wt_dermal_merged, selection.method = "vst", nfeatures = 2000)
wt_dermal_merged <- ScaleData(wt_dermal_merged, vars.to.regress = c("nCount_RNA", "sample"))

# PCA
wt_dermal_merged <- RunPCA(wt_dermal_merged)

wt_dermal_merged <- FindNeighbors(wt_dermal_merged, dims = 1:10)
wt_dermal_merged <- FindClusters(wt_dermal_merged, resolution = 0.1)
wt_dermal_merged <- RunTSNE(wt_dermal_merged, dims = 1:10)
wt_dermal_merged <- RunUMAP(wt_dermal_merged, dims = 1:10)

DimPlot(wt_dermal_merged, group.by = "sample", reduction = "tsne", label = TRUE)
DimPlot(wt_dermal_merged, group.by = "seurat_clusters", reduction = "tsne", label = TRUE)


library(dplyr)

wt_dermal_merged <- JoinLayers(wt_dermal_merged)

markers <- FindAllMarkers(
  object = wt_dermal_merged,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"  
)

markers_sorted <- markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj)

head(markers_sorted, 20)

marker_list_by_cluster <- split(markers_sorted, markers_sorted$cluster)

write.csv(markers_sorted, file = "dermal_cluster_markers_sorted.csv", row.names = FALSE)


top20_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  arrange(cluster, desc(avg_log2FC))


print(top20_per_cluster, n = Inf)

wt_dermal_merged$condition <- ifelse(
  wt_dermal_merged$sample %in% c("E13.5_WT_1", "E13.5_WT_2"),
  "E13.5",
  "E14.5"
)


##Since five subgroups were determined based on the resolution, but there were no marker annotations, 
#and we only needed to know the changes in the number of subgroups, 
#we used manual naming as dermal_subcell_1, 2,3,4,5

library(dplyr)

cluster_vec <- as.character(wt_dermal_merged$seurat_clusters)

wt_dermal_merged$subcelltype <- recode(cluster_vec,
  "0" = "dermal_subcell_1",
  "1" = "dermal_subcell_2",
  "2" = "dermal_subcell_3",
  "3" = "dermal_subcell_4",
  "4" = "dermal_subcell_5"
)

table(wt_dermal_merged$subcelltype)
head(wt_dermal_merged@meta.data)

saveRDS(wt_dermal_merged, file = "/Volumes/XiaoYingxue/DA_datasets/10x_datasets/GSE122043_dermal_cell_obj_annotated.rds")






################################################################################ 
################################    GSE129788     ##############################


library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

file_paths <- list.files("/Volumes/XiaoYingxue/DA_datasets/raw_data/10x_datasets/GSE129788_RAW",
                         full.names = TRUE, pattern = "*.txt.gz")

# 样本名映射
sample_rename <- c(
  "GSM3722100" = "Young_01", "GSM3722101" = "Young_02",
  "GSM3722102" = "Young_03", "GSM3722103" = "Young_04",
  "GSM3722104" = "Young_05", "GSM3722105" = "Young_06",
  "GSM3722106" = "Young_07", "GSM3722107" = "Young_08",
  "GSM3722108" = "Old_01",   "GSM3722109" = "Old_02",
  "GSM3722110" = "Old_03",   "GSM3722111" = "Old_04",
  "GSM3722112" = "Old_05",   "GSM3722113" = "Old_06",
  "GSM3722114" = "Old_07",   "GSM3722115" = "Old_08"
)


seurat_list <- list()

for (file in file_paths) {
  gsm_id <- sub("_.*", "", basename(file))     
  sample_name <- sample_rename[gsm_id]        

  message("read: ", basename(file), " -> ", sample_name)

  mat <- read.table(file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  colnames(mat) <- paste(sample_name, colnames(mat), sep = "_")  

  seu <- CreateSeuratObject(counts = mat, project = sample_name)
  seurat_list[[sample_name]] <- seu
}

seu_merged <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:16],
  add.cell.ids = names(seurat_list),
  project = "GSE129788"
)



library(stringr)
seu_merged$sample <- str_extract(colnames(seu_merged), "^[^_]+_[^_]+")

seu_merged$age_group <- ifelse(grepl("^Young", seu_merged$sample), "Young", "Old")

print(table(seu_merged$sample))

print(table(seu_merged$age_group))


annot <- read.table("/Volumes/XiaoYingxue/DA_datasets/raw_data/10x_datasets/GSE129788_Supplementary_meta_data_Cell_Types_Etc.txt.gz",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

annot <- annot[-1, ]

annot$cell_suffix <- sub("^.*_(\\d+_\\w+)$", "\\1", annot$NAME)

head(annot$cell_suffix)
table(nchar(annot$cell_suffix))  

library(stringr)

seu_suffix <- str_extract(colnames(seu_merged), "\\d+_\\w{15,16}$")
length(seu_suffix) 
head(seu_suffix)



seu_cells_df <- data.frame(
  cell_name = colnames(seu_merged),
  cell_suffix = seu_suffix,
  stringsAsFactors = FALSE
)

annot_reduced <- annot[, c("cell_suffix", "cluster", "animal_type", "cell_classes", "cell_type_age")]

merged_annot <- merge(seu_cells_df, annot_reduced, by = "cell_suffix", all.x = TRUE, sort = FALSE)

merged_annot <- merged_annot[match(seu_cells_df$cell_name, merged_annot$cell_name), ]

all(merged_annot$cell_name == seu_cells_df$cell_name)  # 应该为 TRUE

seu_merged@meta.data$cluster <- merged_annot$cluster
seu_merged@meta.data$animal_type <- merged_annot$animal_type
seu_merged@meta.data$cell_classes <- merged_annot$cell_classes
seu_merged@meta.data$cell_type_age <- merged_annot$cell_type_age

seu_merged <- NormalizeData(seu_merged, normalization.method = "LogNormalize", scale.factor = 10000)
seu_merged <- FindVariableFeatures(seu_merged, selection.method = "vst", nfeatures = 2000)

VariableFeaturePlot(seu_merged)
seu_merged <- ScaleData(seu_merged, features = VariableFeatures(seu_merged))
seu_merged <- RunPCA(seu_merged, features = VariableFeatures(seu_merged), npcs = 40)

ElbowPlot(seu_merged, ndims = 40)  # 检查前多少个 PC 解释较多变异

seu_merged <- FindNeighbors(seu_merged, dims = 1:20)
seu_merged <- FindClusters(seu_merged, resolution = 2.0)
table(Idents(seu_merged))
seu_merged <- RunUMAP(seu_merged, dims = 1:20)

DimPlot(seu_merged, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
  ggtitle("UMAP clustering (resolution = 2.0)")

seu_merged <- RunTSNE(seu_merged, dims = 1:20)

saveRDS(seu_merged, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/GSE129788_obj_annotation.rds")



################################################################################ 
############################      GSE135337    #################################
library(Seurat)
library(readr)
library(readxl)

file_paths <- list.files("/Users/xiaoying/Downloads/GSE135337_RAW/", pattern = "gene_cell_exprs_table", full.names = TRUE)

seurat_list <- list()


for (file_path in file_paths) {
  if (grepl(".txt.gz$", file_path)) {
    data <- read_delim(file_path, delim = "\t", col_names = TRUE)
  } else if (grepl(".xls$", file_path)) {
    data <- read_excel(file_path)
  }
  
  symbol_counts <- table(data$Symbol)
  duplicated_genes <- symbol_counts[symbol_counts > 1]  
  cat("Sample:", gsub(".*/(.*)_gene_cell_exprs_table.*", "\\1", file_path), "\n")
  cat("Number of duplicated genes:", length(duplicated_genes), "\n")
  cat("Duplicated genes:", names(duplicated_genes), "\n")
  
 
  data_unique <- data[!duplicated(data$Symbol), ]
    cat("After removing duplicates, number of genes:", nrow(data_unique), "\n")
  
  gene_data <- data_unique[, -c(1, 2)]  
  
  rownames(gene_data) <- data_unique$Symbol
  
  gene_data_matrix <- as.matrix(gene_data)
  
  sample_name <- gsub(".*/(.*)_gene_cell_exprs_table.*", "\\1", file_path)  
  seurat_obj <- CreateSeuratObject(counts = gene_data_matrix, project = sample_name)
  
  seurat_list[[sample_name]] <- seurat_obj
  
  print(seurat_obj)
  cat("Sample:", sample_name, "- Number of genes:", nrow(seurat_obj), "Cells:", ncol(seurat_obj), "\n")
}

seurat_list$GSM5329919_BCN <- NULL

cat("Total number of Seurat objects:", length(seurat_list), "\n")


combined_seurat <- Reduce(function(x, y) {
  merge(x, y, project = "Combined_BC")
}, seurat_list)


View(combined_seurat)

sample_info <- gsub("GSM[0-9]+_(.*)", "\\1", combined_seurat@meta.data[["orig.ident"]])

combined_seurat <- AddMetaData(combined_seurat, metadata = sample_info, col.name = "sample")

head(combined_seurat@meta.data[["sample"]])

table(combined_seurat@meta.data[["sample"]])

mt_genes <- grep("^MT-", rownames(combined_seurat), value = TRUE)

combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, features = mt_genes)
head(combined_seurat@meta.data[["percent.mt"]])

combined_seurat <- subset(combined_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

cat("Number of cells after quality control:", ncol(combined_seurat), "\n")
cat("Number of genes after quality control:", nrow(combined_seurat), "\n")


library(Seurat)
library(dplyr)
library(SingleCellExperiment)

combined_seurat <- SCTransform(combined_seurat, verbose = FALSE)

combined_seurat <- RunPCA(combined_seurat, npcs = 30, verbose = FALSE)

ElbowPlot(combined_seurat)

combined_seurat <- RunUMAP(combined_seurat, dims = 1:20) 

combined_seurat <- FindNeighbors(combined_seurat, dims = 1:20)  
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)  

DimPlot(combined_seurat, reduction = "umap", label = TRUE)


interest_marker <- c("EPCAM", "KRT8", "KRT18","CD14", "CSF1R", "AIF1","CD2", "CD3D", "CD3E",
                     "DCN", "PDPN", "TAGLN","PECAM1", "VWF", "CLDN5","PLVAP","SPARCL1","GNG11",
                     "LUM","COL1A1","COL1A2")


library(ggplot2)

dotplot <- DotPlot(combined_seurat, features = interest_marker, group.by = "seurat_clusters") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Gene Expression across Cell Types")

marker_genes <- list(
  urothelial = c("EPCAM", "KRT8", "KRT18"),
  myeloid = c("CD14", "CSF1R", "AIF1"),
  T_cells = c("CD2", "CD3D", "CD3E"),
  fibroblasts = c("DCN", "PDPN", "TAGLN"),
  endothelial = c("PECAM1", "VWF", "CLDN5")
)


cell_type_annotation <- list(
  urothelial = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,16),
  myeloid = c(14),
  T_cells = c(17),
  fibroblasts = c(18),
  endothelial = c(19)
)

cell_type_vector <- rep(NA, length(colnames(combined_seurat)))  

for (celltype in names(cell_type_annotation)) {
  clusters <- cell_type_annotation[[celltype]]
  
  cell_indices <- which(combined_seurat$seurat_clusters %in% clusters)  
  cell_type_vector[cell_indices] <- celltype  
}

combined_seurat$celltype <- cell_type_vector

table(combined_seurat$celltype)  

table(combined_seurat$sample)

DimPlot(combined_seurat, reduction = "umap",group.by="celltype", label = TRUE)


condition_info <- data.frame(
  sample = c("BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "BC7"),
  condition = c("pTa", "pT1", "pT1", "pT2N0M0", "pT3N0M0", "pTa", "pTa")  
)

print(condition_info)

combined_seurat$condition <- condition_info$condition[match(combined_seurat$sample, condition_info$sample)]

head(combined_seurat@meta.data[["condition"]])

table(combined_seurat$condition)

table(combined_seurat$condition,combined_seurat$celltype)

saveRDS(combined_seurat, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE135337_seurat.rds")

DefaultAssay(combined_seurat) 

counts_matrix <- GetAssayData(combined_seurat, assay = "SCT", slot = "counts")
combined_seurat[["RNA"]] <- CreateAssayObject(counts = counts_matrix)
DefaultAssay(combined_seurat) <- "RNA"

combined_seurat[["RNA"]] <- as(combined_seurat[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  combined_seurat,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE135337.h5ad"
)

library(SingleCellExperiment)

combined_seurat[["RNA"]] <- as(combined_seurat[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(combined_seurat)

if ("pca" %in% names(combined_seurat@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(combined_seurat, "pca")
}
if ("umap" %in% names(combined_seurat@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(combined_seurat, "umap")
}

saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE135337_sce.rds")


################################################################################ 
############################      GSE136103    #################################

library(Seurat)
library(ggplot2)
library(pheatmap)
library(grid)
library(dplyr)

data_dir <- "/Users/xiaoying/Downloads/GSE136103_RAW"

sample_folders <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)

sample_folders <- sample_folders[!sample_folders %in% c("", "filtered_feature_bc_matrix")]

seurat_list <- list()

for (sample_name in sample_folders) {
  sample_path <- file.path(data_dir, sample_name) 
  
  data <- Read10X(data.dir = sample_path)
  
  combined <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
  
  combined$sample <- sample_name
  
  seurat_list[[sample_name]] <- combined
  
  message("✅ ：", sample_name)
}

combined <- merge(seurat_list[[1]], y = seurat_list[-1])

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

combined <- subset(combined, subset = nCount_RNA >300 & nCount_RNA < 25000 & percent.mt < 30)

combined <- NormalizeData(combined, normalization.method = "RC", scale.factor = 1e6)
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 1e6)

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

combined <- RunUMAP(combined, dims = 1:50)
combined <- FindNeighbors(combined, dims = 1:15)
combined <- FindClusters(combined, resolution = 0.5)

DimPlot(combined, reduction = "umap", group.by = "sample", label = TRUE)
DimPlot(combined, reduction = "umap", label = TRUE)

combined$condition <- ifelse(grepl("cirrhotic", combined$sample), "cirrhotic", "healthy")

head(combined@meta.data)
DimPlot(combined, reduction = "umap",group.by = "condition")

marker_genes <- list(
  "MP" = c("CD68", "ITGAM", "ITGAX","CSF1R", "CD14"),
  "pDC" = c("LILRA4", "CLEC4C", "GZMB"),
  "ILC" = c("KLRF1", "KLRC1", "GZMA",  "NKG7"),
  "T cell" = c("CD3D", "CD3E", "CD3G", "CD8A"),
  "B cell" = c( "CD79B", "CD19", "MS4A1"),
  "Plasma cell" = c("CD79A", "IGHA2"),
  "Mast cell" = c("KIT", "TPSAB1", "TPSB2","CPA3"),
  "Endothelia" = c("PECAM1", "CDH5", "ICAM2", "KDR", "ERG"),
  "Mesenchyme" = c("PDGFRB", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "DES", "DCN"),
  "Hepatocyte" = c("ALB", "TF", "TTR", "HNF4A", "CYP2A6"),
  "Cholangiocyte" = c("EPCAM", "KRT19", "CD24"),
  "Cycling" = c("MKI67", "CCNA2", "CCNB2", "STMN1")
)

print(dot_plot)

DotPlot(combined, features = marker_genes)

DotPlot(object = combined, features = marker_genes) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


  cell_type_map <- list(
    "MP" = c(4,8,13,18,19),
    "pDC" = c(16),
    "ILC" = c(3,5,6),
    "T cell" = c(0,1,2,7),
    "B cell" = c(9),
    "Plasma cell" = c(17),
    "Mast cell" = c(11),
    "Endothelia" = c(10),
    "Mesenchyme" = c(15),
    "Hepatocyte" = c(12),
    "Cycling" = c(14)
  )
  
  combined$celltype <- "Unknown" 
  
  for (celltype in names(cell_type_map)) {
    clusters <- cell_type_map[[celltype]]
    combined$celltype[combined$seurat_clusters %in% clusters] <- celltype
  }
    head(combined@meta.data)
  
  dot_plot <- DotPlot(combined, 
                      features = unlist(marker_genes),  
                      group.by = "celltype",  
                      cols = c("lightgray", "blue"),  
                      dot.scale = 8) +  
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title = "Marker Gene Expression in Different Cell Types")  
  
  
print(dot_plot)

DimPlot(combined, reduction = "umap",group.by = "celltype",label = TRUE)

table(combined$celltype,combined$sample)
table(combined$sample,combined$condition)


saveRDS(combined, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE136103_seurat.rds")

combined[["RNA"]] <- as(combined[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  combined,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE136103.h5ad"
)

library(SingleCellExperiment)

combined[["RNA"]] <- as(combined[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(combined)

if ("pca" %in% names(combined@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(combined, "pca")
}
if ("umap" %in% names(combined@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(combined, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE136103_sce.rds")



################################################################################ 
############################      GSE137971    #################################

library(data.table)
library(Seurat)
library(AnnotationDbi)
#BiocManager::install("org.Dr.eg.db")
library(org.Dr.eg.db)

folder_path <- "/Users/xiaoying/Downloads/GSE137971_RAW/"

file_paths <- list.files(folder_path, pattern = ".*_DGEmatrix.csv.gz$", full.names = TRUE)

seurat_list <- list()

for (file_path in file_paths) {
  
  data <- fread(file_path)
  
  gene_expression_matrix <- data[, -1, with = FALSE]
  
  gene_expression_matrix <- t(as.matrix(gene_expression_matrix))
  
  colnames(gene_expression_matrix) <- data$Barcodes
  
  seurat_obj <- CreateSeuratObject(counts = gene_expression_matrix)
  
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  seurat_list[[file_name]] <- seurat_obj
}

names(seurat_list)

names(seurat_list) <- gsub(".*_(samp\\d+|\\d+dpa\\d+)_DGEmatrix.*", "\\1", names(seurat_list))

names(seurat_list)


dpa_samples <- names(seurat_list)[grep("dpa", names(seurat_list))]  
samp_samples <- names(seurat_list)[grep("samp", names(seurat_list))]  

dpa_seurat_list <- seurat_list[dpa_samples]
samp_seurat_list <- seurat_list[samp_samples]

dpa_seu <- merge(
  x = dpa_seurat_list[[1]], 
  y = dpa_seurat_list[-1],
  add.cell.ids = names(dpa_seurat_list)
)

samp_seu <- merge(
  x = samp_seurat_list[[1]], 
  y = samp_seurat_list[-1],
  add.cell.ids = names(samp_seurat_list)
)


sample_info <- sapply(strsplit(colnames(dpa_seu), "_"), `[`, 1)
dpa_seu$sample <- sample_info

table(dpa_seu$sample)

sample_info <- sapply(strsplit(colnames(samp_seu), "_"), `[`, 1)
samp_seu$sample <- sample_info

table(samp_seu$sample)



batch_convert_genes <- function(ensembl_ids, batch_size = 1000) {
  batches <- split(ensembl_ids, ceiling(seq_along(ensembl_ids)/batch_size))
  results <- list()
  
  for(i in seq_along(batches)) {
    message("process ", i, "/", length(batches))
    batch_result <- select(org.Dr.eg.db, keys = batches[[i]], 
                           columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                           keytype = "ENSEMBL")
    results[[i]] <- batch_result
  }
  
  do.call(rbind, results)
}

process_seurat <- function(seurat_obj) {
  
  ensembl_ids <- rownames(seurat_obj)
  
  gene_map <- batch_convert_genes(ensembl_ids)
  
  gene_map <- gene_map[!is.na(gene_map$SYMBOL) & gene_map$SYMBOL != "", ]
  gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]
  gene_map$SYMBOL <- make.unique(gene_map$SYMBOL)
  
  id_to_symbol <- setNames(gene_map$SYMBOL, gene_map$ENSEMBL)
  
  mapped_genes <- ensembl_ids[ensembl_ids %in% names(id_to_symbol)]
  message("可映射的基因比例: ", round(length(mapped_genes)/length(ensembl_ids)*100, 2), "%")
  
  seurat_obj <- JoinLayers(seurat_obj)
  counts_matrix <- GetAssayData(seurat_obj, slot = "counts")
  
  counts_mapped <- counts_matrix[mapped_genes, ]
  
  rownames(counts_mapped) <- id_to_symbol[rownames(counts_mapped)]
  
  seurat_symbol <- CreateSeuratObject(counts = counts_mapped, meta.data = seurat_obj@meta.data)
  
  seurat_symbol[["percent.mt"]] <- PercentageFeatureSet(seurat_symbol, pattern = "^mt-")
  seurat_symbol[["percent.ribo"]] <- PercentageFeatureSet(seurat_symbol, pattern = "^rp[sl]")
  
  VlnPlot(seurat_symbol, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", pt.size = 0.1)
  
  filtered_seurat <- subset(seurat_symbol, 
                            subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 6000 & 
                              percent.mt < 20)
  
  filtered_seurat <- NormalizeData(filtered_seurat)
  filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)
  
  filtered_seurat <- ScaleData(filtered_seurat, verbose = FALSE)
  
  filtered_seurat <- RunPCA(filtered_seurat, npcs = 50, verbose = FALSE)
  
  pcs_use <- 1:30
  filtered_seurat <- FindNeighbors(filtered_seurat, dims = pcs_use)
  filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.8)
  
  filtered_seurat <- RunUMAP(filtered_seurat, dims = pcs_use)
  filtered_seurat <- RunTSNE(filtered_seurat, dims = pcs_use)
  
  DimPlot(filtered_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
  DimPlot(filtered_seurat, reduction = "umap", group.by = "sample")
  
  return(filtered_seurat)
}

seurat_objects <- list(dpa_seu = dpa_seu, samp_seu = samp_seu)

processed_seurats <- lapply(seurat_objects, process_seurat)

processed_seurats

dpa_seu <- processed_seurats$dpa_seu
samp_seu <- processed_seurats$samp_seu

library(ggplot2)
key_markers <- c("epcam", "cdh1","krt4", "agr2","foxp1b",
                 "foxp4","tp63","krtt1c19e","mpeg1.1", "cxcr3.2", "msx1b", "twist1a")

FeaturePlot(dpa_seu, features = key_markers, 
            reduction = "umap", ncol = 3)
FeaturePlot(samp_seu, features = key_markers, 
            reduction = "umap", ncol = 3)

DotPlot(dpa_seu, features = key_markers, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(samp_seu, features = key_markers, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DimPlot(dpa_seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(samp_seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE)


cell_type_annotation <- list(
  Mucosal_like =c(19),
  epithelial =c(0,1,2,3,4,6,7,9,10,11,13,14,15,16,18,23),
  hemotopoietic =c(8,21,22),
  Mesenchymal =c(5,12,17,20,24)
)

cell_type_vector <- rep(NA, length(colnames(dpa_seu)))  

for (cell_type in names(cell_type_annotation)) {
  clusters <- cell_type_annotation[[cell_type]]
  
  cell_indices <- which(dpa_seu$seurat_clusters %in% clusters)  
  cell_type_vector[cell_indices] <- cell_type  
}

dpa_seu$celltype <- cell_type_vector

table(dpa_seu$celltype)  

cell_type_annotation <- list(
  Mucosal_like =c(10,9),
  epithelial =c(0,1,2,3,4,5,6,7,8,11,14),
  hemotopoietic =c(15),
  Mesenchymal =c(12,13)
)

cell_type_vector <- rep(NA, length(colnames(samp_seu)))  

for (cell_type in names(cell_type_annotation)) {
  clusters <- cell_type_annotation[[cell_type]]
  
  cell_indices <- which(samp_seu$seurat_clusters %in% clusters)  
  cell_type_vector[cell_indices] <- cell_type  
}

samp_seu$celltype <- cell_type_vector

table(samp_seu$celltype)  

table(samp_seu$celltype,samp_seu$sample)

table(dpa_seu$celltype,dpa_seu$sample)

common_genes <- intersect(rownames(dpa_seu), rownames(samp_seu))

dpa_seu <- dpa_seu[common_genes, ]
samp_seu <- samp_seu[common_genes, ]

merged_seurat <- merge(dpa_seu, samp_seu, add.cell.ids = c("dpa", "samp"))

merged_seurat

merged_seurat$condition <- ifelse(merged_seurat$sample %in% c("samp1", "samp2"), "preinjury", 
                                  ifelse(grepl("^1dpa", merged_seurat$sample), "1 dpa",
                                         ifelse(grepl("^2dpa", merged_seurat$sample), "2 dpa", 
                                                ifelse(grepl("^4dpa", merged_seurat$sample), "4 dpa", NA))))
table(merged_seurat$condition)

table(merged_seurat$condition,merged_seurat$celltype)


merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

merged_seurat <- ScaleData(merged_seurat, verbose = FALSE)

merged_seurat <- RunPCA(merged_seurat, npcs = 50, verbose = FALSE)

pcs_use <- 1:30
merged_seurat <- FindNeighbors(merged_seurat, dims = pcs_use)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.8)

merged_seurat <- RunUMAP(merged_seurat, dims = pcs_use)
merged_seurat <- RunTSNE(merged_seurat, dims = pcs_use)

saveRDS(merged_seurat, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE137971_seurat.rds")

merged_seurat[["RNA"]] <- as(merged_seurat[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  merged_seurat,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE137971.h5ad"
)

library(SingleCellExperiment)

merged_seurat[["RNA"]] <- as(merged_seurat[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(merged_seurat)

if ("pca" %in% names(merged_seurat@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(merged_seurat, "pca")
}
if ("umap" %in% names(merged_seurat@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(merged_seurat, "umap")
}

saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE137971_sce.rds")



################################################################################ 
############################      GSE145980    #################################

nonCM_kit_mutant <- read.csv("Users/xiaoying/Downloads/GSE145980_zebrafish_nonCM_kit_mutant_scRNA_raw_count_matrix.csv.gz",sep=" ", header = TRUE, row.names = 1)

nonCM_wildtype <- read.csv("Users/xiaoying/Downloads/GSE145980_zebrafish_nonCM_wildtype_scRNA_raw_count_matrix.csv.gz", sep=",",header = TRUE, row.names = 1)

library(Seurat)

kit <-  CreateSeuratObject(nonCM_kit_mutant)
wildtype <- CreateSeuratObject(nonCM_wildtype)

sample_info <- sub("^([A-Za-z0-9_]+)_.*", "\\1", colnames(kit))
unique_samples <- unique(sample_info)

print(unique_samples)
kit$sample <- sample_info
head(kit@meta.data)

sample_info_wild <- sub("^([A-Za-z0-9_]+)_.*", "\\1", colnames(wildtype))

unique_samples_wild <- unique(sample_info_wild)
print(unique_samples_wild)
wildtype$sample <- sample_info_wild
head(wildtype@meta.data)


library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

wildtype$condition <- ifelse(grepl("unjuried", wildtype$sample), "unjuried",
                             ifelse(grepl("2dpi", wildtype$sample), "2dpi",
                                    ifelse(grepl("7dpi", wildtype$sample), "7dpi", "14dpi")))

#wildtype$condition <- factor(wildtype$condition, levels = c("unjuried", "2dpi", "7dpi", "14dpi"))

wildtype[["percent.mt"]] <- PercentageFeatureSet(wildtype, pattern = "^mt-")

VlnPlot(wildtype, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample")
plot1 <- FeatureScatter(wildtype, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wildtype, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

wildtype <- subset(wildtype, subset = nFeature_RNA > 200)

cardiomyocyte_genes <- c("ckma", "tnnt2", "nppa")
cardiomyocyte_genes <- cardiomyocyte_genes[cardiomyocyte_genes %in% rownames(wildtype)]
wildtype <- subset(wildtype, cells = colnames(wildtype)[!colnames(wildtype) %in% Cells(wildtype)[apply(GetAssayData(wildtype, slot = "counts")[cardiomyocyte_genes, ] > 0, 2, any)]])

wildtype <- NormalizeData(wildtype, normalization.method = "LogNormalize", scale.factor = 10000)
wildtype <- FindVariableFeatures(wildtype, selection.method = "vst", nfeatures = 2000)


wildtype_list <- SplitObject(wildtype, split.by = "sample")

wildtype_list <- lapply(wildtype_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 2000)
})

anchors <- FindIntegrationAnchors(object.list = wildtype_list, dims = 1:30)

wildtype_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(wildtype_integrated) <- "integrated"

wildtype_integrated <- ScaleData(wildtype_integrated)
wildtype_integrated <- RunPCA(wildtype_integrated, npcs = 50)
wildtype_integrated <- FindNeighbors(wildtype_integrated, dims = 1:30)
wildtype_integrated <- FindClusters(wildtype_integrated, resolution = 0.5)
wildtype_integrated <- RunUMAP(wildtype_integrated, dims = 1:30)


DimPlot(wildtype_integrated,group.by = "seurat_clusters",reduction="umap",label = TRUE)
DimPlot(wildtype_integrated,group.by = "integrated_snn_res.0.5",reduction="umap",label = TRUE)

marker_genes <- list(
  ECs = c("cdh5", "kdrl", "fli1a", "flt1"),
  FBs = c("tcf21", "fn1b", "col1a1a"),
  Mes = c("angptl7", "rspo1", "mgp"),
  MC = c("mpeg1.1", "mfap4", "cd74a","c1qb"),
  Neutro = c("lyz", "mpx"),
  T_NK_B = c("sla2", "irf4b", "ccl36.1", "cxcr4a", "lck", "nkl.2", "zbtb32", "cd79a"),
  Eryth = c("cahz", "slc4a1a"),
  Throm = c("itga2b", "gp1bb")
)

DotPlot(wildtype_integrated, features = marker_genes, dot.scale = 8) + RotatedAxis()

cluster_annotation <- list(
  ECs = c(1,2,7,15,19,20,22),
  FBs = c(3,4,6,12,18),
  Mes = c(5),
  MC = c(0,9,10),
  Neutro = c(14),
  T_NK_B = c(11,13,16),
  Eryth = c(8,21,23),
  Throm = c(17)
)

cell_type_vector <- rep(NA, length(colnames(wildtype_integrated)))  

for (celltype in names(cluster_annotation)) {
  
  clusters <- cluster_annotation[[celltype]]
  
  cell_indices <- which(wildtype_integrated$seurat_clusters %in% clusters) 
  cell_type_vector[cell_indices] <- celltype 
}

wildtype_integrated$celltype <- cell_type_vector

table(wildtype_integrated$celltype) 

table(wildtype_integrated$celltype,wildtype_integrated$sample)

table(wildtype_integrated$celltype,wildtype_integrated$condition)

library(ggplot2)
library(dplyr)

saveRDS(wildtype_integrated, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE145980_wildtype_seurat.rds")

wildtype_integrated[["RNA"]] <- as(wildtype_integrated[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  wildtype_integrated,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE145980_wildtype.h5ad"
)

library(SingleCellExperiment)

wildtype_integrated[["RNA"]] <- as(wildtype_integrated[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(wildtype_integrated)

if ("pca" %in% names(wildtype_integrated@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(wildtype_integrated, "pca")
}
if ("umap" %in% names(wildtype_integrated@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(wildtype_integrated, "umap")
}

saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE145980_wildtype_sce.rds")






###kit
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

kit$condition <- ifelse(grepl("unjuried", kit$sample), "unjuried",
                             ifelse(grepl("2dpi", kit$sample), "2dpi",
                                    ifelse(grepl("7dpi", kit$sample), "7dpi", "14dpi")))

kit[["percent.mt"]] <- PercentageFeatureSet(kit, pattern = "^mt-")

VlnPlot(kit, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample")
plot1 <- FeatureScatter(kit, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(kit, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

kit <- subset(kit, subset = nFeature_RNA > 200)

cardiomyocyte_genes <- c("ckma", "tnnt2", "nppa")
cardiomyocyte_genes <- cardiomyocyte_genes[cardiomyocyte_genes %in% rownames(kit)]
kit <- subset(kit, cells = colnames(kit)[!colnames(kit) %in% Cells(kit)[apply(GetAssayData(kit, slot = "counts")[cardiomyocyte_genes, ] > 0, 2, any)]])

kit <- NormalizeData(kit, normalization.method = "LogNormalize", scale.factor = 10000)
kit <- FindVariableFeatures(kit, selection.method = "vst", nfeatures = 2000)


kit_list <- SplitObject(kit, split.by = "sample")


kit_list <- lapply(kit_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 2000)
})

anchors <- FindIntegrationAnchors(object.list = kit_list, dims = 1:30)

kit_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(kit_integrated) <- "integrated"

kit_integrated <- ScaleData(kit_integrated)
kit_integrated <- RunPCA(kit_integrated, npcs = 50)
kit_integrated <- FindNeighbors(kit_integrated, dims = 1:30)
kit_integrated <- FindClusters(kit_integrated, resolution = 0.5)
kit_integrated <- RunUMAP(kit_integrated, dims = 1:30)


DimPlot(kit_integrated,group.by = "seurat_clusters",reduction="umap",label = TRUE)
DimPlot(kit_integrated,group.by = "integrated_snn_res.0.5",reduction="umap",label = TRUE)

marker_genes <- list(
  ECs = c("cdh5", "kdrl", "fli1a", "flt1"),
  FBs = c("tcf21", "fn1b", "col1a1a"),
  Mes = c("angptl7", "rspo1", "mgp"),
  MC = c("mpeg1.1", "mfap4", "cd74a","c1qb"),
  Neutro = c("lyz", "mpx"),
  T_NK_B = c("sla2", "irf4b", "ccl36.1", "cxcr4a", "lck", "nkl.2", "zbtb32", "cd79a"),
  Eryth = c("cahz", "slc4a1a"),
  Throm = c("itga2b", "gp1bb")
)

DotPlot(kit_integrated, features = marker_genes, dot.scale = 8) + RotatedAxis()

cluster_annotation <- list(
  ECs = c(1,2,3,8,9,13,17,18),
  FBs = c(4,5,15,19),
  Mes = c(7),
  MC = c(0,12,16),
  Neutro = c(14),
  T_NK_B = c(10,21,22),
  Eryth = c(6,20),
  Throm = c(11)
)

cell_type_vector <- rep(NA, length(colnames(kit_integrated)))  

for (celltype in names(cluster_annotation)) {
  clusters <- cluster_annotation[[celltype]]
  

  cell_indices <- which(kit_integrated$seurat_clusters %in% clusters)  
  cell_type_vector[cell_indices] <- celltype  
}

kit_integrated$celltype <- cell_type_vector

table(kit_integrated$celltype)  

table(kit_integrated$celltype,kit_integrated$sample)

table(kit_integrated$celltype,kit_integrated$condition)

saveRDS(kit_integrated, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE145980_kit_seurat.rds")

kit_integrated[["RNA"]] <- as(kit_integrated[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

# 转换为 .h5ad 文件
sceasy::convertFormat(
  kit_integrated,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE145980_kit.h5ad"
)

library(SingleCellExperiment)

kit_integrated[["RNA"]] <- as(kit_integrated[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(kit_integrated)

if ("pca" %in% names(kit_integrated@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(kit_integrated, "pca")
}
if ("umap" %in% names(kit_integrated@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(kit_integrated, "umap")
}

saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE145980_kit_sce.rds")



################################################################################ 
############################      GSE160269    #################################

data <- read.table(gzfile("/media/cqnu/XiaoYingxue/GSE160269_CD45neg_UMIs.txt.gz"),
                   header=TRUE,
                   sep=" ",
                   row.names=1)


seurat_obj<- CreateSeuratObject(counts = data)

cellmeta <- read.table("/media/cqnu/XiaoYingxue/GSE160269_CD45neg_cells.txt.gz",
                       header=TRUE,
                       sep=" ")

##add new cellmeta names which compaired with seurat

cellmeta$cell_new <-  gsub("-",".",cellmeta$cell)


cellmeta_new <- cellmeta[,-which(names(cellmeta) == "cell")]
cellmeta_new<- as.data.frame(cellmeta_new)

seurat_obj1 <- AddMetaData(seurat_obj,metadata=cellmeta_new)

colnames(seurat_obj1)
rownames(seurat_obj1)

seurat_obj1$condition <- ifelse(grepl("T",seurat_obj1@meta.data$sample),"Tumor","Normal")

seurat_obj1$celltype <- seurat_obj1$annotated_type

seurat_obj1 <- NormalizeData(seurat_obj1)

seurat_obj1 <- FindVariableFeatures(seurat_obj1, selection.method = "vst", nfeatures = 2000)

seurat_obj1 <- ScaleData(seurat_obj1)

seurat_obj1 <- RunPCA(seurat_obj1)

print(seurat_obj1[["pca"]], dims = 1:5, nfeatures = 5)

#seurat_obj111 <- subset(seurat_obj1,nFeature_RNA>600)

seurat_obj1 <- FindNeighbors(seurat_obj1,dims=1:10, k.param =30 )  
seurat_obj1 <- FindClusters(seurat_obj1, resolution = 0.6)  

seurat_obj1 <- RunUMAP(seurat_obj1, dims = 1:20)
DimPlot(seurat_obj1,group.by = "seurat_clusters",label = TRUE,reduction = "umap")
DimPlot(seurat_obj1,group.by = "sample",label = TRUE,reduction = "umap")
DimPlot(seurat_obj1,group.by = "celltype",label = TRUE,reduction = "umap")
DimPlot(seurat_obj1,group.by = "condition",label = TRUE,reduction = "umap")

saveRDS(seurat_obj1, file = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE160269_cd45neg_seurat.rds")

seurat_obj1[["RNA"]] <- as(seurat_obj1[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

#use_condaenv("base", required = TRUE)
sceasy::convertFormat(
  seurat_obj1,
  from = "seurat",
  to = "anndata",
  outFile = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE160269_cd45neg.h5ad"
)

library(SingleCellExperiment)

seurat_obj1[["RNA"]] <- as(seurat_obj1[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seurat_obj1)

if ("pca" %in% names(seurat_obj1@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(seurat_obj1, "pca")
}
if ("umap" %in% names(seurat_obj1@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seurat_obj1, "umap")
}
saveRDS(sce_obj, file = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE160269_cd45neg_sce.rds")


data <- read.table(gzfile("/media/cqnu/XiaoYingxue/GSE160269/GSE160269_CD45pos_UMIs.txt.gz"),
                   header=TRUE,
                   sep=" ",
                   row.names=1)


seurat_obj<- CreateSeuratObject(counts = data)

cellmeta <- read.table("/media/cqnu/XiaoYingxue/GSE160269/GSE160269_CD45pos_cells.txt.gz",
                       header=TRUE,
                       sep=" ")

##add new cellmeta names which compaired with seurat

cellmeta$cell_new <-  gsub("-",".",cellmeta$cell)


cellmeta_new <- cellmeta[,-which(names(cellmeta) == "cell")]
cellmeta_new<- as.data.frame(cellmeta_new)

# 添加注释信息到Seurat对象的meta.data

seurat_obj1 <- AddMetaData(seurat_obj,metadata=cellmeta_new)


colnames(seurat_obj1)
rownames(seurat_obj1)

seurat_obj1$condition <- ifelse(grepl("T",seurat_obj1@meta.data$sample),"Tumor","Normal")

seurat_obj1$celltype <- seurat_obj1$annotated_type


seurat_obj1 <- NormalizeData(seurat_obj1)

seurat_obj1 <- FindVariableFeatures(seurat_obj1, selection.method = "vst", nfeatures = 2000)
seurat_obj1 <- ScaleData(seurat_obj1)

seurat_obj1 <- RunPCA(seurat_obj1)

print(seurat_obj1[["pca"]], dims = 1:5, nfeatures = 5)
seurat_obj1 <- FindNeighbors(seurat_obj1,dims=1:10, k.param =30 )  
seurat_obj1 <- FindClusters(seurat_obj1, resolution = 0.6)  

seurat_obj1 <- RunUMAP(seurat_obj1, dims = 1:20)
DimPlot(seurat_obj1,group.by = "seurat_clusters",label = TRUE,reduction = "umap")
DimPlot(seurat_obj1,group.by = "sample",label = TRUE,reduction = "umap")
DimPlot(seurat_obj1,group.by = "celltype",label = TRUE,reduction = "umap")
DimPlot(seurat_obj1,group.by = "condition",label = TRUE,reduction = "umap")

saveRDS(seurat_obj1, file = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE160269_cd45pos_seurat.rds")

seurat_obj1[["RNA"]] <- as(seurat_obj1[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

#use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  seurat_obj1,
  from = "seurat",
  to = "anndata",
  outFile = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE160269_cd45pos.h5ad"
)

library(SingleCellExperiment)

seurat_obj1[["RNA"]] <- as(seurat_obj1[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seurat_obj1)

if ("pca" %in% names(seurat_obj1@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(seurat_obj1, "pca")
}
if ("umap" %in% names(seurat_obj1@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seurat_obj1, "umap")
}
saveRDS(sce_obj, file = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE160269_cd45pos_sce.rds")




################################################################################ 
############################      GSE166766    #################################

library(Seurat)
library(Matrix)
library(dplyr)

file_paths <- list(
  mock = list(
    barcodes = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082289_mock_barcodes.tsv.gz",
    features = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082289_mock_features.tsv.gz",
    matrix = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082289_mock_matrix.mtx.gz"
  ),
  dpi_1 = list(
    barcodes = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082290_1dpi_barcodes.tsv.gz",
    features = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082290_1dpi_features.tsv.gz",
    matrix = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082290_1dpi_matrix.mtx.gz"
  ),
  dpi_2 = list(
    barcodes = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082291_2dpi_barcodes.tsv.gz",
    features = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082291_2dpi_features.tsv.gz",
    matrix = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082291_2dpi_matrix.mtx.gz"
  ),
  dpi_3 = list(
    barcodes = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082292_3dpi_barcodes.tsv.gz",
    features = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082292_3dpi_features.tsv.gz",
    matrix = "/Users/xiaoying/Downloads/GSE166766_RAW/GSM5082292_3dpi_matrix.mtx.gz"
  )
)

read_data <- function(barcodes_file, features_file, matrix_file) {
  barcodes <- read.delim(barcodes_file, header = FALSE, stringsAsFactors = FALSE)
  features <- read.delim(features_file, header = FALSE, stringsAsFactors = FALSE)
  matrix <- readMM(matrix_file)
  
  rownames(matrix) <- features$V1
  colnames(matrix) <- barcodes$V1
  
  seurat_obj <- CreateSeuratObject(counts = matrix, project = "SingleCell_Project")
  return(seurat_obj)
}

mock_seurat <- read_data(file_paths$mock$barcodes, file_paths$mock$features, file_paths$mock$matrix)
dpi_1_seurat <- read_data(file_paths$dpi_1$barcodes, file_paths$dpi_1$features, file_paths$dpi_1$matrix)
dpi_2_seurat <- read_data(file_paths$dpi_2$barcodes, file_paths$dpi_2$features, file_paths$dpi_2$matrix)
dpi_3_seurat <- read_data(file_paths$dpi_3$barcodes, file_paths$dpi_3$features, file_paths$dpi_3$matrix)

seurat_combined <- merge(mock_seurat, y = c(dpi_1_seurat, dpi_2_seurat, dpi_3_seurat), 
                         add.cell.ids = c("mock", "1dpi", "2dpi", "3dpi"))

seurat_combined$sample <- rep(c("mock", "1dpi", "2dpi", "3dpi"), 
                              times = c(ncol(mock_seurat), ncol(dpi_1_seurat), 
                                        ncol(dpi_2_seurat), ncol(dpi_3_seurat)))


seurat_combined <- subset(seurat_combined, subset = nFeature_RNA > 200)

seurat_combined <- subset(seurat_combined, subset = nCount_RNA > 200)

seurat_combined[["percent.mt"]] <- PercentageFeatureSet(seurat_combined, pattern = "^MT-")
seurat_combined <- subset(seurat_combined, subset = percent.mt < 10)

seurat_combined <- NormalizeData(seurat_combined)

seurat_combined <- ScaleData(seurat_combined, features = rownames(seurat_combined))

seurat_combined <- RunPCA(seurat_combined)


library(harmony)
seurat_combined <- RunHarmony(seurat_combined, group.by.vars = "sample")

seurat_combined <- RunUMAP(seurat_combined, dims = 1:30)

seurat_combined <- FindNeighbors(seurat_combined, dims = 1:30)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)

DimPlot(seurat_combined, reduction = "umap", group.by = "sample")

seurat_combined$sample <- rep(c("mock", "1dpi", "2dpi", "3dpi"), 
                              times = c(ncol(mock_seurat), ncol(dpi_1_seurat), 
                                        ncol(dpi_2_seurat), ncol(dpi_3_seurat)))

seurat_combined$condition <- ifelse(seurat_combined$sample == "mock", "mock", "infected")

library(org.Hs.eg.db)
head(rownames(seurat_combined))
ids=select(org.Hs.eg.db,keys = rownames(seurat_combined),
           columns = c('ENSEMBL','SYMBOL'),
           keytype = 'ENSEMBL')
head(ids)

dim(ids) 
ids=na.omit(ids)
dim(ids) 
length(unique(ids$SYMBOL)) 
ids=ids[!duplicated(ids$SYMBOL),]
ids=ids[!duplicated(ids$ENSEMBL),]
pos=match(ids$ENSEMBL,rownames(seurat_combined) )
seurat_combined=seurat_combined[pos,]
seurat_combined

rownames(seurat_combined) = ids$SYMBOL

genes_of_interest <- c("KRT5", "DAPL1", "TP63", "FOXJ1", "CCDC153", "CCDC113", 
                       "SCGB1A1", "KRT15", "CYP2F2", "LYPD2", "CBR2", "KRT4", 
                       "KRT13", "CHG1", "ASCL1", "POU2F3", "AVIL", "GNAT3", "TRPM5", 
                       "FOXI1", "CFTR", "ASCL3", "MUC5AC", "MUC5B", "GP2", "SPDEF")

cell_type_annotation <- list(
  basal_cells = c(4, 5, 12, 20),
  ciliated_cells = c(0, 3, 7, 16, 19),
  club_cells = c(6, 8, 9, 13, 20, 21, 22, 15),
  BC_club_cells = c(1, 10, 11),
  neuroendocrine_cells = c(18),
  tuft_cells = c(17),
  ionocytes = c(14, 23),
  goblet_cells = c(2)
)

cell_type_vector <- rep(NA, length(colnames(seurat_combined)))  

for (cell_type in names(cell_type_annotation)) {

  clusters <- cell_type_annotation[[cell_type]]
  
  cell_indices <- which(seurat_combined$seurat_clusters %in% clusters)  
  cell_type_vector[cell_indices] <- cell_type  
}

seurat_combined$cell_type <- cell_type_vector

table(seurat_combined$cell_type)  

saveRDS(seurat_obj_25d11, file = "/Users/xiaoying/Desktop/GSE223626_seurat.rds")

seurat_obj_25d11[["RNA"]] <- as(seurat_obj_25d11[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  seurat_obj_25d11,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/"
)

library(SingleCellExperiment)

seurat_obj_25d11[["RNA"]] <- as(seurat_obj_25d11[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seurat_obj_25d11)

if ("pca" %in% names(seurat_obj_25d11@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(seurat_obj_25d11, "pca")
}
if ("umap" %in% names(seurat_obj_25d11@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seurat_obj_25d11, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/_sce.rds")



################################################################################ 
############################      GSE168732   #################################

library(stringr)
data_dir <- "/Users/xiaoying/Downloads/GSE168732_RAW/"

sample_dirs <- list.dirs(main_dir, full.names = TRUE, recursive = FALSE)

rename_to_standard <- function(folder) {
  files <- list.files(folder, full.names = TRUE)
  
  barcodes_file <- files[grepl("barcode", files, ignore.case = TRUE)]
  if (length(barcodes_file) == 1) {
    file.rename(barcodes_file, file.path(folder, "barcodes.tsv.gz"))
  }
    features_file <- files[grepl("feature|gene", files, ignore.case = TRUE)]
  if (length(features_file) == 1) {
    file.rename(features_file, file.path(folder, "features.tsv.gz"))
  }
  
  matrix_file <- files[grepl("matrix", files, ignore.case = TRUE)]
  if (length(matrix_file) == 1) {
    file.rename(matrix_file, file.path(folder, "matrix.mtx.gz"))
  }
}

for (dir in sample_dirs) {
  rename_to_standard(dir)
}


sample_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
sample_names <- basename(sample_dirs)

seurat_list <- list()

for (i in seq_along(sample_dirs)) {
  sample <- sample_names[i]
  sample_path <- sample_dirs[i]
  
  counts <- Read10X(data.dir = sample_path)
  seu <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
  
  seu$sample <- sample
  seu$group <- ifelse(grepl("^H", sample), "Healthy",
                      ifelse(grepl("^P[1-4]", sample), "Phase1", "Phase2"))
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  if (sample == "P1_before") {
    seu <- subset(seu, subset = nCount_RNA > 1000 & nCount_RNA < 60000 & percent.mt < 5)
  } else {
    seu <- subset(seu, subset = nCount_RNA > 2000 & nCount_RNA < 60000 & percent.mt < 5)
  }
  
  seurat_list[[sample]] <- seu
}


integration_samples <- sample_names[grepl("^P[1-4]|^H[1-3]", sample_names)]
integration_list <- seurat_list[integration_samples]

for (i in 1:length(integration_list)) {
  integration_list[[i]] <- NormalizeData(integration_list[[i]])
  integration_list[[i]] <- FindVariableFeatures(integration_list[[i]], selection.method = "vst", nfeatures = 2000)
}

anchors <- FindIntegrationAnchors(object.list = integration_list, dims = 1:30)

integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 50)

ElbowPlot(integrated)

integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.1)
DimPlot(integrated, label = TRUE, group.by = "seurat_clusters")

integrated1 <- FindClusters(integrated, resolution = 1.2)
DimPlot(integrated1, group.by = "seurat_clusters", label = TRUE)

query_samples <- sample_names[grepl("^P[5-7]", sample_names)]
query_list <- seurat_list[query_samples]

for (i in seq_along(query_list)) {
  query <- query_list[[i]]
  query <- NormalizeData(query)
  query <- FindVariableFeatures(query)
  
  transfer_anchors <- FindTransferAnchors(reference = integrated, query = query, dims = 1:30)
  predictions <- TransferData(anchorset = transfer_anchors, refdata = integrated$seurat_clusters, dims = 1:30)
  
  query <- AddMetaData(query, metadata = predictions)
  query_list[[i]] <- query
}


library(SingleR)
ref <- HumanPrimaryCellAtlasData()

data_expr <- GetAssayData(integrated[["integrated"]], layer = "data")
data_expr_avg <- AverageExpression(integrated, assays = "RNA", slot = "data")$RNA

singleR_results <- SingleR(test = data_expr_avg, ref = ref, labels = ref$label.main)

clusters <- integrated$seurat_clusters

singleR_results <- SingleR(test = data_expr, ref = ref, labels = ref$label.main, clusters = clusters)

integrated$SingleR_label <- singleR_results$labels[match(clusters, rownames(singleR_results))]

library(ggplot2)
DimPlot(integrated, group.by = "SingleR_label", label = TRUE, repel = TRUE) +
  ggtitle("UMAP of Cells Colored by SingleR Cell Type")

DimPlot(integrated,
                 group.by = "SingleR_label",
                 split.by = "group",
                 label = TRUE,
                 repel = TRUE,
                 ncol = 3) +
       ggtitle("UMAP by SingleR Cell Type Split by Condition")

library(dplyr)
library(ggplot2)

integrated$condition <- ifelse(integrated$sample %in% c("H1", "H2", "H3"), "Healthy",
                        ifelse(grepl("_before$", integrated$sample), "Before",
                        ifelse(grepl("_after$", integrated$sample), "After", NA)))


b_cells <- subset(integrated, subset = SingleR_label == "B_cell")
DefaultAssay(b_cells) <- "RNA"
b_cells <- NormalizeData(b_cells)
b_cells <- FindVariableFeatures(b_cells, selection.method = "vst", nfeatures = 2000)
b_cells <- ScaleData(b_cells)
b_cells <- RunPCA(b_cells, npcs = 30)
ElbowPlot(b_cells)
b_cells <- RunUMAP(b_cells, dims = 1:15)
b_cells <- FindNeighbors(b_cells, dims = 1:15)
b_cells <- FindClusters(b_cells, resolution = 0.1)  # resolution可调整

DimPlot(b_cells, label = TRUE)

FeaturePlot(b_cells, features = c("MS4A1", "TCL1A", "CD27", "IGHA1", "IGHG1", "CD38"))


b_markers <- c("CD19", "CD79A", "MS4A1",     
               "TCL1A", "IGHM", "IGHD",   
               "CD27", "CD38", "IGHA1", "IGHG1")  

b_markers_factor <- factor(b_markers, levels = b_markers)

DotPlot(b_cells, features = b_markers_factor) +
  RotatedAxis() +
  labs(title = "B cell marker gene expression (cluster level)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

monocyte_cells <- subset(integrated, subset = SingleR_label == "Monocyte")
DefaultAssay(monocyte_cells) <- "RNA"
monocyte_cells <- NormalizeData(monocyte_cells)
monocyte_cells <- FindVariableFeatures(monocyte_cells, selection.method = "vst", nfeatures = 2000)
monocyte_cells <- ScaleData(monocyte_cells)
monocyte_cells <- RunPCA(monocyte_cells, npcs = 30)
ElbowPlot(monocyte_cells)
monocyte_cells <- RunUMAP(monocyte_cells, dims = 1:15)
monocyte_cells <- FindNeighbors(monocyte_cells, dims = 1:15)
monocyte_cells <- FindClusters(monocyte_cells, resolution = 0.1)  

DimPlot(monocyte_cells, label = TRUE)

FeaturePlot(monocytes, features = c("CD14", "FCGR3A"), cols = c("lightgrey", "darkred"),
            reduction = "umap", pt.size = 0.3)


FeaturePlot(monocytes, features = "CD14", cols = c("lightgrey", "blue"), pt.size = 0.3)
FeaturePlot(monocytes, features = "FCGR3A", cols = c("lightgrey", "green"), pt.size = 0.3)


mono_markers <- factor(c("CD14", "FCGR3A"), levels = c("CD14", "FCGR3A"))

DotPlot(monocytes, features = mono_markers) +
  RotatedAxis() +
  labs(title = "Monocyte marker expression (by cluster)")



monocyte_cells$monocyte_type <- NA  

monocyte_cells$monocyte_type[monocyte_cells$seurat_clusters %in% c(0, 2, 3)] <- "CD14_monocyte"
monocyte_cells$monocyte_type[monocyte_cells$seurat_clusters == 1] <- "CD16_monocyte"
monocyte_cells$monocyte_type[monocyte_cells$seurat_clusters == 4] <- "Unknown"

monocyte_cells$monocyte_type <- factor(monocyte_cells$monocyte_type, levels = c("CD14_monocyte", "CD16_monocyte", "Unknown"))

table(monocyte_cells$seurat_clusters, monocyte_cells$monocyte_type)

DimPlot(monocyte_cells, group.by = "monocyte_type", label = TRUE, pt.size = 0.5) +
  ggtitle("Monocyte Subtypes (CD14+ vs CD16+ vs Unknown)")

monocytes_known <- subset(monocyte_cells, subset = monocyte_type != "Unknown")

mono_counts <- table(monocytes_known$sample, monocytes_known$monocyte_type)

print(mono_counts)

mono_props <- prop.table(mono_counts, margin = 1)

round(mono_props * 100, 1)




################################################################################ 
############################      GSE169147  #################################

library(Seurat)

data_dir <- "/Users/xiaoying/Downloads/GSE169147_RAW /"
sample_files <- c("GSM5176926_counts_scSAR1.txt.gz", "GSM5176927_counts_scSAR2.txt.gz", "GSM5176928_counts_scSAR3.txt.gz", 
                  "GSM5176929_counts_scHealthy1.txt.gz", "GSM5176930_counts_scHealthy2.txt.gz", "GSM5176931_counts_scHealthy3.txt.gz")
sample_names <- c("SAR1", "SAR2", "SAR3", "Healthy1", "Healthy2", "Healthy3")

seurat_list <- list()

for (i in seq_along(sample_files)) {
  sample_file <- sample_files[i]
  sample_name <- sample_names[i]
  
  cat("Processing sample:", sample_name, "\n")
  
  counts <- read.table(file.path(data_dir, sample_file), header = TRUE, row.names = 1, sep = "\t")
  
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
  
  seurat_obj$sample <- sample_name
  
  seurat_list[[sample_name]] <- seurat_obj
  
  cat("Completed sample:", sample_name, "\n")
}

if (length(seurat_list) == 0) {
  stop("No samples were successfully read, please check file paths and formats")
}

if (length(seurat_list) == 1) {
  combined_seurat <- seurat_list[[1]]
} else {
  cat("Merging all samples...\n")
  
  combined_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], 
                           add.cell.ids = sample_names, 
                           project = "combined_project")
}

cat("\nMerged Seurat object information:\n")
print(combined_seurat)

combined_seurat <- subset(combined_seurat, subset = nFeature_RNA > 100)
combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

combined_seurat <- FindVariableFeatures(combined_seurat)

combined_seurat <- ScaleData(combined_seurat)

combined_seurat <- RunPCA(combined_seurat, npcs = 30)
combined_seurat <- RunUMAP(combined_seurat,dims = 1:30)

ElbowPlot(combined_seurat)

combined_seurat <- FindNeighbors(combined_seurat, dims = 1:20)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.8)

DimPlot(combined_seurat, reduction = "pca", group.by = "seurat_clusters")

DimPlot(combined_seurat, reduction = "umap", group.by = "seurat_clusters")

interest_marker <- c("CD45","CD3E","PTPRC","MKI67","KLRC1","COL1A1","LYVE1",
                     "CD20","MS4A1","TPSAB1","Tryptase","CD68","ACTA2","VWF","DCT","DSG1")

library(ggplot2)

dotplot <- DotPlot(combined_seurat, features = interest_marker, group.by = "seurat_clusters") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Gene Expression across Cell Types")

dotplot

DimPlot(combined_seurat, reduction = "umap", group.by = "seurat_clusters",label=TRUE)

Myeloid <-c(8,13,24)
B_cell <- c(26)
Mast_cell <- c(22)
Fibroblasts <- c(17,18,20,6)
Endothelium <- c(11,12,25,1)
Lymphatic_endo <- c(14)
Myocytes <- c(10)
Melanocytes <- c(7)
Keratinocytes <- c(21)

cell_type_annotation <- list(
  T_cells = c(0,2,3,4,5,9,16,19),
  Tcell_mitotic =c(23),
  NK = c(15),
  Myeloid = c(8,13,24),
  B_cell = c(26),
  Mast_cell = c(22),
  Fibroblasts = c(17,18,20,6),
  Endothelium = c(11,12,25,1),
  Lymphatic_endo = c(14),
  Myocytes = c(10),
  Melanocytes = c(7),
  Keratinocytes = c(21)
)

cell_type_vector <- rep(NA, length(colnames(combined_seurat)))

for (celltype in names(cell_type_annotation)) {
  clusters <- cell_type_annotation[[celltype]]
  
  cell_indices <- which(combined_seurat$seurat_clusters %in% clusters)
  cell_type_vector[cell_indices] <- celltype
}

combined_seurat$celltype <- cell_type_vector

table(combined_seurat$celltype)

combined_seurat$condition <- ifelse(combined_seurat$sample %in% c("Healthy1", "Healthy2", "Healthy3"), "Healthy", "SAR")

table(combined_seurat$condition,combined_seurat$sample)

saveRDS(combined_seurat, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE169147_seurat.rds")

combined_seurat[["RNA"]] <- as(combined_seurat[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  combined_seurat,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE169147.h5ad"
)

library(SingleCellExperiment)

combined_seurat[["RNA"]] <- as(combined_seurat[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(combined_seurat)

if ("pca" %in% names(combined_seurat@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(combined_seurat, "pca")
}
if ("umap" %in% names(combined_seurat@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(combined_seurat, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE169147_sce.rds")


################################################################################ 
############################      GSE176078  ###################################
library(R.utils)

data_dir <- "/media/cqnu/XiaoYingxue/DA_datasets/raw_data/10x_datasets/GSE176078_Wu_etal_2021_BRCA_scRNASeq"

gzip(file.path(data_dir, "matrix.mtx"), overwrite = TRUE)
gzip(file.path(data_dir, "barcodes.tsv"), overwrite = TRUE)
gzip(file.path(data_dir, "features.tsv"), overwrite = TRUE)

library(Seurat)
data <- Read10X(data.dir = data_dir, gene.column = 1)

seurat_obj <- CreateSeuratObject(counts = data)

metadata <- read.csv("/media/cqnu/XiaoYingxue/DA_datasets/raw_data/10x_datasets/GSE176078_Wu_etal_2021_BRCA_scRNASeq/metadata.csv", stringsAsFactors = FALSE)
head(colnames(seurat_obj))
head(metadata$X) 

rownames(metadata) <- metadata$X

common_cells <- intersect(colnames(seurat_obj), rownames(metadata))

seurat_obj <- AddMetaData(seurat_obj, metadata = metadata[common_cells, ])

library(Seurat)

seurat_obj <- NormalizeData(seurat_obj)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

seurat_obj <- ScaleData(seurat_obj)

seurat_obj <- RunPCA(seurat_obj, npcs = 30)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

seurat_obj$condition <- seurat_obj$subtype
seurat_obj$sample <- seurat_obj$orig.ident

saveRDS(seurat_obj,file="/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE176078.rds")
saveRDS(seurat_obj,file="/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasets/GSE176078_seurat.rds")
seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seurat_obj)


################################################################################ 
############################      GSE184880  ###################################

ibrary(R.utils)

raw_dir <- "./GSE184880"
setwd(raw_dir)

all_files <- list.files(raw_dir, pattern = "\\.tsv\\.gz$", full.names = FALSE)

sample_groups <- unique(sub(".*_(.+?)\\..*", "\\1", all_files))

for (group in sample_groups) {
  group_files <- list.files(raw_dir, pattern = paste0("_", group, "\\."), full.names = TRUE)
  
  group_dir <- file.path(raw_dir, group)
  if (!dir.exists(group_dir)) {
    dir.create(group_dir)
  }
  
  for (f in group_files) {
    file.rename(from = f, to = file.path(group_dir, basename(f)))
  }
}

sample_dirs <- list.dirs(raw_dir, recursive = FALSE, full.names = TRUE)

for (sample_dir in sample_dirs) {
  sample_files <- list.files(sample_dir, pattern = "\\.gz$", full.names = TRUE)
  
  for (file_path in sample_files) {
    new_name <- ""
    
    if (grepl("barcodes", basename(file_path))) {
      new_name <- "barcodes.tsv.gz"
    } else if (grepl("genes", basename(file_path))) {
      new_name <- "features.tsv.gz"
    } else if (grepl("matrix", basename(file_path))) {
      new_name <- "matrix.mtx.gz"
    }
    
    if (new_name != "") {
      file.rename(file_path, file.path(sample_dir, new_name))
    }
  }
}

library(Seurat)

sample_dirs <- list.dirs(raw_dir, recursive = FALSE, full.names = FALSE)

seurat_list <- lapply(sample_dirs, function(sample) {
  data_path <- file.path(raw_dir, sample)
  counts <- Read10X(data.dir = data_path)
  
  obj <- CreateSeuratObject(counts = counts, project = sample)
  
  obj$sample <- sample
  
  return(obj)
})

names(seurat_list) <- sample_dirs

seurat_list

library(Seurat)
library(scran)
library(scater)
library(batchelor)
library(GSVA)
library(msigdbr)
library(limma)
library(patchwork)

seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 40)
  
  mito_genes <- grep("^MT-", rownames(obj), value = TRUE)
  obj <- obj[!rownames(obj) %in% mito_genes, ]
  
  return(obj)
})

for (i in seq_along(seurat_list)) {
  obj <- seurat_list[[i]]
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30)
  obj <- RunUMAP(obj, dims = 1:10)
  
  seurat_list[[i]] <- obj
}

anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:10)
combined <- IntegrateData(anchorset = anchors, dims = 1:10)

DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, dims = 1:10)
combined <- FindNeighbors(combined, dims = 1:10)
combined <- FindClusters(combined, resolution = 0.5)

markers <- FindAllMarkers(combined,
                          only.pos = TRUE,
                          min.pct = 0.1,
                          logfc.threshold = 0.25)

marker_list <- list(
  T_cell           = c("CD3D", "CD3E", "CD8A"),
  Epithelia        = c("KRT18", "EPCAM", "CD24", "KRT19"),
  Monocytic        = c("CD14", "C1QA"),
  Endothelia       = c("PECAM1", "CLDN5"),
  Cell_cycle       = c("MKI67", "TOP2A"),
  Fibroblast       = c("DCN", "OGN"),
  B_cell_plasma    = c("CD79A", "JCHAIN"),
  SMC_myoFibroblast = c("ACTA2", "MYH11", "TAGLN")
)

all_markers <- unique(unlist(marker_list))

available_markers <- all_markers[all_markers %in% rownames(combined)]

DotPlot(combined, features = available_markers, group.by = "seurat_clusters") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marker Expression per Cluster")

cluster_annotation <- list(
  T_cell            = c(2,3,5,6),     
  Epithelia         = c(1,4,9,11),
  Monocytic         = c(7,10,17),
  Endothelia        = c(12,18,20),
  Cell_cycle        = c(13,16),
  Fibroblast        = c(0,8),
  B_cell_plasma     = c(15),
  SMC_myoFibroblast = c(14,19)
)

cluster_ids <- as.numeric(as.character(combined$seurat_clusters))

manual_anno <- sapply(cluster_ids, function(clu) {
  matched <- names(Filter(function(x) clu %in% x, cluster_annotation))
  if (length(matched) == 1) {
    return(matched)
  } else {
    return("Unknown")
  }
})

combined$manual_annotation <- manual_anno

DimPlot(combined, group.by = "manual_annotation", label = TRUE, repel = TRUE) +
  ggtitle("Manual Annotation with Marker-Based Cell Types")

DimPlot(combined, group.by = "sample", label = TRUE, repel = TRUE) 

library(dplyr)

meta <- combined@meta.data

meta$condition <- ifelse(meta$sample %in% c("Norm1", "Norm2", "Norm3", "Norm4", "Norm5"),
                         "norm", "cancer")

meta$sub_stage <- case_when(
  meta$sample %in% c("Norm1", "Norm2", "Norm3", "Norm4", "Norm5") ~ "norm",
  meta$sample %in% c("Cancer3", "Cancer4", "Cancer7") ~ "stage1",
  meta$sample %in% c("Cancer2", "Cancer5") ~ "stage2",
  meta$sample %in% c("Cancer1", "Cancer6") ~ "stage3",
  TRUE ~ "Unknown"
)

combined@meta.data <- meta

saveRDS(combined, file = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/10X_processed_datasers/GSE184880_obj_annotation.rds")




################################################################################ 
############################      GSE201765  ###################################
library(Seurat)
library(dplyr)

data_dir <- "/GSE201765"
sample_names <- c("Relapse_1", "Relapse_2", "Remission_1", "Remission_2")

seurat_list <- lapply(sample_names, function(sample) {
  sample_path <- file.path(data_dir, sample)
  data <- Read10X(sample_path)
  seurat_obj <- CreateSeuratObject(counts = data, project = sample, min.cells = 3, min.features = 200)
  seurat_obj$sample <- sample
  return(seurat_obj)
})

combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_names)

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")

VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

combined <- subset(combined, subset = nCount_RNA > 500 & nCount_RNA < 25000 & percent.mt < 5)

combined <- NormalizeData(combined, normalization.method = "RC", scale.factor = 1e6)
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 1e6)

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

combined <- RunUMAP(combined, dims = 1:15)
combined <- FindNeighbors(combined, dims = 1:15)
combined <- FindClusters(combined, resolution = 0.5)

DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

FeaturePlot(combined, features = c("Cd3e", "Cd3d", "Cd8a", "Cd4", "Nkg7", "Klra1", "Xcl1", "Ifng", "Gzmb"))

library(SingleR)
library(celldex)
library(Seurat)
library(SeuratObject)
library(SummarizedExperiment)

data_matrix <- GetAssayData(combined, slot = "data")

cell_metadata <- combined@meta.data

ref <- celldex::ImmGenData()

singleR_results <- SingleR(
  test = data_matrix,
  ref = ref,
  labels = ref$label.main
)

combined1 <-combined
combined1$SingleR_label <- singleR_results$labels

DimPlot(combined1, group.by = "SingleR_label", label = TRUE, repel = TRUE) + NoLegend()

tcell_clusters <- c("5", "8", "9")

table(combined$seurat_clusters)

tcell_obj <- subset(combined, idents = tcell_clusters)

tcell_obj <- NormalizeData(tcell_obj)
tcell_obj <- FindVariableFeatures(tcell_obj)
tcell_obj <- ScaleData(tcell_obj)
tcell_obj <- RunPCA(tcell_obj)
tcell_obj <- RunUMAP(tcell_obj, dims = 1:15)
tcell_obj <- FindNeighbors(tcell_obj, dims = 1:15)
tcell_obj <- FindClusters(tcell_obj, resolution = 0.5)

DimPlot(tcell_obj, reduction = "umap", label = TRUE)

table(tcell_obj$sample)

tcell_obj$condition <- ifelse(tcell_obj$sample %in% c("Relapse_1", "Relapse_2"), "Relapse", "Remission")

table(tcell_obj$condition)

FeaturePlot(tcell_obj, features = "Cd4", split.by = "condition")

FeaturePlot(tcell_obj, features = c("Il7r", "Ccr7", "Sell"), split.by = "condition")

FeaturePlot(tcell_obj, features = c("Gzma", "Gzmb", "Pdcd1", "Ifng"), split.by = "condition")

FeaturePlot(tcell_obj, features = c("Klra1", "Xcl1"), split.by = "condition")

relapse_obj <- subset(tcell_obj, subset = condition == "Relapse")

genes.use <- c("Cd4", "Il7r", "Ccr7", "Sell", 
               "Gzma", "Gzmb", "Pdcd1", "Ifng", 
               "Klra1", "Xcl1")

DotPlot(relapse_obj, features = genes.use, group.by = "seurat_clusters", dot.scale = 8) +
  RotatedAxis() +
  ggtitle("Relapse 条件下 T细胞亚群 marker 表达（按cluster分组）") +
  theme(plot.title = element_text(hjust = 0.5))

library(dplyr)

tcell_obj$celltype <- recode(tcell_obj$seurat_clusters,
                             "0" = "CD8_T1",
                             "1" = "CD4_T",
                             "2" = "CD8_T2",
                             "3" = "CD8_T1",
                             "4" = "CD8_T1",
                             "5" = "NKT")

table(tcell_obj$seurat_clusters, tcell_obj$celltype)

library(dplyr)
cell_counts <- as.data.frame(table(tcell_obj$condition, tcell_obj$celltype))
colnames(cell_counts) <- c("condition", "celltype", "count")

saveRDS(tcell_obj, file = "/Volumes/XiaoYingxue/DA_datasets/10x_datasets/GSE201765_tcell_obj_annotated.rds")




