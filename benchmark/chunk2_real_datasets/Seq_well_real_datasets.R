##Seq-Well_real_datasets:

########################################################################
##################--------         GSE116256     ---------##############

library(data.table)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scater)

data_dir <- "/Users/xiaoying/Downloads/GSE116256_RAW/"

read_dem_file <- function(dem_path){
  dt <- tryCatch(fread(dem_path, data.table = FALSE), error = function(e) NULL)
  if(is.null(dt) || nrow(dt) == 0) stop("error: ", dem_path)
  
  first_col_name <- colnames(dt)[1]
  if(first_col_name %in% c("gene","GENE","V1","Gene")){
    genes <- dt[[1]]
    mat <- dt[,-1, drop=FALSE]
  } else {
    mat <- dt
    genes <- rownames(mat)
  }
  
  mat <- as.matrix(mat)
  
  mode(mat) <- "numeric"
  mat <- round(mat)
  storage.mode(mat) <- "integer"
  
  if(!inherits(mat, "Matrix")) mat <- Matrix(mat, sparse = TRUE)
  
  if(!is.null(genes)) rownames(mat) <- genes
  
  if(any(grepl("^V[0-9]+$", colnames(mat)))) colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
  
  return(mat)
}

samples_to_use <- c(
  "GSM3587996_BM1.dem.txt.gz","GSM3587997_BM2.dem.txt.gz",
  "GSM3587998_BM3.dem.txt.gz","GSM3588000_BM4.dem.txt.gz",
  "GSM3587927_AML314-D0.dem.txt.gz","GSM3587929_AML314-D31.dem.txt.gz",
  "GSM3587946_AML371-D0.dem.txt.gz","GSM3587948_AML371-D34.dem.txt.gz",
  "GSM3587959_AML475-D0.dem.txt.gz","GSM3587961_AML475-D29.dem.txt.gz",
  "GSM3587992_AML997-D0.dem.txt.gz","GSM3587994_AML997-D35.dem.txt.gz",
  "GSM3587953_AML420B-D0.dem.txt.gz","GSM3587957_AML420B-D35.dem.txt.gz",
  "GSM3587963_AML556-D0.dem.txt.gz","GSM3587967_AML556-D31.dem.txt.gz",
  "GSM3587931_AML328-D0.dem.txt.gz","GSM3587937_AML328-D29.dem.txt.gz"
)

seu_list <- list()
for(p in samples_to_use){
  dem_path <- file.path(data_dir, p)
  if(!file.exists(dem_path)) stop("error: ", dem_path)
  
  mat <- read_dem_file(dem_path)
  if(nrow(mat) < ncol(mat)){ mat <- t(mat) }
  
  sample_id <- gsub("\\.dem\\.txt\\.gz$", "", p)
  seu <- CreateSeuratObject(counts = mat, project = sample_id,
                            min.cells = 3, min.features = 200)
  
  mito_genes <- grep("^MT-", rownames(seu), value = TRUE)
  ribo_genes <- grep("^RPS|^RPL", rownames(seu), value = TRUE)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito_genes)
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, features = ribo_genes)
  seu$sample_id <- sample_id
  
  seu_list[[sample_id]] <- seu
  message("complete: ", sample_id, " | cell: ", ncol(seu))
}

sapply(seu_list, ncol)
sapply(seu_list, nrow)

anno_files <- list.files(data_dir, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE)

for(p in samples_to_use){
  core_name <- gsub(".*_(AML.*|BM[0-9]+)-?.*\\.dem\\.txt\\.gz$", "\\1", p)
  sample_id <- gsub("\\.dem\\.txt\\.gz$", "", p)
  
  matched_anno <- anno_files[grepl(core_name, basename(anno_files))]
  
  if(length(matched_anno) == 1){
    adf <- fread(matched_anno, data.table = FALSE)
    cb_col <- colnames(adf)[1]
    rownames(adf) <- adf[[cb_col]]
    
    common_bcs <- intersect(colnames(seu_list[[sample_id]]), rownames(adf))
    
    if(length(common_bcs) > 0){
      seu_list[[sample_id]] <- AddMetaData(seu_list[[sample_id]], metadata = adf[common_bcs, , drop = FALSE])
      message("complete: ", sample_id, " | cell: ", length(common_bcs))
    } else {
      message("⚠️ : ", sample_id)
    }
  } else {
    message("⚠️ ", sample_id)
  }
}

seu_merged <- Reduce(function(x, y) merge(x, y), seu_list)

sample_map <- c(
  "GSM3587927_AML314-D0" = "AML314-D0",
  "GSM3587929_AML314-D31" = "AML314-D31",
  "GSM3587931_AML328-D0" = "AML328-D0",
  "GSM3587937_AML328-D29" = "AML328-D29",
  "GSM3587946_AML371-D0" = "AML371-D0",
  "GSM3587948_AML371-D34" = "AML371-D34",
  "GSM3587953_AML420B-D0" = "AML420B-D0",
  "GSM3587957_AML420B-D35" = "AML420B-D35",
  "GSM3587959_AML475-D0" = "AML475-D0",
  "GSM3587961_AML475-D29" = "AML475-D29",
  "GSM3587963_AML556-D0" = "AML556-D0",
  "GSM3587967_AML556-D31" = "AML556-D31",
  "GSM3587992_AML997-D0" = "AML997-D0",
  "GSM3587994_AML997-D35" = "AML997-D35",
  "GSM3587996_BM1" = "BM1",
  "GSM3587997_BM2" = "BM2",
  "GSM3587998_BM3" = "BM3",
  "GSM3588000_BM4" = "BM4"
)

seu_merged@meta.data$sample <- sample_map[seu_merged@meta.data$sample_id]

table(seu_merged@meta.data$sample)

seu_merged@meta.data$sample_id <- NULL

seu_merged$celltype <- seu_merged$CellType
seu_merged@meta.data$CellType<- NULL

seu_merged@meta.data$condition <- ifelse(
  seu_merged@meta.data$sample %in% c("BM1", "BM2", "BM3", "BM4"),
  "healthy",
  ifelse(grepl("-D0$", seu_merged@meta.data$sample), "D0", "D30")
)

table(seu_merged@meta.data$condition)

drop_cols <- c( "AlignedToGenome", "AlignedToTranscriptome", 
               "TranscriptomeUMIs", "CyclingScore", "CyclingBinary","MutTranscripts" ,
               "WtTranscripts", "PredictionRF2","Score_HSC","Score_Prog", "Score_GMP",             
               "Score_ProMono", "Score_Mono", "Score_cDC", 
               "Score_pDC", "Score_earlyEry" , "Score_lateEry",         
               "Score_ProB",  "Score_B" ,"Score_Plasma" ,         
               "Score_T" ,  "Score_CTL" ,"Score_NK")

seu_merged@meta.data <- seu_merged@meta.data[, !(colnames(seu_merged@meta.data) %in% drop_cols)]

qc_min_umis <- 1000
qc_min_genes <- 500
qc_max_mito  <- 20

seu_merged <- subset(
  seu_merged,
  subset = nCount_RNA >= qc_min_umis &
    nFeature_RNA >= qc_min_genes &
    percent.mt <= qc_max_mito
)

seu_merged <- NormalizeData(seu_merged, normalization.method = "LogNormalize", scale.factor = 1e4)

seu_merged <- FindVariableFeatures(seu_merged, selection.method = "vst", nfeatures = 2000)

seu_merged <- ScaleData(seu_merged, verbose = FALSE)

seu_merged <- RunPCA(seu_merged,  npcs = 50, verbose = FALSE)

seu_merged <- FindNeighbors(seu_merged, dims = 1:30)

seu_merged <- FindClusters(seu_merged, resolution = 0.5)

seu_merged <- RunUMAP(seu_merged, dims = 1:30)

DimPlot(seu_merged, reduction = "umap", group.by = "sample", pt.size = 0.5) + ggtitle("UMAP by sample")
DimPlot(seu_merged, reduction = "umap", group.by = "condition", label = TRUE)
DimPlot(seu_merged, reduction = "umap", group.by = "celltype", label = TRUE)













########################################################################
##################--------         GSE200151     ---------##############


library(Seurat)
library(Matrix)
library(dplyr)
library(stringr)

data_dir <- "/Users/xiaoying/Downloads/GSE200151_RAW"

tenweek_barcodes <- read.table(file.path(data_dir, "GSE200151_10Week_countsbarcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
tenweek_features <- read.table(file.path(data_dir, "GSE200151_10Week_countsfeatures.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
tenweek_matrix <- readMM(file.path(data_dir, "GSE200151_10Week_countsmatrix.mtx.gz"))

tenweek_metadata <- read.delim(file.path(data_dir, "GSE200151_Updated10wk_alexandria_structured_metadata10.txt.gz"), 
                                header = FALSE, stringsAsFactors = FALSE)
colnames(tenweek_metadata) <- tenweek_metadata[1, ]
tenweek_metadata <- tenweek_metadata[-1, ]
tenweek_metadata <- tenweek_metadata[-1, ]

tenweek_seurat <- CreateSeuratObject(counts = tenweek_matrix, 
                                      assay = "RNA", 
                                      project = "GSE200151_10Week", 
                                      meta.data = tenweek_metadata)

library(dplyr)
tenweek_seurat@meta.data <- tenweek_seurat@meta.data %>%
  select(orig.ident, nCount_RNA, 
         nFeature_RNA, CellID, 
         GenericFinal, disease, 
         donor_id, cell_type__ontology_label, 
         GranulomaBurden, GranulomaTiming)

library(stringr)
tenweek_seurat@meta.data$sample <- str_extract(tenweek_seurat@meta.data$CellID, "^Array\\d+_\\d+")

head(tenweek_seurat@meta.data$sample)

tenweek_seurat <- subset(tenweek_seurat, subset = nFeature_RNA > 500 & nCount_RNA > 750)

tenweek_seurat <- NormalizeData(tenweek_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

tenweek_seurat <- ScaleData(tenweek_seurat, features = rownames(tenweek_seurat))

tenweek_seurat <- FindVariableFeatures(tenweek_seurat, selection.method = "vst", nfeatures = 2000)

tenweek_seurat <- RunPCA(tenweek_seurat, features = VariableFeatures(object = tenweek_seurat))

tenweek_seurat <- FindNeighbors(tenweek_seurat, dims = 1:20)
tenweek_seurat <- FindClusters(tenweek_seurat, resolution = 1.00)

tenweek_seurat <- RunUMAP(tenweek_seurat, dims = 1:20)

DimPlot(tenweek_seurat, reduction = "umap", label = TRUE)

DimPlot(tenweek_seurat, reduction = "umap", group.by="celltype",label = TRUE)

saveRDS(tenweek_seurat, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE200151_tenweek_seurat.rds")

tenweek_seurat[["RNA"]] <- as(tenweek_seurat[["RNA"]], "Assay")

library(sceasy)
library(reticulate)
use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  tenweek_seurat,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE200151_tenweek.h5ad"
)

library(SingleCellExperiment)

tenweek_seurat[["RNA"]] <- as(tenweek_seurat[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(tenweek_seurat)

if ("pca" %in% names(tenweek_seurat@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(tenweek_seurat, "pca")
}
if ("umap" %in% names(tenweek_seurat@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(tenweek_seurat, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE200151_tenweek_sce.rds")




########################################################################
##################--------         GSE211191     ---------##############



library(Seurat)
mm <- readRDS("/Users/xiaoying/Downloads/GSE211191_integrated_object.rds")
tnk <- readRDS("/Users/xiaoying/Downloads/GSE211191_tnk_object.rds")
myeloid <- readRDS("/Users/xiaoying/Downloads/GSE211191_myeloid_object.rds")

subset_13_25 <- subset(mm, subset = time %in% c("13 weeks", "25 weeks"))
subset_13_25@meta.data[["orig.ident"]] <- NULL
subset_13_25@meta.data[["unique_identifiers"]] <- NULL
subset_13_25@meta.data[["vaccine_route_unique"]] <- NULL

columns_to_remove <- c(
  "log_dose_vaccine", "pitt_nhp_id", "log_total_cfu", "log_lung_cfu", "log_ln_cfu", "has_scrnaseq",
  "total_cfu_estimated", "mtb_dose", "classification", "dose_group", "date_infected", "date_nx",
  "total_cfu", "lung_cfu", "lymph_cfu", "num_grans", "fdg", "gps", "eps", "burden_group", "days_infected",
  "gender", "age_at_vaccination", "nhp_cohort", "born", "holding", "biosample_id", "barcode", "vrc_nhp_id",
  "percent_ribo", "percent_hsp", "percent_hgb", "binary_group", "tertiary_group"
)

subset_13_25@meta.data <- subset_13_25@meta.data[, !colnames(subset_13_25@meta.data) %in% columns_to_remove]

colnames(subset_13_25@meta.data)

subset_13_25@meta.data[["route"]] <- NULL
subset_13_25@meta.data[["dose"]] <- NULL
subset_13_25@meta.data[["tree.ident"]] <- NULL
subset_13_25@meta.data[["passing_qc"]] <- NULL
subset_13_25@meta.data[["percent_soup"]] <- NULL
subset_13_25@meta.data[["consensus_celltype"]] <- NULL

week_13_subset <- subset(subset_13_25, subset = week == "w13")
week_25_subset <- subset(subset_13_25, subset = week == "w25")

table(week_13_subset@meta.data[["approx_group"]])
table(week_25_subset@meta.data[["approx_group"]])

week_13_subset <- NormalizeData(week_13_subset, normalization.method = "LogNormalize", scale.factor = 10000)
week_13_subset <- FindVariableFeatures(week_13_subset, selection.method = "vst", nfeatures = 2000)
week_13_subset <- ScaleData(week_13_subset, features = rownames(week_13_subset))
week_13_subset <- RunPCA(week_13_subset, features = VariableFeatures(week_13_subset))
week_13_subset <- RunUMAP(week_13_subset, dims = 1:20)
week_13_subset <- FindNeighbors(week_13_subset, dims = 1:20)
week_13_subset <- FindClusters(week_13_subset, resolution = 0.5)
DimPlot(week_13_subset, reduction = "umap", group.by = "seurat_clusters")
DimPlot(week_13_subset, reduction = "umap", group.by = "integrated_celltypel2")

week_13_subset@meta.data[["celltype"]] <- NULL
week_13_subset@meta.data[["array"]] <- NULL
week_13_subset@meta.data[["cluster"]] <- NULL
week_13_subset@meta.data[["study"]] <- NULL
week_13_subset@meta.data[["log_dose"]] <- NULL
week_13_subset@meta.data[["integrated_celltype"]] <- NULL
week_13_subset@meta.data[["integrated_celltypel2"]] <- NULL
week_13_subset@meta.data[["quarternary_group"]] <- NULL
week_13_subset@meta.data[["harmony_leiden"]] <- NULL

week_13_subset$sample <- week_13_subset$vaccine_route_unique_group
week_13_subset$condition <- week_13_subset$approx_group
week_13_subset$celltype <- week_13_subset$integrated_celltypel1

meta <- week_13_subset@meta.data

saveRDS(week_13_subset, "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE211191_13week_seurat.rds")

library(Seurat)
library(sceasy)
library(reticulate)
use_condaenv("base", required = TRUE)
week_13_subset[["RNA"]] <- as(week_13_subset[["RNA"]], "Assay")
sceasy::convertFormat(
  week_13_subset,
  from = "seurat",
  to = "anndata",
  outFile = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE211191_13week.h5ad"
)

sce_obj <- as.SingleCellExperiment(week_13_subset)
library(SingleCellExperiment)
if ("pca" %in% names(week_13_subset@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(week_13_subset, "pca")
}
if ("umap" %in% names(week_13_subset@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(week_13_subset, "umap")
}
saveRDS(sce_obj, file = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE211191_13week_sce.rds")

week_25_subset <- NormalizeData(week_25_subset, normalization.method = "LogNormalize", scale.factor = 10000)
week_25_subset <- FindVariableFeatures(week_25_subset, selection.method = "vst", nfeatures = 2000)
week_25_subset <- ScaleData(week_25_subset, features = rownames(week_25_subset))
week_25_subset <- RunPCA(week_25_subset, features = VariableFeatures(week_25_subset))
week_25_subset <- RunUMAP(week_25_subset, dims = 1:20)
week_25_subset <- FindNeighbors(week_25_subset, dims = 1:20)
week_25_subset <- FindClusters(week_25_subset, resolution = 0.5)
DimPlot(week_25_subset, reduction = "umap", group.by = "seurat_clusters")
DimPlot(week_25_subset, reduction = "umap", group.by = "integrated_celltypel1")

week_25_subset@meta.data[["celltype"]] <- NULL
week_25_subset@meta.data[["array"]] <- NULL
week_25_subset@meta.data[["cluster"]] <- NULL
week_25_subset@meta.data[["study"]] <- NULL
week_25_subset@meta.data[["log_dose"]] <- NULL
week_25_subset@meta.data[["integrated_celltype"]] <- NULL
week_25_subset@meta.data[["integrated_celltypel2"]] <- NULL
week_25_subset@meta.data[["quarternary_group"]] <- NULL
week_25_subset@meta.data[["harmony_leiden"]] <- NULL

week_25_subset$sample <- week_25_subset$vaccine_route_unique_group
week_25_subset$condition <- week_25_subset$approx_group
week_25_subset$celltype <- week_25_subset$integrated_celltypel1

meta <- week_25_subset@meta.data

saveRDS(week_25_subset, "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE211191_25week_seurat.rds")

library(Seurat)
library(sceasy)
library(reticulate)
use_condaenv("base", required = TRUE)
week_25_subset[["RNA"]] <- as(week_25_subset[["RNA"]], "Assay")
sceasy::convertFormat(
  week_25_subset,
  from = "seurat",
  to = "anndata",
  outFile = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE211191_25week.h5ad"
)

sce_obj <- as.SingleCellExperiment(week_25_subset)
library(SingleCellExperiment)
if ("pca" %in% names(week_25_subset@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(week_25_subset, "pca")
}
if ("umap" %in% names(week_25_subset@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(week_25_subset, "umap")
}
saveRDS(sce_obj, file = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE211191_25week_sce.rds")







########################################################################
##################--------         GSE224626     ---------##############


library(Seurat)
setwd("/Users/xiaoying/Downloads/GSE224626_RAW")

files <- list.files(pattern = "_dge.txt.gz")

seurat_objects <- list()

for (file in files) {
  data <- read.table(gzfile(file), header = TRUE, row.names = 1, sep = "\t")
  
  sample_name <- gsub("_dge.txt.gz", "", file)
  
  seurat_obj <- CreateSeuratObject(counts = data)
  
  seurat_obj$sample <- sample_name
  
  seurat_objects[[sample_name]] <- seurat_obj
}

merged_seurat <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)], add.cell.ids = names(seurat_objects))

merged_seurat

merged_seurat@meta.data$orig.ident <- gsub("^GSM[0-9]+_", "", merged_seurat@meta.data$orig.ident)

merged_seurat@meta.data$sample <- gsub("^GSM[0-9]+_", "", merged_seurat@meta.data$sample)

head(merged_seurat@meta.data)

head(merged_seurat@meta.data)

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^mt-")

merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 400 & nFeature_RNA < 10000 & percent.mt < 15)

merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)

merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
merged_seurat <- FindNeighbors(merged_seurat)
merged_seurat <- FindClusters(merged_seurat)
merged_seurat <- RunTSNE(merged_seurat)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:20)

DimPlot(merged_seurat, reduction = "umap", group.by = "seurat_clusters")

cell_markers <- list(
  Neutrophil = c("S100a9", "S100a8", "Niacr1", "Hdc", "Irg1", "Cxcl2"),
  Monocyte = c("Ifit2", "Mx1", "Ly6c2", "Cmpk2", "Chi3l3"),
  Macrophage = c("C1qa", "C1qb", "C1qc", "Cx3cr1", "Apoe", "Mrc1"),
  DC = c("Ccr7", "Tbc1d4", "Ccl22", "Il4i1", "Fscn1", "Ccl5"),
  T = c("Cd3g", "Trbc2", "Cd8a", "Cd3d", "Cd2", "Ikzf2"),
  NK = c("Gzmc", "Gzmd", "Gzma", "Gzmg", "Gzmb", "Gzme"),
  Fibroblast = c("Col3a1", "Col1a2", "Col1a1", "Dcn", "Apod", "Sparc"),
  Osteoclast = c("Ctsk", "Mmp9", "Atp6v0d2", "Acp5", "Slc37a2", "Ckb")
)

marker_genes <- unlist(cell_markers)

head(merged_seurat@meta.data)
library(ggplot2)
dot_plot <- DotPlot(merged_seurat, 
                    features = marker_genes,
                    group.by = "seurat_clusters") + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dot_plot

cluster_to_celltype <- list(
  Neutrophil = c(6),
  Monocyte = c(9, 13),
  Macrophage = c(0,1, 2, 7, 8, 14,17),
  DC = c(11),
  T = c(10,15),
  NK = c(4,12),
  Fibroblast = c(3,5,16,19),
  Osteoclast = c(18)
)

merged_seurat$celltype <- NA

for (celltype in names(cluster_to_celltype)) {
  clusters <- cluster_to_celltype[[celltype]]
  merged_seurat$celltype[merged_seurat$seurat_clusters %in% clusters] <- celltype
}

head(merged_seurat@meta.data)

DimPlot(merged_seurat, group.by = "celltype")
DimPlot(merged_seurat, group.by = "sample")

library(dplyr)

merged_seurat$condition <- case_when(
  grepl("Drug", merged_seurat$sample) ~ "Drug",
  grepl("Polymer", merged_seurat$sample) ~ "Polymer",
  grepl("Saline", merged_seurat$sample) ~ "Saline",
  TRUE ~ "Unknown"
)

table(merged_seurat$condition)

saveRDS(merged_seurat, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE224626_seurat.rds")

merged_seurat[["RNA"]] <- as(merged_seurat[["RNA"]], "Assay")

library(sceasy)
library(reticulate)
use_condaenv("base", required = TRUE)

sceasy::convertFormat(
  merged_seurat,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE224626.h5ad"
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
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE224626_sce.rds")










########################################################################
##################--------         GSE234241    ---------###############

library(Seurat)
library(Matrix)
library(readr)

data_dir <- "//Users/xiaoying/Downloads/"

matrix_data <- readMM(file.path(data_dir, "GSE234241_matrix.mtx.gz"))
features <- read_tsv(file.path(data_dir, "GSE234241_features.tsv.gz"), col_names = FALSE)
barcodes <- read_tsv(file.path(data_dir, "GSE234241_barcodes.tsv.gz"), col_names = FALSE)

rownames(matrix_data) <- features$X1
colnames(matrix_data) <- barcodes$X1

metadata <- read.delim("/Users/xiaoying/Downloads/GSE234241_Metadata.txt.gz", header = TRUE, stringsAsFactors = FALSE)
head(metadata)
colnames(metadata)

seurat_obj <- CreateSeuratObject(counts = matrix_data, project = "GSE234241")

head(seurat_obj@meta.data)
head(metadata)

metadata <- metadata[match(colnames(seurat_obj), metadata$CellID), ]
head(metadata)
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata[, -1])
head(seurat_obj@meta.data)

table(seurat_obj$Treatment)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
head(VariableFeatures(seurat_obj))
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
DimPlot(seurat_obj, reduction = "pca")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

macrophages <- c("CTSL", "C1QB", "MAFB", "C1QA", "EMP1", "RNASE1", "MARCO", "FN1")
DCs <- c("HLA-DQB1", "HLA-DRB1", "HLA-DPB1", "HLA-DPA1", "CST3", "HLA-DRB5", "CD1C", "HLA-DRA", "HLA-DQA1")
B <- c("CD22", "MS4A1", "CD19", "CD79B", "CD79A")
Proliferating <- c("MKI67", "TOP2A")
Monocytes <- c("CLEC12A", "FCN1", "VCAN")
Neutrophils <- c("AQP9", "CXCR2", "FCGR3B")
GD_T <- c("TRGV9", "TRDV2", "TRGC1", "TRDC")
CD8_T <- c("CD8B", "CD8A", "CD3D")
NK_T <- c("PRF1", "IL2RB")
Tregs <- c("CCR4", "CTLA4", "FOXP3", "IL2RA")
CD4_T <- c("CCR7", "RCAN3", "TCF7", "IL7R")
Granulocytes <- c("FCER1A", "GATA2", "CLC")
Endothelial <- c("FLT1", "SPARC", "ENG", "LYVE1")
Hepatocytes <- c("APOA1", "APOC1", "APOC3", "CYP2E1", "ALB")

all_genes <- c(macrophages, DCs, B, Proliferating, Monocytes, Neutrophils, GD_T, CD8_T, NK_T, Tregs, CD4_T, Granulocytes, Endothelial, Hepatocytes)

dotplot <- DotPlot(seurat_obj, 
                   features = all_genes, 
                   group.by = "seurat_clusters") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  theme(axis.text.y = element_text(size = 6)) +  
  scale_size_continuous(range = c(2, 6))

dotplot

cluster_to_celltype <- list(
  "macrophages" = c(7, 20),
  "DCs" = c(8, 16, 17, 28),
  "B" = c(9, 15, 29, 34),
  "Proliferating" = c(24),
  "Monocytes" = c(3, 12, 19, 31, 35),
  "Neutrophils" = c(5, 11, 13, 25, 33),
  "GD_T" = c(6, 21),
  "CD8_T" = c(0, 10, 18),
  "NK_T" = c(4, 14),
  "Tregs" = c(2),
  "CD4_T" = c(1, 23),
  "Granulocytes" = c(27, 30),
  "Endothelial" = c(26),
  "Hepatocytes" = c(22, 32)
)

cell_types <- rep(NA, length = length(seurat_obj$seurat_clusters))

for (cell_type in names(cluster_to_celltype)) {
  clusters <- cluster_to_celltype[[cell_type]]
  cell_types[seurat_obj$seurat_clusters %in% clusters] <- cell_type
}

seurat_obj$celltype <- cell_types

head(seurat_obj@meta.data)

DimPlot(seurat_obj, group.by = "celltype", label = TRUE)

blood_subset <- subset(seurat_obj, subset = Tissue == "Blood")
fna_subset <- subset(seurat_obj, subset = Tissue == "FNA")

head(blood_subset@meta.data)
head(fna_subset@meta.data)

blood_subset$condition <- paste(blood_subset$Disease_State, blood_subset$Treatment, sep = "_")
blood_subset$condition <- gsub("Healthy_Treated", "Healthy", blood_subset$condition)
blood_subset$condition <- gsub("CHB_Treated", "CHB_treated", blood_subset$condition)
blood_subset$condition <- gsub("CHB_Untreated", "CHB_untreated", blood_subset$condition)

saveRDS(blood_subset, "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE234241_blood_seurat.rds")

library(Seurat)
library(sceasy)
library(reticulate)
use_condaenv("base", required = TRUE)
blood_subset[["RNA"]] <- as(blood_subset[["RNA"]], "Assay")
sceasy::convertFormat(
       blood_subset,
       from = "seurat",
       to = "anndata",
       outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE234241_blood.h5ad"
)

sce_obj <- as.SingleCellExperiment(blood_subset)
library(SingleCellExperiment)
if ("pca" %in% names(blood_subset@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(blood_subset, "pca")
}
if ("umap" %in% names(blood_subset@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(blood_subset, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/GSE234241_blood_sce.rds")











########################################################################
##################--------         SCP1289       ---------##############

setwd("/media/cqnu/XiaoYingxue/SCP1289")

RawCounts <- read.table("expression/20210220_NasalSwab_RawCounts.txt", header = TRUE, sep = "\t", row.names = 1)
meta_data <- read.table("metadata/20210701_NasalSwab_MetaData.txt", header = TRUE, sep = "\t", row.names = 1)

head(RawCounts)
head(meta_data)

all(rownames(RawCounts) == rownames(meta_data))
sum(is.na(RawCounts))
sum(is.na(meta_data))

library(Seurat)

seurat_obj <- CreateSeuratObject(counts = RawCounts, meta.data = meta_data)

seurat_subset <- subset(seurat_obj, subset = Cohort_Disease_WHO_Score %in% c("Control_WHO_0", "COVID19_WHO_6-8"))

table(seurat_subset@meta.data[["Cohort_Disease_WHO_Score"]])

columns_to_keep <- c(
  "nCount_RNA", "nFeature_RNA", "donor_id", "Peak_Respiratory_Support_WHO_Score", "orig.ident",
  "Percent_Mitochondrial", "SARSCoV2_PCR_Status_and_WHO_Score","Cohort_Disease_WHO_Score",
  "biosample_id", "sex", "disease","disease__ontology_label","age","Coarse_Cell_Annotations","Detailed_Cell_Annotations"
)

seurat_subset@meta.data <- seurat_subset@meta.data[, columns_to_keep]

seurat_subset@meta.data[["Peak_Respiratory_Support_WHO_Score"]] <- NULL
seurat_subset@meta.data[["SARSCoV2_PCR_Status_and_WHO_Score"]] <- NULL
seurat_subset$condition <- seurat_subset@meta.data[["Cohort_Disease_WHO_Score"]]
seurat_subset$sample <- seurat_subset@meta.data[["donor_id"]]
seurat_subset@meta.data[["donor_id"]] <- NULL

seurat_subset <- NormalizeData(seurat_subset)
row_sums <- rowSums(seurat_subset@assays$RNA$counts > 0)
table(row_sums == 0)
dim(seurat_subset)

seurat_subset <- FindVariableFeatures(seurat_subset, selection.method = "vst", nfeatures = 2000)

head(VariableFeatures(seurat_subset))

seurat_subset <- ScaleData(seurat_subset, features = rownames(seurat_subset))

seurat_subset <- RunPCA(seurat_subset, features = VariableFeatures(seurat_subset))

ElbowPlot(seurat_subset)

seurat_subset <- FindNeighbors(seurat_subset, dims = 1:30)
seurat_subset <- FindClusters(seurat_subset, resolution = 0.4)

table(seurat_subset$seurat_clusters)
seurat_subset <- RunUMAP(seurat_subset, dims = 1:30)

DimPlot(seurat_subset, reduction = "umap", group.by = "seurat_clusters")
DimPlot(seurat_subset, reduction = "umap", group.by = "Detailed_Cell_Annotations")

seurat_subset$celltype<- seurat_subset$Detailed_Cell_Annotations

saveRDS(seurat_subset, file = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/SCP1289_068_seurat.rds")

seurat_subset[["RNA"]] <- as(seurat_subset[["RNA"]], "Assay")

library(sceasy)
library(reticulate)

sceasy::convertFormat(
  seurat_subset,
  from = "seurat",
  to = "anndata",
  outFile = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/SCP1289_068.h5ad"
)

library(SingleCellExperiment)

seurat_subset[["RNA"]] <- as(seurat_subset[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seurat_subset)

if ("pca" %in% names(seurat_subset@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(seurat_subset, "pca")
}
if ("umap" %in% names(seurat_subset@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seurat_subset, "umap")
}
saveRDS(sce_obj, file = "/media/cqnu/XiaoYingxue/DA_datasets/preprocessed_data/seq_well/SCP1289_068_sce.rds")







