#Drop_seq_datasets:

##########################################################################
####################------      GSE174400        ---------################

library(Seurat)
library(dplyr)
library(Matrix)

counts_file <- "/Users/xiaoying/Downloads/GSE174400_Seurat_object.RNAcounts_P14retina_ENOSproject.txt.gz"
metadata_file <- "/Users/xiaoying/Downloads/GSE174400_Seurat_object.metadata_P14retina_ENOSproject.txt.gz"

counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")
counts <- as(as.matrix(counts), "dgCMatrix")

metadata <- read.table(metadata_file, header = TRUE, row.names = 1, sep = "\t")

seu <- CreateSeuratObject(counts = counts, meta.data = metadata, min.features = 100)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

seu <- subset(seu, subset = nFeature_RNA >= 100 & percent.mt <= 10)

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

seu$celltype <- seu$Cell_Type
seu$Cell_Type <- NULL

seu$condition <- seu$Genotype
seu$Genotype <- NULL

seu$sample <- seu$Replicate
seu$Replicate <- NULL

seu_list <- SplitObject(seu, split.by = "sample")

seu_list <- lapply(seu_list, NormalizeData)
seu_list <- lapply(seu_list, FindVariableFeatures)

anchors <- FindIntegrationAnchors(object.list = seu_list, dims = 1:20)

seu_integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(seu_integrated) <- "integrated"

seu_integrated <- ScaleData(seu_integrated, verbose = FALSE)

seu_integrated <- RunPCA(seu_integrated, npcs = 20, verbose = FALSE)

seu_integrated <- FindNeighbors(seu_integrated, dims = 1:20)
seu_integrated <- FindClusters(seu_integrated, resolution = 0.5)

seu_integrated <- RunUMAP(seu_integrated, dims = 1:20)
DimPlot(seu_integrated, reduction = "umap", group.by = "celltype")
DimPlot(seu_integrated, reduction = "umap", group.by = "condition")
DimPlot(seu_integrated, reduction = "umap", group.by = "sample")

saveRDS(seu_integrated, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE174400_seurat.rds")

seu_integrated[["RNA"]] <- as(seu_integrated[["RNA"]], "Assay")

library(sceasy)
library(reticulate)
use_condaenv("base", required = TRUE)
sceasy::convertFormat(
  seu_integrated,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE174400.h5ad"
)

library(SingleCellExperiment)
seu_integrated[["RNA"]] <- as(seu_integrated[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seu_integrated)

if ("pca" %in% names(seu_integrated@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(seu_integrated, "pca")
}
if ("umap" %in% names(seu_integrated@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seu_integrated, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE174400_sce.rds")




##########################################################################
####################------      GSE180862        ---------################

library(Seurat)
library(Matrix)
library(readr)
library(dplyr)

data_dir <- "/Users/xiaoying/Downloads/GSE180862/"

mtx <- readMM(file.path(data_dir, "GSE180862_DropSeq.Cortex.7days.24hrs.TBI.Sham.digital_expression.mtx.gz"))
features <- read_tsv(file.path(data_dir, "GSE180862_DropSeq.Cortex.7days.24hrs.TBI.Sham.features.tsv.gz"),
                     col_names = FALSE)
barcodes <- read_tsv(file.path(data_dir, "GSE180862_DropSeq.Cortex.7days.24hrs.TBI.Sham.barcodes.tsv.gz"),
                     col_names = FALSE)[[1]]
meta <- read_tsv(file.path(data_dir, "GSE180862_DropSeq.Cortex.7days.24hrs.TBI.Sham.metaData.tsv.gz"))

rownames(mtx) <- make.unique(features$X2)
colnames(mtx) <- barcodes

so <- CreateSeuratObject(counts = mtx, project = "Cortex")

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")

rownames(meta) <- meta$barcode
so <- AddMetaData(so, meta)

so <- subset(so, subset = nFeature_RNA >= 200 & percent.mt <= 15)

so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so,npcs = 50, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30)
so <- FindClusters(so, resolution = 0.6)
so <- RunUMAP(so, dims = 1:30)

DimPlot(so, group.by = "NeuronCellType", label = TRUE)

table(so$Condition)
table(so$Animal)

so_7 <- subset(so, subset = Condition == "7days" & Animal %in% c("Sham","TBI"))

table(so_7$Animal)
table(so_7$Tissue)
table(so_7$Animal, so_7$Tissue)

so_7 <- NormalizeData(so_7, normalization.method = "LogNormalize", scale.factor = 1e4)
so_7 <- FindVariableFeatures(so_7, selection.method = "vst", nfeatures = 3000)
so_7 <- ScaleData(so_7, verbose = TRUE)
so_7 <- RunPCA(so_7, npcs = 50, verbose = FALSE)
so_7 <- FindNeighbors(so_7, dims = 1:30)
so_7 <- FindClusters(so_7, resolution = 0.6)
so_7 <- RunUMAP(so_7, dims = 1:30, seed.use = 123)

DimPlot(so_7, group.by = "NeuronCellType", label = TRUE, repel = TRUE)

so_7$celltype <- so_7$NeuronCellType
so_7$NeuronCellType <- NULL
so_7$sample <- so_7$Tissue
so_7$Tissue <- NULL
so_7$condition <- so_7@meta.data[["Animal"]]
so_7@meta.data[["Animal"]] <- NULL

saveRDS(so_7, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE180862_7days_seurat.rds")

so_7[["RNA"]] <- as(so_7[["RNA"]], "Assay")

library(sceasy)
library(reticulate)
use_condaenv("base", required = TRUE)
sceasy::convertFormat(
  so_7,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE180862_7days.h5ad"
)

library(SingleCellExperiment)
so_7[["RNA"]] <- as(so_7[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(so_7)

if ("pca" %in% names(so_7@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(so_7, "pca")
}
if ("umap" %in% names(so_7@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(so_7, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE180862_7days_sce.rds")




##########################################################################
####################------      GSE195986        ---------################

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(stringr)
library(ggplot2)

base_dir <- "/Users/xiaoying/Downloads/GSE195986_RAW/"

files <- c(
  "GSM5857297_HT1.dropseq.dge.txt.gz",
  "GSM5857298_HT2.dropseq.dge.txt.gz",
  "GSM5857299_HT3.dropseq.dge.txt.gz",
  "GSM5857300_HT4.dropseq.dge.txt.gz",
  "GSM5857301_HT5.dropseq.dge.txt.gz",
  "GSM5857302_HT6.dropseq.dge.txt.gz",
  "GSM5857303_HT7.dropseq.dge.txt.gz",
  "GSM5857304_T2D1.dropseq.dge.txt.gz",
  "GSM5857305_T2D2.dropseq.dge.txt.gz",
  "GSM5857306_T2D3.dropseq.dge.txt.gz",
  "GSM5857307_T2D4.dropseq.dge.txt.gz"
)

sample_condition <- tibble(
  sample = c("HT1","HT2","HT3","HT4","HT5","HT6","HT7","T2D1","T2D2","T2D3","T2D4"),
  condition = c(rep("Healthy", 7), rep("T2D", 4))
)

read_dropseq_dge <- function(fp) {
  dt <- fread(file.path(base_dir, fp))
  gene_col <- colnames(dt)[1]
  genes <- dt[[gene_col]]
  dt[[gene_col]] <- NULL
  m <- as.matrix(dt)
  rownames(m) <- genes
  m <- Matrix(m, sparse = TRUE)
  return(m)
}

make_seurat_from_dge <- function(mat, sample_name) {
  so <- CreateSeuratObject(counts = mat, project = sample_name, min.cells = 3, min.features = 100)
  so$sample <- sample_name
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so <- subset(so, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt <= 20)
  return(so)
}

objects <- list()
for (f in files) {
  message("Reading ", f)
  mat <- read_dropseq_dge(f)
  sample_name <- str_match(f, "(HT[1-7]|T2D[1-4])")[,2]
  so <- make_seurat_from_dge(mat, sample_name)
  objects[[sample_name]] <- so
}

seu <- Reduce(function(x, y) merge(x, y), objects)

seu$condition <- sample_condition$condition[match(seu$sample, sample_condition$sample)]

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 500, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

pcs_use <- 1:10

seu <- FindNeighbors(seu, reduction = "pca", dims = pcs_use, k.param = 20, verbose = FALSE)

seu <- FindClusters(seu, resolution = 0.6, cluster.name = "seurat_clusters", algorithm = 3, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = pcs_use, verbose = FALSE)

hormone_genes <- c("INS","GCG","SST","PPY","IAPP","GHRL")

get_hormone_thresholds <- function(seu, hormone_genes, layer = "counts") {
  expr <- LayerData(seu[["RNA"]], layer = layer)
  
  thresholds <- sapply(hormone_genes, function(gene) {
    if (!(gene %in% rownames(expr))) return(NA)
    vals <- expr[gene, ]
    if (sum(vals > 0) == 0) return(0)
    return(median(vals[vals > 0]) * 0.1)
  })
  
  return(thresholds)
}

hormone_thr <- get_hormone_thresholds(seu, hormone_genes, layer = "counts")

count_high_hormones <- function(seu, genes, thr, layer = "counts") {
  mat <- LayerData(seu[["RNA"]], layer = layer)
  
  present <- intersect(genes, rownames(mat))
  if (length(present) == 0) return(rep(0L, ncol(seu)))
  
  m <- mat[present, , drop = FALSE]
  res <- rep(0L, ncol(seu))
  
  for (g in present) {
    res <- res + as.integer(m[g, ] > thr[g])
  }
  return(res)
}

seu$high_hormone_count <- count_high_hormones(seu, hormone_genes, hormone_thr)

seu$doublet_hormone <- seu$high_hormone_count >= 2

message("Doublets (hormone-based) detected: ", sum(seu$doublet_hormone))

seu <- subset(seu, subset = doublet_hormone == FALSE)

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 500, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = pcs_use, k.param = 20, verbose = FALSE)

seu <- FindClusters(seu, resolution = 0.6, algorithm = 3, verbose = FALSE)
seu <- RunUMAP(seu, dims = pcs_use, verbose = FALSE)

library(ggplot2)

DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 5) +
  ggtitle("UMAP of scRNA-seq clusters") +
  theme_minimal()

marker_genes <- list(
  Alpha = "GCG",
  Beta  = "INS",
  Delta = "SST",
  PP    = "PPY",
  Duct  = "KRT19",
  PSC   = "COL1A2",
  Acinar = "REG1A",
  Endothelial = "FLT1"
)

markers_vec <- unlist(marker_genes)

dp <- DotPlot(seu, features = marker_genes, cols = c("lightgrey", "red"), dot.scale = 6) +
  RotatedAxis() +
  ggtitle("Marker gene expression across clusters")

FeaturePlot(seu, features = c("SST","PPY"), cols = c("lightgrey","red"))

cluster_ids <- levels(seu$seurat_clusters)

cluster_to_celltype <- c(
  "0" = "alpha",
  "1" = "beta",
  "2" = "alpha",
  "3" = "PSC",
  "4" = "Duct",
  "5" = "beta",
  "6" = "alpha",
  "7" = "beta",
  "8" = "Duct",
  "9" = "Endothelial",
  "10"="beta",
  "11"="Acinar",
  "12"="beta"
)

seu$celltype <- plyr::mapvalues(
  x = seu$seurat_clusters,
  from = names(cluster_to_celltype),
  to = cluster_to_celltype
)

DimPlot(seu, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 0.5) +
  ggtitle("Cell Type Annotation") +
  theme_minimal()

DimPlot(seu, reduction = "umap", group.by = "condition", label = TRUE, pt.size = 0.5)

saveRDS(seu, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE195986_seurat.rds")

seu[["RNA"]] <- as(seu[["RNA"]], "Assay")

library(sceasy)
library(reticulate)
use_condaenv("base", required = TRUE)
sceasy::convertFormat(
  seu,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE195986.h5ad"
)

library(SingleCellExperiment)
seu[["RNA"]] <- as(seu[["RNA"]], Class = "Assay")
sce_obj <- as.SingleCellExperiment(seu)

if ("pca" %in% names(seu@reductions)) {
  reducedDims(sce_obj)[["PCA"]] <- Embeddings(seu, "pca")
}
if ("umap" %in% names(seu@reductions)) {
  reducedDims(sce_obj)[["UMAP"]] <- Embeddings(seu, "umap")
}
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE195986_sce.rds")




##########################################################################
####################------      GSE225948        ---------################
library(Seurat)
library(data.table)
library(Matrix)
library(celda)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(stringr)
library(ggplot2)

options(stringsAsFactors = FALSE)
set.seed(1234)

DATA_DIR <- "/Users/xiaoying/Downloads/GSE225948_RAW/"

sample_map <- tribble(
  ~gsm,        ~file_prefix,                   ~age_group, ~condition, ~day, ~replicate, ~sex,
  "GSM7060815","Brain_GR180716",               "young",    "Sham",     2,     1,          "male",
  "GSM7060816","Brain_GR181128",               "young",    "Sham",     2,     2,          "male",
  "GSM7060817","Brain_GR181212",               "young",    "Sham",     2,     3,          "male",
  "GSM7060818","Brain_GR190110",               "young",    "Sham",     2,     4,          "male",
  "GSM7060819","Brain_GR180426",               "young",    "Stroke",   2,     1,          "male",
  "GSM7060820","Brain_GR180614",               "young",    "Stroke",   2,     2,          "male",
  "GSM7060821","Brain_GR180919",               "young",    "Stroke",   2,     3,          "male",
  "GSM7060822","Brain_GR181024",               "young",    "Stroke",   2,     4,          "male",
  "GSM7060823","Brain_GR180125",               "young",    "Stroke",   14,    1,          "male",
  "GSM7060824","Brain_GR180613",               "young",    "Stroke",   14,    2,          "male",
  "GSM7060825","Brain_GR180905",               "young",    "Stroke",   14,    3,          "male",
  "GSM7060826","Brain_GR181114",               "young",    "Stroke",   14,    4,          "male",
  "GSM7060827","Brain_aged_GR210708",          "aged",     "Sham",     2,     1,          "male",
  "GSM7060828","Brain_aged_GR200728",          "aged",     "Sham",     2,     1,          "female",
  "GSM7060829","Brain_aged_GR200723",          "aged",     "Stroke",   2,     1,          "male",
  "GSM7060830","Brain_aged_GR200716",          "aged",     "Stroke",   2,     1,          "female",
  "GSM7060831","Brain_aged_GR210225",          "aged",     "Stroke",   2,     2,          "female",
  "GSM7060832","Brain_aged_GR200812",          "aged",     "Stroke",   14,    1,          "female"
)

sample_map <- sample_map %>%
  mutate(
    counts_file = file.path(DATA_DIR, paste0(gsm, "_", file_prefix, "_counts.csv.gz")),
    meta_file   = file.path(DATA_DIR, paste0(gsm, "_", file_prefix, "_metadata.csv.gz")),
    sample_id   = sub(".*_", "", file_prefix)
  )

stopifnot(all(file.exists(sample_map$counts_file)),
          all(file.exists(sample_map$meta_file)))

process_one_sample <- function(counts_path, meta_path, sample_id,
                               umi_min = 200, umi_max = 10000, mt_max = 20,
                               pcs_df = 20) {
  message(">>> Loading sample: ", sample_id)
  
  raw_counts <- fread(counts_path)
  rn <- raw_counts[[1]]
  raw_counts[[1]] <- NULL
  mat <- as.matrix(raw_counts)
  rownames(mat) <- rn
  
  meta <- fread(meta_path)
  
  seu <- CreateSeuratObject(counts = mat, project = sample_id)
  seu$sample <- sample_id
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
  seu <- subset(seu, subset = nCount_RNA >= umi_min & nCount_RNA <= umi_max & percent.mt <= mt_max)
  
  message(">>> DecontX: ", sample_id)
  decont <- decontX(as.matrix(GetAssayData(seu, layer = "counts")))
  seu[["RNA"]] <- CreateAssayObject(counts = decont$decontXcounts)
  DefaultAssay(seu) <- "RNA"
  
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = max(30, pcs_df), verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:pcs_df, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
  
  message(">>> Done sample: ", sample_id, " | Cells kept: ", ncol(seu))
  return(seu)
}

brain_objs <- list()
for (i in seq_len(nrow(sample_map))) {
  r <- sample_map[i, ]
  obj <- process_one_sample(
    counts_path = r$counts_file,
    meta_path   = r$meta_file,
    sample_id   = r$sample_id
  )
  
  obj$age_group <- r$age_group
  obj$condition <- r$condition
  obj$day       <- r$day
  obj$replicate <- r$replicate
  obj$sex       <- r$sex
  
  brain_objs[[r$sample_id]] <- obj
}

seu@meta.data <- data.frame(lapply(seu@meta.data, function(x) {
  if (is.list(x)) return(as.character(unlist(x)))
  else return(x)
}), row.names = rownames(seu@meta.data))

library(DoubletFinder)
library(dplyr)

pcs_df <- 20
doublet_rate <- 0.05

for (sample_id in names(brain_objs)) {
  seu <- brain_objs[[sample_id]]
  message(">>> Processing DoubletFinder: ", sample_id)
  
  sweep.res.list <- paramSweep(seu, PCs = 1:pcs_df, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  bcmvn <- as.data.frame(bcmvn)
  bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
  bcmvn$BCmetric <- as.numeric(as.character(bcmvn$BCmetric))
  
  best_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  
  annotations <- seu$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp <- round(doublet_rate * ncol(seu))
  nExp.adj <- round(nExp * (1 - homotypic.prop))
  
  seu <- doubletFinder(
    seu,
    PCs = 1:pcs_df,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp.adj,
    sct = FALSE
  )
  
  dcols <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)
  stopifnot(length(dcols) >= 1)
  class_col <- tail(dcols, 1)
  
  seu <- subset(seu, subset = !!as.name(class_col) == "Singlet")
  
  message(">>> DoubletFinder done: ", sample_id, " | Cells kept: ", ncol(seu))
  brain_objs[[sample_id]] <- seu
}

library(Seurat)

combined <- merge(
  x = brain_objs[[1]],
  y = brain_objs[-1],
  add.cell.ids = names(brain_objs),
  project = "Brain_Combined"
)

combined

all_pann_df <- grep("^(pANN|DF.classifications)_", colnames(combined@meta.data), value = TRUE)

keep_cols <- c("pANN_0.25_0.01_132", "DF.classifications_0.25_0.01_132")

cols_to_remove <- setdiff(all_pann_df, keep_cols)

combined@meta.data <- combined@meta.data[, !colnames(combined@meta.data) %in% cols_to_remove]

colnames(combined@meta.data)

combined$condition_day <- paste0(combined$condition, "_", combined$day, "d")
table(combined$condition_day)

combined@meta.data[["replicate"]] <- NULL

combined$age <- combined$age_group
combined$age_group <- NULL

combined$treatment <- combined$condition
combined$condition<- NULL

combined$condition <- combined$condition_day
combined$condition_day <- NULL

combined$condition[combined$condition == "Sham_2d"]    <- "Sham"
combined$condition[combined$condition == "Stroke_2d"]  <- "D02"
combined$condition[combined$condition == "Stroke_14d"] <- "D14"

table(combined$condition)

gsm2sample <- c(
  "GR180716" = "Sham_D02_1",
  "GR181128" = "Sham_D02_2",
  "GR181212" = "Sham_D02_3",
  "GR190110" = "Sham_D02_4",
  "GR180426" = "Stroke_D02_1",
  "GR180614" = "Stroke_D02_2",
  "GR180919" = "Stroke_D02_3",
  "GR181024" = "Stroke_D02_4",
  "GR180125" = "Stroke_D14_1",
  "GR180613" = "Stroke_D14_2",
  "GR180905" = "Stroke_D14_3",
  "GR181114" = "Stroke_D14_4",
  "GR210708" = "Sham_D02_aged_1",
  "GR200728" = "Sham_D02_aged_2",
  "GR200723" = "Stroke_D02_aged_1",
  "GR200716" = "Stroke_D02_aged_2",
  "GR210225" = "Stroke_D02_aged_3",
  "GR200812" = "Stroke_D14_aged_1"
)

combined$sample <- unname(gsm2sample[combined$sample])
table(is.na(combined$sample))

table(combined$condition, combined$sample)

library(Seurat)
library(harmony)

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)

combined <- ScaleData(combined, features = VariableFeatures(combined))
combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = 50, verbose = FALSE)

ElbowPlot(combined, ndims = 50)

combined <- RunHarmony(
  object = combined,
  group.by.vars = "sample"
)

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:40)

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:40)
combined <- FindClusters(combined, resolution = 0.7)

DimPlot(combined, reduction = "umap", group.by = "condition")
DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

library(SingleR)
library(celldex)
library(Seurat)
library(dplyr)

immgen_ref <- ImmGenData()

expr_mat <- as.matrix(GetAssayData(combined, slot = "data"))

mouse_rna_seq_ref <- MouseRNAseqData()

singleR_immgen <- SingleR(
  test = expr_mat,
  ref  = immgen_ref,
  labels = immgen_ref$label.main
)

mouse_rna_seq_singleR <- SingleR(
  test = expr_mat,
  ref = mouse_rna_seq_ref,
  labels = colData(mouse_rna_seq_ref)$label.main
)

mouse_rna_seq_singleR <- SingleR(test = combined,
                                 ref = mouse_rna_seq_ref, 
                                 labels = mouse_rna_seq_ref$label.main)

combined$SingleR_ImmGen <- singleR_immgen$labels

DimPlot(combined, reduction = "umap", group.by = "SingleR_ImmGen", label = TRUE)

combined$SingleR_MouseRNAseq <- mouse_rna_seq_singleR$labels
DimPlot(combined, reduction = "umap", group.by = "SingleR_MouseRNAseq", label = TRUE)

p1 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Seurat Clusters") + theme_minimal()

p2 <- DimPlot(combined, reduction = "umap", group.by = "SingleR_MouseRNAseq", label = TRUE) +
  ggtitle("SingleR Annotation") + theme_minimal()

p1 + p2

library(Seurat)
library(ggplot2)
library(dplyr)

markers_list <- list(
  Microglia        = c("Hexb", "P2ry12", "Cx3cr1", "Siglech", "Tmem119"),
  BAMs             = c("Mrc1", "Pf4", "Cbr2", "Cd163", "Lyve1", "H2-Ab1"),
  MdCs             = c("Itgax", "Csf1r", "Cd68"),
  Granulocytes     = c("S100a8", "S100a9", "Cxcr2", "Ly6g"),
  Mast_cells       = c("Cpa3"),
  DCs              = c("Cd209a", "Xcr1", "Ccr7"),
  T_cells          = c("Cd3d", "Cd3e", "Trac", "Trbc2"),
  NK_cells         = c("Nkg7", "Gzma", "Ncr1"),
  B_cells          = c("Ms4a1", "Cd19", "Ly6d"),
  ECs              = c("Cldn5", "Ly6c1", "Flt1"),
  Vascular_mural   = c("Acta2", "Tagln", "Vtn"),
  Epithelial_like  = c("Krt8", "Krt18"),
  Oligodendrocytes = c("Mbp", "Plp1")
)

all_markers <- unlist(markers_list)
marker_classes <- rep(names(markers_list), sapply(markers_list, length))

p <- DotPlot(
  combined, 
  features = all_markers, 
  group.by = "seurat_clusters"
) + 
  RotatedAxis() +
  theme_minimal() +
  labs(title = "Marker Gene Expression per Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

celltype_annotations <- c(
  "0"  = "ECs",
  "1"  = "Microglia",
  "2"  = "Microglia",
  "3"  = "Microglia",
  "4"  = "MdCs",
  "5"  = "MdCs",
  "6"  = "BAMs",
  "7"  = "Microglia",
  "8"  = "T_cells",
  "9"  = "MdCs",
  "10" = "Microglia",
  "11" = "Granulocytes",
  "12" = "Microglia",
  "13" = "BAMs",
  "14" = "BAMs",
  "15" = "NK_cells",
  "16" = "B_cells",
  "17" = "DCs",
  "18" = "Rare_cell",
  "19" = "Oligodendrocytes",
  "20" = "ECs"
)

combined$celltype <- plyr::mapvalues(
  combined$seurat_clusters,
  from = names(celltype_annotations),
  to   = celltype_annotations
)

DimPlot(combined, group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("UMAP with Cell Type Annotations")

saveRDS(combined, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE225948_seurat.rds")

combined[["RNA"]] <- as(combined[["RNA"]], "Assay")

library(sceasy)
library(reticulate)
use_condaenv("base", required = TRUE)
sceasy::convertFormat(
  combined,
  from = "seurat",
  to = "anndata",
  outFile = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE225948.h5ad"
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
saveRDS(sce_obj, file = "/Volumes/XiaoYingxue/DA_datasets/preprocessed_data/drop_seq/GSE225948_sce.rds")






