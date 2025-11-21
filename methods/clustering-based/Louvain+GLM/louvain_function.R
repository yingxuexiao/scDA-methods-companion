library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(scran)
library(stringr)
library(igraph)
library(edgeR)
library(purrr)

run_da_louvain_nbglm <- function(
    sce_obj,
    condition_col,
    sample_col,
    reduced.dim,
    d=30,
    k=15,
    resolution=1,
    batch=NULL,
    norm.method="TMM"
){
  ## Make design matrix
  colnames_needed <- c("sample", "condition")
   if (!is.null(batch)) colnames_needed <- c(colnames_needed, "batch")

  
  # Check if these columns exist in colData
  missing_cols <- setdiff(colnames_needed, colnames(colData(sce_obj)))

  if (length(missing_cols) > 0) {
     stop("Missing column(s) in colData(sce_obj): ", paste(missing_cols, collapse = ", "))
  }
   
   # Extract the existing columns and convert to design_df
  design_df <- as_tibble(colData(sce_obj)[, colnames_needed]) %>%
  distinct()

  condition <- "condition"
   if (is.null(batch)) {
    design <- formula(paste('~', condition, collapse = ' '))  
   } else {
     design <- formula(paste('~', batch, "+", condition, collapse = ' '))
   }


  ## Louvain clustering
  X_red_dim = reducedDim(sce_obj, reduced.dim)[,1:d]

  sce.graph <- buildKNNGraph(t(X_red_dim), k=k)

  louvain.clust <- cluster_louvain(sce.graph, resolution=resolution)

  message(str_c("#cluster: ", length(louvain.clust)))

  louvain.clust.ids <- membership(louvain.clust)
  
  condition_vec <- colData(sce_obj)[["condition"]]

  sample_vec <- colData(sce_obj)[["sample"]]

  clust.df <- data.frame("cell_id"=colnames(sce_obj), "Louvain.Clust"=as.character(louvain.clust.ids))
 
  clust.df$sample <- sample_vec
   
  clust.df$Condition <- condition_vec


  louvain.count <- table(clust.df$Louvain.Clust, clust.df$sample)

  attributes(louvain.count)$class <- "matrix"
  
  ## Test with the same NB-GLM model as the Milo
  if(norm.method %in% c("TMM")){
    message("Using TMM normalisation")
    dge <- DGEList(counts=louvain.count,
                   lib.size=colSums(louvain.count))
    dge <- calcNormFactors(dge, method="TMM")
  } else if(norm.method %in% c("logMS")){
    message("Using logMS normalisation")
    dge <- DGEList(counts=louvain.count,
                   lib.size=colSums(louvain.count))
  }
  
  model <- model.matrix(design, data=design_df)
  rownames(model) <- design_df$sample
  model <- model[colnames(louvain.count), ]
  #model <- model[colnames(louvain.count), , drop = FALSE] 
  
  dge <- estimateDisp(dge, model)
  fit <- glmQLFit(dge, model, robust=TRUE)
  n.coef <- ncol(model)
  louvain.res <- as.data.frame(topTags(glmQLFTest(fit, coef=n.coef), sort.by='none', n=Inf))
  
  clust.df$logFC <- louvain.res[clust.df$Louvain.Clust, 'logFC']

  clust.df$FDR <- louvain.res[clust.df$Louvain.Clust, 'FDR']

  ###### ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  ######           New features as you requested              ######

  # Celltype annotation table
  cell_annotations <- data.frame(
    cell_id = colnames(sce_obj),
    celltype = sce_obj$celltype,
    sample = sce_obj$sample,
    condition = sce_obj$condition
  )

  # Merge annotation
  res_annotated <- clust.df %>%
    left_join(cell_annotations, by = "cell_id")

  # Get all condition combinations
  unique_conditions <- unique(res_annotated$condition)
  condition_pairs <- combn(unique_conditions, 2, simplify = FALSE)

  # Extract aggregated information for each pair of conditions
  res_all_pairs <- map_dfr(condition_pairs, function(cond_pair) {
    cond_1 <- cond_pair[1]
    cond_2 <- cond_pair[2]
    
    res_pair <- res_annotated %>%
      filter(condition %in% c(cond_1, cond_2)) %>%
      distinct(celltype, Louvain.Clust, logFC, FDR) %>%
      mutate(
        condition_1 = cond_1,
        condition_2 = cond_2,
        significant = ifelse(FDR < 0.05, "yes", "no")
      ) %>%
      select(celltype, Louvain.Clust, condition_1, condition_2, logFC, FDR, significant)

    return(res_pair)
  })

  # Return raw and parsed results
  return(list(
    raw_result = clust.df,
    parsed_result = res_all_pairs
  ))
}
