##dawnn_function


# tensorflow=2.13.0;
# keras=2.13.1;
# numpy=1.23.2;
# python=3.8.19;
# wheel=0.45.1;


# Step 1: Install Dawnn package (may need to install `remotes` package first)
#remotes::install_github("george-hall-ucl/dawnn")

# Step 2: Download Dawnn's model
# By default, model stored at ~/.dawnn/dawnn_nn_model.h5
#dawnn::download_model()

#  model_path <- path.expand("~/.dawnn/dawnn_nn_model.h5")
#  print(model_path)  

#  Step 3: Install Tensorflow Python package in Reticulate environment
#  reticulate::py_install("tensorflow") 


library(reticulate)
reticulate::use_condaenv("dawnn_env",required =TRUE)  
py_config()
library(Seurat)
library(dawnn)
library(dplyr)
library(ggplot2)
library(purrr)



run_da_dawnn_analysis <- function(seurat_obj, 
                              label_names ,
                              label_1 ,
                              label_2 ,
                              reduced_dim ,
                              n_dims ,
                              nn_model,
                              da_threshold = 0.1) {

  

  # 1.  run Dawnn
  seurat_obj <- run_dawnn(
    seurat_obj,
    label_names ,
    label_1 = label_1,
    label_2 = label_2,
    reduced_dim,
    n_dims ,
    nn_model 
  )
  
  # 2. Extraction result
  da_results <- seurat_obj@meta.data %>%
     dplyr::select(
         seurat_clusters,
         celltype,
         !!sym(label_names),  
         dawnn_scores,
         dawnn_lfc,
         dawnn_p_vals,
         dawnn_da_verdict
     ) %>%
     tibble::as_tibble()
  

  # 3. Statistical results
  stats_list <- list(
    by_cluster = da_results %>%
      group_by(seurat_clusters, condition) %>%
      summarise(
        n_cells = n(),
        n_da = sum(dawnn_da_verdict),
        prop_da = n_da / n_cells,
        mean_lfc = mean(dawnn_lfc),
        .groups = "drop"
      ),
    
    by_celltype = da_results %>%
      group_by(celltype, condition) %>%
      summarise(
        n_cells = n(),
        n_da = sum(dawnn_da_verdict),
        prop_da = n_da / n_cells,
        mean_lfc = mean(dawnn_lfc),
        .groups = "drop"
      )
  )
  
  return(list(
    seurat_obj = seurat_obj,
    stats = stats_list
  ))
}
















