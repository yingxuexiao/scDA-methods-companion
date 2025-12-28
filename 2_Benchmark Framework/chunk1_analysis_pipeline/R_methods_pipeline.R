# ==== pipeline_final.R ====
options(future.globals.maxSize = 200 * 1024^3)
install.packages("pryr")
install.packages("peakRAM")
library(dplyr)
library(pryr)
library(peakRAM)

# ==== 1. Basic Path Configuration ====
base_dir <- "/app"
base_dataset_dir <- file.path(base_dir, "preprocessed_data")
da_results_dir <- file.path(base_dir, "da_results")
methods_dir <- file.path(base_dir, "methods")
log_dir <- file.path(base_dir, "logs")
tmp_dir <- file.path(base_dir, "tmp")

condition_col <- "condition"
sample_col <- "sample"
celltype_col <- "celltype"

# ==== Function to Automatically Retrieve Dataset Files ====
get_dataset_files <- function(base_dataset_dir) {
  # Recursively search for dataset files in all subfolders
  dataset_files <- list()
  
  # Supported file extension patterns
  patterns <- c("\\_seurat\\.rds$", "\\_sce\\.rds$", "\\.h5ad$")
  
  # Traverse all subfolders
  all_folders <- list.dirs(base_dataset_dir, recursive = TRUE, full.names = TRUE)
  
  for (folder in all_folders) {
    # Look for all supported dataset files in each folder
    for (pattern in patterns) {
      files <- list.files(folder, pattern = pattern, full.names = TRUE)
      
      if (length(files) > 0) {
        # Extract dataset name (remove extension)
        dataset_names <- gsub(pattern, "", basename(files))
        
        # Add file path to corresponding dataset
        for (i in seq_along(dataset_names)) {
          dataset_name <- dataset_names[i]
          file_path <- files[i]
          
          # Classify by file type
          if (grepl("\\_seurat\\.rds$", file_path)) {
            if (is.null(dataset_files[[dataset_name]])) {
              dataset_files[[dataset_name]] <- list()
            }
            dataset_files[[dataset_name]]$seurat <- file_path
          } else if (grepl("\\_sce\\.rds$", file_path)) {
            if (is.null(dataset_files[[dataset_name]])) {
              dataset_files[[dataset_name]] <- list()
            }
            dataset_files[[dataset_name]]$sce <- file_path
          } else if (grepl("\\.h5ad$", file_path)) {
            if (is.null(dataset_files[[dataset_name]])) {
              dataset_files[[dataset_name]] <- list()
            }
            dataset_files[[dataset_name]]$h5ad <- file_path
          }
        }
      }
    }
  }
  
  # Print the found dataset information
  if (length(dataset_files) > 0) {
    cat(sprintf("Found %d datasets:\n", length(dataset_files)))
    for (dataset_name in names(dataset_files)) {
      cat(sprintf("  %s:\n", dataset_name))
      if (!is.null(dataset_files[[dataset_name]]$seurat)) {
        cat(sprintf("    Seurat: %s\n", dataset_files[[dataset_name]]$seurat))
      }
      if (!is.null(dataset_files[[dataset_name]]$sce)) {
        cat(sprintf("    SCE: %s\n", dataset_files[[dataset_name]]$sce))
      }
      if (!is.null(dataset_files[[dataset_name]]$h5ad)) {
        cat(sprintf("    H5AD: %s\n", dataset_files[[dataset_name]]$h5ad))
      }
    }
  } else {
    cat("No dataset files found\n")
  }
  
  return(dataset_files)
}

# Get all dataset file paths
dataset_files <- get_dataset_files(base_dataset_dir)

# Set seurat_paths, sce_paths, etc., based on the retrieved file paths
seurat_paths <- sapply(names(dataset_files), function(x) dataset_files[[x]]$seurat, simplify = FALSE)
sce_paths <- sapply(names(dataset_files), function(x) dataset_files[[x]]$sce, simplify = FALSE)

method_list <- c(
  "DAwnn",
  "milo",
  "DCATS",
  "cna",
  "cydar",
  "Louvain",
  "DirichletReg",
  "propeller",
  "sccomp",
  "DA_seq"
)

method_data_format <- c(
  milo = "SingleCellExperiment",
  DCATS = "Seurat",
  cna = "Seurat",
  cydar = "SingleCellExperiment",
  DAwnn = "Seurat",
  Louvain = "SingleCellExperiment",
  DirichletReg = "Seurat",
  propeller = "Seurat",
  sccomp = "Seurat",
  DA_seq = "SingleCellExperiment"
)

# Whether to keep tmp cache for the method
method_keep_tmp <- c(
  DAwnn=TRUE,
  sccomp=TRUE,
  milo=FALSE,
  DCATS=FALSE,
  CTDS=FALSE,
  cna=FALSE,
  cydar=FALSE,
  Louvain=FALSE,
  DirichletReg=FALSE,
  propeller=FALSE,
  DA_seq=FALSE
)

# ==== 2. Configuration Loading ====
load_config <- function() {
  dataset_names <- names(dataset_files)  # Get dataset names from get_dataset_files
  dataset_da_dirs <- setNames(file.path(da_results_dir, dataset_names), dataset_names)
  dataset_log_dirs <- setNames(file.path(log_dir, dataset_names), dataset_names)

 
  lapply(dataset_da_dirs, dir.create, recursive=TRUE, showWarnings=FALSE)
  lapply(dataset_log_dirs, dir.create, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive=TRUE)

  method_scripts <- setNames(file.path(methods_dir, paste0("run_", method_list, ".R")), method_list)
  list(
    dataset_names=dataset_names,
    seurat_paths=seurat_paths,
    sce_paths=sce_paths,
    method_list=method_list,
    method_scripts=method_scripts,
    method_data_format=method_data_format,
    da_results_dir=dataset_da_dirs,
    dataset_log_dirs=dataset_log_dirs,
    condition_col=condition_col,
    sample_col=sample_col,
    celltype_col=celltype_col,
    tmp_dir=tmp_dir
  )
}

# ==== 3. Resource Monitoring ====
get_cache_size <- function(dir) {
  if (!dir.exists(dir)) return(0)
  tryCatch({
    as.numeric(system(paste("du -sm", shQuote(dir), "| cut -f1"), intern=TRUE))
  }, error=function(e) 0)
}

# ==== 4. Logging System ====
init_logger <- function(log_file) {
  if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive=TRUE)
  file.create(log_file, showWarnings=FALSE)
  list(
    info=function(msg){
      timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
      message <- sprintf("%s %s\n", timestamp, msg)
      cat(message)
      cat(message, file=log_file, append=TRUE)
    },
    error=function(msg){
      timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
      message <- sprintf("%s ERROR: %s\n", timestamp, msg)
      cat(message, file=stderr())
      cat(message, file=log_file, append=TRUE)
    }
  )
}

# ==== 5. Cache Cleaning ====
clean_dataset_cache <- function(keep_tmp=FALSE) {
  rm(list=c("obj","da_result"), envir=.GlobalEnv)
  gc(full=TRUE)
  
  if (!keep_tmp) {
    if (exists("tmp_dir", envir=.GlobalEnv) && dir.exists(tmp_dir)) {
      unlink(file.path(tmp_dir, "*"), recursive=TRUE, force=TRUE)
    }
  }
}

# ==== 6. Running Methods ====
run_method <- function(method, config, summary_csv) {
  global_logger <- init_logger(file.path(log_dir, "global_scheduler.log"))
  global_logger$info(sprintf("==== Starting method execution: %s ====", method))
  
  method_env <- new.env()
  sys.source(config$method_scripts[[method]], envir=method_env)
  run_fun <- get(paste0("run_", method), envir=method_env)

  for(dataset in config$dataset_names){
    log_file <- file.path(config$dataset_log_dirs[[dataset]], paste0("log_", dataset, "_", method, ".txt"))
    logger <- init_logger(log_file)
    status <- "Success"
    cache_before <- get_cache_size(config$tmp_dir)
    cache_after <- cache_before
    
    tryCatch({
      logger$info(sprintf("Starting analysis | Method: %s | Dataset: %s", method, dataset))
      data_format <- config$method_data_format[[method]]
      data_path <- if(data_format=="Seurat") config$seurat_paths[[dataset]] else config$sce_paths[[dataset]]
      if(!file.exists(data_path)) stop("Data file does not exist: ", data_path)
      obj <- readRDS(data_path)

      # ==== Method Execution Time ====
      method_start <- Sys.time()
      mem_before <- pryr::mem_used()
      peak_res <- peakRAM({
        da_result <- run_fun(obj=obj,
                             condition_col=config$condition_col,
                             sample_col=config$sample_col,
                             celltype_col=config$celltype_col)
      })
      mem_after <- pryr::mem_used()
      method_time <- as.numeric(difftime(Sys.time(), method_start, units="secs"))

      peak_mem <- if(!is.na(peak_res$Peak_RAM_Used_MiB)) {
        peak_res$Peak_RAM_Used_MiB
      } else {
        as.numeric(mem_after - mem_before)/1024^2
      }

      # ==== Save Results + Total Time ====
      result_file <- file.path(config$da_results_dir[[dataset]], paste0(dataset,"_",method,"_result.rds"))
      saveRDS(da_result, result_file)
      total_time <- as.numeric(difftime(Sys.time(), method_start, units="secs"))

      cache_after <- get_cache_size(config$tmp_dir)
      logger$info(sprintf("Method runtime: %.1f s | Total time: %.1f s | Peak memory: %.1f MB | Cache: %.1f -> %.1f MB | Dataset %s Method %s completed",
                          method_time, total_time, peak_mem, cache_before, cache_after, dataset, method))

    }, error=function(e){
      logger$error(e$message)
      status <<- "Fail"
      method_time <<- NA
      total_time <<- NA
      peak_mem <<- NA
      cache_after <- get_cache_size(config$tmp_dir)
    })

    # ==== CSV Append ====
    line <- data.frame(
      Dataset=dataset,
      Method=method,
      Status=status,
      Method_Time_sec=ifelse(!is.null(method_time), round(method_time), NA),
      Total_Time_sec=ifelse(!is.null(total_time), round(total_time), NA),
      Peak_Mem_MB=ifelse(!is.null(peak_mem), round(peak_mem,1), NA),
      Cache_Before_MB=cache_before,
      Cache_After_MB=cache_after,
      R_Version=R.version.string,
      System_CPU=as.integer(system("nproc", intern=TRUE)),
      System_Mem=gsub(" kB","",sub("MemTotal:","",system("grep MemTotal /proc/meminfo", intern=TRUE))),
      Timestamp=format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors=FALSE
    )
    write.table(line, file=summary_csv, sep=",", row.names=FALSE,
                col.names=!file.exists(summary_csv), append=TRUE)

    # ==== Cache/Environment Cleanup ====
    keep_tmp_flag <- method_keep_tmp[[method]]
    clean_dataset_cache(keep_tmp = keep_tmp_flag)
    Sys.sleep(2)
  }

  rm(list=c("method_env","run_fun"))
  gc(full=TRUE)
  global_logger$info(sprintf("Method %s execution completed | Total time: %.1f s", method, total_time))
  return(invisible(NULL))
}

# ==== 7. Main Controller ====
main <- function(){
  config <- load_config()
  summary_csv <- file.path(base_dir,"12da_summary.csv")

  global_logger <- init_logger(file.path(log_dir,"pipeline_master.log"))
  global_logger$info("=== Analysis pipeline started ===")
  global_logger$info(sprintf("Available methods: %s", paste(config$method_list,collapse=", ")))
  global_logger$info(sprintf("Datasets to process: %s", paste(config$dataset_names,collapse=", ")))

  for(method in config$method_list){
    method_start <- Sys.time()
    global_logger$info(sprintf("\n>> Starting method: %s", method))
    
    tryCatch({
      run_method(method, config, summary_csv)
      elapsed <- difftime(Sys.time(), method_start, units="mins")
      global_logger$info(sprintf("<< Method %s completed | Total time: %.1f minutes", method, elapsed))
    }, error=function(e){
      global_logger$error(sprintf("Method %s execution failed: %s", method, e$message))
    })
  }

  global_logger$info("\n=== All analyses successfully completed ===")
}

# ==== 8. Execution Entry Point ====
main()
