# ============================================================================
# DA Method Batch Analysis Script - Supports Excel Format
# Function: Batch process 21 datasets (supports xlsx format), calculate performance metrics, and generate visualizations
# ============================================================================

# Load necessary packages
library(pROC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(readr)
library(readxl)  # Add readxl package to read Excel files

# ============================================================================
# 1. Define File Paths
# ============================================================================

# Input folder (containing 21 datasets in Excel format)
input_folder <- "/Users/xiaoying/Desktop/DA_Results_Analysis/21_Datasets_Accuracy_Interpretation/21_Datasets_Results/"

# Output folder
output_folder <- "/Users/xiaoying/Desktop/DA_Results_Analysis/21_Datasets_Accuracy_Interpretation/"
metrics_folder <- file.path(output_folder, "Metrics_Results")
plots_folder <- file.path(output_folder, "Visualizations")
pdf_folder <- file.path(output_folder, "PDF_Graphs")

# Create output folders
dir.create(metrics_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(pdf_folder, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 2. Define Core Calculation Function (Unchanged)
# ============================================================================

calculate_single_dataset_enhanced <- function(df) {
  # Extract true labels (already in 0/1 format)
  truth_difference <- as.numeric(df[1, ])
  
  # Ensure no NA values
  if(any(is.na(truth_difference))) {
    warning("NA values detected in true labels, they will be removed")
    truth_difference <- truth_difference[!is.na(truth_difference)]
  }
  
  results <- list()
  methods <- rownames(df)[-1]
  
  for (method in methods) {
    # Get consistency metric
    consistency <- as.numeric(df[method, ])
    
    # Check data length match
    if(length(consistency) != length(truth_difference)) {
      warning(sprintf("Data length for method %s does not match true labels", method))
      next
    }
    
    # Derive predictions
    predictions <- ifelse(consistency == 1, truth_difference, 1 - truth_difference)
    
    # Compute confusion matrix
    TP <- sum(truth_difference == 1 & predictions == 1)
    TN <- sum(truth_difference == 0 & predictions == 0)
    FP <- sum(truth_difference == 0 & predictions == 1)
    FN <- sum(truth_difference == 1 & predictions == 0)
    
    # Calculate basic metrics
    total <- TP + TN + FP + FN
    
    TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
    TNR <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
    FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), NA)
    FNR <- ifelse((FN + TP) > 0, FN / (FN + TP), NA)
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
    Recall <- TPR
    F1 <- ifelse((Precision + Recall) > 0, 
                 2 * (Precision * Recall) / (Precision + Recall), 
                 NA)
    Accuracy <- ifelse(total > 0, (TP + TN) / total, NA)
    FDR <- ifelse((FP + TP) > 0, FP / (FP + TP), NA)
    NPV <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)
    
    # Compute AUC
    auc_value <- NA
    tryCatch({
      if (length(unique(truth_difference)) > 1 && 
          length(unique(consistency)) > 1 &&
          sum(!is.na(truth_difference)) >= 2 &&
          sum(!is.na(consistency)) >= 2) {
        
        valid_idx <- !is.na(truth_difference) & !is.na(consistency)
        if(sum(valid_idx) >= 2) {
          roc_obj <- roc(truth_difference[valid_idx], consistency[valid_idx])
          auc_value <- auc(roc_obj)
        }
      }
    }, error = function(e) {
      # AUC calculation failed
    })
    
    # If AUC is still NA, use estimated value
    if(is.na(auc_value) && !is.na(TPR) && !is.na(FPR)) {
      auc_value <- (1 + TPR - FPR) / 2
    }
    
    # Calculate other useful metrics
    consistency_rate <- mean(consistency, na.rm = TRUE)
    
    # Matthews Correlation Coefficient (MCC)
    numerator <- TP * TN - FP * FN
    denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    MCC <- ifelse(denominator > 0, numerator / denominator, NA)
    
    # Balanced Accuracy
    Balanced_Accuracy <- ifelse(!is.na(TPR) && !is.na(TNR), (TPR + TNR) / 2, NA)
    
    # Youden Index
    Youden_Index <- ifelse(!is.na(TPR) && !is.na(TNR), TPR + TNR - 1, NA)
    
    # Store results
    results[[method]] <- data.frame(
      Method = method,
      TP = TP, TN = TN, FP = FP, FN = FN,
      Total = total,
      TPR = round(TPR, 4),
      TNR = round(TNR, 4),
      FPR = round(FPR, 4),
      FNR = round(FNR, 4),
      Precision = round(Precision, 4),
      Recall = round(Recall, 4),
      F1 = round(F1, 4),
      Accuracy = round(Accuracy, 4),
      FDR = round(FDR, 4),
      NPV = round(NPV, 4),
      MCC = round(MCC, 4),
      Balanced_Accuracy = round(Balanced_Accuracy, 4),
      Youden_Index = round(Youden_Index, 4),
      AUC = round(auc_value, 4),
      Consistency_Rate = round(consistency_rate, 4),
      Positive_Cells = sum(truth_difference == 1),
      Negative_Cells = sum(truth_difference == 0),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  if(length(results) == 0) {
    warning("No successful results calculated for any methods")
    return(NULL)
  }
  
  result_df <- do.call(rbind, results)
  
  # Sort by F1 score
  result_df <- result_df %>% arrange(desc(F1))
  
  return(result_df)
}

# ============================================================================
# 3. Define Visualization Function (Unchanged)
# ============================================================================

visualize_single_dataset <- function(metrics_df, dataset_name = "Dataset") {
  
  # 1. Performance Bar Chart
  p1 <- metrics_df %>%
    mutate(Method = reorder(Method, F1)) %>%
    ggplot(aes(x = Method, y = F1, fill = Method)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(F1, 3)), hjust = -0.2, size = 3.5) +
    coord_flip() +
    scale_fill_brewer(palette = "Set3") +
    labs(title = paste("F1 Score by Method -", dataset_name),
         x = "Method", y = "F1 Score") +
    theme_minimal() +
    theme(
      legend.position = "none",
      # Remove background grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Add only outer border
      panel.border = element_rect(colour = "gray", fill = NA, size = 0.7)
    )
  
  # 2. TPR-FDR Scatter Plot - Add outer borders, remove background grid lines
  p2 <- ggplot(metrics_df, aes(x = FDR, y = TPR)) +
    geom_point(aes(size = AUC, color = F1), alpha = 0.8) +
    geom_text_repel(aes(label = Method), size = 3.5, max.overlaps = 20) +
    geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.3) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(3, 8)) +
    labs(title = "TPR vs FDR Trade-off",
         x = "False Discovery Rate (FDR)",
         y = "True Positive Rate (TPR)",
         color = "F1 Score", size = "AUC") +
    theme_minimal() +
    theme(
      # Remove background grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Add only outer border
      panel.border = element_rect(colour = "gray", fill = NA, size = 0.7)
    )
  
  # 3. Confusion Matrix Heatmap
  confusion_data <- metrics_df %>%
    select(Method, TP, TN, FP, FN) %>%
    pivot_longer(cols = c(TP, TN, FP, FN), 
                 names_to = "Metric", values_to = "Count") %>%
    mutate(Result_Type = ifelse(Metric %in% c("TP", "TN"), "Correct", "Error"))
  
  p3 <- ggplot(confusion_data, aes(x = Metric, y = Method, fill = Count)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Count), color = "black", size = 4) +
    scale_fill_gradient(low = "white", high = "#2E8B57") +
    facet_grid(~Result_Type, scales = "free", space = "free") +
    labs(title = "Confusion Matrix Elements",
         x = "Confusion Matrix Element",
         y = "Method") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      # Remove background grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Add only outer border
      panel.border = element_rect(colour = "gray", fill = NA, size = 0.7)
    )
  
  # Use patchwork to combine the plots
  combined_plot <- (p1 | p2) / p3 +
    plot_annotation(title = paste("DA Methods Analysis -", dataset_name),
                    subtitle = paste("Total cells:", unique(metrics_df$Total_Cells)[1],
                                     "| Positive cells:", unique(metrics_df$Positive_Cells)[1]))
  
  return(list(p1 = p1, p2 = p2, p3 = p3, combined = combined_plot))
}

# ============================================================================
# 4. Read and Process Functions (Supports Excel)
# ============================================================================

# Function to read Excel files
read_excel_file <- function(file_path) {
  # Try to read Excel file
  tryCatch({
    # Read the first worksheet
    data_raw <- read_excel(file_path, col_names = TRUE)
    
    # Convert to dataframe
    data_df <- as.data.frame(data_raw)
    
    # Check data format
    if(nrow(data_df) < 2 || ncol(data_df) < 2) {
      warning(sprintf("File %s has too few data", basename(file_path)))
      return(NULL)
    }
    
    # Set row names (first column as method names)
    if(ncol(data_df) > 0) {
      # First column as row names
      row_names <- as.character(data_df[, 1])
      # Check for duplicate row names
      if(any(duplicated(row_names))) {
        warning(sprintf("File %s has duplicate row names", basename(file_path)))
        # Add suffix to make unique
        row_names <- make.unique(row_names)
      }
      rownames(data_df) <- row_names
      
      # Remove the first column
      if(ncol(data_df) > 1) {
        data_df <- data_df[, -1, drop = FALSE]
      }
    }
    
    # Check if dataframe is valid
    if(nrow(data_df) < 1 || ncol(data_df) < 1) {
      warning(sprintf("File %s has empty data after processing", basename(file_path)))
      return(NULL)
    }
    
    return(data_df)
    
  }, error = function(e) {
    warning(sprintf("Error reading file %s: %s", basename(file_path), e$message))
    return(NULL)
  })
}

# Smart file reading function (auto detects format)
read_data_file <- function(file_path) {
  cat(sprintf("  Reading file: %s\n", basename(file_path)))
  
  # Check file extension
  file_ext <- tolower(tools::file_ext(file_path))
  
  if(file_ext %in% c("xlsx", "xls")) {
    # Excel file
    cat("  Format: Excel\n")
    return(read_excel_file(file_path))
    
  } else if(file_ext == "csv") {
    # CSV file
    cat("  Format: CSV\n")
    tryCatch({
      data_raw <- read_csv(file_path, show_col_types = FALSE)
      data_df <- as.data.frame(data_raw)
      
      # Set row names and columns (same as Excel handling)
      if(ncol(data_df) > 0) {
        rownames(data_df) <- data_df[, 1]
        data_df <- data_df[, -1, drop = FALSE]
      }
      return(data_df)
    }, error = function(e) {
      warning(sprintf("Failed to read CSV file: %s", e$message))
      return(NULL)
    })
    
  } else if(file_ext %in% c("tsv", "txt")) {
    # TSV or TXT file
    cat("  Format: TSV/TXT\n")
    tryCatch({
      data_raw <- read_delim(file_path, delim = "\t", show_col_types = FALSE)
      data_df <- as.data.frame(data_raw)
      
      # Set row names and columns (same as Excel handling)
      if(ncol(data_df) > 0) {
        rownames(data_df) <- data_df[, 1]
        data_df <- data_df[, -1, drop = FALSE]
      }
      return(data_df)
    }, error = function(e) {
      warning(sprintf("Failed to read TSV file: %s", e$message))
      return(NULL)
    })
    
  } else {
    warning(sprintf("Unsupported file format: %s", file_ext))
    return(NULL)
  }
}

# Process a single dataset
process_single_dataset <- function(file_path, dataset_id) {
  # Clean up memory
  gc()
  
  # Extract dataset name (without extension)
  dataset_name <- tools::file_path_sans_ext(basename(file_path))
  
  cat(sprintf("\n%s Processing dataset: %s (%s) %s\n",
              "=" %>% rep(20) %>% paste(collapse = ""),
              dataset_name,
              file_path,
              "=" %>% rep(20) %>% paste(collapse = "")))
  
  # Read data (using new read function)
  cat("1. Reading data...\n")
  data_df <- read_data_file(file_path)
  
  if(is.null(data_df)) {
    cat("   Reading failed, skipping this dataset\n")
    return(NULL)
  }
  
  cat(sprintf("   Data dimensions: %d rows × %d columns\n", nrow(data_df), ncol(data_df)))
  cat(sprintf("   Number of methods: %d\n", nrow(data_df) - 1))
  
  # Check data structure
  cat("2. Checking data structure...\n")
  
  # Check if first row is the true label (should have 0 and 1)
  first_row <- as.numeric(data_df[1, ])
  unique_values <- unique(na.omit(first_row))
  
  if(!all(unique_values %in% c(0, 1))) {
    cat(sprintf("   Warning: First row contains non-0/1 values: %s\n", paste(unique_values, collapse = ", ")))
    cat("   Attempting to convert T/F to 0/1...\n")
    
    # Attempt to convert T/F
    first_row_char <- as.character(data_df[1, ])
    first_row <- ifelse(first_row_char == "T", 1, 
                        ifelse(first_row_char == "F", 0, NA))
    data_df[1, ] <- first_row
  }
  
  # Compute metrics
  cat("3. Calculating performance metrics...\n")
  metrics_result <- calculate_single_dataset_enhanced(data_df)
  
  if(is.null(metrics_result)) {
    cat("   Metrics calculation failed, skipping this dataset\n")
    return(NULL)
  }
  
  cat(sprintf("   Successfully calculated metrics for %d methods\n", nrow(metrics_result)))
  
  # Save metrics results
  cat("4. Saving metrics results...\n")
  metrics_file <- file.path(metrics_folder, sprintf("%s_metrics.csv", dataset_name))
  write.csv(metrics_result, metrics_file, row.names = FALSE)
  cat(sprintf("   Metrics saved: %s\n", metrics_file))
  
  # Generate visualization
  cat("5. Generating visualization...\n")
  plots <- visualize_single_dataset(metrics_result, dataset_name)
  
  # Save as PNG
  png_file <- file.path(plots_folder, sprintf("%s_analysis.png", dataset_name))
  ggsave(png_file, plots$combined, width = 14, height = 10, dpi = 300, bg = "white")
  cat(sprintf("   PNG plot saved: %s\n", png_file))
  
  # Save as PDF (landscape, 10x8 inches)
  pdf_file <- file.path(pdf_folder, sprintf("%s_analysis.pdf", dataset_name))
  ggsave(pdf_file, plots$combined, 
         width = 10, height = 8,  # custom 10*8 inches
         device = "pdf", bg = "white")
  cat(sprintf("   PDF plot saved: %s (10x8 inches, landscape)\n", pdf_file))
  
  # Generate brief report
  cat("6. Generating brief report...\n")
  report_file <- file.path(metrics_folder, sprintf("%s_report.txt", dataset_name))
  sink(report_file)
  cat(paste("=" %>% rep(50), collapse = ""), "\n")
  cat(sprintf("DA Method Performance Analysis Report - %s\n", dataset_name))
  cat(paste("=" %>% rep(50), collapse = ""), "\n\n")
  
  cat("Method ranking (by F1 score):\n")
  top_n <- min(3, nrow(metrics_result))
  for(i in 1:top_n) {
    m <- metrics_result[i, ]
    rank_symbol <- c("🥇", "🥈", "🥉")[i]
    cat(sprintf("  %s %s: F1=%.3f, Accuracy=%.3f, AUC=%.3f\n",
                rank_symbol, m$Method, m$F1, m$Accuracy, m$AUC))
  }
  
  sink()
  cat(sprintf("   Report saved: %s\n", report_file))
  
  # Clean up all objects related to the current dataset
  rm(data_df, metrics_result)
  gc()
  
  cat(sprintf("✅ Dataset %s processed successfully!\n", dataset_name))
  
  return(list(
    dataset_name = dataset_name,
    file_path = file_path
  ))
}

# ============================================================================
# 5. Main Program: Batch Process All Datasets
# ============================================================================

batch_process_all_datasets <- function(input_folder) {
  cat("\n")
  cat(paste("=" %>% rep(60), collapse = ""), "\n")
  cat("Starting batch processing of datasets (supports Excel format)\n")
  cat(paste("=" %>% rep(60), collapse = ""), "\n\n")
  
  # Get all data files (supports multiple formats)
  data_files <- list.files(input_folder, 
                           pattern = "\\.(xlsx|xls|csv|tsv|txt)$", 
                           full.names = TRUE,
                           ignore.case = TRUE)
  
  if(length(data_files) == 0) {
    cat("Error: No data files found in the specified folder\n")
    cat(sprintf("Folder: %s\n", input_folder))
    cat("Supported file formats: .xlsx, .xls, .csv, .tsv, .txt\n")
    return(NULL)
  }
  
  cat(sprintf("Found %d data files:\n", length(data_files)))
  for(i in 1:length(data_files)) {
    cat(sprintf("  %d. %s\n", i, basename(data_files[i])))
  }
  cat("\n")
  
  # Store all results
  all_results <- list()
  success_count <- 0
  failed_datasets <- c()
  
  # Process each dataset
  for(i in 1:length(data_files)) {
    file_path <- data_files[i]
    
    cat(sprintf("\n>>> Processing dataset %d/%d <<<\n", i, length(data_files)))
    
    # Process a single dataset
    result <- process_single_dataset(file_path, i)
    
    if(!is.null(result)) {
      all_results[[result$dataset_name]] <- result
      success_count <- success_count + 1
    } else {
      failed_datasets <- c(failed_datasets, basename(file_path))
    }
    
    # Clean up memory
    gc()
    
    # Progress bar
    progress <- round(i/length(data_files) * 100)
  }
  
  # Generate comprehensive report
  cat("\n")
  cat(paste("=" %>% rep(60), collapse = ""), "\n")
  cat("Batch processing completed!\n")
  cat(paste("=" %>% rep(60), collapse = ""), "\n\n")
  
  cat(sprintf("Processing results summary:\n"))
  cat(sprintf("  Successfully processed: %d datasets\n", success_count))
  cat(sprintf("  Failed: %d datasets\n", length(failed_datasets)))
  
  if(length(failed_datasets) > 0) {
    cat("  Failed datasets:\n")
    for(ds in failed_datasets) {
      cat(sprintf("    - %s\n", ds))
    }
  }
  
  # If successful datasets were processed, attempt to generate a summary
  if(success_count > 0) {
    cat("\nAttempting to generate comprehensive summary...\n")
    
    # Try to read all successful dataset's metrics files
    all_metrics <- list()
    
    for(ds_name in names(all_results)) {
      metrics_file <- file.path(metrics_folder, sprintf("%s_metrics.csv", ds_name))
      
      if(file.exists(metrics_file)) {
        tryCatch({
          metrics_df <- read.csv(metrics_file)
          metrics_df$Dataset <- ds_name
          all_metrics[[ds_name]] <- metrics_df
          cat(sprintf("   Loaded: %s\n", ds_name))
        }, error = function(e) {
          cat(sprintf("   Failed to load: %s (%s)\n", ds_name, e$message))
        })
      }
    }
    
    if(length(all_metrics) > 0) {
      # Combine all metrics
      combined_metrics <- do.call(rbind, all_metrics)
      
      # Save combined metrics
      combined_file <- file.path(output_folder, "all_datasets_combined_metrics.csv")
      write.csv(combined_metrics, combined_file, row.names = FALSE)
      cat(sprintf("   Combined metrics saved: %s\n", combined_file))
      
      # Generate method performance summary
      method_summary <- combined_metrics %>%
        group_by(Method) %>%
        summarise(
          Mean_F1 = mean(F1, na.rm = TRUE),
          SD_F1 = sd(F1, na.rm = TRUE),
          Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
          SD_Accuracy = sd(Accuracy, na.rm = TRUE),
          Mean_AUC = mean(AUC, na.rm = TRUE),
          SD_AUC = sd(AUC, na.rm = TRUE),
          Mean_Precision = mean(Precision, na.rm = TRUE),
          Mean_Recall = mean(Recall, na.rm = TRUE),
          n_Datasets = n()
        ) %>%
        arrange(desc(Mean_F1))
      
      summary_file <- file.path(output_folder, "method_summary_across_datasets.csv")
      write.csv(method_summary, summary_file, row.names = FALSE)
      cat(sprintf("   Method summary saved: %s\n", summary_file))
      
      # Generate comprehensive report
      report_file <- file.path(output_folder, "batch_processing_report.txt")
      sink(report_file)
      
      cat("DA Method Batch Analysis Comprehensive Report\n")
      cat(paste("=" %>% rep(50), collapse = ""), "\n\n")
      
      cat(sprintf("Analysis time: %s\n", Sys.time()))
      cat(sprintf("Total datasets: %d\n", length(data_files)))
      cat(sprintf("Successfully analyzed: %d\n", success_count))
      cat(sprintf("Failed: %d\n", length(failed_datasets)))
      cat("\n")
      
      cat("Best Methods (Average F1 Score):\n")
      top_n <- min(5, nrow(method_summary))
      for(i in 1:top_n) {
        m <- method_summary[i, ]
        cat(sprintf("  %d. %s: F1=%.3f ± %.3f, Accuracy=%.3f ± %.3f (based on %d datasets)\n",
                    i, m$Method, m$Mean_F1, m$SD_F1, 
                    m$Mean_Accuracy, m$SD_Accuracy, m$n_Datasets))
      }
      
      sink()
      cat(sprintf("   Comprehensive report saved: %s\n", report_file))
    }
  }
  
  cat("\nAll processing completed!\n")
  cat(sprintf("Results saved in: %s\n", output_folder))
  cat(sprintf("  - Metrics Results: %s\n", metrics_folder))
  cat(sprintf("  - Visualizations: %s\n", plots_folder))
  cat(sprintf("  - PDF Graphs: %s\n", pdf_folder))
  
  return(list(
    all_results = all_results,
    success_count = success_count,
    failed_datasets = failed_datasets
  ))
}

# ============================================================================
# 6. Execute Batch Processing
# ============================================================================

# Set random seed for reproducibility
set.seed(123)

# Start batch processing
final_results <- batch_process_all_datasets(input_folder)

# Clean workspace
cat("\nCleaning workspace...\n")
objects_to_keep <- c("final_results", "input_folder", "output_folder", 
                     "metrics_folder", "plots_folder", "pdf_folder",
                     "batch_process_all_datasets", "process_single_dataset")
objects_to_remove <- setdiff(ls(), objects_to_keep)
if(length(objects_to_remove) > 0) {
  rm(list = objects_to_remove)
}
gc()

cat("\nScript execution completed!\n")



######## Comprehensive Interpretation of Results:
# First, read the merged data
combined_data <- read.csv("/Users/xiaoying/Desktop/DA_Results_Analysis/21_Datasets_Accuracy_Interpretation/all_datasets_combined_metrics.csv")

# View method names
cat("Original method names (may have inconsistent casing):\n")
unique(combined_data$Method)

# Standardize method names (convert to lowercase and then capitalize first letter)
standardize_method_names <- function(method_names) {
  # Remove whitespace
  cleaned <- trimws(method_names)
  
  # Convert to lowercase
  lower_case <- tolower(cleaned)
  
  # Common method standardization mapping
  method_mapping <- c(
    "sccomp" = "sccomp",
    "scCODA" = "scCODA", 
    "scanpro" = "scanpro",
    "Proprller" = "Propeller",
    "louvain" = "Louvain",
    "louvain_clustering" = "Louvain",
    "dcats" = "DCATS",
    "dirichletreg" = "DirichletReg",
    "dirichlet_regression" = "DirichletReg"
  )
  
  # Apply mapping
  standardized <- method_mapping[lower_case]
  # If no match, keep original (but standardize case)
  standardized[is.na(standardized)] <- tools::toTitleCase(lower_case[is.na(standardized)])
  
  return(standardized)
}

# Apply standardization
combined_data$Method_Standardized <- standardize_method_names(combined_data$Method)

# View standardized names
cat("\nStandardized method names:\n")
print(table(combined_data$Method_Standardized))

# Save standardized data
write.csv(combined_data, 
          "/Users/xiaoying/Desktop/DA_Results_Analysis/21_Datasets_Accuracy_Interpretation/all_datasets_combined_metrics_standardized.csv",
          row.names = FALSE)






# Comprehensive Evaluation: -------------------------------------------------------------------
# Fixed calculation function
calculate_comprehensive_scores_fixed <- function(data) {
  # Basic statistics summary
  basic_stats <- data %>% 
    group_by(Method_Standardized) %>%
    summarise(
      # Basic statistics
      n_Datasets = n(),
      Mean_F1 = mean(F1, na.rm = TRUE),
      Median_F1 = median(F1, na.rm = TRUE),
      SD_F1 = sd(F1, na.rm = TRUE),
      CV_F1 = ifelse(Mean_F1 > 0, SD_F1 / Mean_F1 * 100, NA),
      
      Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
      Mean_AUC = mean(AUC, na.rm = TRUE),
      Mean_Precision = mean(Precision, na.rm = TRUE),
      Mean_Recall = mean(Recall, na.rm = TRUE),
      Mean_MCC = mean(MCC, na.rm = TRUE),
      Mean_Balanced_Accuracy = mean(Balanced_Accuracy, na.rm = TRUE),
      
      # Stability metrics
      Min_F1 = min(F1, na.rm = TRUE),
      Max_F1 = max(F1, na.rm = TRUE),
      Range_F1 = Max_F1 - Min_F1,
      IQR_F1 = IQR(F1, na.rm = TRUE),
      
      # Success rate
      Success_Rate_F1_0.6 = sum(F1 >= 0.6, na.rm = TRUE) / n_Datasets * 100,
      Success_Rate_F1_0.7 = sum(F1 >= 0.7, na.rm = TRUE) / n_Datasets * 100,
      Success_Rate_F1_0.8 = sum(F1 >= 0.8, na.rm = TRUE) / n_Datasets * 100,
      
      # Consistency
      Mean_Consistency = mean(Consistency_Rate, na.rm = TRUE)
    )
  
  # Calculate ranking statistics
  cat("Calculating ranking statistics...\n")
  
  # Prepare ranking data
  rank_data_list <- list()
  datasets <- unique(data$Dataset)
  
  for(dataset in datasets) {
    dataset_data <- data %>% filter(Dataset == dataset)
    
    if(nrow(dataset_data) > 0) {
      # F1 ranking (higher value is better)
      f1_ranks <- rank(-dataset_data$F1, ties.method = "average", na.last = "keep")
      # Accuracy ranking
      acc_ranks <- rank(-dataset_data$Accuracy, ties.method = "average", na.last = "keep")
      
      temp_df <- data.frame(
        Method_Standardized = dataset_data$Method_Standardized,
        Dataset = dataset,
        Rank_F1 = f1_ranks,
        Rank_Accuracy = acc_ranks,
        stringsAsFactors = FALSE
      )
      rank_data_list[[dataset]] <- temp_df
    }
  }
  
  # Combine all ranking data
  if(length(rank_data_list) > 0) {
    rank_data <- do.call(rbind, rank_data_list)
    
    rank_summary <- rank_data %>%
      group_by(Method_Standardized) %>%
      summarise(
        Mean_Rank_F1 = mean(Rank_F1, na.rm = TRUE),
        Mean_Rank_Accuracy = mean(Rank_Accuracy, na.rm = TRUE),
        Best_Rank_F1 = min(Rank_F1, na.rm = TRUE),
        Worst_Rank_F1 = max(Rank_F1, na.rm = TRUE),
        Rank_Stability_F1 = sd(Rank_F1, na.rm = TRUE)
      )
  } else {
    # If no ranking data, create an empty dataframe
    rank_summary <- data.frame(
      Method_Standardized = unique(data$Method_Standardized),
      Mean_Rank_F1 = NA,
      Mean_Rank_Accuracy = NA,
      Best_Rank_F1 = NA,
      Worst_Rank_F1 = NA,
      Rank_Stability_F1 = NA
    )
  }
  
  # Combine basic statistics and ranking statistics
  combined_stats <- basic_stats %>%
    left_join(rank_summary, by = "Method_Standardized")
  
  # Calculate various composite scores - fixed vectorization issue
  final_result <- combined_stats %>%
    mutate(
      # 1. Weighted average score
      Composite_Score_Weighted = (0.4 * Mean_F1 + 0.2 * Mean_Accuracy + 
                                    0.15 * Mean_AUC + 0.1 * Mean_MCC + 
                                    0.1 * Mean_Balanced_Accuracy + 0.05 * Mean_Consistency),
      
      # 2. Stability-adjusted score - fixed vectorization issue
      Composite_Score_Stable = case_when(
        is.na(CV_F1) ~ Mean_F1,
        CV_F1 == 0 ~ Mean_F1,
        TRUE ~ Mean_F1 * (1 - pmin(CV_F1, 100)/200)
      ),
      
      # 3. Geometric mean score
      Composite_Score_Geometric = apply(
        cbind(
          pmax(Mean_F1, 0.01), 
          pmax(Mean_Accuracy, 0.01), 
          pmax(Mean_AUC, 0.01),
          pmax((Mean_MCC + 1)/2, 0.01)
        ), 
        1, 
        function(x) exp(mean(log(x), na.rm = TRUE))
      ),
      
      # 4. Success rate weighted score
      Composite_Score_Success = (Success_Rate_F1_0.6 * 0.3 + 
                                   Success_Rate_F1_0.7 * 0.4 + 
                                   Success_Rate_F1_0.8 * 0.3) / 100
    )
  
  # 5. Ranking score (considering rank and stability)
  final_result$Composite_Score_Rank <- mapply(
    function(rank, stability) {
      if(is.na(rank) || is.na(stability)) return(NA)
      (1 / rank) * (1 - pmin(stability, 10)/20)
    },
    final_result$Mean_Rank_F1,
    final_result$Rank_Stability_F1
  )
  
  # Standardize scores (so that they are on the same scale)
  # Only standardize non-NA columns
  score_columns <- c("Composite_Score_Weighted", "Composite_Score_Stable", 
                     "Composite_Score_Geometric", "Composite_Score_Success", 
                     "Composite_Score_Rank")
  
  # Check which columns have valid data
  valid_columns <- score_columns[sapply(score_columns, function(col) {
    any(!is.na(final_result[[col]]))
  })]
  
  # If there are valid score columns, standardize them
  if(length(valid_columns) > 0) {
    for(col in valid_columns) {
      scaled_col <- paste0(col, "_Scaled")
      # Only standardize non-NA values
      non_na_idx <- !is.na(final_result[[col]])
      if(sum(non_na_idx) > 1) {
        scaled_values <- rep(NA, nrow(final_result))
        scaled_values[non_na_idx] <- scale(final_result[[col]][non_na_idx])[,1]
        final_result[[scaled_col]] <- scaled_values
      } else {
        final_result[[scaled_col]] <- NA
      }
    }
    
    # Calculate final composite score (using available standardized scores)
    scaled_cols <- paste0(valid_columns, "_Scaled")
    
    # Calculate row mean, ignoring NA values
    final_result$Final_Composite_Score <- apply(
      final_result[, scaled_cols, drop = FALSE], 
      1, 
      function(x) mean(x, na.rm = TRUE)
    )
  } else {
    # If no valid scores, use the weighted score
    final_result$Final_Composite_Score <- final_result$Composite_Score_Weighted
  }
  
  # Sort by final score
  final_result <- final_result %>%
    arrange(desc(Final_Composite_Score))
  
  return(final_result)
}

# Run the fixed function
cat("Calculating comprehensive scores...\n")
comprehensive_scores <- calculate_comprehensive_scores_fixed(combined_data)

# View the results
cat("\n=== Comprehensive Scores Results ===\n")
print(comprehensive_scores %>% 
        select(Method_Standardized, n_Datasets, Mean_F1, SD_F1, CV_F1,
               Final_Composite_Score, Success_Rate_F1_0.7, Mean_Rank_F1) %>%
        arrange(desc(Final_Composite_Score)))

write.csv(comprehensive_scores, 
          "/Users/xiaoying/Desktop/DA结果分析/21个双条件数据集准确性解读结果/comprehensive_scores_final.csv",
          row.names = FALSE)




# --------------------------------------------------------------------

visualize_simple_scores <- function(scores_df) {   
  library(ggplot2)   
  library(ggrepel)   
  library(patchwork)      
  
  # 1. Composite Score Bar Chart   
  p1 <- scores_df %>%     
    mutate(Method = reorder(Method_Standardized, Final_Composite_Score)) %>%     
    ggplot(aes(x = Method, y = Final_Composite_Score, fill = Final_Composite_Score)) +     
    geom_bar(stat = "identity", color = "gray", linewidth = 0.3) +     
    geom_text(aes(label = sprintf("%.3f", Final_Composite_Score)),                
              hjust = -0.1, size = 3.5, fontface = "bold") +     
    coord_flip() +     
    scale_fill_gradient(low = "#FFE5CC", high = "#FF6600", name = "Composite Score") +     
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +     
    labs(title = "DA Methods Composite Score Ranking",          
         subtitle = "Cross-dataset Evaluation based on 21 datasets",          
         x = "Method", y = "Composite Score") +     
    theme_minimal() +     
    theme(       
      panel.grid.major = element_blank(),       
      panel.grid.minor = element_blank(),       
      panel.border = element_rect(fill = NA, color = "gray", linewidth = 0.7),       
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),       
      axis.title = element_text(face = "bold"),       
      legend.position = "none"     
    )      
  
  # 2. Multi-metric Scatter Matrix   
  # Select key metrics   
  scatter_data <- scores_df %>%     
    select(Method_Standardized,             
           "F1 Score" = Mean_F1,            
           "Accuracy" = Mean_Accuracy,            
           "AUC" = Mean_AUC,            
           "MCC" = Mean_MCC,            
           "Success Rate" = Success_Rate_F1_0.7,            
           "Composite Score" = Final_Composite_Score)      
  
  # Create scatter plot matrix   
  p2 <- ggplot(scatter_data, aes(x = `F1 Score`, y = `Accuracy`)) +     
    geom_point(aes(size = `Composite Score`, color = `Success Rate`), alpha = 0.8) +     
    geom_text_repel(aes(label = Method_Standardized), size = 3.2, max.overlaps = 15) +     
    scale_color_gradient(low = "blue", high = "red", name = "Success Rate (%)") +     
    scale_size_continuous(range = c(3, 10), name = "Composite Score") +     
    labs(title = "F1 Score vs Accuracy",          
         subtitle = "Point size represents composite score, color represents success rate",          
         x = "Mean F1 Score", y = "Mean Accuracy") +     
    theme_minimal() +     
    theme(       
      panel.grid.major = element_blank(),       
      panel.grid.minor = element_blank(),       
      panel.border = element_rect(fill = NA, color = "gray", linewidth = 0.7),       
      plot.title = element_text(hjust = 0.5, face = "bold")     
    )      
  
  # 3. Radar Chart Data Preparation   
  radar_data <- scores_df %>%     
    select(Method_Standardized,             
           `F1 Score` = Mean_F1,            
           `Accuracy` = Mean_Accuracy,            
           `AUC` = Mean_AUC,            
           `Balanced Accuracy` = Mean_Balanced_Accuracy,            
           `Consistency` = Mean_Consistency) %>%     
    mutate(across(-Method_Standardized, ~ (. - min(.)) / (max(.) - min(.))))      
  
  # Convert to long format   
  radar_long <- radar_data %>%     
    pivot_longer(cols = -Method_Standardized, names_to = "Metric", values_to = "Value")      
  
  # Only display top 5 methods   
  top_methods <- head(scores_df$Method_Standardized, 5)   
  radar_top <- radar_long %>% filter(Method_Standardized %in% top_methods)      
  
  p3 <- ggplot(radar_top, aes(x = Metric, y = Value, group = Method_Standardized,                                
                              color = Method_Standardized)) +     
    geom_line(linewidth = 1, alpha = 0.7) +     
    geom_point(size = 2) +     
    coord_polar() +     
    scale_y_continuous(limits = c(0, 1)) +     
    labs(title = "Top 5 Methods Performance Radar Chart",          
         subtitle = "Normalized comparison across metrics",          
         x = "", y = "") +     
    theme_minimal() +     
    theme(       
      axis.text.x = element_text(angle = 0, hjust = 1),       
      panel.grid.major = element_line(color = "gray90"),       
      panel.border = element_rect(fill = NA, color = "gray", linewidth = 0.7),       
      plot.title = element_text(hjust = 0.5, face = "bold"),       
      legend.position = "right"     
    )      
  
  # 4. Success Rate Bar Chart   
  p4 <- scores_df %>%     
    mutate(Method = reorder(Method_Standardized, Success_Rate_F1_0.7)) %>%     
    ggplot(aes(x = Method, y = Success_Rate_F1_0.7)) +     
    geom_bar(stat = "identity", fill = "#66C2A5", color = "gray", linewidth = 0.2) +     
    geom_text(aes(label = sprintf("%.1f%%", Success_Rate_F1_0.7)),                
              hjust = -0.1, size = 3.2, fontface = "bold") +     
    coord_flip() +     
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +     
    labs(title = "Success Rate Analysis (F1 ≥ 0.7)",          
         subtitle = "Proportion of methods achieving good performance across 21 datasets",          
         x = "Method", y = "Success Rate (%)") +     
    theme_minimal() +     
    theme(       
      panel.grid.major = element_blank(),       
      panel.grid.minor = element_blank(),       
      panel.border = element_rect(fill = NA, color = "gray", linewidth = 0.7),       
      plot.title = element_text(hjust = 0.5, face = "bold")     
    )      
  
  # Combine Plots   
  combined_plot <- (p1 | p2) / (p3 | p4) +     
    plot_annotation(       
      title = "DA Methods Comprehensive Performance Evaluation",       
      subtitle = "Comprehensive Analysis across 21 datasets",       
      theme = theme(         
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),         
        plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray40")       
      )     
    )      
  
  # Display Plot   
  print(combined_plot)      
  
  # Save Plot   
  ggsave("/Users/xiaoying/Desktop/DA_Results_Analysis/21_Datasets_Accuracy_Interpretation/comprehensive_analysis.png",           
         combined_plot, width = 16, height = 12, dpi = 300, bg = "white")      
  
  ggsave("/Users/xiaoying/Desktop/DA_Results_Analysis/21_Datasets_Accuracy_Interpretation/comprehensive_analysis.pdf",           
         combined_plot, width = 10, height = 8, device = "pdf", bg = "white")      
  
  cat("Comprehensive Evaluation Plot Saved!\n")      
  
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, combined = combined_plot)) 
}

# Generate Visualization 
plots <- visualize_simple_scores(comprehensive_scores)
