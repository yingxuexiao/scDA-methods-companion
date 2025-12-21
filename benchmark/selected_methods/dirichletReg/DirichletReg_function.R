library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
library(DirichletReg)
library(reshape2)
library(tibble)

run_da_dirichletreg <- function(meta_data, sample_col, condition_col, celltype_col) {
  all_results <- list()
  
  # Ensure required columns are present
  stopifnot(all(c(sample_col, condition_col, celltype_col) %in% colnames(meta_data)))
  
  # ==== Construct df table ====
  df <- meta_data %>%
  rename(
    sample    = !!sym(sample_col),
    condition = !!sym(condition_col),
    celltype  = !!sym(celltype_col)
  ) %>%
  mutate(
    test = as.character(celltype)  # Ensure it's a character type
  )

  # ==== Construct df_summary (percentage table) ====
  df_summary <- df %>%
  group_by(sample, condition, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample, condition) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup() %>%
  mutate(percent = as.numeric(percent))  # Ensure percentage is numeric

  # ==== Construct condition pair combinations ====
  unique_conditions <- unique(df$condition)
  if (length(unique_conditions) > 2) {
    pairs_list <- as.data.frame(t(combn(unique_conditions, 2)), stringsAsFactors = FALSE)
    colnames(pairs_list) <- c("state_1", "state_2")
  } else if (length(unique_conditions) == 2) {
    pairs_list <- data.frame(state_1 = unique_conditions[1], state_2 = unique_conditions[2])
  } else {
    stop("Need at least two conditions for comparison.")
  }
  
  # ==== Loop through each comparison pair ====
  for (i in seq_len(nrow(pairs_list))) {
    state_1 <- as.character(pairs_list$state_1[i])
    state_2 <- as.character(pairs_list$state_2[i])
    
    message("\nProcessing comparison: ", state_1, " vs ", state_2)
    
    metaDF_tmp <- df_summary %>% filter(condition %in% c(state_1, state_2))
    if (nrow(metaDF_tmp) == 0) {
      warning("No data for ", state_1, " vs ", state_2, ". Skipping.")
      next
    }
    
    # Construct rowname_ID
    metaDF_tmp <- metaDF_tmp %>% mutate(rowname_ID = paste(sample, condition, sep = "|"))
    
    # Reshape data to wide format
    dummy_mat <- tryCatch(
      {
        reshape2::dcast(
          metaDF_tmp,
          rowname_ID ~ celltype,
          value.var = "percent",
          fill = 0
        ) %>%
          tibble::column_to_rownames("rowname_ID") %>%
          as.data.frame() %>%
          mutate(across(everything(), as.numeric)) %>%
          as.matrix()
      },
      error = function(e) {
        warning("Error reshaping data for ", state_1, " vs ", state_2, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(dummy_mat)) next
    
    # Split rowname_ID
    rowname_ids <- rownames(dummy_mat)
    split_parts <- strsplit(rowname_ids, "\\|")
    
    # Construct model data frame
    model_data <- data.frame(
      rowname_ID = rowname_ids,
      sample = sapply(split_parts, `[`, 1),
      condition = factor(sapply(split_parts, `[`, 2), levels = c(state_1, state_2)),
      dummy_mat,
      check.names = FALSE
    )
    
    # Construct DirichletReg object
    model_data$counts <- DirichletReg::DR_data(dummy_mat)
    
    # Fit Dirichlet regression
    fit <- tryCatch(
      {
        DirichletReg::DirichReg(counts ~ condition, data = model_data, model = "common"
        #,control=list(interlim = 2000,tol1=1e-2,tol2=1e-2)
        )
      },
      error = function(e) {
        warning("Error fitting model for ", state_1, " vs ", state_2, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(fit)) next
    
    # Extract results
    fit_summary <- summary(fit)
    celltypes <- rep(fit_summary$varnames, each = 2)
    param_names <- rownames(fit_summary$coef.mat)
    
    results_table <- data.frame(
      celltype = celltypes,
      parameter = param_names,
      Estimate = fit_summary$coef.mat[, "Estimate"],
      Std.Error = fit_summary$coef.mat[, "Std. Error"],
      z.value = fit_summary$coef.mat[, "z value"],
      p.value = fit_summary$coef.mat[, "Pr(>|z|)"],
      condition_1 = state_1,
      condition_2 = state_2
    )
    
    all_results[[paste(state_1, state_2, sep = "_vs_")]] <- results_table
  }
  
  if (length(all_results) == 0) {
    warning("No successful DA results.")
    da_results <- NULL
  } else {
    da_results <- dplyr::bind_rows(all_results)
  }
  
  return(list(
    df = df,                  # Original table
    df_summary = df_summary,  # Percentage table
    da_results = da_results   # DirichletReg fit results
  ))
}
