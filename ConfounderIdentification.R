# ------------------------------------------------------------------------------
# Script:        ConfounderIdentification.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Functions for the identification of Confounders
# ------------------------------------------------------------------------------


################################################################################
# Confounder identification
################################################################################
confounderIdentificationMenu <- function() {
  # Run preliminary variable setup
  vars <- setup_variables()
  if (is.null(vars)) return(invisible(NULL))
  ptx   <- vars$ptx
  sinfo <- vars$sinfo
  binfo <- vars$binfo
  
  repeat {
    cat("\n=== Confounder Identification Menu: ===\n")
    cat("  1 = Boxplots\n")
    cat("  2 = Covariate-colored PCA-plot\n")
    cat("  3 = Linear Regression\n")
    cat("  4 = Surrogate Variable Analysis\n")
    cat("  0 = Exit Confounder Identification\n")
    confound_choice <- read_choice("Enter your choice [0-4]: ")
    
    if (is.na(confound_choice) || confound_choice == 0) {
      cat("\nExiting Confounder Identification.\n")
      break
    } else if (confound_choice == 1) {
      confounderBoxplots(sinfo)
    } else if (confound_choice == 2) {
      confoundersPCA(ptx, sinfo)
    } else if (confound_choice == 3) {
      prop_shift <- as.numeric(readline(prompt = "Enter cutoff value for large proportional effect:\n"))
      DEAconfoundersX(expr = NULL, sinfo = NULL, prop_shift = prop_shift)
    } else if (confound_choice == 4)
      runSVA()
  }
}

################################################################################
# Function for boxplots
################################################################################
confounderBoxplots <- function(sinfo) { 
  cat("\nAvailable columns in sinfo:\n")
  print(names(sinfo))
  col1 <- trimws(readline("First column: "))
  col2 <- trimws(readline("Second column: "))
  if ((col1 %in% names(sinfo)) && (col2 %in% names(sinfo))) {
    p <- ggplot(sinfo, aes_string(x = col1, y = col2)) +
      geom_boxplot() +
      labs(title = paste("Boxplot of", col1, "vs", col2),
           x = col1, y = col2) +
      theme_minimal()
    print(p)
    cat("\nBoxplot created.\n")
  } else {
    cat("\nInvalid column selection. Please check the names and try again.\n")
  }
}

################################################################################
# Function for pca plot colored by selected variable (column) 
################################################################################
confoundersPCA <- function(ptx,sinfo) {
  cat("\nRunning PCA for Sample Batch Analysis...\n")
  # Run PCA (Assuming data has already been scaled, else change scale. argument to TRUE)
  pca <- prcomp(as.data.frame(ptx), scale. = FALSE)
  pcaX <- as.data.frame(pca$x, row.names = rownames(ptx))
  cat("\nAvailable columns in sinfo:\n")
  print(names(sinfo))
  
  # Prompt for column to use for coloring
  batch_color <- trimws(readline("Please enter the column name to use for color: "))
  
  # Create plot 
  if (batch_color %in% names(sinfo)) {
    color_values <- sinfo[match(rownames(ptx), rownames(sinfo)), batch_color]
    pcaX$Color <- as.factor(color_values)
    p <- ggplot(pcaX, aes(x = PC1, y = PC2, color = Color)) +
      geom_point(size = 3) +
      labs(title = "PCA for Sample Batch Analysis",
           x = paste0("PC1: ", round((pca$sdev[1]^2 / sum(pca$sdev^2))*100, 2), " %"),
           y = paste0("PC2: ", round((pca$sdev[2]^2 / sum(pca$sdev^2))*100, 2), " %")) +
      theme_minimal()
  } else {
    cat("Invalid column for coloring; proceeding without color.\n")
    p <- ggplot(pcaX, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) +
      labs(title = "PCA for Sample Batch Analysis",
           x = paste0("PC1: ", round((pca$sdev[1]^2 / sum(pca$sdev^2))*100, 2), " %"),
           y = paste0("PC2: ", round((pca$sdev[2]^2 / sum(pca$sdev^2))*100, 2), " %")) +
      theme_minimal()
  }
  print(p)
  cat("\nPCA for Sample Batch Analysis completed.\n")
}

################################################################################
# Linear regression for confounder identification
################################################################################
DEAconfoundersX <- function(expr = NULL, sinfo = NULL, prop_shift) {
  # Check required objects
  if (is.null(expr) && !exists("select.ptx", envir = .GlobalEnv))
    stop("Expression data not provided and 'select.ptx' not found.")
  if (is.null(sinfo) && !exists("select.sinfo", envir = .GlobalEnv))
    stop("Sample info not provided and 'select.sinfo' not found.")
  
  if (is.null(expr)) expr <- get("select.ptx", envir = .GlobalEnv)
  if (is.null(sinfo)) sinfo <- get("select.sinfo", envir = .GlobalEnv)
  
  # Transpose expression matrix as limma requires features as rows
  expr <- t(expr)
  
  # Prompt for main variable
  cat("\nAvailable clinical variables:\n")
  for (i in seq_along(colnames(sinfo))) {
    cat(i, ":", colnames(sinfo)[i], "\n")
  }
  main_idx <- as.numeric(readline("Select the main variable of interest: "))
  if (is.na(main_idx) || main_idx < 1 || main_idx > ncol(sinfo)) stop("Invalid selection.")
  main_var <- colnames(sinfo)[main_idx]
  cat("Main variable selected:", main_var, "\n")
  
  # Prompt for candidate confounders 
  cat("\nSelect candidate confounders (comma or space separated indices):\n")
  input <- readline("Enter numbers: ")
  nums <- as.integer(unlist(strsplit(input, "[, ]+")))
  nums <- nums[!is.na(nums) & nums >= 1 & nums <= ncol(sinfo)]
  confounders <- setdiff(colnames(sinfo)[nums], main_var)
  if (length(confounders) == 0) stop("No valid confounders selected.")
  cat("Confounders selected:", paste(confounders, collapse = ", "), "\n")
  
  # Fit base model 
  design_base <- model.matrix(as.formula(paste("~", main_var)), data = sinfo)
  fit_base <- limma::eBayes(limma::lmFit(expr, design_base))
  coef_main <- grep(paste0("^", main_var), colnames(design_base), value = TRUE)
  base_tab <- limma::topTable(fit_base, coef = coef_main, number = Inf, sort.by = "none")
  base_logFC <- base_tab$logFC
  
  # Loop over chosen confounders
  results <- data.frame(
    confounder = character(),
    mean_abs_shift = numeric(),
    prop_large_shift = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (conf in confounders) {
    design_ext <- model.matrix(as.formula(paste("~", main_var, "+", conf)), data = sinfo)
    fit_ext <- limma::eBayes(limma::lmFit(expr, design_ext))
    coef_main_ext <- grep(paste0("^", main_var), colnames(design_ext), value = TRUE)
    ext_tab <- limma::topTable(fit_ext, coef = coef_main_ext, number = Inf, sort.by = "none")
    ext_logFC <- ext_tab$logFC
    
    delta <- ext_logFC - base_logFC
    results <- rbind(results, data.frame(
      confounder = conf,
      mean_abs_shift = mean(abs(delta), na.rm = TRUE),
      prop_large_shift = mean(abs(delta) > prop_shift, na.rm = TRUE)
    ))
  }
  
  # Rank confounders and present summary
  results <- results[order(-results$mean_abs_shift), ]
  cat("\n--- Confounder impact summary ---\n")
  print(results, row.names = FALSE)
  
  invisible(results)
}

################################################################################
# Function for Surrogate Variable Analysis (SVA)
################################################################################
runSVA <- function() {
  # Check that select.ptx exists.
  if (!exists("select.ptx") || is.null(select.ptx)) {
    cat("\nError: 'select.ptx' is not available. Ensure the expression data is loaded before SVA.\n")
    return(invisible(NULL))
  }
  
  # Check that select.sinfo exists.
  if (!exists("select.sinfo") || is.null(select.sinfo)) {
    cat("\nError: 'select.sinfo' is not available. Clinical information is needed for SVA.\n")
    return(invisible(NULL))
  }
  
  # Check that the number of samples (rows) in select.ptx matches that in select.sinfo.
  if (nrow(select.ptx) != nrow(select.sinfo)) {
    cat("\nError: The number of rows in 'select.ptx' (samples) must equal the number of rows in 'select.sinfo'.\n")
    return(invisible(NULL))
  }
  
  cat("\nPerforming Surrogate Variable Analysis (SVA)...\n")
  
  # Build the Full Model Matrix
  # List available columns in select.sinfo in a numbered list.
  cat("Available columns in select.sinfo:\n")
  for (i in 1:ncol(select.sinfo)) {
    cat(i, ": ", colnames(select.sinfo)[i], "\n")
  }
  
  # Prompt for column selection by numbers.
  group_input <- readline(prompt = "Enter the number(s) corresponding to the column(s) to include as known factors (comma separated; leave blank for none): ")
  group_input <- trimws(group_input)
  
  if (group_input == "") {
    cat("No known factors specified. Using only the intercept in the model.\n")
    mod <- model.matrix(~ 1, data = select.sinfo)
  } else {
    # Split and convert to numeric.
    indices <- as.numeric(unlist(strsplit(group_input, split = ",")))
    
    # Validate indices.
    if (any(is.na(indices)) || any(indices < 1) || any(indices > ncol(select.sinfo))) {
      cat("Error: Invalid selection(s).\n")
      return(invisible(NULL))
    }
    
    # Get the corresponding column names.
    selected_vars <- colnames(select.sinfo)[indices]
    
    # Build the design matrix using the selected variables.
    formula_str <- paste("~", paste(selected_vars, collapse = " + "))
    mod <- model.matrix(as.formula(formula_str), data = select.sinfo)
  }
  
  # Construct the null model (only intercept).
  mod0 <- model.matrix(~ 1, data = select.sinfo)
  
  # Prompt for Number of Surrogate Variables
  nsv_input <- readline(prompt = "Enter number of surrogate variables to estimate (leave empty for automatic estimation): ")
  nsv_input <- trimws(nsv_input)
  
  if (nsv_input == "") {
    nsv <- NULL
  } else {
    nsv <- as.numeric(nsv_input)
    if (is.na(nsv) || nsv < 1) {
      cat("Invalid number provided. Using automatic estimation.\n")
      nsv <- NULL
    }
  }
  
  # Run SVA 
  sva_result <- sva(t(select.ptx), mod, mod0, n.sv = nsv)
  cat("\nSVA analysis complete. Estimated number of surrogate variables: ", sva_result$n.sv, "\n")
  
  # Add surrogate variable columns to select.sinfo
  sinfo_sva <- cbind(select.sinfo, sva_result$sv)
  
  # Rename the surrogate variable columns for clarity.
  surrogate_names <- paste0("SV", 1:sva_result$n.sv)
  colnames(sinfo_sva)[(ncol(select.sinfo) + 1):ncol(sinfo_sva)] <- surrogate_names
  
  # Plot the Surrogate Variables (scatter-plots)
  cat("Plotting surrogate variables...\n")
  if (sva_result$n.sv > 1) {
    pairs(sva_result$sv, main = "Pairwise Scatterplots of Surrogate Variables")
  } else {
    hist(sva_result$sv, main = "Distribution of Surrogate Variable", xlab = "SV1")
  }
  
  # Assign results to global environment and return
  assign("sva_result", sva_result, envir = .GlobalEnv)
  assign("sinfo_sva", sinfo_sva, envir = .GlobalEnv)
  bup.sinfo <<- select.sinfo
  select.sinfo <<- sinfo_sva
  
  cat("Results stored in the global environment as 'sva_result' and 'sinfo_sva'.\n")
  cat("Surrogate Variables added to select.sinfo.\n")
  cat("sinfo without surrogate variables stored in bup.sinfo.\n")
  return(list(sva_result = sva_result, sinfo_sva = sinfo_sva))
}


