# ------------------------------------------------------------------------------
# Script:        PreProcessing.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Functions for pre-processing of Protein Abundance Data
# ------------------------------------------------------------------------------

################################################################################
# Pre-processing functions:
################################################################################

################################################################################
# Pre-processing Main menu:
################################################################################
preprocess_main_menu <- function() {
  # First, call pre-existing variable setup function.
  setup_variables()
  cat("\nVariables set up successfully.\n")
  
  # Main interactive menu loop
  repeat {
    cat("\n=== Preprocessing Menu: ===\n")
    cat("  1 = Variance Filtering\n")
    cat("  2 = Normalization and Scaling\n")
    cat("  3 = Correlation Analysis\n")
    cat("  4 = Confounder Identification Menu\n")
    cat("  5 = Adjustment for Confounders and Batch Effects\n")
    cat("  6 = Dimensionality Reduction\n")
    cat("  7 = Inspect Loadings\n")
    cat("  0 = Exit\n")
    
    main_choice <- as.numeric(readline(prompt = "Enter your choice [0-7]:"))
    
    if (main_choice == 1) {
      filter_data(reproduce = NULL)           
    } else if (main_choice == 2) {
      normalize_data(reproduce = NULL)        
    } else if (main_choice == 3) {
      analyze_correlation(reproduce = NULL)   
    } else if (main_choice == 4) {
      confounderIdentificationMenu()   
    } else if (main_choice == 5) {
      AdjustForConfounderAndBatchEffects(reproduce = NULL)
    } else if (main_choice == 6) {
      dimReductionMenu()  
    } else if (main_choice == 7){
      loading_inspector()
    } else if (main_choice == 0) {
      cat("\nExiting the script.\n")
      break
    } else {
      cat("\nInvalid selection. Please try again.\n")
    }
  }
}

################################################################################
# Min-Max normalization:
################################################################################
min_max_normalize <-function(x) {
  return((x-min(x)) / (max(x) - min(x)))
}

################################################################################
# Inverse Normal Transformation:
################################################################################
inverse_normalize_df <- function(df, ties.method = "average") {
  # Check that df is a dataframe
  if (!is.data.frame(df)) {
    stop("Input must be a dataframe.")
  }
  
  # Function to transform a numeric vector using rank-based inverse normal transformation
  inv_norm <- function(x, ties.method = "average") {
    if (!is.numeric(x)) {
      warning("Non-numeric column encountered; returning unchanged column.")
      return(x)
    }
    
    # Count non-missing values
    N <- sum(!is.na(x))
    if (N == 0) return(x)  # If all values are NA, leave the column unchanged
    
    # Compute ranks with NA preserved
    r <- rank(x, na.last = "keep", ties.method = ties.method)
    
    # Map the scaled ranks to standard normal quantiles using qnorm
    transformed <- qnorm((r - 0.5) / N)
    return(transformed)
  }
  
  # Apply the transformation to every column and return as a dataframe
  transformed_df <- data.frame(lapply(df, inv_norm, ties.method = ties.method), 
                               row.names = rownames(df))
  assign("select.ptx", transformed_df, .GlobalEnv)
  return(transformed_df)
}

################################################################################
# Variance filtering:
################################################################################
filter_data <- function(reproduce = NULL) {
  if (!is.null(reproduce)) {
    variance_choice <- reproduce$varChoice
    percentile_input <- reproduce$percInput
    
    if (variance_choice == 1) {
      # Column Variance Filtering
      col_vars <- apply(select.ptx, 2, var)
      cutoff <- quantile(col_vars, percentile_input)
      select.ptx <<- as.data.frame(select.ptx[, col_vars > cutoff, drop = FALSE])
    } else if (variance_choice == 2) {
      # Row Variance Filtering
      row_vars <- apply(select.ptx, 1, var)
      cutoff <- quantile(row_vars, percentile_input)
      select.ptx <<- as.data.frame(select.ptx[row_vars > cutoff, , drop = FALSE])
      select.sinfo <<- select.sinfo[rownames(select.sinfo) %in% rownames(select.ptx), ]
    } else {
      reproduce <- data.frame(
        varChoice = NA,
        percInput = NA,
        stringsAsFactors = FALSE
      )
      
      cat("\n--- Variance Filtering ---\n")
      cat("Choose filtering option:\n")
      cat("  1 = Column Variance Filtering\n")
      cat("  2 = Row Variance Filtering\n")
      cat("  0 = Back/Skip\n")
      
      variance_choice <- as.numeric(readline(prompt = "Enter your choice: "))
      reproduce$varChoice <- variance_choice
      
      if (variance_choice == 0) {
        cat("\nReturning without filtering.\n")
        return(invisible(TRUE))
      }
      
      cat("\nEnter a percentile cutoff (e.g., 0.9, or 0 to skip):\n")
      percentile_input <- as.numeric(readline(prompt = "Percentile cutoff: "))
      reproduce$percInput <- percentile_input
      
      # Check for valid numeric input
      if (is.na(percentile_input)) {
        cat("\nInvalid numeric input. Skipping filtering.\n")
        return(invisible(TRUE))
      }
      
      if (percentile_input == 0) {
        cat("\nSkipping variance filtering as requested.\n")
      } else if (percentile_input > 0 && percentile_input <= 1) {
        if (variance_choice == 1) {
          # Column Variance Filtering
          col_vars <- apply(select.ptx, 2, var)
          cutoff <- quantile(col_vars, percentile_input)
          select.ptx <<- as.data.frame(select.ptx[, col_vars > cutoff, drop = FALSE])
          cat("\nColumn variance filtering applied. Number of columns retained:", ncol(select.ptx), "\n")
        } else if (variance_choice == 2) {
          # Row Variance Filtering
          row_vars <- apply(select.ptx, 1, var)
          cutoff <- quantile(row_vars, percentile_input)
          select.ptx <<- as.data.frame(select.ptx[row_vars > cutoff, , drop = FALSE])
          # Also filter select.sinfo rows to match the filtered select.ptx rows
          if (exists("select.sinfo")) {
            select.sinfo <<- select.sinfo[rownames(select.sinfo) %in% rownames(select.ptx), ]
          }
          cat("\nRow variance filtering applied. Number of rows retained:", nrow(select.ptx), "\n")
        } else {
          cat("\nInvalid filtering option selected.\n")
        }
      } else {
        cat("\nInvalid percentile value. Please enter a value between 0 and 1.\n")
      }
      assign("reproduceVarFilter", reproduce, envir = .GlobalEnv)
    }
  return(invisible(TRUE)) 
  }
  }


################################################################################
# Normalization function:
################################################################################
normalize_data <- function(reproduce = NULL) {
  if (is.null(reproduce)) {
    
    reproduce <- data.frame(
      normChoice = NA,
      minMaxChoice = NA,
      stringsAsFactors = FALSE
    )
    cat("\n--- Normalization / Scaling ---\n")
    cat("Choose a normalization method:\n")
    cat("  1 - Inverse log₂ transformation\n")
    cat("  2 - ProtPQN\n")
    cat("  3 - ProtPQN (multiple panels)\n") # Study-specific
    cat("  4 - Scaling (z-score)\n")
    cat("  5 - Min-Max normalization\n")
    cat("  6 - VSN\n")
    cat("  7 - Inverse Normal Transformation\n")
    cat("  8 - VSN separate panels\n") # Study-specific
    
    norm_choice <- as.numeric(readline(prompt = "Enter your choice: "))
    reproduce$normChoice <- norm_choice
  } else {
    norm_choice <- reproduce$normChoice
  }
  
  if (norm_choice == 1) {
    cat("\nApplying inverse log₂ transformation to undo any prior log₂ normalization...\n")
    select.ptx <<- 2^(select.ptx)
    
  } else if (norm_choice == 2) {
    # ProtPQN normalization (over the whole data):
    if (!exists("apply_protpqn", mode = "function")) {
      cat("\nError: apply_protpqn() function is not available.\n")
      return(invisible(NULL))
    }
    cat("\nApplying ProtPQN normalization...\n")
    select.ptx <<- apply_protpqn(select.ptx, transform = TRUE)
    # Assign normalized data to global environment
    assign("select.ptx", select.ptx, .GlobalEnv)
    
  } else if (norm_choice == 3) {
    cat("\nApplying ProtPQN (multiple panels) normalization...\n")
    # Check that protpqn exists.
    if (!exists("apply_protpqn", mode = "function")) {
      cat("\nError: apply_protpqn() function is not available.\n")
      return(invisible(NULL))
    }
    # Ensure select.binfo exists and contains the panel column.
    if (!exists("select.binfo")) {
      cat("\nError: select.binfo is not available in the environment.\n")
      return(invisible(NULL))
    }
    if (!"panel" %in% colnames(select.binfo)) {
      cat("\nError: 'panel' column is not present in select.binfo.\n")
      return(invisible(NULL))
    }
    
    # Split select.ptx into two subsets based on panel info.
    cam_idx <- which(select.binfo$panel == "CAM")
    imonc_idx <- which(select.binfo$panel == "IMONC")
    
    # Check that both panels are present.
    if (length(cam_idx) == 0 || length(imonc_idx) == 0) {
      cat("\nError: One or both panel groups (CAM, IMONC) not found in select.binfo.\n")
      return(invisible(NULL))
    }
    
    # Create temporary matrices for each panel.
    select.ptx$sample_id <<- rownames(select.ptx)
    temp.CAM.ptx <- select.ptx %>% dplyr::select(sample_id, 
                                                 all_of(colnames(select.ptx)[cam_idx])) 
    temp.IMONC.ptx <- select.ptx %>% dplyr::select(sample_id, 
                                                   all_of(colnames(select.ptx)[imonc_idx]))
    
    # Apply ProtPQN normalization separately per panel
    temp.CAM.ptx <- apply_protpqn(temp.CAM.ptx, transform = TRUE, 
                                  long_format = FALSE,
                                  kitwise = FALSE)
    temp.IMONC.ptx <- apply_protpqn(temp.IMONC.ptx, transform = TRUE,
                                    long_format = FALSE,
                                    kitwise = FALSE)
    
    temp.CAM.ptx <- temp.CAM.ptx %>% dplyr::select(-sample_id)
    temp.IMONC.ptx <- temp.IMONC.ptx %>% dplyr::select(-sample_id)
    
    # Rejoin the two normalized subsets into a new select.ptx, preserving original column order.
    norm_ptx <- select.ptx %>% dplyr::select(-sample_id) 
    norm_ptx[, cam_idx] <- temp.CAM.ptx
    norm_ptx[, imonc_idx] <- temp.IMONC.ptx
    
    # Assign normalized data to global environment
    select.ptx <<- norm_ptx
    assign('select.ptx', select.ptx, .GlobalEnv)
    
  } else if (norm_choice == 4) {
    cat("\nApplying Scaling (z-score) normalization...\n")
    unscaled.ptx <<- select.ptx
    select.ptx <<- as.data.frame(scale(select.ptx, scale = TRUE, center = TRUE ))
    # Assign to global environment
    assign("select.ptx", select.ptx, .GlobalEnv)
    
  } else if (norm_choice == 5) {
    if (is.na(reproduce$minMaxChioice)) {
      cat("\nYou selected Min-Max normalization.\n")
      cat("Choose how to apply Min-Max normalization:\n")
      cat("  1 = By columns\n")
      cat("  2 = By rows\n")
      cat("  3 = Matrix wide (global normalization)\n")
      mm_choice <- as.numeric(readline(prompt = "Enter your choice: "))
      reproduce$minMaxChoice <- mm_choice
    } else {
      mm_choice <- reproduce$minMaxChoice
    }
    
    if (mm_choice == 1) {
      # Normalize by columns
      select.ptx <<- as.data.frame(apply(select.ptx, 2, min_max_normalize))
      
    } else if (mm_choice == 2) {
      # Normalize by rows
      select.ptx <<- as.data.frame(t(apply(select.ptx, 1, min_max_normalize)))
      
    } else if (mm_choice == 3) {
      # Global (rows and columns)
      mm_mat <- as.matrix(select.ptx)
      global_min <- min(mm_mat)
      global_max <- max(mm_mat)
      # Prevent division by zero
      if (global_max == global_min) {
        cat("\nAll values are the same. Skipping global normalization.\n")
      } else {
        norm_mat <- (mm_mat - global_min) / (global_max - global_min)
        select.ptx <<- as.data.frame(norm_mat)
      }
    } else {
      cat("\nInvalid choice for Min-Max normalization. Skipping normalization.\n")
    }
    
  } else if (norm_choice == 6) {
    # VSN normalization: 
    if (!exists("normalizeVSN", mode = "function")) {
      cat("\nError: normalizeVSN() function is not available.\n")
      return(invisible(NULL))
    }
    cat("\nApplying VSN normalization...\n")
    select.ptx <<- t(normalizeVSN(t(select.ptx)))
    
  } else if (norm_choice == 7) {
    # Inverse Normal Transformation using qnorm
    inverse_normalize_df(df = as.data.frame(select.ptx))
    
  } else if (norm_choice == 8) {
    cat("\nApplying vsn (multiple panels) normalization...\n")
    
    # Ensure select.binfo exists and contains the panel column.
    if (!exists("select.binfo")) {
      cat("\nError: select.binfo is not available in the environment.\n")
      return(invisible(NULL))
    }
    if (!"panel" %in% colnames(select.binfo)) {
      cat("\nError: 'panel' column is not present in select.binfo.\n")
      return(invisible(NULL))
    }
    
    # Split select.ptx (expression data) into subsets based on panel info.
    cam_idx <- which(select.binfo$panel == "CAM")
    imonc_idx <- which(select.binfo$panel == "IMONC")
    
    # Check that both panels are present.
    if (length(cam_idx) == 0 || length(imonc_idx) == 0) {
      cat("\nError: One or both panel groups (CAM, IMONC) not found in select.binfo.\n")
      return(invisible(NULL))
    }
    
    # Create temporary matrices for each panel.
    temp.CAM.ptx <- select.ptx[, cam_idx, drop = FALSE]
    temp.IMONC.ptx <- select.ptx[, imonc_idx, drop = FALSE]
    
    # Apply vsn normalization separately per panel
    temp.CAM.ptx <- t(normalizeVSN(t(temp.CAM.ptx)))
    temp.IMONC.ptx <- t(normalizeVSN(t(temp.IMONC.ptx)))
    
    # Rejoin the two normalized subsets into a new select.ptx, preserving the original column order.
    norm_ptx <- select.ptx  # copy original structure
    norm_ptx[, cam_idx] <- temp.CAM.ptx
    norm_ptx[, imonc_idx] <- temp.IMONC.ptx
    select.ptx <<- norm_ptx
    assign("select.ptx", select.ptx, .GlobalEnv)
  } else {
    cat("\nInvalid normalization option. Skipping normalization.\n")
  }
  # Assign normalized data to global environment
  assign("reproduceNormalization", reproduce, envir = .GlobalEnv)
  cat("\nNormalization completed.\n")
  return(invisible(TRUE))
}

################################################################################
# Adjustments for confounders and batch effects
################################################################################
AdjustForConfounderAndBatchEffects <- function(reproduce = NULL) {
  # Ensure required objects exist in the global environment
  if (!exists("select.sinfo", envir = .GlobalEnv))
    stop("Global object 'select.sinfo' not found.")
  if (!exists("select.ptx", envir = .GlobalEnv))
    stop("Global object 'select.ptx' not found.")
  
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      conCols = NA,
      adjMethod = NA,
      stringsAsFactors = FALSE
    )
    # List available columns in select.sinfo
    cat("Available columns in select.sinfo:\n")
    cols <- colnames(select.sinfo)
    for (i in seq_along(cols)) {
      cat(i, ": ", cols[i], "\n", sep = "")
    }
    
    # Prompt for one or more confounder column numbers (option 0 to return)
    input <- readline(prompt = "Enter the numbers corresponding to the confounder columns separated by commas (or 0 to return): ")
    if (input == "0") {
      cat("Returning to previous menu.\n")
      return(invisible(NULL))
    }
    reproduce$conCols <- input
  } else{
    input <- reproduce$conCols
  }
  
  # Split input by commas and convert to numeric indices
  conf_indices <- as.numeric(unlist(strsplit(input, ",")))
  if (any(is.na(conf_indices)) || any(conf_indices < 1) || any(conf_indices > length(cols))) {
    stop("Invalid column number selected.")
  }
  
  confounder_columns <- cols[conf_indices]
  cat("Selected confounder columns: ", paste(confounder_columns, collapse = ", "), "\n\n", sep = "")
  
  if (is.na(reproduce$adjMethod)) {
    # Prompt for choice of adjustment method, with an option to return (option 0)
    cat("Choose adjustment method:\n")
    cat("0 : Return to previous menu\n")
    cat("1 : Adjust through regression using residuals\n")
    cat("2 : Adjust using removeBatchEffect\n")
    method_choice <- as.numeric(readline(prompt = "Enter 0, 1, or 2: "))
    if (is.na(method_choice) || !(method_choice %in% c(0, 1, 2))) {
      stop("Invalid method choice.")
    }
    if (method_choice == 0) {
      cat("Returning to previous menu.\n")
      return(invisible(NULL))
    }
    reproduce$adjMethod <- method_choice
  } else {
    method_choice <- reproduce$adjMethod
  }
  
  if (method_choice == 1) {
    cat("Running adjustment through regression using residuals...\n")
    
    # Build a design matrix using the selected confounder columns.
    # This constructs a formula like: ~ conf1 + conf2 etc...
    formula_str <- paste("~", paste(confounder_columns, collapse = " + "))
    design <- model.matrix(as.formula(formula_str), data = select.sinfo)
    fit <- limma::lmFit(t(select.ptx), design)
    
    # Calculate adjusted expression values by subtracting the confounder effects.
    adjustedExpr <- select.ptx - design %*% t(fit$coefficients)
    
    # Update the global select.ptx with the adjusted expression data.
    unadjusted.ptx <<- select.ptx
    select.ptx <<- adjustedExpr
    cat("Regression adjustment complete. Global variable 'select.ptx' updated.\n")
    
  } else if (method_choice == 2) {
    cat("Running adjustment using removeBatchEffect...\n")
    
    # For removeBatchEffect, convert the selected confounders into a numeric covariate matrix.
    confounder_data <- select.sinfo[, confounder_columns, drop = FALSE]
    covariate_design <- model.matrix(~ . , data = confounder_data)[, -1, drop = FALSE]
    adjustedExpr <- limma::removeBatchEffect(t(select.ptx), covariates = covariate_design)
    
    # Update the global select.ptx after transposing back.
    unadjusted.ptx <<- select.ptx
    select.ptx <<- t(adjustedExpr)
    assign("reproduceAdjustment", reproduce, envir = .GlobalEnv)
    cat("removeBatchEffect adjustment complete. Global variable 'select.ptx' updated.\n")
  }
  
  cat("Adjustment complete.\n")
}


################################################################################
# Function for inspection of loadings (PCA, PLS-DA, sPLS-DA)
################################################################################
loading_inspector <- function() {
  cat("\n--- Loadings Inspection ---\n")
  
  # Choose method:
  cat("Select the method for loadings inspection:\n")
  cat("  1 = PCA\n")
  cat("  2 = PLS-DA\n")
  cat("  3 = sPLS-DA\n")
  method_choice <- as.numeric(readline(prompt = "Enter your choice: "))
  
  if (method_choice == 1) {
    if (!exists("pca.loadings")) {
      cat("\nError: 'pca.loadings' not found in the environment.\n")
      return(invisible(NULL))
    }
    loadings_obj <- get("pca.loadings", envir = .GlobalEnv)
    method_label <- "PCA"
  } else if (method_choice == 2) {
    if (!exists("plsda.loadings")) {
      cat("\nError: 'plsda.loadings' not found in the environment.\n")
      return(invisible(NULL))
    }
    loadings_obj <- get("plsda.loadings", envir = .GlobalEnv)
    method_label <- "PLS-DA"
  } else if (method_choice == 3) {
    if (!exists("splsda.loadings")) {
      cat("\nError: 'splsda.loadings' not found in the environment.\n")
      return(invisible(NULL))
    }
    loadings_obj <- get("splsda.loadings", envir = .GlobalEnv)
    method_label <- "sPLS-DA"
  } else {
    cat("\nInvalid selection.\n")
    return(invisible(NULL))
  }
  
  # Ask which component to inspect (default column 1)
  comp_input <- as.numeric(readline(prompt = "Enter the component (column) number to inspect [default = 1]: "))
  if(is.na(comp_input) || comp_input < 1 || comp_input > ncol(loadings_obj)) {
    comp_num <- 1
    cat("No valid component provided. Defaulting to component 1.\n")
  } else {
    comp_num <- comp_input
  }
  
  # Ask which influence type
  cat("\nSelect influence option:\n")
  cat("  1 = Absolute influence (order proteins by absolute loading value)\n")
  cat("  2 = Positive/Negative influence (display top proteins with positive and negative loadings)\n")
  influence_choice <- as.numeric(readline(prompt = "Enter your choice: "))
  if(is.na(influence_choice) || !(influence_choice %in% c(1, 2))) {
    cat("Invalid input. Defaulting to Absolute influence.\n")
    influence_choice <- 1
  }
  
  # Ask how many proteins to display
  num_proteins <- as.numeric(readline(prompt = "Enter the number of proteins to show: "))
  if(is.na(num_proteins) || num_proteins < 1) {
    cat("Invalid number. Defaulting to 10 proteins.\n")
    num_proteins <- 10
  }
  
  # Extract the loadings for the chosen component
  loadings_vec <- loadings_obj[, comp_num]
  
  if(influence_choice == 1) {
    # Absolute influence option:
    # Order proteins by the absolute value of their loadings
    top_proteins <- sort(abs(loadings_vec), decreasing = TRUE)
    n_to_show <- min(num_proteins, length(top_proteins))
    top_proteins <- head(top_proteins, n_to_show)
    
    cat("\nTop", n_to_show, "proteins by absolute loading (", method_label, " component ", comp_num, "):\n", sep = " ")
    print(top_proteins)
    
    # Retrieve the original loading values using the protein names
    top_names <- names(top_proteins)
    top_values <- loadings_vec[top_names]
    
    barplot(top_values,
            main = paste("Top", n_to_show, "Protein Loadings for", method_label, "Comp", comp_num),
            ylab = "Loading Value", las = 2, col = "skyblue")
    
  } else if (influence_choice == 2) {
    # Positive/Negative influence option:
    n_to_show <- num_proteins  # Show this many proteins for each sign.
    
    # For top positive influences, sort in descending order.
    top_positive <- sort(loadings_vec, decreasing = TRUE)[1:n_to_show]
    
    # For top negative influences, sort in ascending order.
    top_negative <- sort(loadings_vec, decreasing = FALSE)[1:n_to_show]
    
    cat("\nTop", n_to_show, "proteins with the most positive loadings for", method_label, "Comp", comp_num, ":\n")
    print(top_positive)
    
    cat("\nTop", n_to_show, "proteins with the most negative loadings for", method_label, "Comp", comp_num, ":\n")
    print(top_negative)
    
    # Set up side-by-side plotting panels
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(1, 2))
    
    # Plot barplots
    barplot(top_positive,
            main = paste("Positive Loadings\n", method_label, "Comp", comp_num),
            ylab = "Loading Value", las = 2, col = "darkgreen")
    
    barplot(top_negative,
            main = paste("Negative Loadings\n", method_label, "Comp", comp_num),
            ylab = "Loading Value", las = 2, col = "red")
    
    par(old_par)
  }
  
  invisible(list(method = method_label,
                 component = comp_num,
                 influence = if(influence_choice == 1) "Absolute" else "Positive/Negative"))
}
