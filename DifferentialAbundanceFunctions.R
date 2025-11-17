# ------------------------------------------------------------------------------
# Script:        DifferentialAbundanceFunctions.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Functions for Differential Abundance Analysis using limma
# ------------------------------------------------------------------------------

################################################################################
# Functions for Differential Abundance Analysis
################################################################################

################################################################################
# Interactive Design Matrix builder
################################################################################
build_design_matrix_interactive <- function(df, preselect = NULL, save_name = "last_design_select") {
  
  # Load preselection for reproduction or run interactive 
  if (!is.null(preselect)) {
    group_traits <- preselect$group_traits
    contrast_traits <- preselect$contrast_traits
  } else {
    group_traits <- list()
    contrast_traits <- list()
    
    repeat {
      cat("\nCurrent selections:\n")
      cat("Group traits:\n")
      if (length(group_traits) == 0) {
        cat("  (none)\n")
      } else {
        for (i in seq_along(group_traits)) {
          cat(sprintf("  [%d] %s = %s\n", i, group_traits[[i]]$trait, group_traits[[i]]$level))
        }
      }
      cat("Contrast traits:\n")
      if (length(contrast_traits) == 0) {
        cat("  (none)\n")
      } else {
        for (i in seq_along(contrast_traits)) {
          cat(sprintf("  [%d] %s = %s\n", i, contrast_traits[[i]]$trait, contrast_traits[[i]]$level))
        }
      }
      
      # Menu for composing the design matrix
      cat("\nOptions:\n")
      cat("1. Add group trait\n")
      cat("2. Add contrasting trait\n")
      cat("3. Remove group trait\n")
      cat("4. Remove contrasting trait\n")
      cat("5. Finish design\n")
      
      choice <- as.integer(readline(prompt = "Enter choice: "))
      
      if (choice == 1) {
        sel <- choose_filter(df)
        group_traits <- append(group_traits, list(list(trait = sel$filter_trait, level = sel$filter_level)))
        
      } else if (choice == 2) {
        sel <- choose_filter(df)
        contrast_traits <- append(contrast_traits, list(list(trait = sel$filter_trait, level = sel$filter_level)))
        
      } else if (choice == 3) {
        if (length(group_traits) == 0) {
          cat("No group traits to remove.\n")
        } else {
          rem <- as.integer(readline(prompt = "Enter number of group trait to remove: "))
          if (!is.na(rem) && rem >= 1 && rem <= length(group_traits)) {
            group_traits <- group_traits[-rem]
          } else {
            cat("Invalid selection.\n")
          }
        }
        
      } else if (choice == 4) {
        if (length(contrast_traits) == 0) {
          cat("No contrast traits to remove.\n")
        } else {
          rem <- as.integer(readline(prompt = "Enter number of contrast trait to remove: "))
          if (!is.na(rem) && rem >= 1 && rem <= length(contrast_traits)) {
            contrast_traits <- contrast_traits[-rem]
          } else {
            cat("Invalid selection.\n")
          }
        }
        
      } else if (choice == 5) {
        break
      } else {
        cat("Invalid choice.\n")
      }
    }
  }
  
  # Save current selections for reproduction
  preselect <- list(
    group_traits = group_traits,
    contrast_traits = contrast_traits
  )
  assign(save_name, preselect, envir = .GlobalEnv)
  
  # Detect if this is a cluster vs cluster case
  same_trait <- length(group_traits) == 1 &&
    length(contrast_traits) == 1 &&
    group_traits[[1]]$trait == contrast_traits[[1]]$trait
  
  if (same_trait) {
    cat("\nDetected same trait for group and contrast — using full model mode.\n")
    
    factor_var <- factor(df[[group_traits[[1]]$trait]])
    
    # Create valid R names for levels
    safe_levels <- levels(factor_var)
    safe_levels <- gsub(">", "gt", safe_levels)
    safe_levels <- gsub("<", "lt", safe_levels)
    safe_levels <- gsub("=", "", safe_levels)   # remove equals if present
    safe_levels <- gsub("[^[:alnum:]_]", "_", safe_levels) # replace other non-alphanumerics with _
    
    # Apply safe names to factor levels
    levels(factor_var) <- safe_levels
    
    # Build design matrix with safe column names
    design <- model.matrix(~0 + factor_var)
    colnames(design) <- safe_levels
    
    # Sanitize the group/contrast level names for the contrast string
    safe_group <- gsub(">", "gt", group_traits[[1]]$level)
    safe_group <- gsub("<", "lt", safe_group)
    safe_group <- gsub("=", "", safe_group)
    safe_group <- gsub("[^[:alnum:]_]", "_", safe_group)
    
    safe_contrast <- gsub(">", "gt", contrast_traits[[1]]$level)
    safe_contrast <- gsub("<", "lt", safe_contrast)
    safe_contrast <- gsub("=", "", safe_contrast)
    safe_contrast <- gsub("[^[:alnum:]_]", "_", safe_contrast)
    
    contrast_name <- paste(safe_group, "-", safe_contrast)
    
    contrast.matrix <- makeContrasts(contrasts = contrast_name, levels = design)
    
    return(list(
      design = design,
      contrast.matrix = contrast.matrix,
      samples = rownames(df),
      mode = "full_model",
      group_traits = group_traits,
      contrast_traits = contrast_traits
    ))
  }
  else {
    # Filtered two-group mode:
    group_flag <- rep(TRUE, nrow(df))
    if (length(group_traits) > 0) {
      for (gt in group_traits) {
        group_flag <- group_flag & (as.character(df[[gt$trait]]) == gt$level)
      }
    } else {
      group_flag <- rep(FALSE, nrow(df))
    }
    
    if (length(contrast_traits) > 0) {
      contrast_flag <- rep(TRUE, nrow(df))
      for (ct in contrast_traits) {
        contrast_flag <- contrast_flag & (as.character(df[[ct$trait]]) == ct$level)
      }
    } else {
      contrast_flag <- !group_flag
    }
    
    group_factor <- rep(NA, nrow(df))
    group_factor[group_flag] <- "Group"
    group_factor[contrast_flag] <- "Contrast"
    group_factor <- factor(group_factor, levels = c("Group", "Contrast"))
    
    keep <- !is.na(group_factor)
    df <- df[keep, , drop = FALSE]
    group_factor <- droplevels(group_factor[keep])
    
    design <- model.matrix(~0 + group_factor)
    colnames(design) <- levels(group_factor)
    contrast.matrix <- makeContrasts("Group - Contrast", levels = design)
    
    return(list(
      design = design,
      contrast.matrix = contrast.matrix,
      samples = rownames(df),
      mode = "filtered_two_group",
      group_traits = group_traits,
      contrast_traits = contrast_traits
    ))
  }
}

################################################################################
# Differential abundance analysis runner for categorical variables
################################################################################
run_differential_abundance <- function(design_output, 
                                       p_cutoff = 0.05, 
                                       logFC_cutoff = 0.3,
                                       adj_method = "bonferroni") {
  cat("### Differential Abundance Analysis ###\n\n")
  
  # Prompt for which data to use
  cat("Select expression data type to use for analysis:\n")
  cat(" 1 - Scaled data (select.ptx)\n")
  cat(" 2 - Unscaled data (unscaled.ptx)\n")
  cat(" 3 - Unadjusted data (unadjusted.ptx)\n")
  cat(" 4 - Unreduced data (unreduced.ptx)\n")
  data_choice <- as.integer(readline(prompt = "Enter your choice (1, 2, 3, or 4): "))
  
  if (!is.na(data_choice) && data_choice == 2 && exists("unscaled.ptx")) {
    expr_matrix <- unscaled.ptx
    cat("Using unscaled data (unscaled.ptx).\n\n")
  } else if (!is.na(data_choice) && data_choice == 3 && exists("unadjusted.ptx")) {
    expr_matrix <- unadjusted.ptx
    cat("Using unadjusted data (unadjusted.ptx).\n\n")
  } else if (!is.na(data_choice) && data_choice == 4 && exists("unreduced.ptx")) {
    expr_matrix <- unreduced.ptx
    cat("Using unreduced data (unreduced.ptx).\n\n")
  } else {
    expr_matrix <- select.ptx
    cat("Using scaled data (select.ptx).\n\n")
  }
  
  # Extract design and contrast created using the design function
  design <- design_output$design
  contrast.matrix <- design_output$contrast.matrix
  samples <- design_output$samples
  
  # Ensure sample order matches
  if (!all(rownames(expr_matrix) %in% samples)) {
    stop("Sample names in expression matrix do not match design matrix samples.")
  }
  expr_matrix <- expr_matrix[samples, ]
  expr_matrix <- t(expr_matrix)
  
  # Run limma
  fit <- limma::lmFit(expr_matrix, design)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit3 <<- limma::eBayes(fit2)   # <-- keep this as the post-eBayes object
  
  # Extract results for each contrast
  results_list <- list()
  for (contrast in colnames(contrast.matrix)) {
    results_list[[contrast]] <- limma::topTable(
      fit3, 
      coef = contrast, 
      number = Inf, 
      adjust.method = "bonferroni"
    )
  }
  
  # Generate Volcano plots
  VP_plots <- list()
  for (contrast in names(results_list)) {
    VP_results <- results_list[[contrast]]
    plt <- EnhancedVolcano::EnhancedVolcano(
      VP_results,
      lab = rownames(VP_results),
      x = "logFC",
      y = "adj.P.Val",
      title = contrast,
      pCutoff = p_cutoff,
      FCcutoff = logFC_cutoff,
      pointSize = 3,
      labSize = 5,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      colConnectors = "grey30"
    )
    VP_plots[[contrast]] <- plt
    
    # Plot the volcano plot
    print(plt)
    
    # Prompt if to save as variable and for variable name
    save_choice <- readline(prompt = "Do you want to save the plot? y/n: ")
    if (tolower(save_choice) == "y") {
      var_name2 <- readline("Enter a name to save your plot object under: ")
      assign(var_name2, plt, .GlobalEnv)
      cat("Plot saved as variable:", var_name2, "\n")
    } else {
      cat("Plot not saved as a variable\n")
    }
  }
  cat("Results saved in variable 'fit3'\n")
}

################################################################################
# Differential abundance analysis function for continuous variables
################################################################################
DEAcontinuous <- function() {
  # Check required objects in the GlobalEnv
  if (!exists("select.sinfo", envir = .GlobalEnv))
    stop("Global object 'select.sinfo' not found.")
  if (!exists("unreduced.ptx", envir = .GlobalEnv))
    stop("Global object 'unreduced.ptx' not found.")
  
  # Main menu loop
  main_loop <- TRUE
  while(main_loop) {
    cat("\n--- Clinical Variable Selection ---\n")
    clinical_cols <- colnames(select.sinfo)
    for (i in seq_along(clinical_cols)) {
      cat(i, ": ", clinical_cols[i], "\n", sep = "")
    }
    # Prompt for clinical covariate for design
    design_input <- trimws(readline(
      prompt = "Enter the numbers corresponding to the clinical covariates (comma-separated), or 'q' to quit: "
    ))
    if(tolower(design_input) == "q") {
      cat("Exiting analysis.\n")
      return(invisible(NULL))
    }
    design_indices <- as.numeric(unlist(strsplit(design_input, ",")))
    if (any(is.na(design_indices)) || any(design_indices < 1) || any(design_indices > length(clinical_cols))) {
      stop("Invalid clinical covariate selection.")
    }
    design_vars <- clinical_cols[design_indices]
    cat("Selected covariates for the design matrix: ", paste(design_vars, collapse = ", "), "\n")
    
    # Ensure each selected variable is an unordered factor
    for (v in design_vars) {
      if (is.ordered(select.sinfo[[v]])) {
        cat("Converting ordered factor '", v, "' to an unordered factor.\n", sep = "")
        select.sinfo[[v]] <- factor(select.sinfo[[v]], ordered = FALSE)
      }
    }
    # Prompt for variable of interest to use as contrast 
    cat("\nSelect the variable to contrast (from the selected design variables):\n")
    for (i in seq_along(design_vars)) {
      cat(i, ": ", design_vars[i], "\n", sep = "")
    }
    contrast_input <- trimws(readline(prompt = "Enter the number corresponding to the variable of interest: "))
    contrast_index <- as.numeric(contrast_input)
    if (is.na(contrast_index) || contrast_index < 1 || contrast_index > length(design_vars)) {
      stop("Invalid variable of interest selection.")
    }
    contrast_var <- design_vars[contrast_index]
    cat("Selected variable of interest (contrast): ", contrast_var, "\n")
    
    # Handle contrast creation differently based on variable type
    if (is.factor(select.sinfo[[contrast_var]]) || is.character(select.sinfo[[contrast_var]])) {
      # Convert to factor if it is not already
      select.sinfo[[contrast_var]] <- factor(select.sinfo[[contrast_var]], ordered = FALSE)
      contrast_levels <- levels(select.sinfo[[contrast_var]])
      cat("\nLevels for ", contrast_var, ":\n", sep = "")
      for (i in seq_along(contrast_levels)) {
        cat(i, ": ", contrast_levels[i], "\n", sep = "")
      }
      
      # Interpret the two numbers as (baseline, comparison)
      # Prompt for which levels to contrast
      level_input <- trimws(readline(
        prompt = "Enter the two numbers (comma-separated) corresponding to the levels to contrast (baseline,comparison): "
      ))
      level_indices <- as.numeric(unlist(strsplit(level_input, ",")))
      if (length(level_indices) != 2 ||
          any(is.na(level_indices)) ||
          any(level_indices < 1) ||
          any(level_indices > length(contrast_levels))) {
        stop("Invalid level selection for contrast.")
      }
      # Use the first number as the baseline (reference) and the second one as the comparison
      baseline_level <- contrast_levels[level_indices[1]]
      comparison_level <- contrast_levels[level_indices[2]]
      cat("Contrast: ", comparison_level, " vs ", baseline_level, "\n")
      
      # Relevel the factor so that the chosen baseline becomes the reference
      select.sinfo[[contrast_var]] <- relevel(select.sinfo[[contrast_var]], ref = baseline_level)
      
      # Rebuild the design matrix with updated factor coding
      design_formula <- as.formula(paste("~", paste(design_vars, collapse = " + ")))
      design_matrix <- model.matrix(design_formula, data = select.sinfo)
      colnames(design_matrix) <- make.names(colnames(design_matrix))
      
      # Construct the candidate contrast column using the comparison level
      candidate_contrast <- make.names(paste(contrast_var, comparison_level, sep = ""))
      if (!(candidate_contrast %in% colnames(design_matrix))) {
        stop("Candidate contrast column '", candidate_contrast, "' not found in the design matrix. Check factor levels and coding.")
      }
      contrast_matrix <- limma::makeContrasts(contrasts = candidate_contrast, levels = design_matrix)
    } else {
      # For numeric (continous) variables
      design_formula <- as.formula(paste("~", paste(design_vars, collapse = " + ")))
      design_matrix <- model.matrix(design_formula, data = select.sinfo)
      colnames(design_matrix) <- make.names(colnames(design_matrix))
      cat("The selected variable is numeric. Its coefficient will be used for differential expression.\n")
      
      contrast_matrix <- diag(ncol(design_matrix))
      colnames(contrast_matrix) <- colnames(design_matrix)
      contrast_matrix <- contrast_matrix[, contrast_var, drop = FALSE]
    }
    
    # Prepare the protein data for limma
    expression_data <- t(unreduced.ptx)
    
    # Run limma
    fit <- limma::lmFit(expression_data, design_matrix)
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    fit <- limma::eBayes(fit)
    results <<- limma::topTable(fit, number = Inf, adjust.method = "bonferroni", sort.by = "P")
    
    cat("Differential expression analysis complete.\n")
    
    # Inner loop: choose plot type or re-run with different variables
    repeat {
      cat("\n--- Plotting Options ---\n")
      plot_choice <- trimws(readline(
        prompt = "Enter '1' for volcano plot, '2' for heatmap, '3' to choose another clinical variable(s), or '0' to quit: "
      ))
      if (tolower(plot_choice) == "1") {
        # Volcano plot
        volcano_plot <- EnhancedVolcano::EnhancedVolcano(
          results,
          lab = ifelse(results$adj.P.Val < 0.05, rownames(results), ""),
          x = "logFC",
          y = "adj.P.Val",
          pCutoff = 0.05,
          FCcutoff = 0.25,
          title = "Differential Expression Analysis",
          selectLab = rownames(results)[-log10(results$adj.P.Val) > 2.3],
          subtitle = paste("Contrast:", contrast_var),
          xlab = bquote(~Log[2]~ "Fold Change"),
          ylab = bquote(~-Log[10]~ "Adjusted P-value"),
          drawConnectors = TRUE,
          widthConnectors = 0.5,
          colConnectors = "grey30",
          labSize = 3.5
        )
        # Prompt if to save plot and under what variable name
        save_choice <- readline(prompt ="Do you want to save the plot? y/n: \n")
        if (save_choice == "y") {
          var_name2 <- readline("Enter a name to save your plot object under: ")
          assign(var_name2, volcano_plot, .GlobalEnv ) 
          print(volcano_plot) 
        } else {
          cat("Plot not saved as a variable")
          print(volcano_plot)
        }
      } else if (tolower(plot_choice) == "2") {
        # Heatmap:
        signif_proteins <- rownames(results)[results$adj.P.Val < 0.05]
        if (length(signif_proteins) < 5) {
          cat("Fewer than 5 proteins reached significance (adj.P.Val < 0.05). Displaying top 50 proteins by p-value.\n")
          top50 <- head(rownames(results[order(results$adj.P.Val), ]), 50)
          signif_proteins <- top50
        }
        common_proteins <- intersect(signif_proteins, rownames(expression_data))
        if(length(common_proteins) == 0) {
          stop("No matching proteins found between the results and the expression data. Check row names.")
        }
        heatmap_data <- expression_data[common_proteins, , drop = FALSE]
        heatmap_data_scaled <- t(scale(t(heatmap_data)))
        heatmap_plot <- ComplexHeatmap::Heatmap(
          heatmap_data_scaled,
          name = "Z-score Expression",
          column_title = "Samples",
          row_title = "Proteins",
          show_row_names = TRUE,
          show_column_names = TRUE,
          clustering_distance_rows = "euclidean",
          clustering_distance_columns = "euclidean",
          clustering_method_rows = "complete",
          clustering_method_columns = "complete"
        )
        draw(heatmap_plot, heatmap_legend_side = "right")
      } else if (tolower(plot_choice) == "3") {
        break  # Back to clinical variable selection.
      } else if (tolower(plot_choice) == "0") {
        cat("Exiting analysis.\n")
        return(invisible(results))
      } else {
        cat("Invalid input. Please try again.\n")
      }
    }
  }
}

