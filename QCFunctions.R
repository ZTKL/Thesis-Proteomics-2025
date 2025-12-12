# ------------------------------------------------------------------------------
# Script:        QCFunctions.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Functions for the Quality Control Workflow of the Suite
# ------------------------------------------------------------------------------

################################################################################
# QC Workflow main menu function:
################################################################################
QC_main_menu <- function() {
  # Run preliminary variable setup
  vars <- setup_variables()
  if (is.null(vars)) return(invisible(NULL))
  ptx   <- vars$ptx
  sinfo <- vars$sinfo
  binfo <- vars$binfo
  
  repeat {
    cat("\n=== QC Workflow Menu: ===\n")
    cat("  1 = Filter proteins based on frequency below LOD\n")
    cat("  2 = Remove technical duplicates\n")
    cat("  3 = Remove samples flagged with a WARNING\n")
    cat("  4 = Handle NAs\n")
    cat("  5 = Identify sample outliers with OlinkAnalyze\n")
    cat("  6 = Identify sample outliers with PCA\n")
    cat("  7 = Identify protein outliers with PCA\n")
    cat("  8 = Hierarchical Clustering for sample outlier identification\n")
    cat("  9 = Protein and Sample Correlation Heatmaps\n")
    cat("  0 = Exit QC Workflow\n")
    qc_choice <- read_choice("Enter your choice [0-9]: ")
    
    if (is.na(qc_choice) || qc_choice == 0) {
      cat("\nExiting QC workflow.\n")
      break
    } else if (qc_choice == 1) {
      filterfbLOD(ptx, sinfo, binfo, reproduce = NULL)
    } else if (qc_choice == 2) {
      remove_technical_duplicates(reproduce = NULL)
    } else if (qc_choice == 3) {
      remove_warning_samples(select.ptx, select.sinfo)
    } else if (qc_choice == 4) {
      res <- handle_NAs(ptx = select.ptx, sinfo = select.sinfo)
      select.ptx <<- res$ptx
      select.sinfo <<- res$sinfo
    } else if (qc_choice == 5) {
      olink_runner(select.ptx = select.ptx, select.sinfo = select.sinfo
                   , select.binfo = select.binfo)
    } else if (qc_choice == 6) {
      run_pca_filtering_normal(ptx = select.ptx, sinfo = select.sinfo, reproduce = NULL)
    } else if (qc_choice == 7) {
      run_pca_filtering_proteins(ptx = select.ptx, binfo = select.binfo,
                                   sinfo = select.sinfo, reproduce = NULL)
    } else if (qc_choice == 8) {
      res <- run_hierarchical_clustering(ptx, sinfo)
      ptx   <- res$ptx; sinfo <- res$sinfo;
    } else if (qc_choice == 9) {
      analyze_correlation(reproduce = NULL)
    } else {
      cat("\nInvalid QC option. Please select again.\n")
    }
  }
  cat("\nQC analysis complete.\n")
  return(invisible(TRUE))
}

################################################################################
# Outlier Identification Functions
################################################################################

# Hierarchical Clustering for outlier identification 
########################################################################
run_hierarchical_clustering_outliers <- function(ptx, sinfo) {
  # Interactive trait selection
  Traits4color <- choose_traits(sinfo)
  
  
  temp.sinfo <- sinfo
  temp.sinfo[Traits4color] <- lapply(temp.sinfo[Traits4color],
                                     function(x) as.numeric(as.factor(x)))
  cat("Running Hierarchical Clustering...\n")
  htr <- hclust(dist(ptx), method = "average")
  tr.colors <- WGCNA::numbers2colors(temp.sinfo[, Traits4color, drop = FALSE], 
                                     signed = FALSE)
  WGCNA::plotDendroAndColors(htr, tr.colors,
                                        groupLabels = names(temp.sinfo[, Traits4color
                                                                       , drop = FALSE]),
                                        main = "Sample dendrogram and Trait Heatmap")
  # Interactive loop for sample removal:
  repeat {
    cat("\nWhich sample to remove? Enter sample row name or 0 to exit:\n")
    sample_to_remove <- trimws(readline())
    
    if (sample_to_remove == "0") {
      cat("\nReturning to QC Menu.\n")
      break
      
    } else {
      if (sample_to_remove %in% rownames(ptx)) {
        # Remove samples from local ptx and sinfo
        ptx <- ptx[!rownames(ptx) %in% sample_to_remove, , drop = FALSE]
        sinfo <- sinfo[!rownames(sinfo) %in% rownames(ptx), , drop = FALSE]
        
        # Update global variable
        assign("select.ptx", ptx, envir = .GlobalEnv)
        assign("select.sinfo", sinfo, envir = .GlobalEnv)
        cat("\nSample", sample_to_remove, "removed successfully.\n")
      } else{
        cat("\nInvalid sample selection. Please try again.\n")
      }
    }
  }
  return(list(ptx = ptx, sinfo = sinfo))
}

################################################################################
# Removal of warnings and duplicates
################################################################################

################################################################################ 
# Remove samples flagged with WARNING (only CAM and IMONC):
################################################################################
remove_warning_samples <- function(ptx, sinfo) {
  # Check if required variables exist
  if (!exists("select.sinfo") || !exists("select.ptx")) {
    cat("ERROR: Either 'select.sinfo' or 'select.ptx' was not found in workspace\n")
    return(invisible(NULL))
  }
  
  # Check for required warning columns
  required_cols <- c("qc_warn_CAM", "qc_warn_IMONC")
  if (!all(required_cols %in% colnames(select.sinfo))) {
    cat("ERROR: One or both of columns 'qc_warn_CAM' and 'qc_warn_IMONC' are missing from select.sinfo\n ")
    return(invisible(NULL))
  }
  
  warning_samples <- rownames(select.sinfo)[
    select.sinfo$qc_warn_CAM == "Warning" | select.sinfo$qc_warn_IMONC == "Warning"
  ]
  if (length(warning_samples) > 0) {
    cat("Removing the following sample(s) due to Warning(s):\n", 
        paste(warning_samples, collapse = ", "), "\n")
  } else { 
    cat("No samples with 'Warning' found. No changes made.\n")
    return(list(sinfo = select.sinfo, ptx = select.ptx))
  }
  # Remove identified Warning samples
  new_sinfo <- select.sinfo[!rownames(select.sinfo) %in% warning_samples, , drop = FALSE]
  new_ptx <- select.ptx[!rownames(select.ptx) %in% warning_samples, , drop = FALSE]
  
  cat("After removal, select.sinfo has", nrow(new_sinfo), "samples and select.ptx has",
      nrow(new_ptx), "samples\n")
  
  # Update global variables
  assign("select.sinfo", new_sinfo, envir = .GlobalEnv)
  assign("select.ptx", new_ptx, envir = .GlobalEnv)
  
  return(list(sinfo = new_sinfo, ptx = new_ptx))
}

################################################################################ 
# Removal of duplicate samples
################################################################################ 
remove_technical_duplicates <- function(sinfo, ptx, reproduce = NULL) {
  # Check that required dataframes exist in global environment
  if (!exists("select.sinfo", envir = .GlobalEnv))
    stop("Global object 'select.sinfo' not found")
  if (!exists("select.ptx", envir = .GlobalEnv))
    stop("Global object 'select.ptx' not found")
    
  duplicates <- select.sinfo %>%
    group_by(individual_id) %>%
    filter(n() == 2) %>%
    ungroup()
  
  duplicates <- as.data.frame(duplicates)
  # Set rownames to sample_id 
  rownames(duplicates) <- duplicates$sample_id
  
  cat("Found", nrow(duplicates), "duplicate samples (in pairs) based on individual_id.\n")
  
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      option = NA,
      pbdrw = NA,
      fblChoice = NA,
      stringsAsFactors = FALSE
    )
    # Prompt for removal option
    cat("Choose duplicate removal option:\n")
    cat("1 = Remove based on primary_blood_draw.\n")
    cat("2 = Remove based on freq_below_lod.\n")
    option <- trimws(readline(prompt = "Enter 1 or 2: "))
    
    if (option == "1") {
      # Removal based on primary_blood_draw
      cat("Option 1: Remove based on primary_blood_draw.\n")
      cat("1 = Remove samples from primary blood draw.\n")
      cat("2 = Remove samples not from primary blood draw.\n")
      choice <- trimws(readline(prompt = "Enter your choice: \n"))
      if (choice == "1") {
        pbdrw <- 1
      }
      else if (choice == "2") {
        pbdrw <- 0
      } else {
        stop("Invalid input.")
      }
      
      # Keep specified duplicates
      duplicates <- duplicates[duplicates$primary_blood_draw == pbdrw, ]
      cat("Removing", nrow(duplicates), "samples based on primary_blood_draw.\n")
      select.sinfo <<- select.sinfo[!rownames(select.sinfo) %in% rownames(duplicates), ]
      select.ptx <<- select.ptx[rownames(select.ptx) %in% rownames(select.sinfo), ]
      reproduce$option <- "primary_blood_draw"
      reproduce$pbdrw <- pbdrw
      
    } else if (option == "2") {
      # Option 2: removal based on freq_below_lod
      cat("Option 2: Remove based on freq_below_lod.\n")
      cat("1 - Remove samples with the HIGHEST freq_below_lod.\n")
      cat("2 - Remove samples with the LOWEST freq_below_lod.\n")
      removal_choice_input <- trimws(readline(prompt = "Enter your choice: "))
      
      if (removal_choice_input == "1") {
        removal_choice <- 1
      } else if (removal_choice_input == "2") {
        removal_choice <- 2
      } else {
        stop("Invalid Input.")
      }
      ids_to_remove <- c()
      unique_ids <- unique(duplicates$individual_id)
      for (id in unique_ids) {
        temp <- duplicates[duplicates$individual_id == id, ]
        # Make sure there are 2 samples in each group
        if (nrow(temp) != 2) next
        if(removal_choice == 1) {
          removal_sample <- rownames(temp)[which.max(temp$freq_below_lod)]
        } else {
          removal_sample <- rownames(temp)[which.min(temp$freq_below_lod)]
        }
        ids_to_remove <- c(ids_to_remove, removal_sample)
      }
      duplicates_to_remove <- duplicates[rownames(duplicates) %in% ids_to_remove, ]
      cat("Removing", nrow(duplicates_to_remove), "samples based on freq_below_lod.\n")
      select.sinfo <- select.sinfo[rownames(select.sinfo) %in% rownames(duplicates_to_remove), ]
      select.ptx <- select.ptx[rownames(select.ptx) %in% rownames(select.sinfo), ]
      reproduce$option <- "freq_below_lod"
      reproduce$fblChoice <- removal_choice
    } else { 
      cat("Invalid option. Exiting the duplicate removal function.\n")
      return(invisible(NULL))
    }
    assign("reproduceDuplicateRemoval", reproduce, envir = .GlobalEnv)
    
  } else {
    if (reproduce$option == "primary_blood_draw") {
      pbdrw <- reproduce$pbdrw
      duplicates <- duplicates[duplicates$primary_blood_draw == pbdrw, ]
      select.sinfo <<- select.sinfo[!rownames(select.sinfo) %in% rownames(duplicates), ]
      select.ptx <<- select.ptx[rownames(select.ptx) %in% rownames(select.sinfo), ]
    } else if (reproduce$option == "freq_below_lod") {
      removal_choice <- reproduce$fblChoice
      ids_to_remove <- c()
      unique_ids <- unique(duplicates$individual_id)
      for (id in unique_ids) {
        temp <- duplicates[duplicates$individual_id == id, ]
        # Make sure there are 2 samples in each group
        if (nrow(temp) != 2) next
        if(removal_choice == 1) {
          removal_sample <- rownames(temp)[which.max(temp$freq_below_lod)]
        } else if (removal_choice == 2) {
          removal_sample <- rownames(temp)[which.min(temp$freq_below_lod)]
        }
        ids_to_remove <- c(ids_to_remove, removal_sample)
      }
      duplicates_to_remove <- duplicates[rownames(duplicates) %in% ids_to_remove, ]
      select.sinfo <- select.sinfo[!rownames(select.sinfo) %in% rownames(duplicates_to_remove), ]
      select.ptx <- select.ptx[rownames(select.ptx) %in% rownames(select.sinfo), ]
      
    }
  }
  assign("select.ptx", select.ptx, envir = .GlobalEnv)
  assign("select.sinfo", select.sinfo, envir = .GlobalEnv)
}

################################################################################  
# PCA filtering of protein outliers
################################################################################ 
run_pca_filtering_proteins <- function(ptx, sinfo, binfo, reproduce = NULL) {
  # Transpose for PCA
  ptx_trans <- t(ptx)
  pca <- prcomp(as.data.frame(ptx_trans), scale. = TRUE)
  pcaX <- as.data.frame(pca$x, row.names = rownames(ptx_trans))
  pca.var <- pca$sdev^2
  pca.var.percent <- round((pca.var / sum(pca.var)) * 100, 2)
  
  if (is.null(reproduce)) {
    # Interactive branch
    sd_factor <- as.numeric(readline(
      prompt = "Enter the number of standard deviations to use for gridlines and filtering: "
    ))
    if (is.na(sd_factor) || sd_factor <= 0) {
      sd_factor <- 3
      cat("Invalid input or non-positive value. Using default of 3.\n")
    }
    
    xlines <- c(-sd_factor * pca$sdev[1], sd_factor * pca$sdev[1])
    ylines <- c(-sd_factor * pca$sdev[2], sd_factor * pca$sdev[2])
    
    filtering_active <- TRUE
    while (filtering_active) {
      p <- ggplot(pcaX, aes(PC1, PC2)) +
        geom_point() +
        geom_text(aes(label = rownames(pcaX)), hjust = 0.5, vjust = -0.5) +
        geom_vline(xintercept = xlines, linetype = "dashed", color = "red") +
        geom_hline(yintercept = ylines, linetype = "dashed", color = "red") +
        labs(x = paste0("PC1: ", pca.var.percent[1], "%"),
             y = paste0("PC2: ", pca.var.percent[2], "%"),
             title = "PCA plot of proteins with SD gridlines") +
        theme_minimal()
      print(p)
      
      cat("\nPlease specify cutoffs for filtering or enter 'q' to quit:\n")
      
      # PC1 min
      pc1_min_input <- trimws(readline(
        prompt = paste0("Enter lower cutoff for PC1 (suggestion: ", round(xlines[1], 2), "): ")
      ))
      if (pc1_min_input == "q") {
        cat("Returning to PCA menu.\n")
        return(list(ptx = as.data.frame(t(ptx_trans)), sinfo = sinfo, binfo = binfo))
      }
      pc1_min <- as.numeric(pc1_min_input)
      
      # PC1 max
      pc1_max_input <- trimws(readline(
        prompt = paste0("Enter upper cutoff for PC1 (suggestion: ", round(xlines[2], 2), "): ")
      ))
      if (pc1_max_input == "q") {
        cat("Returning to PCA menu.\n")
        return(list(ptx = as.data.frame(t(ptx_trans)), sinfo = sinfo, binfo = binfo))
      }
      pc1_max <- as.numeric(pc1_max_input)
      
      # PC2 min
      pc2_min_input <- trimws(readline(
        prompt = paste0("Enter lower cutoff for PC2 (suggestion: ", round(ylines[1], 2), "): ")
      ))
      if (pc2_min_input == "q") {
        cat("Returning to PCA menu.\n")
        return(list(ptx = as.data.frame(t(ptx_trans)), sinfo = sinfo, binfo = binfo))
      }
      pc2_min <- as.numeric(pc2_min_input)
      
      # PC2 max
      pc2_max_input <- trimws(readline(
        prompt = paste0("Enter upper cutoff for PC2 (suggestion: ", round(ylines[2], 2), "): ")
      ))
      if (pc2_max_input == "q") {
        cat("Returning to PCA menu.\n")
        return(list(ptx = as.data.frame(t(ptx_trans)), sinfo = sinfo, binfo = binfo))
      }
      pc2_max <- as.numeric(pc2_max_input)
      
      if (any(is.na(c(pc1_min, pc1_max, pc2_min, pc2_max)))) {
        cat("Invalid numeric input. Try again.\n")
        next
      }
      
      # Apply filtering
      filtered_indices <- which(pcaX$PC1 >= pc1_min & pcaX$PC1 <= pc1_max &
                                  pcaX$PC2 >= pc2_min & pcaX$PC2 <= pc2_max)
      if (length(filtered_indices) > 0) {
        ptx_trans <- ptx_trans[filtered_indices, , drop = FALSE]
        cat("Filtering applied. Remaining number of proteins: ", nrow(ptx_trans), "\n")
      } else {
        cat("Warning: No proteins met the cutoff criteria. No changes made.\n")
      }
      
      filtering_active <- FALSE
    }
    
    cat("Transposing matrix back to original state...\n")
    ptx <- as.data.frame(t(ptx_trans))
    rownames(ptx) <- rownames(sinfo)
    
    reproduce <- data.frame(
      SD = sd_factor,
      PC1low = pc1_min,
      PC1hi = pc1_max,
      PC2low = pc2_min,
      PC2hi = pc2_max,
      stringsAsFactors = FALSE
    )
    
  } else {
    # Reproduce branch
    sd_factor <- reproduce$SD
    xlines <- c(-sd_factor * pca$sdev[1], sd_factor * pca$sdev[1])
    ylines <- c(-sd_factor * pca$sdev[2], sd_factor * pca$sdev[2])
    
    pca_protein_outlier_plot <<- ggplot(pcaX, aes(PC1, PC2)) +
      geom_point() +
      geom_text(aes(label = rownames(pcaX)), hjust = 0.5, vjust = -0.5) +
      geom_vline(xintercept = xlines, linetype = "dashed", color = "red") +
      geom_hline(yintercept = ylines, linetype = "dashed", color = "red") +
      labs(x = paste0("PC1: ", pca.var.percent[1], "%"),
           y = paste0("PC2: ", pca.var.percent[2], "%"),
           title = "PCA plot of proteins with SD gridlines") +
      theme_minimal()
    
    pc1_min <- reproduce$PC1low
    pc1_max <- reproduce$PC1hi
    pc2_min <- reproduce$PC2low
    pc2_max <- reproduce$PC2hi
    
    filtered_indices <- which(pcaX$PC1 >= pc1_min & pcaX$PC1 <= pc1_max &
                                pcaX$PC2 >= pc2_min & pcaX$PC2 <= pc2_max)
    ptx_trans <- ptx_trans[filtered_indices, , drop = FALSE]
    ptx <- as.data.frame(t(ptx_trans))
    rownames(ptx) <- rownames(sinfo)
  }
  
  # Update binfo and assign to global environment
  binfo <- binfo[binfo$protein_name %in% colnames(ptx), , drop = FALSE]
  assign("select.ptx", ptx, envir = .GlobalEnv)
  assign("select.binfo", binfo, envir = .GlobalEnv)
  assign("reproducePCAproteinFilter", reproduce, envir = .GlobalEnv)
  
  return(list(ptx = ptx, sinfo = sinfo, binfo = binfo))
}

################################################################################ 
# Filter proteins based on frequency below LOD
################################################################################ 
filterfbLOD <- function(ptx, sinfo, binfo, reproduce = NULL) {
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      cutoff = NA,
      stringsAsFactors = FALSE
    )
    cat("\nEnter a cutoff value between 0 and 1 for frequency below LOD \n (proteins with freq_below_lod above this cutoff will be removed):\n")
    cutoff <- as.numeric(readline())
    reproduce$cutoff <- cutoff
    if (is.na(cutoff) || cutoff < 0 || cutoff > 1) {
      cat("Invalid cutoff provided; using default value of 0.9.\n")
      cutoff <- 0.9
    }
    # Here selected_proteins correspond to proteins where freq_below_lod is below the cutoff
    selected_proteins <- binfo$protein_name[binfo$freq_below_lod > cutoff]
    cat("\n", length(selected_proteins), "proteins found with freq_below_lod above the cutoff.\n")
    if (length(selected_proteins) > 0) {
      binfo <- binfo[!(binfo$protein_name %in% selected_proteins), ]
      ptx <- ptx[, !(colnames(ptx) %in% selected_proteins), drop = FALSE]
      cat("Proteins removed from binfo and ptx based on the cutoff.\n")
    } else {
      cat("No proteins were removed based on the cutoff.\n")
    } 
    assign("reproduceLODfilter", reproduce, envir = .GlobalEnv)
  } else {
    cutoff <- reproduce$cutoff
    selected_proteins <- binfo$protein_name[binfo$freq_below_lod > cutoff]
    binfo <- binfo[!(binfo$protein_name %in% selected_proteins), ]
    ptx <- ptx[, !(colnames(ptx) %in% selected_proteins), drop = FALSE]
  }
}

################################################################################
# Protein and sample correlation Heatmaps
################################################################################
analyze_correlation <- function(reproduce = NULL) {
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      corrType = NA,
      title = NA,
      fontSize = NA,
      savePlot = NA,
      plotVarName = NA,
      stringsAsFactors = FALSE
    )
    cat("\n--- Correlation Analysis ---\n")
    cat("Choose correlation type:\n")
    cat("  1 = Proteins (columns)\n")
    cat("  2 = Samples (rows)\n")
    corr_choice <- as.numeric(readline(prompt = "Enter your choice: "))
    reproduce$corrType <- corr_choice
    
    if (corr_choice == 1) {
      cat("\nCalculating correlation matrix for proteins...\n")
      cor_mat <- cor(select.ptx, use = "pairwise.complete.obs")
      
    } else if (corr_choice == 2) {
      cat("\nCalculating correlation matrix for samples...\n")
      cor_mat <- cor(t(select.ptx), use = "pairwise.complete.obs")
      
    } else {
      cat("\nInvalid correlation choice. Skipping correlation analysis.\n")
      return(invisible(NULL))
    }
    
    main_title <- readline(prompt = "Enter Title of Plot: ")
    fontsize <- as.integer(readline(prompt = "Enter Font-size: "))
    reproduce$title <- main_title
    reproduce$fontSize <- fontsize
  } else {
    corr_choice <- reproduce$corrType
    main_title <- reproduce$title
    fontsize <- reproduce$fontSize
    if (corr_choice == 1) {
      cat("\nCalculating correlation matrix for proteins...\n")
      cor_mat <- cor(select.ptx, use = "pairwise.complete.obs")
      
    } else if (corr_choice == 2) {
      cat("\nCalculating correlation matrix for samples...\n")
      cor_mat <- cor(t(select.ptx), use = "pairwise.complete.obs")
      
    } else {
      cat("\nInvalid correlation choice. Skipping correlation analysis.\n")
      return(invisible(NULL))
    }
    
  }
  
  # Generate a symmetric color scale so that 0 is white.
  # Determine the maximum absolute value in the correlation matrix.
  max_val <- max(abs(cor_mat), na.rm = TRUE)
  # Create 51 breaks (resulting in 50 intervals) spanning from -max_val to max_val.
  breaks <- seq(-max_val, max_val, length.out = 51)
  # Create a palette of 50 colors spanning from blue to white to red.
  my_colors <- colorRampPalette(c("blue", "white", "red"))(50)
  
  # Generate the heatmap
  corrmap <- pheatmap::pheatmap(
    cor_mat,
    color = my_colors,
    breaks = breaks,
    main = main_title,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_col = fontsize, 
    fontsize_row = fontsize,
    silent = TRUE
  )
  
  if (is.na(reproduce$savePlot)) {
    save_choice <- readline(prompt ="Do you want to save the plot? y/n: \n")
    reproduce$savePlot <- save_choice
    if (save_choice == "y") {
      var_name2 <- readline("Enter a name to save your plot object under: ")
      assign(var_name2, corrmap, .GlobalEnv ) 
      print(corrmap) 
    } else {
      cat("Plot not saved as a variable")
      print(corrmap)
    }
  } else {
    if (save_choice == "y") {
      var_name2 <- reproduce$plotVarName
      assign(var_name2, corrmap, .GlobalEnv ) 
      print(corrmap) 
    } else {
      cat("Plot not saved as a variable")
      print(corrmap)
    }
    
  }
  assign("reproduceCorrHeatmap", reproduce, envir = .GlobalEnv)
  cat("\nCorrelation heatmap generated successfully.\n")
  return(invisible(TRUE))
}
################################################################################
# Functions using OlinkAnalyze
################################################################################
# Runner function
olink_runner <- function(select.ptx = NULL,
                         select.sinfo = NULL,
                         select.binfo = NULL,
                         olink_data = NULL,
                         ensure_unique_ids = FALSE) {
  
  repeat {
    cat("\n--- Olink Runner Main Menu ---\n")
    cat("1: Create Olink Format\n")
    cat("2: Olink PCA Plot\n")
    cat("3: Olink QC Plot\n")
    cat("4: Olink ptx Distribution\n")
    cat("5: Remove Outliers\n")
    cat("6: Convert back to Wide Format (Optional)\n")
    cat("7: Exit\n")
    
    option <- as.integer(readline(prompt = "Enter your choice (1-7): "))
    
    if (is.na(option)) {
      cat("Invalid input. Please enter a number between 1 and 7.\n")
      next
    }
    
    if (option == 7) {
      cat("Exiting the Olink Runner. Goodbye!\n")
      break
    } else if (option == 1) {
      if (is.null(select.ptx) || is.null(select.sinfo) || is.null(select.binfo)) {
        cat("For option 1, please provide select.ptx, select.sinfo, and select.binfo.\n")
      } else {
        result <- convert_to_olink_long(select.ptx, select.sinfo, select.binfo, ensure_unique_ids)
        olink_data <- result  # update olink_data for subsequent operations
        cat("Olink Format created successfully.\n")
        print(result)
      }
      
    } else if (option == 2) {
      if (is.null(olink_data)) {
        cat("For option 2 (PCA Plot), please run 'Create Olink Format'.\n")
      } else {
        result <- olink_pca_interactive(olink_data)
        cat("PCA Plot function executed interactively.\n")
        print(result)
      }
      
    } else if (option == 3) {
      if (is.null(olink_data)) {
        cat("For option 3 (QC Plot), please run 'Create Olink Format'.\n")
      } else {
        result <- olink_qc_interactive(olink_data)
        cat("QC Plot function executed interactively.\n")
        print(result)
      }
      
    } else if (option == 4) {
      if (is.null(olink_data)) {
        cat("For option 4 (ptx Distribution), please run 'Create Olink Format'.\n")
      } else {
        result <- olink_ptx_distribution(olink_data)
        cat("ptx Distribution function executed.\n")
        print(result)
      }
      
    } else if (option == 5) {
      if (is.null(olink_data)) {
        cat("For option 5 (Remove Outliers), please run 'Create Olink Format'.\n")
      } else {
        result <- remove_outliers_interactive(olink_data, olink_cam_outliers = olink_cam_outliers
                                              , olink_imonc_outliers= olink_imonc_outliers)
        cat("Outlier removal function executed.\n")
        print(result)
      }
      
    } else if (option == 6) {
      if (is.null(olink_data)) {
        cat("For option 6 (Convert back to Wide Format), please run 'Create Olink Format'.\n")
      } else {
        result <- reverse_convert_to_olink_long(olink_data)
        olink_data <- result  # update olink_data
        cat("Data has been converted back to wide format.\n")
        print(result)
      }
      
    } else {
      cat("Invalid option. Please select a number between 1 and 7.\n")
    }
    
    readline(prompt = "\nPress Enter to return to the main menu...")
  }
}
################################################################################ 
# Convert to Long format
################################################################################ 
convert_to_olink_long <- function(select.ptx, select.sinfo, select.binfo, 
                                  ensure_unique_ids = FALSE) {
  
  #########################################################
  # Prepare and Pivot ptx Data (Wide to Long Format)
  #########################################################
  if (!"sample_id" %in% colnames(select.ptx)) {
    select.ptx <- select.ptx %>% tibble::rownames_to_column("sample_id")
  }
  
  ptx_long <- select.ptx %>%
    tidyr::pivot_longer(
      cols = -sample_id,
      names_to = "protein_name",
      values_to = "ptx"
    )
  
  #########################################################
  # Merge Clinical (Sample) Data
  #########################################################
  if (!"sample_id" %in% colnames(select.sinfo)) {
    select.sinfo <- select.sinfo %>% tibble::rownames_to_column("sample_id")
  }
  
  merged_data <- ptx_long %>%
    dplyr::left_join(select.sinfo, by = "sample_id")
  
  #########################################################
  # Merge Protein Metadata (binfo) 
  #########################################################
  if ("region" %in% colnames(select.binfo) && "region" %in% colnames(merged_data)) {
    merged_data <- merged_data %>%
      dplyr::left_join(select.binfo, by = c("protein_name", "region"), suffix = c("", ".binfo"))
  } else {
    merged_data <- merged_data %>%
      dplyr::left_join(select.binfo, by = "protein_name", suffix = c("", ".binfo"))
  }
  
  #########################################################
  # Create a Unified QC_Warning Column
  #########################################################
  qc_cam <- if ("qc_warn_CAM" %in% colnames(merged_data)) {
    merged_data$qc_warn_CAM
  } else if ("QC_Warning_CAM" %in% colnames(merged_data)) {
    merged_data$QC_Warning_CAM
  } else {
    rep(NA_character_, nrow(merged_data))
  }
  
  qc_imonc <- if ("qc_warn_IMONC" %in% colnames(merged_data)) {
    merged_data$qc_warn_IMONC
  } else if ("QC_Warning_IMONC" %in% colnames(merged_data)) {
    merged_data$QC_Warning_IMONC
  } else {
    rep(NA_character_, nrow(merged_data))
  }
  
  merged_data <- merged_data %>%
    dplyr::mutate(
      QC_Warning = dplyr::case_when(
        panel == "CAM"   ~ as.character(qc_cam),
        panel == "IMONC" ~ as.character(qc_imonc),
        TRUE             ~ NA_character_
      )
    )
  
  # Remove individual QC warning columns using any_of() so that non-existent columns are ignored
  merged_data <- merged_data %>%
    dplyr::select(-dplyr::any_of(c("qc_warn_CAM", "qc_warn_IMONC", "QC_Warning_CAM", "QC_Warning_IMONC")))
  
  #########################################################
  # Optionally Ensure Unique Sample IDs
  #########################################################
  if (ensure_unique_ids) {
    merged_data <- merged_data %>% 
      dplyr::mutate(sample_id = paste0(sample_id, "_", panel))
  }
  
  #########################################################
  # Final Renaming to Standardize Column Names and 
  # Convert Specific Columns to Character Datatype
  #########################################################
  final_data <- merged_data %>%
    { 
      if("freq_below_lod.binfo" %in% colnames(.)) {
        dplyr::rename(., 
                      SampleID      = sample_id,
                      Assay         = protein_name,
                      UniProt       = uniprot_id,
                      OlinkID       = olink_id,
                      Panel         = panel,
                      Panel_Version = panel_long,
                      Normalisation = normalisation,
                      Region        = region,
                      LOD           = lod,
                      MissingFreq   = `freq_below_lod.binfo`
        )
      } else if ("freq_below_lod" %in% colnames(.)) {
        dplyr::rename(., 
                      SampleID      = sample_id,
                      Assay         = protein_name,
                      UniProt       = uniprot_id,
                      OlinkID       = olink_id,
                      Panel         = panel,
                      Panel_Version = panel_long,
                      Normalisation = normalisation,
                      Region        = region,
                      LOD           = lod,
                      MissingFreq   = freq_below_lod
        )
      } else {
        dplyr::rename(., 
                      SampleID      = sample_id,
                      Assay         = protein_name,
                      UniProt       = uniprot_id,
                      OlinkID       = olink_id,
                      Panel         = panel,
                      Panel_Version = panel_long,
                      Normalisation = normalisation,
                      Region        = region,
                      LOD           = lod
        ) %>% dplyr::mutate(MissingFreq = NA)
      }
    } %>%
    dplyr::mutate(Index = dplyr::row_number()) %>%

    dplyr::mutate(dplyr::across(
      dplyr::any_of(c("mht_status", "CancerTime", "x_BC", 
                      "menopause_status", "smoking_status", 
                      "plate_id", "row", "column")), as.character))
  LongData <<- final_data
  return(final_data)
}
################################################################################ 
# Olink PCA script (study-specific)
################################################################################ 
olink_pca_interactive <- function(olink_long_data,
                                  outlierDefX = NULL,
                                  outlierDefY = NULL,
                                  outlierLines = TRUE,
                                  label_outliers = TRUE,
                                  quiet = TRUE,
                                  byPanel = TRUE) {
  # Check input: Ensure data is provided.
  if (is.null(olink_long_data)) {
    stop("Please provide the long-format Olink data (olink_long_data).")
  }
  
  # Prompt for outlierDefX if not supplied.
  if (is.null(outlierDefX)) {
    outlierDefX <- as.numeric(readline(prompt = "Enter a value for outlierDefX: "))
    if (is.na(outlierDefX))
      stop("Invalid input for outlierDefX.")
  }
  
  # Prompt for outlierDefY if not supplied.
  if (is.null(outlierDefY)) {
    outlierDefY <- as.numeric(readline(prompt = "Enter a value for outlierDefY: "))
    if (is.na(outlierDefY))
      stop("Invalid input for outlierDefY.")
  }
  
  # Filter out rows that contain "CONTROL" in the SampleID.
  filtered_data <- olink_long_data %>%
    dplyr::filter(!stringr::str_detect(SampleID, "CONTROL"))
  
  # Run the underlying PCA plotting function.
  pca_result <- olink_pca_plot(
    filtered_data,
    outlierDefX = outlierDefX,
    outlierDefY = outlierDefY,
    outlierLines = outlierLines,
    label_outliers = label_outliers,
    quiet = quiet,
    byPanel = byPanel
  )
  options(ggrepel.max.overlaps = 50)
  
  # If pca_result is a list with two panels, combine them.
  if (is.list(pca_result) && length(pca_result) == 2) {
    # Combine the two PCA plots side by side.
    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined_plot <- pca_result[[1]] + pca_result[[2]] + patchwork::plot_layout(ncol = 2)
      grob1 <- ggplotGrob(pca_result[[1]])
      grob2 <- ggplotGrob(pca_result[[2]])
      grid.newpage()
      grid.arrange(grob1, grob2, ncol = 2)
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      combined_plot <- gridExtra::grid.arrange(pca_result[[1]], pca_result[[2]], ncol = 2)
      grob1 <- ggplotGrob(pca_result[[1]])
      grob2 <- ggplotGrob(pca_result[[2]])
      grid.newpage()
      grid.arrange(grob1, grob2, ncol = 2)
    } else {
      warning("Neither patchwork nor gridExtra is available. Plots will be shown separately.")
      combined_plot <- NULL
    }
    
    # Extract outlier SampleIDs from each panel's data.
    cam_outliers <- pca_result[[1]]$data %>%
      dplyr::filter(Outlier == TRUE) %>%
      dplyr::select(SampleID) %>%
      dplyr::distinct()
    
    imonc_outliers <- pca_result[[2]]$data %>%
      dplyr::filter(Outlier == TRUE) %>%
      dplyr::select(SampleID) %>%
      dplyr::distinct()
    
    # Combine the outlier data frames with an indicator for the panel.
    combined_outliers <- dplyr::bind_rows(
      cam_outliers %>% dplyr::mutate(Panel = "CAM"),
      imonc_outliers %>% dplyr::mutate(Panel = "IMONC")
    )
    
  } else {  # If only one panel exists.
    if (is.list(pca_result)) {
      single_plot <- pca_result[[1]]
    } else {
      single_plot <- pca_result
    }
    combined_plot <- single_plot
    print(single_plot)
    
    # Extract outlier SampleIDs from the single panel's data.
    combined_outliers <- single_plot$data %>%
      dplyr::filter(Outlier == TRUE) %>%
      dplyr::select(SampleID) %>%
      dplyr::distinct() %>%
      dplyr::mutate(Panel = NA) 
    cam_outliers <- combined_outliers
    imonc_outliers <- NULL
  }
  
  # Assign the combined outlier result globally.
  olink_pca_res <<- combined_outliers
  olink_cam_outliers <<- cam_outliers
  olink_imonc_outliers <<- imonc_outliers
  
  
  # Return a list with the outlier data.
  return(list(
    cam_outliers    = cam_outliers,
    imonc_outliers  = imonc_outliers,
    global_outliers = combined_outliers 
  ))
}
################################################################################ 
# Olink QC plots (study-specific)
################################################################################ 
olink_qc_interactive <- function(olink_long_data,
                                 median_outlierDef = NULL,
                                 IQR_outlierDef = NULL,
                                 outlierLines = TRUE,
                                 label_outliers = TRUE) {
  # Ensure that the long-format data is provided.
  if (is.null(olink_long_data)) {
    stop("Please provide the long-format data (olink_long_data).")
  }
  
  # Prompt for median_outlierDef if not supplied.
  if (is.null(median_outlierDef)) {
    median_outlierDef <- as.numeric(readline(prompt = "Enter value for median_outlierDef: "))
    if (is.na(median_outlierDef)) {
      stop("Invalid input for median_outlierDef.")
    }
  }
  
  # Prompt for IQR_outlierDef if not supplied.
  if (is.null(IQR_outlierDef)) {
    IQR_outlierDef <- as.numeric(readline(prompt = "Enter value for IQR_outlierDef: "))
    if (is.na(IQR_outlierDef)) {
      stop("Invalid input for IQR_outlierDef.")
    }
  }
  
  # List available columns for color grouping.
  available_cols <- names(olink_long_data)
  cat("Available columns in the data for color grouping:\n")
  for (i in seq_along(available_cols)) {
    cat(sprintf("%d: %s\n", i, available_cols[i]))
  }
  
  # Prompt for selection of column number for color grouping
  col_num <- as.integer(readline(prompt = "Enter the column number for color grouping (enter 0 for none): "))
  if (is.na(col_num) || col_num < 0 || col_num > length(available_cols)) {
    stop("Invalid column number entered.")
  }
  
  if (col_num == 0) {
    chosenTrait <- NULL
    cat("No color grouping will be used.\n")
  } else {
    chosenTrait <- available_cols[col_num]
    cat(sprintf("You have chosen '%s' for color grouping.\n", chosenTrait))
  }
  
  # Call the underlying QC plotting function.
  # If a column was chosen, convert to a symbol.
  if (!is.null(chosenTrait)) {
    qc_result <- olink_qc_plot(
      olink_long_data,
      median_outlierDef = median_outlierDef,
      IQR_outlierDef = IQR_outlierDef,
      outlierLines = outlierLines,
      label_outliers = label_outliers,
      color_g = !!rlang::sym(chosenTrait)
    )
  } else {
    qc_result <- olink_qc_plot(
      olink_long_data,
      median_outlierDef = median_outlierDef,
      IQR_outlierDef = IQR_outlierDef,
      outlierLines = outlierLines,
      label_outliers = label_outliers
    )
  }
  
  # Extract outlier SampleIDs from the QC plot's data.
  qc_outliers <- qc_result$data %>%
    dplyr::filter(Outlier == TRUE) %>%
    dplyr::select(SampleID) %>%
    dplyr::distinct()
  
  # Assign the extracted outlier data globally as 'olink_qc_res'
  olink_qc_res <<- qc_outliers
  
  # Print the QC plot 
  print(qc_result)
  
  # Return a list with the QC plot result and the outlier data.
  return(list(
    qc_plot_result = qc_result,
    qcplot_outliers = qc_outliers
  ))
}
################################################################################ 
# Outlier removal script for the olink functions
################################################################################ 
remove_outliers_interactive <- function(olink_long_data,
                                        olink_cam_outliers,
                                        olink_imonc_outliers) {
  # Combine all outlier data frames and summarize occurrences
  all_outliers <<- dplyr::bind_rows(olink_qc_res, olink_cam_outliers, olink_imonc_outliers) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(occurrences = dplyr::n(), .groups = "drop")
  
  cat("Combined Outliers Summary:\n")
  print(all_outliers)
  cat("\n")
  
  # Interactive removal choices
  cat("Removal options:\n")
  cat("1: Remove samples flagged in exactly 3 plots\n")
  cat("2: Remove samples flagged in at least 2 plots\n")
  cat("3: Remove samples flagged in at least 1 plot\n")
  cat("4: Remove samples manually using SampleIDs\n")
  cat("5: Back to Main (no removal)\n")
  
  removal_choice <- as.integer(readline(prompt = "Enter your choice (1-5): "))
  
  # Determine which samples to remove based on the choice
  if (removal_choice == 1) {
    samples_to_remove <<- all_outliers %>% dplyr::filter(occurrences == 3) %>% dplyr::pull(SampleID)
    select.ptx <<- select.ptx[!rownames(select.ptx) %in% samples_to_remove,]
    select.sinfo <<- select.sinfo[rownames(select.sinfo) %in% rownames(select.ptx), , drop = FALSE]
  } else if (removal_choice == 2) {
    samples_to_remove <<- all_outliers %>% dplyr::filter(occurrences >= 2) %>% dplyr::pull(SampleID)
    select.ptx <<- select.ptx[!rownames(select.ptx) %in% samples_to_remove,]
    select.sinfo <<- select.sinfo[rownames(select.sinfo) %in% rownames(select.ptx), , drop = FALSE]
  } else if (removal_choice == 3) {
    samples_to_remove <<- all_outliers %>% dplyr::filter(occurrences >= 1) %>% dplyr::pull(SampleID)
    select.ptx <<- select.ptx[!rownames(select.ptx) %in% samples_to_remove,]
    select.sinfo <<- select.sinfo[rownames(select.sinfo) %in% rownames(select.ptx), , drop = FALSE]
  } else if (removal_choice == 4) {
    manual_input <- readline(prompt = "Enter SampleIDs to remove (comma-separated): ")
    samples_to_remove <<- unlist(strsplit(manual_input, split = ",\\s*"))
    select.ptx <<- select.ptx[!rownames(select.ptx) %in% samples_to_remove,]
    select.sinfo <<- select.sinfo[rownames(select.sinfo) %in% rownames(select.ptx), , drop = FALSE]
  } else if (removal_choice == 5) {
    cat("Returning to Main Menu without removing samples.\n")
    return(NULL)
  } else {
    cat("Invalid choice. Aborting removal process.\n")
    return(NULL)
  }
  
  # Remove the flagged samples from the long format data frame
  cleaned_data <- olink_long_data %>%
    dplyr::filter(!(SampleID %in% samples_to_remove))
  
  cat("\nRemoved the following samples:\n")
  print(samples_to_remove)
  cat("Number of remaining samples: ", nrow(cleaned_data), "\n")
  
  return(list(
    cleaned_data  = cleaned_data,
    removed_samples = samples_to_remove,
    outliers_summary = all_outliers
  ))
}
################################################################################ 
# Olink ptx distribution plots
################################################################################
olink_ptx_distribution <- function(olink_data) {
  # Check that data is provided
  if (is.null(olink_data)) {
    stop("Please provide the Olink long-format data (olink_data).")
  }
  
  # Prompt if to adjust data by Sample Median?
  adjust_choice <- tolower(readline(prompt = "Do you want to adjust ptx by sample median? (y/n): "))
  if (adjust_choice == "y") {
    # Calculate the median ptx for each SampleID
    median_ptx <- olink_data %>%
      dplyr::group_by(SampleID) %>%
      dplyr::summarise(Median_ptx = median(ptx, na.rm = TRUE), .groups = "drop")
    
    # Adjust ptx by subtracting the sample median
    olink_data <- olink_data %>%
      dplyr::inner_join(median_ptx, by = "SampleID") %>%
      dplyr::mutate(ptx = ptx - Median_ptx)
    
    cat("Data has been adjusted by subtracting the sample median from ptx values.\n")
  } else {
    cat("Using original ptx data.\n")
  }
  
  # Choose a Variable for Coloring the Plot
  available_cols <- names(olink_data)
  cat("Available columns in your data:\n")
  print(available_cols)
  
  chosenTrait <- readline(prompt = "Enter the column name to use for color grouping (or leave blank for no coloring): ")
  if (chosenTrait != "" && !(chosenTrait %in% available_cols)) {
    cat("Column not found. Proceeding with no color grouping.\n")
    chosenTrait <- NULL
  }
  
  # Determine How Many Samples per Plot
  n_samples <- as.integer(readline(prompt = "Enter the number of samples per plot: "))
  if (is.na(n_samples) || n_samples <= 0) {
    stop("Please enter a valid positive integer for the number of samples per plot.")
  }
  
  # Order Data and Partition by SampleID 
  olink_data <- olink_data %>% dplyr::arrange(SampleID)
  unique_ids <- unique(olink_data$SampleID)
  groups <- split(unique_ids, ceiling(seq_along(unique_ids) / n_samples))
  
  # Create Plots for Each Group Using olink_dist_plot
  plot_list <- list()
  for (i in seq_along(groups)) {
    ids <- groups[[i]]
    subset_data <- olink_data %>% dplyr::filter(SampleID %in% ids)
    
    if (!is.null(chosenTrait) && chosenTrait != "") {
      p <- olink_dist_plot(subset_data, color_g = chosenTrait)
    } else {
      p <- olink_dist_plot(subset_data)
    }
    
    plot_list[[i]] <- p
  }
  
  # Interactive Plot Selection 
  repeat {
    cat("\nAvailable plots:\n")
    for (i in seq_along(plot_list)) {
      cat(paste0(i, ": Plot for samples: ", paste(groups[[i]], collapse = ", "), "\n"))
    }
    
    choice <- readline(prompt = "Enter a plot number to view it (or type 'q' to quit): ")
    if (tolower(choice) == "q") {
      cat("Exiting plot selection.\n")
      break
    }
    
    choice_num <- as.integer(choice)
    if (is.na(choice_num) || choice_num < 1 || choice_num > length(plot_list)) {
      cat("Invalid selection. Please enter a valid plot number or 'q' to quit.\n")
    } else {
      # Display the chosen plot
      print(plot_list[[choice_num]])
    }
  }
  
  # Return plots for further use if needed.
  return(plot_list)
}
################################################################################ 
# Conversion back to wide format
################################################################################ 
reverse_convert_to_olink_long <- function(olink_long_data) {
  
  # Reconstruct ptx data in wide format
  cleaned.ptx <<- olink_long_data %>%
    dplyr::select(dplyr::all_of(c("SampleID", "Assay", "ptx"))) %>%
    pivot_wider(names_from = Assay, values_from = ptx)
  
  # Extract sample (clinical) information from sinfo
  cleaned.sinfo <<- olink_long_data %>%
    dplyr::distinct(SampleID, .keep_all = TRUE) %>%
    dplyr::select(-dplyr::all_of(c("Assay", "ptx", "UniProt", "OlinkID",
                                   "Panel", "Panel_Version", "Normalisation",
                                   "LOD", "MissingFreq", "Index")))
  
  # Extract protein metadata from binfo
  # Use the distinct Assay rows and keep only columns derived from binfo.
  cleaned.binfo <<- olink_long_data %>%
    dplyr::distinct(Assay, .keep_all = TRUE) %>%
    dplyr::select(dplyr::all_of(c("Assay", "UniProt", "OlinkID", "Panel",
                                  "Panel_Version", "Normalisation", "LOD", "MissingFreq")))
  
  return(list(cleaned.ptx = cleaned.ptx,
              cleaned.sinfo = cleaned.sinfo,
              cleaned.binfo = cleaned.binfo))
}

################################################################################ 
# Basic PCA filtering function (without OlinkAnalyze)
################################################################################
run_pca_filtering_normal <- function(ptx, sinfo, reproduce = NULL) {
  # Compute PCA
  pca <- prcomp(as.data.frame(ptx), scale. = FALSE)
  pcaX <- as.data.frame(pca$x, row.names = rownames(ptx))
  pca.var <- pca$sdev^2
  pca.var.percent <- round((pca.var / sum(pca.var)) * 100, 2)
  
  if (is.null(reproduce)) {
    # --- Prompt for SD factor ---
    sd_factor <- as.numeric(readline(
      prompt = "Enter the number of standard deviations to use for gridlines and filtering (default is 3): "
    ))
    if (is.na(sd_factor) || sd_factor <= 0) {
      sd_factor <- 3
      cat("Invalid input or non-positive value. Using default value of 3.\n")
    }
    
    xlines <- c(-sd_factor * pca$sdev[1], sd_factor * pca$sdev[1])
    ylines <- c(-sd_factor * pca$sdev[2], sd_factor * pca$sdev[2])
    
    # --- Plot ---
    p <- ggplot(pcaX, aes(PC1, PC2)) +
      geom_point() +
      geom_text(aes(label = rownames(pcaX)), hjust = 0.5, vjust = -0.5) +
      geom_vline(xintercept = xlines, linetype = "dashed", color = "red") +
      geom_hline(yintercept = ylines, linetype = "dashed", color = "red") +
      labs(x = paste0("PC1: ", pca.var.percent[1], " %"),
           y = paste0("PC2: ", pca.var.percent[2], " %"),
           title = "PCA Plot with Filtering") +
      theme_minimal()
    print(p)
    
    # --- Prompt cutoffs ---
    cat("\nPlease specify cutoffs for filtering or enter 'q' to quit filtering:\n")
    
    pc1_min_input <- trimws(readline(
      prompt = paste0("Enter lower cutoff for PC1 (suggestion: ", round(xlines[1], 2), "): ")
    ))
    if (pc1_min_input == "q") {
      cat("Returning to PCA menu.\n")
      return(list(ptx = ptx, sinfo = sinfo))
    }
    pc1_min <- as.numeric(pc1_min_input)
    
    pc1_max_input <- trimws(readline(
      prompt = paste0("Enter upper cutoff for PC1 (suggestion: ", round(xlines[2], 2), "): ")
    ))
    if (pc1_max_input == "q") {
      cat("Returning to PCA menu.\n")
      return(list(ptx = ptx, sinfo = sinfo))
    }
    pc1_max <- as.numeric(pc1_max_input)
    
    pc2_min_input <- trimws(readline(
      prompt = paste0("Enter lower cutoff for PC2 (suggestion: ", round(ylines[1], 2), "): ")
    ))
    if (pc2_min_input == "q") {
      cat("Returning to PCA menu.\n")
      return(list(ptx = ptx, sinfo = sinfo))
    }
    pc2_min <- as.numeric(pc2_min_input)
    
    pc2_max_input <- trimws(readline(
      prompt = paste0("Enter upper cutoff for PC2 (suggestion: ", round(ylines[2], 2), "): ")
    ))
    if (pc2_max_input == "q") {
      cat("Returning to PCA menu.\n")
      return(list(ptx = ptx, sinfo = sinfo))
    }
    pc2_max <- as.numeric(pc2_max_input)
    
    if (any(is.na(c(pc1_min, pc1_max, pc2_min, pc2_max)))) {
      cat("Invalid numeric input. Returning without filtering.\n")
      return(list(ptx = ptx, sinfo = sinfo))
    }
    
    # --- Apply filtering ---
    filtered_indices <- which(pcaX$PC1 >= pc1_min & pcaX$PC1 <= pc1_max &
                                pcaX$PC2 >= pc2_min & pcaX$PC2 <= pc2_max)
    if (length(filtered_indices) > 0) {
      ptx <- ptx[filtered_indices, , drop = FALSE]
      sinfo <- sinfo[rownames(sinfo) %in% rownames(ptx), , drop = FALSE]
      assign("select.ptx", ptx, .GlobalEnv)
      assign("select.sinfo", sinfo, .GlobalEnv)
      cat("Filtering applied. Remaining samples:", nrow(ptx), "\n")
    } else {
      cat("Warning: No samples met the cutoff criteria. No changes made.\n")
    }
    
    # Save reproducibility object
    reproduce <- data.frame(
      SD = sd_factor,
      PC1low = pc1_min,
      PC1hi = pc1_max,
      PC2low = pc2_min,
      PC2hi = pc2_max,
      stringsAsFactors = FALSE
    )
    
  } else {
    # --- Reproduce branch ---
    sd_factor <- reproduce$SD
    xlines <- c(-sd_factor * pca$sdev[1], sd_factor * pca$sdev[1])
    ylines <- c(-sd_factor * pca$sdev[2], sd_factor * pca$sdev[2])
    
    pca_plot_normal <<- ggplot(pcaX, aes(PC1, PC2)) +
      geom_point() +
      geom_text(aes(label = rownames(pcaX)), hjust = 0.5, vjust = -0.5) +
      geom_vline(xintercept = xlines, linetype = "dashed", color = "red") +
      geom_hline(yintercept = ylines, linetype = "dashed", color = "red") +
      labs(x = paste0("PC1: ", pca.var.percent[1], " %"),
           y = paste0("PC2: ", pca.var.percent[2], " %"),
           title = "PCA Plot with Filtering") +
      theme_minimal()
    
    pc1_min <- reproduce$PC1low
    pc1_max <- reproduce$PC1hi
    pc2_min <- reproduce$PC2low
    pc2_max <- reproduce$PC2hi
    
    filtered_indices <- which(pcaX$PC1 >= pc1_min & pcaX$PC1 <= pc1_max &
                                pcaX$PC2 >= pc2_min & pcaX$PC2 <= pc2_max)
    ptx <- ptx[filtered_indices, , drop = FALSE]
    sinfo <- sinfo[rownames(sinfo) %in% rownames(ptx), , drop = FALSE]
  }
  
  assign("reproducePCAfilterNormal", reproduce, envir = .GlobalEnv)
  return(list(ptx = ptx, sinfo = sinfo))
}

###############################
# Hierarchical Clustering
###############################

run_hierarchical_clustering <- function(ptx, sinfo) {
  cat("\nRunning Hierarchical Clustering...\n")
  
  Traits4color <- c("ageCat", "bmiCat", "menopause_status", "mht_status",
                                 "CancerTime", "x_BC", "alcoholCat",
                                 "smoking_status")
  temp.sinfo <- sinfo
  temp.sinfo[Traits4color] <- lapply(temp.sinfo[Traits4color],
                                     function(x) as.numeric(as.factor(x)))
  # Perform hierarchical clustering on the numeric data in ptx
  htr <- hclust(dist(ptx), method = "average")
  
  if (!requireNamespace("WGCNA", quietly = TRUE))
    stop("WGCNA package is required. Please install it using install.packages('WGCNA').")
  
  
  tr.colors <- WGCNA::numbers2colors(temp.sinfo[, Traits4color, drop = FALSE], signed = FALSE)
  WGCNA::plotDendroAndColors(htr, tr.colors,
                             groupLabels = names(temp.sinfo[, Traits4color, drop = FALSE]),
                             main = "Sample Dendrogram and Trait Heatmap")
  
  # Interactive loop for sample removal; removal will update the global variables ptx and sinfo
  repeat {
    cat("\nWhich sample to remove? Enter sample row name or 0 to exit removal:\n")
    sample_to_remove <- trimws(readline())
    
    if (sample_to_remove == "0") {
      cat("\nReturning to QC menu.\n")
      break
    } else {
      if (sample_to_remove %in% rownames(ptx)) {
        # Remove the sample(s) from local ptx and sinfo
        ptx <- ptx[!rownames(ptx) %in% sample_to_remove, , drop = FALSE]
        sinfo <- sinfo[rownames(sinfo) %in% rownames(ptx), , drop = FALSE]
        
        # Explicitly update the global variables using assign()
        assign("select.ptx", ptx, envir = .GlobalEnv)
        assign("select.sinfo", sinfo, envir = .GlobalEnv)
        
        cat("\nSample", sample_to_remove, "removed successfully.\n")
      } else {
        cat("\nInvalid sample selection. Please try again.\n")
      }
    }
  }
  
  return(list(ptx = ptx, sinfo = sinfo))
}
