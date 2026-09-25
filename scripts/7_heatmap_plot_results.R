# Repo root (run scripts from repo root, or set COPY_NUM_ROOT)
ROOT <- Sys.getenv("COPY_NUM_ROOT", unset = "")
if (!nzchar(ROOT)) {
  ROOT <- if (dir.exists("scripts") && dir.exists("data")) {
    normalizePath(".")
  } else if (dir.exists("../scripts") && dir.exists("../data")) {
    normalizePath("..")
  } else {
    normalizePath(".")
  }
}

## Katia Renault
## Visualizing results from copy number PGLMM analysis

source(file.path(ROOT, "scripts", "FUN_color_mappings.R"))
color_mapping <- create_color_mapping(species_data$order)
########################################################################################################
# 1. Top pathways/heatmap 
########################################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

###########################################
# A. Prepare data
###########################################
load_pathway_data <- function() {
  signatures_dir <- "file.path(ROOT, "data", "gsea")"
  pgls_dir <- file.path(ROOT, "results")
  signatures_files <- list.files(signatures_dir, pattern = "\\.csv$", full.names = TRUE)
  pgls_files <- list.files(pgls_dir, pattern = "\\pathways.csv$", full.names = TRUE)
  all_files <- c(signatures_files, pgls_files)
  cat("Found", length(all_files), "files for analysis\n")
  all_pathways_data <- data.frame(
    pathway = character(),
    dataset = character(),
    NES = numeric(),
    pval = numeric(),
    FDR = numeric(),
    stringsAsFactors = FALSE
  )
  for (file in all_files) {
    file_name <- gsub("\\.csv$", "", basename(file))
    if (file %in% signatures_files) file_name <- gsub("_", " ", file_name)
    pgls_labels <- c("Copy Number ML", "Copy Number MLres")
    if (file %in% pgls_files) file_name <- pgls_labels[match(file, pgls_files)]
    data <- read.csv(file)
    if ("pathway" %in% colnames(data)) {
      path_col <- "pathway"
      nes_col <- "NES"
      pval_col <- "pval"
      fdr_col <- if ("padj" %in% colnames(data)) "padj" else "pval"
    } else {
      potential_path_cols <- c("pathway", "Pathway", "term", "Term", "gene_set", "Gene_set", "description", "Description")
      potential_nes_cols <- c("NES", "nes", "enrichment_score", "EnrichmentScore", "coefficient", "Coefficient")
      potential_pval_cols <- c("pval", "p_value", "pvalue", "p.value", "P.Value", "P_Value")
      potential_fdr_cols <- c("padj", "adj.P.Val", "FDR", "fdr", "qvalue", "q_value", "adjusted_pvalue")
      path_col <- potential_path_cols[potential_path_cols %in% colnames(data)][1]
      nes_col <- potential_nes_cols[potential_nes_cols %in% colnames(data)][1]
      pval_col <- potential_pval_cols[potential_pval_cols %in% colnames(data)][1]
      fdr_col <- potential_fdr_cols[potential_fdr_cols %in% colnames(data)][1]
      if (is.na(path_col) || is.na(nes_col) || is.na(pval_col)) {
        warning(paste("Could not determine columns for file:", file))
        next
      }
      if (is.na(fdr_col)) {
        fdr_col <- pval_col
      }
    }
    temp_df <- data.frame(
      pathway = data[[path_col]],
      dataset = file_name,
      NES = data[[nes_col]],
      pval = data[[pval_col]],
      FDR = data[[fdr_col]],
      stringsAsFactors = FALSE
    )

    all_pathways_data <- rbind(all_pathways_data, temp_df)
  }
  return(all_pathways_data)
}

###########################################
# B. Identify top pathways
###########################################
select_top_pathways_all_files <- function(pathway_data, top_n = 20, selection_method = "significance") {
  cat("Selecting top", top_n, "pathways present across ALL files using method:", selection_method, "\n")
  total_datasets <- length(unique(pathway_data$dataset))
  cat("Total datasets:", total_datasets, "\n")
  pathways_in_all_files <- pathway_data %>%
    group_by(pathway) %>%
    summarise(n_datasets = n_distinct(dataset), .groups = "drop") %>%
    filter(n_datasets == total_datasets) %>%
    pull(pathway)
  cat("Pathways present in all", total_datasets, "files:", length(pathways_in_all_files), "\n")
  if (length(pathways_in_all_files) == 0) {
    stop("No pathways found that are present in all files!")
  }
  if (length(pathways_in_all_files) < top_n) {
    cat("Warning: Only", length(pathways_in_all_files), "pathways available, adjusting top_n\n")
    top_n <- length(pathways_in_all_files)
  }
  all_files_data <- pathway_data %>%
    filter(pathway %in% pathways_in_all_files)
  if (selection_method == "significance") {
    top_pathways <- all_files_data %>%
      group_by(pathway) %>%
      summarise(
        min_fdr = min(FDR, na.rm = TRUE),
        max_abs_nes = max(abs(NES), na.rm = TRUE),
        mean_abs_nes = mean(abs(NES), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(min_fdr, desc(max_abs_nes)) %>%
      head(top_n) %>%
      pull(pathway)
    
  } else if (selection_method == "effect_size") {
    top_pathways <- all_files_data %>%
      group_by(pathway) %>%
      summarise(
        mean_abs_nes = mean(abs(NES), na.rm = TRUE),
        max_abs_nes = max(abs(NES), na.rm = TRUE),
        min_fdr = min(FDR, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_abs_nes), desc(max_abs_nes), min_fdr) %>%
      head(top_n) %>%
      pull(pathway)
  } else if (selection_method == "consistency") {
    top_pathways <- all_files_data %>%
      group_by(pathway) %>%
      summarise(
        n_significant_datasets = sum(FDR <= 0.05, na.rm = TRUE),
        prop_significant = n_significant_datasets / total_datasets,
        mean_abs_nes = mean(abs(NES), na.rm = TRUE),
        min_fdr = min(FDR, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(prop_significant), desc(n_significant_datasets), desc(mean_abs_nes), min_fdr) %>%
      head(top_n) %>%
      pull(pathway)
  } else if (selection_method == "combined") {
    top_pathways <- all_files_data %>%
      group_by(pathway) %>%
      summarise(
        min_fdr = min(FDR, na.rm = TRUE),
        mean_abs_nes = mean(abs(NES), na.rm = TRUE),
        max_abs_nes = max(abs(NES), na.rm = TRUE),
        n_significant = sum(FDR <= 0.05, na.rm = TRUE),
        prop_significant = n_significant / total_datasets,
        .groups = "drop"
      ) %>%
      mutate(
        norm_fdr = 1 - (min_fdr - min(min_fdr, na.rm = TRUE)) / (max(min_fdr, na.rm = TRUE) - min(min_fdr, na.rm = TRUE) + 1e-10),
        norm_nes = (mean_abs_nes - min(mean_abs_nes, na.rm = TRUE)) / (max(mean_abs_nes, na.rm = TRUE) - min(mean_abs_nes, na.rm = TRUE) + 1e-10),
        norm_consistency = prop_significant,
        combined_score = norm_fdr + norm_nes + norm_consistency
      ) %>%
      arrange(desc(combined_score)) %>%
      head(top_n) %>%
      pull(pathway)
    
  } else if (selection_method == "variance") {
    top_pathways <- all_files_data %>%
      group_by(pathway) %>%
      summarise(
        nes_variance = var(NES, na.rm = TRUE),
        mean_abs_nes = mean(abs(NES), na.rm = TRUE),
        min_fdr = min(FDR, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(nes_variance), desc(mean_abs_nes), min_fdr) %>%
      head(top_n) %>%
      pull(pathway)
  }
  cat("Selected pathways (present in all files):\n")
  for(i in 1:length(top_pathways)) {
    cat(i, ". ", top_pathways[i], "\n", sep = "")
  }
  cat("\nSummary statistics for selected pathways:\n")
  summary_stats <- all_files_data %>%
    filter(pathway %in% top_pathways) %>%
    group_by(pathway) %>%
    summarise(
      min_fdr = min(FDR, na.rm = TRUE),
      mean_abs_nes = mean(abs(NES), na.rm = TRUE),
      n_significant = sum(FDR <= 0.05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(match(pathway, top_pathways))
  
  print(summary_stats)
  
  return(top_pathways)
}

###########################################
# C. Make heatmap
###########################################
create_top_pathways_heatmap <- function(pathway_data, top_pathways) {
  filtered_data <- pathway_data %>%
    filter(pathway %in% top_pathways)
  all_datasets <- unique(pathway_data$dataset)
  complete_dataset <- expand.grid(
    pathway = top_pathways,
    dataset = all_datasets,
    stringsAsFactors = FALSE
  )
  filtered_data <- left_join(complete_dataset, filtered_data, by = c("pathway", "dataset"))
  filtered_data$display_name <- gsub("_", " ", filtered_data$pathway)
  filtered_data$significant <- case_when(
    is.na(filtered_data$FDR) ~ "",
    filtered_data$FDR <= 0.001 ~ "***",
    filtered_data$FDR <= 0.01 ~ "**",
    filtered_data$FDR <= 0.05 ~ "*",
    TRUE ~ ""
  )
  pathway_order <- filtered_data %>%
    group_by(pathway) %>%
    summarise(mean_nes = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_nes)) %>%
    pull(pathway)
  filtered_data$display_name <- factor(filtered_data$display_name, 
                                       levels = gsub("_", " ", pathway_order))
  p <- ggplot(filtered_data, aes(x = dataset, y = display_name, fill = NES)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = significant), size = 3, fontface = "bold", color = "black") +
    scale_fill_gradient2(
      low = "#4575b4", 
      mid = "white", 
      high = "#d73027",
      midpoint = 0,
      na.value = "grey90",
      name = "Normalized\nEnrichment\nScore"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(title = paste("Top", length(top_pathways), "pathways across datasets"),
         subtitle = "Significance levels: * FDR ≤ 0.05, ** FDR ≤ 0.01, *** FDR ≤ 0.001")
  
  return(p)
}

###########################################
# Main runner
###########################################
run_top_pathways_analysis <- function(top_n = 20, selection_method = "significance") {
  cat("=== TOP", top_n, "PATHWAYS ANALYSIS (PRESENT IN ALL FILES) ===\n")
  cat("Selection method:", selection_method, "\n\n")
  cat("Loading pathway data...\n")
  pathway_data <- load_pathway_data()
  if (is.null(pathway_data) || nrow(pathway_data) == 0) {
    stop("No pathway data found.")
  }
  total_datasets <- length(unique(pathway_data$dataset))
  total_pathways <- length(unique(pathway_data$pathway))
  cat("Total pathways found:", total_pathways, "\n")
  cat("Total datasets:", total_datasets, "\n")
  pathways_coverage <- pathway_data %>%
    group_by(pathway) %>%
    summarise(n_datasets = n_distinct(dataset), .groups = "drop") %>%
    count(n_datasets, name = "n_pathways")
  cat("\nPathway coverage across datasets:\n")
  print(pathways_coverage)
  pathways_in_all <- pathway_data %>%
    group_by(pathway) %>%
    summarise(n_datasets = n_distinct(dataset), .groups = "drop") %>%
    filter(n_datasets == total_datasets) %>%
    nrow()
  cat("Pathways present in ALL", total_datasets, "files:", pathways_in_all, "\n\n")
  if (pathways_in_all == 0) {
    stop("No pathways are present in all files! Consider using a lower coverage requirement.")
  }
  if (pathways_in_all < top_n) {
    cat("⚠️  Warning: Only", pathways_in_all, "pathways available in all files, but requested", top_n, "\n")
    cat("Will show all", pathways_in_all, "available pathways.\n\n")
  }
  top_pathways <- select_top_pathways_all_files(pathway_data, top_n, selection_method)
  cat("\nCreating heatmap...\n")
  heatmap <- create_top_pathways_heatmap(pathway_data, top_pathways)
  return(list(
    heatmap = heatmap, 
    data = pathway_data, 
    top_pathways = top_pathways,
    total_available_in_all_files = pathways_in_all
  ))
}
cat("🔍 Diagnostic: Checking pathway coverage across files\n")
diagnostic_data <- load_pathway_data()
coverage_summary <- diagnostic_data %>%
  group_by(pathway) %>%
  summarise(n_datasets = n_distinct(dataset), .groups = "drop") %>%
  count(n_datasets, name = "n_pathways") %>%
  arrange(desc(n_datasets))
cat("Pathway coverage distribution:\n")
print(coverage_summary)
total_files <- length(unique(diagnostic_data$dataset))
pathways_in_all <- sum(coverage_summary$n_pathways[coverage_summary$n_datasets == total_files])
cat("\nTotal files:", total_files, "\n")
cat("Pathways present in ALL files:", pathways_in_all, "\n\n")
results_effect <- run_top_pathways_analysis(top_n = 20, selection_method = "effect_size")
print(results_effect$heatmap)
results_combined <- run_top_pathways_analysis(top_n = 30, selection_method = "combined")
print(results_combined$heatmap)

########################################################################################################
# 2. Overall correlation between duplication and expression data
########################################################################################################

create_all_pathways_correlation_heatmap <- function(pathway_data, correlation_method = "pearson") {
  cat("Creating dataset correlation heatmap using ALL pathways present in ALL files...\n")
  cat("Method:", correlation_method, "\n\n")
  all_datasets <- unique(pathway_data$dataset)
  n_datasets <- length(all_datasets)
  cat("Total datasets:", n_datasets, "\n")
  if (n_datasets < 2) {
    stop("Need at least 2 datasets to compute correlations!")
  }
  pathways_in_all_files <- pathway_data %>%
    group_by(pathway) %>%
    summarise(n_datasets = n_distinct(dataset), .groups = "drop") %>%
    filter(n_datasets == n_datasets) %>%
    pull(pathway)
  cat("Pathways present in ALL", n_datasets, "files:", length(pathways_in_all_files), "\n")
  if (length(pathways_in_all_files) == 0) {
    stop("No pathways found that are present in all files!")
  }
  filtered_data <- pathway_data %>%
    filter(pathway %in% pathways_in_all_files)
  correlation_data <- filtered_data %>%
    select(pathway, dataset, NES) %>%
    pivot_wider(names_from = dataset, values_from = NES) %>%
    column_to_rownames("pathway")
  correlation_data_complete <- correlation_data[complete.cases(correlation_data), ]
  cat("Pathways after removing incomplete cases:", nrow(correlation_data_complete), "\n")
  cat("Original pathways in all files:", nrow(correlation_data), "\n")
  if (nrow(correlation_data_complete) < 3) {
    stop("Not enough complete pathway data for reliable correlation analysis!")
  }
  correlation_matrix <- cor(correlation_data_complete, method = correlation_method, use = "complete.obs")
  n_pathways <- nrow(correlation_data_complete)
  p_values <- matrix(NA, nrow = ncol(correlation_matrix), ncol = ncol(correlation_matrix))
  colnames(p_values) <- colnames(correlation_matrix)
  rownames(p_values) <- rownames(correlation_matrix)
  
  for (i in 1:ncol(correlation_matrix)) {
    for (j in 1:ncol(correlation_matrix)) {
      if (i != j) {
        test_result <- cor.test(correlation_data_complete[, i], correlation_data_complete[, j], method = correlation_method)
        p_values[i, j] <- test_result$p.value
      } else {
        p_values[i, j] <- 0
      }
    }
  }
  upper_tri_cors <- correlation_matrix[upper.tri(correlation_matrix)]
  cat("\n=== CORRELATION SUMMARY ===\n")
  cat("Pathways used for correlation:", nrow(correlation_data_complete), "\n")
  cat("Mean correlation:", round(mean(upper_tri_cors), 3), "\n")
  cat("Median correlation:", round(median(upper_tri_cors), 3), "\n")
  cat("Min correlation:", round(min(upper_tri_cors), 3), "\n")
  cat("Max correlation:", round(max(upper_tri_cors), 3), "\n\n")
  cat("Individual dataset correlations:\n")
  for (i in 1:(ncol(correlation_matrix)-1)) {
    for (j in (i+1):ncol(correlation_matrix)) {
      dataset1 <- colnames(correlation_matrix)[i]
      dataset2 <- colnames(correlation_matrix)[j]
      cor_value <- correlation_matrix[i, j]
      
      cat("  ", dataset1, " vs ", dataset2, ": ", round(cor_value, 3), "\n", sep = "")
    }
  }
  correlation_df <- as.data.frame(correlation_matrix) %>%
    rownames_to_column("Dataset1") %>%
    pivot_longer(cols = -Dataset1, names_to = "Dataset2", values_to = "Correlation") %>%
    mutate(
      Dataset1 = factor(Dataset1, levels = colnames(correlation_matrix)),
      Dataset2 = factor(Dataset2, levels = colnames(correlation_matrix))
    )
  p_values_df <- as.data.frame(p_values) %>%
    rownames_to_column("Dataset1") %>%
    pivot_longer(cols = -Dataset1, names_to = "Dataset2", values_to = "P_Value") %>%
    mutate(
      Dataset1 = factor(Dataset1, levels = colnames(correlation_matrix)),
      Dataset2 = factor(Dataset2, levels = colnames(correlation_matrix))
    )
  plot_data <- left_join(correlation_df, p_values_df, by = c("Dataset1", "Dataset2"))
  plot_data$significance <- case_when(
    plot_data$Dataset1 == plot_data$Dataset2 ~ "",
    plot_data$P_Value <= 0.001 ~ "***",
    plot_data$P_Value <= 0.01 ~ "**",
    plot_data$P_Value <= 0.05 ~ "*",
    TRUE ~ ""
  )
  p <- ggplot(plot_data, aes(x = Dataset1, y = Dataset2, fill = Correlation)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Correlation)),
              size = 3, fontface = "bold") +
    geom_text(aes(label = significance), 
              size = 4, fontface = "bold", color = "black", vjust = -0.8) +
    scale_fill_gradient2(
      low = "#4575b4", 
      mid = "white", 
      high = "#d73027",
      midpoint = 0,
      na.value = "grey90",
      name = "Correlation",
      limits = c(-1, 1)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      aspect.ratio = 1
    ) +
    labs(title = "Dataset Correlation Matrix",
         subtitle = paste("All pathways present in ALL files (n=", nrow(correlation_data_complete), ") - ",
                          str_to_title(correlation_method), " correlation\nSignificance levels: * p ≤ 0.05, ** p ≤ 0.01, *** p ≤ 0.001", sep = ""))
  return(list(
    plot = p,
    correlation_matrix = correlation_matrix,
    p_values = p_values,
    pathways_used = pathways_in_all_files,
    correlation_data = correlation_data_complete,
    n_pathways = nrow(correlation_data_complete),
    summary_stats = data.frame(
      n_pathways = nrow(correlation_data_complete),
      mean_correlation = mean(upper_tri_cors),
      median_correlation = median(upper_tri_cors),
      min_correlation = min(upper_tri_cors),
      max_correlation = max(upper_tri_cors),
      sd_correlation = sd(upper_tri_cors)
    )
  ))
}

run_all_pathways_correlation_analysis <- function(correlation_method = "pearson") {
  cat("=== DATASET CORRELATION ANALYSIS (ALL PATHWAYS IN ALL FILES) ===\n")
  cat("Method:", correlation_method, "\n\n")
  pathway_data <- load_pathway_data()
  if (is.null(pathway_data) || nrow(pathway_data) == 0) {
    stop("No pathway data found.")
  }
  total_datasets <- length(unique(pathway_data$dataset))
  pathway_coverage <- pathway_data %>%
    group_by(pathway) %>%
    summarise(n_datasets = n_distinct(dataset), .groups = "drop") %>%
    count(n_datasets, name = "n_pathways") %>%
    arrange(desc(n_datasets))
  cat("Pathway coverage distribution:\n")
  print(pathway_coverage)
  pathways_in_all <- pathway_coverage %>%
    filter(n_datasets == total_datasets) %>%
    pull(n_pathways)
  if (length(pathways_in_all) == 0) pathways_in_all <- 0
  cat("\nPathways present in ALL", total_datasets, "files:", pathways_in_all, "\n\n")
  if (pathways_in_all == 0) {
    stop("No pathways are present in all files! Check your data.")
  }
  correlation_results <- create_all_pathways_correlation_heatmap(pathway_data, correlation_method)
  return(correlation_results)
}

cat("🔗 Running all-pathways correlation analysis...\n\n")
all_pathways_results <- run_all_pathways_correlation_analysis("spearman")
print(all_pathways_results$plot)
all_pathways_results_pearson <- run_all_pathways_correlation_analysis("pearson")
print(all_pathways_results_pearson$plot)
