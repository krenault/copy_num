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

##################################
## Ranked method gene highlight ##
##################################

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(ggrepel)

# Jewel palette function
get_jewel_palette <- function(n) {
  sophisticated_jewels_palette <- c(
    "#9b383a", "#64a590", "#AA4839", "#AF6C36", "#945a87",
    "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
    "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
    "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
    "#c07a7a", "#8a7c6d", "#5c4d8b", "#4878a0", "#3D7B68"
  )
  if (missing(n)) {
    return(sophisticated_jewels_palette)
  } else {
    if (n > length(sophisticated_jewels_palette)) {
      return(sophisticated_jewels_palette[1:length(sophisticated_jewels_palette)])
    } else {
      return(sophisticated_jewels_palette[1:n])
    }
  }
}

# Function to identify top N outliers
identify_top_outliers <- function(x, y, n_outliers = 10) {
  # Method 1: Outliers based on individual variables (Z-score)
  x_scores <- abs(scale(x))
  y_scores <- abs(scale(y))
  
  # Method 2: Outliers based on regression residuals
  model <- lm(y ~ x)
  residuals_abs <- abs(residuals(model))
  residuals_scaled <- residuals_abs / sd(residuals_abs)
  
  # Method 3: Outliers based on leverage 
  leverage_values <- hatvalues(model)
  leverage_scaled <- leverage_values / mean(leverage_values)
  
  # Create composite outlier score for each point
  composite_score <- x_scores + y_scores + residuals_scaled + leverage_scaled
  
  # Get indices of top N outliers based on composite score
  top_outliers <- order(composite_score, decreasing = TRUE)[1:min(n_outliers, length(x))]
  
  return(list(
    top_outliers = top_outliers,
    composite_scores = composite_score
  ))
}

# Load data
# PGLS results to get gene list
file_path <- file.path(ROOT, "results", "max_longevity_simple_phylo_ranking_results_w_or.csv")
pgls_results <- read.csv(file_path)

# DEBUG: Print column names to identify the correct gene column
cat("Column names in pgls_results:\n")
print(colnames(pgls_results))
cat("\n")

gene_copy_data <- read.csv(file.path(ROOT, "data", "All_Species_Orthologous_CopyNumber_Annotated.tsv"), row.names = "t_gene", sep = '\t')
rownames(gene_copy_data) <- gene_copy_data$t_symbol
gene_copy_data <- gene_copy_data %>% dplyr::select(-t_symbol)

# Filter out OR genes
# cat("Filtering out OR genes...\n")
# or_genes <- rownames(gene_copy_data)[grepl("^OR", rownames(gene_copy_data))]
# cat("  Removing", length(or_genes), "OR genes\n")
# gene_copy_data <- gene_copy_data[!grepl("^OR", rownames(gene_copy_data)), ]

# Read metadata
metadata <- read.csv(file.path(ROOT, "data", "raxml_final_metadata_revised.csv"))

# Clean up gene copy data if needed
if("X" %in% colnames(gene_copy_data)) {
  gene_copy_data$X <- NULL
}

# Prepare metadata
metadata <- metadata %>%
  filter(!is.na(MLres)) %>%
  mutate(log_lifespan = log(maximum_longevity_y))

# Check if 'order' column exists, if not try 'Order'
if (!"order" %in% names(metadata)) {
  if ("Order" %in% names(metadata)) {
    metadata$order <- metadata$Order
  } else {
    cat("Warning: Neither 'order' nor 'Order' column found in metadata. Creating dummy order.\n")
    metadata$order <- "Unknown"
  }
}

# Specify genes you want to plot
specific_genes <- head(pgls_results %>% arrange(pgls_results$adj_p_value))$Gene

# Alternative: Use top significant genes from PGLS results
# Get top 5 most significant genes
if (nrow(pgls_results) > 0) {
  # First, let's check what the first few rows look like
  cat("First few rows of pgls_results:\n")
  print(head(pgls_results, 3))
  cat("\n")
  
  # Try to identify the gene column - common names are: gene, Gene, gene_name, symbol, etc.
  possible_gene_columns <- c("gene", "Gene", "gene_name", "symbol", "t_symbol", "gene_symbol")
  gene_column <- intersect(possible_gene_columns, colnames(pgls_results))[1]
  
  if (is.na(gene_column)) {
    # If no standard gene column found, use the first column (often contains gene names)
    gene_column <- colnames(pgls_results)[1]
    cat("Warning: No standard gene column found. Using first column:", gene_column, "\n")
  } else {
    cat("Using gene column:", gene_column, "\n")
  }
  
  # Sort by p-value and get top genes
  top_genes <- pgls_results %>%
    filter(!is.na(p_value)) %>%
    arrange(p_value) %>%
    slice_head(n = 5) %>%
    pull(!!sym(gene_column))  # Use !!sym() to handle column name as string
  
  # Use top genes if specific genes aren't in the data
  genes_in_data <- intersect(specific_genes, rownames(gene_copy_data))
  if (length(genes_in_data) == 0) {
    cat("Specified genes not found. Using top significant genes instead.\n")
    specific_genes <- top_genes
  } else {
    specific_genes <- genes_in_data
  }
  
  cat("Analyzing genes:", paste(specific_genes, collapse = ", "), "\n")
}

# Function to create gene plot
create_gene_plot <- function(gene_name, n_outliers = 10) {
  # Check if gene exists in data
  if (!gene_name %in% rownames(gene_copy_data)) {
    cat("Gene", gene_name, "not found in copy number data\n")
    return(NULL)
  }
  
  # Get copy number data for this gene
  copy_data <- as.numeric(gene_copy_data[gene_name, ])
  species_names <- colnames(gene_copy_data)
  
  # Create data frame for visualization
  plot_df <- data.frame(
    species = species_names,
    copy_number = copy_data
  ) %>%
    filter(!is.na(copy_number)) %>%
    # Match with metadata
    inner_join(metadata, by = c("species" = "Scientific_name"))
  
  # Skip if no data
  if (nrow(plot_df) == 0) {
    cat("No data available for gene:", gene_name, "\n")
    return(NULL)
  }
  
  # Clean species names (remove underscores)
  plot_df$species_clean <- gsub("_", " ", plot_df$species)
  
  # Set up color mapping for orders
  unique_orders <- unique(plot_df$order)
  color_palette <- get_jewel_palette(length(unique_orders))
  names(color_palette) <- unique_orders
  
  # Identify top N outliers
  outlier_results <- identify_top_outliers(plot_df$copy_number, plot_df$log_lifespan, n_outliers)
  outlier_indices <- outlier_results$top_outliers
  
  # Get PGLS results for this gene
  gene_pgls_result <- pgls_results[pgls_results[[gene_column]] == gene_name, ]
  
  # Extract coefficient and robust p-value if available
  if (nrow(gene_pgls_result) > 0) {
    # Safely extract coefficient
    if ("coefficient" %in% colnames(gene_pgls_result)) {
      coefficient <- gene_pgls_result$coefficient[1]
    } else {
      coefficient <- NA
    }
    
    # Safely extract p-values
    p_value_robust <- NULL
    if ("p_value_robust" %in% colnames(gene_pgls_result)) {
      p_value_robust <- gene_pgls_result$p_value_robust[1]
    }
    
    p_value_regular <- NULL
    if ("p_value" %in% colnames(gene_pgls_result)) {
      p_value_regular <- gene_pgls_result$p_value[1]
    }
    
    # Use robust p-value if available and not NA, otherwise use regular p-value
    if (!is.null(p_value_robust) && length(p_value_robust) > 0 && !is.na(p_value_robust)) {
      p_value_to_show <- p_value_robust
    } else if (!is.null(p_value_regular) && length(p_value_regular) > 0 && !is.na(p_value_regular)) {
      p_value_to_show <- p_value_regular
    } else {
      p_value_to_show <- NA
    }
  } else {
    # Fallback to correlation if PGLS results not found
    cor_test <- cor.test(plot_df$copy_number, plot_df$log_lifespan)
    coefficient <- cor_test$estimate
    p_value_to_show <- cor_test$p.value
  }
  
  # Create the plot
  p <- ggplot(plot_df, aes(x = copy_number, y = log_lifespan)) +
    # Points colored by order
    geom_point(aes(color = order, fill = order),
               size = 3, alpha = 0.7, shape = 21, stroke = 0.6) +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = color_palette) +
    
    # Add regression line
    geom_smooth(method = "lm", se = TRUE, color = "gray40", alpha = 0.3, linetype = 'dashed') +
    
    # Label only top N outliers (using cleaned species names)
    geom_text_repel(data = plot_df[outlier_indices, ],
                    aes(label = species_clean),
                    size = 4,
                    box.padding = 0.5,
                    point.padding = 0.2,
                    segment.color = "gray60",
                    segment.alpha = 0.7,
                    max.overlaps = Inf,
                    force = 2) +
    
    # Labels and theme
    labs(title = paste("Gene:", gene_name),
         subtitle = paste("Coefficient =", 
                          ifelse(is.na(coefficient), "NA", round(coefficient, 3)),
                          "p =", 
                          ifelse(is.na(p_value_to_show), "NA", formatC(p_value_to_show, format = "e", digits = 2))),
         x = "Copy number",
         y = "Mass-corrected longevity",
         color = "Order",
         fill = "Order") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 13),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  
  # Print outlier summary
  cat("\nTop", n_outliers, "outliers for", gene_name, ":\n")
  if (length(outlier_indices) > 0) {
    outlier_summary <- plot_df[outlier_indices, c("species_clean", "order")] %>%
      mutate(copy_number = plot_df$copy_number[outlier_indices],
             log_lifespan = plot_df$log_lifespan[outlier_indices]) %>%
      arrange(desc(outlier_results$composite_scores[outlier_indices]))
    colnames(outlier_summary)[1] <- "species"
    print(outlier_summary)
  }
  
  return(p)
}

# Create plots for specified genes
gene_plots <- list()
for (gene in specific_genes) {
  cat("Processing gene:", gene, "\n")
  plot <- create_gene_plot(gene, n_outliers = 10)
  if (!is.null(plot)) {
    gene_plots[[gene]] <- plot
  }
}

# Display plots
for (i in 1:length(gene_plots)) {
  print(gene_plots[[i]])
}

# Save individual plots (uncomment to save)
# for (gene in names(gene_plots)) {
#   filename <- paste0("plots/",
#                      gene, "_copy_vs_lifespan.pdf")
#   ggsave(filename, gene_plots[[gene]], width = 12, height = 8)
# }

# Create combined plot if you have multiple genes
if (length(gene_plots) > 1) {
  library(patchwork)
  combined_plot <- wrap_plots(gene_plots, ncol = 2) +
    plot_annotation(
      title = "Gene Copy Number vs Maximum Longevity",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
  
  print(combined_plot)
  
  # Save combined plot (uncomment to save)
  # ggsave("plots/genes_combined_plots.pdf",
  #        combined_plot, width = 16, height = 12)
}


##########################
##### PGLMM method #######
##########################

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(ggrepel)

# Jewel palette function
get_jewel_palette <- function(n) {
  sophisticated_jewels_palette <- c(
    "#9b383a", "#64a590", "#AA4839", "#AF6C36", "#945a87",
    "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
    "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
    "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
    "#c07a7a", "#8a7c6d", "#5c4d8b", "#4878a0", "#3D7B68"
  )
  if (missing(n)) {
    return(sophisticated_jewels_palette)
  } else {
    if (n > length(sophisticated_jewels_palette)) {
      return(sophisticated_jewels_palette[1:length(sophisticated_jewels_palette)])
    } else {
      return(sophisticated_jewels_palette[1:n])
    }
  }
}

# Function to identify top N outliers
identify_top_outliers <- function(x, y, n_outliers = 5) {
  # Method 1: Outliers based on individual variables (Z-score)
  x_scores <- abs(scale(x))
  y_scores <- abs(scale(y))
  
  # Method 2: Outliers based on regression residuals
  model <- lm(y ~ x)
  residuals_abs <- abs(residuals(model))
  residuals_scaled <- residuals_abs / sd(residuals_abs)
  
  # Method 3: Outliers based on leverage 
  leverage_values <- hatvalues(model)
  leverage_scaled <- leverage_values / mean(leverage_values)
  
  # Create composite outlier score for each point
  composite_score <- x_scores + y_scores + residuals_scaled + leverage_scaled
  
  # Get indices of top N outliers based on composite score
  top_outliers <- order(composite_score, decreasing = TRUE)[1:min(n_outliers, length(x))]
  
  return(list(
    top_outliers = top_outliers,
    composite_scores = composite_score
  ))
}

# Load data
# PGLS results to get gene list
file_path <- file.path(ROOT, "results", "max_longevity_zero_filtered_poisson_pglmm_20250702_110253.csv")
pgls_results <- read.csv(file_path)

# DEBUG: Print column names to identify the correct gene column
cat("Column names in pgls_results:\n")
print(colnames(pgls_results))
cat("\n")

gene_copy_data <- read.csv(file.path(ROOT, "data", "All_Species_Orthologous_CopyNumber_Annotated.tsv"),
                           row.names = "t_symbol", sep = '\t')
gene_copy_data <- gene_copy_data %>% dplyr::select(-t_gene)

# Filter out OR genes
cat("Filtering out OR genes...\n")
or_genes <- rownames(gene_copy_data)[grepl("^OR", rownames(gene_copy_data))]
cat("  Removing", length(or_genes), "OR genes\n")
gene_copy_data <- gene_copy_data[!grepl("^OR", rownames(gene_copy_data)), ]

or_genes <- rownames(gene_copy_data)[grepl("^ZNF", rownames(gene_copy_data))]
cat("  Removing", length(or_genes), "OR genes\n")
gene_copy_data <- gene_copy_data[!grepl("^ZNF", rownames(gene_copy_data)), ]

# Read metadata
metadata <- read.csv(file.path(ROOT, "data", "raxml_final_metadata_revised.csv"))

# Clean up gene copy data if needed
if("X" %in% colnames(gene_copy_data)) {
  gene_copy_data$X <- NULL
}

# Prepare metadata
metadata <- metadata %>%
  filter(!is.na(MLres)) %>%
  mutate(log_lifespan = log(maximum_longevity_y))

# Check if 'order' column exists, if not try 'Order'
if (!"order" %in% names(metadata)) {
  if ("Order" %in% names(metadata)) {
    metadata$order <- metadata$Order
  } else {
    cat("Warning: Neither 'order' nor 'Order' column found in metadata. Creating dummy order.\n")
    metadata$order <- "Unknown"
  }
}

# Specify genes you want to plot
#specific_genes <- head(pgls_results %>% arrange(pgls_results$adjusted_p_longevity), 10)$gene
specific_genes <- c("FBXO31")
#specific_genes <- common_high_genes[29:33]

# Alternative: Use top significant genes from PGLS results
# Get top 5 most significant genes
if (nrow(pgls_results) > 0) {
  # First, let's check what the first few rows look like
  cat("First few rows of pgls_results:\n")
  print(head(pgls_results, 3))
  cat("\n")
  
  # Try to identify the gene column - common names are: gene, Gene, gene_name, symbol, etc.
  possible_gene_columns <- c("gene", "Gene", "gene_name", "symbol", "t_symbol", "gene_symbol")
  gene_column <- intersect(possible_gene_columns, colnames(pgls_results))[1]
  
  if (is.na(gene_column)) {
    # If no standard gene column found, use the first column (often contains gene names)
    gene_column <- colnames(pgls_results)[1]
    cat("Warning: No standard gene column found. Using first column:", gene_column, "\n")
  } else {
    cat("Using gene column:", gene_column, "\n")
  }
  
  # Sort by p-value and get top genes
  top_genes <- pgls_results %>%
    filter(!is.na(adjusted_p_longevity)) %>%
    arrange(adjusted_p_longevity) %>%
    slice_head(n = 5) %>%
    pull(!!sym(gene_column))  # Use !!sym() to handle column name as string
  
  # Use top genes if specific genes aren't in the data
  genes_in_data <- intersect(specific_genes, rownames(gene_copy_data))
  if (length(genes_in_data) == 0) {
    cat("Specified genes not found. Using top significant genes instead.\n")
    specific_genes <- top_genes
  } else {
    specific_genes <- genes_in_data
  }
  
  cat("Analyzing genes:", paste(specific_genes, collapse = ", "), "\n")
}

# Function to create gene plot
create_gene_plot <- function(gene_name, n_outliers = 5) {
  # Check if gene exists in data
  if (!gene_name %in% rownames(gene_copy_data)) {
    cat("Gene", gene_name, "not found in copy number data\n")
    return(NULL)
  }
  
  # Get copy number data for this gene
  copy_data <- as.numeric(gene_copy_data[gene_name, ])
  species_names <- colnames(gene_copy_data)
  
  # Create data frame for visualization
  plot_df <- data.frame(
    species = species_names,
    copy_number = copy_data
  ) %>%
    filter(!is.na(copy_number)) %>%
    # Match with metadata
    inner_join(metadata, by = c("species" = "Scientific_name"))
  
  # Skip if no data
  if (nrow(plot_df) == 0) {
    cat("No data available for gene:", gene_name, "\n")
    return(NULL)
  }
  
  # Clean species names (remove underscores)
  plot_df$species_clean <- gsub("_", " ", plot_df$species)
  
  # Set up color mapping for orders
  unique_orders <- unique(plot_df$order)
  color_palette <- get_jewel_palette(length(unique_orders))
  names(color_palette) <- unique_orders
  
  # Identify top N outliers
  outlier_results <- identify_top_outliers(plot_df$copy_number, plot_df$log_lifespan, n_outliers)
  outlier_indices <- outlier_results$top_outliers
  
  # Get PGLS results for this gene
  gene_pgls_result <- pgls_results[pgls_results[[gene_column]] == gene_name, ]
  
  # Extract estimate_longevity and robust p-value if available
  if (nrow(gene_pgls_result) > 0) {
    # Safely extract estimate_longevity
    if ("estimate_longevity" %in% colnames(gene_pgls_result)) {
      estimate_longevity <- gene_pgls_result$estimate_longevity[1]
    } else {
      estimate_longevity <- NA
    }
    
    # Safely extract p-values
    adjusted_p_longevity_robust <- NULL
    if ("adjusted_p_longevity_robust" %in% colnames(gene_pgls_result)) {
      adjusted_p_longevity_robust <- gene_pgls_result$adjusted_p_longevity_robust[1]
    }
    
    adjusted_p_longevity_regular <- NULL
    if ("adjusted_p_longevity" %in% colnames(gene_pgls_result)) {
      adjusted_p_longevity_regular <- gene_pgls_result$adjusted_p_longevity[1]
    }
    
    # Use robust p-value if available and not NA, otherwise use regular p-value
    if (!is.null(adjusted_p_longevity_robust) && length(adjusted_p_longevity_robust) > 0 && !is.na(adjusted_p_longevity_robust)) {
      adjusted_p_longevity_to_show <- adjusted_p_longevity_robust
    } else if (!is.null(adjusted_p_longevity_regular) && length(adjusted_p_longevity_regular) > 0 && !is.na(adjusted_p_longevity_regular)) {
      adjusted_p_longevity_to_show <- adjusted_p_longevity_regular
    } else {
      adjusted_p_longevity_to_show <- NA
    }
  } else {
    # Fallback to correlation if PGLS results not found
    cor_test <- cor.test(plot_df$copy_number, plot_df$log_lifespan)
    estimate_longevity <- cor_test$estimate
    adjusted_p_longevity_to_show <- cor_test$p.value
  }
  
  # Create the plot
  p <- ggplot(plot_df, aes(x = copy_number, y = log_lifespan)) +
    # Points colored by order
    geom_point(aes(color = order, fill = order),
               size = 3, alpha = 0.7, shape = 21, stroke = 0.6) +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = color_palette) +
    
    # Add regression line
    geom_smooth(method = "lm", se = TRUE, color = "gray40", alpha = 0.3, linetype = 'dashed') +
    
    # Label only top N outliers (using cleaned species names)
    geom_text_repel(data = plot_df[outlier_indices, ],
                    aes(label = species_clean),
                    size = 4,
                    box.padding = 0.5,
                    point.padding = 0.2,
                    segment.color = "gray60",
                    segment.alpha = 0.7,
                    max.overlaps = Inf,
                    force = 2) +
    
    # Labels and theme
    labs(title = paste("Gene:", gene_name),
         subtitle = paste("coefficient =", 
                          ifelse(is.na(estimate_longevity), "NA", round(estimate_longevity, 3)),
                          "adjusted p value =", 
                          ifelse(is.na(adjusted_p_longevity_to_show), "NA", formatC(adjusted_p_longevity_to_show, format = "e", digits = 2))),
         x = "Copy number",
         y = "Maximum longevity (log)",
         color = "Order",
         fill = "Order") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 13),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  
  # Print outlier summary
  cat("\nTop", n_outliers, "outliers for", gene_name, ":\n")
  if (length(outlier_indices) > 0) {
    outlier_summary <- plot_df[outlier_indices, c("species_clean", "order")] %>%
      mutate(copy_number = plot_df$copy_number[outlier_indices],
             log_lifespan = plot_df$log_lifespan[outlier_indices]) %>%
      arrange(desc(outlier_results$composite_scores[outlier_indices]))
    colnames(outlier_summary)[1] <- "species"
    print(outlier_summary)
  }
  
  return(p)
}

# Create plots for specified genes
gene_plots <- list()
for (gene in specific_genes) {
  cat("Processing gene:", gene, "\n")
  plot <- create_gene_plot(gene, n_outliers = 5)
  if (!is.null(plot)) {
    gene_plots[[gene]] <- plot
  }
}

# Display plots
for (i in 1:length(gene_plots)) {
  print(gene_plots[[i]])
}

# Save individual plots (uncomment to save)
# for (gene in names(gene_plots)) {
#   filename <- paste0("plots/",
#                      gene, "_copy_vs_lifespan.pdf")
#   ggsave(filename, gene_plots[[gene]], width = 12, height = 8)
# }

# Create combined plot if you have multiple genes
if (length(gene_plots) > 1) {
  library(patchwork)
  combined_plot <- wrap_plots(gene_plots, ncol = 2) +
    plot_annotation(
      title = "Gene Copy Number vs Maximum Longevity",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
  
  print(combined_plot)
  
  # Save combined plot (uncomment to save)
  # ggsave("plots/genes_combined_plots.pdf",
  #        combined_plot, width = 16, height = 12)
}