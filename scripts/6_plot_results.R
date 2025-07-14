## Katia Renault
## Visualizing results from copy number PGLMM analysis

source("/Users/katiarenault/Documents/Github/copy_num/scripts/FUN_color_mappings.R")
color_mapping <- create_color_mapping(species_data$order)
########################################################################################################
# 1. Gene result plot
########################################################################################################

library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

#volcano_data <- read.csv('/Users/katiarenault/Documents/GitHub/copy_num/results/max_longevity_zero_filtered_poisson_pglmm_20250702_110253.csv')
volcano_data <- read.csv("/Users/katiarenault/Documents/GitHub/copy_num/results/mlres_zero_filtered_poisson_pglmm_20250707_075917.csv")

volcano_data$p_value <- as.numeric(volcano_data$p_value_longevity)
volcano_data$FDR <- volcano_data$adjusted_p_longevity
volcano_data$estimate <- as.numeric(volcano_data$estimate_longevity)

volcano_data <- volcano_data %>%
  # drop_na() %>%
  mutate(
    significance_level = case_when(
      FDR < 0.05 | abs(estimate) >= 2.2 ~ "FDR_significant",
      p_value < 0.05 ~ "p_significant", 
      TRUE ~ "not_significant"
    ),
    point_color = case_when(
      significance_level == "FDR_significant" ~ estimate,
      significance_level == "p_significant" ~ estimate * 0.5,  # Lighter version
      TRUE ~ NA_real_
    ),
    log_p = -log(p_value)
  ) %>%
  mutate(significance_rank = rank(-(abs(estimate) * log_p))) %>%
  mutate(to_label = ifelse(significance_level == "FDR_significant" & significance_rank <= 15, gene, NA))

volcano_plot <- ggplot(volcano_data, aes(x = estimate, y = log_p)) +
  geom_point(data = filter(volcano_data, significance_level == "not_significant"),
             size = 4, alpha = 0.3, shape = 21, 
             fill = "grey80", color = "grey60") +
  geom_point(data = filter(volcano_data, significance_level == "p_significant"),
             aes(fill = point_color),
             size = 4, alpha = 0.5, shape = 21, 
             color = "grey70", stroke = 0.8) +
  geom_point(data = filter(volcano_data, significance_level == "FDR_significant"),
             aes(fill = point_color, 
                 color = after_scale(scales::alpha(fill, 0.4))),
             size = 4, alpha = 0.8, shape = 21, stroke = 1.2) +
  geom_text_repel(aes(label = to_label),
                  size = 5,
                  box.padding = 0.5,
                  point.padding = 0.2,
                  min.segment.length = 0.2,
                  segment.color = "grey50",
                  segment.alpha = 0.7,
                  max.overlaps = Inf) +
  scale_fill_gradient2(low = "#86acb9", mid = "white", high = "#9b383a", 
                       midpoint = 0, name = "Estimated coefficient",
                       guide = guide_colorbar(barheight = unit(5, "cm"))) +
  scale_color_identity() +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "gray40", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", alpha = 0.7) +
  labs(title = "Volcano Plot of GO Pathway PGLS Results",
       x = "Estimated coefficient",
       y = "-log(p-value)") +
  
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12)
  )

print(volcano_plot)
table(volcano_data$significance_level, useNA = "ifany")
ggsave("/Users/katiarenault/Documents/GitHub/copy_num/plots/mlres_volcano.png",
       volcano_plot, width = 12, height = 8, dpi = 300)


########################################################################################################
# 2. Pathways result plot
########################################################################################################

cors_pathways <- read.csv("/Users/katiarenault/Documents/GitHub/copy_num/results/mlres_zero_filtered_poisson_pglmm_20250707_075917.csv")
#cors_pathways <- read.csv("/Users/katiarenault/Documents/GitHub/copy_num/results/max_longevity_zero_filtered_poisson_pglmm_20250702_110253.csv")
pathways.reactome <- gmtPathways("/Users/katiarenault/Desktop/PhD/Human_GSEA/h.all.v2023.2.Hs.symbols.gmt")
pathways.hallmark <- gmtPathways("/Users/katiarenault/Desktop/PhD/Human_GSEA/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
pathways.kegg <- gmtPathways("/Users/katiarenault/Desktop/PhD/Human_GSEA/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
#pathways.go <- gmtPathways("/Nori_1/krenault/hibernator_rer_converge/data/c5.go.v2023.2.Hs.symbols.gmt")
pathways.hallmark <- c(pathways.kegg, pathways.hallmark, pathways.reactome)

cors_pathways$estimate_longevity <- as.numeric(cors_pathways$estimate_longevity)
cors_pathways <- cors_pathways %>% tidyr::drop_na()
positive_significant_genes <- cors_pathways %>% 
  # filter(Rho < 0) %>%
  as_tibble() %>%
  arrange(p_value_longevity)

positive_stats_vector <- setNames(positive_significant_genes$estimate_longevity, positive_significant_genes$gene)
positive_fgsea_results <- fgsea(pathways = pathways.hallmark, stats = positive_stats_vector)
positive_fgsea_results <- data.frame(positive_fgsea_results)
positive_fgsea_results <- positive_fgsea_results %>% select(-leadingEdge)
#write.csv(positive_fgsea_results, "/Users/katiarenault/Documents/GitHub/copy_num/results/mlres_zero_filtered_poisson_pglmm_20250707_075917_pathways.csv")

volcano_data <- read.csv("/Users/katiarenault/Documents/GitHub/copy_num/results/mlres_zero_filtered_poisson_pglmm_20250707_075917_pathways.csv")
#volcano_data <- read.csv("/Users/katiarenault/Documents/GitHub/copy_num/results/max_longevity_zero_filtered_poisson_pglmm_20250702_110253_pathways.csv")
volcano_data <- volcano_data %>% mutate(
  pathway = stringr::str_replace_all(pathway, "_", " "))
volcano_data$p_value <- as.numeric(volcano_data$pval)
volcano_data$FDR <- volcano_data$padj
volcano_data$estimate <- as.numeric(volcano_data$NES)

volcano_data <- volcano_data %>%
  mutate(
    significance_level = case_when(
      FDR < 0.05 | abs(estimate) >= 2.2 ~ "FDR_significant",
      p_value < 0.05 ~ "p_significant", 
      TRUE ~ "not_significant"
    ),
    point_color = case_when(
      significance_level == "FDR_significant" ~ estimate,
      significance_level == "p_significant" ~ estimate * 0.5,  # Lighter version
      TRUE ~ NA_real_
    ),
    log_p = -log(p_value),
    log_p_transformed = case_when(
      log_p <= 10 ~ log_p,                    # Keep 0-10 as is
      log_p > 10 ~ 10 + (log_p - 10) * 0.3   # Compress 10+ range
    ),
    direction = case_when(
      estimate > 0 ~ "positive",
      estimate < 0 ~ "negative",
      TRUE ~ "neutral"
    )
  ) %>%
  group_by(direction) %>%
  mutate(
    direction_significance_rank = rank(-(abs(estimate) * log_p), ties.method = "first")
  ) %>%
  ungroup() %>%
  mutate(
    to_label = case_when(
      direction == "positive" & 
        significance_level %in% c("FDR_significant", "p_significant") & 
        direction_significance_rank <= 5 ~ pathway,
      direction == "negative" & 
        significance_level %in% c("FDR_significant", "p_significant") & 
        direction_significance_rank <= 5 ~ pathway,
      TRUE ~ NA_character_
    )
  )

volcano_plot <- ggplot(volcano_data, aes(x = estimate, y = log_p_transformed)) +
  geom_point(data = filter(volcano_data, significance_level == "not_significant"),
             size = 4, alpha = 0.3, shape = 21, 
             fill = "grey80", color = "grey60") +
  geom_point(data = filter(volcano_data, significance_level == "p_significant"),
             aes(fill = point_color),
             size = 4, alpha = 0.5, shape = 21, 
             color = "grey70", stroke = 0.8) +
  geom_point(data = filter(volcano_data, significance_level == "FDR_significant"),
             aes(fill = point_color, 
                 color = after_scale(scales::alpha(fill, 0.4))),
             size = 4, alpha = 0.8, shape = 21, stroke = 1.2) +
  geom_text_repel(aes(label = to_label),
                  size = 3,
                  box.padding = 0.5,
                  point.padding = 0.2,
                  min.segment.length = 0.2,
                  segment.color = "grey50",
                  segment.alpha = 0.7,
                  max.overlaps = 12, 
                  force = 2,  
                  seed = 123) +  
  scale_fill_gradient2(low = "#86acb9", mid = "white", high = "#9b383a", 
                       midpoint = 0, name = "Estimated coefficient",
                       guide = guide_colorbar(barheight = unit(5, "cm"))) +
  scale_color_identity() +
  scale_y_continuous(
    breaks = c(0, 2, 4, 6, 8, 10, 11, 12, 13),
    labels = c("0", "2", "4", "6", "8", "10", "15", "20", "25"),
    minor_breaks = seq(0, 10, by = 1),
    expand = c(0.02, 0.02)
  ) +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "gray40", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", alpha = 0.7) +
  labs(title = "Volcano Plot of GO Pathway PGLS Results",
       subtitle = "Y-axis: 0-10 shown in detail, 10+ compressed (top 5 positive & negative pathways labeled)",
       x = "Estimated coefficient",
       y = "-log(p-value)") +
  
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_line(color = "grey95", size = 0.3),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12)
  )

print(volcano_plot)
cat("Summary of labeled pathways:\n")
labeled_summary <- volcano_data %>% 
  filter(!is.na(to_label)) %>% 
  group_by(direction) %>%
  summarise(
    count = n(),
    .groups = "drop"
  )
print(labeled_summary)
cat("\nDetailed labeled pathways:\n")
labeled_pathways <- volcano_data %>% 
  filter(!is.na(to_label)) %>% 
  arrange(direction, direction_significance_rank) %>%
  select(pathway, direction, significance_level, p_value, FDR, estimate, direction_significance_rank)
print(labeled_pathways)

cat("\nSummary of significance levels:\n")
table(volcano_data$significance_level, useNA = "ifany")
ggsave("/Users/katiarenault/Documents/GitHub/copy_num/plots/mlres_pathway_volcano.png",
       volcano_plot, width = 12, height = 8, dpi = 300)


###########
## FGSEA ##
###########
# 
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(scales)
# library(stringr)
# 
# go_data_file <- read.csv("/Nori_1/krenault/copy_num/results/phyr_pglmm_zero_filtered_poisson/pathway_results/zero_filtered_poisson_pglmm_20250707_075917_pathways.csv")
# 
# plot_go_terms <- function(go_data_file, keywords = NULL, n_pathways = 15) {
#   go_data <- as.data.frame(go_data_file)
#   go_data$pathway <- gsub("_", " ", go_data$pathway)
#   
#   # Apply keyword filtering if provided
#   if (!is.null(keywords) && length(keywords) > 0) {
#     # Create a regex pattern to match any of the keywords (case insensitive)
#     pattern <- paste(keywords, collapse = "|")
#     go_data <- go_data %>%
#       filter(str_detect(tolower(pathway), tolower(pattern)))
#     
#     # If no pathways match the keywords, warn the user
#     if (nrow(go_data) == 0) {
#       warning("No pathways match the provided keywords. Showing all pathways instead.")
#       go_data <- as.data.frame(go_data_file)
#       go_data$pathway <- gsub("_", " ", go_data$pathway)
#     }
#   }
#   
#   # Prepare data for plotting
#   plot_data <- go_data %>%
#     arrange(padj) %>%  
#     head(n_pathways) %>%  # Use the parameter for number of pathways      
#     mutate(
#       pathway = factor(pathway, levels = pathway),
#       log_p_adj = -log(padj)
#     )
#   
#   # Create the plot
#   p <- ggplot(plot_data, aes(x = NES, y = reorder(pathway, NES))) +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", alpha = 0.5) +
#     geom_point(aes(size = log_p_adj, color = NES), alpha = 0.9) +
#     scale_color_gradient2(
#       low = "#86acb9",   
#       mid = "#ffffff",  
#       high = "#9b383a",  
#       midpoint = 0
#     ) +
#     scale_size_continuous(range = c(3, 10), name = "-log(adj.P)") +
#     labs(
#       title = "Pathway enrichment analysis",
#       x = "NES",
#       y = "Term",
#       color = "NES"
#     ) +
#     theme_minimal() +
#     theme(
#       plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
#       axis.text.y = element_text(size = 14),
#       axis.title = element_text(size = 14),
#       legend.position = "right",
#       panel.grid.major.y = element_line(color = "gray90"),
#       panel.grid.minor = element_blank(),
#       plot.margin = margin(1, 1, 1, 1, "cm")
#     )
#   
#   # Add subtitle if keywords were used
#   if (!is.null(keywords) && length(keywords) > 0) {
#     p <- p + labs(subtitle = paste0("Filtered by keywords: ", paste(keywords, collapse = ", ")))
#   }
#   
#   return(p)
# }
# 
# # Example usage:
# # Default (no keywords)
# p1 <- plot_go_terms(go_data_file)
# 
# # With keywords (e.g., to focus on pathways related to DNA repair)
# #p2 <- plot_go_terms(go_data_file, keywords = c("dna repair", "mitochondrial translation", "parkinsons", "extracellular matrix organization", "alzheimers", "Kegg spliceosome", "peroxisome", "Reactome translation", "ubiquitin-mediated proteolysis", "fatty acid metabolism", "citrate cycle", "growth hormone", "base excision repair", "kegg ribosome", "oxidative phosphorylation", "spliceosome", "telomere maintenance", "insulin processing", "respiratory electron transport"))
# p2 <- plot_go_terms(go_data_file, keywords = c("Hallmark dna repair", "citric acid", "growth hormone", "base excision repair", "telomere maintenance", "insulin processing", "dna repair"))
# p2 <- plot_go_terms(go_data_file, 
#                     keywords = c("growth hormone", 
#                                  "reactome base excision repair", "rna polymerase ii transcription",
#                                  "ubiquitin specific processing proteases", 
#                                  "deubiquitination", "neuronal system", "hcmv late", "ub specific",
#                                  "orc complex assembly", "sirt1", "reactome dna methylation", "b wich"),
#                     n_pathways = 20)
########################################################################################################
# 3. Top and bottom gene contribution
########################################################################################################

library(dplyr)
library(ggplot2)
library(stringr)
plot_gene_results_top_bottom <- function(gene_data_file, keywords = NULL, n_genes = 20, 
                                         metric = "estimate_longevity", p_threshold = 0.05) {
  gene_data <- as.data.frame(gene_data_file)
  if ("gene" %in% names(gene_data)) {
    gene_data$gene <- gsub("_", " ", gene_data$gene)
    gene_column <- "gene"
  } else if ("Gene" %in% names(gene_data)) {
    gene_data$Gene <- gsub("_", " ", gene_data$Gene)
    gene_column <- "Gene"
  } else {
    # Try to find a likely gene column
    possible_gene_cols <- c("gene_name", "symbol", "gene_symbol", "t_symbol")
    gene_column <- intersect(possible_gene_cols, names(gene_data))[1]
    if (is.na(gene_column)) {
      gene_column <- names(gene_data)[1]  # Use first column as fallback
      warning("Gene column not clearly identified. Using first column: ", gene_column)
    }
    gene_data[[gene_column]] <- gsub("_", " ", gene_data[[gene_column]])
  }
  if (!is.null(keywords) && length(keywords) > 0) {
    pattern <- paste(keywords, collapse = "|")
    gene_data <- gene_data %>%
      filter(str_detect(tolower(.data[[gene_column]]), tolower(pattern)))
    if (nrow(gene_data) == 0) {
      warning("No genes match the provided keywords. Showing all genes instead.")
      gene_data <- as.data.frame(gene_data_file)
      if (gene_column %in% names(gene_data)) {
        gene_data[[gene_column]] <- gsub("_", " ", gene_data[[gene_column]])
      }
    }
  }
  p_value_cols <- c("p_adj", "padj", "FDR", "p_value_longevity", "p_value", "pval", "P.Value")
  available_p_cols <- intersect(p_value_cols, names(gene_data))
  
  if (length(available_p_cols) > 0) {
    p_col <- available_p_cols[1]  # Use the first available p-value column (prioritizes adjusted)
    gene_data <- gene_data %>%
      filter(.data[[p_col]] <= p_threshold)
    if (nrow(gene_data) == 0) {
      warning(paste("No genes pass the p-value threshold of", p_threshold, 
                    ". Showing all genes instead."))
      gene_data <- as.data.frame(gene_data_file)
      if (gene_column %in% names(gene_data)) {
        gene_data[[gene_column]] <- gsub("_", " ", gene_data[[gene_column]])
      }
    }
  } else {
    p_col <- NULL
    warning("No p-value column found. Showing results without p-value filtering.")
  }
  if (!metric %in% names(gene_data)) {
    stop(paste("Column '", metric, "' not found in the data. Available columns: ", 
               paste(names(gene_data), collapse = ", ")))
  }
  gene_data <- gene_data %>%
    filter(!is.na(.data[[metric]]))
  top_genes <- gene_data %>%
    filter(.data[[metric]] > 0) %>%  # Only positive values
    arrange(.data[[if (!is.null(p_col)) p_col else metric]], desc(.data[[metric]])) %>%
    head(n_genes) %>%
    mutate(category = "Top genes")
  bottom_genes <- gene_data %>%
    filter(.data[[metric]] < 0) %>%  # Only negative values
    arrange(.data[[if (!is.null(p_col)) p_col else metric]], .data[[metric]]) %>%
    head(n_genes) %>%
    mutate(category = "Bottom genes")
  plot_data <- bind_rows(top_genes, bottom_genes)
  if (nrow(plot_data) == 0) {
    stop("No genes found for plotting. Check your data and filtering criteria.")
  }
  positive_data <- plot_data %>%
    filter(category == "Top genes") %>%
    arrange(if (!is.null(p_col)) .data[[p_col]] else .data[[metric]])
  negative_data <- plot_data %>%
    filter(category == "Bottom genes") %>%
    arrange(if (!is.null(p_col)) desc(.data[[p_col]]) else desc(.data[[metric]]))
  plot_data <- bind_rows(positive_data, negative_data) %>%
    mutate(
      gene_factor = factor(.data[[gene_column]], levels = .data[[gene_column]]),
      metric_value = .data[[metric]],
      direction = ifelse(metric_value > 0, "Positive", "Negative")
    )
  if (!is.null(p_col)) {
    plot_data <- plot_data %>%
      mutate(
        log_p = -log(.data[[p_col]]),
        p_intensity = (log_p - min(log_p, na.rm = TRUE)) / (max(log_p, na.rm = TRUE) - min(log_p, na.rm = TRUE))
      )
    plot_data$p_intensity <- pmax(plot_data$p_intensity, 0.3)
  } else {
    plot_data$p_intensity <- 1
  }
  direction_colors <- c("Positive" = "#9b383a", "Negative" = "#86acb9")
  p <- ggplot(plot_data, aes(x = gene_factor, y = metric_value, fill = direction, alpha = p_intensity)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
    geom_text(aes(label = .data[[gene_column]], 
                  y = ifelse(metric_value > 0, metric_value * 0.5, metric_value * 0.5)),
              hjust = 0.5, size = 3, color = "black", fontface = "bold") +
    scale_fill_manual(values = direction_colors, name = "Effect direction") +
    scale_alpha_identity() +  # Use the actual alpha values we calculated
    coord_flip() +
    labs(
     # title = paste("Top", n_genes, "and Bottom", n_genes, "Genes by", metric),
      x = "",  # Remove x-axis label since gene names are on bars
      y = "Estimated coefficient",
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.y = element_blank(),  # Remove y-axis text since labels are on bars
      axis.ticks.y = element_blank(), # Remove y-axis ticks
      axis.text.x = element_text("Estimated coefficient", size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  return(p)
}
gene_data_file <- read.csv("/Users/katiarenault/Documents/GitHub/copy_num/results/max_longevity_zero_filtered_poisson_pglmm_20250702_110253.csv")
p_combined <- plot_gene_results_top_bottom(gene_data_file, n_genes = 20)
print(p_combined)
ggsave("/Users/katiarenault/Documents/GitHub/copy_num/plots/max_longevity_top_genes.png",
       p_combined, width = 12, height = 8, dpi = 300)
