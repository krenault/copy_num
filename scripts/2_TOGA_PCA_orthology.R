
###############
##### PCA #####
###############
# 
# species_data <- read.csv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber.tsv", sep = '\t')
# gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/loss/Mammalia_all_species_copynumber.csv",
#                            row.names = "Gene")
# gene_copy_data <- gene_copy_data %>% dplyr::select(-X)
# 
# available_species <- species_data$Scientific_name
# gene_copy_data <- gene_copy_data[, colnames(gene_copy_data) %in% available_species]
# gene_copy_data_species <- colnames(gene_copy_data)
# all_final_data <- species_data[species_data$Scientific_name %in% gene_copy_data_species, ]

#############################
### Tree pruning section ###
#############################

# library(ape)
# original_tree <- read.tree("/Users/katiarenault/PhD/Databases/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")
# cat("\nOriginal tree information:\n")
# cat("Number of tips in original tree:", length(original_tree$tip.label), "\n")
# modified_tip_names <- sapply(strsplit(original_tree$tip.label, "_"), function(x) {
#   paste(head(x, -2), collapse = "_")
# })
# original_tree$tip.label <- modified_tip_names
# tree_species <- original_tree$tip.label
# data_species <- all_final_data$Scientific_name
# 
# common_species <- intersect(tree_species, data_species)
# cat("\nMatching species information:\n")
# cat("Number of species in original tree:", length(tree_species), "\n") # 4099
# cat("Number of species in data frame:", length(data_species), "\n") # 410
# cat("Number of species in common:", length(common_species), "\n") # 391
# 
# trimmed_tree <- drop.tip(original_tree, setdiff(tree_species, common_species))
# trimmed_data <- all_final_data[all_final_data$Scientific_name %in% common_species, ]
# 
# cat("\nAfter trimming:\n")
# cat("Number of species in pruned tree:", length(trimmed_tree$tip.label), "\n") # 366
# cat("Number of species in filtered data frame:", nrow(trimmed_data), "\n") # 366
# 
# final_tree_species <- trimmed_tree$tip.label
# final_metadata <- all_final_data %>%
#   filter(Scientific_name %in% final_tree_species)
# 
# # Verify the results
# cat("\nFinal matching verification:\n")
# cat("Number of species in pruned tree:", length(final_tree_species), "\n")
# cat("Number of species in final metadata:", nrow(final_metadata), "\n")
# write.tree(trimmed_tree, "/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_species_tree_revised.nwk")
# write.csv(final_metadata, "/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")

############
### PCA ####
############
gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber_Annotated.tsv", row.names = "t_gene", sep = '\t')
rownames(gene_copy_data) <- gene_copy_data$t_symbol
gene_copy_data <- gene_copy_data %>% select(-t_symbol)
gene_copy_data <- gene_copy_data[apply(gene_copy_data, 1, function(x) sum(x == 0, na.rm=TRUE)/sum(!is.na(x))) <= 0.5, ]

#gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/loss/Mammalia_all_species_partially_intact_copynumber.csv",
#                         row.names = "Gene")
#gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/pseudogenes/Mammalia_all_species_pseudo.csv", row.names = "Gene")
#gene_copy_data <- gene_copy_data %>% dplyr::select(-X)

###################################################################################
final_metadata <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")
# species_to_order <- setNames(final_metadata$order, final_metadata$Scientific_name)
# is_primate_specific <- function(gene_row, species_names, species_to_order) {
#   species_with_gene <- species_names[gene_row > 0]
#   orders_with_gene <- unique(species_to_order[species_with_gene])
#   orders_with_gene <- orders_with_gene[!is.na(orders_with_gene)]
#   if(length(orders_with_gene) == 1 && orders_with_gene == "Primates") {
#     return(TRUE)
#   } else {
#     return(FALSE)
#   }
# }
# species_names <- colnames(gene_copy_data)
# primate_specific_genes <- apply(gene_copy_data, 1, function(row) {
#   is_primate_specific(row, species_names, species_to_order)
# })
# 
# gene_copy_data <- gene_copy_data[!primate_specific_genes, ]
# 
# pan_column <- NULL
# possible_names <- c("Pan_troglodytes", "Pan troglodytes")
# for (name in possible_names) {
#   if (name %in% colnames(gene_copy_data)) {
#     pan_column <- name
#     cat("Found Pan troglodytes as:", pan_column, "\n")
#     break
#   }
# }
# if (is.null(pan_column)) {
#   cat("Warning: Pan_troglodytes not found in data columns!\n")
#   cat("Available species (first 10):", head(colnames(gene_copy_data), 10), "...\n")
#   possible_match <- colnames(gene_copy_data)[grep("Pan", colnames(gene_copy_data))]
#   if (length(possible_match) > 0) {
#     pan_column <- possible_match[1]
#     cat("Using closest match instead:", pan_column, "\n")
#   } else {
#     stop("Cannot proceed: Pan_troglodytes not found in dataset")
#   }
# }
# has_higher_copies_than_pan <- function(gene_row, pan_column) {
#   pan_copies <- gene_row[pan_column]
#   any_higher <- any(gene_row > pan_copies, na.rm = TRUE)
#   return(any_higher)
# }
# 
# higher_than_pan_genes <- apply(gene_copy_data, 1, function(row) {
#   has_higher_copies_than_pan(row, pan_column)
# })
# gene_copy_data <- gene_copy_data[higher_than_pan_genes, ]

###################################################################################

# gene_copy_data <- log(gene_copy_data + 0.01)
# gene_variances <- apply(gene_copy_data, 1, var, na.rm = TRUE)
# sorted_genes <- names(sort(gene_variances, decreasing = TRUE))
# num_genes_to_keep <- ceiling(length(sorted_genes) * 0.25)
# gene_copy_data <- gene_copy_data[sorted_genes[1:num_genes_to_keep], ]

available_species <- final_metadata$Scientific_name
gene_copy_data <- gene_copy_data[, colnames(gene_copy_data) %in% available_species]
write.csv(gene_copy_data, "/Users/katiarenault/PhD/TOGA/revised_results/data/gene_copy_data_revised.csv")
t2t_species <- read.csv('/Users/katiarenault/Downloads/overview.table (16).tsv', sep = '\t')
#t2t_species <- t2t_species %>% filter(contig.N50..bp. > 1000000)
t2t_species$Species <- gsub(" ", "_", t2t_species$Species)
t2t_species_available <- t2t_species$Species
gene_copy_data <- gene_copy_data[, colnames(gene_copy_data) %in% t2t_species_available]

# Simplified PCA of gene copy data with matching aesthetics
library(tidyverse)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(ggrepel)

################################################################################
sophisticated_jewels_palette <- c(
  "#9b383a", "#64a590", "#4878a0", "#AF6C36", "#945a87",
  "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
  "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
  "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
  "#c07a7a", "#8a7c6d", "#5c4d8b", "#4878a0", "#3D7B68"
)
################################################################################

# Load data (gene_copy_data is already filtered)
gene_copy_data <- t(gene_copy_data)  # Transpose to have species as rows
final_metadata <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")
final_metadata <- final_metadata %>%
  filter(Scientific_name %in% t2t_species_available)
final_metadata <- final_metadata %>% filter(order == "Primates" | order == "Rodentia" | order == "Chiroptera" | order == "Artiodactyla" | order == "Carnivora"  )
final_metadata <- final_metadata %>% filter(order == "Primates" | order == "Rodentia" | order == "Chiroptera" | order == "Artiodactyla")


# Match species in gene data to metadata
species_names <- rownames(gene_copy_data)
matched_data <- data.frame()
matched_orders <- c()

for (species in species_names) {
  # Try different formats for matching
  species_space <- gsub("_", " ", species)
  idx <- which(final_metadata$Scientific_name == species | 
                 final_metadata$Scientific_name == species_space)
  
  if (length(idx) > 0) {
    matched_data <- rbind(matched_data, gene_copy_data[species, ])
    matched_orders <- c(matched_orders, final_metadata$order[idx[1]])
    rownames(matched_data)[nrow(matched_data)] <- species
  }
}

# Clean data (replace NAs with column means)
for (col in 1:ncol(matched_data)) {
  matched_data[is.na(matched_data[, col]), col] <- mean(matched_data[, col], na.rm = TRUE)
}
#matched_data <- log1p(matched_data)
#matched_data <- t(scale(t(matched_data)))

# Run PCA
pca_result <- PCA(matched_data, graph = FALSE)

# Extract PCA results
pca_data <- data.frame(
  PC1 = pca_result$ind$coord[, 1],
  PC2 = pca_result$ind$coord[, 2],
  PC3 = pca_result$ind$coord[, 3],
  PC4 = pca_result$ind$coord[, 4],
  Species = rownames(matched_data),
  Order = matched_orders
)

# Get variance explained
var_explained <- pca_result$eig[, "percentage of variance"]

# Set up color mapping
orders <- unique(pca_data$Order)
if (length(orders) > length(sophisticated_jewels_palette)) {
  color_mapping <- setNames(
    sophisticated_jewels_palette[1:length(sophisticated_jewels_palette)], 
    orders[1:length(sophisticated_jewels_palette)]
  )
  warning("More orders than available colors. Some colors will be recycled.")
} else {
  color_mapping <- setNames(sophisticated_jewels_palette[1:length(orders)], orders)
}

# Create PCA plot with matching aesthetics
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Order)) +
  # Use shape = 21 for points with fill and color
  geom_point(size = 4, alpha = 0.7, stroke = 0.8, shape = 21, aes(fill = Order)) +
  scale_fill_manual(values = color_mapping) +
  scale_color_manual(values = color_mapping) +
  labs(
    title = "PCA of Gene Copy Number Variation",
    subtitle = paste0("PC1: ", round(var_explained[1], 1), "% variance, PC2: ", 
                      round(var_explained[2], 1), "% variance"),
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    guides(fill = guide_legend(title = "Order"),
           color = guide_legend(title = "Order")),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12)
  )

# ggsave("/Users/katiarenault/PhD/TOGA/results/preliminary_april/pca_gene_copy_by_order.png",
#        pca_plot, width = 12, height = 8, dpi = 300)

# Get top genes contributing to PC1 and PC2
gene_contrib <- as.data.frame(pca_result$var$contrib[, 1:2])
gene_contrib$Gene <- colnames(matched_data)
top_genes_pc1 <- gene_contrib[order(-gene_contrib[, 1]), ][1:15, ]
top_genes_pc2 <- gene_contrib[order(-gene_contrib[, 2]), ][1:15, ]

# Create combined dataframe for top genes
top_genes <- rbind(
  data.frame(Gene = top_genes_pc1$Gene, 
             Contribution = top_genes_pc1[, 1], 
             PC = "PC1"),
  data.frame(Gene = top_genes_pc2$Gene, 
             Contribution = top_genes_pc2[, 2], 
             PC = "PC2")
)

# Direction for plotting
top_genes$Direction <- ifelse(top_genes$PC == "PC1", 
                              "PC1 (species differentiation)", 
                              "PC2 (secondary variation)")

# Palette for PC1 and PC2
pc_colors <- c("PC1 (species differentiation)" = "#9b383a", 
               "PC2 (secondary variation)" = "#86acb9")

# Plot top genes contributions
genes_plot <- ggplot(top_genes, aes(x = reorder(Gene, Contribution), 
                                    y = Contribution, 
                                    fill = Direction)) +
  geom_col() +
  scale_fill_manual(values = pc_colors) +
  coord_flip() +
  facet_wrap(~PC, scales = "free_y") +
  labs(
    title = "Top genes contributing to principal components",
    x = "Gene",
    y = "Contribution (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

pca_plot

#########
library(plotly)
# 
# # Use your existing data
# # (assuming matched_data, pca_result, matched_orders, and color_mapping are already set up)
# 
# # Create 3D PCA plot data
pca_3d_data <- data.frame(
  PC1 = pca_result$ind$coord[, 1],
  PC2 = pca_result$ind$coord[, 2],
  PC3 = pca_result$ind$coord[, 3],
  Species = rownames(matched_data),
  Order = matched_orders
)

# Create interactive 3D PCA plot
# Create interactive 3D PCA plot
plot_3d <- plot_ly(pca_3d_data,
                   x = ~PC1,
                   y = ~PC2,
                   z = ~PC3,
                   color = ~Order,
                   colors = color_mapping,
                   text = ~Species,
                   hoverinfo = "text+x+y+z",
                   type = "scatter3d",
                   mode = "markers",
                   marker = list(size = 5, opacity = 0.8)) %>%
  plotly::layout(  # Explicitly use plotly's layout function
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(var_explained[1], 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(var_explained[2], 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(var_explained[3], 1), "%)")),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    ),
    title = "3D PCA of Gene Copy Number Variation"
  )

# Display the 3D plot
plot_3d

###########################################
# UMAP

library(tidyverse)
library(umap)
library(ggrepel)

# Load and prepare data
copy_num_data <- read_tsv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber_Annotated.tsv")
metadata <- read_csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")
#metadata <- metadata %>% filter(order == "Primates" | order == "Carnivora" | order == "Chiroptera" | order == "Artiodactyla")

# Prepare gene matrix
gene_matrix <- copy_num_data %>%
  select(-t_symbol) %>%
  column_to_rownames("t_gene") %>%
  as.matrix()
# Log-transform and scale
gene_matrix <- log1p(gene_matrix)
gene_matrix <- t(scale(t(gene_matrix)))

# Define your color palette
sophisticated_jewels_palette <- c(
  "#9b383a", "#64a590", "#4878a0", "#AF6C36", "#945a87",
  "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
  "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
  "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
  "#c07a7a", "#8a7c6d", "#5c4d8b", "#4878a0", "#3D7B68"
)

# Run UMAP
#set.seed(42) # For reproducibility
set.seed(123)
umap_result <- umap(t(gene_matrix))

# Prepare UMAP data with metadata
umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  Species = colnames(gene_matrix)
) %>%
  mutate(Species = str_replace(Species, "\\.", " ")) %>%
  left_join(metadata, by = c("Species" = "Scientific_name")) %>%
  filter(!is.na(order)) # Remove species without order classification

# Get unique orders and assign colors
unique_orders <- unique(umap_df$order)
color_mapping <- setNames(
  sophisticated_jewels_palette[1:length(unique_orders)],
  unique_orders
)

# Create the plot
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = order, label = Species)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = color_mapping) +
  labs(
    title = "UMAP Projection of Gene Copy Number Variation",
    subtitle = "Colored by Taxonomic Order",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Taxonomic Order"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# For better visualization of dense clusters, you can add:
# geom_text_repel(data = umap_df %>% group_by(order) %>% slice(1), 
#                aes(label = order), size = 3, show.legend = FALSE)

# Display the plot
print(umap_plot)

#################################################
# Load required libraries
library(Rtsne)
library(tidyverse)

# Load and prepare data
copy_num_data <- read_tsv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber_Annotated.tsv")
metadata <- read_csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")
#metadata <- metadata %>% filter(order == "Primates" | order == "Carnivora" | order == "Chiroptera" | order == "Artiodactyla")

# Prepare gene matrix
gene_matrix <- copy_num_data %>%
  select(-t_symbol) %>%
  column_to_rownames("t_gene") %>%
  as.matrix()
gene_matrix <- log1p(gene_matrix)
gene_matrix <- t(scale(t(gene_matrix)))

# Define your color palette
sophisticated_jewels_palette <- c(
  "#9b383a", "#64a590", "#4878a0", "#AF6C36", "#945a87",
  "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
  "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
  "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
  "#c07a7a", "#8a7c6d", "#5c4d8b", "#4878a0", "#3D7B68"
)

# Run t-SNE
set.seed(42) # For reproducibility
tsne_result <- Rtsne(t(gene_matrix), perplexity = 30, dims = 2, check_duplicates = FALSE)

# Prepare t-SNE data with metadata
tsne_df <- data.frame(
  tSNE1 = tsne_result$Y[,1],
  tSNE2 = tsne_result$Y[,2],
  Species = colnames(gene_matrix)
) %>%
  mutate(Species = str_replace(Species, "\\.", " ")) %>%
  left_join(metadata, by = c("Species" = "Scientific_name")) %>%
  filter(!is.na(order)) # Remove species without order classification

# Get unique orders and assign colors
unique_orders <- unique(tsne_df$order)
color_mapping <- setNames(
  sophisticated_jewels_palette[1:length(unique_orders)],
  unique_orders
)

# Create the t-SNE plot
ggplot(tsne_df, aes(tSNE1, tSNE2, color = order)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = color_mapping) +
  labs(title = "t-SNE of All Gene Copy Numbers",
       x = "t-SNE 1",
       y = "t-SNE 2",
       color = "Order") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
