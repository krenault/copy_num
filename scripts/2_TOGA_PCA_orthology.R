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
## Visualizing data using various methods

source(file.path(ROOT, "scripts", "FUN_color_mappings.R"))
color_mapping <- create_color_mapping(species_data$order)
########################################################################################################
# 1. PCA
########################################################################################################
library(tidyverse)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(ggrepel)
gene_copy_data <- read.csv(file.path(ROOT, "data", "All_Species_Orthologous_CopyNumber_Annotated.tsv"), row.names = "t_gene", sep = '\t')
rownames(gene_copy_data) <- gene_copy_data$t_symbol
gene_copy_data <- gene_copy_data %>% select(-t_symbol)

#available_species <- final_metadata$Scientific_name
#gene_copy_data <- gene_copy_data[, colnames(gene_copy_data) %in% available_species]
#write.csv(gene_copy_data, "file.path(ROOT, "data")/gene_copy_data_revised.csv")
#t2t_species <- read.csv('"path/to/local/file"', sep = '\t')
#t2t_species <- t2t_species %>% filter(contig.N50..bp. > 1000000)
#t2t_species$Species <- gsub(" ", "_", t2t_species$Species)
#t2t_species_available <- t2t_species$Species
#gene_copy_data <- gene_copy_data[, colnames(gene_copy_data) %in% t2t_species_available]

gene_copy_data <- t(gene_copy_data)  
final_metadata <- read.csv(file.path(ROOT, "data", "raxml_final_metadata_revised.csv"))
final_metadata <- final_metadata %>% filter(order == "Primates" | order == "Rodentia" | order == "Chiroptera" | order == "Artiodactyla" | order == "Carnivora"  )
final_metadata <- final_metadata %>% filter(order == "Primates" | order == "Rodentia" | order == "Chiroptera" | order == "Artiodactyla")

species_names <- rownames(gene_copy_data)
matched_data <- data.frame()
matched_orders <- c()

for (species in species_names) {
  species_space <- gsub("_", " ", species)
  idx <- which(final_metadata$Scientific_name == species | 
                 final_metadata$Scientific_name == species_space)
  
  if (length(idx) > 0) {
    matched_data <- rbind(matched_data, gene_copy_data[species, ])
    matched_orders <- c(matched_orders, final_metadata$order[idx[1]])
    rownames(matched_data)[nrow(matched_data)] <- species
  }
}

for (col in 1:ncol(matched_data)) {
  matched_data[is.na(matched_data[, col]), col] <- mean(matched_data[, col], na.rm = TRUE)
}
pca_result <- PCA(matched_data, graph = FALSE)
pca_data <- data.frame(
  PC1 = pca_result$ind$coord[, 1],
  PC2 = pca_result$ind$coord[, 2],
  PC3 = pca_result$ind$coord[, 3],
  PC4 = pca_result$ind$coord[, 4],
  Species = rownames(matched_data),
  Order = matched_orders
)
var_explained <- pca_result$eig[, "percentage of variance"]

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Order)) +
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


gene_contrib <- as.data.frame(pca_result$var$contrib[, 1:2])
gene_contrib$Gene <- colnames(matched_data)
top_genes_pc1 <- gene_contrib[order(-gene_contrib[, 1]), ][1:15, ]
top_genes_pc2 <- gene_contrib[order(-gene_contrib[, 2]), ][1:15, ]
top_genes <- rbind(
  data.frame(Gene = top_genes_pc1$Gene, 
             Contribution = top_genes_pc1[, 1], 
             PC = "PC1"),
  data.frame(Gene = top_genes_pc2$Gene, 
             Contribution = top_genes_pc2[, 2], 
             PC = "PC2")
)
top_genes$Direction <- ifelse(top_genes$PC == "PC1", 
                              "PC1 (species differentiation)", 
                              "PC2 (secondary variation)")
pc_colors <- c("PC1 (species differentiation)" = "#9b383a", 
               "PC2 (secondary variation)" = "#86acb9")
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
pca_3d_data <- data.frame(
  PC1 = pca_result$ind$coord[, 1],
  PC2 = pca_result$ind$coord[, 2],
  PC3 = pca_result$ind$coord[, 3],
  Species = rownames(matched_data),
  Order = matched_orders
)

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

########################################################################################################
# 2. UMAP
########################################################################################################

library(tidyverse)
library(umap)
library(ggrepel)

# Load and prepare data
gene_copy_data <- read.csv(file.path(ROOT, "data", "All_Species_Orthologous_CopyNumber_Annotated.tsv"),  sep = '\t')
gene_copy_data <- gene_copy_data %>% select(-t_symbol)
metadata <- read_csv(file.path(ROOT, "data", "raxml_final_metadata_revised.csv"))
#metadata <- metadata %>% filter(order == "Primates" | order == "Carnivora" | order == "Chiroptera" | order == "Artiodactyla")

gene_matrix <- gene_copy_data %>%
  `rownames<-`(NULL) %>%  
  column_to_rownames("t_gene") %>%
  as.matrix()
gene_matrix <- log1p(gene_matrix)
gene_matrix <- t(scale(t(gene_matrix)))

# Run UMAP
set.seed(123)
umap_result <- umap(t(gene_matrix))

# Prepare UMAP data with metadata
umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  Species = colnames(gene_matrix)
) %>%
  mutate(Species = str_replace(Species, "\\.", " ")) %>%
  left_join(metadata, by = c("Species" = "Scientific_name")) 

umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = order, label = Species)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = color_mapping) +
  labs(
    title = "UMAP of gene copy number",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Order"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

# geom_text_repel(data = umap_df %>% group_by(order) %>% slice(1), 
#                aes(label = order), size = 3, show.legend = FALSE)
print(umap_plot)
ggsave(file.path(ROOT, "plots", "umap_gene_copy_by_order.png"),
       umap_plot, width = 12, height = 8, dpi = 300)
####### Interactive ####### 
library(plotly)
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = order, 
                                 text = paste("Species:", Species, "<br>",
                                              "Order:", order, "<br>",
                                              "UMAP1:", round(UMAP1, 2), "<br>",
                                              "UMAP2:", round(UMAP2, 2)))) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = color_mapping) +
  labs(
    title = "UMAP of gene copy number",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Order"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

# Convert to interactive plot
interactive_plot <- ggplotly(umap_plot, tooltip = "text") %>%
  layout(
    hoverlabel = list(
      bgcolor = "white",
      font = list(size = 12)
    )
  )
interactive_plot

########################################################################################################
# 3. tSNE
########################################################################################################
library(Rtsne)
library(tidyverse)

# Load and prepare data
gene_copy_data <- read.csv(file.path(ROOT, "data", "All_Species_Orthologous_CopyNumber_Annotated.tsv"),  sep = '\t')
gene_copy_data <- gene_copy_data %>% select(-t_symbol)
metadata <- read_csv(file.path(ROOT, "data", "raxml_final_metadata_revised.csv"))
#metadata <- metadata %>% filter(order == "Primates" | order == "Carnivora" | order == "Chiroptera" | order == "Artiodactyla")

gene_matrix <- gene_copy_data %>%
  `rownames<-`(NULL) %>%  
  column_to_rownames("t_gene") %>%
  as.matrix()
gene_matrix <- log1p(gene_matrix)
gene_matrix <- t(scale(t(gene_matrix)))

set.seed(42)
tsne_result <- Rtsne(t(gene_matrix), perplexity = 30, dims = 2, check_duplicates = FALSE)

tsne_df <- data.frame(
  tSNE1 = tsne_result$Y[,1],
  tSNE2 = tsne_result$Y[,2],
  Species = colnames(gene_matrix)
) %>%
  mutate(Species = str_replace(Species, "\\.", " ")) %>%
  left_join(metadata, by = c("Species" = "Scientific_name")) %>%
  filter(!is.na(order)) # Remove species without order classification

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
