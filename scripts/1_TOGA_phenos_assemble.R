###################################
## Phenos assembly for the TOGA ##
##################################
###############
### Amniote ###
###############


# amniote_df <- read.delim("/Users/katiarenault/PhD/Databases/Amniote-Life-History-Database.tsv")
# amniote_df$scientific_name <- paste(amniote_df$genus, amniote_df$species)
# amniote_df$scientific_name <- gsub(" ", "_", amniote_df$scientific_name)

library(dplyr)
copy_number_animal_data <- read.csv("/Users/katiarenault/PhD/TOGA/unique_species_list.csv", header = TRUE)
# amniote_df <- amniote_df %>%
#   dplyr::select(scientific_name, common_name, maximum_longevity_y, adult_body_mass_g, order)
# columns_to_process <- c("maximum_longevity_y", "adult_body_mass_g")
# 
# for (col in columns_to_process) {
#   # Make sure the column exists
#   if (col %in% names(amniote_df)) {
#     # Convert column to character if it isn't already
#     amniote_df[[col]] <- as.character(amniote_df[[col]])
#     
#     # Replace "Missing" with NA
#     amniote_df[[col]][amniote_df[[col]] == "Missing"] <- NA
#     
#     # Process remaining values to extract number between "[" and ","
#     for (i in 1:nrow(amniote_df)) {
#       if (!is.na(amniote_df[[col]][i])) {
#         # Check if the value contains both "[" and ","
#         if (grepl("\\[.*,", amniote_df[[col]][i])) {
#           # Extract the number between "[" and ","
#           extracted_number <- gsub(".*\\[(.*),.*", "\\1", amniote_df[[col]][i])
#           amniote_df[[col]][i] <- extracted_number
#         }
#       }
#     }
#     
#     # Convert the column back to numeric
#     amniote_df[[col]] <- as.numeric(amniote_df[[col]])
#   } else {
#     warning(paste("Column", col, "not found in the dataframe"))
#   }
# }
# write.csv(amniote_df, "/Users/katiarenault/PhD/Databases/minimal_amniote_df.csv")

amniote_df <- read.csv("/Users/katiarenault/PhD/Databases/minimal_amniote_df.csv")
merged_data <- copy_number_animal_data %>%
  left_join(amniote_df %>% dplyr::select(scientific_name, maximum_longevity_y, adult_body_mass_g, order, common_name), 
            by = c("Scientific_name" = "scientific_name"))

#############
### AnAge ###
#############
anage_df <- read.delim("/Users/katiarenault/PhD/Databases/anage_data.txt")
anage_df <- anage_df %>% 
  mutate(Genus_species = paste(Genus, Species, sep = "_"))
merged_data_anage <- copy_number_animal_data %>%
  left_join(anage_df %>% dplyr::select(`Genus_species`, Maximum.longevity..yrs., Adult.weight..g., Order, Common.name), 
            by = c("Scientific_name" = "Genus_species"))

merged_data$maximum_longevity_y[merged_data$maximum_longevity_y < 0] <- NA
print("Number of species with known age in Amniote:") 
print(sum(!is.na(merged_data$maximum_longevity_y))) ### 404
print("Number of species with unknown age in Amniote:")
print(sum(is.na(merged_data$maximum_longevity_y))) ### 146

merged_data_anage$Maximum.longevity..yrs.[merged_data_anage$Maximum.longevity..yrs. < 0] <- NA
print("Number of species with known age in Amniote:") 
print(sum(!is.na(merged_data_anage$Maximum.longevity..yrs.))) ### 325
print("Number of species with unknown age in Amniote:")
print(sum(is.na(merged_data_anage$Maximum.longevity..yrs.))) ### 225

# Merge the two data frames
final_merged_data <- merged_data %>%
  full_join(merged_data_anage, by = "Scientific_name") %>%
  mutate(
    # Default to anage values if available, otherwise keep amniote values
    maximum_longevity_y = coalesce(Maximum.longevity..yrs., maximum_longevity_y),
    adult_body_mass_g = coalesce(Adult.weight..g., adult_body_mass_g),
    order = coalesce(Order, order),
    common_name = coalesce(Common.name, common_name)
  ) %>%
  # Drop redundant columns from the second data frame
  dplyr::select(-Maximum.longevity..yrs., -Adult.weight..g., -Order, -Common_name.x, -Common_name.y, -Common.name)

# Count the number of updated values
updated_longevity <- sum(!is.na(merged_data_anage$Maximum.longevity..yrs.) & is.na(merged_data$maximum_longevity_y))
updated_mass <- sum(!is.na(merged_data_anage$Adult.weight..g.) & is.na(merged_data$adult_body_mass_g))

# Print the counts
cat("Number of maximum longevity values updated from anage:", updated_longevity, "\n") # 11 
cat("Number of adult body mass values updated from anage:", updated_mass, "\n") # 29

unmatched_data <- final_merged_data %>%
  filter(is.na(maximum_longevity_y)) %>%
  dplyr::select(-maximum_longevity_y)

new_matches <- unmatched_data %>%
  left_join(
    amniote_df %>%
      dplyr::select(common_name, maximum_longevity_y) %>%
      mutate(common_name = tolower(common_name)),  # Convert to lowercase in txt_data
    by = c("common_name" = "common_name")  
  ) %>%
  mutate(common_name  = tolower(common_name))

new_matches$common_name <- gsub(" ", "_", new_matches$common_name)
new_matches_unique <- new_matches %>%
  distinct(common_name, .keep_all = TRUE)

final_data <- final_merged_data %>%
  rows_update(new_matches_unique, by = c("Scientific_name", "common_name"), unmatched = "ignore")
final_data$maximum_longevity_y[final_data$maximum_longevity_y < 0] <- NA
print("Number of species with known age:")
print(sum(!is.na(final_data$maximum_longevity_y))) ### 415
print("Number of species with unknown age:")
print(sum(is.na(final_data$maximum_longevity_y))) ### 135
print("Number of species with unknown body mass:")
final_data$adult_body_mass_g[final_data$adult_body_mass_g < 0] <- NA
print(sum(is.na(final_data$adult_body_mass_g))) ### 108
analysis_data <- final_data[!is.na(final_data$maximum_longevity_y) & 
                              !is.na(final_data$adult_body_mass_g) &
                              final_data$adult_body_mass_g > 0 &
                              final_data$maximum_longevity_y > 0, ]
analysis_data$log_longevity <- log10(analysis_data$maximum_longevity_y)
analysis_data$log_mass <- log10(analysis_data$adult_body_mass_g)
lm_model <- lm(log_longevity ~ log_mass, data = analysis_data)
analysis_data$MLres <- residuals(lm_model)
final_data$MLres <- NA
final_data$MLres[match(rownames(analysis_data), rownames(final_data))] <- analysis_data$MLres

unique_species_data <- final_data %>%
  filter(!is.na(maximum_longevity_y)) %>%
  group_by(Scientific_name) %>%
  slice_tail(n = 1) %>%
  ungroup()
all_final_data <- final_data %>% filter(!is.na(maximum_longevity_y)) %>% filter(!is.na(adult_body_mass_g))

write.csv(unique_species_data, "/Users/katiarenault/PhD/TOGA/all_max_lifespan_unique_mammalian_species_metadata.csv")
write.csv(all_final_data, "/Users/katiarenault/PhD/TOGA/all_max_mass_lifespan_mammalian_species_metadata.csv")

### 410 unique species at the end with both lifespan and body mass info
### 415 unique species at the end with lifespan info

# library(ape)  # For tree manipulation
# tree <- read.tree("/Users/katiarenault/PhD/TOGA/unique_species_list.nwk")
# tree_species <- tree$tip.label
# data_species <- all_final_data$Scientific_name
# common_species <- intersect(tree_species, data_species)
# 
# cat("Number of species in tree:", length(tree_species), "\n")
# cat("Number of species in data frame:", length(data_species), "\n")
# cat("Number of species in common:", length(common_species), "\n")
# trimmed_tree <- drop.tip(tree, setdiff(tree_species, common_species))
# trimmed_data <- all_final_data[all_final_data$Scientific_name %in% common_species, ]
# cat("\nAfter trimming:\n")
# cat("Number of species in trimmed tree:", length(trimmed_tree$tip.label), "\n")
# cat("Number of species in trimmed data frame:", nrow(trimmed_data), "\n")
# 
# write.tree(trimmed_tree, "/Users/katiarenault/PhD/TOGA/final_species_tree.nwk")
# write.csv(trimmed_data, "/Users/katiarenault/PhD/TOGA/final_metadata.csv", row.names = FALSE)
# 


#############################
### Tree pruning section ###
#############################

library(ape)

# Load the original tree
original_tree <- read.tree("/Users/katiarenault/PhD/Databases/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")
cat("\nOriginal tree information:\n")
cat("Number of tips in original tree:", length(original_tree$tip.label), "\n")

# Modify tip names by removing last two elements separated by "_"
modified_tip_names <- sapply(strsplit(original_tree$tip.label, "_"), function(x) {
  paste(head(x, -2), collapse = "_")
})
original_tree$tip.label <- modified_tip_names

# Get species lists
tree_species <- original_tree$tip.label
data_species <- all_final_data$Scientific_name

# Find matching species
common_species <- intersect(tree_species, data_species)
cat("\nMatching species information:\n")
cat("Number of species in original tree:", length(tree_species), "\n") # 4099
cat("Number of species in data frame:", length(data_species), "\n") # 410
cat("Number of species in common:", length(common_species), "\n") # 391

# Prune tree and filter data
trimmed_tree <- drop.tip(original_tree, setdiff(tree_species, common_species))
trimmed_data <- all_final_data[all_final_data$Scientific_name %in% common_species, ]

cat("\nAfter trimming:\n")
cat("Number of species in pruned tree:", length(trimmed_tree$tip.label), "\n") # 391
cat("Number of species in filtered data frame:", nrow(trimmed_data), "\n") # 391

final_tree_species <- trimmed_tree$tip.label

# Filter the metadata to keep only species present in the pruned tree
final_metadata <- all_final_data %>%
  filter(Scientific_name %in% final_tree_species)

# Verify the results
cat("\nFinal matching verification:\n")
cat("Number of species in pruned tree:", length(final_tree_species), "\n")
cat("Number of species in final metadata:", nrow(final_metadata), "\n")
#write.tree(trimmed_tree, "/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_species_tree_revised.nwk")
#write.csv(final_metadata, "/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv", row.names = FALSE)

#############################
## Body mass and longevity ##
#############################
sophisticated_jewels_palette <- c(
  "#9b383a", "#64a590", "#AA4839", "#AF6C36", "#945a87",
  "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
  "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
  "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
  "#c07a7a", "#8a7c6d", "#5c4d8b", "#4878a0", "#3D7B68"
)

if (length(orders) > length(sophisticated_jewels_palette)) {
  # If more orders than colors, recycle colors
  color_mapping <- setNames(
    sophisticated_jewels_palette[1:length(sophisticated_jewels_palette)], 
    orders[1:length(sophisticated_jewels_palette)]
  )
  warning("More orders than available colors. Some colors will be recycled.")
} else {
  # Use just enough colors from the palette
  color_mapping <- setNames(sophisticated_jewels_palette[1:length(orders)], orders)
}
color_df <- data.frame(
  order = names(color_mapping),
  color = color_mapping
)

species_data <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")
species_data$log_longevity <- log10(species_data$maximum_longevity_y)
species_data$log_mass <- log10(species_data$adult_body_mass_g)
mass_longevity_df <- data.frame(
  log_mass = species_data$log_mass,
  log_longevity = species_data$log_longevity,
  order = species_data$order
)

# calculate correlation
mass_longevity_cor <- cor(mass_longevity_df$log_mass, mass_longevity_df$log_longevity)
mass_longevity_rho <- cor(mass_longevity_df$log_mass, mass_longevity_df$log_longevity, 
                          method = "spearman")

p_mass_longevity <- ggplot(mass_longevity_df, aes(x = log_mass, y = log_longevity, color = order)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  geom_smooth(method = "lm",  linetype = "dashed", color = "gray40", se = TRUE, alpha = 0.2) +
  labs(
    title = "Relationship between body mass and longevity",
    subtitle = paste0("Pearson's r = ", round(mass_longevity_cor, 3), 
                      ", Spearman's ρ = ", round(mass_longevity_rho, 3)),
    x = "Body mass (g, log10)",
    y = "Maximum longevity (yrs, log10)"
  ) +
  theme_bw() + labs(fill = "Order", color = "Order") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

print(p_mass_longevity)
#ggsave("/Users/katiarenault/PhD/TOGA/revised_results/plots/body_mass_longevity.png", p_mass_longevity, width = 8, height = 10)

####################
# Load your data first
species_data <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")

# Define your palette
sophisticated_jewels_palette <- c(
  "#9b383a", "#64a590", "#AA4839", "#AF6C36", "#945a87",
  "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
  "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
  "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
  "#c07a7a", "#8a7c6d", "#5c4d8b", "#4878a0", "#3D7B68"
)

# FIX: Get the actual orders from your data
orders <- unique(species_data$order)
orders <- orders[!is.na(orders)]  # Remove any NA values

cat("Number of unique orders in your data:", length(orders), "\n")
cat("Orders found:", paste(orders, collapse = ", "), "\n")

# Create color mapping with the correct orders
if (length(orders) > length(sophisticated_jewels_palette)) {
  # If more orders than colors, expand the palette by recycling
  expanded_palette <- rep(sophisticated_jewels_palette, ceiling(length(orders) / length(sophisticated_jewels_palette)))
  color_mapping <- setNames(expanded_palette[1:length(orders)], orders)
  warning("More orders than available colors. Colors will be recycled.")
} else {
  # Use just enough colors from the palette
  color_mapping <- setNames(sophisticated_jewels_palette[1:length(orders)], orders)
}

# Verify the color mapping
cat("Color mapping created for", length(color_mapping), "orders\n")
print(color_mapping)

# Create the data for plotting
species_data$log_longevity <- log10(species_data$maximum_longevity_y)
species_data$log_mass <- log10(species_data$adult_body_mass_g)

mass_longevity_df <- data.frame(
  log_mass = species_data$log_mass,
  log_longevity = species_data$log_longevity,
  order = species_data$order
) %>%
  filter(!is.na(order) & !is.na(log_mass) & !is.na(log_longevity))  # Remove any rows with NA values

# Check how many species per order
order_counts <- table(mass_longevity_df$order)
cat("\nSpecies count per order:\n")
print(order_counts)

# Calculate correlations
mass_longevity_cor <- cor(mass_longevity_df$log_mass, mass_longevity_df$log_longevity)
mass_longevity_rho <- cor(mass_longevity_df$log_mass, mass_longevity_df$log_longevity, method = "spearman")

# Create the plot
p_mass_longevity <- ggplot(mass_longevity_df, aes(x = log_mass, y = log_longevity, color = order)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  geom_smooth(method = "lm", linetype = "dashed", color = "gray40", se = TRUE, alpha = 0.2) +
  labs(
    title = "Relationship between body mass and longevity",
    subtitle = paste0("Pearson's r = ", round(mass_longevity_cor, 3),
                      ", Spearman's ρ = ", round(mass_longevity_rho, 3)),
    x = "Body mass (g, log10)",
    y = "Maximum longevity (yrs, log10)"
  ) +
  theme_bw() + 
  labs(fill = "Order", color = "Order") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

print(p_mass_longevity)
ggsave("/Users/katiarenault/PhD/TOGA/revised_results/plots/body_mass_longevity.png", p_mass_longevity, width = 11, height = 8)


######################
###### Heatmaps ######

library(tidyverse)
library(viridis)
library(fgsea)

# Load your data
gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber_Annotated.tsv", 
                           row.names = "t_gene", sep = '\t')
rownames(gene_copy_data) <- gene_copy_data$t_symbol
gene_copy_data <- gene_copy_data %>% dplyr::select(-t_symbol)
gene_copy_data <- gene_copy_data[apply(gene_copy_data, 1, function(x) sum(x == 0, na.rm=TRUE)/sum(!is.na(x))) <= 0.5, ]

metadata <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")

# Load pathway databases
pathways.reactome <- gmtPathways("/Users/katiarenault/Desktop/PhD/Human_GSEA/h.all.v2023.2.Hs.symbols.gmt")
pathways.kegg <- gmtPathways("/Users/katiarenault/Desktop/PhD/Human_GSEA/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
pathways.hallmark <- gmtPathways("/Users/katiarenault/Desktop/PhD/Human_GSEA/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
all_pathways <- c(pathways.kegg, pathways.hallmark, pathways.reactome)

cat("Loaded", length(all_pathways), "pathways total\n")

# Your sophisticated jewel palette
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
  }
  rep(sophisticated_jewels_palette, length.out = n)
}

cat("Gene copy data dimensions:", dim(gene_copy_data), "\n")
cat("Metadata dimensions:", dim(metadata), "\n")

# Function to create pathway-based heatmap
create_pathway_clustered_heatmap <- function(gene_copy_matrix, metadata_df, pathway_list, top_n_pathways = 5) {
  
  # Find longevity and order columns
  longevity_cols <- grep("lifespan|longevity|age", colnames(metadata_df), ignore.case = TRUE, value = TRUE)
  species_cols <- grep("species|name", colnames(metadata_df), ignore.case = TRUE, value = TRUE)
  order_cols <- grep("order", colnames(metadata_df), ignore.case = TRUE, value = TRUE)
  
  longevity_col <- longevity_cols[1]
  species_col <- species_cols[1]
  order_col <- if(length(order_cols) > 0) order_cols[1] else NULL
  
  cat("Using longevity column:", longevity_col, "\n")
  cat("Using species column:", species_col, "\n")
  if(!is.null(order_col)) cat("Using order column:", order_col, "\n")
  
  # Prepare metadata
  metadata_clean <- metadata_df %>%
    filter(!is.na(.data[[longevity_col]]) & .data[[longevity_col]] > 0) %>%
    mutate(
      Species_clean = gsub("[^A-Za-z0-9_]", "_", .data[[species_col]]),
      Species_display = gsub("_", " ", .data[[species_col]])
    ) %>%
    arrange(.data[[longevity_col]])
  
  # Add order information
  if(!is.null(order_col)) {
    metadata_clean <- metadata_clean %>%
      mutate(Order = .data[[order_col]])
  } else {
    metadata_clean <- metadata_clean %>%
      mutate(Order = "Unknown")
  }
  
  # Match species
  common_species <- intersect(colnames(gene_copy_matrix), metadata_clean$Species_clean)
  if(length(common_species) < 5) {
    common_species <- intersect(colnames(gene_copy_matrix), metadata_clean[[species_col]])
    metadata_clean$Species_clean <- metadata_clean[[species_col]]
    metadata_clean$Species_display <- gsub("_", " ", metadata_clean[[species_col]])
  }
  
  cat("Found", length(common_species), "overlapping species\n")
  
  # Filter data for all genes
  filtered_copy <- gene_copy_matrix[, common_species, drop = FALSE]
  metadata_final <- metadata_clean[metadata_clean$Species_clean %in% common_species, ]
  
  # Calculate correlations for ALL genes
  longevity_values <- metadata_final[[longevity_col]][match(common_species, metadata_final$Species_clean)]
  
  cat("Calculating correlations for ALL genes...\n")
  correlations <- apply(filtered_copy, 1, function(gene_counts) {
    if(sum(!is.na(gene_counts)) < 10) return(NA)
    cor(gene_counts, longevity_values, method = "spearman", use = "complete.obs")
  })
  
  # Remove genes with no correlation data
  valid_genes <- names(correlations)[!is.na(correlations)]
  correlations_clean <- correlations[!is.na(correlations)]
  filtered_copy_clean <- filtered_copy[valid_genes, ]
  
  cat("Found", length(valid_genes), "genes with valid correlations\n")
  gene_variance <- apply(filtered_copy, 1, function(gene_counts) {
    var(gene_counts, na.rm = TRUE)
  })
  
  # Filter genes based on variance - keep top 50% most variable genes
  variance_threshold <- quantile(gene_variance, 0.5, na.rm = TRUE)
  high_variance_genes <- names(gene_variance)[gene_variance >= variance_threshold]
  
  # Remove genes with no correlation data or low variance
  valid_genes <- intersect(names(correlations)[!is.na(correlations)], high_variance_genes)
  correlations_clean <- correlations[valid_genes]
  filtered_copy_clean <- filtered_copy[valid_genes, ]
  
  cat("Found", length(valid_genes), "genes with valid correlations and top 50% variance\n")
  
  
  # Find pathway overlaps with our genes
  cat("Analyzing pathway overlaps...\n")
  pathway_overlaps <- lapply(pathway_list, function(pathway_genes) {
    overlap <- intersect(pathway_genes, valid_genes)
    list(genes = overlap, count = length(overlap))
  })
  
  # Get pathway statistics
  pathway_stats <- data.frame(
    Pathway = names(pathway_overlaps),
    Gene_Count = sapply(pathway_overlaps, function(x) x$count),
    stringsAsFactors = FALSE
  ) %>%
    filter(Gene_Count >= 10) %>%  # Only pathways with 10+ genes
    arrange(desc(Gene_Count))
  
  cat("Found", nrow(pathway_stats), "pathways with 10+ genes\n")
  
  # Calculate mean correlation for each pathway
  pathway_stats$Mean_Abs_Correlation <- sapply(pathway_stats$Pathway, function(pw) {
    pathway_genes <- pathway_overlaps[[pw]]$genes
    if(length(pathway_genes) > 0) {
      mean(abs(correlations_clean[pathway_genes]), na.rm = TRUE)
    } else {
      0
    }
  })
  
  # Select top pathways by gene count and correlation
  pathway_stats <- pathway_stats %>%
    arrange(desc(Gene_Count), desc(Mean_Abs_Correlation))
  
  top_pathways <- head(pathway_stats, top_n_pathways)
  
  cat("Top", top_n_pathways, "pathways selected:\n")
  for(i in 1:nrow(top_pathways)) {
    cat(sprintf("%d. %s: %d genes, mean |r| = %.3f\n", 
                i, top_pathways$Pathway[i], top_pathways$Gene_Count[i], 
                top_pathways$Mean_Abs_Correlation[i]))
  }
  
  # Create gene-pathway mapping with unique assignments
  gene_pathway_map <- data.frame(
    Gene = character(),
    Pathway = character(),
    Correlation = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Track which genes have been assigned
  assigned_genes <- character()
  
  # Assign genes to pathways, starting with smallest pathway
  # This ensures genes go to their most specific pathway
  pathway_order <- top_pathways$Pathway[order(top_pathways$Gene_Count)]
  
  cat("Assigning genes to pathways (smallest first to avoid duplicates):\n")
  
  for(pw in pathway_order) {
    pathway_genes <- pathway_overlaps[[pw]]$genes
    # Only take genes not already assigned
    new_genes <- setdiff(pathway_genes, assigned_genes)
    
    if(length(new_genes) > 0) {
      gene_pathway_map <- rbind(gene_pathway_map, data.frame(
        Gene = new_genes,
        Pathway = pw,
        Correlation = correlations_clean[new_genes],
        stringsAsFactors = FALSE
      ))
      
      assigned_genes <- c(assigned_genes, new_genes)
      cat(sprintf("  %s: assigned %d new genes (total: %d)\n", 
                  gsub("REACTOME_|KEGG_|HALLMARK_", "", pw), 
                  length(new_genes), length(assigned_genes)))
    }
  }
  
  cat("Total unique gene assignments to top pathways:", nrow(gene_pathway_map), "\n")
  
  # Verify no duplicates
  if(any(duplicated(gene_pathway_map$Gene))) {
    stop("Error: Duplicate genes found after assignment!")
  }
  
  # Order genes by pathway, then by correlation within pathway
  # Use the order we want pathways to appear (largest to smallest for display)
  display_pathway_order <- top_pathways$Pathway[order(-top_pathways$Gene_Count)]
  
  gene_pathway_final <- gene_pathway_map %>%
    mutate(Abs_Correlation = abs(Correlation)) %>%
    arrange(match(Pathway, display_pathway_order), desc(Abs_Correlation))
  
  # Limit genes per pathway for visualization
  max_genes_per_pathway <- 100
  gene_pathway_final <- gene_pathway_final %>%
    group_by(Pathway) %>%
    slice_head(n = max_genes_per_pathway) %>%
    ungroup()
  
  final_genes <- gene_pathway_final$Gene
  
  cat("Final selection:", length(final_genes), "unique genes from", length(unique(gene_pathway_final$Pathway)), "pathways\n")
  
  # Double-check no duplicates in final selection
  if(any(duplicated(final_genes))) {
    cat("Warning: Found duplicates in final genes, removing...\n")
    gene_pathway_final <- gene_pathway_final[!duplicated(gene_pathway_final$Gene), ]
    final_genes <- gene_pathway_final$Gene
  }
  
  # Order species by longevity
  species_order <- metadata_final$Species_clean[order(metadata_final[[longevity_col]])]
  filtered_copy_final <- filtered_copy_clean[final_genes, species_order]
  
  # Apply row-wise scaling (each gene normalized to its own baseline)
  cat("Applying row-wise scaling for relative gene expression visualization...\n")
  
  # First, replace 0s with NA to handle missing gene copies
  filtered_copy_clean_na <- filtered_copy_final
  filtered_copy_clean_na[filtered_copy_clean_na == 0] <- NA
  
  # Calculate row-wise scaling: each value as fold-change from gene's minimum value
  scaled_copy_final <- t(apply(filtered_copy_clean_na, 1, function(gene_row) {
    # Remove NA values for calculation
    valid_values <- gene_row[!is.na(gene_row)]
    
    if(length(valid_values) == 0) {
      return(gene_row)  # Return original if all NA
    }
    
    if(length(valid_values) == 1) {
      # If only one valid value, set it as baseline (1×)
      gene_row[!is.na(gene_row)] <- 1
      return(gene_row)
    }
    
    # Get minimum value as baseline
    min_val <- min(valid_values)
    
    # Calculate fold change from baseline
    fold_change <- gene_row / min_val
    
    # Cap extreme values for better visualization
    fold_change[fold_change > 10 & !is.na(fold_change)] <- 10  # Cap at 10-fold increase
    
    return(fold_change)
  }))
  
  cat("Row-wise scaling complete. 0 copy numbers converted to NA (grey tiles).\n")
  cat("Values now represent fold-change from each gene's baseline.\n")
  
  # Convert to long format
  plot_data <- scaled_copy_final %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(cols = -Gene, names_to = "Species", values_to = "Fold_Change") %>%
    mutate(
      Gene_clean = gsub("_", " ", Gene),
      Species_display = gsub("_", " ", Species)
    ) %>%
    left_join(
      metadata_final %>% 
        select(Species_clean, Order, all_of(longevity_col)) %>% 
        rename(Species = Species_clean, Longevity = all_of(longevity_col)),
      by = "Species"
    ) %>%
    left_join(
      gene_pathway_final %>% select(Gene, Pathway, Correlation),
      by = "Gene"
    )
  
  # Create proper factor ordering maintaining pathway grouping
  gene_display_order <- gsub("_", " ", gene_pathway_final$Gene)
  
  plot_data$Species_display <- factor(plot_data$Species_display, 
                                      levels = gsub("_", " ", species_order))
  plot_data$Gene_clean <- factor(plot_data$Gene_clean, 
                                 levels = gene_display_order)
  
  # Create pathway break positions for visual separation
  pathway_breaks <- gene_pathway_final %>%
    mutate(Position = match(Gene, final_genes)) %>%
    group_by(Pathway) %>%
    summarise(
      Start = min(Position),
      End = max(Position),
      Count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      Mid = (Start + End) / 2,
      Label = gsub("HALLMARK_|KEGG_|REACTOME_", "", Pathway),
      Label = gsub("_", " ", Label)  # Replace underscores with spaces
    )
  
  cat("Creating main heatmap with", nrow(plot_data), "data points...\n")
  
  # Create main heatmap
  main_plot <- ggplot(plot_data, aes(x = Species_display, y = Gene_clean)) +
    geom_tile(aes(fill = Fold_Change), color = "white", size = 0.05) +
    scale_fill_gradient2(
      low = "#4575b4", 
      mid = "#f9e1e1", 
      high = "#9b383a",
      midpoint = 1,  # 1 = baseline (no change)
      na.value = "grey90",
      name = "Fold change\nfrom gene\nnaseline",
      breaks = c(0.1, 1, 2, 4, 8),
      labels = c("missing", "1×", "2×", "4×", "8×"),
      limits = c(0.1, 10)  # Set reasonable limits to avoid infinite values
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 0),
      axis.text.y = element_text(size = 3),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "left",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.margin = margin(20, 210, 40, 40)  # Extra right margin for pathway labels
    ) +
    labs(
      title = "Relative Gene Copy Number Changes Across Species Longevity Spectrum",
      subtitle = paste("Top", top_n_pathways, "pathways • Each gene scaled to its own baseline • Grey = gene absent •", 
                       length(unique(plot_data$Species)), "species ordered by longevity")
    )
  
  # Add pathway labels on the right
  if(nrow(pathway_breaks) > 0) {
    main_plot <- main_plot +
      annotate("text", 
               x = length(unique(plot_data$Species_display)) + 2,  # Move labels closer
               y = pathway_breaks$Mid,
               label = pathway_breaks$Label,
               hjust = 0, 
               size = 3,  # Slightly larger text
               color = "#333333", 
               fontface = "bold") +
      coord_cartesian(clip = "off", 
                      xlim = c(0.5, length(unique(plot_data$Species_display))), 
                      expand = FALSE)  # Prevent expansion beyond data
  }
  
  cat("Main plot created successfully\n")
  
  # Create order color bar
  unique_orders <- unique(metadata_final$Order)
  order_colors <- get_jewel_palette(length(unique_orders))
  names(order_colors) <- unique_orders
  
  order_data <- metadata_final %>%
    filter(Species_clean %in% species_order) %>%
    arrange(match(Species_clean, species_order)) %>%
    mutate(
      Species_display = gsub("_", " ", Species_clean),
      y_value = 1
    )
  
  order_data$Species_display <- factor(order_data$Species_display, 
                                       levels = levels(plot_data$Species_display))
  
  cat("Creating order bar...\n")
  order_plot <- ggplot(order_data, aes(x = Species_display, y = y_value)) +
    geom_tile(aes(fill = Order), color = "white", size = 0.1, height = 0.8) +
    scale_fill_manual(values = order_colors, name = "Order") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    ) +
    guides(fill = guide_legend(nrow = 3))
  
  cat("Order plot created successfully\n")
  
  return(list(
    main_plot = main_plot,
    order_plot = order_plot,
    correlations = correlations_clean,
    pathway_stats = pathway_stats,
    pathway_breaks = pathway_breaks,
    order_colors = order_colors,
    plot_data = plot_data,
    gene_pathway_map = gene_pathway_final
  ))
}

# Create the comprehensive pathway-clustered heatmap
cat("Creating comprehensive pathway-clustered longevity heatmap...\n")
result <- create_pathway_clustered_heatmap(gene_copy_data, metadata, all_pathways, top_n_pathways = 5)

# Display results
if(!is.null(result$main_plot)) {
  cat("\n=== DISPLAYING PATHWAY-CLUSTERED HEATMAP ===\n")
  print(result$main_plot)
} else {
  cat("Failed to create main plot\n")
}

if(!is.null(result$order_plot)) {
  cat("\n=== DISPLAYING TAXONOMIC ORDER BAR ===\n")
  print(result$order_plot)
} else {
  cat("Failed to create order plot\n")
}

# Create combined plot if both exist
if(!is.null(result$main_plot) && !is.null(result$order_plot)) {
  if(require(patchwork, quietly = TRUE)) {
    cat("\n=== DISPLAYING COMBINED PLOT ===\n")
    combined <- result$main_plot / result$order_plot + 
      plot_layout(heights = c(12, 1))
    print(combined)
  }
}

# Summary information
cat("\n=== PATHWAY ANALYSIS SUMMARY ===\n")
cat("Top pathways by gene count and correlation:\n")
top_pathway_summary <- result$pathway_stats %>% 
  arrange(desc(Gene_Count)) %>%
  head(10) %>%
  mutate(
    Clean_Name = gsub("HALLMARK_|KEGG_|REACTOME_", "", Pathway),
    Clean_Name = str_wrap(Clean_Name, 50)
  )
print(top_pathway_summary %>% select(Clean_Name, Gene_Count, Mean_Abs_Correlation))

cat("\nPathways included in heatmap:\n")
if(nrow(result$pathway_breaks) > 0) {
  pathway_summary <- result$pathway_breaks %>%
    mutate(Clean_Label = gsub("\n.*", "", Label)) %>%
    select(Clean_Label, Count)
  print(pathway_summary)
} else {
  cat("No pathway breaks found\n")
}

cat("\nGenes per pathway in visualization:\n")
pathway_gene_counts <- table(result$gene_pathway_map$Pathway)
print(sort(pathway_gene_counts, decreasing = TRUE))

cat("\nOverall correlation distribution:\n")
cor_summary <- summary(abs(result$correlations))
print(cor_summary)

cat("\nTotal genes analyzed:", length(result$correlations), "\n")
cat("Total genes in visualization:", nrow(result$gene_pathway_map), "\n")
cat("Genes with |r| > 0.3:", sum(abs(result$correlations) > 0.3, na.rm = TRUE), "\n")
cat("Genes with |r| > 0.5:", sum(abs(result$correlations) > 0.5, na.rm = TRUE), "\n")

# Show some example high-correlation genes from each pathway
cat("\nTop correlated genes by pathway:\n")
top_genes_by_pathway <- result$gene_pathway_map %>%
  group_by(Pathway) %>%
  slice_max(order_by = abs(Correlation), n = 3) %>%
  ungroup() %>%
  arrange(Pathway, desc(abs(Correlation)))

for(pw in unique(top_genes_by_pathway$Pathway)) {
  cat("\n", gsub("REACTOME_", "", pw), ":\n")
  pathway_genes <- top_genes_by_pathway %>% filter(Pathway == pw)
  for(i in 1:nrow(pathway_genes)) {
    cat(sprintf("  %s: r = %.3f\n", pathway_genes$Gene[i], pathway_genes$Correlation[i]))
  }
}

cat("\nFold-change visualization complete!\n")
cat("Each gene row now shows relative changes from its own baseline.\n")
cat("White = baseline (1×), Red = higher copies, Blue = lower copies\n")