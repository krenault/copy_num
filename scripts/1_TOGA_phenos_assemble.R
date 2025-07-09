## Katia Renault


########################################################################################################
# 1. Phenos assembly for species that have TOGA
########################################################################################################
# generating a metadata file for the species for which I have copy number data and lifespan information 
# uses amniote and anage data
########################################################################################################
###############
### Amniote ###
###############
library(dplyr)
copy_number_animal_data <- read.csv("/Users/katiarenault/PhD/TOGA/unique_species_list.csv", header = TRUE)
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
# merge the two data frames
final_merged_data <- merged_data %>%
  full_join(merged_data_anage, by = "Scientific_name") %>%
  mutate(
    # default to anage values if available, otherwise keep amniote values
    maximum_longevity_y = coalesce(Maximum.longevity..yrs., maximum_longevity_y),
    adult_body_mass_g = coalesce(Adult.weight..g., adult_body_mass_g),
    order = coalesce(Order, order),
    common_name = coalesce(Common.name, common_name)
  ) %>%
  dplyr::select(-Maximum.longevity..yrs., -Adult.weight..g., -Order, -Common_name.x, -Common_name.y, -Common.name)
# count the number of updated values
updated_longevity <- sum(!is.na(merged_data_anage$Maximum.longevity..yrs.) & is.na(merged_data$maximum_longevity_y))
updated_mass <- sum(!is.na(merged_data_anage$Adult.weight..g.) & is.na(merged_data$adult_body_mass_g))
cat("Number of maximum longevity values updated from anage:", updated_longevity, "\n") # 11 
cat("Number of adult body mass values updated from anage:", updated_mass, "\n") # 29
unmatched_data <- final_merged_data %>%
  filter(is.na(maximum_longevity_y)) %>%
  dplyr::select(-maximum_longevity_y)
new_matches <- unmatched_data %>%
  left_join(
    amniote_df %>%
      dplyr::select(common_name, maximum_longevity_y) %>%
      mutate(common_name = tolower(common_name)), 
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

########################################################################################################
# 2. Tree pruning section 
########################################################################################################
# generating a tree file for the species for which I have metadata for phylogenetically informed analysis
# output: raxml_final_species_tree_revised.nwk and raxml_final_metadata_revised.csv
########################################################################################################

library(ape)
original_tree <- read.tree("/Users/katiarenault/PhD/Databases/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")
cat("\nOriginal tree information:\n")
cat("Number of tips in original tree:", length(original_tree$tip.label), "\n")
modified_tip_names <- sapply(strsplit(original_tree$tip.label, "_"), function(x) {
  paste(head(x, -2), collapse = "_")
})
original_tree$tip.label <- modified_tip_names
# get species lists
tree_species <- original_tree$tip.label
data_species <- all_final_data$Scientific_name
# matching species
common_species <- intersect(tree_species, data_species)
cat("\nMatching species information:\n")
cat("Number of species in original tree:", length(tree_species), "\n") # 4099
cat("Number of species in data frame:", length(data_species), "\n") # 410
cat("Number of species in common:", length(common_species), "\n") # 391
trimmed_tree <- drop.tip(original_tree, setdiff(tree_species, common_species))
trimmed_data <- all_final_data[all_final_data$Scientific_name %in% common_species, ]
cat("\nAfter trimming:\n")
cat("Number of species in pruned tree:", length(trimmed_tree$tip.label), "\n") # 391
cat("Number of species in filtered data frame:", nrow(trimmed_data), "\n") # 391
final_tree_species <- trimmed_tree$tip.label

# filter the metadata to keep only species present in the pruned tree
final_metadata <- all_final_data %>%
  filter(Scientific_name %in% final_tree_species)
cat("\nFinal matching verification:\n")
cat("Number of species in pruned tree:", length(final_tree_species), "\n")
cat("Number of species in final metadata:", nrow(final_metadata), "\n")
#write.tree(trimmed_tree, "/Users/katiarenault/Documents/Github/copy_num/data/raxml_final_species_tree_revised.nwk")
#write.csv(final_metadata, "/Users/katiarenault/Documents/Github/copy_num/data/raxml_final_metadata_revised.csv", row.names = FALSE)


########################################################################################################
# 3. Body mass and longevity 
########################################################################################################
# generating plot illustrating body mass and longevity correlation between species
########################################################################################################

species_data <- read.csv("/Users/katiarenault/Documents/Github/copy_num/data/raxml_final_metadata_revised.csv")
source("/Users/katiarenault/Documents/Github/copy_num/scripts/FUN_color_mappings.R")
color_mapping <- create_color_mapping(species_data$order)
print(color_mapping)

species_data$log_longevity <- log10(species_data$maximum_longevity_y)
species_data$log_mass <- log10(species_data$adult_body_mass_g)

mass_longevity_df <- data.frame(
  log_mass = species_data$log_mass,
  log_longevity = species_data$log_longevity,
  order = species_data$order
) %>%
  filter(!is.na(order) & !is.na(log_mass) & !is.na(log_longevity))  # Remove any rows with NA values

order_counts <- table(mass_longevity_df$order)
cat("\nSpecies count per order:\n")
print(order_counts)
mass_longevity_cor <- cor(mass_longevity_df$log_mass, mass_longevity_df$log_longevity)
mass_longevity_rho <- cor(mass_longevity_df$log_mass, mass_longevity_df$log_longevity, method = "spearman")

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
ggsave("/Users/katiarenault/Documents/GitHub/copy_num/plots/body_mass_longevity.png", p_mass_longevity, width = 11, height = 8)
