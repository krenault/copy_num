## Katia Renault


########################################################################################################
# 1. Phylogeny of available species
########################################################################################################
library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(ape)
library(phytools)

source("/Users/katiarenault/Documents/Github/copy_num/scripts/FUN_color_mappings.R")
color_mapping <- create_color_mapping(species_data$order)

# Load data
newick_tree <- "/Users/katiarenault/Documents/Github/copy_num/data/raxml_final_species_tree_revised.nwk"
tree <- read.tree(newick_tree)
metadata <- read.csv("/Users/katiarenault/Documents/Github/copy_num/data/raxml_final_metadata_revised.csv")
gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber_Annotated.tsv", 
                           row.names = 1, sep = '\t', check.names = FALSE)
species_in_gene_data <- colnames(gene_copy_data)
tree_species <- tree$tip.label
species_to_keep <- intersect(tree_species, species_in_gene_data)
pruned_tree <- drop.tip(tree, setdiff(tree_species, species_to_keep))

plot_data <- metadata %>%
  filter(Scientific_name %in% species_to_keep) %>%
  dplyr::select(Scientific_name, order, maximum_longevity_y)
order_colors <- color_mapping
tips_data <- plot_data %>%
  dplyr::select(Scientific_name, order) %>%
  dplyr::rename(label = `Scientific_name`)

tree_data <- as_tibble(pruned_tree)
tree_with_orders <- full_join(tree_data, tips_data, by = "label")
get_node_orders <- function(tree, tip_orders) {
  node_orders <- rep(NA, tree$Nnode + length(tree$tip.label))
  for (i in 1:length(tree$tip.label)) {
    tip_name <- tree$tip.label[i]
    if (tip_name %in% names(tip_orders)) {
      node_orders[i] <- tip_orders[tip_name]
    }
  }
  return(node_orders)
}
tip_order_vector <- setNames(plot_data$order, plot_data$Scientific_name)
all_orders <- get_node_orders(pruned_tree, tip_order_vector)

# Create the fan-shaped tree plot with equal branch lengths and colored branches
p <- ggtree(pruned_tree, layout = 'fan', branch.length = "none") +
  geom_tree(aes(color = all_orders[node])) +
  scale_color_manual(
    name = "Order",
    values = order_colors,
    na.value = "#504f4f",
    guide = 'none'
  )

# Add the longevity bars as fruit
p_with_bars <- p +
  geom_fruit(
    data = plot_data,
    geom = geom_col,
    mapping = aes(
      x = maximum_longevity_y,
      y = Scientific_name,
      fill = order
    ),
    orientation = "y",
    width = 0.6,
    color = "#504f4f",
    size = 0.5,
    offset = 0.1
  ) +
  scale_fill_manual(
    name = "Order",
    values = order_colors
  ) +
  theme(
    legend.position = "right",
    plot.margin = margin(100, 100, 100, 100)
  )

print(p_with_bars)
write.tree(pruned_tree, file = "/Users/katiarenault/Documents/Github/copy_num/data/pruned_species_tree.nwk")
cat("Original tree had", length(tree$tip.label), "species\n")
cat("Pruned tree has", length(pruned_tree$tip.label), "species\n")
cat("Number of orders represented:", n_orders, "\n")
ggsave("/Users/katiarenault/Documents/Github/copy_num/plots/phylogenetic_tree_with_longevity.pdf", p_with_bars, width = 25, height = 25)