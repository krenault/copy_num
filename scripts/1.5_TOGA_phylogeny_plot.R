library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)

# Load the color palette function
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
    return(sophisticated_jewels_palette[1:min(n, length(sophisticated_jewels_palette))])
  }
}

# Load data
newick_tree <- "/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_species_tree_revised.nwk"
tree <- read.tree(newick_tree)
metadata <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")
gene_copy_data <- read.csv("/Users/katiarenault/Desktop/PhD/TOGA/results/Mammalian_All_species_copynumber_updated.tsv", 
                           row.names = 1, sep = '\t', check.names = FALSE)

# Get species names from gene copy data (excluding the first column which is likely gene names)
species_in_gene_data <- colnames(gene_copy_data)

# Prune the tree to keep only species present in gene copy data
tree_species <- tree$tip.label
species_to_keep <- intersect(tree_species, species_in_gene_data)

# Prune the tree
pruned_tree <- drop.tip(tree, setdiff(tree_species, species_to_keep))

# Prepare data for the barplot
plot_data <- metadata %>%
  filter(Scientific_name %in% species_to_keep) %>%
  dplyr::select(Scientific_name, order, maximum_longevity_y)

# Create a mapping of orders to colors
unique_orders <- unique(plot_data$order)
n_orders <- length(unique_orders)
order_colors <- setNames(get_jewel_palette(n_orders), unique_orders)

# Map orders to internal nodes in the tree
# First get the tips and their orders
tips_data <- plot_data %>%
  dplyr::select(Scientific_name, order) %>%
  dplyr::rename(label = `Scientific_name`)

# Create tree data with order information
tree_data <- as_tibble(pruned_tree)
tree_with_orders <- full_join(tree_data, tips_data, by = "label")

# Function to get most common order for internal nodes
get_node_orders <- function(tree, tip_orders) {
  node_orders <- rep(NA, tree$Nnode + length(tree$tip.label))
  
  # Set tip orders
  for (i in 1:length(tree$tip.label)) {
    tip_name <- tree$tip.label[i]
    if (tip_name %in% names(tip_orders)) {
      node_orders[i] <- tip_orders[tip_name]
    }
  }
  
  return(node_orders)
}

# Create a named vector of tip orders
tip_order_vector <- setNames(plot_data$order, plot_data$Scientific_name)
all_orders <- get_node_orders(pruned_tree, tip_order_vector)

# Plot the tree with longevity bars
# Create the base circular tree
p <- ggtree(pruned_tree, layout = 'fan') +
  # Add colored branches based on orders
  geom_tree(aes(color = all_orders[node])) +
  # Apply custom color palette for branches
  scale_color_manual(
    name = "Order",
    values = order_colors,
    na.value = "#504f4f" , # Color for internal branches
    guide = 'none'
  )

# Add bars showing maximum longevity years
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
    color = "#504f4f",   # Dark grey outline for all bars
    size = 0.5,           # Outline thickness
    offset = 0.1
  ) +
  scale_fill_manual(
    name = "Order",
    values = order_colors
  ) +
  theme(
    legend.position = "right",
    plot.margin = margin(100, 100, 100, 100)  # Increase margins (top, right, bottom, left)
  )

# Display the plot
print(p_with_bars)

# Save the pruned tree
write.tree(pruned_tree, file = "pruned_species_tree.nwk")

# Print summary
cat("Original tree had", length(tree$tip.label), "species\n")
cat("Pruned tree has", length(pruned_tree$tip.label), "species\n")
cat("Number of orders represented:", n_orders, "\n")

# Export the plot (uncomment to save)
ggsave("/Users/katiarenault/Desktop/PhD/TOGA/results/phylogenetic_tree_with_longevity.pdf", p_with_bars, width = 25, height = 25)