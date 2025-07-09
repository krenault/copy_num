# color_mappings.R
# This file contains consistent color mappings for mammalian orders

#' Sophisticated Jewels Color Palette
#' 
#' A vector of 25 carefully selected jewel-toned colors
sophisticated_jewels_palette <- c(
  "#9b383a", "#64a590", "#AA4839", "#AF6C36", "#945a87",
  "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
  "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
  "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
  "#c07a7a", "#8a7c6d", "#5c4d8b", "#4878a0", "#3D7B68"
)

#' Fixed Color Mapping for Mammalian Orders
#' 
#' A named vector where each mammalian order is assigned a specific color
fixed_color_mapping <- c(
  "Carnivora" = "#9b383a",
  "Rodentia" = "#64a590",
  "Artiodactyla" = "#AA4839",
  "Primates" = "#AF6C36",
  "Dasyuromorphia" = "#945a87",
  "Chiroptera" = "#b68ab9",
  "Cetacea" = "#9c6e5a",
  "Pilosa" = "#9fbd8b",
  "Afrosoricida" = "#a85472",
  "Soricomorpha" = "#86acb9",
  "Perissodactyla" = "#20854c",
  "Didelphimorphia" = "#9c534a",
  "Microbiotheria" = "#2c6e5f",
  "Sirenia" = "#31a6ad",
  "Macroscelidea" = "#f0a0a3",
  "Proboscidea" = "#c57251",
  "Erinaceomorpha" = "#856890",
  "Diprotodontia" = "#7e9960",
  "Hyracoidea" = "#7d9a9e",
  "Lagomorpha" = "#B65B32",
  "Pholidota" = "#c07a7a",
  "Monotremata" = "#8a7c6d",
  "Scandentia" = "#5c4d8b"
)

#' Create Color Mapping for Orders in a Dataset
#'
#' @param order_vector A vector containing order names (e.g., species_data$order)
#' @return A named vector of colors for the orders
#' @examples
#' color_mapping <- create_color_mapping(species_data$order)
create_color_mapping <- function(order_vector) {
  orders <- unique(order_vector)
  orders <- orders[!is.na(orders)] 
  
  cat("Number of unique orders in your data:", length(orders), "\n")
  cat("Orders found:", paste(orders, collapse = ", "), "\n")
  
  # Start with fixed colors for known orders
  mapping <- fixed_color_mapping[names(fixed_color_mapping) %in% orders]
  
  # Handle any new orders not in our fixed mapping
  new_orders <- setdiff(orders, names(fixed_color_mapping))
  if (length(new_orders) > 0) {
    remaining_colors <- setdiff(sophisticated_jewels_palette, mapping)
    if (length(new_orders) > length(remaining_colors)) {
      warning("Not enough unique colors for all orders. Some colors will be recycled.")
      remaining_colors <- rep(remaining_colors, length.out = length(new_orders))
    }
    new_mapping <- setNames(remaining_colors[1:length(new_orders)], new_orders)
    mapping <- c(mapping, new_mapping)
  }
  
  cat("Color mapping created for", length(mapping), "orders\n")
  return(mapping)
}