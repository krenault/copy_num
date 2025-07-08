##############################################################################
# Species Longevity Prediction Analysis
# This script analyzes species-specific longevity predictions and determines
# which genes contribute most to differences in predicted lifespans.
##############################################################################

###### run 3_Elastic_net and with seed 928 for naked mole rat and mouse
## - seed 982 for Molossus_molossus and Myotis brandtii
library(ggplot2)
library(dplyr)

#####################################################
# 1. FUNCTION TO ANALYZE SPECIES PREDICTIONS
#####################################################

analyze_species_predictions <- function(species_list, model, x_test_scaled, test_species, 
                                        var_importance, top_n = 20) {
  
  # Create a list to store results
  results <- list()
  
  # Find indices of requested species in the test set
  species_indices <- match(species_list, test_species)
  
  # Check if any species were not found
  not_found <- species_list[is.na(species_indices)]
  if (length(not_found) > 0) {
    warning(paste("Species not found in test set:", paste(not_found, collapse=", ")))
  }
  
  # Remove NA values (species not found)
  species_indices <- species_indices[!is.na(species_indices)]
  
  if (length(species_indices) == 0) {
    stop("None of the requested species were found in the test set")
  }
  
  # Extract feature values for the selected species
  species_features <- x_test_scaled[species_indices, , drop=FALSE]
  
  # Get predictions for these species
  species_preds <- predict(model, newx=species_features, s=best_lambda)
  
  # Calculate gene contributions for each species
  gene_contributions <- list()
  
  for (i in 1:length(species_indices)) {
    idx <- species_indices[i]
    species_name <- test_species[idx]
    
    # Calculate contribution of each gene to the prediction
    # (feature value * coefficient)
    contributions <- species_features[i, ] * var_importance
    
    # Sort by absolute contribution
    abs_contributions <- abs(contributions)
    sorted_idx <- order(abs_contributions, decreasing=TRUE)
    
    # Store top contributing genes
    top_genes <- data.frame(
      Gene = names(contributions)[sorted_idx[1:top_n]],
      Contribution = contributions[sorted_idx[1:top_n]],
      Feature_Value = species_features[i, sorted_idx[1:top_n]],
      Coefficient = var_importance[sorted_idx[1:top_n]],
      stringsAsFactors = FALSE
    )
    
    gene_contributions[[species_name]] <- top_genes
    
    # Store basic prediction info
    results[[species_name]] <- list(
      prediction = species_preds[i],
      actual = y_test[idx],
      top_genes = top_genes
    )
  }
  
  # Compare between species
  if (length(species_indices) >= 2) {
    species_pairs <- combn(names(results), 2, simplify=FALSE)
    
    comparisons <- list()
    
    for (pair in species_pairs) {
      sp1 <- pair[1]
      sp2 <- pair[2]
      
      # Calculate prediction difference
      pred_diff <- results[[sp1]]$prediction - results[[sp2]]$prediction
      
      # Find genes that differ most in their contribution
      all_genes <- unique(c(results[[sp1]]$top_genes$Gene, 
                            results[[sp2]]$top_genes$Gene))
      
      diff_contributions <- data.frame(
        Gene = all_genes,
        Contribution_Diff = NA,
        Feature_Value_Diff = NA,
        stringsAsFactors = FALSE
      )
      
      for (g in all_genes) {
        # Find gene in each species
        idx1 <- match(g, results[[sp1]]$top_genes$Gene)
        idx2 <- match(g, results[[sp2]]$top_genes$Gene)
        
        # Get contribution for each species (0 if not in top genes)
        contrib1 <- if(!is.na(idx1)) results[[sp1]]$top_genes$Contribution[idx1] else 0
        contrib2 <- if(!is.na(idx2)) results[[sp2]]$top_genes$Contribution[idx2] else 0
        
        # Calculate feature value difference
        feat_idx <- which(names(var_importance) == g)
        feat_val1 <- species_features[which(test_species[species_indices] == sp1), feat_idx]
        feat_val2 <- species_features[which(test_species[species_indices] == sp2), feat_idx]
        
        # Store differences
        diff_contributions$Contribution_Diff[diff_contributions$Gene == g] <- contrib1 - contrib2
        diff_contributions$Feature_Value_Diff[diff_contributions$Gene == g] <- feat_val1 - feat_val2
      }
      
      # Sort by absolute contribution difference
      diff_contributions <- diff_contributions[order(abs(diff_contributions$Contribution_Diff), 
                                                     decreasing=TRUE), ]
      
      comparisons[[paste(sp1, "vs", sp2)]] <- list(
        prediction_diff = pred_diff,
        actual_diff = results[[sp1]]$actual - results[[sp2]]$actual,
        gene_contributions = diff_contributions[1:min(top_n, nrow(diff_contributions)), ]
      )
    }
    
    results$comparisons <- comparisons
  }
  
  return(results)
}

#####################################################
# 2. FUNCTION TO VISUALIZE SINGLE SPECIES PREDICTIONS
#####################################################

visualize_species_genes <- function(analysis_results, species_name) {
  # Extract species data
  species_data <- analysis_results[[species_name]]
  
  if (is.null(species_data)) {
    stop(paste("Species not found:", species_name))
  }
  
  # Prepare data for plotting
  genes_df <- species_data$top_genes
  genes_df$Gene <- factor(genes_df$Gene, levels = genes_df$Gene[order(genes_df$Contribution)])
  
  # Add direction column
  genes_df$Direction <- ifelse(genes_df$Contribution > 0, 
                               "Positive impact on longevity", 
                               "Negative impact on longevity")
  
  # Colors
  direction_colors <- c("Positive impact on longevity" = "#9b383a", 
                        "Negative impact on longevity" = "#86acb9")
  
  # Create plot
  p <- ggplot(genes_df, aes(x = Gene, y = Contribution, fill = Direction)) +
    geom_col() +
    scale_fill_manual(values = direction_colors) +
    coord_flip() +
    labs(
      title = paste("Top genes contributing to longevity in", species_name),
      subtitle = paste("Predicted log10 lifespan:", round(species_data$prediction, 3),
                       "| Actual log10 lifespan:", round(species_data$actual, 3)),
      x = "Gene",
      y = "Contribution to Longevity Prediction"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 12)
    )
  
  return(p)
}

#####################################################
# 3. FUNCTION TO VISUALIZE SPECIES COMPARISONS
#####################################################

visualize_species_comparison <- function(analysis_results, comparison_name) {
  # Extract comparison data
  comp_data <- analysis_results$comparisons[[comparison_name]]
  
  if (is.null(comp_data)) {
    stop(paste("Comparison not found:", comparison_name))
  }
  
  # Prepare data for plotting
  genes_df <- comp_data$gene_contributions
  genes_df$Gene <- factor(genes_df$Gene, levels = genes_df$Gene[order(genes_df$Contribution_Diff)])
  
  # Extract species names from comparison name
  species_names <- strsplit(comparison_name, " vs ")[[1]]
  species1 <- species_names[1]
  species2 <- species_names[2]
  
  # Add direction column for clearer legend
  genes_df$Direction <- ifelse(genes_df$Contribution_Diff > 0, 
                               paste("Higher in", species1), 
                               paste("Higher in", species2))
  
  # Colors
  direction_names <- c(paste("Higher in", species1), paste("Higher in", species2))
  direction_colors <- c("#64a590", "#31a6ad")
  names(direction_colors) <- direction_names
  
  # Create plot
  p <- ggplot(genes_df, aes(x = Gene, y = Contribution_Diff, fill = Direction)) +
    geom_col() +
    scale_fill_manual(values = direction_colors, name = "Contribution") +
    coord_flip() +
    labs(
      title = paste("Gene contribution differences:", comparison_name),
      subtitle = paste("Prediction difference:", round(comp_data$prediction_diff, 3),
                       "| Actual difference:", round(comp_data$actual_diff, 3)),
      x = "Gene",
      y = "Difference in Contribution to Longevity Prediction"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 12)
    )
  
  return(p)
}

#####################################################
# 4. FUNCTION TO PRINT SPECIES PREDICTION REPORT
#####################################################

print_species_report <- function(analysis_results, species_name) {
  # Extract species data
  species_data <- analysis_results[[species_name]]
  
  if (is.null(species_data)) {
    stop(paste("Species not found:", species_name))
  }
  
  # Print header
  cat("=====================================================\n")
  cat(paste("LONGEVITY PREDICTION REPORT FOR:", species_name, "\n"))
  cat("=====================================================\n\n")
  
  # Print prediction info
  cat(paste("Predicted log10 lifespan:", round(species_data$prediction, 3), "\n"))
  cat(paste("Actual log10 lifespan:", round(species_data$actual, 3), "\n"))
  cat(paste("Difference:", round(species_data$prediction - species_data$actual, 3), "\n\n"))
  
  # Print top contributing genes
  cat("TOP CONTRIBUTING GENES:\n")
  cat("-----------------------\n\n")
  
  top_genes <- species_data$top_genes
  for (i in 1:nrow(top_genes)) {
    direction <- ifelse(top_genes$Contribution[i] > 0, "POSITIVE", "NEGATIVE")
    cat(sprintf("%2d. %-20s: %+.4f (%s impact)\n", 
                i, 
                top_genes$Gene[i],
                top_genes$Contribution[i],
                direction))
  }
  
  cat("\n")
}

#####################################################
# 5. FUNCTION TO PRINT COMPARISON REPORT
#####################################################

print_comparison_report <- function(analysis_results, comparison_name) {
  # Extract comparison data
  comp_data <- analysis_results$comparisons[[comparison_name]]
  
  if (is.null(comp_data)) {
    stop(paste("Comparison not found:", comparison_name))
  }
  
  # Extract species names
  species_names <- strsplit(comparison_name, " vs ")[[1]]
  species1 <- species_names[1]
  species2 <- species_names[2]
  
  # Print header
  cat("=====================================================\n")
  cat(paste("LONGEVITY COMPARISON:", species1, "vs", species2, "\n"))
  cat("=====================================================\n\n")
  
  # Print prediction info
  cat(paste("Prediction difference:", round(comp_data$prediction_diff, 3),
            "(", species1, "-", species2, ")\n"))
  cat(paste("Actual difference:", round(comp_data$actual_diff, 3),
            "(", species1, "-", species2, ")\n\n"))
  
  # Print top differentiating genes
  cat("TOP DIFFERENTIATING GENES:\n")
  cat("--------------------------\n\n")
  
  genes_df <- comp_data$gene_contributions
  for (i in 1:nrow(genes_df)) {
    favors <- ifelse(genes_df$Contribution_Diff[i] > 0, species1, species2)
    cat(sprintf("%2d. %-20s: %+.4f (favors %s)\n", 
                i, 
                genes_df$Gene[i],
                genes_df$Contribution_Diff[i],
                favors))
  }
  
  cat("\n")
}

#####################################################
# 6. FUNCTION TO ANALYZE GENE-LONGEVITY CORRELATIONS
#####################################################

analyze_longevity_gene_correlations <- function(gene_data, longevity_data, top_n = 20) {
  # Create matrix of gene values
  gene_matrix <- as.matrix(gene_data)
  
  # Calculate correlation for each gene
  cors <- numeric(ncol(gene_matrix))
  names(cors) <- colnames(gene_matrix)
  
  for (i in 1:ncol(gene_matrix)) {
    cors[i] <- cor(gene_matrix[, i], longevity_data, use = "pairwise.complete.obs")
  }
  
  # Sort by absolute correlation
  abs_cors <- abs(cors)
  sorted_idx <- order(abs_cors, decreasing = TRUE)
  
  # Create dataframe of top correlated genes
  top_genes <- data.frame(
    Gene = names(cors)[sorted_idx[1:min(top_n, length(sorted_idx))]],
    Correlation = cors[sorted_idx[1:min(top_n, length(sorted_idx))]],
    stringsAsFactors = FALSE
  )
  
  # Add direction
  top_genes$Direction <- ifelse(top_genes$Correlation > 0, 
                                "Positive correlation with longevity", 
                                "Negative correlation with longevity")
  
  return(top_genes)
}

# Function to visualize gene correlations
visualize_gene_correlations <- function(correlation_results) {
  # Prepare data for plotting
  genes_df <- correlation_results
  genes_df$Gene <- factor(genes_df$Gene, levels = genes_df$Gene[order(genes_df$Correlation)])
  
  # Colors
  direction_colors <- c("Positive correlation with longevity" = "#9b383a", 
                        "Negative correlation with longevity" = "#86acb9")
  
  # Create plot
  p <- ggplot(genes_df, aes(x = Gene, y = Correlation, fill = Direction)) +
    geom_col() +
    scale_fill_manual(values = direction_colors) +
    coord_flip() +
    labs(
      title = "Genes most correlated with longevity",
      subtitle = paste("Based on", nrow(gene_data), "species"),
      x = "Gene",
      y = "Correlation with log10 lifespan"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 12)
    )
  
  return(p)
}

#####################################################
# 7. FUNCTION TO ANALYZE A SPECIES AGAINST ALL OTHERS
#####################################################

compare_species_against_all <- function(target_species, model, x_test_scaled, test_species, 
                                        var_importance, top_n = 20) {
  
  # Check if target species is in test set
  target_idx <- match(target_species, test_species)
  
  if (is.na(target_idx)) {
    stop(paste("Target species not found in test set:", target_species))
  }
  
  # Get feature values for target species
  target_features <- x_test_scaled[target_idx, , drop=FALSE]
  
  # Calculate average feature values for all other species
  other_indices <- setdiff(1:length(test_species), target_idx)
  other_features <- colMeans(x_test_scaled[other_indices, , drop=FALSE])
  
  # Calculate contributions for target species
  target_contributions <- target_features[1, ] * var_importance
  
  # Calculate contributions for average of other species
  other_contributions <- other_features * var_importance
  
  # Find differences in contributions
  contribution_diffs <- target_contributions - other_contributions
  
  # Sort by absolute difference
  sorted_idx <- order(abs(contribution_diffs), decreasing=TRUE)
  
  # Create results dataframe
  diff_genes <- data.frame(
    Gene = names(contribution_diffs)[sorted_idx[1:min(top_n, length(sorted_idx))]],
    Contribution_Diff = contribution_diffs[sorted_idx[1:min(top_n, length(sorted_idx))]],
    Target_Feature = target_features[1, sorted_idx[1:min(top_n, length(sorted_idx))]],
    Others_Feature = other_features[sorted_idx[1:min(top_n, length(sorted_idx))]],
    Coefficient = var_importance[sorted_idx[1:min(top_n, length(sorted_idx))]],
    stringsAsFactors = FALSE
  )
  
  # Add direction
  diff_genes$Direction <- ifelse(diff_genes$Contribution_Diff > 0, 
                                 "Higher in target species", 
                                 "Lower in target species")
  
  # Get predictions
  target_pred <- predict(model, newx=target_features, s=best_lambda)
  others_pred <- mean(predict(model, newx=x_test_scaled[other_indices, , drop=FALSE], s=best_lambda))
  
  # Create results list
  results <- list(
    target_species = target_species,
    target_prediction = target_pred,
    others_prediction = others_pred,
    prediction_diff = target_pred - others_pred,
    diff_genes = diff_genes
  )
  
  return(results)
}

# Function to visualize species vs all others
visualize_species_vs_all <- function(comparison_results) {
  # Extract data
  diff_genes <- comparison_results$diff_genes
  target_species <- comparison_results$target_species
  
  # Prepare data for plotting
  diff_genes$Gene <- factor(diff_genes$Gene, levels = diff_genes$Gene[order(diff_genes$Contribution_Diff)])
  
  # Direction colors
  direction_colors <- c("Higher in target species" = "#64a590", "Lower in target species" = "#31a6ad")
  
  # Create plot
  p <- ggplot(diff_genes, aes(x = Gene, y = Contribution_Diff, fill = Direction)) +
    geom_col() +
    scale_fill_manual(values = direction_colors) +
    coord_flip() +
    labs(
      title = paste(target_species, "vs Average of All Other Species"),
      subtitle = paste("Prediction difference:", 
                       round(comparison_results$prediction_diff, 3),
                       "| Target prediction:", 
                       round(comparison_results$target_prediction, 3)),
      x = "Gene",
      y = "Difference in Contribution to Longevity Prediction"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 12)
    )
  
  return(p)
}

#####################################################
# EXAMPLE CODE FOR PRIMATES VS RODENTS COMPARISON
#####################################################

# Choose species to analyze
species_of_interest <- c("Molossus_molossus", "Myotis_brandtii")

# Run the analysis
species_analysis <- analyze_species_predictions(
  species_list = species_of_interest,
  model = final_model,
  x_test_scaled = x_test_scaled, 
  test_species = test_species,
  var_importance = var_importance
)

# Generate reports
print_species_report(species_analysis, "Molossus_molossus")
print_comparison_report(species_analysis, "Molossus_molossus vs Myotis_brandtii")

# Visualize results
p1 <- visualize_species_genes(species_analysis, "Myotis_brandtii")
p2 <- visualize_species_comparison(species_analysis, "Molossus_molossus vs Myotis_brandtii")
print(p1)
print(p2)

ggsave("/Users/katiarenault/PhD/TOGA/revised_results/plots/myotis_brandtii_elastic_net_v_molmol.png", p1)
ggsave("/Users/katiarenault/PhD/TOGA/revised_results/plots/myotis_brandtii_v_molmol_elastic_net.png", p2)
