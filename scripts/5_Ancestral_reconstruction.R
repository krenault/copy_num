

# Gene Copy Number Ancestral Reconstruction using Marginal MCMC
# Based on discrete character evolution models - FIXED VERSION

# Load required libraries
if (!require("ape")) install.packages("ape")
if (!require("phytools")) install.packages("phytools")
if (!require("dplyr")) install.packages("dplyr")
if (!require("parallel")) install.packages("parallel")
if (!require("foreach")) install.packages("foreach")
if (!require("doParallel")) install.packages("doParallel")

library(ape)
library(phytools)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)

# Load your data
metadata <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")
gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber_Annotated.tsv",
                           row.names = "t_symbol", sep = '\t')

# Filter genes (remove those with >50% zeros)
gene_copy_data <- gene_copy_data[apply(gene_copy_data, 1, function(x) sum(x == 0, na.rm=TRUE)/sum(!is.na(x))) <= 0.5, ]

# Log transform longevity
metadata$maximum_longevity_y <- log(metadata$maximum_longevity_y)

# Load tree
tree <- read.tree("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_species_tree_revised.nwk")

# Function to prepare data for ancestral reconstruction
prepare_gene_data <- function(gene_copy_data, tree, gene_name) {
  # Extract copy numbers for the specific gene
  gene_copies <- gene_copy_data[gene_name, ]
  
  # Match species names between tree and data
  common_species <- intersect(tree$tip.label, names(gene_copies))
  
  if (length(common_species) < 10) {
    warning(paste("Too few species overlap for gene", gene_name, "- skipping"))
    return(NULL)
  }
  
  # Prune tree to match available data
  pruned_tree <- keep.tip(tree, common_species)
  
  # Extract copy numbers for pruned tree
  gene_vector <- gene_copies[pruned_tree$tip.label]
  gene_vector <- gene_vector[!is.na(gene_vector)]
  
  # Convert to factor for discrete analysis
  gene_factor <- as.factor(gene_vector)
  names(gene_factor) <- names(gene_vector)
  
  return(list(tree = pruned_tree, data = gene_factor, numeric_data = gene_vector))
}

# Function for ancestral state reconstruction using multiple methods
ancestral_reconstruction <- function(prepared_data, gene_name, n_simulations = 1000) {
  
  if (is.null(prepared_data)) return(NULL)
  
  tree <- prepared_data$tree
  gene_data <- prepared_data$data
  numeric_data <- prepared_data$numeric_data
  
  cat(paste("Processing gene:", gene_name, "\n"))
  cat(paste("Tree size:", Ntip(tree), "species\n"))
  cat(paste("Copy number range:", min(numeric_data), "-", max(numeric_data), "\n"))
  
  results <- list()
  results$gene_name <- gene_name
  results$n_species <- Ntip(tree)
  results$copy_range <- range(numeric_data)
  results$tree <- tree
  
  # 1. Maximum Likelihood Ancestral State Reconstruction
  cat("Running ML ancestral state reconstruction...\n")
  
  # Equal rates model (ER)
  fit_ER <- try(ace(gene_data, tree, model = "ER", type = "discrete"), silent = TRUE)
  if (!inherits(fit_ER, "try-error")) {
    results$ML_ER <- fit_ER
    results$ML_ER_loglik <- fit_ER$loglik
  }
  
  # All rates different model (ARD) - only if reasonable number of states
  if (length(levels(gene_data)) <= 6) {
    fit_ARD <- try(ace(gene_data, tree, model = "ARD", type = "discrete"), silent = TRUE)
    if (!inherits(fit_ARD, "try-error")) {
      results$ML_ARD <- fit_ARD
      results$ML_ARD_loglik <- fit_ARD$loglik
    }
  }
  
  # 2. Stochastic Character Mapping (MCMC approach)
  cat("Running stochastic character mapping...\n")
  
  # Use the better ML model for stochastic mapping
  if (exists("fit_ARD") && !inherits(fit_ARD, "try-error") && 
      exists("fit_ER") && !inherits(fit_ER, "try-error")) {
    best_fit <- if (fit_ARD$loglik > fit_ER$loglik) fit_ARD else fit_ER
    model_type <- if (fit_ARD$loglik > fit_ER$loglik) "ARD" else "ER"
  } else if (exists("fit_ER") && !inherits(fit_ER, "try-error")) {
    best_fit <- fit_ER
    model_type <- "ER"
  } else {
    cat("No successful ML fit - skipping stochastic mapping\n")
    return(results)
  }
  
  # Extract transition matrix from best fit
  Q_matrix <- best_fit$rates
  
  # Perform stochastic character mapping
  stoch_maps <- try(make.simmap(tree, gene_data, model = model_type, 
                                nsim = n_simulations, Q = Q_matrix), silent = TRUE)
  
  if (!inherits(stoch_maps, "try-error")) {
    results$stochastic_maps <- stoch_maps
    
    # Summarize stochastic maps
    summary_maps <- summary(stoch_maps)
    results$stoch_summary <- summary_maps
    
    # Extract marginal ancestral state probabilities
    results$marginal_ancestral_states <- summary_maps$ace
    
    # Calculate transition counts and rates
    results$transition_summary <- summary_maps$count
  }
  
  # 3. Model comparison using AIC
  if (exists("fit_ER") && exists("fit_ARD") && 
      !inherits(fit_ER, "try-error") && !inherits(fit_ARD, "try-error")) {
    
    # Calculate AIC for model comparison
    aic_ER <- -2 * fit_ER$loglik + 2 * length(fit_ER$rates)
    aic_ARD <- -2 * fit_ARD$loglik + 2 * length(fit_ARD$rates)
    
    results$model_comparison <- data.frame(
      Model = c("ER", "ARD"),
      LogLik = c(fit_ER$loglik, fit_ARD$loglik),
      AIC = c(aic_ER, aic_ARD),
      Delta_AIC = c(aic_ER - min(aic_ER, aic_ARD), aic_ARD - min(aic_ER, aic_ARD))
    )
  }
  
  return(results)
}

# Function to visualize results
plot_ancestral_reconstruction <- function(results, output_dir = "ancestral_plots") {
  
  if (is.null(results) || is.null(results$tree)) return()
  
  dir.create(output_dir, showWarnings = FALSE)
  
  gene_name <- results$gene_name
  tree <- results$tree
  
  # Create filename
  filename <- file.path(output_dir, paste0(gene_name, "_ancestral_reconstruction.pdf"))
  
  pdf(filename, width = 12, height = 8)
  
  # Set up multiple plots
  par(mfrow = c(2, 2))
  
  # 1. Plot tree with ML ancestral states (if available)
  if (!is.null(results$ML_ER)) {
    # Create color palette for copy numbers
    all_states <- levels(as.factor(c(as.character(results$ML_ER$lik.anc), 
                                     names(results$ML_ER$lik.anc))))
    cols <- rainbow(length(all_states))
    names(cols) <- all_states
    
    plot(tree, type = "phylogram", show.tip.label = TRUE, cex = 0.6)
    title(paste("ML Ancestral States -", gene_name))
    
    # Add node labels with ancestral state probabilities
    nodelabels(node = 1:tree$Nnode + Ntip(tree), 
               pie = results$ML_ER$lik.anc, 
               piecol = cols, cex = 0.5)
  }
  
  # 2. Plot stochastic mapping summary (if available)
  if (!is.null(results$stochastic_maps)) {
    plot(summary(results$stochastic_maps), fsize = 0.6, ftype = "i")
    title(paste("Stochastic Character Mapping -", gene_name))
  }
  
  # 3. Plot model comparison (if available)
  if (!is.null(results$model_comparison)) {
    barplot(results$model_comparison$AIC, 
            names.arg = results$model_comparison$Model,
            main = paste("Model Comparison -", gene_name),
            ylab = "AIC")
  }
  
  # 4. Plot copy number distribution
  if (!is.null(results$copy_range)) {
    hist(as.numeric(names(results$tree$tip.label)), 
         main = paste("Copy Number Distribution -", gene_name),
         xlab = "Copy Number", 
         breaks = seq(results$copy_range[1] - 0.5, results$copy_range[2] + 0.5, 1))
  }
  
  dev.off()
  
  cat(paste("Plot saved to:", filename, "\n"))
}

# FIXED: Function to run analysis on multiple genes in parallel
run_ancestral_analysis <- function(gene_copy_data, tree, selected_genes = NULL, 
                                   n_cores = 4, n_simulations = 1000) {
  
  # Use all genes if none specified
  if (is.null(selected_genes)) {
    selected_genes <- rownames(gene_copy_data)
  }
  
  # Limit to available genes
  selected_genes <- intersect(selected_genes, rownames(gene_copy_data))
  
  cat(paste("Running ancestral reconstruction for", length(selected_genes), "genes\n"))
  
  # OPTION 1: Use lapply instead of parallel processing (most reliable)
  cat("Running analysis sequentially for reliability...\n")
  
  results_list <- lapply(selected_genes, function(gene) {
    cat(paste("Processing gene:", gene, "\n"))
    prepared_data <- prepare_gene_data(gene_copy_data, tree, gene)
    ancestral_reconstruction(prepared_data, gene, n_simulations)
  })
  
  names(results_list) <- selected_genes
  
  # Remove NULL results
  results_list <- results_list[!sapply(results_list, is.null)]
  
  return(results_list)
  
  # OPTION 2: Fixed parallel version (uncomment if you want to use parallel processing)
  # cat("Setting up parallel processing...\n")
  # 
  # # Create a wrapper function that includes n_simulations as a parameter
  # process_gene_wrapper <- function(gene, gene_copy_data, tree, n_simulations) {
  #   prepared_data <- prepare_gene_data(gene_copy_data, tree, gene)
  #   ancestral_reconstruction(prepared_data, gene, n_simulations)
  # }
  # 
  # # Set up parallel processing
  # cl <- makeCluster(n_cores)
  # registerDoParallel(cl)
  # 
  # # Export ALL necessary objects and functions to workers
  # clusterExport(cl, c("gene_copy_data", "tree", "n_simulations",
  #                     "prepare_gene_data", "ancestral_reconstruction", 
  #                     "process_gene_wrapper"))
  # 
  # # Load required packages on workers
  # clusterEvalQ(cl, {
  #   library(ape)
  #   library(phytools)
  #   library(dplyr)
  # })
  # 
  # # Run analysis in parallel with explicit parameter passing
  # results_list <- foreach(gene = selected_genes, 
  #                         .packages = c("ape", "phytools", "dplyr"),
  #                         .export = c("gene_copy_data", "tree", "n_simulations",
  #                                     "prepare_gene_data", "ancestral_reconstruction")) %dopar% {
  #   process_gene_wrapper(gene, gene_copy_data, tree, n_simulations)
  # }
  # 
  # # Stop cluster
  # stopCluster(cl)
  # 
  # names(results_list) <- selected_genes
  # 
  # # Remove NULL results
  # results_list <- results_list[!sapply(results_list, is.null)]
  # 
  # return(results_list)
}

# Function to summarize results across genes
summarize_results <- function(results_list) {
  
  summary_df <- data.frame(
    Gene = character(),
    N_Species = integer(),
    Copy_Range_Min = integer(),
    Copy_Range_Max = integer(),
    ML_ER_LogLik = numeric(),
    ML_ARD_LogLik = numeric(),
    Best_Model = character(),
    N_Transitions = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (gene in names(results_list)) {
    result <- results_list[[gene]]
    
    if (is.null(result)) next
    
    best_model <- "None"
    if (!is.null(result$model_comparison)) {
      best_model <- result$model_comparison$Model[which.min(result$model_comparison$AIC)]
    } else if (!is.null(result$ML_ER)) {
      best_model <- "ER"
    }
    
    n_transitions <- NA
    if (!is.null(result$transition_summary)) {
      n_transitions <- sum(result$transition_summary)
    }
    
    summary_df <- rbind(summary_df, data.frame(
      Gene = gene,
      N_Species = result$n_species,
      Copy_Range_Min = result$copy_range[1],
      Copy_Range_Max = result$copy_range[2],
      ML_ER_LogLik = ifelse(is.null(result$ML_ER_loglik), NA, result$ML_ER_loglik),
      ML_ARD_LogLik = ifelse(is.null(result$ML_ARD_loglik), NA, result$ML_ARD_loglik),
      Best_Model = best_model,
      N_Transitions = n_transitions,
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_df)
}

# MAIN ANALYSIS EXECUTION
# ======================

# Test with a small subset first (remove this line to run on all genes)
test_genes <- c("TAF11L10", "NBPF12", "BPY2")

cat("Starting ancestral reconstruction analysis...\n")
cat("Tree has", Ntip(tree), "tips\n")
cat("Gene copy data has", nrow(gene_copy_data), "genes and", ncol(gene_copy_data), "species\n")

# Run the analysis (now using fixed function)
results <- run_ancestral_analysis(
  gene_copy_data = gene_copy_data,
  tree = tree,
  selected_genes = test_genes,  # Use all genes: selected_genes = NULL
  n_cores = min(4, detectCores() - 1),  # Use available cores (not used in sequential version)
  n_simulations = 500  # Increase for final analysis
)

# Summarize results
summary_table <- summarize_results(results)
print(summary_table)

# Generate plots for each gene
cat("Generating plots...\n")
for (gene in names(results)) {
  if (!is.null(results[[gene]])) {
    plot_ancestral_reconstruction(results[[gene]])
  }
}

# Save results
save(results, summary_table, file = "ancestral_reconstruction_results.RData")

cat("Analysis complete! Results saved to ancestral_reconstruction_results.RData\n")
cat("Plots saved in ancestral_plots/ directory\n")

# Example: Access results for a specific gene
if (length(results) > 0) {
  example_gene <- names(results)[1]
  cat(paste("\nExample results for", example_gene, ":\n"))
  
  if (!is.null(results[[example_gene]]$marginal_ancestral_states)) {
    cat("Marginal ancestral state probabilities (first 5 nodes):\n")
    print(head(results[[example_gene]]$marginal_ancestral_states, 5))
  }
  
  if (!is.null(results[[example_gene]]$model_comparison)) {
    cat("\nModel comparison:\n")
    print(results[[example_gene]]$model_comparison)
  }
}

# ALTERNATIVE: If you want to enable parallel processing later, uncomment the parallel
# version in the run_ancestral_analysis function and comment out the lapply version
###########################################################################################################################################################################
# Mk MODEL-BASED EVOLUTIONARY CHANGE ANALYSIS
# Uses your ancestral reconstruction pie chart data for sophisticated change analysis

library(ape)
library(phytools)
library(dplyr)

# Function to extract evolutionary changes from Mk model results
extract_mk_evolutionary_changes <- function(ancestral_results, gene_name) {
  
  result <- ancestral_results[[gene_name]]
  if (is.null(result)) return(NULL)
  
  # Get the best model (ARD preferred)
  best_model <- if (!is.null(result$ML_ARD)) result$ML_ARD else result$ML_ER
  if (is.null(best_model)) return(NULL)
  
  tree <- result$tree
  ancestral_probs <- best_model$lik.anc  # This is your pie chart data!
  
  # Get tip states (actual observed copy numbers)
  tip_states <- rep(NA, Ntip(tree))
  names(tip_states) <- tree$tip.label
  
  # Extract tip states from the original gene data that was used
  # Check what's available in the results structure
  if (!is.null(result$numeric_data)) {
    gene_data_subset <- result$numeric_data
    tip_states[names(gene_data_subset)] <- gene_data_subset
  } else {
    # Fallback: try to get from global gene_copy_data
    if (exists("gene_copy_data")) {
      gene_row <- gene_copy_data[gene_name, ]
      common_tips <- intersect(names(gene_row), tree$tip.label)
      tip_states[common_tips] <- as.numeric(gene_row[common_tips])
    } else {
      cat("  Warning: Cannot find tip state data for", gene_name, "\n")
      return(NULL)
    }
  }
  
  # Calculate expected ancestral states (weighted means from probabilities)
  state_names <- as.numeric(colnames(ancestral_probs))
  expected_ancestral <- apply(ancestral_probs, 1, function(prob_row) {
    sum(prob_row * state_names, na.rm = TRUE)
  })
  
  # Initialize empty changes data frame with correct structure
  changes <- data.frame(
    gene = character(0),
    node_parent = numeric(0),
    node_child = numeric(0),
    parent_state = numeric(0),
    child_state = numeric(0),
    change_magnitude = numeric(0),
    branch_length = numeric(0),
    change_rate = numeric(0),
    change_type = character(0),
    stringsAsFactors = FALSE
  )
  
  # Process each edge in the tree
  for (i in 1:nrow(tree$edge)) {
    parent_node <- tree$edge[i, 1]
    child_node <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get parent state
    if (parent_node <= Ntip(tree)) {
      # Parent is a tip (shouldn't happen in rooted tree, but just in case)
      parent_state <- tip_states[parent_node]
    } else {
      # Parent is internal node
      internal_index <- parent_node - Ntip(tree)
      parent_state <- expected_ancestral[internal_index]
    }
    
    # Get child state
    if (child_node <= Ntip(tree)) {
      # Child is a tip
      child_state <- tip_states[child_node]
    } else {
      # Child is internal node
      internal_index <- child_node - Ntip(tree)
      child_state <- expected_ancestral[internal_index]
    }
    
    # Skip if either state is missing
    if (is.na(parent_state) || is.na(child_state)) next
    
    # Calculate change
    change_magnitude <- child_state - parent_state
    change_rate <- change_magnitude / branch_length
    
    # Classify change type
    if (abs(change_magnitude) < 0.1) {
      change_type <- "No change"
    } else if (change_magnitude > 0) {
      change_type <- "Gene gain"
    } else {
      change_type <- "Gene loss"
    }
    
    # Add row to changes data frame
    new_row <- data.frame(
      gene = gene_name,
      node_parent = parent_node,
      node_child = child_node,
      parent_state = parent_state,
      child_state = child_state,
      change_magnitude = change_magnitude,
      branch_length = branch_length,
      change_rate = change_rate,
      change_type = change_type,
      stringsAsFactors = FALSE
    )
    changes <- rbind(changes, new_row)
  }
  
  return(list(
    changes = changes,
    ancestral_states = expected_ancestral,
    tip_states = tip_states,
    tree = tree
  ))
}

# Function to map trait evolution onto the same tree
map_trait_evolution <- function(tree, metadata, trait_column = "maximum_longevity_y") {
  
  # Find species column
  species_col <- find_species_column(metadata, tree)
  if (is.null(species_col)) {
    stop("Cannot match species between metadata and tree")
  }
  
  # Get trait data
  trait_data <- metadata[[trait_column]]
  names(trait_data) <- metadata[[species_col]]
  
  # Match to tree
  tip_traits <- trait_data[tree$tip.label]
  names(tip_traits) <- tree$tip.label
  
  # Reconstruct ancestral trait values using ML
  complete_cases <- !is.na(tip_traits)
  if (sum(complete_cases) < 10) {
    stop("Insufficient trait data")
  }
  
  tree_subset <- keep.tip(tree, tree$tip.label[complete_cases])
  trait_subset <- tip_traits[complete_cases][tree_subset$tip.label]
  
  # Ancestral trait reconstruction using fastAnc (ML under Brownian motion)
  ancestral_traits <- fastAnc(tree_subset, trait_subset)
  
  # Calculate trait changes along branches
  trait_changes <- data.frame(
    node_parent = numeric(),
    node_child = numeric(),
    parent_trait = numeric(),
    child_trait = numeric(),
    trait_change = numeric(),
    branch_length = numeric(),
    trait_change_rate = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Create combined tip + ancestral trait vector
  all_traits <- c(trait_subset, ancestral_traits)
  names(all_traits) <- c(tree_subset$tip.label, 
                         paste0("Node", (Ntip(tree_subset) + 1):(Ntip(tree_subset) + Nnode(tree_subset))))
  
  # Process each edge
  for (i in 1:nrow(tree_subset$edge)) {
    parent_node <- tree_subset$edge[i, 1]
    child_node <- tree_subset$edge[i, 2]
    branch_length <- tree_subset$edge.length[i]
    
    # Get trait values
    if (parent_node <= Ntip(tree_subset)) {
      parent_trait <- all_traits[tree_subset$tip.label[parent_node]]
    } else {
      parent_trait <- all_traits[paste0("Node", parent_node)]
    }
    
    if (child_node <= Ntip(tree_subset)) {
      child_trait <- all_traits[tree_subset$tip.label[child_node]]
    } else {
      child_trait <- all_traits[paste0("Node", child_node)]
    }
    
    trait_change <- child_trait - parent_trait
    trait_change_rate <- trait_change / branch_length
    
    trait_changes <- rbind(trait_changes, data.frame(
      node_parent = parent_node,
      node_child = child_node,
      parent_trait = parent_trait,
      child_trait = child_trait,
      trait_change = trait_change,
      branch_length = branch_length,
      trait_change_rate = trait_change_rate,
      stringsAsFactors = FALSE
    ))
  }
  
  return(list(
    trait_changes = trait_changes,
    ancestral_traits = ancestral_traits,
    tip_traits = trait_subset,
    tree = tree_subset
  ))
}

# Main function to analyze Mk-based evolutionary correlations
analyze_mk_evolutionary_correlations <- function(ancestral_results, metadata, 
                                                 trait_column = "maximum_longevity_y") {
  
  cat("=== Mk MODEL-BASED EVOLUTIONARY CORRELATION ANALYSIS ===\n")
  cat("Using sophisticated ancestral reconstruction data (your pie charts!)\n")
  cat("to analyze evolutionary correlations\n\n")
  
  correlation_results <- data.frame(
    Gene = character(),
    Branch_Change_Correlation = numeric(),
    Branch_Change_P_value = numeric(),
    Rate_Correlation = numeric(),
    Rate_P_value = numeric(),
    N_Branches = integer(),
    Interpretation = character(),
    Method = character(),
    stringsAsFactors = FALSE
  )
  
  for (gene_name in names(ancestral_results)) {
    
    cat("Processing gene:", gene_name, "\n")
    
    # Extract gene copy number changes using Mk model results
    gene_changes <- extract_mk_evolutionary_changes(ancestral_results, gene_name)
    if (is.null(gene_changes)) {
      cat("  Skipping: no Mk model results available\n")
      next
    }
    
    # Map trait evolution onto the same tree
    tryCatch({
      trait_evolution <- map_trait_evolution(gene_changes$tree, metadata, trait_column)
      
      # Match branches between gene and trait changes
      gene_branch_data <- gene_changes$changes
      trait_branch_data <- trait_evolution$trait_changes
      
      # Create matching keys for branches
      gene_branch_data$branch_key <- paste(gene_branch_data$node_parent, 
                                           gene_branch_data$node_child, sep = "_")
      trait_branch_data$branch_key <- paste(trait_branch_data$node_parent, 
                                            trait_branch_data$node_child, sep = "_")
      
      # Merge branch data
      merged_data <- merge(gene_branch_data, trait_branch_data, by = "branch_key", all = FALSE)
      
      if (nrow(merged_data) < 5) {
        cat("  Skipping: insufficient matching branches (", nrow(merged_data), ")\n")
        next
      }
      
      # Test correlation between evolutionary changes
      change_test <- cor.test(merged_data$change_magnitude, merged_data$trait_change)
      rate_test <- cor.test(merged_data$change_rate, merged_data$trait_change_rate)
      
      # Interpretation
      if (change_test$p.value < 0.05) {
        if (change_test$estimate > 0) {
          interpretation <- paste("Positive correlation: branches with copy number gains",
                                  "also show increases in", trait_column)
        } else {
          interpretation <- paste("Negative correlation: branches with copy number gains", 
                                  "show decreases in", trait_column)
        }
      } else {
        interpretation <- "No significant correlation between evolutionary changes"
      }
      
      correlation_results <- rbind(correlation_results, data.frame(
        Gene = gene_name,
        Branch_Change_Correlation = change_test$estimate,
        Branch_Change_P_value = change_test$p.value,
        Rate_Correlation = rate_test$estimate,
        Rate_P_value = rate_test$p.value,
        N_Branches = nrow(merged_data),
        Interpretation = interpretation,
        Method = "Mk model + ML trait reconstruction",
        stringsAsFactors = FALSE
      ))
      
      cat("  ✓ Analysis successful:", nrow(merged_data), "branches analyzed\n")
      
    }, error = function(e) {
      cat("  ✗ Analysis failed:", e$message, "\n")
    })
  }
  
  return(correlation_results)
}

# Enhanced visualization function
plot_mk_evolutionary_correlations <- function(correlation_results, ancestral_results, 
                                              top_genes = 3) {
  
  # Get top significant results
  sig_results <- correlation_results[!is.na(correlation_results$Branch_Change_P_value) & 
                                       correlation_results$Branch_Change_P_value < 0.05, ]
  
  if (nrow(sig_results) == 0) {
    cat("No significant correlations to plot\n")
    return()
  }
  
  sig_results <- sig_results[order(sig_results$Branch_Change_P_value), ]
  plot_genes <- head(sig_results$Gene, top_genes)
  
  pdf("mk_evolutionary_correlations.pdf", width = 15, height = 10)
  
  par(mfrow = c(2, 3))
  
  for (gene in plot_genes) {
    # Extract data for plotting
    gene_changes <- extract_mk_evolutionary_changes(ancestral_results, gene)
    
    if (!is.null(gene_changes)) {
      # Plot 1: Ancestral tree with copy number evolution
      tree <- gene_changes$tree
      plot(tree, show.tip.label = FALSE, main = paste("Copy Number Evolution:", gene))
      
      # Add ancestral states as node labels
      expected_states <- round(gene_changes$ancestral_states, 1)
      nodelabels(expected_states, cex = 0.8, bg = "lightblue")
      
      # Plot 2: Change magnitude vs branch length
      changes_data <- gene_changes$changes
      plot(changes_data$branch_length, abs(changes_data$change_magnitude),
           main = paste("Change Magnitude:", gene),
           xlab = "Branch Length", ylab = "Abs(Change in Copy Number)",
           pch = 19, col = ifelse(changes_data$change_magnitude > 0, "red", "blue"))
      legend("topright", c("Gains", "Losses"), col = c("red", "blue"), pch = 19)
    }
  }
  
  dev.off()
  
  cat("Mk evolutionary correlation plots saved to 'mk_evolutionary_correlations.pdf'\n")
}

# Function to integrate with your existing code
run_mk_change_analysis_with_your_data <- function() {
  
  cat("=== RUNNING Mk MODEL CHANGE ANALYSIS WITH YOUR DATA ===\n")
  
  # Check if your results are loaded
  if (!exists("results") || !exists("metadata") || !exists("tree")) {
    cat("Please ensure your data is loaded:\n")
    cat("- results: from ancestral reconstruction\n")
    cat("- metadata: your species metadata\n") 
    cat("- tree: your phylogenetic tree\n")
    return()
  }
  
  # Debug: Check the structure of your results
  cat("Checking results structure...\n")
  cat("Number of genes in results:", length(results), "\n")
  cat("First gene structure:\n")
  if (length(results) > 0) {
    first_gene <- results[[1]]
    cat("Available elements:", names(first_gene), "\n")
    if (!is.null(first_gene$tree)) {
      cat("Tree has", Ntip(first_gene$tree), "tips\n")
    }
  }
  
  # Run the analysis
  mk_correlations <- analyze_mk_evolutionary_correlations(
    ancestral_results = results,
    metadata = metadata,
    trait_column = "maximum_longevity_y"
  )
  
  # Only proceed if we have results
  if (nrow(mk_correlations) == 0) {
    cat("No results obtained. Check for errors above.\n")
    return(mk_correlations)
  }
  
  # Add multiple testing correction
  mk_correlations$Branch_Change_P_adjusted <- p.adjust(mk_correlations$Branch_Change_P_value, method = "fdr")
  mk_correlations$Rate_P_adjusted <- p.adjust(mk_correlations$Rate_P_value, method = "fdr")
  
  # Summary
  cat("\n=== Mk MODEL CORRELATION RESULTS ===\n")
  cat("Total genes analyzed:", nrow(mk_correlations), "\n")
  cat("Successful analyses:", sum(!is.na(mk_correlations$Branch_Change_P_value)), "\n")
  cat("Significant correlations (p < 0.05):", 
      sum(mk_correlations$Branch_Change_P_value < 0.05, na.rm = TRUE), "\n")
  cat("Significant after FDR correction:", 
      sum(mk_correlations$Branch_Change_P_adjusted < 0.05, na.rm = TRUE), "\n")
  
  # Show top results
  if (nrow(mk_correlations) > 0) {
    cat("\nTop results:\n")
    valid_results <- mk_correlations[!is.na(mk_correlations$Branch_Change_P_value), ]
    if (nrow(valid_results) > 0) {
      top_results <- head(valid_results[order(valid_results$Branch_Change_P_value), ], 5)
      print(top_results[, c("Gene", "Branch_Change_Correlation", "Branch_Change_P_value", 
                            "N_Branches", "Method")])
    }
  }
  
  # Create plots if we have significant results
  sig_results <- mk_correlations[!is.na(mk_correlations$Branch_Change_P_value) & 
                                   mk_correlations$Branch_Change_P_value < 0.05, ]
  if (nrow(sig_results) > 0) {
    plot_mk_evolutionary_correlations(mk_correlations, results)
  } else {
    cat("No significant results to plot.\n")
  }
  
  return(mk_correlations)
}

# Helper function (keep the same as before)
find_species_column <- function(metadata, tree) {
  possible_cols <- c("Scientific_name", "species_name", "species", "Species", 
                     "scientific_name", "tip_label", "name")
  
  for (col in possible_cols) {
    if (col %in% colnames(metadata)) {
      overlap <- length(intersect(metadata[[col]], tree$tip.label))
      if (overlap > 0.5 * min(length(metadata[[col]]), length(tree$tip.label))) {
        return(col)
      }
    }
  }
  return(NULL)
}

cat("=== Mk MODEL EVOLUTIONARY CHANGE ANALYSIS LOADED ===\n")
cat("This uses your sophisticated ancestral reconstruction results!\n")
cat("\nTo run with your data:\n")
cat("mk_results <- run_mk_change_analysis_with_your_data()\n")
cat("\nThis approach is superior to PIC because it:\n")
cat("1. Uses full probabilistic ancestral states (your pie charts)\n")
cat("2. Accounts for uncertainty in reconstructions\n")
cat("3. Analyzes changes along specific branches\n")
cat("4. Uses the complex evolutionary models you already fitted\n")

mk_results <- run_mk_change_analysis_with_your_data()
