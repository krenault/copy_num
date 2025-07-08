########################################### 
# SIMPLE PHYLOGENETIC RANKING - ROBUST ALTERNATIVE
# This approach will analyze ALL your genes with phylogenetic correction
###########################################

cat("Starting SIMPLE phylogenetic ranking analysis...\n")

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ape)
library(fgsea)

# Set memory options
if(.Platform$OS.type == "windows") {
  memory.limit(size = 16000)
}
options(expressions = 500000)
gc()

# Read data (same as before)
cat("Reading input files...\n")
gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber_Annotated.tsv",
                           row.names = "t_symbol", sep = '\t')

# For testing: Use subset (remove this line for full analysis)
# IMPORTANT: Uncomment the next line to use only 1000 genes for testing
#gene_copy_data <- head(gene_copy_data, 1000)
gene_copy_data <- gene_copy_data %>% dplyr::select(-t_gene)

metadata <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")

# Data cleaning and preparation (same as before)
cat("Cleaning species names and aggregating individuals...\n")
clean_species_names <- function(names) {
  str_replace(names, "_\\d+$", "")
}

species_cols <- colnames(gene_copy_data)
cleaned_cols <- clean_species_names(species_cols)
colname_mapping <- data.frame(
  original = species_cols,
  cleaned = cleaned_cols,
  stringsAsFactors = FALSE
)

# Memory-efficient aggregation
unique_species <- unique(colname_mapping$cleaned)
unique_genes <- rownames(gene_copy_data)

aggregated_matrix <- matrix(NA, nrow = length(unique_genes), ncol = length(unique_species))
rownames(aggregated_matrix) <- unique_genes
colnames(aggregated_matrix) <- unique_species

for(species in unique_species) {
  original_cols <- colname_mapping$original[colname_mapping$cleaned == species]
  original_cols <- original_cols[original_cols %in% colnames(gene_copy_data)]
  
  if(length(original_cols) > 0) {
    if(length(original_cols) == 1) {
      aggregated_matrix[, species] <- gene_copy_data[, original_cols]
    } else {
      species_data <- gene_copy_data[, original_cols, drop = FALSE]
      aggregated_matrix[, species] <- rowMeans(species_data, na.rm = TRUE)
    }
  }
}

gene_copy_data_aggregated <- as.data.frame(aggregated_matrix)
rm(aggregated_matrix, gene_copy_data)
gc()

# Process metadata
cat("Processing metadata...\n")
metadata_processed <- metadata %>%
  dplyr::select(Scientific_name, MLres, order) %>%
  distinct() %>%
  filter(!is.na(MLres)) %>%
  mutate(Scientific_name = str_replace_all(Scientific_name, " ", "_"))

common_species <- intersect(colnames(gene_copy_data_aggregated),
                            metadata_processed$Scientific_name)
gene_copy_data <- gene_copy_data_aggregated[, common_species]
rm(gene_copy_data_aggregated)
gc()

cat(sprintf("Analyzing %d species that appear in both datasets\n", length(common_species)))

# Load phylogenetic tree
cat("Loading phylogenetic tree...\n")
phylo_tree <- tryCatch({
  read.tree("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_species_tree_revised.nwk")
}, error = function(e) {
  cat("WARNING: Could not load tree file. Using equal weights.\n")
  NULL
})

# Identify outliers
identify_outliers <- function(lifespans, sd_threshold = 1.5) {
  mean_ls <- mean(lifespans, na.rm = TRUE)
  sd_ls <- sd(lifespans, na.rm = TRUE)
  threshold <- mean_ls + (sd_threshold * sd_ls)
  lifespans > threshold
}

cat("Identifying outliers...\n")
global_outliers <- metadata_processed %>%
  mutate(is_global_outlier = identify_outliers(MLres)) %>%
  filter(is_global_outlier)

cat("DEBUG: Number of global outliers identified:", nrow(global_outliers), "\n")

outlier_species <- global_outliers$Scientific_name
normal_species <- setdiff(metadata_processed$Scientific_name, outlier_species)

cat("Outlier species:", length(outlier_species), "\n")
cat("Normal species:", length(normal_species), "\n")

# Gene effects calculation (same as before)
calculate_geneffects_gains <- function(species_data) {
  copy_numbers <- as.numeric(species_data[, 1])
  names(copy_numbers) <- rownames(species_data)
  
  copy_numbers[copy_numbers < 1] <- NA
  
  diff_from_baseline <- copy_numbers - 1
  absolute_diff <- abs(diff_from_baseline)
  
  ranked_genes <- rank(-absolute_diff, ties.method = "min", na.last = "keep")
  
  return(list(
    genes = names(copy_numbers),
    rank = ranked_genes,
    direction = ifelse(diff_from_baseline > 0, "gain", NA),
    magnitude = diff_from_baseline
  ))
}

# Process species data (same as before but streamlined)
cat("Calculating gene effects across species...\n")
all_species_results <- list()
species_counter <- 1

for(species in common_species) {
  species_data <- gene_copy_data[, species, drop = FALSE]
  
  if(!all(is.na(species_data))) {
    effects <- calculate_geneffects_gains(species_data)
    
    if(!is.null(effects)) {
      all_species_results[[species_counter]] <- data.frame(
        Scientific_name = species,
        Gene = effects$genes,
        Rank = effects$rank,
        Direction = effects$direction,
        Magnitude = effects$magnitude,
        stringsAsFactors = FALSE
      )
    }
  }
  species_counter <- species_counter + 1
}

# Combine results
all_species_results <- all_species_results[!sapply(all_species_results, is.null)]
gene_effect_df <- bind_rows(all_species_results)
rm(all_species_results)
gc()

cat("DEBUG: Gene effect data dimensions:", dim(gene_effect_df), "\n")

############################################ 
# SIMPLE PHYLOGENETIC RANKING METHOD
###########################################

cat("Implementing SIMPLE phylogenetic ranking...\n")

# Function to calculate phylogenetic diversity score
calculate_phylo_diversity <- function(species_list, phylo_tree) {
  if(is.null(phylo_tree) || length(species_list) < 2) {
    return(1)  # Default weight
  }
  
  common_species <- intersect(species_list, phylo_tree$tip.label)
  if(length(common_species) < 2) {
    return(1)
  }
  
  tryCatch({
    tree_subset <- keep.tip(phylo_tree, common_species)
    distances <- cophenetic.phylo(tree_subset)
    
    # Return mean pairwise distance as diversity measure
    mean(distances, na.rm = TRUE)
  }, error = function(e) {
    return(1)
  })
}

# Simple phylogenetic ranking function - MEMORY OPTIMIZED
simple_phylo_ranking <- function(gene_effect_df, phylo_tree, outlier_species, normal_species) {
  
  cat("Processing genes for simple phylogenetic ranking...\n")
  
  # Get unique genes
  unique_genes <- unique(gene_effect_df$Gene)
  total_genes <- length(unique_genes)
  
  cat("Total genes to process:", total_genes, "\n")
  
  # MEMORY CHECK: If too many genes, process in batches
  if(total_genes > 5000) {
    cat("Large dataset detected. Processing in batches of 1000 genes.\n")
    batch_size <- 1000
  } else {
    batch_size <- total_genes
  }
  
  # Storage for results
  all_results <- list()
  
  # Process in batches
  for(batch_start in seq(1, total_genes, batch_size)) {
    batch_end <- min(batch_start + batch_size - 1, total_genes)
    batch_genes <- unique_genes[batch_start:batch_end]
    
    cat(sprintf("Processing batch %d-%d of %d genes\n", batch_start, batch_end, total_genes))
    
    # Process current batch
    batch_results <- list()
    
    for(gene_idx in 1:length(batch_genes)) {
      gene <- batch_genes[gene_idx]
      
      if(gene_idx %% 100 == 0) {
        cat(sprintf("  Batch gene %d of %d: %s\n", gene_idx, length(batch_genes), gene))
      }
      
      # Use base R filtering to avoid dplyr memory issues
      gene_indices <- which(gene_effect_df$Gene == gene & 
                              gene_effect_df$Direction == "gain" & 
                              !is.na(gene_effect_df$Rank))
      
      if(length(gene_indices) < 2) {
        next  # Skip genes with insufficient data
      }
      
      # Extract gene data using base R
      gene_data <- gene_effect_df[gene_indices, ]
      
      # Separate outliers and normals using base R
      outlier_indices <- gene_data$Scientific_name %in% outlier_species
      normal_indices <- gene_data$Scientific_name %in% normal_species
      
      outlier_data <- gene_data[outlier_indices, ]
      normal_data <- gene_data[normal_indices, ]
      
      # Need at least TWO of each for reliable t-tests
      if(nrow(outlier_data) < 2 || nrow(normal_data) < 2) {
        next
      }
      
      # Calculate basic statistics
      outlier_mean_rank <- mean(outlier_data$Rank, na.rm = TRUE)
      normal_mean_rank <- mean(normal_data$Rank, na.rm = TRUE)
      outlier_mean_magnitude <- mean(outlier_data$Magnitude, na.rm = TRUE)
      normal_mean_magnitude <- mean(normal_data$Magnitude, na.rm = TRUE)
      
      # Basic effect size
      effect_size <- outlier_mean_rank - normal_mean_rank
      magnitude_effect <- outlier_mean_magnitude - normal_mean_magnitude
      
      # Phylogenetic diversity weighting
      all_species_in_gene <- gene_data$Scientific_name
      phylo_diversity <- calculate_phylo_diversity(all_species_in_gene, phylo_tree)
      
      # Phylogenetically-weighted effect size
      phylo_weighted_effect <- effect_size * phylo_diversity
      phylo_weighted_magnitude <- magnitude_effect * phylo_diversity
      
      # Simple t-test for p-value
      t_test_result <- tryCatch({
        t.test(outlier_data$Rank, normal_data$Rank)
      }, error = function(e) NULL)
      
      # Store results
      batch_results[[gene]] <- list(
        Gene = gene,
        n_outliers = nrow(outlier_data),
        n_normals = nrow(normal_data),
        outlier_mean_rank = outlier_mean_rank,
        normal_mean_rank = normal_mean_rank,
        outlier_mean_magnitude = outlier_mean_magnitude,
        normal_mean_magnitude = normal_mean_magnitude,
        standard_effect_size = effect_size,
        magnitude_effect_size = magnitude_effect,
        phylo_diversity_score = phylo_diversity,
        phylo_weighted_effect = phylo_weighted_effect,
        phylo_weighted_magnitude = phylo_weighted_magnitude,
        p_value = if(!is.null(t_test_result)) t_test_result$p.value else NA,
        species_in_analysis = length(all_species_in_gene)
      )
    }
    
    # Add batch results to main results
    all_results <- c(all_results, batch_results)
    
    # Force garbage collection after each batch
    rm(batch_results)
    gc()
    
    cat(sprintf("Completed batch. Total results so far: %d\n", length(all_results)))
  }
  
  return(all_results)
}

# Run the simple phylogenetic ranking
cat("Running simple phylogenetic ranking analysis...\n")
simple_results <- simple_phylo_ranking(gene_effect_df, phylo_tree, outlier_species, normal_species)

cat("Successfully analyzed", length(simple_results), "genes\n")

# Convert results to data frame
cat("Converting results to data frame...\n")
simple_results_df <- do.call(rbind, lapply(simple_results, function(x) {
  data.frame(
    Gene = x$Gene,
    n_outliers = x$n_outliers,
    n_normals = x$n_normals,
    outlier_mean_rank = x$outlier_mean_rank,
    normal_mean_rank = x$normal_mean_rank,
    outlier_mean_magnitude = x$outlier_mean_magnitude,
    normal_mean_magnitude = x$normal_mean_magnitude,
    standard_effect_size = x$standard_effect_size,
    magnitude_effect_size = x$magnitude_effect_size,
    phylo_diversity_score = x$phylo_diversity_score,
    phylo_weighted_effect = x$phylo_weighted_effect,
    phylo_weighted_magnitude = x$phylo_weighted_magnitude,
    p_value = x$p_value,
    species_in_analysis = x$species_in_analysis,
    stringsAsFactors = FALSE
  )
}))

# Add FDR correction
simple_results_df$adj_p_value <- p.adjust(simple_results_df$p_value, method = "BH")

# Sort by phylogenetically weighted effect
simple_results_df <- simple_results_df[order(simple_results_df$phylo_weighted_effect, decreasing = TRUE), ]

cat("Final results dimensions:", dim(simple_results_df), "\n")

# CREATE MULTIPLE RANKING STRATEGIES - FIXED
###########################################

cat("Creating multiple ranking strategies...\n")

# Create ranking vectors for GSEA - with better data cleaning
rankings <- list()

# Strategy 1: Standard effect size
valid_standard <- !is.na(simple_results_df$standard_effect_size) & 
  is.finite(simple_results_df$standard_effect_size)
if(sum(valid_standard) > 10) {
  rankings[["standard_effect"]] <- setNames(simple_results_df$standard_effect_size[valid_standard], 
                                            simple_results_df$Gene[valid_standard])
  # Remove any remaining invalid values
  rankings[["standard_effect"]] <- rankings[["standard_effect"]][is.finite(rankings[["standard_effect"]])]
}

# Strategy 2: Phylogenetically weighted effect
valid_phylo <- !is.na(simple_results_df$phylo_weighted_effect) & 
  is.finite(simple_results_df$phylo_weighted_effect)
if(sum(valid_phylo) > 10) {
  rankings[["phylo_weighted_effect"]] <- setNames(simple_results_df$phylo_weighted_effect[valid_phylo], 
                                                  simple_results_df$Gene[valid_phylo])
  # Remove any remaining invalid values
  rankings[["phylo_weighted_effect"]] <- rankings[["phylo_weighted_effect"]][is.finite(rankings[["phylo_weighted_effect"]])]
}

# Strategy 3: Magnitude-based phylogenetic weighting
valid_magnitude <- !is.na(simple_results_df$phylo_weighted_magnitude) & 
  is.finite(simple_results_df$phylo_weighted_magnitude)
if(sum(valid_magnitude) > 10) {
  rankings[["phylo_weighted_magnitude"]] <- setNames(simple_results_df$phylo_weighted_magnitude[valid_magnitude], 
                                                     simple_results_df$Gene[valid_magnitude])
  # Remove any remaining invalid values
  rankings[["phylo_weighted_magnitude"]] <- rankings[["phylo_weighted_magnitude"]][is.finite(rankings[["phylo_weighted_magnitude"]])]
}

# Strategy 4: -log10(p) * sign(effect) - more careful filtering
valid_p_combo <- !is.na(simple_results_df$p_value) & 
  !is.na(simple_results_df$phylo_weighted_effect) &
  simple_results_df$p_value > 0 & 
  simple_results_df$p_value <= 1 &
  is.finite(simple_results_df$phylo_weighted_effect)

if(sum(valid_p_combo) > 10) {
  log_p_values <- -log10(simple_results_df$p_value[valid_p_combo] + 1e-300)
  effect_signs <- sign(simple_results_df$phylo_weighted_effect[valid_p_combo])
  log_p_signed <- log_p_values * effect_signs
  
  # Remove any invalid values
  log_p_signed <- log_p_signed[is.finite(log_p_signed)]
  
  if(length(log_p_signed) > 10) {
    rankings[["log_p_signed_phylo"]] <- setNames(log_p_signed, 
                                                 simple_results_df$Gene[valid_p_combo][is.finite(log_p_signed)])
  }
}

cat("Created", length(rankings), "ranking strategies:\n")
for(name in names(rankings)) {
  valid_count <- sum(is.finite(rankings[[name]]))
  cat("  -", name, ":", valid_count, "valid genes\n")
  
  # Show summary statistics
  if(valid_count > 0) {
    cat("    Range:", round(min(rankings[[name]], na.rm=TRUE), 3), "to", 
        round(max(rankings[[name]], na.rm=TRUE), 3), "\n")
  }
}

############################################ 
# ENHANCED GSEA WITH SIMPLE PHYLOGENETIC RANKINGS
###########################################

# Enhanced GSEA function - FIXED for data issues
run_enhanced_gsea <- function(rankings_list, pathway_files, output_dir, 
                              nperm = 1000, min_size = 15, max_size = 500) {  # Reduced nperm for speed
  
  gsea_results <- list()
  
  for(ranking_name in names(rankings_list)) {
    cat("Running GSEA with", ranking_name, "ranking...\n")
    
    stats_vector <- rankings_list[[ranking_name]]
    
    # Thorough data cleaning
    stats_vector <- stats_vector[!is.na(stats_vector) & is.finite(stats_vector)]
    
    # Remove duplicate gene names (take first occurrence)
    if(any(duplicated(names(stats_vector)))) {
      cat("  Removing", sum(duplicated(names(stats_vector))), "duplicate genes\n")
      stats_vector <- stats_vector[!duplicated(names(stats_vector))]
    }
    
    # Check minimum genes
    if(length(stats_vector) < 10) {
      cat("  Skipping", ranking_name, "- insufficient genes (", length(stats_vector), ")\n")
      next
    }
    
    # Sort in descending order
    stats_vector <- sort(stats_vector, decreasing = TRUE)
    
    cat("  Using", length(stats_vector), "genes for GSEA\n")
    cat("  Score range:", round(min(stats_vector), 3), "to", round(max(stats_vector), 3), "\n")
    
    tryCatch({
      # Run GSEA - use fgseaMultilevel (recommended)
      fgsea_results <- fgsea(
        pathways = pathway_files,
        stats = stats_vector,
        minSize = min_size,
        maxSize = max_size  # Removed nperm to use fgseaMultilevel
      )
      
      # Check if results are valid
      if(is.null(fgsea_results) || nrow(fgsea_results) == 0) {
        cat("  Warning: No pathways tested for", ranking_name, "\n")
        next
      }
      
      # Convert to data frame and clean
      fgsea_df <- data.frame(fgsea_results)
      
      # Convert list columns to strings safely
      list_cols <- sapply(fgsea_df, is.list)
      for (col in names(fgsea_df)[list_cols]) {
        fgsea_df[[col]] <- sapply(fgsea_df[[col]], function(x) {
          if(is.null(x) || length(x) == 0) return("")
          paste(x, collapse=";")
        })
      }
      
      # Add ranking method
      fgsea_df$ranking_method <- ranking_name
      
      # Save individual results
      output_file <- file.path(output_dir, paste0("simple_phylo_fgsea_", ranking_name, "_results.csv"))
      write.csv(fgsea_df, output_file, row.names = FALSE)
      
      gsea_results[[ranking_name]] <- fgsea_df
      
      # Report significant pathways
      sig_pathways <- sum(fgsea_df$padj < 0.05, na.rm = TRUE)
      cat("  ", ranking_name, ": ", sig_pathways, " significant pathways (FDR < 0.05)\n")
      
      if(sig_pathways > 0) {
        top_pathway <- fgsea_df[which.min(fgsea_df$padj), ]
        cat("    Top pathway:", top_pathway$pathway, "(padj =", 
            format(top_pathway$padj, scientific = TRUE), ")\n")
      }
      
    }, error = function(e) {
      cat("  Error with", ranking_name, ":", e$message, "\n")
      cat("  Stats vector summary: length =", length(stats_vector), 
          ", range =", paste(range(stats_vector, na.rm=TRUE), collapse=" to "), "\n")
    })
  }
  
  return(gsea_results)
}

# Save results
write_output <- function(data, filename) {
  if(is.data.frame(data) && nrow(data) > 0) {
    write.csv(data, filename, row.names = FALSE)
    cat("SUCCESS: Wrote", nrow(data), "results to", filename, "\n")
    return(TRUE)
  } else {
    cat("WARNING: No results to write to", filename, "\n")
    return(FALSE)
  }
}

# Save main results
results_written <- write_output(
  simple_results_df,
  "/Users/katiarenault/PhD/TOGA/revised_results/results/simple_phylo_ranking_results_w_or.csv"
)
# 
# # Run GSEA if we have results
# if(results_written && length(rankings) > 0) {
#   cat("Proceeding with GSEA analysis...\n")
#   
#   # Load pathways
#   cat("Loading pathway databases...\n")
#   pathway_files <- list()
#   
#   tryCatch({
#     pathways.hallmark <- gmtPathways("/Users/katiarenault/PhD/GSEA/h.all.v2023.2.Hs.symbols.gmt")
#     pathway_files <- c(pathway_files, pathways.hallmark)
#     cat("  Loaded Hallmark pathways:", length(pathways.hallmark), "\n")
#   }, error = function(e) cat("  Could not load Hallmark pathways\n"))
#   
#   tryCatch({
#     pathways.kegg <- gmtPathways("/Users/katiarenault/PhD/GSEA/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
#     pathway_files <- c(pathway_files, pathways.kegg)
#     cat("  Loaded KEGG pathways:", length(pathways.kegg), "\n")
#   }, error = function(e) cat("  Could not load KEGG pathways\n"))
#   
#   # Add these pathway databases for better coverage:
#   tryCatch({
#     pathways.go_bp <- gmtPathways("/Users/katiarenault/PhD/GSEA/c5.go.bp.v2023.2.Hs.symbols.gmt")
#     pathway_files <- c(pathway_files, pathways.go_bp)
#     cat("  Loaded GO Biological Process pathways:", length(pathways.go_bp), "\n")
#   }, error = function(e) cat("  Could not load GO BP pathways\n"))
#   
#   tryCatch({
#     pathways.reactome <- gmtPathways("/Users/katiarenault/PhD/GSEA/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
#     pathway_files <- c(pathway_files, pathways.reactome)
#     cat("  Loaded Reactome pathways:", length(pathways.reactome), "\n")
#   }, error = function(e) cat("  Could not load Reactome pathways\n"))
#   
#   if(length(pathway_files) > 0) {
#     cat("Total pathways loaded:", length(pathway_files), "\n")
#     
#     # Run GSEA
#     output_dir <- "/Users/katiarenault/PhD/TOGA/revised_results/results/"
#     gsea_results <- run_enhanced_gsea(rankings, pathway_files, output_dir, nperm = 10000)
#     
#     # Summary
#     cat("\n=== SIMPLE PHYLOGENETIC GSEA SUMMARY ===\n")
#     for(ranking_name in names(gsea_results)) {
#       result <- gsea_results[[ranking_name]]
#       sig_count <- sum(result$padj < 0.05, na.rm = TRUE)
#       cat(ranking_name, ": ", sig_count, " significant pathways\n")
#       
#       if(sig_count > 0) {
#         top_3 <- result[order(result$padj)[1:min(3, sig_count)], ]
#         for(i in 1:nrow(top_3)) {
#           cat("  ", i, ". ", top_3$pathway[i], " (padj = ", 
#               format(top_3$padj[i], scientific = TRUE), ")\n", sep = "")
#         }
#       }
#     }
#   }
# }

############################################ 
# ANALYSIS SUMMARY
###########################################

cat("\n=== SIMPLE PHYLOGENETIC RANKING SUMMARY ===\n")
cat("Total genes analyzed:", nrow(simple_results_df), "\n")
cat("Genes with phylogenetic weighting:", sum(!is.na(simple_results_df$phylo_weighted_effect)), "\n")
cat("Genes with p-values:", sum(!is.na(simple_results_df$p_value)), "\n")

if(nrow(simple_results_df) > 0) {
  # Effect size comparison
  valid_both <- !is.na(simple_results_df$standard_effect_size) & !is.na(simple_results_df$phylo_weighted_effect)
  if(sum(valid_both) > 1) {
    correlation <- cor(simple_results_df$standard_effect_size[valid_both], 
                       simple_results_df$phylo_weighted_effect[valid_both])
    cat("Correlation between standard and phylogenetic effects:", round(correlation, 3), "\n")
  }
  
  # Significance
  sig_genes <- sum(simple_results_df$adj_p_value < 0.05, na.rm = TRUE)
  cat("Significant genes (FDR < 0.05):", sig_genes, "\n")
  
  if(sig_genes > 0) {
    top_genes <- simple_results_df[simple_results_df$adj_p_value < 0.05 & !is.na(simple_results_df$adj_p_value), ]
    top_genes <- head(top_genes[order(top_genes$adj_p_value), ], 5)
    cat("Top significant genes:\n")
    for(i in 1:nrow(top_genes)) {
      cat("  ", top_genes$Gene[i], " (adj.p =", format(top_genes$adj_p_value[i], scientific = TRUE), ")\n")
    }
  }
  
  # Phylogenetic diversity
  phylo_scores <- simple_results_df$phylo_diversity_score[!is.na(simple_results_df$phylo_diversity_score)]
  if(length(phylo_scores) > 0) {
    cat("Phylogenetic diversity scores:\n")
    cat("  Mean:", round(mean(phylo_scores), 3), "\n")
    cat("  Range:", round(min(phylo_scores), 3), "to", round(max(phylo_scores), 3), "\n")
  }
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Enhanced analysis with improved phylogenetic weighting complete!\n")
cat("Key advantages of this approach:\n")
cat("  ✓ Analyzes genes with ≥2 species per group (outliers AND normals)\n")
cat("  ✓ Consistent statistical testing (t-tests) for all analyzed genes\n")
cat("  ✓ Simple and robust phylogenetic correction\n")
cat("  ✓ Multiple ranking strategies for GSEA\n")
cat("  ✓ Works well with count data\n")
cat("  ✓ No incomplete p-value results\n")
cat("\nFiles saved:\n")
cat("  - Main results: simple_phylo_ranking_results.csv\n")
cat("  - GSEA results: simple_phylo_fgsea_*_results.csv\n")
cat("\nDataframes available:\n")
cat("  - simple_results_df: Main phylogenetic ranking results\n")
cat("  - rankings: Multiple ranking strategies\n")
cat("  - gsea_results: GSEA results for each ranking\n")

###########
## GSEA ###
###########

cors_pathways <- read.csv('/Users/katiarenault/PhD/TOGA/revised_results/results/simple_phylo_ranking_results_w_or.csv')
pathways.reactome <- gmtPathways("/Users/katiarenault/Desktop/PhD/Human_GSEA/h.all.v2023.2.Hs.symbols.gmt")
pathways.hallmark <- gmtPathways("/Users/katiarenault/Desktop/PhD/Human_GSEA/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
pathways.kegg <- gmtPathways("/Users/katiarenault/Desktop/PhD/Human_GSEA/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
#pathways.go <- gmtPathways("/Nori_1/krenault/hibernator_rer_converge/data/c5.go.v2023.2.Hs.symbols.gmt")
pathways.hallmark <- c(pathways.kegg, pathways.hallmark, pathways.reactome)

cors_pathways$phylo_weighted_effect <- as.numeric(cors_pathways$phylo_weighted_effect)
cors_pathways <- cors_pathways %>% tidyr::drop_na()
positive_significant_genes <- cors_pathways %>%
  # filter(Rho < 0) %>%
  as_tibble() %>%
  arrange(adj_p_value)

positive_stats_vector <- setNames(positive_significant_genes$phylo_weighted_effect, positive_significant_genes$Gene)
positive_fgsea_results <- fgsea(pathways = pathways.hallmark, stats = positive_stats_vector)
positive_fgsea_results <- data.frame(positive_fgsea_results)
positive_fgsea_results <- positive_fgsea_results %>% select(-leadingEdge)
write.csv(positive_fgsea_results, '/Users/katiarenault/PhD/TOGA/revised_results/results/simple_phylo_ranking_results_pathways_w_or.csv')
########################
## Plotting results ###
########################

# Load libraries
library(ggplot2)
library(dplyr)

# Read the data
file_path <- "/Users/katiarenault/PhD/TOGA/revised_results/results/simple_phylo_ranking_results_w_or.csv"
gene_data <- read.csv(file_path)

# Highlight significant genes (adj_p_value < 0.05) and scale size by n_outliers
p <- ggplot(simple_results_df, aes(x = normal_mean_rank, y = outlier_mean_rank)) +
  # Non-significant genes (faint, small)
  geom_point(data = subset(simple_results_df, adj_p_value >= 0.05),
             alpha = 0.3, color = "grey70", size = 1.5) +
  # Significant genes (colored by p-value, sized by n_outliers)
  geom_point(data = subset(simple_results_df, adj_p_value < 0.05),
             aes(color = -log10(adj_p_value), size = n_outliers),
             alpha = 0.7) +
  # Diagonal reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  scale_color_gradient(low = "#f1baaa", high = "#9b383a", name = "-log10(adj p-value)") +
  # Size scale for n_outliers (adjust range as needed)
  scale_size_continuous(range = c(3, 8), name = "Number of outliers") +
  # Labels and theme
  labs(
    title = "Outlier vs. normal mean ranks",
  #  subtitle = "Significant genes (adj p < 0.05) colored by significance and sized by outlier count",
    x = "Normal Mean Rank",
    y = "Outlier Mean Rank"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9))

##################
## Volcano data ##
##################

library(ggplot2)
library(ggrepel)
library(dplyr)

# Prepare the data with similar logic to your example
volcano_data <- gene_data %>%
  mutate(
    # Create significance categories
    significance_level = case_when(
      adj_p_value < 0.05 | abs(standard_effect_size) >= 2.2 ~ "FDR_significant",
      p_value < 0.05 ~ "p_significant",
      TRUE ~ "not_significant"
    ),
    # Create color values based on effect size
    point_color = case_when(
      significance_level == "FDR_significant" ~ standard_effect_size,
      significance_level == "p_significant" ~ standard_effect_size * 0.5,
      TRUE ~ NA_real_
    ),
    log_p = -log10(p_value)
  ) %>%
  # Rank by significance (combination of effect size and p-value)
  mutate(significance_rank = rank(-(abs(standard_effect_size) * log_p))) %>%
  # Flag top genes for labeling (only from FDR significant points)
  mutate(to_label = ifelse(significance_level == "FDR_significant" & significance_rank <= 15, Gene, NA))

# Create the volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = standard_effect_size, y = log_p)) +
  # Non-significant points
  geom_point(data = filter(volcano_data, significance_level == "not_significant"),
             size = 3, alpha = 0.3, shape = 21,
             fill = "grey80", color = "grey60") +

  # p < 0.05 significant points
  geom_point(data = filter(volcano_data, significance_level == "p_significant"),
             aes(fill = point_color),
             size = 3, alpha = 0.5, shape = 21,
             color = "grey70", stroke = 0.8) +

  # FDR < 0.05 significant points
  geom_point(data = filter(volcano_data, significance_level == "FDR_significant"),
             aes(fill = point_color,
                 color = after_scale(scales::alpha(fill, 0.4))),
             size = 3, alpha = 0.8, shape = 21, stroke = 1.2) +

  # Label top genes
  geom_text_repel(aes(label = to_label),
                  size = 3.5,
                  box.padding = 0.5,
                  point.padding = 0.2,
                  min.segment.length = 0.2,
                  segment.color = "grey50",
                  segment.alpha = 0.7,
                  max.overlaps = Inf) +

  # Color gradient
  scale_fill_gradient2(low = "#86acb9", mid = "white", high = "#9b383a",
                       midpoint = 0, name = "Effect Size",
                       guide = guide_colorbar(barheight = unit(4, "cm"))) +
  scale_color_identity() +

  # Reference lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", alpha = 0.7) +

  # Labels and title
  labs(title = "Volcano Plot of Gene Expression Differences",
       subtitle = "Outliers vs Normal Samples",
       x = "Standard Effect Size",
       y = "-log10(p-value)") +

  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12)
  )

# Display the plot
print(volcano_plot)

# Print summary of significance levels
cat("Summary of significance levels:\n")
table(volcano_data$significance_level)

