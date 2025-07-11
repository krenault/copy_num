## Katia Renault
## PGLMM to identify genes associated with longevity   

########################################################################################################
# 1. Zero filtered poisson PGLMM pipeline
########################################################################################################
# Excludes species with 0 copies and uses Poisson family for count data
# Model: copy_number ~ longevity + (1|species__) [only for species with copy_number > 0]

library(ape)
library(dplyr)
library(phyr)
library(tidyr)
library(parallel)
library(foreach)
library(doParallel)

OUTPUT_DIR <- "/Nori_1/krenault/copy_num/results/phyr_pglmm_zero_filtered_poisson"
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")
log_progress <- function(message, timestamp = TRUE) {
  if(timestamp) {
    cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n", sep = "")
  } else {
    cat(message, "\n")
  }
  flush.console()
}
log_gene_result <- function(gene_name, trait, result, current_count, total_count) {
  status <- "Unknown"
  p_val <- NA
  estimate <- NA
  n_species <- 0
  n_total_species <- 0
  n_nonzero <- 0
  n_zeros_excluded <- 0
  mean_copy <- NA
  if(!is.null(result) && is.data.frame(result) && nrow(result) > 0) {
    status <- as.character(result$status[1])
    if("p_value_longevity" %in% colnames(result)) p_val <- result$p_value_longevity[1]
    if("estimate_longevity" %in% colnames(result)) estimate <- result$estimate_longevity[1]
    if("n_species" %in% colnames(result)) n_species <- result$n_species[1]
    if("n_total_species" %in% colnames(result)) n_total_species <- result$n_total_species[1]
    if("n_nonzero_species" %in% colnames(result)) n_nonzero <- result$n_nonzero_species[1]
    if("n_zeros_excluded" %in% colnames(result)) n_zeros_excluded <- result$n_zeros_excluded[1]
    if("mean_copy_number" %in% colnames(result)) mean_copy <- result$mean_copy_number[1]
  }
  if(!is.na(status) && status == "success" && !is.na(p_val)) {
    significance <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else if(p_val < 0.05) "*" else ""
    cat(sprintf("[%s] %d/%d ✓ %s vs %s | β=%.3f, p=%.2e%s | Mean=%.2f | N=%d/%d (-%d zeros)\n",
                format(Sys.time(), "%H:%M:%S"), current_count, total_count,
                gene_name, trait, 
                ifelse(is.na(estimate), 0, estimate), 
                ifelse(is.na(p_val), 1, p_val), significance,
                ifelse(is.na(mean_copy), 0, mean_copy),
                n_species, n_total_species, n_zeros_excluded))
  } else {
    cat(sprintf("[%s] %d/%d ✗ %s vs %s | %s | N=%d/%d (-%d zeros)\n",
                format(Sys.time(), "%H:%M:%S"), current_count, total_count,
                gene_name, trait, substr(status, 1, 50), n_species, n_total_species, n_zeros_excluded))
  }
  flush.console()
}

############################################
# A. Data preparation 
############################################

prepare_zero_filtered_pglmm_data <- function(gene_row, metadata, tree, trait_column = "MLres") {
  species_in_all <- intersect(
    intersect(names(gene_row), metadata$Scientific_name),
    tree$tip.label
  )
  total_shared_species <- length(species_in_all)
  gene_df <- data.frame(
    species = species_in_all,
    copy_number = as.integer(gene_row[species_in_all]),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(copy_number), copy_number >= 0)  
  n_zeros <- sum(gene_df$copy_number == 0)
  # FILTER OUT ZEROS 
  gene_df_nonzero <- gene_df %>% filter(copy_number > 0)
  nonzero_species_count <- nrow(gene_df_nonzero)
  if(nonzero_species_count < 10) {
    return(list(
      data = data.frame(),
      tree = NULL,
      total_species = total_shared_species,
      nonzero_species = nonzero_species_count,
      zeros_excluded = n_zeros,
      status = "Insufficient nonzero species"
    ))
  }
  if(var(gene_df_nonzero$copy_number) == 0) {
    return(list(
      data = data.frame(),
      tree = NULL,
      total_species = total_shared_species,
      nonzero_species = nonzero_species_count,
      zeros_excluded = n_zeros,
      status = "No variation in nonzero copy numbers"
    ))
  }
  merged <- merge(gene_df_nonzero, metadata, by.x = "species", by.y = "Scientific_name")
  merged <- merged[!is.na(merged[[trait_column]]), ]
  if(nrow(merged) < 10) {
    return(list(
      data = data.frame(),
      tree = NULL,
      total_species = total_shared_species,
      nonzero_species = nrow(merged),
      zeros_excluded = n_zeros,
      status = "Insufficient species after metadata merge"
    ))
  }
  # prune tree to final species set (only nonzero species)
  final_tree <- keep.tip(tree, merged$species)
  return(list(
    data = merged,
    tree = final_tree,
    total_species = total_shared_species,
    nonzero_species = nrow(merged),
    zeros_excluded = n_zeros,
    status = "success"
  ))
}

############################################
# B. PGLMM function
############################################

run_zero_filtered_poisson_pglmm <- function(gene_name, gene_row, metadata, tree, trait_column = "MLres") {
  pglmm_input <- prepare_zero_filtered_pglmm_data(gene_row, metadata, tree, trait_column)
  if (pglmm_input$status != "success" || nrow(pglmm_input$data) < 10) {
    return(data.frame(
      gene = gene_name,
      trait = trait_column,
      estimate_longevity = NA,
      std_error_longevity = NA,
      p_value_longevity = NA,
      n_species = ifelse(is.null(pglmm_input$data) || nrow(pglmm_input$data) == 0, 0, nrow(pglmm_input$data)),
      n_nonzero_species = pglmm_input$nonzero_species,
      n_zeros_excluded = pglmm_input$zeros_excluded,
      n_total_species = pglmm_input$total_species,
      mean_copy_number = NA,
      copy_number_range = NA,
      status = pglmm_input$status
    ))
  }
  tryCatch({
    formula_str <- paste("copy_number ~", trait_column, "+ (1|species__)")
    model <- pglmm(
      as.formula(formula_str),
      data = pglmm_input$data,
      cov_ranef = list(species = pglmm_input$tree),
      family = "poisson",  
      REML = FALSE,        
      verbose = FALSE
    )
    if ("B" %in% names(model)) {
      fixed_effects <- model$B
      if (length(fixed_effects) >= 2) {
        estimate_longevity <- fixed_effects[2]  
        std_error_longevity <- if ("B.se" %in% names(model)) model$B.se[2] else NA
        p_value_longevity <- if ("B.pvalue" %in% names(model)) model$B.pvalue[2] else NA
        mean_copy <- mean(pglmm_input$data$copy_number, na.rm = TRUE)
        copy_range <- paste(min(pglmm_input$data$copy_number), "-", max(pglmm_input$data$copy_number))
        
        return(data.frame(
          gene = gene_name,
          trait = trait_column,
          estimate_longevity = estimate_longevity,
          std_error_longevity = std_error_longevity,
          p_value_longevity = p_value_longevity,
          n_species = nrow(pglmm_input$data),
          n_nonzero_species = pglmm_input$nonzero_species,
          n_zeros_excluded = pglmm_input$zeros_excluded,
          n_total_species = pglmm_input$total_species,
          mean_copy_number = mean_copy,
          copy_number_range = copy_range,
          status = "success"
        ))
      }
    }
    return(data.frame(
      gene = gene_name,
      trait = trait_column,
      estimate_longevity = NA,
      std_error_longevity = NA,
      p_value_longevity = NA,
      n_species = nrow(pglmm_input$data),
      n_nonzero_species = pglmm_input$nonzero_species,
      n_zeros_excluded = pglmm_input$zeros_excluded,
      n_total_species = pglmm_input$total_species,
      mean_copy_number = mean(pglmm_input$data$copy_number, na.rm = TRUE),
      copy_number_range = paste(min(pglmm_input$data$copy_number), "-", max(pglmm_input$data$copy_number)),
      status = "Could not extract coefficients"
    ))
    
  }, error = function(e) {
    return(data.frame(
      gene = gene_name,
      trait = trait_column,
      estimate_longevity = NA,
      std_error_longevity = NA,
      p_value_longevity = NA,
      n_species = nrow(pglmm_input$data),
      n_nonzero_species = pglmm_input$nonzero_species,
      n_zeros_excluded = pglmm_input$zeros_excluded,
      n_total_species = pglmm_input$total_species,
      mean_copy_number = NA,
      copy_number_range = NA,
      status = as.character(e$message)
    ))
  })
}

############################################
# C. Parallel runner
############################################
run_parallel_zero_filtered_poisson <- function(gene_copy_data, metadata, tree,
                                               trait_columns = c("MLres"),
                                               n_cores = 10) {
  
  log_progress("=== RUNNING ZERO-FILTERED POISSON PGLMM ANALYSIS ===")
  log_progress("MODEL: copy_number ~ longevity + (1|species__)")
  log_progress("FAMILY: Poisson (appropriate for count data)")
  log_progress("FILTERING: EXCLUDES species with copy_number = 0")
  log_progress("RESPONSE: Copy numbers 1, 2, 3, ... (no zeros)")
  log_progress("INTERPRETATION: β > 0 means longer-lived species have more copies (among species with gene)")
  log_progress(paste("Analyzing", nrow(gene_copy_data), "genes"))
  log_progress(paste("Traits:", paste(trait_columns, collapse = ", ")))
  log_progress(paste("Using", n_cores, "cores"))

  gene_names <- rownames(gene_copy_data)
  gene_trait_combinations <- expand.grid(
    gene = gene_names,
    trait = trait_columns,
    stringsAsFactors = FALSE
  )
  total_combinations <- nrow(gene_trait_combinations)
  log_progress(paste("Total combinations:", total_combinations))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterExport(cl, c("gene_copy_data", "metadata", "tree", "trait_columns",
                      "prepare_zero_filtered_pglmm_data", "run_zero_filtered_poisson_pglmm"),
                envir = environment())
  clusterEvalQ(cl, {
    library(ape)
    library(dplyr)
    library(phyr)
  })
  batch_size <- 50
  n_batches <- ceiling(total_combinations / batch_size)
  log_progress(paste("Processing in", n_batches, "batches"))
  log_progress("Real-time results (β = effect of longevity on copy number, excluding zeros):")
  log_progress(strrep("-", 100))
  all_results <- list()
  current_count <- 0
  analysis_start_time <- Sys.time()
  for(batch_num in 1:n_batches) {
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx <- min(batch_num * batch_size, total_combinations)
    batch_combinations <- gene_trait_combinations[start_idx:end_idx, ]
    batch_results <- foreach(i = 1:nrow(batch_combinations),
                             .combine = 'rbind',
                             .packages = c('ape', 'dplyr', 'phyr'),
                             .errorhandling = 'pass') %dopar% {
                               gene_name <- batch_combinations$gene[i]
                               trait_col <- batch_combinations$trait[i]
                               gene_row <- gene_copy_data[gene_name, , drop = FALSE]
                               run_zero_filtered_poisson_pglmm(gene_name, gene_row, metadata, tree, trait_col)
                             }
    if(is.data.frame(batch_results)) {
      for(i in 1:nrow(batch_results)) {
        current_count <- current_count + 1
        result <- batch_results[i, ]
        log_gene_result(result$gene, result$trait, result, current_count, total_combinations)
      }
      all_results[[batch_num]] <- batch_results
      if(batch_num %% 10 == 0 || batch_num == n_batches) {
        elapsed <- round(as.numeric(Sys.time() - analysis_start_time, units = "mins"), 1)
        rate <- round(current_count / elapsed, 1)
        remaining <- round((total_combinations - current_count) / rate, 1)
        log_progress(sprintf("PROGRESS: %d/%d (%.1f%%) | %.1f/min | ETA: %.1f min",
                             current_count, total_combinations,
                             100 * current_count / total_combinations, rate, remaining))
      }
      if(batch_num %% 20 == 0) {
        temp_results <- do.call(rbind, all_results)
        if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
        if(!dir.exists(CHECKPOINT_DIR)) dir.create(CHECKPOINT_DIR, recursive = TRUE)
        checkpoint_file <- file.path(CHECKPOINT_DIR, paste0("progress_batch_", batch_num, "_of_", n_batches, ".csv"))
        write.csv(temp_results, checkpoint_file, row.names = FALSE)
        latest_file <- file.path(OUTPUT_DIR, "latest_results.csv")
        write.csv(temp_results, latest_file, row.names = FALSE)
        log_progress(paste("✓ Progress saved:", nrow(temp_results), "results so far"))
        successful_so_far <- sum(temp_results$status == "success", na.rm = TRUE)
        if(successful_so_far > 0) {
          significant_so_far <- sum(temp_results$status == "success" & temp_results$p_value_longevity < 0.05, na.rm = TRUE)
          log_progress(sprintf("  Current stats: %d successful, %d significant (%.1f%%)", 
                               successful_so_far, significant_so_far,
                               100 * significant_so_far / successful_so_far))
          avg_species <- round(mean(temp_results$n_species, na.rm = TRUE), 1)
          avg_total <- round(mean(temp_results$n_total_species, na.rm = TRUE), 1)
          avg_zeros_excluded <- round(mean(temp_results$n_zeros_excluded, na.rm = TRUE), 1)
          log_progress(sprintf("    Average per gene: %.1f analyzed / %.1f total (%.1f zeros excluded)", 
                               avg_species, avg_total, avg_zeros_excluded))
        }
      }
    }
  }
  stopCluster(cl)
  if(length(all_results) > 0) {
    final_results <- do.call(rbind, all_results)
    successful_results <- final_results[final_results$status == "success" & 
                                          !is.na(final_results$p_value_longevity), ]
    if(nrow(successful_results) > 0) {
      final_results$adjusted_p_longevity <- NA
      final_results[final_results$status == "success" & !is.na(final_results$p_value_longevity), "adjusted_p_longevity"] <-
        p.adjust(successful_results$p_value_longevity, method = "BH")
    }
    return(final_results)
  } else {
    return(data.frame())
  }
}

############################################
# D. Results summary
############################################
summarize_zero_filtered_results <- function(results) {
  log_progress("=== RESULTS SUMMARY: ZERO-FILTERED POISSON PGLMM ===")
  if(nrow(results) == 0) {
    log_progress("No results to summarize")
    return(NULL)
  }
  successful <- results[results$status == "success" & 
                          !is.na(results$p_value_longevity), ]
  if(nrow(successful) > 0) {
    log_progress(paste("Successful analyses:", nrow(successful)))
    log_progress(paste("Significant (p < 0.05):", sum(successful$p_value_longevity < 0.05, na.rm = TRUE)))
    log_progress(paste("FDR significant (q < 0.05):", sum(successful$adjusted_p_longevity < 0.05, na.rm = TRUE)))
    positive_effects <- sum(successful$estimate_longevity > 0, na.rm = TRUE)
    negative_effects <- sum(successful$estimate_longevity < 0, na.rm = TRUE)
    log_progress(sprintf("Effect directions: %d positive, %d negative", positive_effects, negative_effects))
    avg_species <- round(mean(successful$n_species, na.rm = TRUE), 1)
    avg_total <- round(mean(successful$n_total_species, na.rm = TRUE), 1)
    avg_zeros_excluded <- round(mean(successful$n_zeros_excluded, na.rm = TRUE), 1)
    log_progress(sprintf("Average per gene: %.1f analyzed / %.1f total (%.1f zeros excluded)", 
                         avg_species, avg_total, avg_zeros_excluded))
    avg_mean_copy <- round(mean(successful$mean_copy_number, na.rm = TRUE), 2)
    log_progress(sprintf("Average copy number among nonzero species: %.2f", avg_mean_copy))
    top_hits <- head(successful[order(successful$p_value_longevity), ], 10)
    log_progress("\nTop 10 associations:")
    for(i in 1:nrow(top_hits)) {
      direction <- if(top_hits$estimate_longevity[i] > 0) "↑" else "↓"
      log_progress(sprintf("  %2d. %s vs %s | β=%.3f %s, p=%.2e | Mean copies=%.2f | N=%d (-%d zeros)",
                           i, top_hits$gene[i], top_hits$trait[i],
                           top_hits$estimate_longevity[i], direction,
                           top_hits$p_value_longevity[i],
                           top_hits$mean_copy_number[i],
                           top_hits$n_species[i],
                           top_hits$n_zeros_excluded[i]))
    }
    # FDR significant results
    fdr_sig <- successful[successful$adjusted_p_longevity < 0.05, ]
    if(nrow(fdr_sig) > 0) {
      log_progress(paste("\n*** FDR SIGNIFICANT RESULTS ***:", nrow(fdr_sig)))
      for(i in 1:min(5, nrow(fdr_sig))) {
        direction <- if(fdr_sig$estimate_longevity[i] > 0) "↑" else "↓"
        log_progress(sprintf("  ★ %s vs %s | β=%.3f %s, q=%.2e | Mean copies=%.2f | N=%d",
                             fdr_sig$gene[i], fdr_sig$trait[i],
                             fdr_sig$estimate_longevity[i], direction,
                             fdr_sig$adjusted_p_longevity[i],
                             fdr_sig$mean_copy_number[i],
                             fdr_sig$n_species[i]))
      }
    }
  }
  
  return(successful)
}

########################################################################################################
# 2. Main analysis 
########################################################################################################
############################################
# A. Data loading
############################################
log_progress("=== LOADING DATA ===")
metadata <- read.csv("/Nori_1/krenault/copy_num/data/raxml_final_metadata_revised.csv", sep = '\t')
log_progress("✓ Metadata loaded")
gene_copy_data <- read.csv("/Nori_1/krenault/copy_num/data/All_Species_Orthologous_CopyNumber_Annotated.tsv",
                           row.names = "t_gene", sep = '\t')
log_progress("✓ Gene copy data loaded")

gene_copy_raw <- read.csv("/Nori_1/krenault/copy_num/data/All_Species_Orthologous_CopyNumber_Annotated.tsv", sep = '\t')
if("t_symbol" %in% colnames(gene_copy_raw)) {
  gene_name_mapping <- setNames(gene_copy_raw$t_symbol, gene_copy_raw$t_gene)
  valid_symbols <- gene_name_mapping[rownames(gene_copy_data)]
  valid_symbols <- valid_symbols[!is.na(valid_symbols) & valid_symbols != ""]
  
  if(length(valid_symbols) > nrow(gene_copy_data) * 0.8) {
    rownames(gene_copy_data) <- ifelse(rownames(gene_copy_data) %in% names(valid_symbols),
                                       valid_symbols[rownames(gene_copy_data)],
                                       rownames(gene_copy_data))
    log_progress("✓ Using gene symbols as identifiers")
  }
}

#original_count <- nrow(gene_copy_data)
#gene_copy_data <- gene_copy_data[apply(gene_copy_data, 1, function(x) sum(x == 0, na.rm=TRUE)/sum(!is.na(x))) <= 0.8, ]
#gene_copy_data <- head(gene_copy_data)
#log_progress(paste("Filtered genes:", original_count, "→", nrow(gene_copy_data)))
tree <- read.tree("/Nori_1/krenault/copy_num/data/raxml_final_species_tree_revised.nwk")
log_progress("✓ Tree loaded")
metadata_processed <- metadata %>%
  dplyr::select(Scientific_name, MLres, order) %>%
  distinct() %>%
  filter(!is.na(MLres)) %>%
  mutate(
    Scientific_name = stringr::str_replace_all(Scientific_name, " ", "_"),
    MLres = MLres  # Log transform
  )
log_progress("=== DATA LOADED ===")
log_progress(paste("Dataset:", nrow(gene_copy_data), "genes x", ncol(gene_copy_data), "species"))
log_progress("Testing with first gene...")
test_gene <- gene_copy_data[1, , drop = FALSE]
test_result <- prepare_zero_filtered_pglmm_data(test_gene, metadata_processed, tree, "MLres")
log_progress(paste("Test result:", test_result$status))
if(test_result$status == "success") {
  log_progress(paste("  -", test_result$nonzero_species, "nonzero species analyzed"))
  log_progress(paste("  -", test_result$zeros_excluded, "zeros excluded"))
  log_progress(paste("  - Copy range:", min(test_result$data$copy_number), "-", max(test_result$data$copy_number)))
}

############################################
# B. Run analysis
############################################
log_progress(strrep("=", 80))
log_progress("STARTING ZERO-FILTERED POISSON PGLMM ANALYSIS")
log_progress("Model: copy_number ~ longevity + (1|species__)")
log_progress("Family: Poisson (appropriate for count data)")
log_progress("Filtering: EXCLUDES species with copy_number = 0")
log_progress("Focus: Copy number variation among species that retain the gene")
log_progress("Interpretation:")
log_progress("  β_longevity > 0: Longer-lived species have MORE copies (among species with gene)")
log_progress("  β_longevity < 0: Longer-lived species have FEWER copies (among species with gene)")
log_progress("Conservative approach: Avoids uncertain gene loss calls")
log_progress(strrep("=", 80))

results <- run_parallel_zero_filtered_poisson(
  gene_copy_data = gene_copy_data,
  metadata = metadata_processed,
  tree = tree,
  trait_columns = c("MLres"),
  n_cores = 20
)
if(nrow(results) > 0) {
  summary_results <- summarize_zero_filtered_results(results)
  for(dir in c(OUTPUT_DIR)) {
    if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  results_file <- file.path(OUTPUT_DIR, paste0("zero_filtered_poisson_pglmm_", timestamp, ".csv"))
  write.csv(results, results_file, row.names = FALSE)
  if(!is.null(summary_results)) {
    fdr_sig <- summary_results[summary_results$adjusted_p_longevity < 0.05, ]
    if(nrow(fdr_sig) > 0) {
      fdr_file <- file.path(OUTPUT_DIR, "fdr_significant_zero_filtered.csv")
      write.csv(fdr_sig, fdr_file, row.names = FALSE)
      log_progress(paste("✓ FDR significant results saved:", fdr_file))
    }
  }
  log_progress(paste("✓ All results saved:", results_file))
}

log_progress(strrep("=", 80))
log_progress("ZERO-FILTERED POISSON ANALYSIS COMPLETE")
log_progress(strrep("=", 80))

############################################
# C. Interpretation guide
############################################
log_progress("\n=== INTERPRETATION GUIDE ===")
log_progress("Model: copy_number ~ longevity + (1|species__) [Poisson family, zeros excluded]")
log_progress("")
log_progress("ANALYSIS FOCUS:")
log_progress("  - Only analyzes species with copy_number > 0 (gene present)")
log_progress("  - Tests copy number variation among gene-retaining species")
log_progress("  - Avoids uncertain gene loss calls (copy_number = 0)")
log_progress("  - Conservative approach focusing on high-confidence data")
log_progress("")
log_progress("LONGEVITY EFFECTS:")
log_progress("  - estimate_longevity > 0: Longer-lived species have MORE copies (among species with gene)")
log_progress("  - estimate_longevity < 0: Longer-lived species have FEWER copies (among species with gene)")
log_progress("  - p_value_longevity: Significance of longevity effect")
log_progress("  - adjusted_p_longevity: FDR-corrected significance")
log_progress("")
log_progress("POISSON MODEL INTERPRETATION:")
log_progress("  - Uses log-link: exp(β) = multiplicative effect on copy number")
log_progress("  - β = 0.1 means ~10% more copies per unit longevity increase")
log_progress("  - β = -0.1 means ~10% fewer copies per unit longevity increase")
log_progress("  - Naturally handles count data structure")
log_progress("")
log_progress("SAMPLE SIZE INFORMATION:")
log_progress("  - n_species: Number of species analyzed (nonzero copies + complete data)")
log_progress("  - n_total_species: Total species available in all datasets")
log_progress("  - n_zeros_excluded: Number of species excluded due to zero copies")
log_progress("  - Higher n_zeros_excluded = more gene losses (potentially informative)")
log_progress("")
log_progress("BIOLOGICAL INTERPRETATION:")
log_progress("  - Tests gene dosage effects among species that retain the gene")
log_progress("  - Positive effects: Gene duplication benefits longevity")
log_progress("  - Negative effects: Single copy is optimal for longevity")
log_progress("  - Focuses on evolutionary fine-tuning rather than presence/absence")
log_progress("")
log_progress("COMPARISON TO OTHER APPROACHES:")
log_progress("  - More conservative than including zeros (avoids loss uncertainty)")
log_progress("  - More powerful than threshold methods (uses full count range)")
log_progress("  - More appropriate than Gaussian (respects count data nature)")
log_progress("  - Complements gene presence/absence analyses")

############################################
# D. Additional functions
############################################
examine_high_loss_genes <- function(results, loss_threshold = 10) {
  log_progress(paste("\n=== GENES WITH HIGH ZERO EXCLUSIONS (>", loss_threshold, ") ==="))
  high_loss <- results[results$n_zeros_excluded > loss_threshold & !is.na(results$n_zeros_excluded), ]
  if(nrow(high_loss) > 0) {
    high_loss <- high_loss[order(high_loss$n_zeros_excluded, decreasing = TRUE), ]
    log_progress(paste("Found", nrow(high_loss), "genes with >", loss_threshold, "zero exclusions"))
    for(i in 1:min(10, nrow(high_loss))) {
      gene_info <- high_loss[i, ]
      status_note <- if(gene_info$status == "success") {
        if(!is.na(gene_info$p_value_longevity) && gene_info$p_value_longevity < 0.05) "SIGNIFICANT" else "non-sig"
      } else "failed"
      log_progress(sprintf("  %2d. %s: %d zeros excluded, %d analyzed, %s", 
                           i, gene_info$gene, gene_info$n_zeros_excluded, 
                           gene_info$n_species, status_note))
    }
    successful_high_loss <- high_loss[high_loss$status == "success", ]
    if(nrow(successful_high_loss) > 0) {
      sig_rate <- mean(successful_high_loss$p_value_longevity < 0.05, na.rm = TRUE)
      log_progress(sprintf("\nHigh-loss gene stats: %.1f%% significant among successful analyses", 
                           100 * sig_rate))
    }
  } else {
    log_progress("No genes found with high zero exclusions")
  }
  return(high_loss)
}

analyze_copy_number_effects <- function(results) {
  log_progress("\n=== COPY NUMBER RANGE ANALYSIS ===")
  successful <- results[results$status == "success" & !is.na(results$mean_copy_number), ]
  if(nrow(successful) > 0) {
    # Categorize genes by mean copy number
    successful$copy_category <- cut(successful$mean_copy_number,
                                    breaks = c(0, 1.5, 2.5, 5, Inf),
                                    labels = c("Low (1-1.5)", "Medium (1.5-2.5)", "High (2.5-5)", "Very High (5+)"),
                                    include.lowest = TRUE)
    category_summary <- successful %>%
      group_by(copy_category) %>%
      summarise(
        n_genes = n(),
        mean_estimate = round(mean(estimate_longevity, na.rm = TRUE), 4),
        median_estimate = round(median(estimate_longevity, na.rm = TRUE), 4),
        significant = sum(p_value_longevity < 0.05, na.rm = TRUE),
        sig_rate = round(100 * significant / n_genes, 1),
        avg_sample_size = round(mean(n_species), 1),
        .groups = 'drop'
      )
    log_progress("Effect sizes by copy number range:")
    print(category_summary)
    if(nrow(successful) > 20) {
      tryCatch({
        aov_result <- aov(estimate_longevity ~ copy_category, data = successful)
        p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
        if(!is.na(p_value)) {
          log_progress(sprintf("ANOVA p-value for category differences: %.4f %s", 
                               p_value, if(p_value < 0.05) "(significant)" else "(not significant)"))
        }
      }, error = function(e) {
        log_progress("Could not perform ANOVA on categories")
      })
    }
  }
  return(successful)
}

identify_gene_patterns <- function(results, pattern_analysis = TRUE) {
  log_progress("\n=== GENE PATTERN ANALYSIS ===")
  successful <- results[results$status == "success", ]
  significant <- successful[successful$p_value_longevity < 0.05, ]
  if(nrow(significant) > 0) {
    log_progress(paste("Analyzing", nrow(significant), "significant genes for patterns"))
    positive_genes <- significant[significant$estimate_longevity > 0, ]
    negative_genes <- significant[significant$estimate_longevity < 0, ]
    log_progress(sprintf("Direction split: %d positive effects, %d negative effects", 
                         nrow(positive_genes), nrow(negative_genes)))
    if(pattern_analysis && nrow(significant) > 5) {
      gene_names <- significant$gene
      prefixes <- substr(gene_names, 1, 3)
      prefix_counts <- table(prefixes)
      common_prefixes <- prefix_counts[prefix_counts >= 2]
      if(length(common_prefixes) > 0) {
        log_progress("\nPotential gene family enrichment (≥2 genes with same 3-letter prefix):")
        for(i in 1:min(5, length(common_prefixes))) {
          prefix <- names(common_prefixes)[i]
          count <- common_prefixes[i]
          genes_with_prefix <- gene_names[substr(gene_names, 1, 3) == prefix]
          log_progress(sprintf("  %s*: %d genes (%s)", prefix, count, 
                               paste(head(genes_with_prefix, 3), collapse = ", ")))
        }
      }
      numeric_genes <- gene_names[grepl("\\d", gene_names)]
      if(length(numeric_genes) > 2) {
        log_progress(sprintf("\nGenes with numbers (potential clusters): %d", length(numeric_genes)))
        log_progress(paste("  Examples:", paste(head(numeric_genes, 5), collapse = ", ")))
      }
    }
    effect_quartiles <- quantile(significant$estimate_longevity, c(0.25, 0.5, 0.75))
    log_progress(sprintf("\nEffect size distribution (significant genes):"))
    log_progress(sprintf("  Q1: %.4f, Median: %.4f, Q3: %.4f", 
                         effect_quartiles[1], effect_quartiles[2], effect_quartiles[3]))
    strongest_positive <- head(positive_genes[order(positive_genes$estimate_longevity, decreasing = TRUE), ], 3)
    strongest_negative <- head(negative_genes[order(negative_genes$estimate_longevity), ], 3)
    
    if(nrow(strongest_positive) > 0) {
      log_progress("\nStrongest positive effects:")
      for(i in 1:nrow(strongest_positive)) {
        gene_info <- strongest_positive[i, ]
        log_progress(sprintf("  %s: β=%.4f, p=%.2e, N=%d", 
                             gene_info$gene, gene_info$estimate_longevity, 
                             gene_info$p_value_longevity, gene_info$n_species))
      }
    }
    
    if(nrow(strongest_negative) > 0) {
      log_progress("\nStrongest negative effects:")
      for(i in 1:nrow(strongest_negative)) {
        gene_info <- strongest_negative[i, ]
        log_progress(sprintf("  %s: β=%.4f, p=%.2e, N=%d", 
                             gene_info$gene, gene_info$estimate_longevity, 
                             gene_info$p_value_longevity, gene_info$n_species))
      }
    }
  }
  
  return(significant)
}

############################################
# E. Complete analysis wrapper
############################################

run_complete_zero_filtered_analysis <- function(gene_copy_data, metadata_processed, tree, n_cores = 20) {
  log_progress("=== STARTING COMPLETE ZERO-FILTERED POISSON ANALYSIS ===")
  results <- run_parallel_zero_filtered_poisson(
    gene_copy_data = gene_copy_data,
    metadata = metadata_processed,
    tree = tree,
    trait_columns = c("MLres"),
    n_cores = n_cores
  )
  if(nrow(results) > 0) {
    summary_results <- summarize_zero_filtered_results(results)
    high_loss_genes <- examine_high_loss_genes(results, loss_threshold = 10)
    copy_analysis <- analyze_copy_number_effects(results)
    pattern_analysis <- identify_gene_patterns(results)
    log_progress("\n=== FINAL RECOMMENDATIONS ===")
    successful <- results[results$status == "success", ]
    if(nrow(successful) > 0) {
      success_rate <- 100 * nrow(successful) / nrow(results)
      sig_rate <- 100 * sum(successful$p_value_longevity < 0.05, na.rm = TRUE) / nrow(successful)
      log_progress(sprintf("Analysis completed successfully: %.1f%% success rate", success_rate))
      log_progress(sprintf("Significance rate: %.1f%% (reasonable for multiple testing)", sig_rate))
      avg_zeros_excluded <- mean(successful$n_zeros_excluded, na.rm = TRUE)
      log_progress(sprintf("Average zeros excluded per gene: %.1f (gene loss frequency)", avg_zeros_excluded))
      if(exists("summary_results") && !is.null(summary_results)) {
        fdr_sig_count <- sum(summary_results$adjusted_p_longevity < 0.05, na.rm = TRUE)
        if(fdr_sig_count > 0) {
          log_progress(sprintf("FDR significant results: %d genes (high confidence)", fdr_sig_count))
          log_progress("→ Focus on these for follow-up studies")
        }
      }
      log_progress("\nNext steps:")
      log_progress("1. Examine FDR significant genes for biological relevance")
      log_progress("2. Consider pathway/GO enrichment analysis")
      log_progress("3. Validate key findings with independent data")
      log_progress("4. Compare with gene presence/absence analysis")
    }
    return(list(
      results = results,
      summary = summary_results,
      high_loss_genes = high_loss_genes,
      copy_analysis = copy_analysis,
      patterns = pattern_analysis
    ))
  } else {
    log_progress("Analysis failed - no results obtained")
    return(NULL)
  }
}

############################################
# COMPLETE ANALYSIS
############################################
log_progress("\n=== ZERO-FILTERED POISSON PGLMM PIPELINE READY ===")
log_progress("Complete pipeline loaded with all analysis functions")
log_progress("")
log_progress("To run the full analysis:")
log_progress("  complete_analysis <- run_complete_zero_filtered_analysis(gene_copy_data, metadata_processed, tree)")
log_progress("")
log_progress("Or run components separately:")
log_progress("  results <- run_parallel_zero_filtered_poisson(gene_copy_data, metadata_processed, tree)")
log_progress("  summary <- summarize_zero_filtered_results(results)")
log_progress("  patterns <- identify_gene_patterns(results)")
metadata <- read.csv("/Nori_1/krenault/copy_num/data/raxml_final_metadata_revised.csv", sep = '\t')
log_progress("✓ Metadata loaded")
gene_copy_data <- read.csv("/Nori_1/krenault/copy_num/data/All_Species_Orthologous_CopyNumber_Annotated.tsv",
                           row.names = "t_gene", sep = '\t')
log_progress("✓ Gene copy data loaded")
if("t_symbol" %in% colnames(gene_copy_data)) {
  rownames(gene_copy_data) <- gene_copy_data$t_symbol
  gene_copy_data <- gene_copy_data %>% dplyr::select(-t_symbol)
}
tree <- read.tree("/Nori_1/krenault/copy_num/data/raxml_final_species_tree_revised.nwk")
log_progress("✓ Tree loaded")
metadata <- metadata %>%
  mutate(
    adult_body_mass_g = log10(adult_body_mass_g),
    maximum_longevity_y = log10(maximum_longevity_y)
  )
log_progress("=== DATA LOADED ===")
log_progress(paste("Dataset:", nrow(gene_copy_data), "genes x", ncol(gene_copy_data), "species"))
log_progress("Testing with first gene...")
test_gene <- gene_copy_data[1, , drop = FALSE]
test_result <- prepare_zero_filtered_pglmm_data(test_gene, metadata_processed, tree, "MLres")  # ✅ CORRECT FUNCTION
log_progress(paste("Test result:", test_result$status))
if(test_result$status == "success") {
  log_progress(paste("  -", test_result$nonzero_species, "nonzero species analyzed"))  # ✅ CORRECT FIELD
  log_progress(paste("  -", test_result$zeros_excluded, "zeros excluded"))  # ✅ CORRECT FIELD
  log_progress(paste("  - Copy range:", min(test_result$data$copy_number), "-", max(test_result$data$copy_number)))
}

results <- run_parallel_zero_filtered_poisson(
  gene_copy_data = gene_copy_data,
  metadata = metadata_processed,  # ✅ Use processed metadata
  tree = tree,
  trait_columns = c("MLres"),
  n_cores = 10  
)
summary <- summarize_zero_filtered_results(results)
significant_genes <- results[results$adjusted_p_longevity < 0.05 & !is.na(results$adjusted_p_longevity), ]
