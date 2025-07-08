############################################
# # Load required libraries
# library(glmnet)
# library(caret)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(patchwork) # For combining plots
# library(viridis)
# library(RColorBrewer)
# library(doParallel) # For parallel processing (optional)
# 
# ##########################################
# # Elastic net model based on copy number #
# ##########################################

###################
### New version ###
###################

# Load required libraries
library(glmnet)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork) # For combining plots
library(viridis)
library(RColorBrewer)
library(doParallel) # For parallel processing (optional)

########################################### 
## Elastic net model based on copy number ##
###########################################
#seed <- 2024
# 1258
# 3411
# 9207
seed = 13
set.seed(seed)
library(glmnet)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(tibble)     
library(viridis)    
library(RColorBrewer) 

gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber_Annotated.tsv", row.names = "t_gene", sep = '\t')
rownames(gene_copy_data) <- gene_copy_data$t_symbol
gene_copy_data <- gene_copy_data %>% dplyr::select(-t_symbol)
gene_copy_data <- gene_copy_data[apply(gene_copy_data, 1, function(x) sum(x == 0, na.rm=TRUE)/sum(!is.na(x))) <= 0.5, ]

metadata <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")
gene_copy_t <- as.data.frame(t(gene_copy_data))

############ HERE IM TAKING THE FIRST SPECIES, I CAN TAKE THE AVERAGE IN THE FUTURE
gene_copy_t$species <- sub("_\\d+$", "", rownames(gene_copy_t))
gene_copy_by_species <- gene_copy_t %>%
  group_by(species) %>%
  filter(row_number() == 1) %>%
  ungroup()

species_data <- gene_copy_by_species %>%
  inner_join(metadata[, c("Scientific_name", "maximum_longevity_y", "order")], # Using lowercase order
             by = c("species" = "Scientific_name")) %>%
  filter(!is.na(maximum_longevity_y)) %>%
  mutate(log_longevity = log10(maximum_longevity_y))

print(paste("Number of species after merging:", nrow(species_data)))

# Create bins for stratified sampling
species_data$longevity_bin <- cut(species_data$log_longevity,
                                  breaks = quantile(species_data$log_longevity,
                                                    probs = seq(0, 1, 0.2)),
                                  include.lowest = TRUE)

# feature selection
valid_genes <- setdiff(colnames(species_data),
                       c("species", "maximum_longevity_y", "log_longevity", "longevity_bin", "order"))
print(paste("Total number of genes:", length(valid_genes)))

# remove near-zero variance predictors
nzv_results <- nearZeroVar(species_data[, valid_genes],
                           saveMetrics = TRUE,
                           freqCut = 99/1,        # More permissive ratio
                           uniqueCut = 3)         # More permissive uniqueness
nzv_genes <- rownames(nzv_results)[nzv_results$nzv]
valid_genes_filtered <- setdiff(valid_genes, nzv_genes)
print(paste("Genes after removing extreme near-zero variance:", length(valid_genes_filtered)))

# calculate correlation with longevity (only for filtering extreme outliers)
cors <- numeric(length(valid_genes_filtered))
names(cors) <- valid_genes_filtered

for (i in seq_along(valid_genes_filtered)) {
  g <- valid_genes_filtered[i]
  tryCatch({
    cors[i] <- cor(species_data[[g]], species_data$log_longevity, use = "pairwise.complete.obs")
    if (i %% 5000 == 0) {
      print(paste("Processed", i, "of", length(valid_genes_filtered), "genes"))
    }
  }, error = function(e) {
    cors[i] <- NA
    print(paste("Error with gene", g, ":", e$message))
  })
}

# *** CHANGED: Select top 500 correlated genes ***
cors <- cors[!is.na(cors)]
abs_cors <- abs(cors)
ranked_genes <- names(sort(abs_cors, decreasing = TRUE))
top_genes <- ranked_genes[1:min(25, length(ranked_genes))]  # CHANGED FROM 36 TO 500
print(paste("Selected top", length(top_genes), "genes by correlation"))

# create matrix with selected genes
feature_matrix <- species_data %>%
  dplyr::select(all_of(top_genes)) %>%
  as.matrix()
rownames(feature_matrix) <- species_data$species
print(paste("Final matrix dimensions:", nrow(feature_matrix), "x", ncol(feature_matrix)))

######################################## 
## stratified Train-Test Split (80/20) ##
########################################

set.seed(seed)
train_index <- createDataPartition(species_data$longevity_bin,
                                   p = 0.8,
                                   list = FALSE)

x_train <- feature_matrix[train_index, ]
y_train <- species_data$log_longevity[train_index]
x_test <- feature_matrix[-train_index, ]
y_test <- species_data$log_longevity[-train_index]

train_species <- species_data$species[train_index]
test_species <- species_data$species[-train_index]
train_orders <- species_data$order[train_index] 
test_orders <- species_data$order[-train_index] 

print(paste("Training set size:", length(y_train)))
print(paste("Test set size:", length(y_test)))

set.seed(seed)
# center and scale features 
x_train_scaled <- scale(x_train)
x_test_scaled <- scale(x_test, center = attr(x_train_scaled, "scaled:center"),
                       scale = attr(x_train_scaled, "scaled:scale"))

# *** OPTIMIZED: More comprehensive alpha search for better L1 regularization ***
alpha_values <- seq(0, 1, by = 0.05)  # More fine-grained search
results <- list()

# loop through alpha values
for (a in alpha_values) {
  print(paste("Fitting model with alpha =", a))
  
  # fit model with cross-validation
  set.seed(seed)  # Ensure reproducible CV folds
  cv_fit <- cv.glmnet(
    x = x_train_scaled,
    y = y_train,
    alpha = a,
    nfolds = 10,
    type.measure = "mse",
    # *** ADDED: Let glmnet choose lambda range automatically for better regularization ***
    nlambda = 100  # More lambda values to test
  )
  
  # store results
  results[[as.character(a)]] <- list(
    alpha = a,
    cv_fit = cv_fit,
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se,
    cvm_min = min(cv_fit$cvm)  # Minimum cross-validated error
  )
  
  print(paste("  Lambda min:", cv_fit$lambda.min))
  print(paste("  Lambda 1se:", cv_fit$lambda.1se))
  print(paste("  Minimum CV error:", min(cv_fit$cvm)))
}

# find best alpha based on minimum cross-validated error
best_alpha <- NULL
best_error <- Inf
best_lambda <- NULL

for (a in names(results)) {
  if (results[[a]]$cvm_min < best_error) {
    best_error <- results[[a]]$cvm_min
    best_alpha <- results[[a]]$alpha
    best_lambda <- results[[a]]$lambda_min
  }
}

print(paste("Best alpha:", best_alpha))
print(paste("Best lambda:", best_lambda))
print(paste("Best CV error:", best_error))

set.seed(seed)
# fit final model with best parameters
final_model <- glmnet(
  x = x_train_scaled,
  y = y_train,
  alpha = best_alpha,
  lambda = best_lambda
)

# get predictions
train_preds <- predict(final_model, newx = x_train_scaled, s = best_lambda)
test_preds <- predict(final_model, newx = x_test_scaled, s = best_lambda)

# calculate performance 
train_perf <- postResample(pred = train_preds, obs = y_train)
test_perf <- postResample(pred = test_preds, obs = y_test)
test_r2 <- test_perf["Rsquared"]
test_mse <- test_perf["RMSE"]^2
test_mae <- test_perf["MAE"]
test_rho <- cor(y_test, as.numeric(test_preds), method = "spearman")

# *** ADDED: Check for overfitting ***
train_r2 <- train_perf["Rsquared"]
overfitting_check <- train_r2 - test_r2
print(paste("Train R²:", round(train_r2, 4)))
print(paste("Test R²:", round(test_r2, 4)))
print(paste("Overfitting (difference):", round(overfitting_check, 4)))
if (overfitting_check > 0.1) {
  warning("Potential overfitting detected! Consider using lambda.1se instead of lambda.min")
}

# create performance dataframe
performance_df <- data.frame(
  Metric = c("RMSE", "R-squared", "MAE"),
  Train = c(train_perf["RMSE"], train_perf["Rsquared"], train_perf["MAE"]),
  Test = c(test_perf["RMSE"], test_perf["Rsquared"], test_perf["MAE"])
)

print("Model Performance:")
print(performance_df)

################################### 
## Most predictive genes in model ##
###################################

coef_matrix <- as.matrix(coef(final_model, s = best_lambda))
var_importance <- coef_matrix[-1, 1] 
names(var_importance) <- rownames(coef_matrix)[-1]  

# *** IMPROVED: Better handling of feature selection results ***
# Count features actually used by L1 regularization
selected_features <- names(var_importance)[var_importance != 0]
n_features_used <- length(selected_features)
n_features_input <- ncol(x_train)

print(paste("=== L1 REGULARIZATION RESULTS ==="))
print(paste("Input features:", n_features_input))
print(paste("Features actually used:", n_features_used))
print(paste("Features zeroed out by L1:", n_features_input - n_features_used))
print(paste("Feature usage rate:", round(n_features_used/n_features_input * 100, 1), "%"))

# *** CHANGED: Show all non-zero features, not just top 19 ***
non_zero_features <- var_importance[var_importance != 0]
top_features_to_show <- min(20, length(non_zero_features))
top_genes_indices <- order(abs(var_importance), decreasing = TRUE)[1:top_features_to_show]

top_important_genes <- data.frame(
  Gene = names(var_importance)[top_genes_indices],
  Effect = var_importance[top_genes_indices],
  Direction = ifelse(var_importance[top_genes_indices] > 0, "Positive impact on longevity", "Negative impact on longevity")
)

# Create comprehensive gene results
all_genes_df <- data.frame(
  Gene = names(var_importance),
  Coefficient = var_importance,
  Abs_Coefficient = abs(var_importance),
  Direction = ifelse(var_importance > 0, "Positive impact on longevity",
                     ifelse(var_importance < 0, "Negative impact on longevity", "No impact"))
)

all_genes_sorted <- all_genes_df[order(all_genes_df$Abs_Coefficient, decreasing = TRUE), ]
all_genes_sorted$Rank <- 1:nrow(all_genes_sorted)

# *** ADDED: Save results ***
write.csv(all_genes_sorted, paste0("all_predictive_genes_500_features_seed_", seed, ".csv"), row.names = FALSE)

# *** Rest of plotting code remains the same but with updated titles ***
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
    if (n > length(sophisticated_jewels_palette)) {
      return(sophisticated_jewels_palette[1:length(sophisticated_jewels_palette)])
    } else {
      return(sophisticated_jewels_palette[1:n])
    }
  }
}

test_plot_data <- data.frame(
  actual = y_test,
  predicted = as.numeric(test_preds),
  order = test_orders,
  Species = test_species
)

#### Color mapping ####
orders <- unique(test_plot_data$order)
sophisticated_jewels_palette <- c(
  "#9b383a", "#64a590", "#AA4839", "#AF6C36", "#945a87",
  "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
  "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
  "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
  "#c07a7a", "#8a7c6d", "#5c4d8b", "#4878a0", "#3D7B68"
)

if (length(orders) > length(sophisticated_jewels_palette)) {
  color_mapping <- setNames(
    sophisticated_jewels_palette[1:length(sophisticated_jewels_palette)], 
    orders[1:length(sophisticated_jewels_palette)]
  )
  warning("More orders than available colors. Some colors will be recycled.")
} else {
  color_mapping <- setNames(sophisticated_jewels_palette[1:length(orders)], orders)
}

colnames(test_plot_data)[colnames(test_plot_data) == "order"] <- "Order"

# create the test set plot with updated subtitle
p_test <- ggplot(test_plot_data, aes(x = actual, y = predicted, color = Order)) +
  geom_point(size = 4, alpha = 0.7, stroke = 0.8, shape = 21, aes(fill = Order)) +
  scale_fill_manual(values = color_mapping) +
  scale_color_manual(values = color_mapping) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  labs(title = "Test Set Performance",
       subtitle = paste0("R² = ", round(test_r2, 3),
                         ", ρ = ", round(test_rho, 3),
                         ", Features used = ", n_features_used, "/", n_features_input),
       x = "Actual maximum lifespan (log10)",
       y = "Predicted maximum lifespan (log10)") +
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

top_important_genes$Direction <- ifelse(top_important_genes$Effect > 0,
                                        "Positive impact on longevity",
                                        "Negative impact on longevity")
top_important_genes <- top_important_genes %>%
  arrange(Effect) 

direction_colors <- c("Positive impact on longevity" = "#9b383a", "Negative impact on longevity" = "#86acb9")

# create enhanced bar plot with directionality
p_genes <- ggplot(top_important_genes, aes(x = reorder(Gene, Effect), y = Effect, fill = Direction)) +
  geom_col() +
  scale_fill_manual(values = direction_colors) +
  coord_flip() +
  labs(
    title = paste0("Top ", nrow(top_important_genes), " predictive genes"),
    x = "Gene",
    y = "Effect on Longevity (model coefficient)"
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

print(p_test)
print(p_genes)

# *** ADDED: Summary of top selected features ***
cat("\n=== TOP SELECTED FEATURES ===\n")
top_selected <- sort(abs(non_zero_features), decreasing = TRUE)[1:min(10, length(non_zero_features))]
for (i in 1:length(top_selected)) {
  gene <- names(top_selected)[i]
  coef <- non_zero_features[gene]
  direction <- ifelse(coef > 0, "positive", "negative")
  cat(sprintf("%2d. %s: %.4f (%s effect)\n", i, gene, coef, direction))
}

################################################################################
# library(fgsea)
# pathways.hallmark <- gmtPathways( "/Users/katiarenault/Desktop/PhD/Human_GSEA/h.all.v2023.2.Hs.symbols.gmt")
# pathways.kegg <- gmtPathways("/Users/katiarenault/Desktop/PhD/Human_GSEA/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
# pathways.reactome <- gmtPathways(  "/Users/katiarenault/Desktop/PhD/Human_GSEA/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
# pathways.go <- gmtPathways("/Users/katiarenault/Downloads/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/c5.go.v2023.2.Hs.symbols.gmt")
# pathway_files <- c(pathways.go, pathways.hallmark, pathways.kegg, pathways.reactome)
# coef_matrix <- coef_matrix[rownames(coef_matrix) != "(Intercept)", , drop=FALSE]
# coef_matrix <- as.data.frame(coef_matrix)
# coef_matrix$gene <- rownames(coef_matrix)
# 
# fgsea_results <- coef_matrix %>%
#   as_tibble() %>%
#   # filter(coefficient > 0) %>%
#   drop_na() 
# 
# stats_vector <- setNames(fgsea_results$s1, fgsea_results$gene)
# pgls_fgsea_results <- fgsea(pathways = pathway_files, stats = stats_vector)
# pgls_fgsea_results <- data.frame(pgls_fgsea_results)

# # Save the plots
# ggsave("/Users/katiarenault/PhD/TOGA/revised_results/plots/test_prediction_plot_20_features.png", p_test, width = 10, height = 8)
# ggsave("/Users/katiarenault/PhD/TOGA/revised_results/plots/gene_importance_with_direction_20_features.png", p_genes, width = 8, height = 10)
# 
# # Save with directional information
# write.csv(performance_df, "/Users/katiarenault/PhD/TOGA/revised_results/data/model_performance_metrics_20_features.csv", row.names = FALSE)
# write.csv(top_important_genes, "/Users/katiarenault/PhD/TOGA/revised_results/data/top_predictive_genes_with_direction_20_features.csv", row.names = FALSE)

############################################
## Including body mass + gene predictions ##
############################################

set.seed(seed)

library(glmnet)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(viridis)
library(RColorBrewer)
library(purrr)

# Use the same subsetting method as your first model
gene_copy_by_species <- gene_copy_t %>%
  group_by(species) %>%
  filter(row_number() == 1) %>% # Using filter like in your first model
  ungroup()

# Create a dataset with ALL species that have longevity data (for fair comparison)
species_data <- gene_copy_by_species %>%
  inner_join(metadata[, c("Scientific_name", "maximum_longevity_y", "order", "adult_body_mass_g")],
             by = c("species" = "Scientific_name")) %>%
  filter(!is.na(maximum_longevity_y)) %>%
  mutate(log_longevity = log10(maximum_longevity_y),
         log_mass = log10(adult_body_mass_g)) # CRITICAL: Create log_mass here!

# Create bins for stratified sampling
species_data$longevity_bin <- cut(species_data$log_longevity,
                                  breaks = quantile(species_data$log_longevity,
                                                    probs = seq(0, 1, 0.2)),
                                  include.lowest = TRUE)

# Feature selection
valid_genes <- setdiff(colnames(species_data),
                       c("species", "maximum_longevity_y", "log_longevity", "longevity_bin",
                         "order", "adult_body_mass_g", "log_mass"))

print(paste("Total number of genes:", length(valid_genes)))

# Remove near-zero variance predictors
set.seed(seed) # Ensure reproducibility
nzv_results <- nearZeroVar(species_data[, valid_genes],
                           saveMetrics = TRUE,
                           freqCut = 99/1,
                           uniqueCut = 3)
nzv_genes <- rownames(nzv_results)[nzv_results$nzv]
valid_genes_filtered <- setdiff(valid_genes, nzv_genes)

print(paste("Genes after removing extreme near-zero variance:", length(valid_genes_filtered)))

# Calculate correlation with longevity
cors <- numeric(length(valid_genes_filtered))
names(cors) <- valid_genes_filtered

for (i in seq_along(valid_genes_filtered)) {
  g <- valid_genes_filtered[i]
  tryCatch({
    cors[i] <- cor(species_data[[g]], species_data$log_longevity, use = "pairwise.complete.obs")
    # Print progress every 5000 genes
    if (i %% 5000 == 0) {
      print(paste("Processed", i, "of", length(valid_genes_filtered), "genes"))
    }
  }, error = function(e) {
    cors[i] <- NA
    print(paste("Error with gene", g, ":", e$message))
  })
}

# Select top 1000 correlated genes
cors <- cors[!is.na(cors)]
abs_cors <- abs(cors)
ranked_genes <- names(sort(abs_cors, decreasing = TRUE))
top_genes <- ranked_genes[1:min(1000, length(ranked_genes))]

print(paste("Selected top", length(top_genes), "genes by correlation"))
