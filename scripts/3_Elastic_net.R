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

##########################################
# Elastic net model based on copy number #
##########################################
#set.seed(7451)
#set.seed(2024)
#seed = 3411 for 1000
#seed = 1258 for 20
#seed = 7451
#seed = 2024
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

gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/loss/Mammalia_all_species_copynumber.csv",
                           row.names = "Gene")

#gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/gene_copy_data_revised.csv")
#rownames(gene_copy_data) <- gene_copy_data$X
gene_copy_data <- gene_copy_data %>% dplyr::select(-X)
#gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/pseudogenes/Mammalia_all_species_pseudo.csv", row.names = "Gene")
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

# select top 1000 correlated genes (more features than previous 300)
cors <- cors[!is.na(cors)]
abs_cors <- abs(cors)
ranked_genes <- names(sort(abs_cors, decreasing = TRUE))
top_genes <- ranked_genes[1:min(1000, length(ranked_genes))]
print(paste("Selected top", length(top_genes), "genes by correlation"))

# create matrix with selected genes
feature_matrix <- species_data %>%
  dplyr::select(all_of(top_genes)) %>%
  as.matrix()
rownames(feature_matrix) <- species_data$species
print(paste("Final matrix dimensions:", nrow(feature_matrix), "x", ncol(feature_matrix)))

#######################################
# stratified Train-Test Split (80/20) #
#######################################

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

# define alpha values 
alpha_values <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,  0.7, 0.8, 0.9, 1)
results <- list()

# loop through alpha values
for (a in alpha_values) {
  print(paste("Fitting model with alpha =", a))
  
  # fit model with cross-validation
  cv_fit <- cv.glmnet(
    x = x_train_scaled,
    y = y_train,
    alpha = a,
    nfolds = 10,
    type.measure = "mse"
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
  lambda = best_lambda)

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

# create performance dataframe
performance_df <- data.frame(
  Metric = c("RMSE", "R-squared", "MAE"),
  Train = c(train_perf["RMSE"], train_perf["Rsquared"], train_perf["MAE"]),
  Test = c(test_perf["RMSE"], test_perf["Rsquared"], test_perf["MAE"]))

print("Model Performance:")
print(performance_df)

##################################
# Most predictive genes in model #
##################################

coef_matrix <- as.matrix(coef(final_model, s = best_lambda))
var_importance <- coef_matrix[-1, 1] 
names(var_importance) <- rownames(coef_matrix)[-1]  

top_genes_indices <- order(abs(var_importance), decreasing = TRUE)[1:20]
top_important_genes <- data.frame(
  Gene = names(var_importance)[top_genes_indices],
  Effect = var_importance[top_genes_indices],
  Direction = ifelse(var_importance[top_genes_indices] > 0, "Positive impact on longevity", "Negative impact on longevity")
)

var_importance <- coef_matrix[-1, 1]  # Skip the intercept
names(var_importance) <- rownames(coef_matrix)[-1]
non_zero_genes <- var_importance[var_importance != 0]
all_genes_df <- data.frame(
  Gene = names(var_importance),
  Coefficient = var_importance,
  Abs_Coefficient = abs(var_importance),
  Direction = ifelse(var_importance > 0, "Positive impact on longevity", 
                     ifelse(var_importance < 0, "Negative impact on longevity", "No impact"))
)
all_genes_sorted <- all_genes_df[order(all_genes_df$Abs_Coefficient, decreasing = TRUE), ]
all_genes_sorted$Rank <- 1:nrow(all_genes_sorted)
 # write.csv(all_genes_sorted, "/Users/katiarenault/PhD/TOGA/revised_results/results/all_predictive_genes_with_direction_20_features.csv",
 #           row.names = FALSE)

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
    # Handle case where n exceeds the palette length by recycling colors
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
####################################################################################
orders <- unique(test_plot_data$order)
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
####################################################################################

colnames(test_plot_data)[colnames(test_plot_data) == "order"] <- "Order" # Capitalize column name
# create the test set plot with R², rho, and MSE metrics
p_test <- ggplot(test_plot_data, aes(x = actual, y = predicted, color = Order)) +
  # Use shape = 21 for points with fill and color
  geom_point(size = 4, alpha = 0.7, stroke = 0.8, shape = 21, aes(fill = Order)) +
  scale_fill_manual(values = color_mapping) +
  scale_color_manual(values = color_mapping) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  labs(title = "Test Set Performance",
       subtitle = paste0("R² = ", round(test_r2, 3), 
                         ", ρ = ", round(test_rho, 3),
                         ", MSE = ", round(test_mse, 3)),
       x = "Actual maximum lifespan (log10)",
       y = "Predicted maximum lifespan (log10)") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),  # Style
    guides(fill = guide_legend(title = "Order"),            # Title text for fill
           color = guide_legend(title = "Order")),                 # Title text for color
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
    title = "Top predictive genes",
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

# Define seed explicitly at the beginning
# seed <- 3411
# set.seed(seed)
# 
# library(glmnet)
# library(caret)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(ggrepel)
# library(tibble)
# library(viridis)
# library(RColorBrewer)
# library(purrr)
# 
# gene_copy_data <- read.csv("/Users/katiarenault/PhD/TOGA/loss/Mammalia_all_species_copynumber.csv",
#                            row.names = "Gene")
# metadata <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")
# 
# gene_copy_t <- as.data.frame(t(gene_copy_data))
# gene_copy_t$species <- sub("_\\d+$", "", rownames(gene_copy_t))
# 
# # Use the same subsetting method as your first model
# gene_copy_by_species <- gene_copy_t %>%
#   group_by(species) %>%
#   filter(row_number() == 1) %>% # Using filter like in your first model
#   ungroup()
# 
# # Create a dataset with ALL species that have longevity data (for fair comparison)
# species_data <- gene_copy_by_species %>%
#   inner_join(metadata[, c("Scientific_name", "maximum_longevity_y", "order", "adult_body_mass_g")],
#              by = c("species" = "Scientific_name")) %>%
#   filter(!is.na(maximum_longevity_y)) %>%
#   mutate(log_longevity = log10(maximum_longevity_y),
#          log_mass = log10(adult_body_mass_g)) # CRITICAL: Create log_mass here!
# 
# # Create bins for stratified sampling
# species_data$longevity_bin <- cut(species_data$log_longevity,
#                                   breaks = quantile(species_data$log_longevity,
#                                                     probs = seq(0, 1, 0.2)),
#                                   include.lowest = TRUE)
# 
# # Feature selection
# valid_genes <- setdiff(colnames(species_data),
#                        c("species", "maximum_longevity_y", "log_longevity", "longevity_bin",
#                          "order", "adult_body_mass_g", "log_mass"))
# 
# print(paste("Total number of genes:", length(valid_genes)))
# 
# # Remove near-zero variance predictors
# set.seed(seed) # Ensure reproducibility
# nzv_results <- nearZeroVar(species_data[, valid_genes],
#                            saveMetrics = TRUE,
#                            freqCut = 99/1,
#                            uniqueCut = 3)
# nzv_genes <- rownames(nzv_results)[nzv_results$nzv]
# valid_genes_filtered <- setdiff(valid_genes, nzv_genes)
# 
# print(paste("Genes after removing extreme near-zero variance:", length(valid_genes_filtered)))
# 
# # Calculate correlation with longevity
# cors <- numeric(length(valid_genes_filtered))
# names(cors) <- valid_genes_filtered
# 
# for (i in seq_along(valid_genes_filtered)) {
#   g <- valid_genes_filtered[i]
#   tryCatch({
#     cors[i] <- cor(species_data[[g]], species_data$log_longevity, use = "pairwise.complete.obs")
#     # Print progress every 5000 genes
#     if (i %% 5000 == 0) {
#       print(paste("Processed", i, "of", length(valid_genes_filtered), "genes"))
#     }
#   }, error = function(e) {
#     cors[i] <- NA
#     print(paste("Error with gene", g, ":", e$message))
#   })
# }
# 
# # Select top 1000 correlated genes
# cors <- cors[!is.na(cors)]
# abs_cors <- abs(cors)
# ranked_genes <- names(sort(abs_cors, decreasing = TRUE))
# top_genes <- ranked_genes[1:min(1000, length(ranked_genes))]
# 
# print(paste("Selected top", length(top_genes), "genes by correlation"))
# 
# ##############################
# ## Compare the two models   ##
# ##############################
# 
# # For fair comparison, create a dataset that ONLY includes species with both longevity AND mass data
# species_with_mass <- species_data %>%
#   filter(!is.na(adult_body_mass_g))
# 
# print(paste("Species with both longevity and mass data:", nrow(species_with_mass)))
# 
# # MODEL 1: Genes only
# feature_matrix_genes_only <- species_with_mass %>%
#   dplyr::select(all_of(top_genes)) %>%
#   as.matrix()
# rownames(feature_matrix_genes_only) <- species_with_mass$species
# 
# # MODEL 2: Genes + body mass
# # Now log_mass exists, so this will work properly
# feature_matrix_with_mass <- cbind("Mass (log10)" = species_with_mass$log_mass, 
#                                   feature_matrix_genes_only)
# rownames(feature_matrix_with_mass) <- species_with_mass$species
# 
# print(paste("Genes-only matrix dimensions:", nrow(feature_matrix_genes_only), "x", ncol(feature_matrix_genes_only)))
# print(paste("Genes+mass matrix dimensions:", nrow(feature_matrix_with_mass), "x", ncol(feature_matrix_with_mass)))
# 
# # Create train/test split (using the same split for both models)
# set.seed(seed) # Ensure reproducibility
# train_index <- createDataPartition(species_with_mass$longevity_bin,
#                                    p = 0.8,
#                                    list = FALSE)
# 
# # Create train and test sets for GENES ONLY model
# x_train_genes <- feature_matrix_genes_only[train_index, ]
# y_train <- species_with_mass$log_longevity[train_index]
# x_test_genes <- feature_matrix_genes_only[-train_index, ]
# y_test <- species_with_mass$log_longevity[-train_index]
# 
# # Create train and test sets for GENES+MASS model
# x_train_mass <- feature_matrix_with_mass[train_index, ]
# x_test_mass <- feature_matrix_with_mass[-train_index, ]
# 
# # Store additional info
# train_species <- species_with_mass$species[train_index]
# test_species <- species_with_mass$species[-train_index]
# train_orders <- species_with_mass$order[train_index]
# test_orders <- species_with_mass$order[-train_index]
# train_mass <- species_with_mass$log_mass[train_index]
# test_mass <- species_with_mass$log_mass[-train_index]
# 
# print(paste("Training set size:", length(y_train)))
# print(paste("Test set size:", length(y_test)))
# 
# ##############################
# # MODEL 1: GENES-ONLY MODEL  #
# ##############################
# 
# # Center and scale features
# set.seed(seed) # Ensure reproducibility
# x_train_genes_scaled <- scale(x_train_genes)
# x_test_genes_scaled <- scale(x_test_genes,
#                              center = attr(x_train_genes_scaled, "scaled:center"),
#                              scale = attr(x_train_genes_scaled, "scaled:scale"))
# 
# # Define alpha values to try
# alpha_values <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# results_genes <- list()
# 
# # Loop through alpha values
# for (a in alpha_values) {
#   print(paste("GENES-ONLY Model: Fitting with alpha =", a))
#   
#   # Set seed before cv.glmnet which has randomization
#   set.seed(seed)
#   
#   # Fit model with cross-validation
#   cv_fit <- cv.glmnet(
#     x = x_train_genes_scaled,
#     y = y_train,
#     alpha = a,
#     nfolds = 10,
#     type.measure = "mse"
#   )
#   
#   # Store results
#   results_genes[[as.character(a)]] <- list(
#     alpha = a,
#     cv_fit = cv_fit,
#     lambda_min = cv_fit$lambda.min,
#     lambda_1se = cv_fit$lambda.1se,
#     cvm_min = min(cv_fit$cvm)
#   )
#   
#   print(paste("  Lambda min:", cv_fit$lambda.min))
#   print(paste("  Lambda 1se:", cv_fit$lambda.1se))
#   print(paste("  Minimum CV error:", min(cv_fit$cvm)))
# }
# 
# # Find best alpha
# best_alpha_genes <- NULL
# best_error_genes <- Inf
# best_lambda_genes <- NULL
# 
# for (a in names(results_genes)) {
#   if (results_genes[[a]]$cvm_min < best_error_genes) {
#     best_error_genes <- results_genes[[a]]$cvm_min
#     best_alpha_genes <- results_genes[[a]]$alpha
#     best_lambda_genes <- results_genes[[a]]$lambda_min
#   }
# }
# 
# print(paste("GENES-ONLY Model - Best alpha:", best_alpha_genes))
# print(paste("GENES-ONLY Model - Best lambda:", best_lambda_genes))
# print(paste("GENES-ONLY Model - Best CV error:", best_error_genes))
# 
# # Fit final model
# set.seed(seed) # Ensure reproducibility
# final_model_genes <- glmnet(
#   x = x_train_genes_scaled,
#   y = y_train,
#   alpha = best_alpha_genes,
#   lambda = best_lambda_genes
# )
# 
# # Get predictions
# train_preds_genes <- predict(final_model_genes, newx = x_train_genes_scaled, s = best_lambda_genes)
# test_preds_genes <- predict(final_model_genes, newx = x_test_genes_scaled, s = best_lambda_genes)
# 
# # Calculate performance metrics
# train_perf_genes <- postResample(pred = train_preds_genes, obs = y_train)
# test_perf_genes <- postResample(pred = test_preds_genes, obs = y_test)
# test_r2_genes <- test_perf_genes["Rsquared"]
# test_mse_genes <- test_perf_genes["RMSE"]^2
# test_mae_genes <- test_perf_genes["MAE"]
# test_rho_genes <- cor(y_test, as.numeric(test_preds_genes), method = "spearman")
# 
# # Create performance dataframe
# performance_df_genes <- data.frame(
#   Metric = c("RMSE", "R-squared", "MAE"),
#   Train = c(train_perf_genes["RMSE"], train_perf_genes["Rsquared"], train_perf_genes["MAE"]),
#   Test = c(test_perf_genes["RMSE"], test_perf_genes["Rsquared"], test_perf_genes["MAE"])
# )
# 
# print("GENES-ONLY Model Performance:")
# print(performance_df_genes)
# 
# ##############################
# # MODEL 2: GENES+MASS MODEL  #
# ##############################
# 
# # Center and scale features
# set.seed(seed) # Ensure reproducibility
# x_train_mass_scaled <- scale(x_train_mass)
# x_test_mass_scaled <- scale(x_test_mass,
#                             center = attr(x_train_mass_scaled, "scaled:center"),
#                             scale = attr(x_train_mass_scaled, "scaled:scale"))
# 
# # Results for GENES+MASS model
# results_mass <- list()
# 
# for (a in alpha_values) {
#   print(paste("GENES+MASS Model: Fitting with alpha =", a))
#   
#   # Set seed before cv.glmnet
#   set.seed(seed)
#   
#   # Fit model with cross-validation
#   cv_fit <- cv.glmnet(
#     x = x_train_mass_scaled,
#     y = y_train,
#     alpha = a,
#     nfolds = 10,
#     type.measure = "mse"
#   )
#   
#   # Store results
#   results_mass[[as.character(a)]] <- list(
#     alpha = a,
#     cv_fit = cv_fit,
#     lambda_min = cv_fit$lambda.min,
#     lambda_1se = cv_fit$lambda.1se,
#     cvm_min = min(cv_fit$cvm)
#   )
#   
#   print(paste("  Lambda min:", cv_fit$lambda.min))
#   print(paste("  Lambda 1se:", cv_fit$lambda.1se))
#   print(paste("  Minimum CV error:", min(cv_fit$cvm)))
# }
# 
# # Find best alpha
# best_alpha_mass <- NULL
# best_error_mass <- Inf
# best_lambda_mass <- NULL
# 
# for (a in names(results_mass)) {
#   if (results_mass[[a]]$cvm_min < best_error_mass) {
#     best_error_mass <- results_mass[[a]]$cvm_min
#     best_alpha_mass <- results_mass[[a]]$alpha
#     best_lambda_mass <- results_mass[[a]]$lambda_min
#   }
# }
# 
# print(paste("GENES+MASS Model - Best alpha:", best_alpha_mass))
# print(paste("GENES+MASS Model - Best lambda:", best_lambda_mass))
# print(paste("GENES+MASS Model - Best CV error:", best_error_mass))
# 
# # Fit final model
# set.seed(seed) # Ensure reproducibility
# final_model_mass <- glmnet(
#   x = x_train_mass_scaled,
#   y = y_train,
#   alpha = best_alpha_mass,
#   lambda = best_lambda_mass
# )
# 
# # Get predictions
# train_preds_mass <- predict(final_model_mass, newx = x_train_mass_scaled, s = best_lambda_mass)
# test_preds_mass <- predict(final_model_mass, newx = x_test_mass_scaled, s = best_lambda_mass)
# 
# # Calculate performance metrics
# train_perf_mass <- postResample(pred = train_preds_mass, obs = y_train)
# test_perf_mass <- postResample(pred = test_preds_mass, obs = y_test)
# test_r2_mass <- test_perf_mass["Rsquared"]
# test_mse_mass <- test_perf_mass["RMSE"]^2
# test_mae_mass <- test_perf_mass["MAE"]
# test_rho_mass <- cor(y_test, as.numeric(test_preds_mass), method = "spearman")
# 
# # Create performance dataframe
# performance_df_mass <- data.frame(
#   Metric = c("RMSE", "R-squared", "MAE"),
#   Train = c(train_perf_mass["RMSE"], train_perf_mass["Rsquared"], train_perf_mass["MAE"]),
#   Test = c(test_perf_mass["RMSE"], test_perf_mass["Rsquared"], test_perf_mass["MAE"])
# )
# 
# print("GENES+MASS Model Performance:")
# print(performance_df_mass)
# 
# ################################
# # Compare model performances   #
# ################################
# 
# comparison_df <- data.frame(
#   Metric = c("RMSE", "R-squared", "MAE", "Spearman's ρ"),
#   Genes_Only = c(test_perf_genes["RMSE"], test_perf_genes["Rsquared"],
#                  test_perf_genes["MAE"], test_rho_genes),
#   Genes_Plus_Mass = c(test_perf_mass["RMSE"], test_perf_mass["Rsquared"],
#                       test_perf_mass["MAE"], test_rho_mass),
#   Improvement = c(
#     test_perf_genes["RMSE"] - test_perf_mass["RMSE"],
#     test_perf_mass["Rsquared"] - test_perf_genes["Rsquared"],
#     test_perf_genes["MAE"] - test_perf_mass["MAE"],
#     test_rho_mass - test_rho_genes
#   ),
#   Percent_Change = c(
#     (test_perf_genes["RMSE"] - test_perf_mass["RMSE"]) / test_perf_genes["RMSE"] * 100,
#     (test_perf_mass["Rsquared"] - test_perf_genes["Rsquared"]) / test_perf_genes["Rsquared"] * 100,
#     (test_perf_genes["MAE"] - test_perf_mass["MAE"]) / test_perf_genes["MAE"] * 100,
#     (test_rho_mass - test_rho_genes) / test_rho_genes * 100
#   )
# )
# 
# print("Model Performance Comparison:")
# print(comparison_df)
# 
# # Feature importance = genes and mass
# coef_matrix_mass <- as.matrix(coef(final_model_mass, s = best_lambda_mass))
# var_importance_mass <- coef_matrix_mass[-1, 1]
# names(var_importance_mass) <- rownames(coef_matrix_mass)[-1]  
# 
# # Find mass importance and rank
# mass_importance <- var_importance_mass["Mass (log10)"]
# mass_rank <- which(names(sort(abs(var_importance_mass), decreasing = TRUE)) == "Mass (log10)")
# print(paste("Body mass coefficient:", mass_importance))
# print(paste("Body mass rank among all predictors:", mass_rank, "out of", length(var_importance_mass)))
# 
# # Sort by absolute importance and get top predictors
# top_predictors_indices <- order(abs(var_importance_mass), decreasing = TRUE)[1:20]
# top_important_predictors <- data.frame(
#   Predictor = names(var_importance_mass)[top_predictors_indices],
#   Effect = var_importance_mass[top_predictors_indices],
#   Direction = ifelse(var_importance_mass[top_predictors_indices] > 0,
#                      "Positive impact on longevity", "Negative impact on longevity")
# )
# 
# # Test set performance comparison
# test_plot_data <- data.frame(
#   actual = rep(y_test, 2),
#   predicted = c(as.numeric(test_preds_genes), as.numeric(test_preds_mass)),
#   order = rep(test_orders, 2),
#   Species = rep(test_species, 2),
#   Model = rep(c("Genes", "Genes + Body Mass"), each = length(y_test))
# )
# 
# # Create the comparison plot
# p_comparison <- ggplot(test_plot_data, aes(x = actual, y = predicted, color = order)) +
#   geom_point(size = 3.5, alpha = 0.7, stroke = 0.8, shape = 21, aes(fill = order)) +
#   scale_fill_manual(values = color_mapping) +
#   scale_color_manual(values = color_mapping) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
#   facet_wrap(~Model) +
#   labs(title = "Model Performance Comparison - Test Set",
#        subtitle = paste0("Genes copy number only: R² = ", round(test_r2_genes, 3),
#                          ", ρ = ", round(test_rho_genes, 3), "\n",
#                          "Genes copy number + body mass: R² = ", round(test_r2_mass, 3),
#                          ", ρ = ", round(test_rho_mass, 3)),
#        x = "Actual maximum lifespan (log10)",
#        y = "Predicted maximum lifespan (log10)") +
#   theme_bw() + labs(fill = "Order", color = "Order") + 
#   theme(
#     legend.position = "right",
#     legend.title = element_text(face = "bold"),
#     panel.grid.minor = element_blank(),
#     plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(size = 11),
#     axis.title = element_text(size = 12)
#   )
# 
# # Coefficient plot for GENES+MASS model
# direction_colors <- c("Positive impact on longevity" = "#9b383a", 
#                       "Negative impact on longevity" = "#86acb9")
# 
# # Create enhanced bar plot with directionality
# p_predictors <- ggplot(top_important_predictors,
#                        aes(x = reorder(Predictor, Effect), y = Effect,
#                            fill = Direction)) +
#   geom_col() +
#   scale_fill_manual(values = direction_colors) +
#   coord_flip() +
#   labs(
#     title = "Top predictors in combined copy number and mass model",
#     x = "Predictor",
#     y = "Effect on longevity (model coefficient)"
#   ) +
#   theme_bw() +
#   theme(
#     legend.position = "top",
#     legend.title = element_text(face = "bold"),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(size = 11),
#     axis.title = element_text(size = 12)
#   )
# 
# print(p_comparison)
# print(p_predictors)
# 
# ggsave("/Users/katiarenault/PhD/TOGA/revised_results/plots/model_comparison_plot.png", p_comparison, width = 12, height = 6)
# ggsave("/Users/katiarenault/PhD/TOGA/revised_results/plots/top_predictors_with_mass.png", p_predictors, width = 8, height = 10)
# 
# # Save the model comparison results
# write.csv(comparison_df, "/Users/katiarenault/PhD/TOGA/revised_results/results/model_comparison_metrics.csv", row.names = FALSE)
# write.csv(top_important_predictors, "/Users/katiarenault/PhD/TOGA/revised_results/results/top_predictors_with_mass.csv", row.names = FALSE)
# 
# 
# ##############################################
# ## Model predicting based on body mass only ##
# ##############################################
# 
# library(caret)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(ggrepel)
# library(tibble)     
# library(viridis)  
# library(RColorBrewer)
# 
# metadata <- read.csv("/Users/katiarenault/PhD/TOGA/revised_results/data/raxml_final_metadata_revised.csv")
# 
# species_data <- metadata %>%
#   filter(!is.na(maximum_longevity_y), !is.na(adult_body_mass_g)) %>%
#   mutate(log_longevity = log10(maximum_longevity_y),
#          log_mass = log10(adult_body_mass_g))
# print(paste("Number of species with complete data:", nrow(species_data)))
# 
# species_data$longevity_bin <- cut(species_data$log_longevity,
#                                   breaks = quantile(species_data$log_longevity,
#                                                     probs = seq(0, 1, 0.2)),
#                                   include.lowest = TRUE)
# bin_counts <- table(species_data$longevity_bin)
# print("Distribution of species across longevity bins:")
# print(bin_counts)
# 
# # calculate correlation between log mass and log longevity
# mass_longevity_cor <- cor(species_data$log_mass, species_data$log_longevity)
# mass_longevity_rho <- cor(species_data$log_mass, species_data$log_longevity, 
#                           method = "spearman")
# 
# print(paste("Pearson correlation (r) between log mass and log longevity:", 
#             round(mass_longevity_cor, 3)))
# print(paste("Spearman correlation (rho) between log mass and log longevity:", 
#             round(mass_longevity_rho, 3)))
# 
# # create a simple linear model to see the relationship
# initial_lm <- lm(log_longevity ~ log_mass, data = species_data)
# summary(initial_lm)
# 
# #########################################################
# # Set up train-test split with same seed as gene model ##
# #########################################################
# set.seed(3411) 
# # create a stratified train-test split (80/20)
# train_index <- createDataPartition(species_data$longevity_bin,
#                                    p = 0.8,
#                                    list = FALSE)
# 
# # create train and test sets
# train_data <- species_data[train_index, ]
# test_data <- species_data[-train_index, ]
# 
# print(paste("Training set size:", nrow(train_data)))
# print(paste("Test set size:", nrow(test_data)))
# 
# # create mass variables in both train and test sets
# train_data <- train_data %>%
#   mutate(mass_linear = adult_body_mass_g,
#          mass_log = log_mass)
# 
# test_data <- test_data %>%
#   mutate(mass_linear = adult_body_mass_g,
#          mass_log = log_mass)
# 
# # create model matrix with both predictors
# x_train <- model.matrix(~ mass_linear + mass_log, train_data)[,-1] # Remove intercept
# y_train <- train_data$log_longevity
# x_test <- model.matrix(~ mass_linear + mass_log, test_data)[,-1]
# y_test <- test_data$log_longevity
# 
# # scale features
# x_train_scaled <- scale(x_train)
# x_test_scaled <- scale(x_test, 
#                        center = attr(x_train_scaled, "scaled:center"),
#                        scale = attr(x_train_scaled, "scaled:scale"))
# 
# # elastic net tuning
# alphas <- seq(0, 1, 0.1)
# tune_results <- map_df(alphas, function(a) {
#   cv <- cv.glmnet(x_train_scaled, y_train, alpha = a, nfolds = 10)
#   tibble(alpha = a,
#          lambda_min = cv$lambda.min,
#          mse_min = min(cv$cvm))
# })
# 
# # get best model
# best_alpha <- tune_results$alpha[which.min(tune_results$mse_min)]
# best_model <- glmnet(x_train_scaled, y_train, 
#                      alpha = best_alpha,
#                      lambda = tune_results$lambda_min[which.min(tune_results$mse_min)])
# 
# # predictions
# test_data$predicted <- predict(best_model, newx = x_test_scaled)[,1]  # Extract column vector
# 
# # performance metrics
# metrics <- list(
#   MAE = MAE(test_data$predicted, test_data$log_longevity),
#   R2 = R2(test_data$predicted, test_data$log_longevity),
#   Pearson = cor(test_data$predicted, test_data$log_longevity)
# )
# 
# # Create color mapping (assuming you have this)
# orders <- unique(test_data$order)
# sophisticated_jewels_palette <- c(
#   "#9b383a", "#64a590", "#AA4839", "#AF6C36", "#945a87",
#   "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
#   "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
#   "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
#   "#c07a7a", "#8a7c6d", "#5c4d8b", "#4878a0", "#3D7B68"
# )
# 
# if (length(orders) > length(sophisticated_jewels_palette)) {
#   # If more orders than colors, recycle colors
#   color_mapping <- setNames(
#     sophisticated_jewels_palette[1:length(sophisticated_jewels_palette)], 
#     orders[1:length(sophisticated_jewels_palette)]
#   )
#   warning("More orders than available colors. Some colors will be recycled.")
# } else {
#   # Use just enough colors from the palette
#   color_mapping <- setNames(sophisticated_jewels_palette[1:length(orders)], orders)
# }
# color_df <- data.frame(
#   order = names(color_mapping),
#   color = color_mapping
# )
# 
# p_test <- ggplot(test_data, aes(x = log_longevity, y = predicted, color = order)) +
#   # Use shape = 21 for points with fill and color
#   geom_point(size = 4, alpha = 0.7, stroke = 0.8, shape = 21, aes(fill = order)) +
#   scale_fill_manual(values = color_mapping) +
#   scale_color_manual(values = color_mapping) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
#   labs(title = "Test Set Performance",
#        subtitle = paste0("R² = ", round(metrics$R2, 3), 
#                          ", ρ = ", round(metrics$Pearson, 3),
#                          ", MAE = ", round(metrics$MA, 3)),
#        x = "Actual maximum lifespan (log10)",
#        y = "Predicted maximum lifespan (log10)") +
#   theme_bw() + labs(fill = "Order", color = "Order") +
#   theme(
#     legend.position = "right",
#     legend.title = element_text(face = "bold", size = 11),  # Style
#     guides(fill = guide_legend(title = "Order"),            # Title text for fill
#            color = guide_legend(title = "Order")),                 # Title text for color
#     panel.grid.minor = element_blank(),
#     plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(size = 11, hjust = 0.5),
#     axis.title = element_text(size = 12)
#   )
# 
# p_test
# 
# ggsave("/Users/katiarenault/PhD/TOGA/revised_results/plots/only_mass_model_longevity_prediction.png", p_test, width = 10, height = 8)

# ###################################################################################
# pan_column <- NULL
# possible_names <- c("Pan_troglodytes", "Pan troglodytes")
# for (name in possible_names) {
#   if (name %in% colnames(gene_copy_data)) {
#     pan_column <- name
#     cat("Found Pan troglodytes as:", pan_column, "\n")
#     break
#   }
# }
# if (is.null(pan_column)) {
#   cat("Warning: Pan_troglodytes not found in data columns!\n")
#   cat("Available species (first 10):", head(colnames(gene_copy_data), 10), "...\n")
#   possible_match <- colnames(gene_copy_data)[grep("Pan", colnames(gene_copy_data))]
#   if (length(possible_match) > 0) {
#     pan_column <- possible_match[1]
#     cat("Using closest match instead:", pan_column, "\n")
#   } else {
#     stop("Cannot proceed: Pan_troglodytes not found in dataset")
#   }
# }
# has_higher_copies_than_pan <- function(gene_row, pan_column) {
#   pan_copies <- gene_row[pan_column]
#   any_higher <- any(gene_row > pan_copies, na.rm = TRUE)
#   return(any_higher)
# }
# 
# higher_than_pan_genes <- apply(gene_copy_data, 1, function(row) {
#   has_higher_copies_than_pan(row, pan_column)
# })
# gene_copy_data <- gene_copy_data[higher_than_pan_genes, ]
#gene_copy_data <- log(gene_copy_data + 0.01)
###################################################################################
