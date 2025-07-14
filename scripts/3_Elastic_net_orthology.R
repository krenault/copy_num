## Katia Renault
## Prediction of species lifespan based on gene copy number  

seed = 13
########################################################################################################
# 1. Elastic net model based on copy number
########################################################################################################

library(glmnet)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork) 
library(viridis)
library(RColorBrewer)
library(doParallel)
library(ggrepel)
library(tibble) 
#seed <- 2024
# 1258
# 3411
# 9207
set.seed(seed)

gene_copy_data <- read.csv("/Users/katiarenault/Documents/GitHub/copy_num/data/All_Species_Orthologous_CopyNumber_Annotated.tsv", row.names = "t_gene", sep = '\t')
rownames(gene_copy_data) <- gene_copy_data$t_symbol
gene_copy_data <- gene_copy_data[!grepl("^ZNF", rownames(gene_copy_data)), ]
gene_copy_data <- gene_copy_data[!grepl("^SPAN", rownames(gene_copy_data)), ]
gene_copy_data <- gene_copy_data %>% dplyr::select(-t_symbol)

metadata <- read.csv("/Users/katiarenault/Documents/GitHub/copy_num/data/raxml_final_metadata_revised.csv")
gene_copy_t <- as.data.frame(t(gene_copy_data))
gene_copy_t$species <- sub("_\\d+$", "", rownames(gene_copy_t))
gene_copy_by_species <- gene_copy_t %>%
  group_by(species) %>%
  filter(row_number() == 1) %>%
  ungroup()

species_data <- gene_copy_by_species %>%
  inner_join(metadata[, c("Scientific_name", "maximum_longevity_y", "order")], 
             by = c("species" = "Scientific_name")) %>%
  filter(!is.na(maximum_longevity_y)) %>%
  mutate(log_longevity = log(maximum_longevity_y))

print(paste("Number of species after merging:", nrow(species_data)))
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

###########################
#### Number of genes ######
###########################

cors <- cors[!is.na(cors)]
abs_cors <- abs(cors)
ranked_genes <- names(sort(abs_cors, decreasing = TRUE))
top_genes <- ranked_genes[1:min(20, length(ranked_genes))]  
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
alpha_values <- seq(0, 1, by = 0.05)  # More fine-grained search
results <- list()

# loop through alpha values
for (a in alpha_values) {
  print(paste("Fitting model with alpha =", a))
  
  # fit model with cross-validation
  set.seed(seed)  
  cv_fit <- cv.glmnet(
    x = x_train_scaled,
    y = y_train,
    alpha = a,
    nfolds = 10,
    type.measure = "mse",
    nlambda = 100  
  )
  # store results
  results[[as.character(a)]] <- list(
    alpha = a,
    cv_fit = cv_fit,
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se,
    cvm_min = min(cv_fit$cvm) 
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

########################################################################################################
# 2. Most predictive genes in the model
########################################################################################################

coef_matrix <- as.matrix(coef(final_model, s = best_lambda))
var_importance <- coef_matrix[-1, 1] 
names(var_importance) <- rownames(coef_matrix)[-1]  
selected_features <- names(var_importance)[var_importance != 0]
n_features_used <- length(selected_features)
n_features_input <- ncol(x_train)
print(paste("=== L1 REGULARIZATION RESULTS ==="))
print(paste("Input features:", n_features_input))
print(paste("Features actually used:", n_features_used))
print(paste("Features zeroed out by L1:", n_features_input - n_features_used))
print(paste("Feature usage rate:", round(n_features_used/n_features_input * 100, 1), "%"))
non_zero_features <- var_importance[var_importance != 0]
top_features_to_show <- min(17, length(non_zero_features))
top_genes_indices <- order(abs(var_importance), decreasing = TRUE)[1:top_features_to_show]

top_important_genes <- data.frame(
  Gene = names(var_importance)[top_genes_indices],
  Effect = var_importance[top_genes_indices],
  Direction = ifelse(var_importance[top_genes_indices] > 0, "Positive impact on longevity", "Negative impact on longevity")
)
all_genes_df <- data.frame(
  Gene = names(var_importance),
  Coefficient = var_importance,
  Abs_Coefficient = abs(var_importance),
  Direction = ifelse(var_importance > 0, "Positive impact on longevity",
                     ifelse(var_importance < 0, "Negative impact on longevity", "No impact"))
)

all_genes_sorted <- all_genes_df[order(all_genes_df$Abs_Coefficient, decreasing = TRUE), ]
all_genes_sorted$Rank <- 1:nrow(all_genes_sorted)
#write.csv(all_genes_sorted, paste0("/Users/katiarenault/Documents/GitHub/copy_num/results/all_predictive_genes_25_features_seed_", seed, ".csv"), row.names = FALSE)


test_plot_data <- data.frame(
  actual = y_test,
  predicted = as.numeric(test_preds),
  order = test_orders,
  Species = test_species
)

#### Color mapping ####
source("/Users/katiarenault/Documents/Github/copy_num/scripts/FUN_color_mappings.R")
color_mapping <- create_color_mapping(species_data$order)
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
       x = "Actual maximum lifespan (log)",
       y = "Predicted maximum lifespan (log)") +
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
cat("\n=== TOP SELECTED FEATURES ===\n")
top_selected <- sort(abs(non_zero_features), decreasing = TRUE)[1:min(20, length(non_zero_features))]
for (i in 1:length(top_selected)) {
  gene <- names(top_selected)[i]
  coef <- non_zero_features[gene]
  direction <- ifelse(coef > 0, "positive", "negative")
  cat(sprintf("%2d. %s: %.4f (%s effect)\n", i, gene, coef, direction))
}
#ggsave("/Users/katiarenault/Documents/GitHub/copy_num/plots/test_set_performance.png",
 #      p_test, width = 12, height = 8, dpi = 300)
#ggsave("/Users/katiarenault/Documents/GitHub/copy_num/plots/elastic_net_predictive_genes.png",
 #      p_genes, width = 12, height = 8, dpi = 300)

########################################################################################################
# 3. Including body mass as a predictor for lifespan in model
########################################################################################################

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

gene_copy_by_species <- gene_copy_t %>%
  group_by(species) %>%
  filter(row_number() == 1) %>% 
  ungroup()
species_data <- gene_copy_by_species %>%
  inner_join(metadata[, c("Scientific_name", "maximum_longevity_y", "order", "adult_body_mass_g")],
             by = c("species" = "Scientific_name")) %>%
  filter(!is.na(maximum_longevity_y)) %>%
  mutate(log_longevity = log(maximum_longevity_y),
         log_mass = log(adult_body_mass_g))
species_data$longevity_bin <- cut(species_data$log_longevity,
                                  breaks = quantile(species_data$log_longevity,
                                                    probs = seq(0, 1, 0.2)),
                                  include.lowest = TRUE)
valid_genes <- setdiff(colnames(species_data),
                       c("species", "maximum_longevity_y", "log_longevity", "longevity_bin",
                         "order", "adult_body_mass_g", "log_mass"))
print(paste("Total number of genes:", length(valid_genes)))
set.seed(seed)
nzv_results <- nearZeroVar(species_data[, valid_genes],
                           saveMetrics = TRUE,
                           freqCut = 99/1,
                           uniqueCut = 3)
nzv_genes <- rownames(nzv_results)[nzv_results$nzv]
valid_genes_filtered <- setdiff(valid_genes, nzv_genes)
print(paste("Genes after removing extreme near-zero variance:", length(valid_genes_filtered)))
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
cors <- cors[!is.na(cors)]
abs_cors <- abs(cors)
ranked_genes <- names(sort(abs_cors, decreasing = TRUE))
top_genes <- ranked_genes[1:min(20, length(ranked_genes))]
print(paste("Selected top", length(top_genes), "genes by correlation"))

########################################################################################################
# 4. Comparing genes only model and genes+mass model
########################################################################################################

species_with_mass <- species_data %>%
  filter(!is.na(adult_body_mass_g))
print(paste("Species with both longevity and mass data:", nrow(species_with_mass)))

# MODEL 1: Genes only
feature_matrix_genes_only <- species_with_mass %>%
  dplyr::select(all_of(top_genes)) %>%
  as.matrix()
rownames(feature_matrix_genes_only) <- species_with_mass$species

# MODEL 2: Genes + body mass
feature_matrix_with_mass <- cbind("Mass (log)" = species_with_mass$log_mass,
                                  feature_matrix_genes_only)
rownames(feature_matrix_with_mass) <- species_with_mass$species
print(paste("Genes-only matrix dimensions:", nrow(feature_matrix_genes_only), "x", ncol(feature_matrix_genes_only)))
print(paste("Genes+mass matrix dimensions:", nrow(feature_matrix_with_mass), "x", ncol(feature_matrix_with_mass)))

set.seed(seed) 
train_index <- createDataPartition(species_with_mass$longevity_bin,
                                   p = 0.8,
                                   list = FALSE)
x_train_genes <- feature_matrix_genes_only[train_index, ]
y_train <- species_with_mass$log_longevity[train_index]
x_test_genes <- feature_matrix_genes_only[-train_index, ]
y_test <- species_with_mass$log_longevity[-train_index]


x_train_mass <- feature_matrix_with_mass[train_index, ]
x_test_mass <- feature_matrix_with_mass[-train_index, ]


train_species <- species_with_mass$species[train_index]
test_species <- species_with_mass$species[-train_index]
train_orders <- species_with_mass$order[train_index]
test_orders <- species_with_mass$order[-train_index]
train_mass <- species_with_mass$log_mass[train_index]
test_mass <- species_with_mass$log_mass[-train_index]

print(paste("Training set size:", length(y_train)))
print(paste("Test set size:", length(y_test)))

##############################
# MODEL 1: GENES-ONLY MODEL  #
##############################
set.seed(seed) 
x_train_genes_scaled <- scale(x_train_genes)
x_test_genes_scaled <- scale(x_test_genes,
                             center = attr(x_train_genes_scaled, "scaled:center"),
                             scale = attr(x_train_genes_scaled, "scaled:scale"))
alpha_values <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
results_genes <- list()
for (a in alpha_values) {
  print(paste("GENES-ONLY Model: Fitting with alpha =", a))
  set.seed(seed)
  cv_fit <- cv.glmnet(
    x = x_train_genes_scaled,
    y = y_train,
    alpha = a,
    nfolds = 10,
    type.measure = "mse"
  )
  results_genes[[as.character(a)]] <- list(
    alpha = a,
    cv_fit = cv_fit,
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se,
    cvm_min = min(cv_fit$cvm)
  )
  print(paste("  Lambda min:", cv_fit$lambda.min))
  print(paste("  Lambda 1se:", cv_fit$lambda.1se))
  print(paste("  Minimum CV error:", min(cv_fit$cvm)))
}
best_alpha_genes <- NULL
best_error_genes <- Inf
best_lambda_genes <- NULL

for (a in names(results_genes)) {
  if (results_genes[[a]]$cvm_min < best_error_genes) {
    best_error_genes <- results_genes[[a]]$cvm_min
    best_alpha_genes <- results_genes[[a]]$alpha
    best_lambda_genes <- results_genes[[a]]$lambda_min
  }
}
print(paste("GENES-ONLY Model - Best alpha:", best_alpha_genes))
print(paste("GENES-ONLY Model - Best lambda:", best_lambda_genes))
print(paste("GENES-ONLY Model - Best CV error:", best_error_genes))
set.seed(seed) 
final_model_genes <- glmnet(
  x = x_train_genes_scaled,
  y = y_train,
  alpha = best_alpha_genes,
  lambda = best_lambda_genes
)
train_preds_genes <- predict(final_model_genes, newx = x_train_genes_scaled, s = best_lambda_genes)
test_preds_genes <- predict(final_model_genes, newx = x_test_genes_scaled, s = best_lambda_genes)

train_perf_genes <- postResample(pred = train_preds_genes, obs = y_train)
test_perf_genes <- postResample(pred = test_preds_genes, obs = y_test)
test_r2_genes <- test_perf_genes["Rsquared"]
test_mse_genes <- test_perf_genes["RMSE"]^2
test_mae_genes <- test_perf_genes["MAE"]
test_rho_genes <- cor(y_test, as.numeric(test_preds_genes), method = "spearman")

performance_df_genes <- data.frame(
  Metric = c("RMSE", "R-squared", "MAE"),
  Train = c(train_perf_genes["RMSE"], train_perf_genes["Rsquared"], train_perf_genes["MAE"]),
  Test = c(test_perf_genes["RMSE"], test_perf_genes["Rsquared"], test_perf_genes["MAE"])
)

print("GENES-ONLY Model Performance:")
print(performance_df_genes)

##############################
# MODEL 2: GENES+MASS MODEL  #
##############################

set.seed(seed)
x_train_mass_scaled <- scale(x_train_mass)
x_test_mass_scaled <- scale(x_test_mass,
                            center = attr(x_train_mass_scaled, "scaled:center"),
                            scale = attr(x_train_mass_scaled, "scaled:scale"))
results_mass <- list()

for (a in alpha_values) {
  print(paste("GENES+MASS Model: Fitting with alpha =", a))
  set.seed(seed)
  cv_fit <- cv.glmnet(
    x = x_train_mass_scaled,
    y = y_train,
    alpha = a,
    nfolds = 10,
    type.measure = "mse"
  )
  results_mass[[as.character(a)]] <- list(
    alpha = a,
    cv_fit = cv_fit,
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se,
    cvm_min = min(cv_fit$cvm)
  )
  print(paste("  Lambda min:", cv_fit$lambda.min))
  print(paste("  Lambda 1se:", cv_fit$lambda.1se))
  print(paste("  Minimum CV error:", min(cv_fit$cvm)))
}
best_alpha_mass <- NULL
best_error_mass <- Inf
best_lambda_mass <- NULL

for (a in names(results_mass)) {
  if (results_mass[[a]]$cvm_min < best_error_mass) {
    best_error_mass <- results_mass[[a]]$cvm_min
    best_alpha_mass <- results_mass[[a]]$alpha
    best_lambda_mass <- results_mass[[a]]$lambda_min
  }
}
print(paste("GENES+MASS Model - Best alpha:", best_alpha_mass))
print(paste("GENES+MASS Model - Best lambda:", best_lambda_mass))
print(paste("GENES+MASS Model - Best CV error:", best_error_mass))
set.seed(seed)
final_model_mass <- glmnet(
  x = x_train_mass_scaled,
  y = y_train,
  alpha = best_alpha_mass,
  lambda = best_lambda_mass
)
train_preds_mass <- predict(final_model_mass, newx = x_train_mass_scaled, s = best_lambda_mass)
test_preds_mass <- predict(final_model_mass, newx = x_test_mass_scaled, s = best_lambda_mass)

train_perf_mass <- postResample(pred = train_preds_mass, obs = y_train)
test_perf_mass <- postResample(pred = test_preds_mass, obs = y_test)
test_r2_mass <- test_perf_mass["Rsquared"]
test_mse_mass <- test_perf_mass["RMSE"]^2
test_mae_mass <- test_perf_mass["MAE"]
test_rho_mass <- cor(y_test, as.numeric(test_preds_mass), method = "spearman")

performance_df_mass <- data.frame(
  Metric = c("RMSE", "R-squared", "MAE"),
  Train = c(train_perf_mass["RMSE"], train_perf_mass["Rsquared"], train_perf_mass["MAE"]),
  Test = c(test_perf_mass["RMSE"], test_perf_mass["Rsquared"], test_perf_mass["MAE"])
)
print("GENES+MASS Model Performance:")
print(performance_df_mass)

################################
# Compare model performances   #
################################

comparison_df <- data.frame(
  Metric = c("RMSE", "R-squared", "MAE", "Spearman's ρ"),
  Genes_Only = c(test_perf_genes["RMSE"], test_perf_genes["Rsquared"],
                 test_perf_genes["MAE"], test_rho_genes),
  Genes_Plus_Mass = c(test_perf_mass["RMSE"], test_perf_mass["Rsquared"],
                      test_perf_mass["MAE"], test_rho_mass),
  Improvement = c(
    test_perf_genes["RMSE"] - test_perf_mass["RMSE"],
    test_perf_mass["Rsquared"] - test_perf_genes["Rsquared"],
    test_perf_genes["MAE"] - test_perf_mass["MAE"],
    test_rho_mass - test_rho_genes
  ),
  Percent_Change = c(
    (test_perf_genes["RMSE"] - test_perf_mass["RMSE"]) / test_perf_genes["RMSE"] * 100,
    (test_perf_mass["Rsquared"] - test_perf_genes["Rsquared"]) / test_perf_genes["Rsquared"] * 100,
    (test_perf_genes["MAE"] - test_perf_mass["MAE"]) / test_perf_genes["MAE"] * 100,
    (test_rho_mass - test_rho_genes) / test_rho_genes * 100
  )
)

print("Model Performance Comparison:")
print(comparison_df)

coef_matrix_mass <- as.matrix(coef(final_model_mass, s = best_lambda_mass))
var_importance_mass <- coef_matrix_mass[-1, 1]
names(var_importance_mass) <- rownames(coef_matrix_mass)[-1]
mass_importance <- var_importance_mass["Mass (log)"]
mass_rank <- which(names(sort(abs(var_importance_mass), decreasing = TRUE)) == "Mass (log)")
print(paste("Body mass coefficient:", mass_importance))
print(paste("Body mass rank among all predictors:", mass_rank, "out of", length(var_importance_mass)))
top_predictors_indices <- order(abs(var_importance_mass), decreasing = TRUE)[1:10]
top_important_predictors <- data.frame(
  Predictor = names(var_importance_mass)[top_predictors_indices],
  Effect = var_importance_mass[top_predictors_indices],
  Direction = ifelse(var_importance_mass[top_predictors_indices] > 0,
                     "Positive impact on longevity", "Negative impact on longevity")
)

test_plot_data <- data.frame(
  actual = rep(y_test, 2),
  predicted = c(as.numeric(test_preds_genes), as.numeric(test_preds_mass)),
  order = rep(test_orders, 2),
  Species = rep(test_species, 2),
  Model = rep(c("Genes", "Genes + Body Mass"), each = length(y_test))
)

p_comparison <- ggplot(test_plot_data, aes(x = actual, y = predicted, color = order)) +
  geom_point(size = 3.5, alpha = 0.7, stroke = 0.8, shape = 21, aes(fill = order)) +
  scale_fill_manual(values = color_mapping) +
  scale_color_manual(values = color_mapping) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  facet_wrap(~Model) +
  labs(title = "Model performance comparison - Test set",
       subtitle = paste0("Genes copy number only: R² = ", round(test_r2_genes, 3),
                         ", ρ = ", round(test_rho_genes, 3), "\n",
                         "Genes copy number + body mass: R² = ", round(test_r2_mass, 3),
                         ", ρ = ", round(test_rho_mass, 3)),
       x = "Actual maximum lifespan (log)",
       y = "Predicted maximum lifespan (log)") +
  theme_bw() + labs(fill = "Order", color = "Order") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

direction_colors <- c("Positive impact on longevity" = "#9b383a",
                      "Negative impact on longevity" = "#86acb9")
p_predictors <- ggplot(top_important_predictors,
                       aes(x = reorder(Predictor, Effect), y = Effect,
                           fill = Direction)) +
  geom_col() +
  scale_fill_manual(values = direction_colors) +
  coord_flip() +
  labs(
    title = "Top predictors in combined copy number and mass model",
    x = "Predictor",
    y = "Effect on longevity (model coefficient)"
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

print(p_comparison)
print(p_predictors)

#ggsave("/Users/katiarenault/Documents/GitHub/copy_num/plots/model_comparison_plot.png", p_comparison, width = 12, height = 6)
#ggsave("/Users/katiarenault/Documents/GitHub/copy_num/plots/top_predictors_with_mass.png", p_predictors, width = 8, height = 10)
#write.csv(comparison_df, "/Users/katiarenault/Documents/GitHub/copy_num/results/model_comparison_metrics.csv", row.names = FALSE)
#write.csv(top_important_predictors, "/Users/katiarenault/Documents/GitHub/copy_num/results/top_predictors_with_mass.csv", row.names = FALSE)
