library(dplyr)
library(tidyr)
library(R.utils)

base_dir <- "/Users/katiarenault/PhD/TOGA/orthology"

# Function to find orthologsClassification.tsv.gz files
find_ortholog_files <- function(dir) {
  list.files(dir, pattern = "orthologsClassification.tsv.gz$", full.names = TRUE, recursive = TRUE)
}

# Extract species name from file path using "__" convention
extract_species_name <- function(file_path) {
  parts <- unlist(strsplit(file_path, "/"))
  for (i in length(parts):1) {
    if (grepl("__", parts[i])) {
      return(unlist(strsplit(parts[i], "__"))[1])
    }
  }
  for (i in length(parts):1) {
    if (grepl("_", parts[i]) && !grepl("\\.tsv", parts[i])) {
      return(parts[i])
    }
  }
  return("unknown_species")
}

# Get list of relevant files
ortholog_files <- find_ortholog_files(base_dir)
cat("Found ortholog files:\n")
print(ortholog_files)

species_count <- list()
sp <- character()
ortholog_dfs <- list()

for (file_path in ortholog_files) {
  scientific_name <- extract_species_name(file_path)
  print(paste("File:", basename(file_path), "Species:", scientific_name))
  
  if (scientific_name %in% names(species_count)) {
    species_count[[scientific_name]] <- species_count[[scientific_name]] + 1
    species_with_count <- paste0(scientific_name, "_", species_count[[scientific_name]])
  } else {
    species_count[[scientific_name]] <- 1
    species_with_count <- scientific_name
  }
  sp <- c(sp, species_with_count)
  
  temp_file <- tempfile(fileext = ".tsv")
  tryCatch({
    R.utils::gunzip(file_path, destname = temp_file, remove = FALSE)
    tab <- read.table(temp_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(tab) <- tolower(colnames(tab))
    
    print(paste("Read file with", nrow(tab), "rows and", ncol(tab), "columns"))
    print("Columns:")
    print(colnames(tab))
    
    counts <- tab %>%
      group_by(t_gene) %>%
      summarize(!!species_with_count := n_distinct(q_gene))
    
    print(paste("Number of rows in counts for", species_with_count, ":", nrow(counts)))
    
    if (nrow(counts) > 0) {
      ortholog_dfs[[species_with_count]] <- counts
      print(paste("Processed:", species_with_count))
    } else {
      warning(paste("No data to add for species:", species_with_count))
    }
    
    file.remove(temp_file)
  }, error = function(e) {
    warning(paste("Error processing file:", file_path, "- Error:", e$message))
    if (file.exists(temp_file)) file.remove(temp_file)
  })
}

####################################################################
###################### Add gene symbols ############################
####################################################################

# Re-scan ortholog files to extract gene symbol mapping
ortholog_files <- list.files(base_dir, pattern = "orthologsClassification.tsv.gz$", full.names = TRUE, recursive = TRUE)

gene_symbol_map <- data.frame()

for (file_path in ortholog_files) {
  temp_file <- tempfile(fileext = ".tsv")
  tryCatch({
    R.utils::gunzip(file_path, destname = temp_file, remove = FALSE)
    tab <- read.table(temp_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(tab) <- tolower(colnames(tab))
    
    # Extract gene symbol from t_transcript
    if (all(c("t_gene", "t_transcript") %in% colnames(tab))) {
      tab$t_symbol <- sub(".*\\.", "", tab$t_transcript)  # extract portion after last "."
      
      # Reduce to unique mapping
      gene_map <- tab[, c("t_gene", "t_symbol")]
      gene_map <- gene_map[!duplicated(gene_map$t_gene), ]
      
      gene_symbol_map <- rbind(gene_symbol_map, gene_map)
    }
    
    file.remove(temp_file)
  }, error = function(e) {
    warning(paste("Failed to process for mapping:", file_path, "-", e$message))
    if (file.exists(temp_file)) file.remove(temp_file)
  })
}

# Remove duplicates across files by keeping the first occurrence
gene_symbol_map <- gene_symbol_map[!duplicated(gene_symbol_map$t_gene), ]

# Load the original merged table
merged_df_path <- file.path(base_dir, "All_Species_Orthologous_CopyNumber.tsv")
merged_df <- read.table(merged_df_path, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Merge in the gene symbols
annotated_df <- merge(gene_symbol_map, merged_df, by = "t_gene", all.y = TRUE)

# Optional: move gene symbol column up front
annotated_df <- annotated_df[, c("t_gene", "t_symbol", setdiff(colnames(annotated_df), c("t_gene", "t_symbol")))]

# Write out the annotated version
output_file <- file.path(base_dir, "All_Species_Orthologous_CopyNumber_Annotated.tsv")
write.table(annotated_df, output_file, sep = '\t', row.names = FALSE, quote = FALSE)

cat("✅ Added gene symbols. Output written to:", output_file, "\n")


if (length(ortholog_dfs) > 0) {
  merged_df <- Reduce(function(x, y) full_join(x, y, by = "t_gene"), ortholog_dfs)
  merged_df[is.na(merged_df)] <- 0
  output_file <- file.path(base_dir, "All_Species_Orthologous_CopyNumber.tsv")
  write.table(merged_df, output_file, sep = '\t', row.names = FALSE, quote = FALSE)
  print(paste("Wrote orthologous gene copy counts to:", output_file))
} else {
  warning("No ortholog data frames available to merge.")
}

######################################################################################
################### Only keeping the orthologs with intact transcripts ###############
######################################################################################
library(dplyr)
library(readr)

ortholog <- read.csv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber_Annotated.tsv",
                     sep="\t", stringsAsFactors=FALSE)
copynumber <- read.csv("/Users/katiarenault/PhD/TOGA/loss/Mammalia_all_species_copynumber.csv",
                       stringsAsFactors=FALSE)

# Rows: match ortholog$t_symbol to copynumber$Gene
common_genes <- intersect(ortholog$t_symbol, copynumber$Gene)
ortholog_sub <- ortholog[ortholog$t_symbol %in% common_genes, ]
copynumber_sub <- copynumber[copynumber$Gene %in% common_genes, ]

# Reorder rows to be identical (based on gene order)
ortholog_sub <- ortholog_sub[match(common_genes, ortholog_sub$t_symbol), ]
copynumber_sub <- copynumber_sub[match(common_genes, copynumber_sub$Gene), ]

# Species columns (exclude gene columns)
ortholog_species_cols <- setdiff(colnames(ortholog_sub), c("t_symbol", "t_gene", "t_transcript"))
copynumber_species_cols <- setdiff(colnames(copynumber_sub), "Gene")

# Find common species columns
common_species <- intersect(ortholog_species_cols, copynumber_species_cols)

# For each species, replace ortholog value with 0 where copynumber == 0
for (species in common_species) {
  zero_mask <- copynumber_sub[[species]] == 0
  ortholog_sub[[species]][zero_mask] <- 0
}

ortholog_sub
output_file <- file.path(base_dir, "All_Species_Orthologous_CopyNumber_Annotated.tsv")
write.table(ortholog_sub, output_file, sep = '\t', row.names = FALSE, quote = FALSE)
