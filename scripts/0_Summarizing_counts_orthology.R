## Katia Renault
## Summarizing counts of genes to identify genes whose copy number correlates with longevity of species

########################################################################################################
# 1. Summarizing the ortholog counts
########################################################################################################
# Using the downloaded files from the base directory to assemble a df of gene counts for each species ##
# by counting the number of unique q_gene corresponding to a single t_gene
# Output: produces a list of dfs of counts of unique t_genes from each species
########################################################################################################

library(dplyr)
library(tidyr)
library(R.utils)
base_dir <- "/Users/katiarenault/PhD/TOGA/orthology"
find_ortholog_files <- function(dir) {
  list.files(dir, pattern = "orthologsClassification.tsv.gz$", full.names = TRUE, recursive = TRUE)
}
# extract species name from file path using "__" convention
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


########################################################################################################
# 2. Add gene symbols
########################################################################################################
# merging in the gene symbols for each species
# Output: gene counts with symbols for each species (All_Species_Orthologous_CopyNumber_Annotated.tsv) 
########################################################################################################

# re-scan ortholog files to extract gene symbol mapping
ortholog_files <- list.files(base_dir, pattern = "orthologsClassification.tsv.gz$", full.names = TRUE, recursive = TRUE)
gene_symbol_map <- data.frame()
for (file_path in ortholog_files) {
  temp_file <- tempfile(fileext = ".tsv")
  tryCatch({
    R.utils::gunzip(file_path, destname = temp_file, remove = FALSE)
    tab <- read.table(temp_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(tab) <- tolower(colnames(tab))
    # extract gene symbol from t_transcript
    if (all(c("t_gene", "t_transcript") %in% colnames(tab))) {
      tab$t_symbol <- sub(".*\\.", "", tab$t_transcript)  # extract portion after last "."
      # reduce to unique mapping
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
# remove duplicates across files by keeping the first occurrence
gene_symbol_map <- gene_symbol_map[!duplicated(gene_symbol_map$t_gene), ]
# load the original merged table
merged_df_path <- file.path(base_dir, "All_Species_Orthologous_CopyNumber.tsv")
merged_df <- read.table(merged_df_path, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# merge in the gene symbols
annotated_df <- merge(gene_symbol_map, merged_df, by = "t_gene", all.y = TRUE)
annotated_df <- annotated_df[, c("t_gene", "t_symbol", setdiff(colnames(annotated_df), c("t_gene", "t_symbol")))]
# write out the annotated version
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


#########################################################################################################
# 3. Identifying intact transcripts
#########################################################################################################
# not all genes may have intact transcripts, so cross-checking with loss_sum files to identify if there #
# are at least one intact transcript or partially intact transcript per gene
# Output: Mammalia_all_species_copynumber.csv with count of intact or partially intact transcripts from #
# a given gene
#########################################################################################################

library(dplyr)
library(tidyr)
library(ape)
library(phangorn)
library(R.utils)
base_dir <- "/Users/katiarenault/PhD/TOGA/loss"
# function to find gzipped TSV files
find_tsv_files <- function(dir) {
  all_files <- list.files(dir, pattern = "\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)
  return(all_files)
}
# function to extract species name from directory path
extract_species_name <- function(file_path) {
  parts <- unlist(strsplit(file_path, "/"))
  # find the part containing the species name (should be before the file name)
  for (i in length(parts):1) {
    # Check if this part contains the species name pattern (typically has underscores)
    if (grepl("__", parts[i])) {
      # Extract the species name from pattern like "Uropsilus_gracilis__gracile_shrew_mole__HLuroGra1"
      species_name <- unlist(strsplit(parts[i], "__"))[1]
      return(species_name)
    }
  }
  # fallback
  return("unknown_species")
}
tsv_files <- find_tsv_files(base_dir)
loss_data_files <- tsv_files[grepl("loss_summ_data.tsv.gz$", tsv_files)]
# initialize empty list of species
sp <- character()
species_count <- list()
# process each file
for(file_path in loss_data_files) {
  # extract species name from directory path
  scientific_name <- extract_species_name(file_path)
  # handle multiple files for the same species
  if (scientific_name %in% names(species_count)) {
    species_count[[scientific_name]] <- species_count[[scientific_name]] + 1
    species_with_count <- paste0(scientific_name, "_", species_count[[scientific_name]])
  } else {
    species_count[[scientific_name]] <- 1
    species_with_count <- scientific_name
  }
  sp <- c(sp, species_with_count)
  # create a temporary file to unzip
  temp_file <- tempfile(fileext = ".tsv")
  # unzip the file
  tryCatch({
    R.utils::gunzip(file_path, destname = temp_file, remove = FALSE)
    # read the unzipped file
    tab <- read.table(temp_file)
    colnames(tab) <- c("Type", "ID", "State")
    assign(species_with_count, tab)
    print(paste("Processed:", species_with_count, "from file:", basename(file_path)))
    #  clean up the temporary file
    file.remove(temp_file)
  }, error = function(e) {
    warning(paste("Error processing file:", file_path, "- Error:", e$message))
    if (file.exists(temp_file)) {
      file.remove(temp_file)
    }
  })
}
# create a list of transcripts for each species
transcript_dfs <- list()
for(i in 1:length(sp)) {
  print(paste("Processing transcripts for:", sp[i]))
  if (exists(sp[i])) {
    tab1 <- get(sp[i])
    transcript_only <- tab1[tab1$Type=="TRANSCRIPT", c(2,3)]
    colnames(transcript_only) <- c("ID", sp[i])
    assign(paste0(sp[i], "_transcript"), transcript_only)
    # add to the list of data frames
    transcript_dfs[[i]] <- transcript_only
  } else {
    warning(paste("Data frame for", sp[i], "does not exist. Skipping."))
  }
}
# create the list of transcript dataframes
dfs <- lapply(sp, function(species) {
  transcript_df_name <- paste0(species, "_transcript")
  if (exists(transcript_df_name)) {
    return(get(transcript_df_name))
  } else {
    return(NULL)
  }
})

# remove NULL entries from dfs
dfs <- dfs[!sapply(dfs, is.null)]
# merge all transcript data frames
if (length(dfs) > 0) {
  all.species.transcripts <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), dfs)
  all.sp.transcripts <- all.species.transcripts[, -1]
  row.names(all.sp.transcripts) <- all.species.transcripts$ID
  all.sp.transcripts[is.na(all.sp.transcripts)] <- "NA"
  selected_transcripts <- all.sp.transcripts[apply(all.sp.transcripts, 1, function(x) all(x %in% c("I")))] #, "PI"))), ]
  # create output file
  #output_file <- file.path(base_dir, "Mammalia_Intact_transcripts.tsv")
  #write.table(selected_transcripts, output_file, sep = '\t')
  print(paste("Wrote intact transcripts to:", output_file))
} else {
  warning("No transcript data frames available to merge.")
}

projection_dfs <- list()
for(i in 1:length(sp)) {
  print(paste("Processing projections for:", sp[i]))
  # Check if the species data frame exists
  if (exists(sp[i])) {
    tab1 <- get(sp[i])
    projections_only <- tab1[tab1$Type=="PROJECTION", c(2,3)]
    colnames(projections_only) <- c("ID", sp[i])
    assign(paste0(sp[i], "_projections"), projections_only)
    projection_dfs[[i]] <- projections_only
  } else {
    warning(paste("Data frame for", sp[i], "does not exist. Skipping."))
  }
}
# create the list of projection dataframes
dfp <- lapply(sp, function(species) {
  projection_df_name <- paste0(species, "_projections")
  if (exists(projection_df_name)) {
    return(get(projection_df_name))
  } else {
    return(NULL)
  }
})
# remove NULL entries from dfp
dfp <- dfp[!sapply(dfp, is.null)]
# process each species
copynumber_dfs <- list()
for(i in 1:length(dfp)) {
  species_name <- sp[i]
  print(paste("Processing copy numbers for:", species_name))
  df_projections <- dfp[[i]]
  df.sp <- df_projections %>% tidyr::separate(ID, c('TranscriptID', 'Gene', 'Chain'), sep = "[.]")
  # keep only I and PI
  df_projections_intact <- df.sp[df.sp[, ncol(df.sp)]=="PI", ]# | df.sp[, ncol(df.sp)]=="PI", ]
  copy_number <- as.data.frame(table(df_projections_intact$Gene))
  colnames(copy_number) <- c("Gene", species_name)
  assign(paste0(species_name, "_copynumber"), copy_number)
  copynumber_dfs[[i]] <- copy_number
}
# create the list of copy number dataframes
dfcn <- lapply(sp, function(species) {
  copynumber_df_name <- paste0(species, "_copynumber")
  if (exists(copynumber_df_name)) {
    return(get(copynumber_df_name))
  } else {
    return(NULL)
  }
})
# remove NULL entries from dfcn
dfcn <- dfcn[!sapply(dfcn, is.null)]
# merge all data frames
if (length(dfcn) > 0) {
  all.species.copynumber <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), dfcn)
  all.species.copynumber[is.na(all.species.copynumber)] <- 0
  output_dir <- file.path(base_dir)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  # save the copy number table
  #output_file <- file.path(output_dir, "All_Species_Orthologous_Intactness_Annotated.csv")
  #write.csv(all.species.copynumber, output_file)
  print(paste("Wrote copy number table to:", output_file))
} else {
  warning("No copy number data frames available to merge.")
}



#########################################################################################################
# 4. Only keeping the orthologs with intact transcripts
#########################################################################################################
# replacing counts from gene copy number table with zero if there are no intact or partially 
# intact transcripts for that gene in the query species by using the intactness csv
# produces a revised All_Species_Orthologous_CopyNumber_Annotated.tsv file
#########################################################################################################

library(dplyr)
library(readr)
ortholog <- read.csv("/Users/katiarenault/PhD/TOGA/orthology/All_Species_Orthologous_CopyNumber_Annotated.tsv",
                     sep="\t", stringsAsFactors=FALSE)
copynumber <- read.csv("/Users/katiarenault/PhD/TOGA/loss/All_Species_Orthologous_Intactness_Annotated.csv",
                       stringsAsFactors=FALSE)
# rows: match ortholog$t_symbol to copynumber$Gene
common_genes <- intersect(ortholog$t_symbol, copynumber$Gene)
ortholog_sub <- ortholog[ortholog$t_symbol %in% common_genes, ]
copynumber_sub <- copynumber[copynumber$Gene %in% common_genes, ]
# reorder rows to be identical (based on gene order)
ortholog_sub <- ortholog_sub[match(common_genes, ortholog_sub$t_symbol), ]
copynumber_sub <- copynumber_sub[match(common_genes, copynumber_sub$Gene), ]
# species columns (exclude gene columns)
ortholog_species_cols <- setdiff(colnames(ortholog_sub), c("t_symbol", "t_gene", "t_transcript"))
copynumber_species_cols <- setdiff(colnames(copynumber_sub), "Gene")
# find common species columns
common_species <- intersect(ortholog_species_cols, copynumber_species_cols)
# for each species, replace ortholog value with 0 where intactness copy number == 0
for (species in common_species) {
  zero_mask <- copynumber_sub[[species]] == 0
  ortholog_sub[[species]][zero_mask] <- 0
}
ortholog_sub
output_file <- file.path(base_dir, "All_Species_Orthologous_CopyNumber_Annotated.tsv")
write.table(ortholog_sub, output_file, sep = '\t', row.names = FALSE, quote = FALSE)
