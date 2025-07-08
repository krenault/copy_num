############################
## TOGA copy number count ##
############################


#####################
# Transcript count ##
#####################

library(dplyr)
library(tidyr)
library(ape)
library(phangorn)
library(R.utils)  # For handling gzipped files

base_dir <- "/Users/katiarenault/PhD/TOGA/loss"

# Function to find gzipped TSV files
find_tsv_files <- function(dir) {
  all_files <- list.files(dir, pattern = "\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)
  return(all_files)
}

# Function to extract species name from directory path
extract_species_name <- function(file_path) {
  # Split the path by "/"
  parts <- unlist(strsplit(file_path, "/"))

  # Find the part containing the species name (should be before the file name)
  for (i in length(parts):1) {
    # Check if this part contains the species name pattern (typically has underscores)
    if (grepl("__", parts[i])) {
      # Extract the species name from pattern like "Uropsilus_gracilis__gracile_shrew_mole__HLuroGra1"
      species_name <- unlist(strsplit(parts[i], "__"))[1]
      return(species_name)
    }
  }

  # If not found using the pattern above, try another approach
  for (i in length(parts):1) {
    if (grepl("_", parts[i]) && !grepl("\\.tsv", parts[i])) {
      return(parts[i])
    }
  }

  # Default fallback
  return("unknown_species")
}

tsv_files <- find_tsv_files(base_dir)
loss_data_files <- tsv_files[grepl("loss_summ_data.tsv.gz$", tsv_files)]

# Initialize empty list of species
sp <- character()
species_count <- list()

# Process each file
for(file_path in loss_data_files) {
  # Extract species name from directory path
  scientific_name <- extract_species_name(file_path)

  # Handle multiple files for the same species
  if (scientific_name %in% names(species_count)) {
    species_count[[scientific_name]] <- species_count[[scientific_name]] + 1
    species_with_count <- paste0(scientific_name, "_", species_count[[scientific_name]])
  } else {
    species_count[[scientific_name]] <- 1
    species_with_count <- scientific_name
  }

  # Add to species list
  sp <- c(sp, species_with_count)

  # Create a temporary file to unzip
  temp_file <- tempfile(fileext = ".tsv")

  # Unzip the file
  tryCatch({
    R.utils::gunzip(file_path, destname = temp_file, remove = FALSE)

    # Read the unzipped file
    tab <- read.table(temp_file)
    colnames(tab) <- c("Type", "ID", "State")
    assign(species_with_count, tab)
    print(paste("Processed:", species_with_count, "from file:", basename(file_path)))

    # Clean up the temporary file
    file.remove(temp_file)
  }, error = function(e) {
    warning(paste("Error processing file:", file_path, "- Error:", e$message))
    if (file.exists(temp_file)) {
      file.remove(temp_file)
    }
  })
}

# Create a list of transcripts for each species
transcript_dfs <- list()
for(i in 1:length(sp)) {
  print(paste("Processing transcripts for:", sp[i]))
  # Check if the species data frame exists
  if (exists(sp[i])) {
    tab1 <- get(sp[i])
    transcript_only <- tab1[tab1$Type=="TRANSCRIPT", c(2,3)]
    colnames(transcript_only) <- c("ID", sp[i])
    assign(paste0(sp[i], "_transcript"), transcript_only)
    # Add to the list of data frames
    transcript_dfs[[i]] <- transcript_only
  } else {
    warning(paste("Data frame for", sp[i], "does not exist. Skipping."))
  }
}

# Create the list of transcript dataframes
dfs <- lapply(sp, function(species) {
  transcript_df_name <- paste0(species, "_transcript")
  if (exists(transcript_df_name)) {
    return(get(transcript_df_name))
  } else {
    return(NULL)
  }
})

# Remove NULL entries from dfs
dfs <- dfs[!sapply(dfs, is.null)]

# Merge all transcript data frames
if (length(dfs) > 0) {
  all.species.transcripts <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), dfs)
  all.sp.transcripts <- all.species.transcripts[, -1]
  row.names(all.sp.transcripts) <- all.species.transcripts$ID

  # Which of these transcripts are intact or partially intact in all species?
  # Replace NAs with a placeholder that won't match the intact conditions
  all.sp.transcripts[is.na(all.sp.transcripts)] <- "NA"
  selected_transcripts <- all.sp.transcripts[apply(all.sp.transcripts, 1, function(x) all(x %in% c("PI")))] #, "PI"))), ]

  # Create output file
  output_file <- file.path(base_dir, "Mammalia_Intact_transcripts.tsv")
  write.table(selected_transcripts, output_file, sep = '\t')
  print(paste("Wrote intact transcripts to:", output_file))
} else {
  warning("No transcript data frames available to merge.")
}

######################### Copy number count #########################
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

# Create the list of projection dataframes
dfp <- lapply(sp, function(species) {
  projection_df_name <- paste0(species, "_projections")
  if (exists(projection_df_name)) {
    return(get(projection_df_name))
  } else {
    return(NULL)
  }
})

# Remove NULL entries from dfp
dfp <- dfp[!sapply(dfp, is.null)]

# Process each species
copynumber_dfs <- list()
for(i in 1:length(dfp)) {
  species_name <- sp[i]
  print(paste("Processing copy numbers for:", species_name))

  df_projections <- dfp[[i]]
  df.sp <- df_projections %>% tidyr::separate(ID, c('TranscriptID', 'Gene', 'Chain'), sep = "[.]")

  # Keep only I and PI
  df_projections_intact <- df.sp[df.sp[, ncol(df.sp)]=="I" | df.sp[, ncol(df.sp)]=="PI", ]
  copy_number <- as.data.frame(table(df_projections_intact$Gene))
  colnames(copy_number) <- c("Gene", species_name)

  assign(paste0(species_name, "_copynumber"), copy_number)
  copynumber_dfs[[i]] <- copy_number
}

# Create the list of copy number dataframes
dfcn <- lapply(sp, function(species) {
  copynumber_df_name <- paste0(species, "_copynumber")
  if (exists(copynumber_df_name)) {
    return(get(copynumber_df_name))
  } else {
    return(NULL)
  }
})

# Remove NULL entries from dfcn
dfcn <- dfcn[!sapply(dfcn, is.null)]

# Merge all data frames
if (length(dfcn) > 0) {
  all.species.copynumber <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), dfcn)
  all.species.copynumber[is.na(all.species.copynumber)] <- 0

  output_dir <- file.path(base_dir)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Save the copy number table
  output_file <- file.path(output_dir, "Mammalia_all_species_gene_copynumber.csv")
  write.csv(all.species.copynumber, output_file)
  print(paste("Wrote copy number table to:", output_file))
} else {
  warning("No copy number data frames available to merge.")
}
###############################
# # Assuming your source dataframe is called 'source_df' and target is 'target_df'
# # Define the columns you want to extract
# columns_to_extract <- c(
#   "Ceratotherium_simum_cottoni", "Ceratotherium_simum_cottoni_2",
#   "Ceratotherium_simum_simum", "Dicerorhinus_sumatrensis_sumatrensis",
#   "Diceros_bicornis", "Diceros_bicornis_minor", "Equus_asinus",
#   "Equus_asinus_asinus", "Equus_burchellii_boehmi", "Equus_burchellii_quagga",
#   "Equus_caballus", "Equus_przewalskii", "Equus_zebra",
#   "Rhinoceros_unicornis", "Tapirus_indicus", "Tapirus_indicus_2",
#   "Tapirus_terrestris", "Tapirus_terrestris_2"
# )
# source_df <- read.csv("/Users/katiarenault/Desktop/PhD/TOGA/results/Aves_All_species_copynumber.tsv", sep = '\t')
# target_df <- read.csv("/Users/katiarenault/Desktop/PhD/TOGA/results/Mammalian_All_species_copynumber.tsv", sep = '\t')
# # Check which columns actually exist in the source dataframe
# existing_columns <- columns_to_extract[columns_to_extract %in% colnames(source_df)]
# 
# if (length(existing_columns) > 0) {
#   # Extract the subset of data
#   extracted_df <- source_df[, existing_columns, drop = FALSE]
#   
#   # Add rownames as a column for merging
#   extracted_df$Gene <- source_df$Gene
#   
#   # Convert target_df's rownames to a column for merging
#   
#   # Merge the dataframes
#   target_df <- merge(target_df, extracted_df, by = "Gene", all.x = TRUE)
#   source_df <- source_df[, !colnames(source_df) %in% existing_columns]
#   
#   cat("Transferred", length(existing_columns), "columns to the target dataframe\n")
# } else {
#   cat("None of the specified columns were found in the source dataframe\n")
# }
# write.table(target_df, "/Users/katiarenault/Desktop/PhD/TOGA/results/Mammalian_All_species_copynumber_updated.tsv", sep = '\t', row.names = FALSE)
# write.table(source_df, "/Users/katiarenault/Desktop/PhD/TOGA/results/Aves_All_species_copynumber_updated.tsv", sep = '\t', row.names = FALSE)

###############################

###############################
#### Phylogenetic analysis ####
###############################

#   mat <- as.matrix(t(cbind(all.species.copynumber[,-1], human=1)))
#   dist <- dist(mat)
#   
#   # NJ tree
#   nj_tree <- nj(dist)
#   pdf(file.path(output_dir, "nj_tree.pdf"))
#   plot(nj_tree, main="Neighbor-joining tree")
#   dev.off()
#   
#   # UPGMA tree
#   upgma <- upgma(dist)
#   pdf(file.path(output_dir, "upgma_tree.pdf"))
#   plot(upgma, main="UPGMA tree")
#   dev.off()
#   
#   # Maximum parsimony analysis
#   dat <- phyDat(as.matrix(t(cbind(all.species.copynumber[,-1], human=1))), type="USER", levels = 0:max(all.species.copynumber[,-1], na.rm=TRUE))
#   
#   # Try to build the parsimony tree, with error handling
#   tryCatch({
#     mm_start <- random.addition(dat)
#     mm_tree <- pratchet(dat, start = mm_start, minit = 1000, maxit = 10000,
#                         all = TRUE, trace = 0)
#     mm_tree <- acctran(mm_tree, dat)
#     mm_cons <- consensus(mm_tree)
#     
#     pdf(file.path(output_dir, "parsimony_tree.pdf"))
#     plot(mm_cons, main="Unrooted pratchet consensus tree")
#     dev.off()
#   }, error = function(e) {
#     warning("Error in maximum parsimony analysis: ", e$message)
#   })
# } else {
#   warning("No copy number data frames available to merge.")
# }

# print("Analysis complete!")

###################
## TOGA paralogs ##
###################
# 
# library(dplyr)
# library(tidyr)
# library(ape)
# library(phangorn)
# library(R.utils)  # For handling gzipped files
# 
# base_dir <- "/Users/katiarenault/PhD/TOGA/loss"
# 
# # Function to find gzipped TSV files
# find_tsv_files <- function(dir) {
#   all_files <- list.files(dir, pattern = "\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)
#   return(all_files)
# }
# 
# # Function to extract species name from directory path
# extract_species_name <- function(file_path) {
#   # Split the path by "/"
#   parts <- unlist(strsplit(file_path, "/"))
#   
#   # Find the part containing the species name (should be before the file name)
#   for (i in length(parts):1) {
#     # Check if this part contains the species name pattern (typically has underscores)
#     if (grepl("__", parts[i])) {
#       # Extract the species name from pattern like "Uropsilus_gracilis__gracile_shrew_mole__HLuroGra1"
#       species_name <- unlist(strsplit(parts[i], "__"))[1]
#       return(species_name)
#     }
#   }
#   
#   # If not found using the pattern above, try another approach
#   for (i in length(parts):1) {
#     if (grepl("_", parts[i]) && !grepl("\\.tsv", parts[i])) {
#       return(parts[i])
#     }
#   }
#   
#   # Default fallback
#   return("unknown_species")
# }
# 
# tsv_files <- find_tsv_files(base_dir)
# loss_data_files <- tsv_files[grepl("loss_summ_data.tsv.gz$", tsv_files)]
# 
# # Initialize empty list of species
# sp <- character()
# species_count <- list()
# 
# # Process each file
# for(file_path in loss_data_files) {
#   # Extract species name from directory path
#   scientific_name <- extract_species_name(file_path)
#   
#   # Handle multiple files for the same species
#   if (scientific_name %in% names(species_count)) {
#     species_count[[scientific_name]] <- species_count[[scientific_name]] + 1
#     species_with_count <- paste0(scientific_name, "_", species_count[[scientific_name]])
#   } else {
#     species_count[[scientific_name]] <- 1
#     species_with_count <- scientific_name
#   }
#   
#   # Add to species list
#   sp <- c(sp, species_with_count)
#   
#   # Create a temporary file to unzip
#   temp_file <- tempfile(fileext = ".tsv")
#   
#   # Unzip the file
#   tryCatch({
#     R.utils::gunzip(file_path, destname = temp_file, remove = FALSE)
#     
#     # Read the unzipped file
#     tab <- read.table(temp_file)
#     colnames(tab) <- c("Type", "ID", "State")
#     assign(species_with_count, tab)
#     print(paste("Processed:", species_with_count, "from file:", basename(file_path)))
#     
#     # Clean up the temporary file
#     file.remove(temp_file)
#   }, error = function(e) {
#     warning(paste("Error processing file:", file_path, "- Error:", e$message))
#     if (file.exists(temp_file)) {
#       file.remove(temp_file)
#     }
#   })
# }
# 
# # Create a list of transcripts for each species
# transcript_dfs <- list()
# for(i in 1:length(sp)) {
#   print(paste("Processing transcripts for:", sp[i]))
#   # Check if the species data frame exists
#   if (exists(sp[i])) {
#     tab1 <- get(sp[i])
#     transcript_only <- tab1[tab1$Type=="TRANSCRIPT", c(2,3)]
#     colnames(transcript_only) <- c("ID", sp[i])
#     assign(paste0(sp[i], "_transcript"), transcript_only)
#     # Add to the list of data frames
#     transcript_dfs[[i]] <- transcript_only
#   } else {
#     warning(paste("Data frame for", sp[i], "does not exist. Skipping."))
#   }
# }
# 
# # Create the list of transcript dataframes
# dfs <- lapply(sp, function(species) {
#   transcript_df_name <- paste0(species, "_transcript")
#   if (exists(transcript_df_name)) {
#     return(get(transcript_df_name))
#   } else {
#     return(NULL)
#   }
# })
# 
# # Remove NULL entries from dfs
# dfs <- dfs[!sapply(dfs, is.null)]
# 
# # Merge all transcript data frames
# if (length(dfs) > 0) {
#   all.species.transcripts <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), dfs)
#   all.sp.transcripts <- all.species.transcripts[, -1]
#   row.names(all.sp.transcripts) <- all.species.transcripts$ID
#   
#   # Which of these transcripts are intact or partially intact in all species?
#   # Replace NAs with a placeholder that won't match the intact conditions
#   all.sp.transcripts[is.na(all.sp.transcripts)] <- "NA"
#   selected_transcripts <- all.sp.transcripts[apply(all.sp.transcripts, 1, function(x) all(x %in% c("I"))), ]
#   
#   # Create output file
#   output_file <- file.path(base_dir, "Mammalia_Intact_transcripts.tsv")
#   write.table(selected_transcripts, output_file, sep = '\t')
#   print(paste("Wrote intact transcripts to:", output_file))
# } else {
#   warning("No transcript data frames available to merge.")
# }
# 
# ######################### Copy number count #########################
# projection_dfs <- list()
# for(i in 1:length(sp)) {
#   print(paste("Processing projections for:", sp[i]))
#   # Check if the species data frame exists
#   if (exists(sp[i])) {
#     tab1 <- get(sp[i])
#     projections_only <- tab1[tab1$Type=="PROJECTION", c(2,3)]
#     colnames(projections_only) <- c("ID", sp[i])
#     assign(paste0(sp[i], "_projections"), projections_only)
#     projection_dfs[[i]] <- projections_only
#   } else {
#     warning(paste("Data frame for", sp[i], "does not exist. Skipping."))
#   }
# }
# 
# # Create the list of projection dataframes
# dfp <- lapply(sp, function(species) {
#   projection_df_name <- paste0(species, "_projections")
#   if (exists(projection_df_name)) {
#     return(get(projection_df_name))
#   } else {
#     return(NULL)
#   }
# })
# 
# # Remove NULL entries from dfp
# dfp <- dfp[!sapply(dfp, is.null)]
# 
# # Process each species
# copynumber_dfs <- list()
# for(i in 1:length(dfp)) {
#   species_name <- sp[i]
#   print(paste("Processing copy numbers for:", species_name))
#   
#   df_projections <- dfp[[i]]
#   df.sp <- df_projections %>% tidyr::separate(ID, c('TranscriptID', 'Gene', 'Chain'), sep = "[.]")
#   
#   # Keep only I and PI
#   df_projections_intact <- df.sp[df.sp[, ncol(df.sp)]=="I", ] #| df.sp[, ncol(df.sp)]=="PI", ]
#   copy_number <- as.data.frame(table(df_projections_intact$Gene))
#   colnames(copy_number) <- c("Gene", species_name)
#   
#   assign(paste0(species_name, "_copynumber"), copy_number) 
#   copynumber_dfs[[i]] <- copy_number
# }
# 
# # Create the list of copy number dataframes
# dfcn <- lapply(sp, function(species) {
#   copynumber_df_name <- paste0(species, "_copynumber")
#   if (exists(copynumber_df_name)) {
#     return(get(copynumber_df_name))
#   } else {
#     return(NULL)
#   }
# })
# 
# # Remove NULL entries from dfcn
# dfcn <- dfcn[!sapply(dfcn, is.null)]
# 
# # Merge all data frames
# if (length(dfcn) > 0) {
#   all.species.copynumber <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), dfcn)
#   all.species.copynumber[is.na(all.species.copynumber)] <- 0
#   
#   output_dir <- file.path(base_dir)
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir)
#   }
#   
#   # Save the copy number table
#   output_file <- file.path(output_dir, "Mammalia_all_species_intact_copynumber.csv")
#   write.csv(all.species.copynumber, output_file)
#   print(paste("Wrote copy number table to:", output_file))
# } else {
#   warning("No copy number data frames available to merge.")
# }