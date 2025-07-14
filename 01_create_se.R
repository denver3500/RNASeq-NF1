library(dplyr)
library(tidyr)
library(SummarizedExperiment)

# Find all .tsv files in the root folder
tsv_files <- list.files(pattern = "*.tsv", full.names = TRUE)
print(paste("Found", length(tsv_files), "TSV files:"))
print(tsv_files)

# Read the first file to initialize the merged data
merged_data <- read.table(tsv_files[1], header = TRUE, sep = "\t")
print(paste("Initialized with file:", tsv_files[1]))
print(paste("Dimensions:", nrow(merged_data), "x", ncol(merged_data)))

# Read the second file to get the correct gene names (if needed)
if (length(tsv_files) > 1) {
  correct_gene_names <- read.table(tsv_files[2], header = TRUE, sep = "\t")
  
  # Replace the gene_name column in the merged data with the correct gene names
  merged_data$gene_name <- correct_gene_names$gene_name
  print("Updated gene names from second file")
}

# Iterate over the remaining files and merge them
if (length(tsv_files) > 1) {
  for (i in 2:length(tsv_files)) {
    file <- tsv_files[i]
    print(paste("Processing file:", file))
    
    data <- read.table(file, header = TRUE, sep = "\t")
    
    # Check if gene_id columns match
    if (!all(merged_data$gene_id == data$gene_id)) {
      mismatched_indices <- which(merged_data$gene_id != data$gene_id)
      mismatched_gene_ids <- merged_data$gene_id[mismatched_indices]
      stop(paste("Gene IDs do not match across files. Mismatched gene IDs:", paste(mismatched_gene_ids, collapse = ", ")))
    }
    
    # Check if gene_name columns match
    if (!all(merged_data$gene_name == data$gene_name)) {
      mismatched_indices <- which(merged_data$gene_name != data$gene_name)
      mismatched_gene_names <- merged_data$gene_name[mismatched_indices]
      stop(paste("Gene names do not match across files. Mismatched gene names:", paste(mismatched_gene_names, collapse = ", ")))
    }
    
    # Merge the data by adding the sample columns (excluding gene_id and gene_name columns)
    merged_data <- cbind(merged_data, data[, -c(1, 2)])
    print(paste("Merged data now has dimensions:", nrow(merged_data), "x", ncol(merged_data)))
  }
}

print("Successfully merged all TSV files")

# Read the metadata file
metadata <- read.csv("Simple_sample_map_10062023.csv")
print(paste("Metadata dimensions:", nrow(metadata), "x", ncol(metadata)))

# Create a new column by merging specimenID and aliquotID
metadata <- metadata %>%
  mutate(merged_id = paste0(gsub("-", ".", specimenID), ".", aliquotID))

print("Created merged_id column in metadata")

# Filter metadata for the specified conditions
filtered_metadata <- metadata %>%
  filter(tumorType %in% c("Malignant Peripheral Nerve Sheath Tumor", "Plexiform Neurofibroma"),
         assay == "Bulk RNA sequencing",
         tissue == "primary tumor")

print(paste("Filtered metadata dimensions:", nrow(filtered_metadata), "x", ncol(filtered_metadata)))
print("Sample distribution:")
print(table(filtered_metadata$tumorType))

# Get the list of merged_ids to keep
merged_ids_to_keep <- filtered_metadata$merged_id

# Subset the merged data to keep only the columns corresponding to the filtered merged_ids
columns_to_keep <- c("gene_id", "gene_name", intersect(merged_ids_to_keep, colnames(merged_data)))
filtered_expression_data <- merged_data[, columns_to_keep, drop = FALSE]

print(paste("Filtered expression data dimensions:", nrow(filtered_expression_data), "x", ncol(filtered_expression_data)))

# Store the gene annotation
gene_annotation <- filtered_expression_data %>% select(gene_id, gene_name)

# Remove the gene_name column from the expression matrix
expression_matrix <- filtered_expression_data %>% select(-gene_name)

# Set gene_id as row names
rownames(expression_matrix) <- expression_matrix$gene_id
expression_matrix <- expression_matrix %>% select(-gene_id)

# Ensure sample names match between the data and metadata
colnames(expression_matrix) <- gsub("-", ".", colnames(expression_matrix))
filtered_metadata$merged_id <- gsub("-", ".", filtered_metadata$merged_id)

# Filter metadata to include only the samples present in the expression matrix
filtered_metadata <- filtered_metadata %>% 
  filter(merged_id %in% colnames(expression_matrix)) %>%
  arrange(merged_id)

# Reorder expression matrix columns to match metadata order
expression_matrix <- expression_matrix[, filtered_metadata$merged_id, drop = FALSE]

# Ensure merged_id is the first column
filtered_metadata <- filtered_metadata %>% select(merged_id, everything())

# Validation checks
stopifnot(rownames(expression_matrix) == gene_annotation$gene_id)
stopifnot(colnames(expression_matrix) == filtered_metadata$merged_id)

print("Data validation passed")

# Create SummarizedExperiment object
se <- SummarizedExperiment(assays = list(counts = expression_matrix),
                           colData = filtered_metadata,
                           rowData = gene_annotation)

print("SummarizedExperiment object created:")
print(se)

# Save the SummarizedExperiment object
saveRDS(se, "summarized_experiment.rds")
print("SummarizedExperiment object saved to 'summarized_experiment.rds'")

# Print summary statistics
print("Summary:")
print(paste("Total genes:", nrow(se)))
print(paste("Total samples:", ncol(se)))
print("Sample distribution by tumor type:")
print(table(se$tumorType))
print("Sample distribution by batch:")
print(table(se$batch))