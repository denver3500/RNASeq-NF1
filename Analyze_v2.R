library(dplyr)
library(tidyr)
library(SummarizedExperiment)

# List of .tsv files I downloaded using Downnload.py. Those are private :)
files <- c("salmon.merged.gene_tpm_jh_batch1.tsv", "salmon.merged.gene_tpm_wu_batch1.tsv", "salmon.merged.gene_tpm_wu_batch1_additional.tsv", "salmon.merged.gene_tpm_wu_batch2.tsv")

# Read the first file to initialize the merged data
merged_data <- read.table(files[1], header = TRUE, sep = "\t")

# Read the second file to get the correct gene names. Because in one file gene names are not correct. AGD3 is AGD 3.00
correct_gene_names <- read.table(files[2], header = TRUE, sep = "\t")

# Replace the gene_name column in the merged data with the correct gene names
merged_data$gene_name <- correct_gene_names$gene_name

# Iterate over the remaining files and merge them
for (file in files[-1]) {
  data <- read.table(file, header = TRUE, sep = "\t")
  
  # Check if gene_id columns match
  if (!all(merged_data$gene_id == data$gene_id)) {
    mismatched_gene_ids <- merged_data$gene_id[merged_data$gene_id != data$gene_id]
    stop(paste("Gene IDs do not match across files. Mismatched gene IDs:", paste(mismatched_gene_ids, collapse = ", ")))
  }
  
  # Merge the data by adding the sample columns
  merged_data <- cbind(merged_data, data[,-c(1,2)])
}

# Read the metadata file
metadata <- read.csv("Simple_sample_map_10062023.csv")

# Filter the metadata for MPNST + Primary tumor only
filtered_metadata_mpnst <- metadata %>%
  filter(tumorType == "Malignant Peripheral Nerve Sheath Tumor",
         assay == "Bulk RNA sequencing",
         tissue == "primary tumor")

# Filter the metadata for PN + Primary tumor only
filtered_metadata_pn <- metadata %>%
  filter(tumorType == "Plexiform Neurofibroma",
         assay == "Bulk RNA sequencing",
         tissue == "primary tumor")

# Create a new column by merging specimenID and aliquotID for both filtered metadata. We need it because in datafile they are using this combined format
filtered_metadata_mpnst <- filtered_metadata_mpnst %>%
  mutate(merged_id = paste0(gsub("-", ".", specimenID), ".", aliquotID))

filtered_metadata_pn <- filtered_metadata_pn %>%
  mutate(merged_id = paste0(gsub("-", ".", specimenID), ".", aliquotID))

# Get the list of merged_ids to keep for both filtered metadata
merged_ids_to_keep_mpnst <- filtered_metadata_mpnst$merged_id
merged_ids_to_keep_pn <- filtered_metadata_pn$merged_id

# Subset the merged data to keep only the columns corresponding to the filtered merged_ids that are present. e.g - if a sample is not present in the data, we don't need it in the metadata
columns_to_keep_mpnst <- c("gene_id", "gene_name", intersect(merged_ids_to_keep_mpnst, colnames(merged_data)))
filtered_data_mpnst <- merged_data[, columns_to_keep_mpnst, drop = FALSE]

columns_to_keep_pn <- c("gene_id", "gene_name", intersect(merged_ids_to_keep_pn, colnames(merged_data)))
filtered_data_pn <- merged_data[, columns_to_keep_pn, drop = FALSE]

# Write the filtered data to new .tsv files
write.table(filtered_data_mpnst, "filtered_data_mpnst.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(filtered_data_pn, "filtered_data_pn.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

print("Filtered data written successfully to filtered_data_mpnst.tsv and filtered_data_pn.tsv")

# Read the expression matrices
filtered_data_mpnst <- read.table("filtered_data_mpnst.tsv", header = TRUE, sep = "\t")
filtered_data_pn <- read.table("filtered_data_pn.tsv", header = TRUE, sep = "\t")

# Join the expression matrices together
expression_matrix <- full_join(filtered_data_mpnst, filtered_data_pn, by = c("gene_id", "gene_name"))

# Store the gene_name column for later
gene_annotation <- expression_matrix %>% select(gene_id, gene_name)

# Remove the gene_name column from the expression matrix
expression_matrix <- expression_matrix %>% select(-gene_name)

# Set gene_id as row names
rownames(expression_matrix) <- expression_matrix$gene_id
expression_matrix <- expression_matrix %>% select(-gene_id)

# Save the expression matrix as a CSV file
write.csv(expression_matrix, "expression_matrix.csv", row.names = TRUE)

print("Expression matrix saved successfully to expression_matrix.csv")

# Read the metadata file
allmetadata <- read.csv("Simple_sample_map_10062023.csv")

# Filter the metadata for "Malignant Peripheral Nerve Sheath Tumor"
filtered_metadata_mpnst <- allmetadata %>%
  filter(tumorType == "Malignant Peripheral Nerve Sheath Tumor",
         assay == "Bulk RNA sequencing",
         tissue == "primary tumor")

# Filter the metadata for "Plexiform Neurofibroma"
filtered_metadata_pn <- allmetadata %>%
  filter(tumorType == "Plexiform Neurofibroma",
         assay == "Bulk RNA sequencing",
         tissue == "primary tumor")

# Create a new column by merging specimenID and aliquotID for both filtered metadata
filtered_metadata_mpnst <- filtered_metadata_mpnst %>%
  mutate(merged_id = paste0(gsub("-", ".", specimenID), ".", aliquotID))

filtered_metadata_pn <- filtered_metadata_pn %>%
  mutate(merged_id = paste0(gsub("-", ".", specimenID), ".", aliquotID))

# Join the metadata tables together
metadata <- bind_rows(filtered_metadata_mpnst, filtered_metadata_pn)

# Ensure merged_id is the first column
metadata <- metadata %>% select(merged_id, everything())

# Ensure the sample names match between the data and metadata
colnames(expression_matrix) <- gsub("-", ".", colnames(expression_matrix))
metadata$merged_id <- gsub("-", ".", metadata$merged_id)

# Filter metadata to include only the samples present in the expression matrix
metadata <- metadata %>% filter(merged_id %in% colnames(expression_matrix))

stopifnot(rownames(expression_matrix) == gene_annotation$gene_id)
stopifnot(colnames(expression_matrix) == metadata$merged_id)

se <- SummarizedExperiment(assays = list(TPM = expression_matrix),
                           colData = metadata,
                           rowData = gene_annotation)
se

# Save the SummarizedExperiment object that contains metadata about PN and MPNST samples, expression in TPM and gene annotations
saveRDS(se, "summarized_experiment.rds")