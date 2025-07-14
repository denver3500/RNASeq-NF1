# Load necessary libraries
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)

# Load the SummarizedExperiment object
summarized_experiment <- readRDS("summarized_experiment.rds")

# Define the list of gene sets
gene_sets <- list(
  TRAP1 = c("TRAP1", "COL6A1", "GAPDH")
)

# Initialize an empty data frame to store the results
final_table <- data.frame()

# Loop through each gene set
for (set_name in names(gene_sets)) {
  genes_to_plot <- gene_sets[[set_name]]
  
  # Access the rowData to find the gene_ids corresponding to the genes to plot
  gene_annotation <- as.data.frame(rowData(summarized_experiment))  # Convert rowData to a data frame if necessary
  gene_ids_to_plot <- gene_annotation %>% filter(gene_name %in% genes_to_plot) %>% select(gene_id, gene_name)
  
  # Access the assay data (TPM) and subset for the genes to plot
  tpm_data <- assay(summarized_experiment, "TPM")
  subset_tpm_data <- tpm_data[rownames(tpm_data) %in% gene_ids_to_plot$gene_id, , drop = FALSE]
  
  # Convert subset_tpm_data to a long data frame
  long_tpm_data <- as.data.frame(t(subset_tpm_data))
  long_tpm_data$merged_id <- rownames(long_tpm_data)
  long_tpm_data <- pivot_longer(long_tpm_data, cols = -merged_id, names_to = "gene_id", values_to = "expression")
  
  # Merge the long TPM data with the gene annotation to get gene names
  long_tpm_data <- long_tpm_data %>%
    left_join(gene_ids_to_plot, by = "gene_id")
  
  # Extract the sample metadata
  sample_metadata <- as.data.frame(colData(summarized_experiment))
  
  # Merge the sample metadata with the expression data
  plot_data <- long_tpm_data %>%
    left_join(sample_metadata, by = "merged_id")
  
  # Filter the data to include only the relevant tumor types and modify tumorType values for custom legend labels
  plot_data <- plot_data %>%
    select(merged_id, expression, gene_name, tumorType) %>%
    mutate(tumorType = case_when(
      tumorType == "Plexiform Neurofibroma" ~ "PN",
      tumorType == "Malignant Peripheral Nerve Sheath Tumor" ~ "MPNST",
      TRUE ~ tumorType
    )) %>%
    filter(tumorType %in% c("PN", "MPNST")) %>%
    mutate(tumorType = factor(tumorType, levels = c("PN", "MPNST")))
  
  # Append the plot_data to the final_table
  final_table <- bind_rows(final_table, plot_data)
}

# Filter GAPDH data
gapdh_data <- final_table %>%
  filter(gene_name == "GAPDH") %>%
  select(merged_id, gapdh_expression = expression)

# Normalize each sample's expression data by the GAPDH expression of that sample
normalized_table <- final_table %>%
  left_join(gapdh_data, by = "merged_id") %>%
  mutate(normalized_expression = expression / gapdh_expression) %>%
  select(-gapdh_expression) %>%
  filter(gene_name != "GAPDH")

# Filter data for each tumor type
pn_data <- normalized_table %>% filter(tumorType == "PN")
mpnst_data <- normalized_table %>% filter(tumorType == "MPNST")

library(ggplot2)
library(dplyr)

# Filter data for COL6A1 and TRAP1 for PN dataset
col6a1_pn <- pn_data %>% filter(gene_name == "COL6A1") %>% select(merged_id, col6a1_expression = normalized_expression)
trap1_pn <- pn_data %>% filter(gene_name == "TRAP1") %>% select(merged_id, trap1_expression = normalized_expression)

# Merge the data for the two genes
merged_pn <- col6a1_pn %>%
  inner_join(trap1_pn, by = "merged_id") %>%
  mutate(tumorType = "PN")

# Perform Pearson correlation test for PN dataset
cor_test_pn <- cor.test(merged_pn$col6a1_expression, merged_pn$trap1_expression)
r_squared_pn <- cor_test_pn$estimate^2
p_value_pn <- cor_test_pn$p.value

# Filter data for COL6A1 and TRAP1 for MPNST dataset
col6a1_mpnst <- mpnst_data %>% filter(gene_name == "COL6A1") %>% select(merged_id, col6a1_expression = normalized_expression)
trap1_mpnst <- mpnst_data %>% filter(gene_name == "TRAP1") %>% select(merged_id, trap1_expression = normalized_expression)

# Merge the data for the two genes
merged_mpnst <- col6a1_mpnst %>%
  inner_join(trap1_mpnst, by = "merged_id") %>%
  mutate(tumorType = "MPNST")

# Perform Pearson correlation test for MPNST dataset
cor_test_mpnst <- cor.test(merged_mpnst$col6a1_expression, merged_mpnst$trap1_expression)
r_squared_mpnst <- cor_test_mpnst$estimate^2
p_value_mpnst <- cor_test_mpnst$p.value

# Combine the data for both tumor types
combined_data <- bind_rows(merged_pn, merged_mpnst)

# Create annotation data frames
annotation_pn <- data.frame(
  tumorType = "PN",
  label = paste("R² =", round(r_squared_pn, 2), "\nP =", format(p_value_pn, digits = 2)),
  x = Inf,
  y = Inf
)

annotation_mpnst <- data.frame(
  tumorType = "MPNST",
  label = paste("R² =", round(r_squared_mpnst, 2), "\nP =", format(p_value_mpnst, digits = 2)),
  x = Inf,
  y = Inf
)

# Create the combined correlation plot
correlation_plot <- ggplot(combined_data, aes(x = col6a1_expression, y = trap1_expression)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~tumorType, scales = "free_x") +
  labs(title = "Correlation Plot for COL6A1 and TRAP1",
       x = "COL6A1 Expression",
       y = "TRAP1 Expression") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        strip.placement = "outside",
        strip.background = element_blank()) +
  geom_text(data = annotation_pn, aes(x = x, y = y, label = label), hjust = 1.1, vjust = 2, size = 5, color = "blue") +
  geom_text(data = annotation_mpnst, aes(x = x, y = y, label = label), hjust = 1.1, vjust = 2, size = 5, color = "blue")

# Display the plot
print(correlation_plot)

# Save the plot to a file
ggsave("COL6A1_TRAP1_correlation.png", plot = correlation_plot, width = 10, height = 6)