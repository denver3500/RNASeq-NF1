library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)

# Load the SummarizedExperiment object
summarized_experiment <- readRDS("summarized_experiment.rds")

# Define the list of gene sets
gene_sets <- list(
  urea_cycle = c("CPS1", "OTC", "ASS1", "ASL", "ARG1", "ARG2"),
  glutamine_to_proline = c("GLS", "ALDH18A1", "PYCR1"),
  Collagen_I = c("COL1A1", "COL1A2"),
  Collagen_II = c("COL2A1"),
  Collagen_III = c("COL3A1"),
  Collagen_IV = c("COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6"),
  Collagen_V = c("COL5A1", "COL5A2", "COL5A3"),
  Collagen_VI = c("COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6"),
  Collagen_VII = c("COL7A1"),
  Collagen_VIII = c("COL8A1", "COL8A2"),
  Collagen_IX = c("COL9A1", "COL9A2", "COL9A3"),
  Collagen_X = c("COL10A1"),
  Collagen_XI = c("COL11A1", "COL11A2", "COL2A1"),
  Collagen_XII = c("COL12A1"),
  Collagen_XIII = c("COL13A1"),
  Collagen_XIV = c("COL14A1"),
  Collagen_XV = c("COL15A1"),
  Collagen_XVI = c("COL16A1"),
  Collagen_XVII = c("COL17A1"),
  Collagen_XVIII = c("COL18A1"),
  Collagen_XIX = c("COL19A1"),
  Collagen_XX = c("COL20A1"),
  Collagen_XXI = c("COL21A1"),
  Collagen_XXII = c("COL22A1"),
  Collagen_XXIII = c("COL23A1"),
  Collagen_XXIV = c("COL24A1"),
  Collagen_XXV = c("COL25A1"),
  Collagen_XXVI = c("COL26A1"),
  Collagen_XXVII = c("COL27A1"),
  Collagen_XXVIII = c("COL28A1"),
  Epithelial = c("CDH1", "TJP1", "DSP"),
  EMT = c("SNAI1", "SNAI2", "TWIST1", "LEF1"),
  Mesenchymal = c("CDH2", "VIM", "CTNNB1", "ACTA2"),
  Laminin_alpha = c("LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5"),
  Laminin_beta = c("LAMB1", "LAMB2", "LAMB3", "LAMB4"),
  Laminin_gamma = c("LAMC1", "LAMC2", "LAMC3"),
  Intergrin_alpha = c("ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGA7", "ITGA8", "ITGA9", "ITGA10", "ITGA11", "ITGAD", "ITGAE", "ITGAL", "ITGAM", "ITGAV", "ITGA2B", "ITGAX"),
  Integrin_beta = c("ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8"),
  Piezo_channels = c("PIEZO1", "PIEZO2")
)

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
  
  # Ensure the genes are plotted in the specified order
  plot_data$gene_name <- factor(plot_data$gene_name, levels = genes_to_plot)
  
  # Create the plot with significance
  p <- ggplot(plot_data, aes(x = gene_name, y = expression, fill = tumorType)) +
    geom_boxplot() +
    stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test", label.y = max(plot_data$expression) + 1) +
    labs(title = paste("Gene expression -", set_name),
         x = "Gene",
         y = "TPM",
         fill = NULL) +  # Remove the legend title
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot
  ggsave(filename = paste0("gene_expression_", set_name, ".png"), plot = p, width = 10, height = 6)
}
