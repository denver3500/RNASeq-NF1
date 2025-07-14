library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(reshape2)

# Load data
se <- readRDS("summarized_experiment.rds")

# Define collagen genes in the ORDER we want them displayed
collagen_genes <- c(
  "COL1A1", "COL1A2", "COL2A1", "COL3A1", 
  "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
  "COL5A1", "COL5A2", "COL5A3", 
  "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6",
  "COL7A1", "COL8A1", "COL8A2", 
  "COL9A1", "COL9A2", "COL9A3", "COL10A1", 
  "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", 
  "COL16A1", "COL17A1", "COL18A1", "COL19A1", "COL20A1", "COL21A1", 
  "COL22A1", "COL23A1", "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1"
)

print("Desired gene order:")
print(collagen_genes)

# Set up DESeq2
counts_matrix <- round(as.matrix(assay(se, "counts")))
colData(se)$tumorType <- factor(colData(se)$tumorType, 
                                levels = c("Plexiform Neurofibroma", "Malignant Peripheral Nerve Sheath Tumor"))

# Create a new SE with corrected counts
se_corrected <- SummarizedExperiment(
  assays = list(counts = counts_matrix),
  colData = colData(se),
  rowData = rowData(se)
)

dds <- DESeqDataSet(se_corrected, design = ~ tumorType)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# Get results and normalized counts
results <- results(dds, contrast = c("tumorType", "Malignant Peripheral Nerve Sheath Tumor", "Plexiform Neurofibroma"))
normalized_counts <- counts(dds, normalized = TRUE)

# Find collagen genes in our data
gene_annotation <- as.data.frame(rowData(se))
collagen_found <- gene_annotation %>%
  filter(gene_name %in% collagen_genes)

print("Collagen genes found in data:")
print(collagen_found$gene_name)

# Get normalized counts for collagen genes
collagen_counts <- normalized_counts[collagen_found$gene_id, ]
rownames(collagen_counts) <- collagen_found$gene_name

# Calculate mean expression by tumor type
sample_metadata <- as.data.frame(colData(se))
pn_samples <- sample_metadata$merged_id[sample_metadata$tumorType == "Plexiform Neurofibroma"]
mpnst_samples <- sample_metadata$merged_id[sample_metadata$tumorType == "Malignant Peripheral Nerve Sheath Tumor"]

# Calculate means
pn_means <- rowMeans(collagen_counts[, pn_samples, drop = FALSE])
mpnst_means <- rowMeans(collagen_counts[, mpnst_samples, drop = FALSE])

# Create a simple data frame with the exact order we want
heatmap_data <- data.frame(
  gene = collagen_genes,
  PN = NA,
  MPNST = NA,
  stringsAsFactors = FALSE
)

# Fill in the values for genes we found
for (i in 1:nrow(heatmap_data)) {
  gene <- heatmap_data$gene[i]
  if (gene %in% names(pn_means)) {
    heatmap_data$PN[i] <- pn_means[gene]
    heatmap_data$MPNST[i] <- mpnst_means[gene]
  }
}

# Remove genes we don't have data for
heatmap_data <- heatmap_data[!is.na(heatmap_data$PN), ]

print("Final gene order for heatmap:")
print(heatmap_data$gene)

# Log transform
heatmap_data$PN <- log2(heatmap_data$PN + 1)
heatmap_data$MPNST <- log2(heatmap_data$MPNST + 1)

# Get significance information
sig_data <- data.frame(
  gene = character(),
  padj = numeric(),
  log2FoldChange = numeric(),
  stringsAsFactors = FALSE
)

for (gene in heatmap_data$gene) {
  gene_id <- collagen_found$gene_id[collagen_found$gene_name == gene]
  if (length(gene_id) > 0) {
    sig_data <- rbind(sig_data, data.frame(
      gene = gene,
      padj = results[gene_id, "padj"],
      log2FoldChange = results[gene_id, "log2FoldChange"],
      stringsAsFactors = FALSE
    ))
  }
}

# Create long format for ggplot
heatmap_long <- reshape2::melt(heatmap_data, id.vars = "gene", variable.name = "tumor_type", value.name = "expression")

# Add significance stars
heatmap_long <- heatmap_long %>%
  left_join(sig_data, by = "gene") %>%
  mutate(
    stars = case_when(
      is.na(padj) ~ "",
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    ),
    display_stars = case_when(
      stars == "" ~ "",
      log2FoldChange > 0 & tumor_type == "MPNST" ~ stars,
      log2FoldChange < 0 & tumor_type == "PN" ~ stars,
      TRUE ~ ""
    )
  )

# Set factor levels to maintain order (reverse for ggplot2)
heatmap_long$gene <- factor(heatmap_long$gene, levels = rev(heatmap_data$gene))

# Create heatmap
create_heatmap <- function(transparent_bg = TRUE) {
  p <- ggplot(heatmap_long, aes(x = tumor_type, y = gene, fill = expression)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = display_stars), color = "black", size = 4, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "red", 
                        name = "Log2\nExpression") +
    scale_x_discrete(labels = c("PN" = "Plexiform Neurofibroma", "MPNST" = "Malignant Peripheral Nerve Sheath Tumor")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(5, 5, 5, 5, "pt")
    ) +
    labs(title = "Collagen Gene Expression") +
    coord_fixed(ratio = 0.5)
  
  if (!transparent_bg) {
    p <- p + theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA)
    )
  }
  
  return(p)
}

# Save heatmaps
p_transparent <- create_heatmap(transparent_bg = TRUE)
ggsave("collagen_heatmap_transparent.png", p_transparent, width = 6, height = 12, dpi = 300, bg = "transparent")

p_white <- create_heatmap(transparent_bg = FALSE)
ggsave("collagen_heatmap_white.png", p_white, width = 6, height = 12, dpi = 300, bg = "white")

print("Heatmaps saved!")
print("- collagen_heatmap_transparent.png")
print("- collagen_heatmap_white.png")

# Save results
write.csv(sig_data, "collagen_deseq2_results.csv", row.names = FALSE)
significant_genes <- sig_data[sig_data$padj < 0.05 & !is.na(sig_data$padj), ]
write.csv(significant_genes, "significant_collagen_genes.csv", row.names = FALSE)

print(paste("Significant genes:", nrow(significant_genes)))
print(paste("Total genes:", nrow(sig_data)))
