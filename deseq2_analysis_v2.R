#!/usr/bin/env Rscript

# ==========================
# Load libraries
# ==========================
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(dplyr))

# ==========================
# User inputs
# ==========================
args <- commandArgs(trailingOnly = TRUE)
# args[1] = gene count CSV
# args[2] = transcript count CSV
# args[3] = sample metadata CSV
# args[4] = output folder

gene_file <- args[1]
transcript_file <- args[2]
sample_info_file <- args[3]
output_dir <- args[4]

logfc_cutoff <- as.numeric(args[5])
padj_cutoff <- as.numeric(args[6])

# Create output folders
gene_out <- file.path(output_dir, "gene_deseq2_results")
tx_out <- file.path(output_dir, "transcript_deseq2_results")
dir.create(gene_out, recursive = TRUE, showWarnings = FALSE)
dir.create(tx_out, recursive = TRUE, showWarnings = FALSE)

# ==========================
# Load sample metadata
# ==========================
sample_info <- read.csv(sample_info_file, row.names = 1)
# Ensure rownames match count matrix columns
print("Sample metadata:")
print(sample_info)

# ==========================
# Gene-level DESeq2
# ==========================
gene_counts <- read.csv(gene_file, row.names = 1)
# Ensure column order matches sample info
gene_counts <- gene_counts[, rownames(sample_info)]
# Quick check

cat("Count matrix dimensions:", dim(gene_counts), "\n")
cat("First few values (should be integers):\n")
print(gene_counts[1:3, 1:4])

dds_gene <- DESeqDataSetFromMatrix(countData = gene_counts,
                                   colData = sample_info,
                                   design = ~ condition)



#  Pre-filter low-count genes ────────────────────────────

# Keep only genes with at least 10 counts in at least 3 samples
# This removes noise and speeds up computation significantly
keep <- rowSums(counts(dds_gene) >= 10) >= 3
dds_gene  <- dds_gene[keep, ]
cat("Genes remaining after filtering:", nrow(dds_gene), "\n")

#Run DESeq2 ─────────────────────────────────────────────
# This runs: size factor estimation → dispersion estimation → GLM fitting
dds_gene <- DESeq(dds_gene)

# Quick sanity checks
cat("\nSize factors (should be close to 1.0):\n")


print(sizeFactors(dds_gene))

print(nrow(dds_gene))


# Plot dispersion estimates — should show shrinkage toward fitted line
pdf(file.path(gene_out, "Dispersion_estimates.pdf"), width = 7, height = 6)
plotDispEsts(dds_gene, main = "Dispersion estimates")
dev.off()

# ── 6. PCA quality check ─────────────────────────────────────
# Use VST-normalized data for PCA (blind=TRUE for QC purposes)
vst <- vst(dds_gene, blind = TRUE)

# PCA plot — samples should cluster by condition
pca_data <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

# condition levels from metadata
conditions <- unique(pca_data$condition)

# generate automatic colors
colors <- setNames(
  scales::hue_pal()(length(conditions)),
  conditions
)

p<-ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  scale_color_manual(values = colors) +
  theme_bw() +
  labs(
    title = paste(
      "PCA —",
      paste(conditions, collapse = " vs ")
    )
  )
  # create automatic filename
pca_filename <- paste0(
  "PCA_",
  paste(conditions, collapse = "_vs_"),
  ".pdf"
)

  # save plot
ggsave(
  file.path(gene_out,"PCA_plot.pdf"),
  p,
  width=7,
  height=6
)

# Sample-to-sample distance heatmap
sample_dists <- dist(t(assay(vst)))
dist_matrix  <- as.matrix(sample_dists)
colors       <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf(file.path(gene_out, "Sample_distance_heatmap.pdf"), width = 8, height = 7)

pheatmap(
  dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  col  = colors,
  main = "Sample-to-sample distances"
)
dev.off()
#extract results
res_gene <- results(dds_gene)
cat("\nDESeq2 results summary:\n")
summary(res_gene)

#LFC shrinkage (strongly recommended) ──────────────────
# apeglm shrinkage stabilizes LFC estimates for low-count genes
# This does NOT change which genes are significant — only the LFC values
coef_name <- resultsNames(dds_gene)
print(coef_name)




# Build full annotated results table
res_df <- as.data.frame(res_gene) %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>%
  mutate(
    direction = case_when(
      log2FoldChange >  logfc_cutoff & padj < padj_cutoff ~ "upregulated",
      log2FoldChange < -logfc_cutoff & padj < padj_cutoff ~ "downregulated",
      TRUE ~ "not_significant"
    )
  )
  print(mcols(res_gene)$description)
 cat("\nDEG summary:\n")
print(table(res_df$direction))
write.csv(
  res_df,
  file.path(gene_out, "DESeq2_all_genes_results.csv"),
  row.names = FALSE
)
 # Filter significant DEGs ──────────────────────────────
# Standard cutoffs: adjusted p-value < 0.05 AND |log2FC| > cutoff

degs <- res_df %>%
  filter(padj < padj_cutoff,
         abs(log2FoldChange) > logfc_cutoff) %>%
  arrange(desc(abs(log2FoldChange)))

# Summary
cat("\nTotal significant DEGs:", nrow(degs), "\n")
cat("Upregulated:   ", sum(degs$log2FoldChange > 0), "\n")
cat("Downregulated: ", sum(degs$log2FoldChange < 0), "\n")

# Save DEG list
write.csv(
  degs,
  file.path(gene_out, "DESeq2_significant_DEGs.csv"),
  row.names = FALSE
)
# Top upregulated and downregulated genes ─────────────────

top_up <- degs %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(20)

top_down <- degs %>%
  filter(log2FoldChange < 0) %>%
  arrange(log2FoldChange) %>%
  head(20)

cat("\nTop 10 upregulated genes:\n")
print(top_up[1:10, c("gene","log2FoldChange","padj")])

cat("\nTop 10 downregulated genes:\n")
print(top_down[1:10, c("gene","log2FoldChange","padj")])

genes_to_label <- bind_rows(
  top_up[1:15, ],
  top_down[1:15, ]
)


# ==========================
# Transcript-level DESeq2
# ==========================
transcript_counts <- read.csv(transcript_file, row.names = 1)
transcript_counts <- transcript_counts[, rownames(sample_info)]

dds_tx <- DESeqDataSetFromMatrix(countData = transcript_counts,
                                 colData = sample_info,
                                 design = ~ condition)
dds_tx <- DESeq(dds_tx)
res_tx <- results(dds_tx)

# Save transcript results
write.csv(as.data.frame(res_tx),
          file.path(tx_out, "transcript_deseq2_results.csv"))
    
    
cat("Transcript-level DESeq2 complete. Results saved in:", tx_out, "\n")
saveRDS(
  dds_gene,
  file.path(gene_out, "dds_gene_object.rds")
)

saveRDS(
 dds_tx,
  file.path(tx_out, "dds_tx_object.rds")
)

# ==========================
# Ending message
# ==========================
cat("\n====================================================\n")
cat("✅ Pipeline executed successfully!\n")
cat("Gene-level and transcript-level DESeq2 results are ready.\n")
cat("Gene results folder:", gene_out, "\n")
cat("Transcript results folder:", tx_out, "\n")
cat("====================================================\n\n")

