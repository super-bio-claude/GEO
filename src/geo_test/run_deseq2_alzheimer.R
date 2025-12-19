#!/usr/bin/env Rscript
# DESeq2 analysis for Alzheimer's disease data (GSE125583)

library(DESeq2)
library(apeglm)

cat("Loading data...\n")

# Load prepared data
counts <- read.csv("data/geo_test/GSE125583/counts_matrix.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/geo_test/GSE125583/sample_metadata.csv", row.names = 1)

cat(sprintf("Counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("Metadata: %d samples\n", nrow(metadata)))

# Ensure sample order matches
common <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common]
metadata <- metadata[common, ]

# Convert to integer
counts <- round(counts)

# Create DESeq2 dataset
cat("\nCreating DESeq2 dataset...\n")
metadata$condition <- factor(metadata$condition, levels = c("control", "AD"))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# Filter low counts
cat("Filtering low count genes...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("Genes after filtering: %d\n", nrow(dds)))

# Run DESeq2
cat("\nRunning DESeq2 (this may take a few minutes for 289 samples)...\n")
dds <- DESeq(dds)

# Get results with apeglm shrinkage
cat("\nApplying apeglm shrinkage...\n")
resultsNames(dds)
res_shrunk <- lfcShrink(dds, coef = "condition_AD_vs_control", type = "apeglm")

# Convert to data frame
res_df <- as.data.frame(res_shrunk)
res_df$gene <- rownames(res_df)

# Summary
cat("\nSummary of differential expression:\n")
print(summary(res_shrunk, alpha = 0.05))

# Count significant genes
sig_genes <- sum(res_df$padj < 0.05, na.rm = TRUE)
sig_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0.5, na.rm = TRUE)
sig_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -0.5, na.rm = TRUE)

cat(sprintf("\nSignificant genes (padj < 0.05): %d\n", sig_genes))
cat(sprintf("Upregulated (log2FC > 0.5): %d\n", sig_up))
cat(sprintf("Downregulated (log2FC < -0.5): %d\n", sig_down))

# Top genes
cat("\nTop 10 upregulated genes:\n")
top_up <- head(res_df[order(-res_df$log2FoldChange), ], 10)
print(top_up[, c("gene", "log2FoldChange", "padj")])

cat("\nTop 10 downregulated genes:\n")
top_down <- head(res_df[order(res_df$log2FoldChange), ], 10)
print(top_down[, c("gene", "log2FoldChange", "padj")])

# Save results
output_dir <- "results/geo_test/GSE125583"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file.path(output_dir, "deseq2_results.csv"), row.names = FALSE)

cat(sprintf("\nResults saved to %s/deseq2_results.csv\n", output_dir))
