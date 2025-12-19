#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 Analysis with apeglm Shrinkage
# =============================================================================
#
# This script performs:
# 1. DESeq2 differential expression analysis
# 2. apeglm log fold change shrinkage
# 3. Multiple contrast comparisons
# 4. Export results for downstream analysis
#
# Usage:
#   Rscript deseq2_analysis.R <counts_file> <metadata_file> <output_dir>
#
# =============================================================================

# Suppress package startup messages
suppressPackageStartupMessages({
    library(DESeq2)
    library(apeglm)
    library(BiocParallel)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

# Set up parallel processing
register(MulticoreParam(4))

# =============================================================================
# Functions
# =============================================================================

run_deseq2_analysis <- function(counts_file, metadata_file, output_dir,
                                 design_formula = "~ condition",
                                 contrast_var = "condition",
                                 numerator = "treatment",
                                 denominator = "control") {

    cat("=== DESeq2 Analysis with apeglm ===\n")

    # Create output directory
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # ---------------------------------------------------------------------------
    # 1. Load Data
    # ---------------------------------------------------------------------------
    cat("\n[1/6] Loading data...\n")

    # Load counts matrix (check.names = FALSE to preserve original column names)
    counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    cat(sprintf("  Loaded counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))

    # Load metadata
    metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
    cat(sprintf("  Loaded metadata: %d samples\n", nrow(metadata)))

    # Ensure sample order matches
    common_samples <- intersect(colnames(counts), rownames(metadata))
    counts <- counts[, common_samples]
    metadata <- metadata[common_samples, ]

    # Convert condition to factor with control as reference
    metadata[[contrast_var]] <- factor(metadata[[contrast_var]],
                                        levels = c(denominator, numerator))

    cat(sprintf("  Matched samples: %d\n", length(common_samples)))
    cat(sprintf("  Condition levels: %s (ref) vs %s\n", denominator, numerator))

    # ---------------------------------------------------------------------------
    # 2. Create DESeq2 Dataset
    # ---------------------------------------------------------------------------
    cat("\n[2/6] Creating DESeq2 dataset...\n")

    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(
        countData = as.matrix(counts),
        colData = metadata,
        design = as.formula(design_formula)
    )

    # Pre-filtering: remove genes with very low counts
    keep <- rowSums(counts(dds) >= 10) >= 3
    dds <- dds[keep, ]
    cat(sprintf("  Genes after filtering: %d\n", nrow(dds)))

    # ---------------------------------------------------------------------------
    # 3. Run DESeq2
    # ---------------------------------------------------------------------------
    cat("\n[3/6] Running DESeq2...\n")

    dds <- DESeq(dds, parallel = TRUE)

    # Get normalized counts
    normalized_counts <- counts(dds, normalized = TRUE)
    vst_counts <- assay(vst(dds, blind = FALSE))
    rlog_counts <- assay(rlog(dds, blind = FALSE))

    # Save normalized counts
    write.csv(normalized_counts, file.path(output_dir, "normalized_counts.csv"))
    write.csv(vst_counts, file.path(output_dir, "vst_counts.csv"))
    write.csv(rlog_counts, file.path(output_dir, "rlog_counts.csv"))

    cat("  Saved normalized counts (normalized, VST, rlog)\n")

    # ---------------------------------------------------------------------------
    # 4. Extract Results with apeglm Shrinkage
    # ---------------------------------------------------------------------------
    cat("\n[4/6] Applying apeglm LFC shrinkage...\n")

    # Get coefficient name
    coef_name <- paste0(contrast_var, "_", numerator, "_vs_", denominator)
    cat(sprintf("  Coefficient: %s\n", coef_name))

    # Results without shrinkage
    res_noshrink <- results(dds,
                           contrast = c(contrast_var, numerator, denominator),
                           alpha = 0.05)

    # Results with apeglm shrinkage
    res_apeglm <- lfcShrink(dds,
                            coef = coef_name,
                            type = "apeglm",
                            parallel = TRUE)

    # Results with ashr shrinkage (alternative)
    # res_ashr <- lfcShrink(dds,
    #                       coef = coef_name,
    #                       type = "ashr")

    # Summary
    cat("\n  Results Summary (padj < 0.05):\n")
    cat(sprintf("    Without shrinkage: %d DE genes\n",
                sum(res_noshrink$padj < 0.05, na.rm = TRUE)))
    cat(sprintf("    With apeglm: %d DE genes\n",
                sum(res_apeglm$padj < 0.05, na.rm = TRUE)))

    # ---------------------------------------------------------------------------
    # 5. Save Results
    # ---------------------------------------------------------------------------
    cat("\n[5/6] Saving results...\n")

    # Convert to data frames
    res_noshrink_df <- as.data.frame(res_noshrink)
    res_noshrink_df$gene <- rownames(res_noshrink_df)
    res_noshrink_df <- res_noshrink_df[order(res_noshrink_df$padj), ]

    res_apeglm_df <- as.data.frame(res_apeglm)
    res_apeglm_df$gene <- rownames(res_apeglm_df)
    res_apeglm_df <- res_apeglm_df[order(res_apeglm_df$padj), ]

    # Save to CSV
    write.csv(res_noshrink_df,
              file.path(output_dir, "deseq2_results_noshrink.csv"),
              row.names = FALSE)
    write.csv(res_apeglm_df,
              file.path(output_dir, "deseq2_results_apeglm.csv"),
              row.names = FALSE)

    # Save significant genes
    sig_genes <- res_apeglm_df[!is.na(res_apeglm_df$padj) &
                               res_apeglm_df$padj < 0.05, ]
    write.csv(sig_genes,
              file.path(output_dir, "significant_genes.csv"),
              row.names = FALSE)

    # Save upregulated and downregulated separately
    up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
    down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]

    write.csv(up_genes,
              file.path(output_dir, "upregulated_genes.csv"),
              row.names = FALSE)
    write.csv(down_genes,
              file.path(output_dir, "downregulated_genes.csv"),
              row.names = FALSE)

    cat(sprintf("  Significant genes: %d (up: %d, down: %d)\n",
                nrow(sig_genes), nrow(up_genes), nrow(down_genes)))

    # ---------------------------------------------------------------------------
    # 6. Generate Plots
    # ---------------------------------------------------------------------------
    cat("\n[6/6] Generating plots...\n")

    # MA plot
    pdf(file.path(output_dir, "MA_plot.pdf"), width = 10, height = 8)
    plotMA(res_apeglm, main = "MA Plot (apeglm shrinkage)", ylim = c(-5, 5))
    dev.off()

    # Volcano plot
    pdf(file.path(output_dir, "volcano_plot.pdf"), width = 10, height = 8)

    volcano_data <- res_apeglm_df
    volcano_data$significance <- "Not Significant"
    volcano_data$significance[volcano_data$padj < 0.05 &
                              volcano_data$log2FoldChange > 1] <- "Up"
    volcano_data$significance[volcano_data$padj < 0.05 &
                              volcano_data$log2FoldChange < -1] <- "Down"

    p <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue),
                                   color = significance)) +
        geom_point(alpha = 0.5, size = 1) +
        scale_color_manual(values = c("Up" = "red", "Down" = "blue",
                                      "Not Significant" = "gray")) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
        labs(title = paste("Volcano Plot:", numerator, "vs", denominator),
             x = "Log2 Fold Change",
             y = "-Log10(p-value)") +
        theme_minimal()
    print(p)
    dev.off()

    # PCA plot
    pdf(file.path(output_dir, "PCA_plot.pdf"), width = 10, height = 8)
    vsd <- vst(dds, blind = FALSE)
    pca_data <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
    percentVar <- round(100 * attr(pca_data, "percentVar"))

    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
        geom_point(size = 4) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        labs(title = "PCA Plot") +
        theme_minimal()
    print(p)
    dev.off()

    # Heatmap of top DE genes
    pdf(file.path(output_dir, "heatmap_top50.pdf"), width = 12, height = 10)
    top_genes <- head(sig_genes$gene, 50)
    if (length(top_genes) > 0) {
        mat <- vst_counts[top_genes, ]
        mat <- t(scale(t(mat)))

        annotation_col <- data.frame(
            Condition = metadata$condition,
            row.names = rownames(metadata)
        )

        pheatmap(mat,
                 annotation_col = annotation_col,
                 show_rownames = TRUE,
                 show_colnames = TRUE,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 main = "Top 50 DE Genes Heatmap")
    }
    dev.off()

    # Shrinkage comparison plot
    pdf(file.path(output_dir, "shrinkage_comparison.pdf"), width = 12, height = 5)
    par(mfrow = c(1, 2))
    plotMA(res_noshrink, main = "Without Shrinkage", ylim = c(-5, 5))
    plotMA(res_apeglm, main = "With apeglm Shrinkage", ylim = c(-5, 5))
    dev.off()

    cat("  Saved all plots\n")

    # ---------------------------------------------------------------------------
    # Return results summary
    # ---------------------------------------------------------------------------
    cat("\n=== Analysis Complete ===\n")
    cat(sprintf("Results saved to: %s\n", output_dir))

    return(list(
        dds = dds,
        res_noshrink = res_noshrink,
        res_apeglm = res_apeglm,
        normalized_counts = normalized_counts,
        vst_counts = vst_counts
    ))
}


# =============================================================================
# Main Execution
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 3) {
    counts_file <- args[1]
    metadata_file <- args[2]
    output_dir <- args[3]

    results <- run_deseq2_analysis(
        counts_file = counts_file,
        metadata_file = metadata_file,
        output_dir = output_dir
    )
} else if (interactive()) {
    # Default paths for interactive use
    cat("Running in interactive mode with default paths...\n")

    project_root <- dirname(dirname(dirname(normalizePath("."))))
    counts_file <- file.path(project_root, "data/raw/GSE313799_counts_matrix_iN.txt")
    metadata_file <- file.path(project_root, "results/sample_metadata.csv")
    output_dir <- file.path(project_root, "results/deseq2")

    if (file.exists(counts_file) && file.exists(metadata_file)) {
        results <- run_deseq2_analysis(
            counts_file = counts_file,
            metadata_file = metadata_file,
            output_dir = output_dir
        )
    } else {
        cat("Required files not found. Please provide paths as arguments.\n")
        cat("Usage: Rscript deseq2_analysis.R <counts_file> <metadata_file> <output_dir>\n")
    }
} else {
    cat("Usage: Rscript deseq2_analysis.R <counts_file> <metadata_file> <output_dir>\n")
}
