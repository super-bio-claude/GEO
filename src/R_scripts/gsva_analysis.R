#!/usr/bin/env Rscript
# =============================================================================
# GSVA: Gene Set Variation Analysis
# =============================================================================
#
# This script performs:
# 1. GSVA enrichment score calculation
# 2. ssGSEA (single-sample GSEA)
# 3. Gene set scoring with various methods
# 4. Differential pathway activity analysis
# 5. Integration with MSigDB gene sets
#
# Usage:
#   Rscript gsva_analysis.R <expression_file> <metadata_file> <output_dir>
#
# =============================================================================

suppressPackageStartupMessages({
    library(GSVA)
    library(GSEABase)
    library(msigdbr)
    library(limma)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(reshape2)
})

# =============================================================================
# Functions
# =============================================================================

get_msigdb_gene_sets <- function(species = "Homo sapiens",
                                  collection = "H",
                                  subcollection = NULL) {
    # Get gene sets from MSigDB (updated API for msigdbr >= 10.0.0)
    cat(sprintf("  Downloading MSigDB gene sets (collection: %s)...\n", collection))

    # Get gene sets using new API
    tryCatch({
        if (is.null(subcollection)) {
            msigdb_df <- msigdbr(species = species, collection = collection)
        } else {
            msigdb_df <- msigdbr(species = species, collection = collection,
                                subcollection = subcollection)
        }
    }, error = function(e) {
        # Fallback to old API for older versions
        if (is.null(subcollection)) {
            msigdb_df <<- msigdbr(species = species, category = collection)
        } else {
            msigdb_df <<- msigdbr(species = species, category = collection,
                                subcategory = subcollection)
        }
    })

    # Convert to list format for GSVA
    gene_sets <- split(msigdb_df$gene_symbol, msigdb_df$gs_name)

    cat(sprintf("    Retrieved %d gene sets\n", length(gene_sets)))

    return(gene_sets)
}


run_gsva <- function(expr_matrix, gene_sets, method = "gsva",
                     min_size = 10, max_size = 500, parallel_sz = 4) {
    # Run GSVA analysis
    cat(sprintf("\n  Running %s...\n", toupper(method)))

    # Filter gene sets by size
    gs_sizes <- sapply(gene_sets, length)
    gene_sets_filtered <- gene_sets[gs_sizes >= min_size & gs_sizes <= max_size]
    cat(sprintf("    Gene sets after size filter: %d\n", length(gene_sets_filtered)))

    # Run GSVA
    gsva_params <- gsvaParam(
        exprData = as.matrix(expr_matrix),
        geneSets = gene_sets_filtered,
        minSize = min_size,
        maxSize = max_size
    )

    if (method == "gsva") {
        gsva_es <- gsva(gsva_params, verbose = FALSE)
    } else if (method == "ssgsea") {
        ssgsea_params <- ssgseaParam(
            exprData = as.matrix(expr_matrix),
            geneSets = gene_sets_filtered,
            minSize = min_size,
            maxSize = max_size
        )
        gsva_es <- gsva(ssgsea_params, verbose = FALSE)
    } else if (method == "zscore") {
        zscore_params <- zscoreParam(
            exprData = as.matrix(expr_matrix),
            geneSets = gene_sets_filtered,
            minSize = min_size,
            maxSize = max_size
        )
        gsva_es <- gsva(zscore_params, verbose = FALSE)
    } else if (method == "plage") {
        plage_params <- plageParam(
            exprData = as.matrix(expr_matrix),
            geneSets = gene_sets_filtered,
            minSize = min_size,
            maxSize = max_size
        )
        gsva_es <- gsva(plage_params, verbose = FALSE)
    }

    cat(sprintf("    Computed scores for %d gene sets x %d samples\n",
                nrow(gsva_es), ncol(gsva_es)))

    return(gsva_es)
}


differential_pathway_analysis <- function(gsva_scores, metadata,
                                          condition_col = "condition",
                                          contrast = c("compress", "control")) {
    # Perform differential pathway activity analysis using limma
    cat("\n  Running differential pathway analysis...\n")

    # Create design matrix
    condition <- factor(metadata[[condition_col]], levels = rev(contrast))
    design <- model.matrix(~ condition)
    colnames(design) <- c("Intercept", paste0(contrast[1], "_vs_", contrast[2]))

    # Fit linear model
    fit <- lmFit(gsva_scores, design)
    fit <- eBayes(fit)

    # Get results
    results <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
    results$pathway <- rownames(results)
    results <- results[order(results$adj.P.Val), ]

    cat(sprintf("    Significant pathways (adj.P < 0.05): %d\n",
                sum(results$adj.P.Val < 0.05, na.rm = TRUE)))

    return(results)
}


# =============================================================================
# Main Analysis Function
# =============================================================================

run_gsva_analysis <- function(expr_file, metadata_file, output_dir) {

    cat("=== GSVA: Gene Set Variation Analysis ===\n")

    # Create output directory
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # ---------------------------------------------------------------------------
    # 1. Load Data
    # ---------------------------------------------------------------------------
    cat("\n[1/5] Loading data...\n")

    # Load expression matrix (VST or log-normalized)
    expr_matrix <- read.csv(expr_file, row.names = 1, check.names = FALSE)
    cat(sprintf("  Loaded expression: %d genes x %d samples\n",
                nrow(expr_matrix), ncol(expr_matrix)))

    # Load metadata
    metadata <- read.csv(metadata_file, check.names = FALSE)

    # Handle row names
    if ("sample_id" %in% colnames(metadata)) {
        rownames(metadata) <- metadata$sample_id
    } else if ("" %in% colnames(metadata) || is.null(colnames(metadata)[1]) || colnames(metadata)[1] == "") {
        rownames(metadata) <- metadata[, 1]
        metadata <- metadata[, -1]
    } else {
        rownames(metadata) <- metadata[, 1]
    }
    cat(sprintf("  Loaded metadata: %d samples\n", nrow(metadata)))

    # Match samples
    common_samples <- intersect(colnames(expr_matrix), rownames(metadata))
    cat(sprintf("  Common samples: %d\n", length(common_samples)))

    if (length(common_samples) == 0) {
        cat("  WARNING: No matching samples found!\n")
        cat("  Expression columns: ", paste(head(colnames(expr_matrix), 3), collapse = ", "), "...\n")
        cat("  Metadata rownames: ", paste(head(rownames(metadata), 3), collapse = ", "), "...\n")
        stop("No common samples between expression and metadata")
    }

    expr_matrix <- expr_matrix[, common_samples, drop = FALSE]
    metadata <- metadata[common_samples, , drop = FALSE]

    # ---------------------------------------------------------------------------
    # 2. Get Gene Sets from MSigDB
    # ---------------------------------------------------------------------------
    cat("\n[2/5] Retrieving gene sets from MSigDB...\n")

    gene_sets_list <- list()

    # Hallmark gene sets (H)
    gene_sets_list[["hallmark"]] <- get_msigdb_gene_sets(collection = "H")

    # KEGG pathways (C2:CP:KEGG_MEDICUS or KEGG_LEGACY)
    tryCatch({
        gene_sets_list[["kegg"]] <- get_msigdb_gene_sets(collection = "C2",
                                                          subcollection = "KEGG_MEDICUS")
    }, error = function(e) {
        cat("    KEGG_MEDICUS not found, trying KEGG...\n")
        gene_sets_list[["kegg"]] <<- get_msigdb_gene_sets(collection = "C2",
                                                           subcollection = "CP:KEGG")
    })

    # Reactome pathways (C2:CP:REACTOME)
    tryCatch({
        gene_sets_list[["reactome"]] <- get_msigdb_gene_sets(collection = "C2",
                                                              subcollection = "REACTOME")
    }, error = function(e) {
        cat("    REACTOME not found, trying CP:REACTOME...\n")
        gene_sets_list[["reactome"]] <<- get_msigdb_gene_sets(collection = "C2",
                                                               subcollection = "CP:REACTOME")
    })

    # GO Biological Process (C5:GO:BP)
    tryCatch({
        gene_sets_list[["go_bp"]] <- get_msigdb_gene_sets(collection = "C5",
                                                           subcollection = "BP")
    }, error = function(e) {
        cat("    BP not found, trying GO:BP...\n")
        gene_sets_list[["go_bp"]] <<- get_msigdb_gene_sets(collection = "C5",
                                                            subcollection = "GO:BP")
    })

    # ---------------------------------------------------------------------------
    # 3. Run GSVA for Each Gene Set Collection
    # ---------------------------------------------------------------------------
    cat("\n[3/5] Computing GSVA enrichment scores...\n")

    gsva_results <- list()

    for (gs_name in names(gene_sets_list)) {
        cat(sprintf("\n  Processing %s gene sets...\n", gs_name))

        # Run GSVA
        gsva_es <- run_gsva(
            expr_matrix = expr_matrix,
            gene_sets = gene_sets_list[[gs_name]],
            method = "gsva"
        )

        gsva_results[[gs_name]] <- gsva_es

        # Save scores
        write.csv(gsva_es,
                  file.path(output_dir, sprintf("gsva_scores_%s.csv", gs_name)))
    }

    # ---------------------------------------------------------------------------
    # 4. Differential Pathway Analysis
    # ---------------------------------------------------------------------------
    cat("\n[4/5] Running differential pathway analysis...\n")

    diff_results <- list()

    for (gs_name in names(gsva_results)) {
        cat(sprintf("\n  Analyzing %s...\n", gs_name))

        diff_res <- differential_pathway_analysis(
            gsva_scores = gsva_results[[gs_name]],
            metadata = metadata,
            condition_col = "condition",
            contrast = c("compress", "control")
        )

        diff_results[[gs_name]] <- diff_res

        # Save results
        write.csv(diff_res,
                  file.path(output_dir, sprintf("differential_%s.csv", gs_name)),
                  row.names = FALSE)
    }

    # ---------------------------------------------------------------------------
    # 5. Generate Visualizations
    # ---------------------------------------------------------------------------
    cat("\n[5/5] Generating visualizations...\n")

    for (gs_name in names(gsva_results)) {
        gsva_es <- gsva_results[[gs_name]]
        diff_res <- diff_results[[gs_name]]

        # Heatmap of top differential pathways
        sig_pathways <- diff_res$pathway[diff_res$adj.P.Val < 0.1]
        if (length(sig_pathways) >= 2) {
            # Need at least 2 pathways for clustering
            top_pathways <- head(sig_pathways, 30)
            mat <- gsva_es[top_pathways, , drop = FALSE]

            # Scale for better visualization
            mat_scaled <- t(scale(t(mat)))

            annotation_col <- data.frame(
                Condition = metadata$condition,
                row.names = rownames(metadata)
            )

            pdf(file.path(output_dir, sprintf("heatmap_%s.pdf", gs_name)),
                width = 12, height = max(8, length(top_pathways) * 0.3))
            pheatmap(mat_scaled,
                     annotation_col = annotation_col,
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     cluster_rows = length(top_pathways) > 1,
                     cluster_cols = TRUE,
                     fontsize_row = 8,
                     main = sprintf("Top Differential %s Pathways", toupper(gs_name)))
            dev.off()
        } else if (length(sig_pathways) == 1) {
            cat(sprintf("    Only 1 significant pathway for %s - skipping heatmap\n", gs_name))
        }

        # Boxplots of top pathways
        top_5 <- head(diff_res$pathway[order(diff_res$adj.P.Val)], 5)
        if (length(top_5) > 0) {
            pdf(file.path(output_dir, sprintf("boxplots_%s.pdf", gs_name)),
                width = 14, height = 10)
            par(mfrow = c(2, 3))

            for (pathway in top_5) {
                scores <- gsva_es[pathway, ]
                df <- data.frame(
                    score = scores,
                    condition = metadata$condition
                )

                boxplot(score ~ condition, data = df,
                       main = substr(pathway, 1, 50),
                       ylab = "GSVA Score",
                       col = c("steelblue", "coral"))
            }
            dev.off()
        }

        # Volcano plot
        pdf(file.path(output_dir, sprintf("volcano_%s.pdf", gs_name)),
            width = 10, height = 8)
        diff_res$significance <- "NS"
        diff_res$significance[diff_res$adj.P.Val < 0.05 & diff_res$logFC > 0] <- "Up"
        diff_res$significance[diff_res$adj.P.Val < 0.05 & diff_res$logFC < 0] <- "Down"

        p <- ggplot(diff_res, aes(x = logFC, y = -log10(P.Value),
                                   color = significance)) +
            geom_point(alpha = 0.6) +
            scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
            labs(title = sprintf("Differential %s Pathway Activity", toupper(gs_name)),
                 x = "Log2 Fold Change",
                 y = "-Log10(P-value)") +
            theme_minimal()
        print(p)
        dev.off()
    }

    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------
    cat("\n=== Analysis Complete ===\n")

    # Create summary
    summary_df <- data.frame(
        Gene_Set = names(diff_results),
        Total_Pathways = sapply(diff_results, nrow),
        Significant_005 = sapply(diff_results, function(x) sum(x$adj.P.Val < 0.05, na.rm = TRUE)),
        Significant_010 = sapply(diff_results, function(x) sum(x$adj.P.Val < 0.10, na.rm = TRUE))
    )

    write.csv(summary_df, file.path(output_dir, "gsva_summary.csv"), row.names = FALSE)
    print(summary_df)

    cat(sprintf("\nResults saved to: %s\n", output_dir))

    return(list(
        gsva_scores = gsva_results,
        differential = diff_results,
        summary = summary_df
    ))
}


# =============================================================================
# Additional: ssGSEA Analysis
# =============================================================================

run_ssgsea_analysis <- function(expr_file, metadata_file, output_dir) {
    # Run ssGSEA analysis as an alternative to GSVA
    cat("=== ssGSEA: Single-Sample GSEA ===\n")

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # Load data
    expr_matrix <- read.csv(expr_file, row.names = 1)
    metadata <- read.csv(metadata_file, row.names = 1)

    common_samples <- intersect(colnames(expr_matrix), rownames(metadata))
    expr_matrix <- expr_matrix[, common_samples]
    metadata <- metadata[common_samples, ]

    # Get Hallmark gene sets
    gene_sets <- get_msigdb_gene_sets(category = "H")

    # Run ssGSEA
    cat("\n  Computing ssGSEA scores...\n")

    ssgsea_params <- ssgseaParam(
        exprData = as.matrix(expr_matrix),
        geneSets = gene_sets,
        minSize = 10,
        maxSize = 500,
        normalize = TRUE
    )

    ssgsea_es <- gsva(ssgsea_params, verbose = FALSE)

    # Save
    write.csv(ssgsea_es, file.path(output_dir, "ssgsea_scores_hallmark.csv"))

    # Differential analysis
    diff_res <- differential_pathway_analysis(ssgsea_es, metadata)
    write.csv(diff_res, file.path(output_dir, "ssgsea_differential.csv"),
              row.names = FALSE)

    cat("\nssGSEA analysis complete.\n")

    return(ssgsea_es)
}


# =============================================================================
# Main Execution
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 3) {
    expr_file <- args[1]
    metadata_file <- args[2]
    output_dir <- args[3]

    results <- run_gsva_analysis(
        expr_file = expr_file,
        metadata_file = metadata_file,
        output_dir = output_dir
    )
} else {
    cat("Usage: Rscript gsva_analysis.R <expression_file> <metadata_file> <output_dir>\n")
    cat("\nExample:\n")
    cat("  Rscript gsva_analysis.R results/deseq2/vst_counts.csv results/sample_metadata.csv results/gsva\n")
}
