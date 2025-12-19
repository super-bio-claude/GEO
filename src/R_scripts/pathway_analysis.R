#!/usr/bin/env Rscript
# =============================================================================
# Pathway Analysis: KEGG, GO, Reactome
# =============================================================================
#
# This script performs:
# 1. Gene ID conversion (Symbol -> Entrez)
# 2. KEGG pathway enrichment analysis
# 3. GO enrichment analysis (BP, MF, CC)
# 4. Reactome pathway analysis
# 5. Gene Set Enrichment Analysis (GSEA)
# 6. Visualization of results
#
# Usage:
#   Rscript pathway_analysis.R <de_results_file> <output_dir> [organism]
#
# =============================================================================

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(DOSE)
    library(enrichplot)
    library(ggplot2)
    library(pathview)
})

# Check for ReactomePA
if (requireNamespace("ReactomePA", quietly = TRUE)) {
    library(ReactomePA)
    REACTOME_AVAILABLE <- TRUE
} else {
    REACTOME_AVAILABLE <- FALSE
    cat("Note: ReactomePA not installed. Reactome analysis will be skipped.\n")
}

# =============================================================================
# Functions
# =============================================================================

convert_gene_ids <- function(gene_symbols, from_type = "SYMBOL", to_type = "ENTREZID") {
    # Convert gene symbols to Entrez IDs
    # Use bitr for conversion
    converted <- bitr(gene_symbols,
                     fromType = from_type,
                     toType = to_type,
                     OrgDb = org.Hs.eg.db)

    cat(sprintf("  Converted %d/%d genes to %s\n",
                nrow(converted), length(gene_symbols), to_type))

    return(converted)
}


run_kegg_enrichment <- function(gene_list, universe = NULL, organism = "hsa",
                                 pvalue_cutoff = 0.05, qvalue_cutoff = 0.1) {
    # Run KEGG pathway enrichment analysis
    cat("\n  Running KEGG enrichment...\n")

    kegg_result <- enrichKEGG(
        gene = gene_list,
        universe = universe,
        organism = organism,
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
        minGSSize = 10,
        maxGSSize = 500
    )

    if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
        cat(sprintf("    Found %d enriched pathways\n", nrow(kegg_result)))
    } else {
        cat("    No enriched pathways found\n")
    }

    return(kegg_result)
}


run_go_enrichment <- function(gene_list, universe = NULL,
                               ont = "ALL",  # BP, MF, CC, or ALL
                               pvalue_cutoff = 0.05, qvalue_cutoff = 0.1) {
    # Run GO enrichment analysis
    cat(sprintf("\n  Running GO enrichment (ont=%s)...\n", ont))

    go_result <- enrichGO(
        gene = gene_list,
        universe = universe,
        OrgDb = org.Hs.eg.db,
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
        readable = TRUE,
        minGSSize = 10,
        maxGSSize = 500
    )

    if (!is.null(go_result) && nrow(go_result) > 0) {
        cat(sprintf("    Found %d enriched GO terms\n", nrow(go_result)))
    } else {
        cat("    No enriched GO terms found\n")
    }

    return(go_result)
}


run_gsea <- function(gene_list_ranked, gene_sets = "C2",
                     pvalue_cutoff = 0.05) {
    # Run Gene Set Enrichment Analysis with KEGG
    cat("\n  Running GSEA...\n")

    # Sort gene list by fold change (descending)
    gene_list_ranked <- sort(gene_list_ranked, decreasing = TRUE)

    # GSEA with KEGG
    gsea_result <- gseKEGG(
        geneList = gene_list_ranked,
        organism = "hsa",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = pvalue_cutoff,
        verbose = FALSE
    )

    if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
        cat(sprintf("    Found %d enriched gene sets\n", nrow(gsea_result)))
    } else {
        cat("    No enriched gene sets found\n")
    }

    return(gsea_result)
}


run_reactome_enrichment <- function(gene_list, universe = NULL,
                                     organism = "human",
                                     pvalue_cutoff = 0.05) {
    # Run Reactome pathway enrichment
    if (!REACTOME_AVAILABLE) {
        cat("\n  Skipping Reactome analysis (package not installed)\n")
        return(NULL)
    }

    cat("\n  Running Reactome enrichment...\n")

    reactome_result <- enrichPathway(
        gene = gene_list,
        universe = universe,
        organism = organism,
        pvalueCutoff = pvalue_cutoff,
        readable = TRUE
    )

    if (!is.null(reactome_result) && nrow(reactome_result) > 0) {
        cat(sprintf("    Found %d enriched pathways\n", nrow(reactome_result)))
    } else {
        cat("    No enriched pathways found\n")
    }

    return(reactome_result)
}


save_enrichment_results <- function(result, output_file) {
    # Save enrichment results to CSV
    if (!is.null(result) && nrow(result) > 0) {
        write.csv(as.data.frame(result), output_file, row.names = FALSE)
        cat(sprintf("    Saved to: %s\n", output_file))
    }
}


plot_enrichment <- function(result, output_prefix, type = "kegg", top_n = 20) {
    # Generate enrichment plots
    if (is.null(result) || nrow(result) == 0) {
        return()
    }

    # Dotplot
    pdf(paste0(output_prefix, "_dotplot.pdf"), width = 12, height = 10)
    p <- dotplot(result, showCategory = top_n) +
         ggtitle(paste(toupper(type), "Enrichment Dotplot"))
    print(p)
    dev.off()

    # Barplot
    pdf(paste0(output_prefix, "_barplot.pdf"), width = 12, height = 10)
    p <- barplot(result, showCategory = top_n) +
         ggtitle(paste(toupper(type), "Enrichment Barplot"))
    print(p)
    dev.off()

    # Cnetplot (gene-concept network)
    if (nrow(result) >= 5) {
        tryCatch({
            pdf(paste0(output_prefix, "_cnetplot.pdf"), width = 14, height = 12)
            p <- cnetplot(result, categorySize = "pvalue", foldChange = NULL,
                         showCategory = min(5, nrow(result)))
            print(p)
            dev.off()
        }, error = function(e) {
            cat(sprintf("    Warning: Could not generate cnetplot: %s\n", e$message))
        })
    }

    # Enrichment map
    if (nrow(result) >= 3) {
        tryCatch({
            result_pairwise <- pairwise_termsim(result)
            pdf(paste0(output_prefix, "_emap.pdf"), width = 12, height = 10)
            p <- emapplot(result_pairwise, showCategory = min(30, nrow(result)))
            print(p)
            dev.off()
        }, error = function(e) {
            cat(sprintf("    Warning: Could not generate emap: %s\n", e$message))
        })
    }

    cat(sprintf("    Saved plots to: %s_*.pdf\n", output_prefix))
}


# =============================================================================
# Main Analysis Function
# =============================================================================

run_pathway_analysis <- function(de_results_file, output_dir, organism = "human") {

    cat("=== Pathway Analysis: KEGG, GO, Reactome ===\n")

    # Create output directory
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # ---------------------------------------------------------------------------
    # 1. Load DE Results
    # ---------------------------------------------------------------------------
    cat("\n[1/5] Loading DE results...\n")

    de_results <- read.csv(de_results_file)
    cat(sprintf("  Loaded %d genes\n", nrow(de_results)))

    # Get significant genes
    sig_genes <- de_results[!is.na(de_results$padj) & de_results$padj < 0.05, ]
    up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
    down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]

    cat(sprintf("  Significant genes: %d (up: %d, down: %d)\n",
                nrow(sig_genes), nrow(up_genes), nrow(down_genes)))

    # ---------------------------------------------------------------------------
    # 2. Gene ID Conversion
    # ---------------------------------------------------------------------------
    cat("\n[2/5] Converting gene IDs...\n")

    # All genes (for universe)
    all_converted <- convert_gene_ids(de_results$gene)

    # Significant genes
    sig_converted <- convert_gene_ids(sig_genes$gene)
    up_converted <- convert_gene_ids(up_genes$gene)
    down_converted <- convert_gene_ids(down_genes$gene)

    # Create ranked gene list for GSEA (by log2FC)
    gene_fc <- de_results$log2FoldChange
    names(gene_fc) <- de_results$gene

    fc_converted <- merge(
        data.frame(SYMBOL = names(gene_fc), log2FC = gene_fc),
        all_converted,
        by = "SYMBOL"
    )
    gene_list_ranked <- fc_converted$log2FC
    names(gene_list_ranked) <- fc_converted$ENTREZID
    gene_list_ranked <- sort(gene_list_ranked, decreasing = TRUE)

    # ---------------------------------------------------------------------------
    # 3. KEGG Analysis
    # ---------------------------------------------------------------------------
    cat("\n[3/5] Running KEGG pathway analysis...\n")

    # All significant genes
    kegg_all <- run_kegg_enrichment(
        gene_list = sig_converted$ENTREZID,
        universe = all_converted$ENTREZID
    )
    save_enrichment_results(kegg_all, file.path(output_dir, "kegg_all_significant.csv"))
    plot_enrichment(kegg_all, file.path(output_dir, "kegg_all"), "kegg")

    # Upregulated genes
    kegg_up <- run_kegg_enrichment(
        gene_list = up_converted$ENTREZID,
        universe = all_converted$ENTREZID
    )
    save_enrichment_results(kegg_up, file.path(output_dir, "kegg_upregulated.csv"))
    plot_enrichment(kegg_up, file.path(output_dir, "kegg_up"), "kegg")

    # Downregulated genes
    kegg_down <- run_kegg_enrichment(
        gene_list = down_converted$ENTREZID,
        universe = all_converted$ENTREZID
    )
    save_enrichment_results(kegg_down, file.path(output_dir, "kegg_downregulated.csv"))
    plot_enrichment(kegg_down, file.path(output_dir, "kegg_down"), "kegg")

    # GSEA KEGG
    gsea_kegg <- run_gsea(gene_list_ranked)
    save_enrichment_results(gsea_kegg, file.path(output_dir, "gsea_kegg.csv"))

    if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
        pdf(file.path(output_dir, "gsea_kegg_plot.pdf"), width = 12, height = 10)
        p <- gseaplot2(gsea_kegg, geneSetID = 1:min(3, nrow(gsea_kegg)),
                      pvalue_table = TRUE)
        print(p)
        dev.off()
    }

    # ---------------------------------------------------------------------------
    # 4. GO Analysis
    # ---------------------------------------------------------------------------
    cat("\n[4/5] Running GO enrichment analysis...\n")

    # Biological Process
    go_bp <- run_go_enrichment(
        gene_list = sig_converted$ENTREZID,
        universe = all_converted$ENTREZID,
        ont = "BP"
    )
    save_enrichment_results(go_bp, file.path(output_dir, "go_bp.csv"))
    plot_enrichment(go_bp, file.path(output_dir, "go_bp"), "GO-BP")

    # Molecular Function
    go_mf <- run_go_enrichment(
        gene_list = sig_converted$ENTREZID,
        universe = all_converted$ENTREZID,
        ont = "MF"
    )
    save_enrichment_results(go_mf, file.path(output_dir, "go_mf.csv"))
    plot_enrichment(go_mf, file.path(output_dir, "go_mf"), "GO-MF")

    # Cellular Component
    go_cc <- run_go_enrichment(
        gene_list = sig_converted$ENTREZID,
        universe = all_converted$ENTREZID,
        ont = "CC"
    )
    save_enrichment_results(go_cc, file.path(output_dir, "go_cc.csv"))
    plot_enrichment(go_cc, file.path(output_dir, "go_cc"), "GO-CC")

    # ---------------------------------------------------------------------------
    # 5. Reactome Analysis
    # ---------------------------------------------------------------------------
    cat("\n[5/5] Running Reactome pathway analysis...\n")

    reactome_all <- run_reactome_enrichment(
        gene_list = sig_converted$ENTREZID,
        universe = all_converted$ENTREZID
    )
    save_enrichment_results(reactome_all, file.path(output_dir, "reactome.csv"))
    plot_enrichment(reactome_all, file.path(output_dir, "reactome"), "Reactome")

    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------
    cat("\n=== Analysis Complete ===\n")
    cat(sprintf("Results saved to: %s\n", output_dir))

    # Create summary
    summary_df <- data.frame(
        Analysis = c("KEGG (all)", "KEGG (up)", "KEGG (down)", "GSEA-KEGG",
                    "GO-BP", "GO-MF", "GO-CC", "Reactome"),
        Enriched_Terms = c(
            ifelse(!is.null(kegg_all), nrow(kegg_all), 0),
            ifelse(!is.null(kegg_up), nrow(kegg_up), 0),
            ifelse(!is.null(kegg_down), nrow(kegg_down), 0),
            ifelse(!is.null(gsea_kegg), nrow(gsea_kegg), 0),
            ifelse(!is.null(go_bp), nrow(go_bp), 0),
            ifelse(!is.null(go_mf), nrow(go_mf), 0),
            ifelse(!is.null(go_cc), nrow(go_cc), 0),
            ifelse(!is.null(reactome_all), nrow(reactome_all), 0)
        )
    )

    write.csv(summary_df, file.path(output_dir, "pathway_analysis_summary.csv"),
              row.names = FALSE)
    print(summary_df)

    return(list(
        kegg_all = kegg_all,
        kegg_up = kegg_up,
        kegg_down = kegg_down,
        gsea_kegg = gsea_kegg,
        go_bp = go_bp,
        go_mf = go_mf,
        go_cc = go_cc,
        reactome = reactome_all
    ))
}


# =============================================================================
# Main Execution
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
    de_results_file <- args[1]
    output_dir <- args[2]
    organism <- ifelse(length(args) >= 3, args[3], "human")

    results <- run_pathway_analysis(
        de_results_file = de_results_file,
        output_dir = output_dir,
        organism = organism
    )
} else {
    cat("Usage: Rscript pathway_analysis.R <de_results_file> <output_dir> [organism]\n")
    cat("\nExample:\n")
    cat("  Rscript pathway_analysis.R results/deseq2/deseq2_results_apeglm.csv results/pathway\n")
}
