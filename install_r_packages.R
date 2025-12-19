#!/usr/bin/env Rscript
# =============================================================================
# Install Required R Packages for RNA-seq Analysis
# =============================================================================
#
# This script installs all required Bioconductor and CRAN packages
#
# Usage:
#   Rscript install_r_packages.R
#
# =============================================================================

cat("=== Installing R Packages for RNA-seq Analysis ===\n\n")

# Install BiocManager if not available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Bioconductor packages
bioc_packages <- c(
    # DESeq2 and related
    "DESeq2",
    "apeglm",
    "ashr",

    # Pathway analysis
    "clusterProfiler",
    "enrichplot",
    "org.Hs.eg.db",
    "DOSE",
    "pathview",
    "ReactomePA",

    # GSVA
    "GSVA",
    "GSEABase",
    "msigdbr",

    # General utilities
    "BiocParallel",
    "limma",
    "SummarizedExperiment"
)

# CRAN packages
cran_packages <- c(
    "ggplot2",
    "pheatmap",
    "RColorBrewer",
    "reshape2",
    "dplyr",
    "tidyr",
    "ggrepel"
)

# Install CRAN packages
cat("\n--- Installing CRAN packages ---\n")
for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("Installing %s...\n", pkg))
        install.packages(pkg, repos = "https://cloud.r-project.org")
    } else {
        cat(sprintf("%s: already installed\n", pkg))
    }
}

# Install Bioconductor packages
cat("\n--- Installing Bioconductor packages ---\n")
for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("Installing %s...\n", pkg))
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
        cat(sprintf("%s: already installed\n", pkg))
    }
}

# Verify installation
cat("\n\n=== Verification ===\n")
all_packages <- c(bioc_packages, cran_packages)
installed <- sapply(all_packages, function(p) requireNamespace(p, quietly = TRUE))

cat("\nInstallation Status:\n")
for (i in seq_along(all_packages)) {
    status <- ifelse(installed[i], "OK", "FAILED")
    cat(sprintf("  %s: %s\n", all_packages[i], status))
}

n_installed <- sum(installed)
n_total <- length(all_packages)

cat(sprintf("\n\nSummary: %d/%d packages installed successfully\n", n_installed, n_total))

if (all(installed)) {
    cat("\nAll packages installed successfully!\n")
    cat("You can now run the RNA-seq analysis pipeline.\n")
} else {
    cat("\nSome packages failed to install. Please check the error messages above.\n")
    failed <- all_packages[!installed]
    cat(sprintf("Failed packages: %s\n", paste(failed, collapse = ", ")))
}
