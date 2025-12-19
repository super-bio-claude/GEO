#!/usr/bin/env Rscript
# Prepare Alzheimer's disease data from GSE125583

library(GEOquery)
library(data.table)
library(org.Hs.eg.db)

cat("Loading GSE125583 metadata...\n")

# Get metadata
gse <- getGEO("GSE125583", destdir = "data/geo_test/GSE125583", GSEMatrix = TRUE)
pdata <- pData(gse[[1]])

# Extract sample info
sample_info <- data.frame(
  geo_accession = pdata$geo_accession,
  diagnosis = pdata$`diagnosis:ch1`,
  age = pdata$`age:ch1`,
  sex = pdata$`Sex:ch1`,
  stringsAsFactors = FALSE
)

# Create simple condition: AD vs control
sample_info$condition <- ifelse(sample_info$diagnosis == "control", "control", "AD")

cat("Sample distribution:\n")
print(table(sample_info$condition))

# Read and combine counts
data_dir <- "data/geo_test/GSE125583"
sample_files <- list.files(data_dir, pattern = "GSM.*tsv.gz$", full.names = TRUE)

cat(sprintf("\nReading %d sample files...\n", length(sample_files)))

counts_list <- list()
for (i in seq_along(sample_files)) {
  sample_id <- gsub("_sample.*", "", basename(sample_files[i]))
  df <- fread(cmd = sprintf("gunzip -c '%s'", sample_files[i]), skip = 2)
  colnames(df) <- c("gene_id", "nrpkm", "count")
  counts_list[[sample_id]] <- df[, c("gene_id", "count")]

  if (i %% 50 == 0) cat(sprintf("  Read %d/%d samples\n", i, length(sample_files)))
}

cat("Merging counts...\n")

# Merge all counts
counts_df <- counts_list[[1]]
colnames(counts_df)[2] <- names(counts_list)[1]

for (i in 2:length(counts_list)) {
  df <- counts_list[[i]]
  colnames(df)[2] <- names(counts_list)[i]
  counts_df <- merge(counts_df, df, by = "gene_id", all = TRUE)
}

# Set rownames
counts_matrix <- as.matrix(counts_df[, -1])
rownames(counts_matrix) <- counts_df$gene_id
counts_matrix[is.na(counts_matrix)] <- 0

cat(sprintf("Counts matrix: %d genes x %d samples\n", nrow(counts_matrix), ncol(counts_matrix)))

# Map Entrez IDs to Gene Symbols
cat("Mapping Entrez IDs to gene symbols...\n")
entrez_ids <- as.character(rownames(counts_matrix))
symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL",
                  keytype = "ENTREZID", multiVals = "first")

# Remove unmapped genes
valid_idx <- !is.na(symbols)
counts_matrix <- counts_matrix[valid_idx, ]
gene_symbols <- symbols[valid_idx]

# Handle duplicate symbols by summing counts
cat("Handling duplicate gene symbols...\n")
unique_symbols <- unique(gene_symbols)
aggregated_counts <- matrix(0, nrow = length(unique_symbols), ncol = ncol(counts_matrix))
rownames(aggregated_counts) <- unique_symbols
colnames(aggregated_counts) <- colnames(counts_matrix)

for (i in seq_along(gene_symbols)) {
  sym <- gene_symbols[i]
  aggregated_counts[sym, ] <- aggregated_counts[sym, ] + counts_matrix[i, ]
}

counts_matrix <- aggregated_counts
cat(sprintf("After symbol mapping: %d genes x %d samples\n", nrow(counts_matrix), ncol(counts_matrix)))

# Match samples with metadata
common_samples <- intersect(colnames(counts_matrix), sample_info$geo_accession)
cat(sprintf("Common samples: %d\n", length(common_samples)))

counts_matrix <- counts_matrix[, common_samples]
sample_info <- sample_info[sample_info$geo_accession %in% common_samples, ]
rownames(sample_info) <- sample_info$geo_accession
sample_info <- sample_info[colnames(counts_matrix), ]

# Save for DESeq2
write.csv(counts_matrix, "data/geo_test/GSE125583/counts_matrix.csv")
write.csv(sample_info, "data/geo_test/GSE125583/sample_metadata.csv")

cat("\nData saved! Ready for DESeq2 analysis.\n")
