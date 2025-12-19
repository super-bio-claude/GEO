-- RNA-seq Analysis Database Schema
-- GSE313799: Effect of chronic mechanical compression on human induced neurons

-- Sample metadata table
CREATE TABLE IF NOT EXISTS samples (
    sample_id TEXT PRIMARY KEY,
    cell_line TEXT NOT NULL,
    condition TEXT NOT NULL,  -- 'compress' or 'control'
    experiment_date TEXT,
    sequencing_date TEXT,
    batch TEXT,
    total_reads INTEGER,
    mapped_reads INTEGER,
    mapping_rate REAL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Gene annotation table
CREATE TABLE IF NOT EXISTS genes (
    gene_id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol TEXT UNIQUE NOT NULL,
    gene_type TEXT,
    chromosome TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    description TEXT
);

-- Raw counts table
CREATE TABLE IF NOT EXISTS raw_counts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol TEXT NOT NULL,
    sample_id TEXT NOT NULL,
    count INTEGER NOT NULL,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    UNIQUE(gene_symbol, sample_id)
);

-- Normalized counts table
CREATE TABLE IF NOT EXISTS normalized_counts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol TEXT NOT NULL,
    sample_id TEXT NOT NULL,
    tpm REAL,
    cpm REAL,
    log2_cpm REAL,
    vst REAL,  -- variance stabilizing transformation
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    UNIQUE(gene_symbol, sample_id)
);

-- Batch corrected counts
CREATE TABLE IF NOT EXISTS batch_corrected_counts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol TEXT NOT NULL,
    sample_id TEXT NOT NULL,
    corrected_value REAL NOT NULL,
    method TEXT NOT NULL,  -- 'ComBat', 'ComBat-seq', 'limma'
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    UNIQUE(gene_symbol, sample_id, method)
);

-- Differential expression results
CREATE TABLE IF NOT EXISTS de_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol TEXT NOT NULL,
    log2_fold_change REAL,
    pvalue REAL,
    padj REAL,  -- adjusted p-value (FDR)
    base_mean REAL,
    stat REAL,
    comparison TEXT NOT NULL,  -- e.g., 'compress_vs_control'
    method TEXT NOT NULL,  -- 'DESeq2', 'edgeR', 'limma-voom'
    cell_line TEXT,  -- NULL for combined analysis
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(gene_symbol, comparison, method, cell_line)
);

-- Feature engineering results
CREATE TABLE IF NOT EXISTS feature_sets (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    feature_name TEXT NOT NULL,
    feature_type TEXT NOT NULL,  -- 'gene', 'pathway', 'module'
    selection_method TEXT,  -- 'variance', 'de', 'lasso', 'random_forest'
    importance_score REAL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ML model results
CREATE TABLE IF NOT EXISTS ml_models (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    model_name TEXT NOT NULL,
    model_type TEXT NOT NULL,  -- 'RandomForest', 'SVM', 'XGBoost', etc.
    parameters TEXT,  -- JSON string of hyperparameters
    train_accuracy REAL,
    test_accuracy REAL,
    cv_score_mean REAL,
    cv_score_std REAL,
    auc_roc REAL,
    f1_score REAL,
    feature_set_id INTEGER,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (feature_set_id) REFERENCES feature_sets(id)
);

-- QC metrics table
CREATE TABLE IF NOT EXISTS qc_metrics (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    total_counts INTEGER,
    detected_genes INTEGER,
    mitochondrial_pct REAL,
    ribosomal_pct REAL,
    passed_qc BOOLEAN DEFAULT TRUE,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

-- Analysis run history
CREATE TABLE IF NOT EXISTS analysis_runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    run_type TEXT NOT NULL,
    parameters TEXT,
    status TEXT DEFAULT 'started',
    started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP
);

-- Create indexes for performance
CREATE INDEX IF NOT EXISTS idx_raw_counts_gene ON raw_counts(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_raw_counts_sample ON raw_counts(sample_id);
CREATE INDEX IF NOT EXISTS idx_normalized_gene ON normalized_counts(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_de_results_gene ON de_results(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_de_results_padj ON de_results(padj);
