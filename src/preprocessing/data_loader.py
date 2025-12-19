"""
RNA-seq Data Loader and Initial Processing
==========================================
GSE313799: Effect of chronic mechanical compression on human induced neurons

This module handles:
1. Loading raw counts matrix
2. Extracting sample metadata from column names
3. Initial data validation and QC
4. Database initialization and data storage
"""

import pandas as pd
import numpy as np
import sqlite3
import re
from pathlib import Path
from typing import Tuple, Dict, Optional
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class RNAseqDataLoader:
    """Load and parse RNA-seq counts data from GEO."""

    def __init__(self, counts_file: str, db_path: str):
        self.counts_file = Path(counts_file)
        self.db_path = Path(db_path)
        self.counts_df: Optional[pd.DataFrame] = None
        self.metadata_df: Optional[pd.DataFrame] = None

    def load_counts(self) -> pd.DataFrame:
        """Load counts matrix from file."""
        logger.info(f"Loading counts from {self.counts_file}")

        # Read the counts matrix
        self.counts_df = pd.read_csv(
            self.counts_file,
            sep='\t',
            index_col=0
        )

        logger.info(f"Loaded {self.counts_df.shape[0]} genes x {self.counts_df.shape[1]} samples")
        return self.counts_df

    def parse_sample_metadata(self) -> pd.DataFrame:
        """
        Extract metadata from sample names.

        Sample name format: {cell_line}_{condition}_{exp_date}_Seq_{seq_date}_{sample_num}
        Example: 1019_compress_04-25_Seq_06-26_S2
        """
        if self.counts_df is None:
            self.load_counts()

        metadata = []
        pattern = r'(.+?)_(compress|control)_(\d{2}-\d{2})_Seq_(\d{2}-\d{2})_S(\d+)'

        for sample_id in self.counts_df.columns:
            match = re.match(pattern, sample_id)
            if match:
                cell_line, condition, exp_date, seq_date, sample_num = match.groups()
                metadata.append({
                    'sample_id': sample_id,
                    'cell_line': cell_line,
                    'condition': condition,
                    'experiment_date': f"2024-{exp_date}",  # Assuming 2024
                    'sequencing_date': f"2024-{seq_date}",
                    'batch': seq_date,  # Use sequencing date as batch
                    'sample_number': int(sample_num)
                })
            else:
                logger.warning(f"Could not parse sample name: {sample_id}")

        self.metadata_df = pd.DataFrame(metadata)
        logger.info(f"Parsed metadata for {len(metadata)} samples")

        # Print summary
        print("\n=== Sample Metadata Summary ===")
        print(f"Cell lines: {self.metadata_df['cell_line'].unique().tolist()}")
        print(f"Conditions: {self.metadata_df['condition'].unique().tolist()}")
        print(f"Batches (seq dates): {self.metadata_df['batch'].unique().tolist()}")
        print("\nSamples per condition:")
        print(self.metadata_df.groupby(['cell_line', 'condition']).size())

        return self.metadata_df

    def calculate_qc_metrics(self) -> pd.DataFrame:
        """Calculate QC metrics for each sample."""
        if self.counts_df is None:
            self.load_counts()

        qc_metrics = []

        # Define mitochondrial and ribosomal gene patterns
        mito_pattern = r'^MT-'
        ribo_pattern = r'^RP[SL]\d+'

        mito_genes = self.counts_df.index[
            self.counts_df.index.str.match(mito_pattern, case=False, na=False)
        ]
        ribo_genes = self.counts_df.index[
            self.counts_df.index.str.match(ribo_pattern, case=False, na=False)
        ]

        for sample_id in self.counts_df.columns:
            counts = self.counts_df[sample_id]
            total_counts = counts.sum()

            qc_metrics.append({
                'sample_id': sample_id,
                'total_counts': int(total_counts),
                'detected_genes': int((counts > 0).sum()),
                'mitochondrial_pct': float(counts[mito_genes].sum() / total_counts * 100) if total_counts > 0 else 0,
                'ribosomal_pct': float(counts[ribo_genes].sum() / total_counts * 100) if total_counts > 0 else 0,
            })

        qc_df = pd.DataFrame(qc_metrics)

        print("\n=== QC Metrics Summary ===")
        print(qc_df.describe())

        return qc_df

    def filter_low_counts(
        self,
        min_counts: int = 10,
        min_samples: int = 3
    ) -> pd.DataFrame:
        """
        Filter genes with low counts.

        Keep genes that have at least `min_counts` in at least `min_samples` samples.
        """
        if self.counts_df is None:
            self.load_counts()

        n_genes_before = self.counts_df.shape[0]

        # Filter criterion: gene must have >= min_counts in >= min_samples samples
        keep = (self.counts_df >= min_counts).sum(axis=1) >= min_samples
        filtered_df = self.counts_df.loc[keep]

        n_genes_after = filtered_df.shape[0]
        logger.info(f"Filtered genes: {n_genes_before} -> {n_genes_after} "
                   f"(removed {n_genes_before - n_genes_after})")

        return filtered_df

    def init_database(self):
        """Initialize the SQLite database with schema."""
        schema_path = self.db_path.parent / 'schema.sql'

        with open(schema_path, 'r') as f:
            schema_sql = f.read()

        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.executescript(schema_sql)
        conn.commit()
        conn.close()

        logger.info(f"Database initialized at {self.db_path}")

    def save_to_database(self):
        """Save counts and metadata to database."""
        if self.counts_df is None:
            self.load_counts()
        if self.metadata_df is None:
            self.parse_sample_metadata()

        conn = sqlite3.connect(self.db_path)

        # Save metadata
        self.metadata_df.to_sql('samples', conn, if_exists='replace', index=False)
        logger.info("Saved sample metadata to database")

        # Save QC metrics
        qc_df = self.calculate_qc_metrics()
        qc_df['passed_qc'] = True  # Will be updated after QC filtering
        qc_df.to_sql('qc_metrics', conn, if_exists='replace', index=False)
        logger.info("Saved QC metrics to database")

        # Save raw counts (long format for database)
        logger.info("Saving raw counts to database (this may take a moment)...")
        counts_long = self.counts_df.reset_index().melt(
            id_vars=[self.counts_df.index.name or 'index'],
            var_name='sample_id',
            value_name='count'
        )
        counts_long.columns = ['gene_symbol', 'sample_id', 'count']

        # Use chunked insert for better performance
        counts_long.to_sql('raw_counts', conn, if_exists='replace', index=False)
        logger.info("Saved raw counts to database")

        conn.close()

    def get_design_matrix(self) -> pd.DataFrame:
        """Create design matrix for differential expression analysis."""
        if self.metadata_df is None:
            self.parse_sample_metadata()

        design = self.metadata_df[['sample_id', 'condition', 'cell_line', 'batch']].copy()

        # Create dummy variables
        design = pd.get_dummies(design, columns=['cell_line', 'batch'], drop_first=True)

        return design


def main():
    """Main function to load and process data."""
    # Paths
    project_root = Path(__file__).parent.parent.parent
    counts_file = project_root / "data/raw/GSE313799_counts_matrix_iN.txt"
    db_path = project_root / "db/rnaseq_analysis.db"

    # Initialize loader
    loader = RNAseqDataLoader(str(counts_file), str(db_path))

    # Load and process
    counts_df = loader.load_counts()
    metadata_df = loader.parse_sample_metadata()
    qc_df = loader.calculate_qc_metrics()

    # Initialize and save to database
    loader.init_database()
    loader.save_to_database()

    # Filter and save processed counts
    filtered_df = loader.filter_low_counts(min_counts=10, min_samples=3)
    filtered_df.to_csv(
        project_root / "data/processed/filtered_counts.csv"
    )

    print("\n=== Data Loading Complete ===")
    print(f"Raw counts: {counts_df.shape}")
    print(f"Filtered counts: {filtered_df.shape}")
    print(f"Database saved to: {db_path}")


if __name__ == "__main__":
    main()
