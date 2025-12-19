#!/usr/bin/env python3
"""
RNA-seq Disease Analysis Pipeline
==================================

Unified pipeline for comprehensive disease analysis from RNA-seq data.

Usage:
    python pipeline.py --counts <counts_file> --output <output_dir>
    python pipeline.py --de_results <de_results_file> --output <output_dir>

The pipeline performs:
1. Data preprocessing and normalization
2. Differential expression analysis
3. Disease similarity search
4. Disease stage prediction
5. Comprehensive report generation
"""

import os
import sys
import json
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from disease_analysis.disease_database import DiseaseDatabaseBuilder
from disease_analysis.similarity_engine import DiseaseSimilarityEngine
from disease_analysis.stage_predictor import DiseaseStagePredictor
from disease_analysis.report_generator import DiseaseAnalysisReportGenerator


class RNAseqDiseaseAnalysisPipeline:
    """
    Unified RNA-seq Disease Analysis Pipeline

    Analyzes gene expression data to:
    1. Identify similar disease states
    2. Predict disease progression stage
    3. Generate comprehensive clinical reports
    """

    def __init__(self,
                 data_dir: str = "data",
                 output_dir: str = "results/disease_analysis"):

        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize components
        print("\n" + "="*60)
        print("Initializing RNA-seq Disease Analysis Pipeline")
        print("="*60)

        self.db_builder = DiseaseDatabaseBuilder(str(self.data_dir / "disease_db"))
        self.similarity_engine = DiseaseSimilarityEngine(
            db_path=str(self.data_dir / "chromadb_disease")
        )
        self.stage_predictor = DiseaseStagePredictor()
        self.report_generator = DiseaseAnalysisReportGenerator(str(self.output_dir))

        self._initialized = False

    def initialize(self, rebuild_db: bool = False):
        """Initialize or load the disease database"""

        print("\n[Step 1] Initializing disease database...")

        db_file = self.data_dir / "disease_db" / "disease_signatures_full.json"

        if rebuild_db or not db_file.exists():
            print("  Building disease database...")
            signatures = self.db_builder.build_database()
        else:
            print("  Loading existing database...")
            signatures = self.db_builder.load_database()

        # Build or load similarity index
        if rebuild_db or not self.similarity_engine.load_index():
            print("\n  Building similarity index...")
            self.similarity_engine.build_index(signatures)

        self._initialized = True
        print("\n  Initialization complete!")

    def analyze_de_results(self,
                           de_file: str,
                           sample_id: str = None,
                           padj_threshold: float = 0.05,
                           lfc_threshold: float = 0.5) -> Dict:
        """
        Analyze differential expression results

        Parameters:
        -----------
        de_file : str
            Path to DE results CSV (must have 'gene', 'log2FoldChange', 'padj' columns)
        sample_id : str
            Sample identifier for the report
        padj_threshold : float
            Adjusted p-value threshold for significance
        lfc_threshold : float
            Log2 fold change threshold

        Returns:
        --------
        dict : Analysis results
        """

        if not self._initialized:
            self.initialize()

        print("\n" + "="*60)
        print("Running Disease Analysis")
        print("="*60)

        # Load DE results
        print("\n[Step 2] Loading DE results...")
        de_df = pd.read_csv(de_file)
        print(f"  Loaded {len(de_df)} genes")

        # Generate sample ID if not provided
        if not sample_id:
            sample_id = Path(de_file).stem

        # Extract significant genes
        if 'padj' in de_df.columns:
            sig_mask = (de_df['padj'] < padj_threshold) & (de_df['log2FoldChange'].abs() > lfc_threshold)
            sig_df = de_df[sig_mask]
        else:
            sig_df = de_df.nlargest(100, de_df.columns[1])  # Use top 100 by first numeric column

        print(f"  Significant genes: {len(sig_df)}")

        # Get gene lists
        all_genes = de_df['gene'].tolist()
        sig_genes = sig_df['gene'].tolist() if len(sig_df) > 0 else all_genes[:100]

        # Get fold changes
        fold_changes = dict(zip(de_df['gene'], de_df['log2FoldChange']))

        # Up/down regulated
        up_genes = de_df[de_df['log2FoldChange'] > lfc_threshold]['gene'].tolist()
        down_genes = de_df[de_df['log2FoldChange'] < -lfc_threshold]['gene'].tolist()

        # If few significant genes, use top by fold change
        if len(sig_genes) < 20:
            print("  Using top genes by |log2FC| due to few significant genes...")
            de_df['abs_lfc'] = de_df['log2FoldChange'].abs()
            top_genes = de_df.nlargest(200, 'abs_lfc')['gene'].tolist()
            query_genes = list(set(sig_genes + top_genes))
        else:
            query_genes = sig_genes

        # DE summary
        de_summary = {
            'n_total': len(de_df),
            'n_significant': len(sig_genes),
            'n_upregulated': len([g for g in sig_genes if fold_changes.get(g, 0) > 0]),
            'n_downregulated': len([g for g in sig_genes if fold_changes.get(g, 0) < 0]),
            'top_genes': sig_genes[:20] if sig_genes else de_df.nlargest(20, 'log2FoldChange')['gene'].tolist()
        }

        print(f"  Up-regulated: {de_summary['n_upregulated']}")
        print(f"  Down-regulated: {de_summary['n_downregulated']}")

        # Run similarity search
        print("\n[Step 3] Searching for similar disease signatures...")
        similarity_results = self.similarity_engine.search(
            query_genes=query_genes,
            query_up_genes=up_genes[:100],
            query_down_genes=down_genes[:100],
            fold_changes=fold_changes,
            top_k=50
        )

        print(f"  Found {len(similarity_results)} matching signatures")

        # Show top matches
        print("\n  Top 5 matches:")
        for i, r in enumerate(similarity_results[:5], 1):
            print(f"    {i}. {r.signature_name[:50]} ({r.category}) - Score: {r.combined_score:.3f}")

        # Predict disease stage
        print("\n[Step 4] Predicting disease stage...")
        stage_assessments = self.stage_predictor.assess_multiple_categories(
            expressed_genes=query_genes,
            fold_changes=fold_changes,
            similarity_results=similarity_results
        )

        primary_cat, primary_assessment = self.stage_predictor.get_primary_assessment(stage_assessments)
        print(f"  Primary category: {primary_cat.upper()}")
        print(f"  Predicted stage: {primary_assessment.predicted_stage}")
        print(f"  Severity score: {primary_assessment.severity_score:.2f}")

        # Generate report
        print("\n[Step 5] Generating analysis report...")
        report_file = self.report_generator.generate_report(
            sample_id=sample_id,
            similarity_results=similarity_results,
            stage_assessments=stage_assessments,
            de_summary=de_summary
        )

        # Save detailed results
        self._save_detailed_results(
            sample_id=sample_id,
            similarity_results=similarity_results,
            stage_assessments=stage_assessments,
            de_summary=de_summary
        )

        print("\n" + "="*60)
        print("Analysis Complete!")
        print("="*60)

        results = {
            'sample_id': sample_id,
            'similarity_results': similarity_results,
            'stage_assessments': stage_assessments,
            'primary_category': primary_cat,
            'primary_assessment': primary_assessment,
            'de_summary': de_summary,
            'report_file': report_file
        }

        self._print_summary(results)

        return results

    def analyze_expression_matrix(self,
                                  counts_file: str,
                                  metadata_file: Optional[str] = None,
                                  condition_column: str = "condition",
                                  control_value: str = "control",
                                  treatment_value: str = "treatment") -> Dict:
        """
        Analyze raw expression matrix (runs DE analysis first)

        Parameters:
        -----------
        counts_file : str
            Path to counts matrix (genes x samples)
        metadata_file : str
            Path to sample metadata CSV
        condition_column : str
            Column name for condition in metadata
        control_value : str
            Value indicating control samples
        treatment_value : str
            Value indicating treatment samples

        Returns:
        --------
        dict : Analysis results
        """

        print("\n[Step 0] Running differential expression analysis...")

        # Run DE analysis using existing R scripts
        import subprocess

        de_output = self.output_dir / "deseq2_results.csv"

        # Create temporary R script for DE analysis
        r_script = f'''
        library(DESeq2)

        # Load data
        counts <- read.csv("{counts_file}", row.names=1, check.names=FALSE)
        metadata <- read.csv("{metadata_file}", row.names=1)

        # Match samples
        common <- intersect(colnames(counts), rownames(metadata))
        counts <- counts[, common]
        metadata <- metadata[common, , drop=FALSE]

        # Create DESeq2 object
        dds <- DESeqDataSetFromMatrix(
            countData = round(counts),
            colData = metadata,
            design = ~ {condition_column}
        )

        # Filter low counts
        dds <- dds[rowSums(counts(dds)) >= 10, ]

        # Run DESeq2
        dds <- DESeq(dds)

        # Get results
        res <- results(dds, contrast=c("{condition_column}", "{treatment_value}", "{control_value}"))
        res <- as.data.frame(res)
        res$gene <- rownames(res)

        # Save
        write.csv(res, "{de_output}", row.names=FALSE)
        '''

        # Run R script
        result = subprocess.run(
            ['Rscript', '-e', r_script],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"  Warning: R script failed: {result.stderr}")
            raise RuntimeError("DE analysis failed")

        # Continue with DE results
        return self.analyze_de_results(str(de_output))

    def _save_detailed_results(self,
                               sample_id: str,
                               similarity_results: List,
                               stage_assessments: Dict,
                               de_summary: Dict):
        """Save detailed results to CSV files"""

        # Similarity results
        sim_df = pd.DataFrame([{
            'rank': i+1,
            'signature_id': r.signature_id,
            'signature_name': r.signature_name,
            'category': r.category,
            'stage': r.stage,
            'severity_score': r.severity_score,
            'embedding_similarity': r.embedding_similarity,
            'gene_overlap': r.gene_overlap,
            'jaccard_index': r.jaccard_index,
            'direction_score': r.direction_score,
            'combined_score': r.combined_score,
            'overlapping_genes': ','.join(r.overlapping_genes)
        } for i, r in enumerate(similarity_results)])

        sim_df.to_csv(self.output_dir / f"similarity_results_{sample_id}.csv", index=False)

        # Stage assessments
        stage_data = []
        for cat, assessment in stage_assessments.items():
            stage_data.append({
                'category': cat,
                'predicted_stage': assessment.predicted_stage,
                'confidence': assessment.confidence,
                'severity_score': assessment.severity_score,
                'biomarkers_present': ','.join(assessment.biomarkers_present),
                'biomarkers_absent': ','.join(assessment.biomarkers_absent),
                'recommendations': ' | '.join(assessment.recommendations)
            })

        stage_df = pd.DataFrame(stage_data)
        stage_df.to_csv(self.output_dir / f"stage_assessment_{sample_id}.csv", index=False)

    def _print_summary(self, results: Dict):
        """Print analysis summary"""

        print(f"\n{'='*60}")
        print("ANALYSIS SUMMARY")
        print(f"{'='*60}")

        print(f"\nSample: {results['sample_id']}")
        print(f"\nDE Analysis:")
        print(f"  Total genes: {results['de_summary']['n_total']}")
        print(f"  Significant: {results['de_summary']['n_significant']}")

        print(f"\nDisease Similarity:")
        print(f"  Top match: {results['similarity_results'][0].signature_name}")
        print(f"  Score: {results['similarity_results'][0].combined_score:.3f}")

        print(f"\nDisease Stage Assessment:")
        print(f"  Primary category: {results['primary_category'].upper()}")
        print(f"  Predicted stage: {results['primary_assessment'].predicted_stage.upper()}")
        print(f"  Severity: {results['primary_assessment'].severity_score:.2f}")
        print(f"  Confidence: {results['primary_assessment'].confidence:.2f}")

        if results['primary_assessment'].recommendations:
            print(f"\nRecommendations:")
            for rec in results['primary_assessment'].recommendations[:3]:
                print(f"  - {rec}")

        print(f"\nOutput files:")
        print(f"  Report: {results['report_file']}")
        print(f"  Results: {self.output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='RNA-seq Disease Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze DE results
  python pipeline.py --de_results results/deseq2/deseq2_results.csv --output results/disease_analysis

  # Analyze with custom thresholds
  python pipeline.py --de_results de_results.csv --padj 0.01 --lfc 1.0

  # Rebuild database
  python pipeline.py --de_results de_results.csv --rebuild_db
        """
    )

    parser.add_argument('--de_results', type=str, required=True,
                       help='Path to DE results CSV file')
    parser.add_argument('--output', type=str, default='results/disease_analysis',
                       help='Output directory')
    parser.add_argument('--data_dir', type=str, default='data',
                       help='Data directory')
    parser.add_argument('--sample_id', type=str, default=None,
                       help='Sample identifier')
    parser.add_argument('--padj', type=float, default=0.05,
                       help='Adjusted p-value threshold')
    parser.add_argument('--lfc', type=float, default=0.5,
                       help='Log2 fold change threshold')
    parser.add_argument('--rebuild_db', action='store_true',
                       help='Force rebuild of disease database')

    args = parser.parse_args()

    # Run pipeline
    pipeline = RNAseqDiseaseAnalysisPipeline(
        data_dir=args.data_dir,
        output_dir=args.output
    )

    pipeline.initialize(rebuild_db=args.rebuild_db)

    results = pipeline.analyze_de_results(
        de_file=args.de_results,
        sample_id=args.sample_id,
        padj_threshold=args.padj,
        lfc_threshold=args.lfc
    )

    return results


if __name__ == "__main__":
    main()
