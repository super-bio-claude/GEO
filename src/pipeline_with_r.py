"""
RNA-seq Analysis Pipeline with R Integration
=============================================

Complete pipeline including:
1. Data loading and preprocessing (Python)
2. Normalization (Python)
3. Batch correction (Python/ComBat)
4. DESeq2 + apeglm (R)
5. KEGG/GO Pathway Analysis (R)
6. GSVA (R)
7. Feature engineering (Python)
8. ML model training (Python)
9. Validation (Python)

Usage:
    python pipeline_with_r.py --config configs/config.yaml
    python pipeline_with_r.py --step deseq2
    python pipeline_with_r.py --step pathway
    python pipeline_with_r.py --step gsva
"""

import argparse
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from datetime import datetime

# Python modules
from preprocessing.data_loader import RNAseqDataLoader
from preprocessing.normalization import RNAseqNormalizer
from batch_correction.combat import ComBat, BatchEffectAnalyzer
from de_analysis.differential_expression import DEAnalysis, filter_significant
from feature_engineering.feature_selection import FeatureSelector, prepare_ml_data
from ml_models.classifiers import RNAseqClassifier, nested_cross_validation
from validation.model_evaluation import ModelEvaluator, generate_evaluation_report

# R integration
from r_integration.r_runner import (
    DESeq2Runner,
    PathwayAnalysisRunner,
    GSVARunner,
    IntegratedAnalysis,
    check_r_packages
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class RNAseqPipelineWithR:
    """Complete RNA-seq analysis pipeline with R integration."""

    def __init__(self, config_path: str):
        """Initialize pipeline."""
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)

        self.project_root = Path(config_path).parent.parent
        self.results_dir = self.project_root / "results"
        self.results_dir.mkdir(exist_ok=True)

        # Data containers
        self.counts_raw = None
        self.counts_filtered = None
        self.metadata = None

        # R runners
        self.deseq2_runner = DESeq2Runner(str(self.project_root))
        self.pathway_runner = PathwayAnalysisRunner(str(self.project_root))
        self.gsva_runner = GSVARunner(str(self.project_root))

        # Results containers
        self.deseq2_results = None
        self.pathway_results = None
        self.gsva_results = None
        self.ml_results = None

        logger.info(f"Initialized pipeline: {self.config['project']['name']}")

    # =========================================================================
    # Python-based preprocessing steps
    # =========================================================================

    def step1_load_data(self):
        """Load raw data."""
        logger.info("=" * 60)
        logger.info("Step 1: Loading Data")
        logger.info("=" * 60)

        counts_file = self.project_root / self.config['data']['raw_counts']
        db_path = self.project_root / self.config['data']['database']

        loader = RNAseqDataLoader(str(counts_file), str(db_path))
        self.counts_raw = loader.load_counts()
        self.metadata = loader.parse_sample_metadata()

        # Save metadata for R scripts
        self.metadata.to_csv(self.results_dir / "sample_metadata.csv", index=False)

        # Initialize database
        loader.init_database()

        return self.counts_raw, self.metadata

    def step2_filter_data(self):
        """Filter low-count genes."""
        logger.info("=" * 60)
        logger.info("Step 2: Filtering Low-Count Genes")
        logger.info("=" * 60)

        if self.counts_raw is None:
            self.step1_load_data()

        min_counts = self.config['preprocessing']['min_counts_per_gene']
        min_samples = self.config['preprocessing']['min_samples_per_gene']

        keep = (self.counts_raw >= min_counts).sum(axis=1) >= min_samples
        self.counts_filtered = self.counts_raw.loc[keep]

        logger.info(f"Filtered: {self.counts_raw.shape[0]} -> {self.counts_filtered.shape[0]} genes")

        # Save filtered counts
        self.counts_filtered.to_csv(
            self.project_root / "data/processed/filtered_counts.csv"
        )

        return self.counts_filtered

    # =========================================================================
    # R-based analysis steps
    # =========================================================================

    def step3_deseq2_analysis(self):
        """Run DESeq2 with apeglm (R)."""
        logger.info("=" * 60)
        logger.info("Step 3: DESeq2 + apeglm Analysis (R)")
        logger.info("=" * 60)

        if self.counts_raw is None:
            self.step1_load_data()

        # Prepare metadata with sample_id as index
        metadata_for_r = self.metadata.set_index('sample_id')
        metadata_file = str(self.results_dir / "sample_metadata.csv")
        metadata_for_r.to_csv(metadata_file)

        counts_file = str(self.project_root / self.config['data']['raw_counts'])
        output_dir = str(self.results_dir / "deseq2")

        self.deseq2_results = self.deseq2_runner.run_deseq2(
            counts_file=counts_file,
            metadata_file=metadata_file,
            output_dir=output_dir
        )

        if self.deseq2_results:
            logger.info("DESeq2 analysis completed successfully")
            if 'significant' in self.deseq2_results:
                n_sig = len(self.deseq2_results['significant'])
                logger.info(f"Significant genes (padj < 0.05): {n_sig}")
        else:
            logger.warning("DESeq2 analysis failed. Check R installation.")

        return self.deseq2_results

    def step4_pathway_analysis(self):
        """Run KEGG/GO pathway analysis (R)."""
        logger.info("=" * 60)
        logger.info("Step 4: KEGG/GO Pathway Analysis (R)")
        logger.info("=" * 60)

        # Check if DESeq2 results exist
        de_results_file = self.results_dir / "deseq2" / "deseq2_results_apeglm.csv"

        if not de_results_file.exists():
            logger.warning("DESeq2 results not found. Running DESeq2 first...")
            self.step3_deseq2_analysis()

        if not de_results_file.exists():
            logger.error("Cannot run pathway analysis without DE results")
            return None

        output_dir = str(self.results_dir / "pathway")

        self.pathway_results = self.pathway_runner.run_pathway_analysis(
            de_results_file=str(de_results_file),
            output_dir=output_dir
        )

        if self.pathway_results:
            logger.info("Pathway analysis completed successfully")
            if 'summary' in self.pathway_results:
                print("\nPathway Analysis Summary:")
                print(self.pathway_results['summary'])
        else:
            logger.warning("Pathway analysis failed")

        return self.pathway_results

    def step5_gsva_analysis(self):
        """Run GSVA analysis (R)."""
        logger.info("=" * 60)
        logger.info("Step 5: GSVA Analysis (R)")
        logger.info("=" * 60)

        # Check if VST counts exist (from DESeq2)
        vst_file = self.results_dir / "deseq2" / "vst_counts.csv"

        if not vst_file.exists():
            logger.warning("VST counts not found. Running DESeq2 first...")
            self.step3_deseq2_analysis()

        if not vst_file.exists():
            logger.error("Cannot run GSVA without normalized counts")
            return None

        metadata_file = str(self.results_dir / "sample_metadata.csv")
        output_dir = str(self.results_dir / "gsva")

        self.gsva_results = self.gsva_runner.run_gsva(
            expression_file=str(vst_file),
            metadata_file=metadata_file,
            output_dir=output_dir
        )

        if self.gsva_results:
            logger.info("GSVA analysis completed successfully")
            if 'summary' in self.gsva_results:
                print("\nGSVA Summary:")
                print(self.gsva_results['summary'])
        else:
            logger.warning("GSVA analysis failed")

        return self.gsva_results

    # =========================================================================
    # ML Pipeline (using GSVA scores as features)
    # =========================================================================

    def step6_ml_with_gsva_features(self):
        """Train ML models using GSVA pathway scores as features."""
        logger.info("=" * 60)
        logger.info("Step 6: ML with GSVA Pathway Features")
        logger.info("=" * 60)

        # Try to load GSVA scores
        gsva_file = self.results_dir / "gsva" / "gsva_scores_hallmark.csv"

        if not gsva_file.exists():
            logger.warning("GSVA scores not found. Running GSVA first...")
            self.step5_gsva_analysis()

        if not gsva_file.exists():
            logger.error("Cannot run ML without GSVA features")
            return None

        # Load GSVA scores (pathways x samples -> samples x pathways)
        gsva_scores = pd.read_csv(gsva_file, index_col=0).T

        # Load metadata
        metadata = pd.read_csv(self.results_dir / "sample_metadata.csv", index_col=0)

        # Prepare data
        common_samples = list(set(gsva_scores.index) & set(metadata.index))
        X = gsva_scores.loc[common_samples]
        y = metadata.loc[common_samples, 'condition']

        logger.info(f"ML Features: {X.shape[1]} pathways, {X.shape[0]} samples")

        # Train models
        classifier = RNAseqClassifier(X, y, test_size=0.2)
        self.ml_results = classifier.train_all_models(
            tune_hyperparameters=True,
            cv_folds=5
        )

        print("\nML Model Comparison (GSVA features):")
        print(self.ml_results)

        # Save results
        self.ml_results.to_csv(self.results_dir / "ml_gsva_results.csv", index=False)

        # Get feature importance (pathway importance)
        if 'RandomForest' in classifier.models:
            importance = classifier.get_feature_importance('RandomForest')
            importance.to_csv(self.results_dir / "pathway_importance.csv", index=False)
            print("\nTop 10 Important Pathways:")
            print(importance.head(10))

        return self.ml_results

    def step7_ml_with_gene_features(self):
        """Train ML models using gene expression features."""
        logger.info("=" * 60)
        logger.info("Step 7: ML with Gene Expression Features")
        logger.info("=" * 60)

        # Load VST counts
        vst_file = self.results_dir / "deseq2" / "vst_counts.csv"

        if not vst_file.exists():
            logger.warning("VST counts not found. Running DESeq2 first...")
            self.step3_deseq2_analysis()

        if not vst_file.exists():
            logger.error("Cannot run ML without expression data")
            return None

        # Load data
        vst_counts = pd.read_csv(vst_file, index_col=0)
        metadata = pd.read_csv(self.results_dir / "sample_metadata.csv", index_col=0)

        # Prepare data (transpose: samples x genes)
        X = vst_counts.T
        common_samples = list(set(X.index) & set(metadata.index))
        X = X.loc[common_samples]
        y = metadata.loc[common_samples, 'condition']

        # Feature selection using DE genes
        sig_genes_file = self.results_dir / "deseq2" / "significant_genes.csv"
        if sig_genes_file.exists():
            sig_genes = pd.read_csv(sig_genes_file)
            # Filter to top genes by adjusted p-value
            top_genes = sig_genes.nsmallest(500, 'padj')['gene'].tolist()
            top_genes = [g for g in top_genes if g in X.columns]
            X_selected = X[top_genes]
            logger.info(f"Selected {len(top_genes)} DE genes for ML")
        else:
            # Variance-based selection
            selector = FeatureSelector(X, y)
            selector.variance_filter(top_n=500)
            X_selected, _ = selector.get_selected_data('variance')
            logger.info(f"Selected {X_selected.shape[1]} high-variance genes")

        # Train models
        classifier = RNAseqClassifier(X_selected, y, test_size=0.2)
        gene_ml_results = classifier.train_all_models(
            tune_hyperparameters=True,
            cv_folds=5
        )

        print("\nML Model Comparison (Gene features):")
        print(gene_ml_results)

        # Save results
        gene_ml_results.to_csv(self.results_dir / "ml_gene_results.csv", index=False)

        # Feature importance (gene importance)
        if 'RandomForest' in classifier.models:
            importance = classifier.get_feature_importance('RandomForest')
            importance.to_csv(self.results_dir / "gene_importance.csv", index=False)
            print("\nTop 10 Important Genes:")
            print(importance.head(10))

        return gene_ml_results

    # =========================================================================
    # Full Pipeline
    # =========================================================================

    def run_full_pipeline(self):
        """Run complete analysis pipeline."""
        logger.info("=" * 60)
        logger.info("Starting Full RNA-seq Analysis Pipeline with R")
        logger.info("=" * 60)

        start_time = datetime.now()

        # Python preprocessing
        self.step1_load_data()
        self.step2_filter_data()

        # R-based analyses
        self.step3_deseq2_analysis()
        self.step4_pathway_analysis()
        self.step5_gsva_analysis()

        # ML analyses
        self.step6_ml_with_gsva_features()
        self.step7_ml_with_gene_features()

        end_time = datetime.now()
        duration = end_time - start_time

        logger.info("=" * 60)
        logger.info(f"Pipeline completed in {duration}")
        logger.info(f"Results saved to: {self.results_dir}")
        logger.info("=" * 60)

        self._generate_final_report()

    def _generate_final_report(self):
        """Generate final analysis report."""
        report = {
            'project': self.config['project']['name'],
            'date': datetime.now().isoformat(),
            'data': {
                'raw_genes': self.counts_raw.shape[0] if self.counts_raw is not None else None,
                'filtered_genes': self.counts_filtered.shape[0] if self.counts_filtered is not None else None,
                'samples': self.counts_raw.shape[1] if self.counts_raw is not None else None
            }
        }

        # DESeq2 results
        if self.deseq2_results and 'significant' in self.deseq2_results:
            report['deseq2'] = {
                'significant_genes': len(self.deseq2_results['significant']),
                'upregulated': len(self.deseq2_results.get('upregulated', [])),
                'downregulated': len(self.deseq2_results.get('downregulated', []))
            }

        # Pathway results
        if self.pathway_results and 'summary' in self.pathway_results:
            report['pathway'] = self.pathway_results['summary'].to_dict('records')

        # GSVA results
        if self.gsva_results and 'summary' in self.gsva_results:
            report['gsva'] = self.gsva_results['summary'].to_dict('records')

        # ML results
        if self.ml_results is not None:
            report['ml_gsva'] = self.ml_results.to_dict('records')

        import json
        with open(self.results_dir / "final_report.json", 'w') as f:
            json.dump(report, f, indent=2, default=str)

        logger.info(f"Final report saved to {self.results_dir / 'final_report.json'}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='RNA-seq Analysis Pipeline with R Integration'
    )
    parser.add_argument(
        '--config',
        type=str,
        default='configs/config.yaml',
        help='Path to configuration file'
    )
    parser.add_argument(
        '--step',
        type=str,
        choices=['all', 'load', 'filter', 'deseq2', 'pathway', 'gsva',
                 'ml-gsva', 'ml-gene', 'check-r'],
        default='all',
        help='Pipeline step to run'
    )

    args = parser.parse_args()

    # Check R packages
    if args.step == 'check-r':
        print("Checking R packages...")
        check_r_packages()
        return

    # Resolve config path
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = Path(__file__).parent.parent / config_path

    # Initialize pipeline
    pipeline = RNAseqPipelineWithR(str(config_path))

    # Run requested step
    step_map = {
        'all': pipeline.run_full_pipeline,
        'load': pipeline.step1_load_data,
        'filter': pipeline.step2_filter_data,
        'deseq2': pipeline.step3_deseq2_analysis,
        'pathway': pipeline.step4_pathway_analysis,
        'gsva': pipeline.step5_gsva_analysis,
        'ml-gsva': pipeline.step6_ml_with_gsva_features,
        'ml-gene': pipeline.step7_ml_with_gene_features
    }

    step_map[args.step]()


if __name__ == "__main__":
    main()
