"""
RNA-seq Analysis Pipeline
=========================

Main pipeline script that orchestrates:
1. Data loading and preprocessing
2. Normalization
3. Batch correction
4. Differential expression analysis
5. Feature engineering
6. ML model training
7. Validation

Usage:
    python pipeline.py --config configs/config.yaml
"""

import argparse
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sqlite3
from datetime import datetime

# Import modules
from preprocessing.data_loader import RNAseqDataLoader
from preprocessing.normalization import RNAseqNormalizer
from batch_correction.combat import ComBat, BatchEffectAnalyzer
from de_analysis.differential_expression import DEAnalysis, filter_significant
from feature_engineering.feature_selection import FeatureSelector, prepare_ml_data
from ml_models.classifiers import RNAseqClassifier, nested_cross_validation
from validation.model_evaluation import ModelEvaluator, generate_evaluation_report

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class RNAseqPipeline:
    """Complete RNA-seq analysis pipeline."""

    def __init__(self, config_path: str):
        """
        Initialize pipeline with configuration.

        Parameters
        ----------
        config_path : str
            Path to YAML configuration file
        """
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)

        self.project_root = Path(config_path).parent.parent
        self.results_dir = self.project_root / "results"
        self.results_dir.mkdir(exist_ok=True)

        # Initialize containers
        self.counts_raw = None
        self.counts_filtered = None
        self.counts_normalized = None
        self.counts_batch_corrected = None
        self.metadata = None
        self.de_results = None
        self.ml_results = None

        logger.info(f"Initialized pipeline for: {self.config['project']['name']}")

    def step1_load_data(self):
        """Load and parse raw data."""
        logger.info("=== Step 1: Loading Data ===")

        counts_file = self.project_root / self.config['data']['raw_counts']
        db_path = self.project_root / self.config['data']['database']

        loader = RNAseqDataLoader(str(counts_file), str(db_path))
        self.counts_raw = loader.load_counts()
        self.metadata = loader.parse_sample_metadata()

        # Initialize database
        loader.init_database()
        loader.save_to_database()

        # Calculate QC metrics
        qc_df = loader.calculate_qc_metrics()

        logger.info(f"Loaded {self.counts_raw.shape[0]} genes x {self.counts_raw.shape[1]} samples")

        return self.counts_raw, self.metadata

    def step2_filter_and_normalize(self):
        """Filter low counts and normalize."""
        logger.info("=== Step 2: Filtering and Normalization ===")

        if self.counts_raw is None:
            self.step1_load_data()

        # Filter low count genes
        min_counts = self.config['preprocessing']['min_counts_per_gene']
        min_samples = self.config['preprocessing']['min_samples_per_gene']

        keep = (self.counts_raw >= min_counts).sum(axis=1) >= min_samples
        self.counts_filtered = self.counts_raw.loc[keep]

        logger.info(f"Filtered: {self.counts_raw.shape[0]} -> {self.counts_filtered.shape[0]} genes")

        # Normalize
        normalizer = RNAseqNormalizer(self.counts_filtered)
        normalized = normalizer.normalize_all()

        # Use log-CPM as default
        self.counts_normalized = normalized['cpm']

        # Save normalized data
        for method, df in normalized.items():
            output_path = self.project_root / f"data/normalized/{method}_normalized.csv"
            output_path.parent.mkdir(exist_ok=True)
            df.to_csv(output_path)

        logger.info("Normalization complete")

        return self.counts_normalized

    def step3_batch_correction(self):
        """Apply batch correction."""
        logger.info("=== Step 3: Batch Correction ===")

        if self.counts_normalized is None:
            self.step2_filter_and_normalize()

        # Get batch information
        batch = self.metadata.set_index('sample_id')['batch']

        # Apply ComBat
        combat = ComBat(
            data=self.counts_normalized,
            batch=batch,
            covariates=None,
            parametric=True
        )
        self.counts_batch_corrected = combat.fit_transform()

        # Analyze batch effects before and after
        metadata_indexed = self.metadata.set_index('sample_id')

        analyzer_before = BatchEffectAnalyzer(self.counts_normalized, metadata_indexed)
        pca_before, var_before = analyzer_before.pca_analysis()

        analyzer_after = BatchEffectAnalyzer(self.counts_batch_corrected, metadata_indexed)
        pca_after, var_after = analyzer_after.pca_analysis()

        logger.info(f"PCA variance explained before: {var_before[:3]}")
        logger.info(f"PCA variance explained after: {var_after[:3]}")

        # Save batch-corrected data
        self.counts_batch_corrected.to_csv(
            self.project_root / "data/processed/batch_corrected_counts.csv"
        )

        return self.counts_batch_corrected

    def step4_differential_expression(self):
        """Run differential expression analysis."""
        logger.info("=== Step 4: Differential Expression Analysis ===")

        if self.counts_filtered is None:
            self.step2_filter_and_normalize()

        metadata_indexed = self.metadata.set_index('sample_id')

        # Use raw filtered counts for DE analysis
        de = DEAnalysis(
            self.counts_filtered,
            metadata_indexed,
            condition_col='condition'
        )

        # Run analysis
        contrast = ('compress', 'control')
        self.de_results = de.run_all_methods(contrast)

        # Compare methods
        comparison = de.compare_methods(self.de_results)
        logger.info(f"\nDE Method Comparison:\n{comparison}")

        # Save results
        for method, df in self.de_results.items():
            output_path = self.results_dir / f"de_{method}.csv"
            df.to_csv(output_path, index=False)

        # Filter significant genes
        padj_thresh = self.config['de_analysis']['thresholds']['padj']
        lfc_thresh = self.config['de_analysis']['thresholds']['log2fc']

        sig_genes = filter_significant(
            self.de_results['deseq2_like'],
            padj_threshold=padj_thresh,
            log2fc_threshold=lfc_thresh
        )

        logger.info(f"Significant genes (padj<{padj_thresh}, |log2FC|>{lfc_thresh}): {len(sig_genes)}")

        return self.de_results

    def step5_feature_engineering(self):
        """Feature selection for ML."""
        logger.info("=== Step 5: Feature Engineering ===")

        if self.counts_batch_corrected is None:
            self.step3_batch_correction()
        if self.de_results is None:
            self.step4_differential_expression()

        # Prepare data for ML
        metadata_indexed = self.metadata.set_index('sample_id')
        X, y = prepare_ml_data(self.counts_batch_corrected, metadata_indexed, 'condition')

        # Feature selection
        selector = FeatureSelector(X, y)

        # Run multiple selection methods
        selector.variance_filter(top_n=1000)
        selector.mutual_information(n_features=200)
        selector.random_forest_importance(n_features=200)

        # DE-based selection
        selector.de_filter(
            self.de_results['deseq2_like'],
            padj_threshold=0.1,
            log2fc_threshold=0.5
        )

        # Ensemble selection
        selector.ensemble_selection(
            methods=['variance', 'mutual_info', 'random_forest'],
            min_votes=2,
            n_features=200
        )

        # Summary
        logger.info(f"\nFeature Selection Summary:\n{selector.get_summary()}")

        # Save selected features
        for method, genes in selector.selected_features.items():
            pd.DataFrame({'gene': genes}).to_csv(
                self.results_dir / f"selected_features_{method}.csv",
                index=False
            )

        return selector

    def step6_train_models(self, feature_method: str = 'ensemble'):
        """Train ML models."""
        logger.info("=== Step 6: ML Model Training ===")

        if self.counts_batch_corrected is None:
            self.step3_batch_correction()

        # Prepare data
        metadata_indexed = self.metadata.set_index('sample_id')
        X, y = prepare_ml_data(self.counts_batch_corrected, metadata_indexed, 'condition')

        # Feature selection
        selector = FeatureSelector(X, y)
        selector.variance_filter(top_n=500)
        X_selected, y_selected = selector.get_selected_data('variance')

        logger.info(f"Training with {X_selected.shape[1]} features")

        # Train models
        classifier = RNAseqClassifier(X_selected, y_selected, test_size=0.2)
        self.ml_results = classifier.train_all_models(
            tune_hyperparameters=True,
            cv_folds=5
        )

        logger.info(f"\nModel Results:\n{self.ml_results}")

        # Save results
        self.ml_results.to_csv(self.results_dir / "ml_model_comparison.csv", index=False)

        # Feature importance from best model
        best_model = self.ml_results.iloc[0]['model']
        if best_model in ['RandomForest', 'GradientBoosting']:
            importance = classifier.get_feature_importance(best_model)
            importance.to_csv(self.results_dir / "feature_importance.csv", index=False)

        # Save best model
        classifier.save_model(
            best_model,
            str(self.results_dir / f"best_model_{best_model}.joblib")
        )

        return classifier

    def step7_validation(self):
        """Comprehensive model validation."""
        logger.info("=== Step 7: Model Validation ===")

        if self.counts_batch_corrected is None:
            self.step3_batch_correction()

        # Prepare data
        metadata_indexed = self.metadata.set_index('sample_id')
        X, y = prepare_ml_data(self.counts_batch_corrected, metadata_indexed, 'condition')

        # Feature selection
        selector = FeatureSelector(X, y)
        selector.variance_filter(top_n=500)
        X_selected, y_selected = selector.get_selected_data('variance')

        # Create model
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler

        model = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', RandomForestClassifier(n_estimators=100, random_state=42))
        ])

        # Evaluate
        evaluator = ModelEvaluator(model, X_selected, y_selected)

        # Cross-validation
        cv_results = evaluator.cross_validation('stratified', n_splits=5)
        logger.info(f"CV Accuracy: {cv_results['mean_score']:.4f} +/- {cv_results['std_score']:.4f}")

        # Permutation test
        perm_results = evaluator.permutation_test(n_permutations=100)
        logger.info(f"Permutation test p-value: {perm_results['pvalue']:.4f}")

        # Generate report
        report = generate_evaluation_report(
            evaluator,
            str(self.results_dir / "evaluation_report.json")
        )

        # Nested CV for unbiased estimate
        nested_results = nested_cross_validation(
            X_selected, y_selected,
            model_name='RandomForest',
            outer_cv=5,
            inner_cv=3
        )
        logger.info(f"Nested CV: {nested_results['mean_score']:.4f} +/- {nested_results['std_score']:.4f}")

        return report

    def run_full_pipeline(self):
        """Run the complete analysis pipeline."""
        logger.info("="*60)
        logger.info("Starting Full RNA-seq Analysis Pipeline")
        logger.info("="*60)

        start_time = datetime.now()

        # Run all steps
        self.step1_load_data()
        self.step2_filter_and_normalize()
        self.step3_batch_correction()
        self.step4_differential_expression()
        self.step5_feature_engineering()
        classifier = self.step6_train_models()
        self.step7_validation()

        end_time = datetime.now()
        duration = end_time - start_time

        logger.info("="*60)
        logger.info(f"Pipeline completed in {duration}")
        logger.info(f"Results saved to: {self.results_dir}")
        logger.info("="*60)

        # Generate summary report
        self._generate_summary_report()

        return self.ml_results

    def _generate_summary_report(self):
        """Generate a summary report of the analysis."""
        summary = {
            'project': self.config['project']['name'],
            'date': datetime.now().isoformat(),
            'data': {
                'raw_genes': self.counts_raw.shape[0] if self.counts_raw is not None else None,
                'filtered_genes': self.counts_filtered.shape[0] if self.counts_filtered is not None else None,
                'samples': self.counts_raw.shape[1] if self.counts_raw is not None else None
            },
            'de_analysis': {
                'significant_genes': len(filter_significant(
                    self.de_results['deseq2_like'],
                    padj_threshold=0.05,
                    log2fc_threshold=1.0
                )) if self.de_results is not None else None
            },
            'ml_results': self.ml_results.to_dict('records') if self.ml_results is not None else None
        }

        import json
        with open(self.results_dir / "pipeline_summary.json", 'w') as f:
            json.dump(summary, f, indent=2, default=str)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='RNA-seq Analysis Pipeline')
    parser.add_argument(
        '--config',
        type=str,
        default='configs/config.yaml',
        help='Path to configuration file'
    )
    parser.add_argument(
        '--step',
        type=str,
        choices=['all', 'load', 'normalize', 'batch', 'de', 'features', 'ml', 'validate'],
        default='all',
        help='Pipeline step to run'
    )

    args = parser.parse_args()

    # Resolve config path
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = Path(__file__).parent.parent / config_path

    # Initialize pipeline
    pipeline = RNAseqPipeline(str(config_path))

    # Run requested step
    if args.step == 'all':
        pipeline.run_full_pipeline()
    elif args.step == 'load':
        pipeline.step1_load_data()
    elif args.step == 'normalize':
        pipeline.step2_filter_and_normalize()
    elif args.step == 'batch':
        pipeline.step3_batch_correction()
    elif args.step == 'de':
        pipeline.step4_differential_expression()
    elif args.step == 'features':
        pipeline.step5_feature_engineering()
    elif args.step == 'ml':
        pipeline.step6_train_models()
    elif args.step == 'validate':
        pipeline.step7_validation()


if __name__ == "__main__":
    main()
