"""
Model Validation and Evaluation
===============================

This module provides comprehensive model evaluation:
1. Cross-validation strategies
2. Performance metrics
3. Visualization of results
4. Statistical comparison of models
5. Learning curves
6. Permutation tests
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import (
    cross_val_score,
    StratifiedKFold,
    LeaveOneOut,
    learning_curve,
    permutation_test_score
)
from sklearn.metrics import (
    accuracy_score,
    precision_recall_fscore_support,
    roc_curve,
    auc,
    precision_recall_curve,
    confusion_matrix,
    classification_report,
    matthews_corrcoef
)
from sklearn.preprocessing import LabelEncoder
from scipy import stats
from typing import Dict, List, Tuple, Optional, Any
import logging
import json

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ModelEvaluator:
    """Comprehensive model evaluation."""

    def __init__(
        self,
        model,
        X: pd.DataFrame,
        y: pd.Series,
        random_state: int = 42
    ):
        """
        Initialize evaluator.

        Parameters
        ----------
        model : sklearn estimator
            Trained model
        X : pd.DataFrame
            Feature matrix
        y : pd.Series
            Labels
        random_state : int
            Random seed
        """
        self.model = model
        self.X = X
        self.y = y
        self.random_state = random_state

        # Encode labels
        self.le = LabelEncoder()
        self.y_encoded = self.le.fit_transform(y)

    def cross_validation(
        self,
        cv_strategy: str = 'stratified',
        n_splits: int = 5,
        scoring: str = 'accuracy'
    ) -> Dict[str, Any]:
        """
        Perform cross-validation.

        Parameters
        ----------
        cv_strategy : str
            'stratified', 'loo' (leave-one-out), or 'kfold'
        n_splits : int
            Number of CV folds
        scoring : str
            Scoring metric

        Returns
        -------
        Dict[str, Any]
            CV results
        """
        logger.info(f"Running {cv_strategy} cross-validation")

        if cv_strategy == 'stratified':
            cv = StratifiedKFold(
                n_splits=n_splits,
                shuffle=True,
                random_state=self.random_state
            )
        elif cv_strategy == 'loo':
            cv = LeaveOneOut()
        else:
            raise ValueError(f"Unknown CV strategy: {cv_strategy}")

        scores = cross_val_score(
            self.model,
            self.X,
            self.y_encoded,
            cv=cv,
            scoring=scoring
        )

        results = {
            'cv_strategy': cv_strategy,
            'n_splits': n_splits if cv_strategy != 'loo' else len(self.y),
            'scoring': scoring,
            'scores': scores.tolist(),
            'mean_score': scores.mean(),
            'std_score': scores.std(),
            'ci_95_lower': scores.mean() - 1.96 * scores.std() / np.sqrt(len(scores)),
            'ci_95_upper': scores.mean() + 1.96 * scores.std() / np.sqrt(len(scores))
        }

        logger.info(f"CV Score: {results['mean_score']:.4f} +/- {results['std_score']:.4f}")

        return results

    def compute_metrics(
        self,
        y_true: np.ndarray,
        y_pred: np.ndarray,
        y_prob: Optional[np.ndarray] = None
    ) -> Dict[str, Any]:
        """
        Compute comprehensive evaluation metrics.

        Parameters
        ----------
        y_true : np.ndarray
            True labels
        y_pred : np.ndarray
            Predicted labels
        y_prob : np.ndarray, optional
            Prediction probabilities

        Returns
        -------
        Dict[str, Any]
            Evaluation metrics
        """
        metrics = {
            'accuracy': accuracy_score(y_true, y_pred),
            'mcc': matthews_corrcoef(y_true, y_pred)
        }

        # Per-class metrics
        precision, recall, f1, support = precision_recall_fscore_support(
            y_true, y_pred, average=None
        )

        metrics['precision_per_class'] = precision.tolist()
        metrics['recall_per_class'] = recall.tolist()
        metrics['f1_per_class'] = f1.tolist()
        metrics['support_per_class'] = support.tolist()

        # Weighted averages
        precision_w, recall_w, f1_w, _ = precision_recall_fscore_support(
            y_true, y_pred, average='weighted'
        )
        metrics['precision_weighted'] = precision_w
        metrics['recall_weighted'] = recall_w
        metrics['f1_weighted'] = f1_w

        # Confusion matrix
        metrics['confusion_matrix'] = confusion_matrix(y_true, y_pred).tolist()

        # ROC-AUC if probabilities available
        if y_prob is not None:
            if len(np.unique(y_true)) == 2:
                fpr, tpr, _ = roc_curve(y_true, y_prob[:, 1])
                metrics['auc_roc'] = auc(fpr, tpr)
                metrics['fpr'] = fpr.tolist()
                metrics['tpr'] = tpr.tolist()

                # Precision-Recall curve
                precision_curve, recall_curve, _ = precision_recall_curve(
                    y_true, y_prob[:, 1]
                )
                metrics['auc_pr'] = auc(recall_curve, precision_curve)

        return metrics

    def learning_curve_analysis(
        self,
        train_sizes: np.ndarray = None,
        cv: int = 5
    ) -> Dict[str, Any]:
        """
        Generate learning curve data.

        Parameters
        ----------
        train_sizes : np.ndarray, optional
            Training set sizes to evaluate
        cv : int
            Cross-validation folds

        Returns
        -------
        Dict[str, Any]
            Learning curve results
        """
        logger.info("Computing learning curve")

        if train_sizes is None:
            train_sizes = np.linspace(0.1, 1.0, 10)

        train_sizes_abs, train_scores, test_scores = learning_curve(
            self.model,
            self.X,
            self.y_encoded,
            train_sizes=train_sizes,
            cv=cv,
            n_jobs=-1,
            random_state=self.random_state
        )

        results = {
            'train_sizes': train_sizes_abs.tolist(),
            'train_scores_mean': train_scores.mean(axis=1).tolist(),
            'train_scores_std': train_scores.std(axis=1).tolist(),
            'test_scores_mean': test_scores.mean(axis=1).tolist(),
            'test_scores_std': test_scores.std(axis=1).tolist()
        }

        return results

    def permutation_test(
        self,
        n_permutations: int = 100,
        cv: int = 5
    ) -> Dict[str, Any]:
        """
        Perform permutation test for statistical significance.

        Parameters
        ----------
        n_permutations : int
            Number of permutations
        cv : int
            Cross-validation folds

        Returns
        -------
        Dict[str, Any]
            Permutation test results
        """
        logger.info(f"Running permutation test ({n_permutations} permutations)")

        score, perm_scores, pvalue = permutation_test_score(
            self.model,
            self.X,
            self.y_encoded,
            cv=cv,
            n_permutations=n_permutations,
            n_jobs=-1,
            random_state=self.random_state
        )

        results = {
            'score': score,
            'permutation_scores_mean': perm_scores.mean(),
            'permutation_scores_std': perm_scores.std(),
            'pvalue': pvalue,
            'significant': pvalue < 0.05
        }

        logger.info(f"Permutation test p-value: {pvalue:.4f}")

        return results


def compare_models(
    models: Dict[str, Any],
    X: pd.DataFrame,
    y: pd.Series,
    cv: int = 5,
    random_state: int = 42
) -> pd.DataFrame:
    """
    Compare multiple models using paired statistical tests.

    Parameters
    ----------
    models : Dict[str, Any]
        Dictionary of model name -> model
    X : pd.DataFrame
        Features
    y : pd.Series
        Labels
    cv : int
        Cross-validation folds
    random_state : int
        Random seed

    Returns
    -------
    pd.DataFrame
        Model comparison results
    """
    logger.info("Comparing models")

    le = LabelEncoder()
    y_encoded = le.fit_transform(y)

    cv_split = StratifiedKFold(
        n_splits=cv,
        shuffle=True,
        random_state=random_state
    )

    # Collect CV scores for each model
    all_scores = {}
    for name, model in models.items():
        scores = cross_val_score(model, X, y_encoded, cv=cv_split, scoring='accuracy')
        all_scores[name] = scores

    # Create results DataFrame
    results = []
    model_names = list(models.keys())

    for name, scores in all_scores.items():
        results.append({
            'model': name,
            'mean_accuracy': scores.mean(),
            'std_accuracy': scores.std(),
            'min_accuracy': scores.min(),
            'max_accuracy': scores.max()
        })

    results_df = pd.DataFrame(results).sort_values('mean_accuracy', ascending=False)

    # Pairwise statistical tests
    print("\n=== Pairwise Model Comparison (Wilcoxon signed-rank test) ===")
    for i, model1 in enumerate(model_names):
        for model2 in model_names[i+1:]:
            stat, pvalue = stats.wilcoxon(all_scores[model1], all_scores[model2])
            sig = "*" if pvalue < 0.05 else ""
            print(f"{model1} vs {model2}: p={pvalue:.4f} {sig}")

    return results_df


def bootstrap_confidence_interval(
    model,
    X: pd.DataFrame,
    y: pd.Series,
    n_bootstrap: int = 1000,
    confidence: float = 0.95,
    random_state: int = 42
) -> Dict[str, float]:
    """
    Calculate bootstrap confidence intervals for model performance.

    Parameters
    ----------
    model : sklearn estimator
        Trained model
    X : pd.DataFrame
        Features
    y : pd.Series
        Labels
    n_bootstrap : int
        Number of bootstrap iterations
    confidence : float
        Confidence level
    random_state : int
        Random seed

    Returns
    -------
    Dict[str, float]
        Bootstrap CI results
    """
    logger.info(f"Computing bootstrap CI ({n_bootstrap} iterations)")

    le = LabelEncoder()
    y_encoded = le.fit_transform(y)

    np.random.seed(random_state)
    n_samples = len(y)

    scores = []
    for _ in range(n_bootstrap):
        # Bootstrap sample
        indices = np.random.choice(n_samples, size=n_samples, replace=True)
        X_boot = X.iloc[indices]
        y_boot = y_encoded[indices]

        # Out-of-bag samples for evaluation
        oob_indices = list(set(range(n_samples)) - set(indices))
        if len(oob_indices) < 2:
            continue

        X_oob = X.iloc[oob_indices]
        y_oob = y_encoded[oob_indices]

        # Fit and evaluate
        model.fit(X_boot, y_boot)
        score = accuracy_score(y_oob, model.predict(X_oob))
        scores.append(score)

    scores = np.array(scores)
    alpha = 1 - confidence

    results = {
        'mean': scores.mean(),
        'std': scores.std(),
        'ci_lower': np.percentile(scores, alpha/2 * 100),
        'ci_upper': np.percentile(scores, (1 - alpha/2) * 100)
    }

    logger.info(f"Bootstrap CI: {results['ci_lower']:.4f} - {results['ci_upper']:.4f}")

    return results


def generate_evaluation_report(
    evaluator: ModelEvaluator,
    output_path: str
) -> Dict[str, Any]:
    """
    Generate comprehensive evaluation report.

    Parameters
    ----------
    evaluator : ModelEvaluator
        Initialized model evaluator
    output_path : str
        Path to save report

    Returns
    -------
    Dict[str, Any]
        Complete evaluation report
    """
    report = {
        'cross_validation': {},
        'learning_curve': {},
        'permutation_test': {}
    }

    # Cross-validation with different strategies
    for cv_type in ['stratified']:
        report['cross_validation'][cv_type] = evaluator.cross_validation(
            cv_strategy=cv_type,
            n_splits=5
        )

    # Learning curve
    report['learning_curve'] = evaluator.learning_curve_analysis()

    # Permutation test
    report['permutation_test'] = evaluator.permutation_test(n_permutations=100)

    # Save report
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

    logger.info(f"Saved evaluation report to {output_path}")

    return report


def main():
    """Test model evaluation."""
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from preprocessing.data_loader import RNAseqDataLoader
    from preprocessing.normalization import RNAseqNormalizer
    from feature_engineering.feature_selection import FeatureSelector, prepare_ml_data
    from ml_models.classifiers import RNAseqClassifier

    project_root = Path(__file__).parent.parent.parent
    counts_file = project_root / "data/raw/GSE313799_counts_matrix_iN.txt"
    db_path = project_root / "db/rnaseq_analysis.db"

    # Load and preprocess
    loader = RNAseqDataLoader(str(counts_file), str(db_path))
    counts = loader.load_counts()
    metadata = loader.parse_sample_metadata()

    # Filter and normalize
    keep = (counts >= 10).sum(axis=1) >= 3
    filtered = counts.loc[keep]
    normalizer = RNAseqNormalizer(filtered)
    log_cpm = normalizer.cpm(log=True)

    # Prepare data
    metadata = metadata.set_index('sample_id')
    X, y = prepare_ml_data(log_cpm, metadata, 'condition')

    # Feature selection
    selector = FeatureSelector(X, y)
    selector.variance_filter(top_n=500)
    X_selected, y_selected = selector.get_selected_data('variance')

    # Train a model
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    model = Pipeline([
        ('scaler', StandardScaler()),
        ('classifier', RandomForestClassifier(n_estimators=100, random_state=42))
    ])

    # Evaluate
    evaluator = ModelEvaluator(model, X_selected, y_selected)

    print("\n=== Cross-Validation Results ===")
    cv_results = evaluator.cross_validation('stratified', n_splits=5)
    print(f"Accuracy: {cv_results['mean_score']:.4f} +/- {cv_results['std_score']:.4f}")

    print("\n=== Permutation Test ===")
    perm_results = evaluator.permutation_test(n_permutations=100)
    print(f"P-value: {perm_results['pvalue']:.4f}")
    print(f"Significant: {perm_results['significant']}")

    # Save report
    report = generate_evaluation_report(
        evaluator,
        str(project_root / "results/evaluation_report.json")
    )


if __name__ == "__main__":
    main()
