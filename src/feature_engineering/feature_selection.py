"""
Feature Engineering and Selection for RNA-seq ML
=================================================

This module implements various feature selection methods:
1. Variance-based filtering
2. Differential expression-based selection
3. Mutual information
4. Lasso/Elastic Net regularization
5. Random Forest importance
6. Recursive Feature Elimination (RFE)
"""

import pandas as pd
import numpy as np
from sklearn.feature_selection import (
    VarianceThreshold,
    SelectKBest,
    mutual_info_classif,
    RFE
)
from sklearn.linear_model import LassoCV, ElasticNetCV, LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder
from typing import Optional, List, Tuple, Dict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class FeatureSelector:
    """Feature selection for RNA-seq data."""

    def __init__(
        self,
        expression_data: pd.DataFrame,
        labels: pd.Series,
        gene_names: Optional[List[str]] = None
    ):
        """
        Initialize feature selector.

        Parameters
        ----------
        expression_data : pd.DataFrame
            Expression matrix (samples x genes)
        labels : pd.Series
            Class labels for each sample
        gene_names : List[str], optional
            Gene names (defaults to column names)
        """
        self.X = expression_data
        self.y = labels

        # Encode labels if needed
        if self.y.dtype == object:
            self.le = LabelEncoder()
            self.y_encoded = pd.Series(
                self.le.fit_transform(self.y),
                index=self.y.index
            )
        else:
            self.y_encoded = self.y
            self.le = None

        self.gene_names = gene_names or expression_data.columns.tolist()
        self.selected_features = {}

    def variance_filter(
        self,
        top_n: Optional[int] = None,
        threshold: Optional[float] = None
    ) -> List[str]:
        """
        Select features based on variance.

        Parameters
        ----------
        top_n : int, optional
            Select top N highest variance genes
        threshold : float, optional
            Minimum variance threshold

        Returns
        -------
        List[str]
            Selected gene names
        """
        logger.info("Running variance-based selection")

        variances = self.X.var(axis=0)

        if top_n is not None:
            selected_genes = variances.nlargest(top_n).index.tolist()
        elif threshold is not None:
            selected_genes = variances[variances > threshold].index.tolist()
        else:
            # Default: top 5000 or all if fewer
            top_n = min(5000, len(self.gene_names))
            selected_genes = variances.nlargest(top_n).index.tolist()

        self.selected_features['variance'] = selected_genes
        logger.info(f"Selected {len(selected_genes)} genes by variance")

        return selected_genes

    def de_filter(
        self,
        de_results: pd.DataFrame,
        padj_threshold: float = 0.05,
        log2fc_threshold: float = 1.0,
        gene_col: str = 'gene'
    ) -> List[str]:
        """
        Select features based on differential expression results.

        Parameters
        ----------
        de_results : pd.DataFrame
            DE results with gene, padj, and log2FC columns
        padj_threshold : float
            Adjusted p-value threshold
        log2fc_threshold : float
            Absolute log2 fold change threshold
        gene_col : str
            Column name for gene identifiers

        Returns
        -------
        List[str]
            Selected gene names
        """
        logger.info("Filtering by DE results")

        # Find log2FC column
        fc_col = [c for c in de_results.columns if 'log' in c.lower() and 'fc' in c.lower()]
        if not fc_col:
            fc_col = [c for c in de_results.columns if 'fold' in c.lower()]
        fc_col = fc_col[0] if fc_col else 'log2FoldChange'

        # Filter significant genes
        mask = (
            (de_results['padj'] < padj_threshold) &
            (abs(de_results[fc_col]) > log2fc_threshold)
        )

        selected_genes = de_results.loc[mask, gene_col].tolist()

        # Keep only genes in our expression data
        selected_genes = [g for g in selected_genes if g in self.X.columns]

        self.selected_features['de'] = selected_genes
        logger.info(f"Selected {len(selected_genes)} DE genes")

        return selected_genes

    def mutual_information(
        self,
        n_features: int = 100,
        random_state: int = 42
    ) -> List[str]:
        """
        Select features using mutual information.

        Parameters
        ----------
        n_features : int
            Number of features to select
        random_state : int
            Random state for reproducibility

        Returns
        -------
        List[str]
            Selected gene names
        """
        logger.info(f"Running mutual information selection (n={n_features})")

        # Scale features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(self.X)

        # Calculate MI scores
        mi_scores = mutual_info_classif(
            X_scaled,
            self.y_encoded,
            random_state=random_state
        )

        # Create DataFrame for selection
        mi_df = pd.DataFrame({
            'gene': self.X.columns,
            'mi_score': mi_scores
        }).sort_values('mi_score', ascending=False)

        selected_genes = mi_df.head(n_features)['gene'].tolist()

        self.selected_features['mutual_info'] = selected_genes
        logger.info(f"Selected {len(selected_genes)} genes by mutual information")

        return selected_genes

    def lasso_selection(
        self,
        n_features: Optional[int] = None,
        cv: int = 5,
        random_state: int = 42
    ) -> List[str]:
        """
        Select features using Lasso regularization.

        Parameters
        ----------
        n_features : int, optional
            Target number of features (adjusts alpha)
        cv : int
            Cross-validation folds
        random_state : int
            Random state

        Returns
        -------
        List[str]
            Selected gene names
        """
        logger.info("Running Lasso feature selection")

        # Scale features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(self.X)

        # Fit Lasso with cross-validation
        lasso = LassoCV(cv=cv, random_state=random_state, max_iter=10000)
        lasso.fit(X_scaled, self.y_encoded)

        # Get non-zero coefficients
        importance = np.abs(lasso.coef_)
        feature_importance = pd.DataFrame({
            'gene': self.X.columns,
            'importance': importance
        }).sort_values('importance', ascending=False)

        # Select features with non-zero coefficients
        if n_features is not None:
            selected_genes = feature_importance.head(n_features)['gene'].tolist()
        else:
            selected_genes = feature_importance[
                feature_importance['importance'] > 0
            ]['gene'].tolist()

        self.selected_features['lasso'] = selected_genes
        logger.info(f"Selected {len(selected_genes)} genes by Lasso")

        return selected_genes

    def random_forest_importance(
        self,
        n_features: int = 100,
        n_estimators: int = 100,
        random_state: int = 42
    ) -> List[str]:
        """
        Select features using Random Forest importance.

        Parameters
        ----------
        n_features : int
            Number of features to select
        n_estimators : int
            Number of trees
        random_state : int
            Random state

        Returns
        -------
        List[str]
            Selected gene names
        """
        logger.info(f"Running Random Forest importance selection (n={n_features})")

        rf = RandomForestClassifier(
            n_estimators=n_estimators,
            random_state=random_state,
            n_jobs=-1
        )
        rf.fit(self.X, self.y_encoded)

        # Get feature importances
        importance = pd.DataFrame({
            'gene': self.X.columns,
            'importance': rf.feature_importances_
        }).sort_values('importance', ascending=False)

        selected_genes = importance.head(n_features)['gene'].tolist()

        self.selected_features['random_forest'] = selected_genes
        logger.info(f"Selected {len(selected_genes)} genes by Random Forest")

        return selected_genes

    def recursive_feature_elimination(
        self,
        n_features: int = 100,
        step: float = 0.1,
        random_state: int = 42
    ) -> List[str]:
        """
        Select features using Recursive Feature Elimination.

        Parameters
        ----------
        n_features : int
            Number of features to select
        step : float
            Fraction of features to remove at each iteration
        random_state : int
            Random state

        Returns
        -------
        List[str]
            Selected gene names
        """
        logger.info(f"Running RFE selection (n={n_features})")

        # Use logistic regression as base estimator
        estimator = LogisticRegression(
            penalty='l2',
            solver='lbfgs',
            max_iter=1000,
            random_state=random_state
        )

        # Scale features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(self.X)

        rfe = RFE(
            estimator=estimator,
            n_features_to_select=n_features,
            step=step
        )
        rfe.fit(X_scaled, self.y_encoded)

        selected_genes = self.X.columns[rfe.support_].tolist()

        self.selected_features['rfe'] = selected_genes
        logger.info(f"Selected {len(selected_genes)} genes by RFE")

        return selected_genes

    def ensemble_selection(
        self,
        methods: List[str] = None,
        min_votes: int = 2,
        n_features: int = 100
    ) -> List[str]:
        """
        Select features that appear in multiple selection methods.

        Parameters
        ----------
        methods : List[str], optional
            Methods to include (defaults to all available)
        min_votes : int
            Minimum number of methods selecting a feature
        n_features : int
            Number of features to select per method

        Returns
        -------
        List[str]
            Selected gene names
        """
        logger.info("Running ensemble feature selection")

        if methods is None:
            methods = ['variance', 'mutual_info', 'random_forest']

        # Run each method if not already done
        for method in methods:
            if method not in self.selected_features:
                if method == 'variance':
                    self.variance_filter(top_n=n_features)
                elif method == 'mutual_info':
                    self.mutual_information(n_features=n_features)
                elif method == 'random_forest':
                    self.random_forest_importance(n_features=n_features)
                elif method == 'lasso':
                    self.lasso_selection(n_features=n_features)
                elif method == 'rfe':
                    self.recursive_feature_elimination(n_features=n_features)

        # Count votes for each gene
        all_genes = set()
        for method in methods:
            all_genes.update(self.selected_features.get(method, []))

        votes = {}
        for gene in all_genes:
            votes[gene] = sum(
                1 for method in methods
                if gene in self.selected_features.get(method, [])
            )

        # Select genes with enough votes
        votes_df = pd.DataFrame({
            'gene': list(votes.keys()),
            'votes': list(votes.values())
        }).sort_values('votes', ascending=False)

        selected_genes = votes_df[votes_df['votes'] >= min_votes]['gene'].tolist()

        self.selected_features['ensemble'] = selected_genes
        logger.info(f"Selected {len(selected_genes)} genes by ensemble (min votes: {min_votes})")

        return selected_genes

    def get_selected_data(
        self,
        method: str = 'ensemble'
    ) -> Tuple[pd.DataFrame, pd.Series]:
        """
        Get expression data with only selected features.

        Parameters
        ----------
        method : str
            Feature selection method to use

        Returns
        -------
        Tuple[pd.DataFrame, pd.Series]
            Selected features and labels
        """
        if method not in self.selected_features:
            raise ValueError(f"Method '{method}' not run yet")

        selected_genes = self.selected_features[method]
        X_selected = self.X[selected_genes]

        return X_selected, self.y

    def get_summary(self) -> pd.DataFrame:
        """Get summary of all feature selection results."""
        summary = []

        for method, genes in self.selected_features.items():
            summary.append({
                'method': method,
                'n_features': len(genes)
            })

        return pd.DataFrame(summary)


def prepare_ml_data(
    expression_data: pd.DataFrame,
    metadata: pd.DataFrame,
    condition_col: str = 'condition'
) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Prepare data for ML: transpose and align with labels.

    Parameters
    ----------
    expression_data : pd.DataFrame
        Expression matrix (genes x samples)
    metadata : pd.DataFrame
        Sample metadata
    condition_col : str
        Column with class labels

    Returns
    -------
    Tuple[pd.DataFrame, pd.Series]
        X (samples x genes) and y (labels)
    """
    # Transpose: samples as rows, genes as columns
    X = expression_data.T

    # Align labels
    y = metadata.loc[X.index, condition_col]

    return X, y


def main():
    """Test feature selection."""
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from preprocessing.data_loader import RNAseqDataLoader
    from preprocessing.normalization import RNAseqNormalizer

    project_root = Path(__file__).parent.parent.parent
    counts_file = project_root / "data/raw/GSE313799_counts_matrix_iN.txt"
    db_path = project_root / "db/rnaseq_analysis.db"

    # Load and preprocess data
    loader = RNAseqDataLoader(str(counts_file), str(db_path))
    counts = loader.load_counts()
    metadata = loader.parse_sample_metadata()

    # Filter
    keep = (counts >= 10).sum(axis=1) >= 3
    filtered = counts.loc[keep]

    # Normalize
    normalizer = RNAseqNormalizer(filtered)
    log_cpm = normalizer.cpm(log=True)

    # Prepare for ML
    metadata = metadata.set_index('sample_id')
    X, y = prepare_ml_data(log_cpm, metadata, 'condition')

    print(f"Data shape: {X.shape}")
    print(f"Labels: {y.value_counts().to_dict()}")

    # Feature selection
    selector = FeatureSelector(X, y)

    # Run various methods
    selector.variance_filter(top_n=1000)
    selector.mutual_information(n_features=100)
    selector.random_forest_importance(n_features=100)

    # Ensemble selection
    selector.ensemble_selection(
        methods=['variance', 'mutual_info', 'random_forest'],
        min_votes=2
    )

    # Summary
    print("\n=== Feature Selection Summary ===")
    print(selector.get_summary())

    # Get selected data
    X_selected, y_selected = selector.get_selected_data('ensemble')
    print(f"\nSelected features shape: {X_selected.shape}")

    # Save selected features
    pd.DataFrame({'gene': selector.selected_features['ensemble']}).to_csv(
        project_root / "results/selected_features.csv",
        index=False
    )


if __name__ == "__main__":
    main()
