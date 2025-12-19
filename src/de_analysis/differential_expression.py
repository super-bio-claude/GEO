"""
Differential Expression Analysis
================================

This module implements methods for differential expression analysis:
1. DESeq2-like analysis (Python approximation using negative binomial)
2. edgeR-like analysis (Python approximation)
3. limma-voom-like analysis
4. Simple t-test/Wilcoxon approaches

Note: For publication-quality DE analysis, use R packages directly.
These Python implementations are for exploration and pipeline integration.
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.special import gammaln
from statsmodels.stats.multitest import multipletests
from typing import Optional, List, Tuple, Dict
import logging
import warnings

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DEAnalysis:
    """Differential Expression Analysis."""

    def __init__(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        condition_col: str = 'condition'
    ):
        """
        Initialize DE analysis.

        Parameters
        ----------
        counts : pd.DataFrame
            Count matrix (genes x samples)
        metadata : pd.DataFrame
            Sample metadata with condition column
        condition_col : str
            Column name for condition/group
        """
        self.counts = counts
        self.metadata = metadata
        self.condition_col = condition_col

        # Ensure sample order matches
        common_samples = list(set(counts.columns) & set(metadata.index))
        self.counts = counts[common_samples]
        self.metadata = metadata.loc[common_samples]

        logger.info(f"Initialized DE analysis with {len(common_samples)} samples")

    def run_deseq2_like(
        self,
        contrast: Tuple[str, str],
        alpha: float = 0.05
    ) -> pd.DataFrame:
        """
        Run DESeq2-like analysis using negative binomial model.

        This is a simplified implementation. For full DESeq2, use R.

        Parameters
        ----------
        contrast : Tuple[str, str]
            (numerator, denominator) for fold change calculation
        alpha : float
            Significance threshold for adjusted p-values

        Returns
        -------
        pd.DataFrame
            DE results with log2FC, pvalue, padj
        """
        logger.info(f"Running DESeq2-like analysis: {contrast[0]} vs {contrast[1]}")

        numerator, denominator = contrast

        # Get sample groups
        group1_samples = self.metadata[self.metadata[self.condition_col] == numerator].index
        group2_samples = self.metadata[self.metadata[self.condition_col] == denominator].index

        group1_counts = self.counts[group1_samples]
        group2_counts = self.counts[group2_samples]

        results = []

        # Calculate size factors (median of ratios method)
        all_counts = self.counts
        geometric_means = np.exp(np.log(all_counts + 1).mean(axis=1))
        size_factors = all_counts.divide(geometric_means, axis=0).median(axis=0)
        size_factors = size_factors / size_factors.mean()

        # Normalize counts
        norm_counts = self.counts.divide(size_factors, axis=1)

        for gene in self.counts.index:
            counts1 = group1_counts.loc[gene].values
            counts2 = group2_counts.loc[gene].values

            # Calculate mean and dispersion
            mean1 = counts1.mean()
            mean2 = counts2.mean()

            # Skip genes with very low counts
            if mean1 < 1 and mean2 < 1:
                continue

            # Base mean
            base_mean = (mean1 + mean2) / 2

            # Log2 fold change
            log2fc = np.log2((mean1 + 0.5) / (mean2 + 0.5))

            # Estimate dispersion (simplified)
            var1 = counts1.var() if len(counts1) > 1 else mean1
            var2 = counts2.var() if len(counts2) > 1 else mean2

            # Negative binomial test approximation using Wald test
            pooled_mean = (mean1 * len(counts1) + mean2 * len(counts2)) / (len(counts1) + len(counts2))
            pooled_var = ((var1 * (len(counts1) - 1) + var2 * (len(counts2) - 1)) /
                         (len(counts1) + len(counts2) - 2))

            # Standard error of log fold change
            se = np.sqrt(1/(mean1 + 0.5) + 1/(mean2 + 0.5) +
                        pooled_var/(mean1 + 0.5)**2/len(counts1) +
                        pooled_var/(mean2 + 0.5)**2/len(counts2))

            if se > 0:
                stat = log2fc / se
                pvalue = 2 * (1 - stats.norm.cdf(abs(stat)))
            else:
                stat = 0
                pvalue = 1.0

            results.append({
                'gene': gene,
                'baseMean': base_mean,
                'log2FoldChange': log2fc,
                'lfcSE': se,
                'stat': stat,
                'pvalue': pvalue
            })

        results_df = pd.DataFrame(results)

        # Multiple testing correction (Benjamini-Hochberg)
        _, padj, _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
        results_df['padj'] = padj

        # Sort by adjusted p-value
        results_df = results_df.sort_values('padj')

        logger.info(f"Found {(results_df['padj'] < alpha).sum()} significant genes (padj < {alpha})")

        return results_df

    def run_limma_voom_like(
        self,
        contrast: Tuple[str, str],
        alpha: float = 0.05
    ) -> pd.DataFrame:
        """
        Run limma-voom-like analysis using weighted linear regression.

        Parameters
        ----------
        contrast : Tuple[str, str]
            (numerator, denominator) conditions
        alpha : float
            Significance threshold

        Returns
        -------
        pd.DataFrame
            DE results
        """
        from sklearn.linear_model import LinearRegression

        logger.info(f"Running limma-voom-like analysis: {contrast[0]} vs {contrast[1]}")

        numerator, denominator = contrast

        # Calculate log-CPM
        lib_sizes = self.counts.sum(axis=0)
        log_cpm = np.log2(self.counts * 1e6 / lib_sizes + 0.5)

        # Create design matrix
        design = (self.metadata[self.condition_col] == numerator).astype(int)

        results = []

        for gene in log_cpm.index:
            y = log_cpm.loc[gene].values
            X = design.values.reshape(-1, 1)

            # Fit linear model
            model = LinearRegression()
            model.fit(X, y)

            coef = model.coef_[0]

            # Calculate statistics
            y_pred = model.predict(X)
            residuals = y - y_pred
            mse = np.sum(residuals**2) / (len(y) - 2)

            # Standard error of coefficient
            X_centered = X - X.mean()
            se = np.sqrt(mse / np.sum(X_centered**2))

            if se > 0:
                t_stat = coef / se
                pvalue = 2 * (1 - stats.t.cdf(abs(t_stat), df=len(y)-2))
            else:
                t_stat = 0
                pvalue = 1.0

            results.append({
                'gene': gene,
                'logFC': coef,
                'AveExpr': y.mean(),
                't': t_stat,
                'pvalue': pvalue
            })

        results_df = pd.DataFrame(results)

        # Multiple testing correction
        _, padj, _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
        results_df['padj'] = padj

        # Add B statistic (log-odds of differential expression) - simplified
        results_df['B'] = np.log((1 - results_df['padj'] + 1e-10) / (results_df['padj'] + 1e-10))

        results_df = results_df.sort_values('padj')

        logger.info(f"Found {(results_df['padj'] < alpha).sum()} significant genes")

        return results_df

    def run_simple_test(
        self,
        contrast: Tuple[str, str],
        test_type: str = 'ttest',
        alpha: float = 0.05
    ) -> pd.DataFrame:
        """
        Run simple statistical test (t-test or Wilcoxon).

        Parameters
        ----------
        contrast : Tuple[str, str]
            (numerator, denominator) conditions
        test_type : str
            'ttest' or 'wilcoxon'
        alpha : float
            Significance threshold

        Returns
        -------
        pd.DataFrame
            Test results
        """
        logger.info(f"Running {test_type}: {contrast[0]} vs {contrast[1]}")

        numerator, denominator = contrast

        group1_samples = self.metadata[self.metadata[self.condition_col] == numerator].index
        group2_samples = self.metadata[self.metadata[self.condition_col] == denominator].index

        # Use log-transformed data
        log_counts = np.log2(self.counts + 1)

        results = []

        for gene in self.counts.index:
            group1 = log_counts.loc[gene, group1_samples].values
            group2 = log_counts.loc[gene, group2_samples].values

            log2fc = group1.mean() - group2.mean()
            mean_expr = (group1.mean() + group2.mean()) / 2

            if test_type == 'ttest':
                stat, pvalue = stats.ttest_ind(group1, group2)
            elif test_type == 'wilcoxon':
                try:
                    stat, pvalue = stats.mannwhitneyu(group1, group2, alternative='two-sided')
                except ValueError:
                    stat, pvalue = 0, 1.0
            else:
                raise ValueError(f"Unknown test type: {test_type}")

            results.append({
                'gene': gene,
                'log2FC': log2fc,
                'meanExpr': mean_expr,
                'statistic': stat,
                'pvalue': pvalue if not np.isnan(pvalue) else 1.0
            })

        results_df = pd.DataFrame(results)

        # Multiple testing correction
        _, padj, _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
        results_df['padj'] = padj

        results_df = results_df.sort_values('padj')

        logger.info(f"Found {(results_df['padj'] < alpha).sum()} significant genes")

        return results_df

    def run_all_methods(
        self,
        contrast: Tuple[str, str],
        alpha: float = 0.05
    ) -> Dict[str, pd.DataFrame]:
        """
        Run all DE methods and return results.

        Parameters
        ----------
        contrast : Tuple[str, str]
            Contrast for comparison
        alpha : float
            Significance threshold

        Returns
        -------
        Dict[str, pd.DataFrame]
            Results from each method
        """
        results = {}

        results['deseq2_like'] = self.run_deseq2_like(contrast, alpha)
        results['limma_voom_like'] = self.run_limma_voom_like(contrast, alpha)
        results['ttest'] = self.run_simple_test(contrast, 'ttest', alpha)
        results['wilcoxon'] = self.run_simple_test(contrast, 'wilcoxon', alpha)

        return results

    def compare_methods(
        self,
        results: Dict[str, pd.DataFrame],
        alpha: float = 0.05
    ) -> pd.DataFrame:
        """
        Compare results from different DE methods.

        Parameters
        ----------
        results : Dict[str, pd.DataFrame]
            Results from different methods
        alpha : float
            Significance threshold

        Returns
        -------
        pd.DataFrame
            Comparison summary
        """
        comparison = []

        for method, df in results.items():
            n_sig = (df['padj'] < alpha).sum()
            comparison.append({
                'method': method,
                'n_significant': n_sig,
                'n_up': ((df['padj'] < alpha) & (df.iloc[:, 1] > 0)).sum(),
                'n_down': ((df['padj'] < alpha) & (df.iloc[:, 1] < 0)).sum()
            })

        return pd.DataFrame(comparison)


def get_top_genes(
    de_results: pd.DataFrame,
    n_top: int = 50,
    by: str = 'padj'
) -> pd.DataFrame:
    """Get top differentially expressed genes."""
    return de_results.nsmallest(n_top, by)


def filter_significant(
    de_results: pd.DataFrame,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 1.0
) -> pd.DataFrame:
    """Filter for significant genes based on padj and log2FC."""
    fc_col = [c for c in de_results.columns if 'log' in c.lower() and 'fc' in c.lower()][0]

    mask = (
        (de_results['padj'] < padj_threshold) &
        (abs(de_results[fc_col]) > log2fc_threshold)
    )

    return de_results[mask]


def main():
    """Test DE analysis."""
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from preprocessing.data_loader import RNAseqDataLoader

    project_root = Path(__file__).parent.parent.parent
    counts_file = project_root / "data/raw/GSE313799_counts_matrix_iN.txt"
    db_path = project_root / "db/rnaseq_analysis.db"

    # Load data
    loader = RNAseqDataLoader(str(counts_file), str(db_path))
    counts = loader.load_counts()
    metadata = loader.parse_sample_metadata()

    # Filter low counts
    keep = (counts >= 10).sum(axis=1) >= 3
    filtered = counts.loc[keep]

    # Set index for metadata
    metadata = metadata.set_index('sample_id')

    # Run DE analysis
    de = DEAnalysis(filtered, metadata, condition_col='condition')

    # Run all methods
    contrast = ('compress', 'control')
    results = de.run_all_methods(contrast)

    # Compare methods
    print("\n=== Method Comparison ===")
    print(de.compare_methods(results))

    # Save results
    for method, df in results.items():
        output_path = project_root / f"results/de_{method}.csv"
        output_path.parent.mkdir(exist_ok=True)
        df.to_csv(output_path, index=False)
        print(f"Saved {method} results to {output_path}")

    # Show top genes
    print("\n=== Top 10 DE Genes (DESeq2-like) ===")
    print(get_top_genes(results['deseq2_like'], n_top=10))


if __name__ == "__main__":
    main()
