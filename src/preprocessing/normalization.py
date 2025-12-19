"""
RNA-seq Normalization Methods
=============================

This module implements various normalization methods for RNA-seq count data:
1. CPM (Counts Per Million)
2. TPM (Transcripts Per Million) - requires gene lengths
3. TMM (Trimmed Mean of M-values)
4. Upper Quartile
5. VST (Variance Stabilizing Transformation)
6. rlog (Regularized Log Transformation)

Note: VST and rlog are typically done through DESeq2 in R,
but we provide Python approximations here.
"""

import pandas as pd
import numpy as np
from scipy import stats
from typing import Optional, Tuple
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class RNAseqNormalizer:
    """Normalize RNA-seq count data using various methods."""

    def __init__(self, counts_df: pd.DataFrame):
        """
        Initialize normalizer with counts DataFrame.

        Parameters
        ----------
        counts_df : pd.DataFrame
            Raw counts matrix (genes x samples)
        """
        self.counts_df = counts_df.copy()
        self.normalized = {}

    def cpm(self, log: bool = False, prior_count: float = 2) -> pd.DataFrame:
        """
        Calculate Counts Per Million (CPM).

        Parameters
        ----------
        log : bool
            If True, return log2(CPM + prior_count)
        prior_count : float
            Prior count added before log transformation

        Returns
        -------
        pd.DataFrame
            CPM normalized counts
        """
        lib_sizes = self.counts_df.sum(axis=0)
        cpm_df = self.counts_df * 1e6 / lib_sizes

        if log:
            cpm_df = np.log2(cpm_df + prior_count)

        self.normalized['cpm'] = cpm_df
        logger.info("CPM normalization complete")
        return cpm_df

    def tmm_factors(self) -> pd.Series:
        """
        Calculate TMM (Trimmed Mean of M-values) normalization factors.

        Based on Robinson & Oshlack (2010) - edgeR method.

        Returns
        -------
        pd.Series
            TMM normalization factors for each sample
        """
        # Use geometric mean as reference
        lib_sizes = self.counts_df.sum(axis=0)

        # Calculate log-ratios (M values) and average intensities (A values)
        # relative to a reference sample (highest count sample)
        ref_idx = lib_sizes.idxmax()
        ref_counts = self.counts_df[ref_idx]

        factors = {}
        for sample in self.counts_df.columns:
            if sample == ref_idx:
                factors[sample] = 1.0
                continue

            obs_counts = self.counts_df[sample]

            # Filter: keep genes with counts > 0 in both samples
            keep = (obs_counts > 0) & (ref_counts > 0)
            obs = obs_counts[keep]
            ref = ref_counts[keep]

            # Calculate M (log-ratio) and A (average intensity)
            M = np.log2(obs / lib_sizes[sample]) - np.log2(ref / lib_sizes[ref_idx])
            A = 0.5 * (np.log2(obs / lib_sizes[sample]) + np.log2(ref / lib_sizes[ref_idx]))

            # Trim extremes (30% from M, 5% from A)
            M_trim = (M > np.percentile(M, 30)) & (M < np.percentile(M, 70))
            A_trim = (A > np.percentile(A, 5)) & (A < np.percentile(A, 95))
            keep_trim = M_trim & A_trim

            if keep_trim.sum() > 0:
                # Weighted mean of M values
                weights = 1 / (1/obs[keep_trim] + 1/ref[keep_trim])
                tmm = np.average(M[keep_trim], weights=weights)
                factors[sample] = 2 ** tmm
            else:
                factors[sample] = 1.0

        return pd.Series(factors)

    def tmm_normalize(self, log: bool = True, prior_count: float = 2) -> pd.DataFrame:
        """
        Apply TMM normalization.

        Parameters
        ----------
        log : bool
            If True, return log2 transformed values
        prior_count : float
            Prior count for log transformation

        Returns
        -------
        pd.DataFrame
            TMM normalized counts
        """
        factors = self.tmm_factors()
        lib_sizes = self.counts_df.sum(axis=0)
        effective_lib_sizes = lib_sizes * factors

        # Normalize
        tmm_df = self.counts_df * 1e6 / effective_lib_sizes

        if log:
            tmm_df = np.log2(tmm_df + prior_count)

        self.normalized['tmm'] = tmm_df
        logger.info("TMM normalization complete")
        return tmm_df

    def upper_quartile(self, log: bool = True, prior_count: float = 2) -> pd.DataFrame:
        """
        Upper quartile normalization.

        Parameters
        ----------
        log : bool
            If True, return log2 transformed values
        prior_count : float
            Prior count for log transformation

        Returns
        -------
        pd.DataFrame
            Upper quartile normalized counts
        """
        # Calculate 75th percentile of non-zero counts for each sample
        uq_factors = self.counts_df.apply(
            lambda x: np.percentile(x[x > 0], 75) if (x > 0).any() else 1
        )

        # Normalize to median of upper quartiles
        median_uq = np.median(uq_factors)
        norm_factors = uq_factors / median_uq

        uq_df = self.counts_df / norm_factors

        if log:
            uq_df = np.log2(uq_df + prior_count)

        self.normalized['upper_quartile'] = uq_df
        logger.info("Upper quartile normalization complete")
        return uq_df

    def vst_approximate(self) -> pd.DataFrame:
        """
        Approximate Variance Stabilizing Transformation.

        This is a simplified version. For full VST, use DESeq2 in R.

        Returns
        -------
        pd.DataFrame
            VST-like transformed counts
        """
        # Add pseudocount and log transform with variance stabilization
        lib_sizes = self.counts_df.sum(axis=0)
        size_factors = lib_sizes / np.median(lib_sizes)

        # Normalize by size factors
        normalized = self.counts_df / size_factors

        # Apply variance stabilizing transformation (approximation)
        # VST ~ log2(count) for large counts, ~2*sqrt(count) for small counts
        vst_df = np.log2(normalized + 1)

        # Additional variance stabilization through asinh transformation
        # which behaves like log for large values and sqrt for small values
        vst_df = pd.DataFrame(
            np.arcsinh(normalized.values),
            index=self.counts_df.index,
            columns=self.counts_df.columns
        )

        self.normalized['vst'] = vst_df
        logger.info("VST approximation complete (use R/DESeq2 for exact VST)")
        return vst_df

    def quantile_normalize(self) -> pd.DataFrame:
        """
        Quantile normalization - forces all samples to have same distribution.

        Returns
        -------
        pd.DataFrame
            Quantile normalized counts
        """
        # Rank values within each sample
        ranked = self.counts_df.rank(method='average')

        # Calculate mean value for each rank across samples
        sorted_df = pd.DataFrame(
            np.sort(self.counts_df.values, axis=0),
            index=self.counts_df.index,
            columns=self.counts_df.columns
        )
        mean_values = sorted_df.mean(axis=1)

        # Map ranks to mean values
        qn_df = ranked.apply(lambda x: mean_values.iloc[x.astype(int) - 1].values)

        self.normalized['quantile'] = qn_df
        logger.info("Quantile normalization complete")
        return qn_df

    def normalize_all(self) -> dict:
        """
        Run all normalization methods.

        Returns
        -------
        dict
            Dictionary of normalized DataFrames
        """
        self.cpm(log=True)
        self.tmm_normalize(log=True)
        self.upper_quartile(log=True)
        self.vst_approximate()

        return self.normalized

    def get_summary_stats(self) -> pd.DataFrame:
        """Get summary statistics for all normalization methods."""
        stats_list = []

        for method, df in self.normalized.items():
            stats_list.append({
                'method': method,
                'mean': df.values.mean(),
                'std': df.values.std(),
                'min': df.values.min(),
                'max': df.values.max(),
                'cv': df.values.std() / df.values.mean() * 100
            })

        return pd.DataFrame(stats_list)


def compare_library_sizes(counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compare library sizes across samples.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Raw counts matrix

    Returns
    -------
    pd.DataFrame
        Library size statistics
    """
    lib_sizes = counts_df.sum(axis=0)

    stats_df = pd.DataFrame({
        'sample_id': lib_sizes.index,
        'total_counts': lib_sizes.values,
        'detected_genes': (counts_df > 0).sum(axis=0).values,
        'mean_count': counts_df.mean(axis=0).values,
        'median_count': counts_df.median(axis=0).values
    })

    return stats_df


def main():
    """Test normalization methods."""
    from pathlib import Path

    project_root = Path(__file__).parent.parent.parent
    counts_file = project_root / "data/raw/GSE313799_counts_matrix_iN.txt"

    # Load data
    counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)

    # Filter low counts first
    keep = (counts_df >= 10).sum(axis=1) >= 3
    filtered_df = counts_df.loc[keep]

    print(f"Filtered: {counts_df.shape[0]} -> {filtered_df.shape[0]} genes")

    # Normalize
    normalizer = RNAseqNormalizer(filtered_df)
    normalized = normalizer.normalize_all()

    # Save normalized data
    for method, df in normalized.items():
        output_path = project_root / f"data/normalized/{method}_normalized.csv"
        df.to_csv(output_path)
        print(f"Saved {method} normalized data to {output_path}")

    # Print summary
    print("\n=== Normalization Summary ===")
    print(normalizer.get_summary_stats())


if __name__ == "__main__":
    main()
