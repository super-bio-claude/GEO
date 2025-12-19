"""
R Integration Module
====================

Python wrapper for running R scripts:
- DESeq2 + apeglm analysis
- KEGG/GO pathway analysis
- GSVA analysis

This module handles:
1. Preparing input files for R
2. Running R scripts via subprocess
3. Loading and processing R output
"""

import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, List
import logging
import json
import shutil

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class RRunner:
    """Execute R scripts and manage results."""

    def __init__(self, project_root: str):
        """
        Initialize R runner.

        Parameters
        ----------
        project_root : str
            Path to project root directory
        """
        self.project_root = Path(project_root)
        self.r_scripts_dir = self.project_root / "src" / "R_scripts"
        self.results_dir = self.project_root / "results"

        # Check if R is available
        self.r_available = self._check_r_installation()

    def _check_r_installation(self) -> bool:
        """Check if R is installed and accessible."""
        try:
            result = subprocess.run(
                ["Rscript", "--version"],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                logger.info("R installation found")
                return True
        except FileNotFoundError:
            logger.warning("R/Rscript not found in PATH")
        return False

    def _run_r_script(
        self,
        script_name: str,
        args: List[str],
        timeout: int = 3600
    ) -> bool:
        """
        Run an R script with arguments.

        Parameters
        ----------
        script_name : str
            Name of the R script
        args : List[str]
            Command line arguments
        timeout : int
            Timeout in seconds

        Returns
        -------
        bool
            True if successful
        """
        if not self.r_available:
            logger.error("R is not available. Please install R and required packages.")
            return False

        script_path = self.r_scripts_dir / script_name

        if not script_path.exists():
            logger.error(f"R script not found: {script_path}")
            return False

        cmd = ["Rscript", str(script_path)] + args
        logger.info(f"Running: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=str(self.project_root)
            )

            # Print R output
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(result.stderr)

            if result.returncode == 0:
                logger.info(f"R script completed successfully")
                return True
            else:
                logger.error(f"R script failed with code {result.returncode}")
                return False

        except subprocess.TimeoutExpired:
            logger.error(f"R script timed out after {timeout} seconds")
            return False
        except Exception as e:
            logger.error(f"Error running R script: {e}")
            return False


class DESeq2Runner(RRunner):
    """Run DESeq2 analysis with apeglm shrinkage."""

    def prepare_metadata(
        self,
        metadata_df: pd.DataFrame,
        output_path: str
    ):
        """Prepare metadata file for R."""
        # Ensure required columns
        required_cols = ['condition']
        for col in required_cols:
            if col not in metadata_df.columns:
                raise ValueError(f"Missing required column: {col}")

        metadata_df.to_csv(output_path)
        logger.info(f"Saved metadata to {output_path}")

    def run_deseq2(
        self,
        counts_file: str,
        metadata_file: str,
        output_dir: str = None
    ) -> Optional[Dict[str, pd.DataFrame]]:
        """
        Run DESeq2 analysis with apeglm.

        Parameters
        ----------
        counts_file : str
            Path to counts matrix
        metadata_file : str
            Path to sample metadata
        output_dir : str, optional
            Output directory

        Returns
        -------
        Dict[str, pd.DataFrame]
            Dictionary of result DataFrames
        """
        if output_dir is None:
            output_dir = str(self.results_dir / "deseq2")

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Run R script
        success = self._run_r_script(
            "deseq2_analysis.R",
            [counts_file, metadata_file, output_dir]
        )

        if not success:
            return None

        # Load results
        results = {}
        result_files = {
            'apeglm': 'deseq2_results_apeglm.csv',
            'noshrink': 'deseq2_results_noshrink.csv',
            'significant': 'significant_genes.csv',
            'upregulated': 'upregulated_genes.csv',
            'downregulated': 'downregulated_genes.csv',
            'normalized_counts': 'normalized_counts.csv',
            'vst_counts': 'vst_counts.csv'
        }

        for key, filename in result_files.items():
            filepath = Path(output_dir) / filename
            if filepath.exists():
                results[key] = pd.read_csv(filepath)
                logger.info(f"Loaded {key}: {len(results[key])} rows")

        return results


class PathwayAnalysisRunner(RRunner):
    """Run KEGG/GO pathway analysis."""

    def run_pathway_analysis(
        self,
        de_results_file: str,
        output_dir: str = None
    ) -> Optional[Dict[str, pd.DataFrame]]:
        """
        Run pathway analysis.

        Parameters
        ----------
        de_results_file : str
            Path to DE results (must have 'gene', 'log2FoldChange', 'padj' columns)
        output_dir : str, optional
            Output directory

        Returns
        -------
        Dict[str, pd.DataFrame]
            Dictionary of pathway analysis results
        """
        if output_dir is None:
            output_dir = str(self.results_dir / "pathway")

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Run R script
        success = self._run_r_script(
            "pathway_analysis.R",
            [de_results_file, output_dir]
        )

        if not success:
            return None

        # Load results
        results = {}
        result_files = {
            'kegg_all': 'kegg_all_significant.csv',
            'kegg_up': 'kegg_upregulated.csv',
            'kegg_down': 'kegg_downregulated.csv',
            'gsea_kegg': 'gsea_kegg.csv',
            'go_bp': 'go_bp.csv',
            'go_mf': 'go_mf.csv',
            'go_cc': 'go_cc.csv',
            'reactome': 'reactome.csv',
            'summary': 'pathway_analysis_summary.csv'
        }

        for key, filename in result_files.items():
            filepath = Path(output_dir) / filename
            if filepath.exists():
                results[key] = pd.read_csv(filepath)
                logger.info(f"Loaded {key}: {len(results[key])} rows")

        return results


class GSVARunner(RRunner):
    """Run GSVA analysis."""

    def run_gsva(
        self,
        expression_file: str,
        metadata_file: str,
        output_dir: str = None
    ) -> Optional[Dict[str, pd.DataFrame]]:
        """
        Run GSVA analysis.

        Parameters
        ----------
        expression_file : str
            Path to normalized expression matrix (VST recommended)
        metadata_file : str
            Path to sample metadata
        output_dir : str, optional
            Output directory

        Returns
        -------
        Dict[str, pd.DataFrame]
            Dictionary of GSVA results
        """
        if output_dir is None:
            output_dir = str(self.results_dir / "gsva")

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Run R script
        success = self._run_r_script(
            "gsva_analysis.R",
            [expression_file, metadata_file, output_dir]
        )

        if not success:
            return None

        # Load results
        results = {}

        # Load GSVA scores
        score_files = ['hallmark', 'kegg', 'reactome', 'go_bp']
        for gs_type in score_files:
            filepath = Path(output_dir) / f"gsva_scores_{gs_type}.csv"
            if filepath.exists():
                results[f'scores_{gs_type}'] = pd.read_csv(filepath, index_col=0)

        # Load differential results
        for gs_type in score_files:
            filepath = Path(output_dir) / f"differential_{gs_type}.csv"
            if filepath.exists():
                results[f'diff_{gs_type}'] = pd.read_csv(filepath)

        # Load summary
        summary_path = Path(output_dir) / "gsva_summary.csv"
        if summary_path.exists():
            results['summary'] = pd.read_csv(summary_path)

        return results


class IntegratedAnalysis:
    """
    Integrated analysis combining Python preprocessing with R-based
    DESeq2, pathway analysis, and GSVA.
    """

    def __init__(self, project_root: str):
        self.project_root = Path(project_root)
        self.deseq2_runner = DESeq2Runner(project_root)
        self.pathway_runner = PathwayAnalysisRunner(project_root)
        self.gsva_runner = GSVARunner(project_root)

    def run_full_r_analysis(
        self,
        counts_file: str,
        metadata_df: pd.DataFrame
    ) -> Dict[str, Any]:
        """
        Run complete R-based analysis pipeline.

        Parameters
        ----------
        counts_file : str
            Path to raw counts matrix
        metadata_df : pd.DataFrame
            Sample metadata

        Returns
        -------
        Dict[str, Any]
            All analysis results
        """
        results = {}
        results_dir = self.project_root / "results"

        # Save metadata for R
        metadata_file = str(results_dir / "sample_metadata.csv")
        metadata_df.to_csv(metadata_file)

        # 1. Run DESeq2
        logger.info("=" * 50)
        logger.info("Running DESeq2 with apeglm")
        logger.info("=" * 50)

        deseq2_results = self.deseq2_runner.run_deseq2(
            counts_file=counts_file,
            metadata_file=metadata_file,
            output_dir=str(results_dir / "deseq2")
        )
        results['deseq2'] = deseq2_results

        if deseq2_results is None:
            logger.warning("DESeq2 analysis failed or R not available")
            return results

        # 2. Run Pathway Analysis
        logger.info("=" * 50)
        logger.info("Running KEGG/GO Pathway Analysis")
        logger.info("=" * 50)

        de_results_file = str(results_dir / "deseq2" / "deseq2_results_apeglm.csv")
        pathway_results = self.pathway_runner.run_pathway_analysis(
            de_results_file=de_results_file,
            output_dir=str(results_dir / "pathway")
        )
        results['pathway'] = pathway_results

        # 3. Run GSVA
        logger.info("=" * 50)
        logger.info("Running GSVA Analysis")
        logger.info("=" * 50)

        vst_file = str(results_dir / "deseq2" / "vst_counts.csv")
        gsva_results = self.gsva_runner.run_gsva(
            expression_file=vst_file,
            metadata_file=metadata_file,
            output_dir=str(results_dir / "gsva")
        )
        results['gsva'] = gsva_results

        # Summary
        logger.info("=" * 50)
        logger.info("Analysis Complete")
        logger.info("=" * 50)

        self._print_summary(results)

        return results

    def _print_summary(self, results: Dict[str, Any]):
        """Print analysis summary."""
        print("\n" + "=" * 50)
        print("ANALYSIS SUMMARY")
        print("=" * 50)

        if results.get('deseq2'):
            deseq2 = results['deseq2']
            if 'significant' in deseq2:
                print(f"\nDESeq2 Results:")
                print(f"  Significant genes (padj < 0.05): {len(deseq2['significant'])}")
                if 'upregulated' in deseq2:
                    print(f"  Upregulated: {len(deseq2['upregulated'])}")
                if 'downregulated' in deseq2:
                    print(f"  Downregulated: {len(deseq2['downregulated'])}")

        if results.get('pathway'):
            pathway = results['pathway']
            if 'summary' in pathway:
                print(f"\nPathway Analysis:")
                print(pathway['summary'].to_string(index=False))

        if results.get('gsva'):
            gsva = results['gsva']
            if 'summary' in gsva:
                print(f"\nGSVA Analysis:")
                print(gsva['summary'].to_string(index=False))


def check_r_packages():
    """Check if required R packages are installed."""
    required_packages = [
        'DESeq2', 'apeglm', 'clusterProfiler', 'GSVA',
        'org.Hs.eg.db', 'msigdbr', 'enrichplot', 'pathview'
    ]

    r_code = f"""
    packages <- c({', '.join([f'"{p}"' for p in required_packages])})
    installed <- sapply(packages, function(p) requireNamespace(p, quietly = TRUE))
    cat("Package Status:\\n")
    for (i in seq_along(packages)) {{
        status <- ifelse(installed[i], "OK", "MISSING")
        cat(sprintf("  %s: %s\\n", packages[i], status))
    }}
    if (!all(installed)) {{
        cat("\\nTo install missing packages, run in R:\\n")
        cat("if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')\\n")
        missing <- packages[!installed]
        cat(sprintf("BiocManager::install(c('%s'))\\n", paste(missing, collapse = "', '")))
    }}
    """

    try:
        result = subprocess.run(
            ["Rscript", "-e", r_code],
            capture_output=True,
            text=True
        )
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
    except FileNotFoundError:
        print("R/Rscript not found. Please install R first.")


def main():
    """Test R integration."""
    from pathlib import Path

    project_root = Path(__file__).parent.parent.parent

    print("Checking R packages...")
    check_r_packages()

    # Example usage
    print("\n" + "=" * 50)
    print("R Integration Module")
    print("=" * 50)
    print(f"\nProject root: {project_root}")
    print(f"R scripts: {project_root / 'src' / 'R_scripts'}")

    print("\nAvailable runners:")
    print("  - DESeq2Runner: DESeq2 + apeglm analysis")
    print("  - PathwayAnalysisRunner: KEGG/GO/Reactome")
    print("  - GSVARunner: Gene Set Variation Analysis")
    print("  - IntegratedAnalysis: Full pipeline")


if __name__ == "__main__":
    main()
