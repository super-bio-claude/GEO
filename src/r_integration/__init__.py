"""
R Integration module for RNA-seq analysis.
"""

from .r_runner import (
    RRunner,
    DESeq2Runner,
    PathwayAnalysisRunner,
    GSVARunner,
    IntegratedAnalysis,
    check_r_packages
)

__all__ = [
    'RRunner',
    'DESeq2Runner',
    'PathwayAnalysisRunner',
    'GSVARunner',
    'IntegratedAnalysis',
    'check_r_packages'
]
