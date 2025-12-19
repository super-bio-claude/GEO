"""
Differential Expression Analysis module.
"""

from .differential_expression import (
    DEAnalysis,
    get_top_genes,
    filter_significant
)

__all__ = ['DEAnalysis', 'get_top_genes', 'filter_significant']
