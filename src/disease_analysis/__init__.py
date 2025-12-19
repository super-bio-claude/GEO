"""
RNA-seq Disease Analysis System
================================

A comprehensive system for analyzing RNA-seq data to:
1. Identify disease similarity based on gene expression patterns
2. Predict disease progression stages
3. Generate clinical interpretation reports

Modules:
- disease_database: Build and manage disease signature databases
- similarity_engine: ChromaDB-based similarity search
- stage_predictor: Disease stage/severity prediction
- report_generator: Comprehensive analysis reports
- pipeline: Unified analysis pipeline
"""

from pathlib import Path

# Package metadata
__version__ = "1.0.0"
__author__ = "RNA-seq Analysis System"

# Default paths
DEFAULT_DATA_DIR = Path(__file__).parent.parent.parent / "data"
DEFAULT_RESULTS_DIR = Path(__file__).parent.parent.parent / "results"
