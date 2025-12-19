"""
Disease Analysis Service
Integrates with the core disease analysis pipeline
"""
import sys
import uuid
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import json

# Add parent project to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent / "src"))

from disease_analysis.disease_database import DiseaseDatabaseBuilder
from disease_analysis.similarity_engine import DiseaseSimilarityEngine
from disease_analysis.stage_predictor import DiseaseStagePredictor
from disease_analysis.report_generator import DiseaseAnalysisReportGenerator

from app.models.schemas import (
    AnalysisResult, SignatureMatch, StageAssessment, BiomarkerAnalysis,
    DiseaseCategory, DiseaseStage, RiskLevel
)


class DiseaseAnalysisService:
    """Service for running disease analysis"""

    _instance = None
    _initialized = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if not self._initialized:
            self._setup()
            DiseaseAnalysisService._initialized = True

    def _setup(self):
        """Initialize analysis components"""
        base_path = Path(__file__).parent.parent.parent.parent.parent

        self.db_builder = DiseaseDatabaseBuilder(
            str(base_path / "data" / "disease_db")
        )
        self.similarity_engine = DiseaseSimilarityEngine(
            db_path=str(base_path / "data" / "chromadb_disease")
        )
        self.stage_predictor = DiseaseStagePredictor()
        self.report_generator = DiseaseAnalysisReportGenerator(
            str(base_path / "web" / "backend" / "reports")
        )

        # Load database and index
        self._load_database()

    def _load_database(self):
        """Load or build disease database"""
        db_file = Path(self.db_builder.data_dir) / "disease_signatures_full.json"

        if db_file.exists():
            self.signatures = self.db_builder.load_database()
        else:
            self.signatures = self.db_builder.build_database()

        if not self.similarity_engine.load_index():
            self.similarity_engine.build_index(self.signatures)

    def get_signature_count(self) -> int:
        """Get total number of disease signatures"""
        return len(self.signatures) if self.signatures else 0

    async def analyze_de_results(
        self,
        de_data: pd.DataFrame,
        sample_id: str,
        padj_threshold: float = 0.05,
        lfc_threshold: float = 0.5
    ) -> AnalysisResult:
        """Run disease analysis on DE results"""

        analysis_id = str(uuid.uuid4())

        # Extract significant genes
        if 'padj' in de_data.columns:
            sig_mask = (
                (de_data['padj'] < padj_threshold) &
                (de_data['log2FoldChange'].abs() > lfc_threshold)
            )
            sig_df = de_data[sig_mask]
        else:
            sig_df = de_data.nlargest(100, 'log2FoldChange')

        # Get gene lists
        all_genes = de_data['gene'].tolist()
        sig_genes = sig_df['gene'].tolist() if len(sig_df) > 0 else all_genes[:100]
        fold_changes = dict(zip(de_data['gene'], de_data['log2FoldChange']))

        up_genes = de_data[de_data['log2FoldChange'] > lfc_threshold]['gene'].tolist()
        down_genes = de_data[de_data['log2FoldChange'] < -lfc_threshold]['gene'].tolist()

        # Query genes
        if len(sig_genes) < 20:
            de_data['abs_lfc'] = de_data['log2FoldChange'].abs()
            top_genes = de_data.nlargest(200, 'abs_lfc')['gene'].tolist()
            query_genes = list(set(sig_genes + top_genes))
        else:
            query_genes = sig_genes

        # DE summary
        de_summary = {
            'n_total': len(de_data),
            'n_significant': len(sig_genes),
            'n_upregulated': len([g for g in sig_genes if fold_changes.get(g, 0) > 0]),
            'n_downregulated': len([g for g in sig_genes if fold_changes.get(g, 0) < 0]),
            'top_genes': sig_genes[:20]
        }

        # Similarity search
        similarity_results = self.similarity_engine.search(
            query_genes=query_genes,
            query_up_genes=up_genes[:100],
            query_down_genes=down_genes[:100],
            fold_changes=fold_changes,
            top_k=50
        )

        # Stage prediction
        stage_assessments = self.stage_predictor.assess_multiple_categories(
            expressed_genes=query_genes,
            fold_changes=fold_changes,
            similarity_results=similarity_results
        )

        primary_cat, primary_assessment = self.stage_predictor.get_primary_assessment(
            stage_assessments
        )

        # Generate report
        report_file = self.report_generator.generate_report(
            sample_id=sample_id,
            similarity_results=similarity_results,
            stage_assessments=stage_assessments,
            de_summary=de_summary
        )

        # Build response
        result = self._build_analysis_result(
            analysis_id=analysis_id,
            sample_id=sample_id,
            similarity_results=similarity_results,
            stage_assessments=stage_assessments,
            primary_cat=primary_cat,
            primary_assessment=primary_assessment,
            de_summary=de_summary,
            report_file=report_file
        )

        return result

    def _build_analysis_result(
        self,
        analysis_id: str,
        sample_id: str,
        similarity_results: List,
        stage_assessments: Dict,
        primary_cat: str,
        primary_assessment,
        de_summary: Dict,
        report_file: str
    ) -> AnalysisResult:
        """Build AnalysisResult from raw results"""

        # Top signatures
        top_signatures = []
        for i, r in enumerate(similarity_results[:20], 1):
            top_signatures.append(SignatureMatch(
                rank=i,
                signature_id=r.signature_id,
                signature_name=r.signature_name,
                category=DiseaseCategory(r.category) if r.category in [e.value for e in DiseaseCategory] else DiseaseCategory.OTHER,
                stage=DiseaseStage(r.stage) if r.stage in [e.value for e in DiseaseStage] else DiseaseStage.UNKNOWN,
                score=round(r.combined_score, 3),
                gene_overlap=r.gene_overlap,
                overlapping_genes=r.overlapping_genes[:10]
            ))

        # Stage assessments
        stage_dict = {}
        for cat, assessment in stage_assessments.items():
            stage_dict[cat] = StageAssessment(
                category=DiseaseCategory(cat) if cat in [e.value for e in DiseaseCategory] else DiseaseCategory.OTHER,
                predicted_stage=DiseaseStage(assessment.predicted_stage) if assessment.predicted_stage in [e.value for e in DiseaseStage] else DiseaseStage.UNKNOWN,
                confidence=round(assessment.confidence, 3),
                severity_score=round(assessment.severity_score, 3)
            )

        # Biomarkers
        biomarkers = BiomarkerAnalysis(
            detected=primary_assessment.biomarkers_present,
            missing=primary_assessment.biomarkers_absent,
            top_de_genes=de_summary.get('top_genes', [])[:10]
        )

        # Risk level
        severity = primary_assessment.severity_score
        if severity > 0.7:
            risk_level = RiskLevel.HIGH
        elif severity > 0.4:
            risk_level = RiskLevel.MODERATE
        else:
            risk_level = RiskLevel.LOW

        return AnalysisResult(
            id=analysis_id,
            sample_id=sample_id,
            created_at=datetime.utcnow(),
            status="completed",
            primary_category=DiseaseCategory(primary_cat) if primary_cat in [e.value for e in DiseaseCategory] else DiseaseCategory.OTHER,
            predicted_stage=DiseaseStage(primary_assessment.predicted_stage) if primary_assessment.predicted_stage in [e.value for e in DiseaseStage] else DiseaseStage.UNKNOWN,
            severity_score=round(primary_assessment.severity_score, 3),
            confidence=round(primary_assessment.confidence, 3),
            risk_level=risk_level,
            total_genes=de_summary['n_total'],
            significant_genes=de_summary['n_significant'],
            upregulated=de_summary['n_upregulated'],
            downregulated=de_summary['n_downregulated'],
            top_signatures=top_signatures,
            stage_assessments=stage_dict,
            biomarkers=biomarkers,
            recommendations=primary_assessment.recommendations,
            supporting_evidence=primary_assessment.supporting_evidence,
            report_url=f"/api/v1/reports/{Path(report_file).name}"
        )

    def get_signatures(
        self,
        category: Optional[str] = None,
        limit: int = 100
    ) -> Tuple[List[Dict], Dict[str, int]]:
        """Get list of disease signatures"""

        signatures_list = []
        category_counts = {}

        for sig_id, sig in self.signatures.items():
            cat = sig.category if hasattr(sig, 'category') else sig.get('category', 'unknown')

            # Count by category
            category_counts[cat] = category_counts.get(cat, 0) + 1

            # Filter if needed
            if category and cat != category:
                continue

            genes = sig.genes if hasattr(sig, 'genes') else sig.get('genes', [])
            name = sig.name if hasattr(sig, 'name') else sig.get('name', sig_id)
            stage = sig.stage if hasattr(sig, 'stage') else sig.get('stage', 'unknown')
            desc = sig.description if hasattr(sig, 'description') else sig.get('description', '')

            signatures_list.append({
                'id': sig_id,
                'name': name,
                'category': cat,
                'stage': stage,
                'gene_count': len(genes),
                'description': desc[:200] if desc else None
            })

            if len(signatures_list) >= limit:
                break

        return signatures_list, category_counts


# Singleton instance
analysis_service = DiseaseAnalysisService()
