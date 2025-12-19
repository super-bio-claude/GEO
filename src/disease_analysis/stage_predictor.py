#!/usr/bin/env python3
"""
Disease Stage Predictor
=======================

Predicts disease progression stage based on:
1. Gene expression pattern matching
2. Severity score aggregation
3. Pathway activity analysis
4. Biomarker presence
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')


@dataclass
class StageAssessment:
    """Disease stage assessment result"""
    disease_category: str
    predicted_stage: str
    confidence: float
    severity_score: float
    supporting_evidence: List[str]
    biomarkers_present: List[str]
    biomarkers_absent: List[str]
    stage_indicators: Dict[str, float]
    recommendations: List[str]


class DiseaseStagePredictor:
    """Predict disease stage from gene expression"""

    def __init__(self):
        self._init_stage_markers()

    def _init_stage_markers(self):
        """Initialize stage-specific biomarkers"""

        # Cancer staging markers
        self.cancer_markers = {
            'early': {
                'positive': ['MKI67_low', 'BCL2', 'CDKN1A'],
                'negative': ['MMP9', 'VEGFA', 'VIM'],
                'pathways': ['cell_cycle_arrest', 'apoptosis']
            },
            'intermediate': {
                'positive': ['MKI67', 'CCND1', 'AURKA'],
                'negative': ['CDH1_loss'],
                'pathways': ['proliferation', 'angiogenesis_early']
            },
            'advanced': {
                'positive': ['MMP2', 'MMP9', 'VEGFA', 'VIM', 'SNAI1', 'TWIST1'],
                'negative': ['CDH1'],
                'pathways': ['emt', 'metastasis', 'angiogenesis']
            }
        }

        # Neurodegeneration staging (including Alzheimer's-specific markers)
        self.neuro_markers = {
            'preclinical': {
                'positive': ['BDNF', 'NTRK2', 'NPAS4', 'ARC', 'EGR1'],
                'negative': ['GFAP', 'AIF1', 'C1QA'],
                'pathways': ['neuroprotection', 'synaptic_plasticity']
            },
            'early': {
                'positive': ['APP', 'MAPT', 'SNCA', 'BACE1', 'PSEN1', 'CLU', 'BIN1'],
                'negative': ['NPAS4', 'ARC'],
                'pathways': ['protein_aggregation', 'amyloid']
            },
            'advanced': {
                'positive': ['GFAP', 'AIF1', 'IL1B', 'IL6', 'S100B', 'TREM2', 'CD68', 'C1QA', 'C3', 'TYROBP'],
                'negative': ['BDNF', 'NPAS4', 'FOSB', 'FOS', 'EGR1', 'EGR2', 'EGR4', 'ARC', 'SYP', 'DLG4'],
                'pathways': ['neuroinflammation', 'neurodegeneration', 'synaptic_loss']
            }
        }

        # Metabolic disease staging
        self.metabolic_markers = {
            'preclinical': {
                'positive': ['ADIPOQ', 'PPARG'],
                'negative': ['IL6', 'TNF'],
                'pathways': ['insulin_sensitivity']
            },
            'early': {
                'positive': ['INS', 'INSR', 'IRS1'],
                'negative': [],
                'pathways': ['glucose_metabolism']
            },
            'advanced': {
                'positive': ['IL6', 'TNF', 'CRP', 'TGFB1'],
                'negative': ['ADIPOQ'],
                'pathways': ['inflammation', 'fibrosis']
            }
        }

        # Aging markers
        self.aging_markers = {
            'early': {
                'positive': ['SIRT1', 'FOXO3', 'TERT'],
                'negative': ['CDKN2A'],
                'pathways': ['longevity']
            },
            'intermediate': {
                'positive': ['CDKN1A', 'TP53'],
                'negative': [],
                'pathways': ['cell_cycle_arrest']
            },
            'advanced': {
                'positive': ['CDKN2A', 'IL6', 'IL1B', 'SERPINE1'],
                'negative': ['SIRT1', 'TERT'],
                'pathways': ['senescence', 'inflammaging']
            }
        }

    def predict_stage(self,
                      expressed_genes: List[str],
                      fold_changes: Optional[Dict[str, float]] = None,
                      disease_category: str = 'cancer',
                      similarity_results: Optional[List] = None) -> StageAssessment:
        """Predict disease stage from expression data"""

        gene_set = set(expressed_genes)

        # Get up/down regulated genes
        up_genes = set()
        down_genes = set()
        if fold_changes:
            up_genes = {g for g, fc in fold_changes.items() if fc > 0.5}
            down_genes = {g for g, fc in fold_changes.items() if fc < -0.5}

        # Select marker set based on category
        if disease_category == 'cancer':
            markers = self.cancer_markers
        elif disease_category == 'neurological':
            markers = self.neuro_markers
        elif disease_category == 'metabolic':
            markers = self.metabolic_markers
        elif disease_category == 'aging':
            markers = self.aging_markers
        else:
            markers = self.cancer_markers  # Default

        # Score each stage
        stage_scores = {}
        stage_evidence = {}

        for stage, stage_markers in markers.items():
            score = 0
            evidence = []

            # Check positive markers (should be expressed, ideally upregulated)
            pos_markers = set(stage_markers['positive'])
            pos_found = gene_set.intersection(pos_markers)
            if pos_found:
                for g in pos_found:
                    if g in up_genes:
                        score += 2
                        evidence.append(f"{g} upregulated (stage marker)")
                    elif g in down_genes:
                        # Gene should be up but is down - negative signal
                        score -= 0.5
                    else:
                        score += 1
                        evidence.append(f"{g} expressed (stage marker)")

            # Check negative markers (should be absent or downregulated)
            neg_markers = set(stage_markers['negative'])

            for g in neg_markers:
                if g in down_genes:
                    # Gene is actively downregulated - strong signal
                    score += 2
                    evidence.append(f"{g} downregulated (expected for this stage)")
                elif g not in gene_set:
                    # Gene is absent - moderate signal
                    score += 0.5
                elif g in up_genes:
                    # Gene should be down but is up - negative signal
                    score -= 1

            # Normalize score
            max_possible = len(pos_markers) * 2 + len(neg_markers) * 2
            stage_scores[stage] = max(0, score / max_possible) if max_possible > 0 else 0
            stage_evidence[stage] = evidence

        # Incorporate similarity results if provided
        if similarity_results:
            for result in similarity_results[:10]:
                if hasattr(result, 'stage') and result.stage != 'unknown':
                    weight = result.combined_score if hasattr(result, 'combined_score') else 0.5
                    if result.stage in stage_scores:
                        stage_scores[result.stage] += weight * 0.3

        # Determine best stage
        best_stage = max(stage_scores, key=stage_scores.get)
        confidence = stage_scores[best_stage]

        # Calculate overall severity (0-1)
        severity = self._calculate_severity(
            stage_scores, best_stage, expressed_genes, fold_changes
        )

        # Get biomarker status
        biomarkers_present = []
        biomarkers_absent = []

        for marker in markers.get(best_stage, {}).get('positive', []):
            if marker in gene_set:
                biomarkers_present.append(marker)
            else:
                biomarkers_absent.append(marker)

        # Generate recommendations
        recommendations = self._generate_recommendations(
            disease_category, best_stage, severity, biomarkers_present
        )

        return StageAssessment(
            disease_category=disease_category,
            predicted_stage=best_stage,
            confidence=min(confidence, 1.0),
            severity_score=severity,
            supporting_evidence=stage_evidence.get(best_stage, []),
            biomarkers_present=biomarkers_present,
            biomarkers_absent=biomarkers_absent,
            stage_indicators=stage_scores,
            recommendations=recommendations
        )

    def _calculate_severity(self,
                           stage_scores: Dict[str, float],
                           best_stage: str,
                           genes: List[str],
                           fold_changes: Optional[Dict[str, float]]) -> float:
        """Calculate overall severity score"""

        # Base severity from stage
        stage_severity = {
            'preclinical': 0.1,
            'early': 0.3,
            'intermediate': 0.5,
            'advanced': 0.8,
            'terminal': 0.95
        }

        base_severity = stage_severity.get(best_stage, 0.5)

        # Adjust by confidence
        severity = base_severity * (0.7 + 0.3 * stage_scores.get(best_stage, 0.5))

        # Adjust by fold change magnitude if available
        if fold_changes:
            avg_abs_fc = np.mean([abs(fc) for fc in fold_changes.values()])
            # Higher fold changes suggest more severe dysregulation
            severity *= (1 + min(avg_abs_fc / 5, 0.3))

        return min(severity, 1.0)

    def _generate_recommendations(self,
                                  category: str,
                                  stage: str,
                                  severity: float,
                                  biomarkers: List[str]) -> List[str]:
        """Generate clinical recommendations"""

        recommendations = []

        if category == 'cancer':
            if stage == 'early':
                recommendations.append("Early stage detected - consider surgical intervention if applicable")
                recommendations.append("Regular monitoring recommended")
            elif stage == 'advanced':
                recommendations.append("Advanced stage markers detected - consider systemic therapy")
                if 'VEGFA' in biomarkers:
                    recommendations.append("Anti-angiogenic therapy may be beneficial (VEGFA positive)")
                if 'VIM' in biomarkers or 'SNAI1' in biomarkers:
                    recommendations.append("EMT signature present - monitor for metastasis")

        elif category == 'neurological':
            if stage == 'early':
                recommendations.append("Early neurodegeneration markers - neuroprotective interventions recommended")
            elif stage == 'advanced':
                recommendations.append("Neuroinflammation markers present - anti-inflammatory therapy consideration")

        elif category == 'metabolic':
            if stage == 'preclinical':
                recommendations.append("Prediabetic indicators - lifestyle modifications strongly recommended")
            elif stage == 'advanced':
                recommendations.append("Advanced metabolic dysfunction - comprehensive treatment needed")

        elif category == 'aging':
            if severity > 0.6:
                recommendations.append("Accelerated aging signature detected")
                recommendations.append("Consider senolytic interventions and lifestyle modifications")

        if severity > 0.7:
            recommendations.append(f"High severity score ({severity:.2f}) - urgent attention recommended")

        return recommendations

    def assess_multiple_categories(self,
                                   expressed_genes: List[str],
                                   fold_changes: Optional[Dict[str, float]] = None,
                                   similarity_results: Optional[List] = None) -> Dict[str, StageAssessment]:
        """Assess disease stage across multiple categories"""

        categories = ['cancer', 'neurological', 'metabolic', 'aging']
        assessments = {}

        for cat in categories:
            assessments[cat] = self.predict_stage(
                expressed_genes=expressed_genes,
                fold_changes=fold_changes,
                disease_category=cat,
                similarity_results=similarity_results
            )

        return assessments

    def get_primary_assessment(self,
                               assessments: Dict[str, StageAssessment]) -> Tuple[str, StageAssessment]:
        """Get the most relevant disease category assessment"""

        # Score by confidence and evidence
        scores = {}
        for cat, assessment in assessments.items():
            score = assessment.confidence * (1 + len(assessment.supporting_evidence) * 0.1)
            scores[cat] = score

        best_cat = max(scores, key=scores.get)
        return best_cat, assessments[best_cat]


if __name__ == "__main__":
    predictor = DiseaseStagePredictor()

    # Test with sample genes
    test_genes = ['BHLHE40', 'DDIT4', 'VEGFA', 'IL6', 'MMP9', 'VIM']
    test_fc = {'BHLHE40': 1.2, 'VEGFA': 0.8, 'IL6': 1.5, 'MMP9': 0.6}

    assessment = predictor.predict_stage(
        expressed_genes=test_genes,
        fold_changes=test_fc,
        disease_category='cancer'
    )

    print(f"Stage: {assessment.predicted_stage}")
    print(f"Confidence: {assessment.confidence:.2f}")
    print(f"Severity: {assessment.severity_score:.2f}")
    print(f"Evidence: {assessment.supporting_evidence}")
