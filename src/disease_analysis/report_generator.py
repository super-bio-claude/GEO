#!/usr/bin/env python3
"""
Disease Analysis Report Generator
==================================

Generates comprehensive analysis reports including:
1. Disease similarity analysis
2. Stage/progression assessment
3. Key biomarkers
4. Clinical interpretations
5. Visualizations
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False


class DiseaseAnalysisReportGenerator:
    """Generate comprehensive disease analysis reports"""

    def __init__(self, output_dir: str = "results/disease_analysis"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_report(self,
                        sample_id: str,
                        similarity_results: List,
                        stage_assessments: Dict,
                        de_summary: Dict,
                        pathway_results: Optional[Dict] = None) -> str:
        """Generate comprehensive analysis report"""

        report_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Create report sections
        sections = []

        # Header
        sections.append(self._generate_header(sample_id, report_time))

        # Executive Summary
        sections.append(self._generate_executive_summary(
            similarity_results, stage_assessments, de_summary
        ))

        # Disease Similarity Analysis
        sections.append(self._generate_similarity_section(similarity_results))

        # Disease Stage Assessment
        sections.append(self._generate_stage_section(stage_assessments))

        # Biomarker Analysis
        sections.append(self._generate_biomarker_section(
            stage_assessments, de_summary
        ))

        # Pathway Analysis
        if pathway_results:
            sections.append(self._generate_pathway_section(pathway_results))

        # Clinical Interpretation
        sections.append(self._generate_clinical_section(
            similarity_results, stage_assessments
        ))

        # Recommendations
        sections.append(self._generate_recommendations(stage_assessments))

        # Combine report
        full_report = "\n\n".join(sections)

        # Save report
        report_file = self.output_dir / f"disease_analysis_report_{sample_id}.md"
        with open(report_file, 'w') as f:
            f.write(full_report)

        # Generate visualizations
        if PLOTTING_AVAILABLE:
            self._generate_visualizations(
                sample_id, similarity_results, stage_assessments
            )

        # Save JSON summary
        self._save_json_summary(
            sample_id, similarity_results, stage_assessments, de_summary
        )

        print(f"  Report saved to: {report_file}")
        return str(report_file)

    def _generate_header(self, sample_id: str, report_time: str) -> str:
        """Generate report header"""
        return f"""# Disease Analysis Report

**Sample ID:** {sample_id}
**Analysis Date:** {report_time}
**Report Version:** 1.0

---
"""

    def _generate_executive_summary(self,
                                    similarity_results: List,
                                    stage_assessments: Dict,
                                    de_summary: Dict) -> str:
        """Generate executive summary"""

        # Get top disease matches
        top_matches = []
        for r in similarity_results[:3]:
            name = r.signature_name if hasattr(r, 'signature_name') else r.get('signature_name', 'Unknown')
            score = r.combined_score if hasattr(r, 'combined_score') else r.get('combined_score', 0)
            top_matches.append(f"- {name} (score: {score:.3f})")

        # Get primary stage assessment
        primary_cat = None
        primary_assessment = None
        best_confidence = 0

        for cat, assessment in stage_assessments.items():
            conf = assessment.confidence if hasattr(assessment, 'confidence') else assessment.get('confidence', 0)
            if conf > best_confidence:
                best_confidence = conf
                primary_cat = cat
                primary_assessment = assessment

        stage = primary_assessment.predicted_stage if hasattr(primary_assessment, 'predicted_stage') else primary_assessment.get('predicted_stage', 'unknown')
        severity = primary_assessment.severity_score if hasattr(primary_assessment, 'severity_score') else primary_assessment.get('severity_score', 0)

        # DE genes summary
        n_sig = de_summary.get('n_significant', 0)
        n_up = de_summary.get('n_upregulated', 0)
        n_down = de_summary.get('n_downregulated', 0)

        return f"""## Executive Summary

### Key Findings

| Metric | Value |
|--------|-------|
| **Primary Disease Category** | {primary_cat.upper() if primary_cat else 'Unknown'} |
| **Predicted Stage** | {stage.upper()} |
| **Severity Score** | {severity:.2f} / 1.00 |
| **Confidence** | {best_confidence:.2f} |
| **Significant DE Genes** | {n_sig} ({n_up} up, {n_down} down) |

### Top Disease Signature Matches
{chr(10).join(top_matches)}

### Risk Assessment
{"**HIGH RISK**: Elevated severity markers detected" if severity > 0.7 else "**MODERATE RISK**: Some disease markers present" if severity > 0.4 else "**LOW RISK**: Minimal disease indicators"}
"""

    def _generate_similarity_section(self, similarity_results: List) -> str:
        """Generate disease similarity section"""

        rows = []
        for i, r in enumerate(similarity_results[:15], 1):
            name = r.signature_name if hasattr(r, 'signature_name') else r.get('signature_name', 'Unknown')
            cat = r.category if hasattr(r, 'category') else r.get('category', 'unknown')
            stage = r.stage if hasattr(r, 'stage') else r.get('stage', 'unknown')
            score = r.combined_score if hasattr(r, 'combined_score') else r.get('combined_score', 0)
            overlap = r.gene_overlap if hasattr(r, 'gene_overlap') else r.get('gene_overlap', 0)

            rows.append(f"| {i} | {name[:50]} | {cat} | {stage} | {score:.3f} | {overlap} |")

        table = "\n".join(rows)

        return f"""## Disease Similarity Analysis

### Top Matching Disease Signatures

| Rank | Signature | Category | Stage | Score | Gene Overlap |
|------|-----------|----------|-------|-------|--------------|
{table}

### Interpretation
The similarity analysis compares your gene expression profile against {len(similarity_results)} known disease signatures.
Higher scores indicate stronger similarity to the disease state.
"""

    def _generate_stage_section(self, stage_assessments: Dict) -> str:
        """Generate disease stage assessment section"""

        rows = []
        for cat, assessment in stage_assessments.items():
            stage = assessment.predicted_stage if hasattr(assessment, 'predicted_stage') else assessment.get('predicted_stage', 'unknown')
            conf = assessment.confidence if hasattr(assessment, 'confidence') else assessment.get('confidence', 0)
            severity = assessment.severity_score if hasattr(assessment, 'severity_score') else assessment.get('severity_score', 0)

            rows.append(f"| {cat.upper()} | {stage} | {conf:.2f} | {severity:.2f} |")

        table = "\n".join(rows)

        # Stage indicators detail
        detail_sections = []
        for cat, assessment in stage_assessments.items():
            indicators = assessment.stage_indicators if hasattr(assessment, 'stage_indicators') else assessment.get('stage_indicators', {})
            if indicators:
                indicator_str = ", ".join([f"{k}: {v:.2f}" for k, v in sorted(indicators.items(), key=lambda x: -x[1])])
                detail_sections.append(f"**{cat.upper()}**: {indicator_str}")

        return f"""## Disease Stage Assessment

### Stage Predictions by Category

| Category | Predicted Stage | Confidence | Severity |
|----------|-----------------|------------|----------|
{table}

### Stage Indicator Scores
{chr(10).join(detail_sections)}

### Methodology
Stage prediction is based on:
- Known stage-specific biomarker expression
- Gene set enrichment patterns
- Pathway activity analysis
"""

    def _generate_biomarker_section(self,
                                    stage_assessments: Dict,
                                    de_summary: Dict) -> str:
        """Generate biomarker analysis section"""

        # Collect all biomarkers
        present_markers = set()
        absent_markers = set()

        for cat, assessment in stage_assessments.items():
            markers = assessment.biomarkers_present if hasattr(assessment, 'biomarkers_present') else assessment.get('biomarkers_present', [])
            absent = assessment.biomarkers_absent if hasattr(assessment, 'biomarkers_absent') else assessment.get('biomarkers_absent', [])
            present_markers.update(markers)
            absent_markers.update(absent)

        # Top DE genes
        top_genes = de_summary.get('top_genes', [])[:10]
        top_genes_str = ", ".join(top_genes) if top_genes else "N/A"

        return f"""## Biomarker Analysis

### Disease Biomarkers Detected
{', '.join(sorted(present_markers)) if present_markers else 'None detected'}

### Expected Biomarkers Not Detected
{', '.join(sorted(absent_markers)) if absent_markers else 'All expected markers present'}

### Top Differentially Expressed Genes
{top_genes_str}

### Biomarker Significance
The presence or absence of specific biomarkers helps determine disease state and progression.
"""

    def _generate_pathway_section(self, pathway_results: Dict) -> str:
        """Generate pathway analysis section"""

        sections = []

        for pathway_type, results in pathway_results.items():
            if isinstance(results, list) and len(results) > 0:
                top_pathways = results[:5]
                pathway_list = "\n".join([f"- {p}" for p in top_pathways])
                sections.append(f"**{pathway_type}:**\n{pathway_list}")

        return f"""## Pathway Analysis

### Enriched Pathways
{chr(10).join(sections) if sections else 'No significant pathway enrichment detected'}
"""

    def _generate_clinical_section(self,
                                   similarity_results: List,
                                   stage_assessments: Dict) -> str:
        """Generate clinical interpretation section"""

        # Analyze patterns
        categories = {}
        for r in similarity_results[:20]:
            cat = r.category if hasattr(r, 'category') else r.get('category', 'unknown')
            if cat not in categories:
                categories[cat] = 0
            categories[cat] += 1

        dominant_category = max(categories, key=categories.get) if categories else 'unknown'

        # Get evidence
        evidence_items = []
        for cat, assessment in stage_assessments.items():
            evidence = assessment.supporting_evidence if hasattr(assessment, 'supporting_evidence') else assessment.get('supporting_evidence', [])
            for e in evidence[:3]:
                evidence_items.append(f"- {e}")

        return f"""## Clinical Interpretation

### Disease Pattern Analysis
The gene expression profile shows strongest similarity to **{dominant_category.upper()}** disease signatures.

Distribution of top matches:
{chr(10).join([f'- {cat}: {count} signatures' for cat, count in sorted(categories.items(), key=lambda x: -x[1])])}

### Supporting Evidence
{chr(10).join(evidence_items[:10]) if evidence_items else '- No specific evidence markers detected'}

### Limitations
- Analysis based on gene expression patterns only
- Clinical correlation recommended
- Results should be interpreted by qualified professionals
"""

    def _generate_recommendations(self, stage_assessments: Dict) -> str:
        """Generate recommendations section"""

        all_recommendations = []

        for cat, assessment in stage_assessments.items():
            recs = assessment.recommendations if hasattr(assessment, 'recommendations') else assessment.get('recommendations', [])
            for rec in recs:
                if rec not in all_recommendations:
                    all_recommendations.append(f"- [{cat.upper()}] {rec}")

        return f"""## Recommendations

### Clinical Recommendations
{chr(10).join(all_recommendations[:10]) if all_recommendations else '- Continue routine monitoring'}

### Follow-up Suggestions
1. Validate findings with clinical assessment
2. Consider additional diagnostic tests for confirmed markers
3. Monitor disease progression with regular expression profiling
4. Consult specialist for treatment planning if indicated

---
*This report is generated automatically and should be reviewed by qualified medical professionals.*
*Report generated using RNA-seq Disease Analysis System v1.0*
"""

    def _generate_visualizations(self,
                                 sample_id: str,
                                 similarity_results: List,
                                 stage_assessments: Dict):
        """Generate visualization plots"""

        if not PLOTTING_AVAILABLE:
            return

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        # 1. Top similarity scores
        ax1 = axes[0, 0]
        names = [r.signature_name[:30] if hasattr(r, 'signature_name') else 'Unknown'
                for r in similarity_results[:10]]
        scores = [r.combined_score if hasattr(r, 'combined_score') else 0
                 for r in similarity_results[:10]]

        colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(names)))
        ax1.barh(range(len(names)), scores, color=colors)
        ax1.set_yticks(range(len(names)))
        ax1.set_yticklabels(names)
        ax1.set_xlabel('Similarity Score')
        ax1.set_title('Top Disease Signature Matches')
        ax1.invert_yaxis()

        # 2. Stage assessment by category
        ax2 = axes[0, 1]
        cats = list(stage_assessments.keys())
        severities = [a.severity_score if hasattr(a, 'severity_score') else a.get('severity_score', 0)
                     for a in stage_assessments.values()]
        confidences = [a.confidence if hasattr(a, 'confidence') else a.get('confidence', 0)
                      for a in stage_assessments.values()]

        x = np.arange(len(cats))
        width = 0.35
        ax2.bar(x - width/2, severities, width, label='Severity', color='coral')
        ax2.bar(x + width/2, confidences, width, label='Confidence', color='steelblue')
        ax2.set_xticks(x)
        ax2.set_xticklabels([c.upper() for c in cats])
        ax2.set_ylabel('Score')
        ax2.set_title('Disease Assessment by Category')
        ax2.legend()
        ax2.set_ylim(0, 1)

        # 3. Category distribution
        ax3 = axes[1, 0]
        category_counts = {}
        for r in similarity_results[:30]:
            cat = r.category if hasattr(r, 'category') else r.get('category', 'unknown')
            category_counts[cat] = category_counts.get(cat, 0) + 1

        if category_counts:
            ax3.pie(category_counts.values(), labels=category_counts.keys(),
                   autopct='%1.1f%%', colors=plt.cm.Set3.colors)
            ax3.set_title('Disease Category Distribution')

        # 4. Gene overlap heatmap
        ax4 = axes[1, 1]
        overlaps = [r.gene_overlap if hasattr(r, 'gene_overlap') else 0
                   for r in similarity_results[:15]]
        emb_sims = [r.embedding_similarity if hasattr(r, 'embedding_similarity') else 0
                   for r in similarity_results[:15]]

        ax4.scatter(emb_sims, overlaps, c=range(len(overlaps)), cmap='viridis', s=100)
        ax4.set_xlabel('Embedding Similarity')
        ax4.set_ylabel('Gene Overlap')
        ax4.set_title('Similarity vs Gene Overlap')

        plt.tight_layout()
        plt.savefig(self.output_dir / f"disease_analysis_{sample_id}.png", dpi=150)
        plt.close()

    def _save_json_summary(self,
                          sample_id: str,
                          similarity_results: List,
                          stage_assessments: Dict,
                          de_summary: Dict):
        """Save JSON summary for programmatic access"""

        summary = {
            'sample_id': sample_id,
            'analysis_time': datetime.now().isoformat(),
            'top_matches': [],
            'stage_assessments': {},
            'de_summary': de_summary
        }

        for r in similarity_results[:20]:
            summary['top_matches'].append({
                'signature_id': r.signature_id if hasattr(r, 'signature_id') else '',
                'signature_name': r.signature_name if hasattr(r, 'signature_name') else '',
                'category': r.category if hasattr(r, 'category') else '',
                'stage': r.stage if hasattr(r, 'stage') else '',
                'combined_score': r.combined_score if hasattr(r, 'combined_score') else 0,
                'gene_overlap': r.gene_overlap if hasattr(r, 'gene_overlap') else 0
            })

        for cat, assessment in stage_assessments.items():
            summary['stage_assessments'][cat] = {
                'predicted_stage': assessment.predicted_stage if hasattr(assessment, 'predicted_stage') else '',
                'confidence': assessment.confidence if hasattr(assessment, 'confidence') else 0,
                'severity_score': assessment.severity_score if hasattr(assessment, 'severity_score') else 0,
                'biomarkers_present': assessment.biomarkers_present if hasattr(assessment, 'biomarkers_present') else [],
                'recommendations': assessment.recommendations if hasattr(assessment, 'recommendations') else []
            }

        json_file = self.output_dir / f"disease_analysis_{sample_id}.json"
        with open(json_file, 'w') as f:
            json.dump(summary, f, indent=2)


if __name__ == "__main__":
    generator = DiseaseAnalysisReportGenerator()
    print("Report generator initialized")
