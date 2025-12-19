#!/usr/bin/env python3
"""
Disease Database Builder
========================

Builds comprehensive disease signature databases including:
1. MSigDB disease/perturbation signatures
2. Cancer staging signatures (from literature)
3. Disease severity markers
4. Drug response signatures

Includes metadata about:
- Disease categories (cancer, neurological, metabolic, etc.)
- Disease stages (early, intermediate, advanced, etc.)
- Severity scores
- Associated pathways
"""

import os
import json
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict
from enum import Enum
import warnings
warnings.filterwarnings('ignore')


class DiseaseCategory(Enum):
    """Disease categories"""
    CANCER = "cancer"
    NEUROLOGICAL = "neurological"
    CARDIOVASCULAR = "cardiovascular"
    METABOLIC = "metabolic"
    INFLAMMATORY = "inflammatory"
    INFECTIOUS = "infectious"
    AUTOIMMUNE = "autoimmune"
    GENETIC = "genetic"
    AGING = "aging"
    OTHER = "other"


class DiseaseStage(Enum):
    """Disease progression stages"""
    PRECLINICAL = "preclinical"
    EARLY = "early"
    INTERMEDIATE = "intermediate"
    ADVANCED = "advanced"
    TERMINAL = "terminal"
    REMISSION = "remission"
    UNKNOWN = "unknown"


@dataclass
class DiseaseSignature:
    """Disease signature with metadata"""
    id: str
    name: str
    genes: List[str]
    category: str
    stage: str
    severity_score: float  # 0-1 scale
    description: str
    source: str
    associated_pathways: List[str]
    up_genes: List[str]
    down_genes: List[str]
    metadata: Dict


class DiseaseDatabaseBuilder:
    """Build comprehensive disease signature database"""

    def __init__(self, data_dir: str = "data/disease_db"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.signatures: Dict[str, DiseaseSignature] = {}

    def build_database(self) -> Dict[str, DiseaseSignature]:
        """Build complete disease database"""
        print("\n" + "="*60)
        print("Building Comprehensive Disease Signature Database")
        print("="*60)

        # 1. Load MSigDB signatures with classification
        self._load_msigdb_signatures()

        # 2. Add cancer staging signatures
        self._add_cancer_staging_signatures()

        # 3. Add neurological disease signatures
        self._add_neurological_signatures()

        # 4. Add metabolic disease signatures
        self._add_metabolic_signatures()

        # 5. Add aging signatures
        self._add_aging_signatures()

        # 6. Save database
        self._save_database()

        print(f"\n  Total signatures: {len(self.signatures)}")
        return self.signatures

    def _load_msigdb_signatures(self):
        """Load and classify MSigDB signatures"""
        print("\n[1/6] Loading MSigDB signatures...")

        # Load existing MSigDB data
        msigdb_file = self.data_dir.parent / "msigdb" / "disease_gene_sets.json"

        if not msigdb_file.exists():
            print("  MSigDB data not found. Fetching...")
            self._fetch_msigdb_data()

        with open(msigdb_file) as f:
            raw_data = json.load(f)

        # Classify signatures
        cancer_keywords = ['CANCER', 'TUMOR', 'CARCINOMA', 'LEUKEMIA', 'LYMPHOMA',
                          'MELANOMA', 'SARCOMA', 'MYELOMA', 'GLIOMA', 'ONCOGENE']
        neuro_keywords = ['ALZHEIMER', 'PARKINSON', 'HUNTINGTON', 'NEURO', 'BRAIN',
                         'SCHIZOPHRENIA', 'AUTISM', 'ALS', 'DEMENTIA']
        metabolic_keywords = ['DIABETES', 'OBESITY', 'METABOLIC', 'INSULIN', 'LIPID',
                             'CHOLESTEROL', 'GLUCOSE']
        inflammatory_keywords = ['INFLAM', 'ARTHRITIS', 'COLITIS', 'ASTHMA', 'ALLERGY']
        aging_keywords = ['AGING', 'SENESCENCE', 'LONGEVITY', 'ELDERLY', 'WERNER']

        for gs_name, gs_data in raw_data.items():
            name_upper = gs_name.upper()

            # Classify category
            if any(kw in name_upper for kw in cancer_keywords):
                category = DiseaseCategory.CANCER.value
            elif any(kw in name_upper for kw in neuro_keywords):
                category = DiseaseCategory.NEUROLOGICAL.value
            elif any(kw in name_upper for kw in metabolic_keywords):
                category = DiseaseCategory.METABOLIC.value
            elif any(kw in name_upper for kw in inflammatory_keywords):
                category = DiseaseCategory.INFLAMMATORY.value
            elif any(kw in name_upper for kw in aging_keywords):
                category = DiseaseCategory.AGING.value
            else:
                category = DiseaseCategory.OTHER.value

            # Determine stage from name
            stage = DiseaseStage.UNKNOWN.value
            if 'EARLY' in name_upper or 'STAGE_I' in name_upper:
                stage = DiseaseStage.EARLY.value
            elif 'ADVANCED' in name_upper or 'STAGE_IV' in name_upper or 'METASTA' in name_upper:
                stage = DiseaseStage.ADVANCED.value
            elif 'INTERMEDIATE' in name_upper or 'STAGE_II' in name_upper:
                stage = DiseaseStage.INTERMEDIATE.value

            # Determine up/down regulation
            genes = gs_data.get('genes', [])
            if isinstance(genes[0], list) if genes else False:
                genes = [g[0] if isinstance(g, list) else g for g in genes]

            up_genes = genes if '_UP' in name_upper else []
            down_genes = genes if '_DN' in name_upper or '_DOWN' in name_upper else []

            name = gs_data.get('name', gs_name)
            if isinstance(name, list):
                name = name[0] if name else gs_name

            sig = DiseaseSignature(
                id=gs_name,
                name=str(name),
                genes=genes,
                category=category,
                stage=stage,
                severity_score=0.5,  # Default
                description=f"MSigDB signature: {name}",
                source="MSigDB",
                associated_pathways=[],
                up_genes=up_genes,
                down_genes=down_genes,
                metadata={'original_category': gs_data.get('category', 'unknown')}
            )
            self.signatures[gs_name] = sig

        print(f"    Loaded {len(self.signatures)} MSigDB signatures")

    def _add_cancer_staging_signatures(self):
        """Add cancer staging gene signatures"""
        print("\n[2/6] Adding cancer staging signatures...")

        # Well-known cancer staging markers
        cancer_stages = {
            # Breast cancer staging markers
            'BREAST_CANCER_STAGE_I': {
                'genes': ['ESR1', 'PGR', 'ERBB2', 'MKI67', 'CCND1', 'BCL2'],
                'stage': DiseaseStage.EARLY.value,
                'severity': 0.2,
                'description': 'Early breast cancer markers - hormone receptor positive'
            },
            'BREAST_CANCER_STAGE_II': {
                'genes': ['ESR1', 'PGR', 'ERBB2', 'MKI67', 'CCND1', 'AURKA', 'TOP2A'],
                'stage': DiseaseStage.INTERMEDIATE.value,
                'severity': 0.4
            },
            'BREAST_CANCER_STAGE_III': {
                'genes': ['ERBB2', 'MKI67', 'AURKA', 'TOP2A', 'BIRC5', 'MMP9', 'VEGFA'],
                'stage': DiseaseStage.INTERMEDIATE.value,
                'severity': 0.6
            },
            'BREAST_CANCER_STAGE_IV_METASTATIC': {
                'genes': ['MMP2', 'MMP9', 'VEGFA', 'TWIST1', 'SNAI1', 'CDH2', 'VIM', 'FN1'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.9,
                'description': 'Metastatic breast cancer - EMT markers'
            },

            # Lung cancer staging
            'LUNG_CANCER_EARLY': {
                'genes': ['EGFR', 'KRAS', 'TP53', 'STK11', 'KEAP1', 'NKX2-1'],
                'stage': DiseaseStage.EARLY.value,
                'severity': 0.3
            },
            'LUNG_CANCER_ADVANCED': {
                'genes': ['EGFR', 'ALK', 'ROS1', 'MET', 'RET', 'BRAF', 'MYC', 'SOX2'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.8
            },

            # Colorectal cancer staging
            'COLORECTAL_CANCER_STAGE_I': {
                'genes': ['APC', 'KRAS', 'CTNNB1', 'MYC', 'CCND1'],
                'stage': DiseaseStage.EARLY.value,
                'severity': 0.2
            },
            'COLORECTAL_CANCER_STAGE_IV': {
                'genes': ['APC', 'KRAS', 'TP53', 'SMAD4', 'PIK3CA', 'BRAF', 'MMP7'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.85
            },

            # General cancer progression markers
            'CANCER_EMT_SIGNATURE': {
                'genes': ['CDH1', 'CDH2', 'VIM', 'SNAI1', 'SNAI2', 'TWIST1', 'ZEB1', 'ZEB2'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.7,
                'description': 'Epithelial-mesenchymal transition markers'
            },
            'CANCER_STEMNESS_SIGNATURE': {
                'genes': ['CD44', 'ALDH1A1', 'NANOG', 'POU5F1', 'SOX2', 'PROM1', 'LGR5'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.8,
                'description': 'Cancer stem cell markers'
            },
            'CANCER_ANGIOGENESIS_SIGNATURE': {
                'genes': ['VEGFA', 'VEGFB', 'VEGFC', 'FLT1', 'KDR', 'ANGPT1', 'ANGPT2'],
                'stage': DiseaseStage.INTERMEDIATE.value,
                'severity': 0.6
            }
        }

        for sig_id, data in cancer_stages.items():
            sig = DiseaseSignature(
                id=sig_id,
                name=sig_id.replace('_', ' ').title(),
                genes=data['genes'],
                category=DiseaseCategory.CANCER.value,
                stage=data['stage'],
                severity_score=data['severity'],
                description=data.get('description', f'Cancer staging signature: {sig_id}'),
                source="Curated",
                associated_pathways=['Cell Cycle', 'Apoptosis', 'Metastasis'],
                up_genes=data['genes'],
                down_genes=[],
                metadata={'curated': True, 'type': 'staging'}
            )
            self.signatures[sig_id] = sig

        print(f"    Added {len(cancer_stages)} cancer staging signatures")

    def _add_neurological_signatures(self):
        """Add neurological disease signatures"""
        print("\n[3/6] Adding neurological disease signatures...")

        neuro_signatures = {
            'ALZHEIMERS_EARLY': {
                'genes': ['APP', 'PSEN1', 'PSEN2', 'APOE', 'MAPT', 'BACE1', 'CLU', 'BIN1', 'PICALM', 'CR1', 'SORL1'],
                'stage': DiseaseStage.EARLY.value,
                'severity': 0.3,
                'description': 'Early Alzheimer\'s disease markers'
            },
            'ALZHEIMERS_ADVANCED': {
                'genes': ['APP', 'MAPT', 'APOE', 'TREM2', 'GFAP', 'AQP4', 'S100B', 'IL1B', 'IL6', 'AIF1', 'CD68', 'ITGAM', 'C3', 'C1QA'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.85,
                'description': 'Advanced Alzheimer\'s with neuroinflammation'
            },
            'ALZHEIMERS_SYNAPTIC_LOSS': {
                'genes': ['SYP', 'DLG4', 'GRIN1', 'GRIN2A', 'GRIN2B', 'SLC17A7', 'SLC17A6', 'GAD1', 'GAD2', 'SNAP25', 'SYT1'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.75,
                'description': 'Synaptic loss in Alzheimer\'s disease',
                'down_genes': ['SYP', 'DLG4', 'GRIN1', 'GRIN2A', 'GRIN2B', 'SLC17A7', 'GAD1', 'GAD2', 'SNAP25', 'SYT1']
            },
            'ALZHEIMERS_NEURONAL_ACTIVITY': {
                'genes': ['NPAS4', 'FOSB', 'FOS', 'EGR1', 'EGR2', 'EGR4', 'ARC', 'BDNF', 'NR4A1', 'NR4A2', 'NR4A3', 'NPTX2'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.7,
                'description': 'Neuronal activity genes downregulated in Alzheimer\'s',
                'down_genes': ['NPAS4', 'FOSB', 'FOS', 'EGR1', 'EGR2', 'EGR4', 'ARC', 'BDNF', 'NPTX2']
            },
            'ALZHEIMERS_BLOOD_SIGNATURE': {
                'genes': ['TOMM40', 'APOE', 'CD33', 'TREM2', 'CLU', 'ABCA7', 'SORL1', 'BIN1', 'PICALM', 'INPP5D', 'MEF2C'],
                'stage': DiseaseStage.UNKNOWN.value,
                'severity': 0.6,
                'description': 'Blood-based Alzheimer\'s biomarkers from GWAS'
            },
            'NEUROINFLAMMATION_SIGNATURE': {
                'genes': ['GFAP', 'AIF1', 'CD68', 'ITGAM', 'TREM2', 'TYROBP', 'C1QA', 'C1QB', 'C3', 'IL1B', 'IL6', 'TNF', 'CCL2'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.7,
                'description': 'Microglial and astrocyte activation markers'
            },
            'PARKINSONS_EARLY': {
                'genes': ['SNCA', 'LRRK2', 'PARK7', 'PINK1', 'PRKN'],
                'stage': DiseaseStage.EARLY.value,
                'severity': 0.3
            },
            'PARKINSONS_ADVANCED': {
                'genes': ['SNCA', 'LRRK2', 'TH', 'DDC', 'SLC6A3', 'GFAP', 'AIF1'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.8
            },
            'ALS_SIGNATURE': {
                'genes': ['SOD1', 'TARDBP', 'FUS', 'C9orf72', 'OPTN', 'VCP'],
                'stage': DiseaseStage.UNKNOWN.value,
                'severity': 0.7
            },
            'HUNTINGTONS_SIGNATURE': {
                'genes': ['HTT', 'BDNF', 'REST', 'PGC1A', 'DARPP32'],
                'stage': DiseaseStage.UNKNOWN.value,
                'severity': 0.6
            }
        }

        for sig_id, data in neuro_signatures.items():
            sig = DiseaseSignature(
                id=sig_id,
                name=sig_id.replace('_', ' ').title(),
                genes=data['genes'],
                category=DiseaseCategory.NEUROLOGICAL.value,
                stage=data['stage'],
                severity_score=data['severity'],
                description=data.get('description', ''),
                source="Curated",
                associated_pathways=['Neurodegeneration', 'Protein Aggregation'],
                up_genes=data['genes'],
                down_genes=[],
                metadata={'curated': True}
            )
            self.signatures[sig_id] = sig

        print(f"    Added {len(neuro_signatures)} neurological signatures")

    def _add_metabolic_signatures(self):
        """Add metabolic disease signatures"""
        print("\n[4/6] Adding metabolic disease signatures...")

        metabolic_signatures = {
            'TYPE2_DIABETES_PREDIABETES': {
                'genes': ['INS', 'INSR', 'IRS1', 'IRS2', 'SLC2A4', 'PPARG', 'ADIPOQ'],
                'stage': DiseaseStage.PRECLINICAL.value,
                'severity': 0.2
            },
            'TYPE2_DIABETES_EARLY': {
                'genes': ['INS', 'INSR', 'IRS1', 'SLC2A4', 'PPARG', 'TCF7L2', 'KCNJ11'],
                'stage': DiseaseStage.EARLY.value,
                'severity': 0.4
            },
            'TYPE2_DIABETES_ADVANCED': {
                'genes': ['INS', 'GCK', 'HNF1A', 'HNF4A', 'PDX1', 'NEUROD1', 'FOXO1'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.7
            },
            'NAFLD_EARLY': {
                'genes': ['PPARA', 'PPARG', 'SREBF1', 'FASN', 'ACACA', 'SCD'],
                'stage': DiseaseStage.EARLY.value,
                'severity': 0.3,
                'description': 'Non-alcoholic fatty liver disease - early stage'
            },
            'NASH_ADVANCED': {
                'genes': ['TGFB1', 'COL1A1', 'ACTA2', 'TIMP1', 'MMP2', 'IL6', 'TNF'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.75,
                'description': 'Non-alcoholic steatohepatitis with fibrosis'
            },
            'OBESITY_SIGNATURE': {
                'genes': ['LEP', 'LEPR', 'ADIPOQ', 'PPARG', 'UCP1', 'ADRB3', 'FTO'],
                'stage': DiseaseStage.UNKNOWN.value,
                'severity': 0.5
            }
        }

        for sig_id, data in metabolic_signatures.items():
            sig = DiseaseSignature(
                id=sig_id,
                name=sig_id.replace('_', ' ').title(),
                genes=data['genes'],
                category=DiseaseCategory.METABOLIC.value,
                stage=data['stage'],
                severity_score=data['severity'],
                description=data.get('description', ''),
                source="Curated",
                associated_pathways=['Metabolism', 'Insulin Signaling'],
                up_genes=data['genes'],
                down_genes=[],
                metadata={'curated': True}
            )
            self.signatures[sig_id] = sig

        print(f"    Added {len(metabolic_signatures)} metabolic signatures")

    def _add_aging_signatures(self):
        """Add aging and senescence signatures"""
        print("\n[5/6] Adding aging signatures...")

        aging_signatures = {
            'CELLULAR_SENESCENCE': {
                'genes': ['CDKN2A', 'CDKN1A', 'TP53', 'RB1', 'SERPINE1', 'IL6', 'IL8', 'MMP3'],
                'stage': DiseaseStage.UNKNOWN.value,
                'severity': 0.5,
                'description': 'Cellular senescence markers (SASP)'
            },
            'BIOLOGICAL_AGING_EARLY': {
                'genes': ['TERT', 'TERC', 'SIRT1', 'SIRT3', 'FOXO3', 'KLOTHO'],
                'stage': DiseaseStage.EARLY.value,
                'severity': 0.3
            },
            'BIOLOGICAL_AGING_ADVANCED': {
                'genes': ['CDKN2A', 'TP53', 'LMNA', 'WRN', 'BLM', 'ERCC1'],
                'stage': DiseaseStage.ADVANCED.value,
                'severity': 0.7
            },
            'INFLAMMAGING': {
                'genes': ['IL6', 'IL1B', 'TNF', 'CRP', 'IL18', 'CXCL8', 'CCL2'],
                'stage': DiseaseStage.UNKNOWN.value,
                'severity': 0.6,
                'description': 'Age-related chronic inflammation'
            },
            'MITOCHONDRIAL_DYSFUNCTION_AGING': {
                'genes': ['MT-ND1', 'MT-CO1', 'MT-ATP6', 'PPARGC1A', 'TFAM', 'NRF1'],
                'stage': DiseaseStage.UNKNOWN.value,
                'severity': 0.5
            }
        }

        for sig_id, data in aging_signatures.items():
            sig = DiseaseSignature(
                id=sig_id,
                name=sig_id.replace('_', ' ').title(),
                genes=data['genes'],
                category=DiseaseCategory.AGING.value,
                stage=data['stage'],
                severity_score=data['severity'],
                description=data.get('description', ''),
                source="Curated",
                associated_pathways=['Senescence', 'DNA Repair', 'Inflammation'],
                up_genes=data['genes'],
                down_genes=[],
                metadata={'curated': True}
            )
            self.signatures[sig_id] = sig

        print(f"    Added {len(aging_signatures)} aging signatures")

    def _fetch_msigdb_data(self):
        """Fetch MSigDB data using R"""
        print("  Fetching MSigDB data via R...")

        msigdb_dir = self.data_dir.parent / "msigdb"
        msigdb_dir.mkdir(parents=True, exist_ok=True)

        r_script = '''
        library(msigdbr)
        library(jsonlite)

        h_sets <- msigdbr(species = "Homo sapiens", collection = "H")
        cgp_sets <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")

        all_data <- rbind(
            transform(h_sets, category = "hallmark"),
            transform(cgp_sets, category = "perturbation")
        )

        gene_sets <- list()
        for (gs_name in unique(all_data$gs_name)) {
            gs_data <- all_data[all_data$gs_name == gs_name, ]
            genes <- unique(gs_data$gene_symbol)
            gene_sets[[gs_name]] <- list(
                name = gs_name,
                genes = as.list(genes),
                category = unique(gs_data$category)[1],
                n_genes = length(genes)
            )
        }

        write(toJSON(gene_sets, auto_unbox = FALSE), "''' + str(msigdb_dir / "disease_gene_sets.json") + '''")
        '''

        subprocess.run(['Rscript', '-e', r_script], capture_output=True)

    def _save_database(self):
        """Save database to file"""
        print("\n[6/6] Saving database...")

        # Convert to serializable format
        db_data = {}
        for sig_id, sig in self.signatures.items():
            db_data[sig_id] = {
                'id': sig.id,
                'name': sig.name,
                'genes': sig.genes,
                'category': sig.category,
                'stage': sig.stage,
                'severity_score': sig.severity_score,
                'description': sig.description,
                'source': sig.source,
                'associated_pathways': sig.associated_pathways,
                'up_genes': sig.up_genes,
                'down_genes': sig.down_genes,
                'metadata': sig.metadata
            }

        # Save main database
        db_file = self.data_dir / "disease_signatures_full.json"
        with open(db_file, 'w') as f:
            json.dump(db_data, f, indent=2)
        print(f"    Saved to {db_file}")

        # Save category index
        categories = {}
        for sig_id, sig in self.signatures.items():
            cat = sig.category
            if cat not in categories:
                categories[cat] = []
            categories[cat].append(sig_id)

        index_file = self.data_dir / "category_index.json"
        with open(index_file, 'w') as f:
            json.dump(categories, f, indent=2)

        # Save statistics
        stats = {
            'total_signatures': len(self.signatures),
            'by_category': {k: len(v) for k, v in categories.items()},
            'by_stage': {},
            'curated_count': sum(1 for s in self.signatures.values() if s.metadata.get('curated'))
        }

        for sig in self.signatures.values():
            stage = sig.stage
            stats['by_stage'][stage] = stats['by_stage'].get(stage, 0) + 1

        stats_file = self.data_dir / "database_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)

        print(f"    Database statistics:")
        print(f"      Total signatures: {stats['total_signatures']}")
        print(f"      Curated: {stats['curated_count']}")
        for cat, count in stats['by_category'].items():
            print(f"      {cat}: {count}")

    def load_database(self) -> Dict[str, DiseaseSignature]:
        """Load existing database"""
        db_file = self.data_dir / "disease_signatures_full.json"

        if not db_file.exists():
            print("  Database not found. Building...")
            return self.build_database()

        with open(db_file) as f:
            db_data = json.load(f)

        for sig_id, data in db_data.items():
            self.signatures[sig_id] = DiseaseSignature(**data)

        print(f"  Loaded {len(self.signatures)} signatures from database")
        return self.signatures


if __name__ == "__main__":
    builder = DiseaseDatabaseBuilder()
    builder.build_database()
