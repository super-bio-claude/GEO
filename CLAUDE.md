# RNA-seq Disease Analysis Platform

## Project Overview

This is a comprehensive RNA-seq analysis platform that performs:
1. **Data Preprocessing** - Normalization, batch correction
2. **Differential Expression Analysis** - DESeq2 with apeglm shrinkage
3. **Pathway Analysis** - KEGG, GO, GSVA
4. **Disease Similarity Analysis** - ChromaDB vector search against MSigDB signatures
5. **Disease Stage Prediction** - Multi-category assessment (Cancer, Neurological, Metabolic, Aging)
6. **Clinical Report Generation** - Comprehensive markdown reports with visualizations

## Tech Stack

### Backend
- **Python 3.x** - Main analysis pipeline
- **R/Bioconductor** - DESeq2, pathway analysis, GSVA
- **ChromaDB** - Vector database for disease signature similarity
- **Sentence Transformers** - Gene set embeddings (all-MiniLM-L6-v2)

### Key Libraries
- pandas, numpy, scikit-learn
- DESeq2, apeglm, clusterProfiler, GSVA, msigdbr
- chromadb, sentence-transformers

## Directory Structure

```
/Users/admin/GEO/
├── src/
│   ├── disease_analysis/          # Core disease analysis module
│   │   ├── pipeline.py            # Main entry point
│   │   ├── disease_database.py    # 3,620 disease signatures
│   │   ├── similarity_engine.py   # ChromaDB vector search
│   │   ├── stage_predictor.py     # Disease stage prediction
│   │   └── report_generator.py    # Report generation
│   ├── preprocessing/
│   │   ├── data_loader.py
│   │   └── normalization.py
│   ├── de_analysis/
│   │   └── differential_expression.py
│   ├── batch_correction/
│   │   └── combat.py
│   ├── R_scripts/
│   │   ├── deseq2_analysis.R
│   │   ├── pathway_analysis.R
│   │   └── gsva_analysis.R
│   ├── ml_models/
│   │   └── classifiers.py
│   └── r_integration/
│       └── r_runner.py
├── data/
│   ├── disease_db/                # Disease signature database
│   ├── chromadb_disease/          # ChromaDB vector index
│   └── geo_test/                  # Test datasets from GEO
├── results/
│   └── disease_analysis/          # Analysis outputs
├── configs/
├── notebooks/
└── requirements.txt
```

## Core Modules

### 1. Disease Analysis Pipeline (`src/disease_analysis/pipeline.py`)

Main entry point for disease analysis.

```bash
# Usage
python src/disease_analysis/pipeline.py \
  --de_results <path_to_deseq2_results.csv> \
  --output <output_directory> \
  --sample_id <sample_name> \
  --padj 0.05 \
  --lfc 0.5
```

**Input Format (DE Results CSV):**
- Required columns: `gene`, `log2FoldChange`, `padj`
- Optional: `baseMean`, `lfcSE`, `pvalue`

### 2. Disease Database (`src/disease_analysis/disease_database.py`)

Contains 3,620 disease signatures:
- 3,588 from MSigDB (Hallmark, C2:CGP perturbation signatures)
- 32 curated signatures (Cancer staging, Neurological, Metabolic, Aging)

**Disease Categories:**
- `cancer` - Oncology signatures with staging
- `neurological` - Alzheimer's, Parkinson's, ALS, etc.
- `metabolic` - Diabetes, NAFLD, Obesity
- `aging` - Senescence, Inflammaging
- `inflammatory` - General inflammation
- `other` - Miscellaneous perturbation signatures

### 3. Similarity Engine (`src/disease_analysis/similarity_engine.py`)

ChromaDB-based vector similarity search.

**Scoring Components:**
- Embedding similarity (40%) - Semantic similarity via sentence-transformers
- Gene overlap / Jaccard index (30%)
- Direction score (30%) - Agreement in up/down regulation

### 4. Stage Predictor (`src/disease_analysis/stage_predictor.py`)

Predicts disease stage based on biomarker expression.

**Neurological Markers (Alzheimer's-optimized):**
- Preclinical: BDNF+, NPAS4+, ARC+, NTRK2+
- Early: APP+, MAPT+, BACE1+, CLU+
- Advanced: GFAP+, AIF1+, IL1B+, NPAS4-, FOSB-, EGR2-

**Output:**
```python
StageAssessment(
    disease_category: str,
    predicted_stage: str,        # preclinical, early, intermediate, advanced
    confidence: float,           # 0-1
    severity_score: float,       # 0-1
    supporting_evidence: List[str],
    biomarkers_present: List[str],
    recommendations: List[str]
)
```

## API Design (for Web Platform)

### Planned Endpoints

```
POST /api/analyze
  - Upload DE results CSV
  - Returns: Analysis report (JSON/Markdown)

GET /api/signatures
  - List available disease signatures
  - Filter by category, stage

GET /api/reports/{sample_id}
  - Get analysis report

POST /api/upload-counts
  - Upload raw counts matrix
  - Triggers DESeq2 + disease analysis
```

### Output Format

```json
{
  "sample_id": "GSE125583_Alzheimer",
  "primary_category": "NEUROLOGICAL",
  "predicted_stage": "ADVANCED",
  "severity_score": 0.70,
  "confidence": 0.47,
  "risk_assessment": "HIGH",
  "top_signatures": [...],
  "biomarkers": {
    "detected": ["GFAP", "FOSB_down", "EGR2_down"],
    "missing": ["BDNF", "TREM2"]
  },
  "recommendations": [...],
  "report_url": "/reports/GSE125583_Alzheimer.md"
}
```

## Development Notes

### Running DESeq2 Analysis

```r
# R script for DE analysis
library(DESeq2)
library(apeglm)

dds <- DESeqDataSetFromMatrix(counts, metadata, ~ condition)
dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "condition_treatment_vs_control", type = "apeglm")
```

### Building Disease Index

```python
from disease_analysis.disease_database import DiseaseDatabaseBuilder
from disease_analysis.similarity_engine import DiseaseSimilarityEngine

# Build database
builder = DiseaseDatabaseBuilder("data/disease_db")
signatures = builder.build_database()

# Build ChromaDB index
engine = DiseaseSimilarityEngine("data/chromadb_disease")
engine.build_index(signatures)
```

### Testing with GEO Data

```bash
# Download and process GEO dataset
Rscript src/geo_test/prepare_alzheimer_data.R

# Run DESeq2
Rscript src/geo_test/run_deseq2_alzheimer.R

# Run disease analysis
python src/disease_analysis/pipeline.py \
  --de_results results/geo_test/GSE125583/deseq2_results.csv \
  --sample_id GSE125583_Alzheimer
```

## Validated Test Cases

### GSE125583 (Alzheimer's Disease Blood)
- 289 samples (219 AD vs 70 controls)
- **Result:** NEUROLOGICAL, ADVANCED, Severity 0.70
- Correctly detected: FOSB/EGR2 downregulation, GFAP upregulation

### GSE313799 (Compression Study)
- Hypoxia signature detected
- BHLHE40 as top DE gene

## Web Platform Requirements (TODO)

### Frontend
- React/Next.js or Vue.js
- File upload component
- Interactive report viewer
- Visualization dashboards

### Backend API
- FastAPI or Flask
- Async job processing (Celery/Redis)
- File storage (S3 or local)
- Database (PostgreSQL for metadata, ChromaDB for signatures)

### Deployment
- Docker containers
- GPU support optional (for embeddings)
- Memory: 8GB+ recommended for ChromaDB

## Known Limitations

1. MSigDB signatures are mostly cancer-focused
2. Blood-based analysis may differ from tissue analysis
3. Confidence scores are relative, not absolute probabilities
4. Requires manual curation for disease-specific signatures

## Contributing

When adding new disease signatures:
1. Add to `disease_database.py` in appropriate category
2. Include stage-specific markers in `stage_predictor.py`
3. Rebuild index with `--rebuild_db` flag
4. Test with known disease datasets
