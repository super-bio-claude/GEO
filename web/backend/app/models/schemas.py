"""
Pydantic Schemas for API Request/Response
"""
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from datetime import datetime
from enum import Enum


class DiseaseCategory(str, Enum):
    CANCER = "cancer"
    NEUROLOGICAL = "neurological"
    METABOLIC = "metabolic"
    AGING = "aging"
    INFLAMMATORY = "inflammatory"
    OTHER = "other"


class DiseaseStage(str, Enum):
    PRECLINICAL = "preclinical"
    EARLY = "early"
    INTERMEDIATE = "intermediate"
    ADVANCED = "advanced"
    UNKNOWN = "unknown"


class RiskLevel(str, Enum):
    LOW = "LOW"
    MODERATE = "MODERATE"
    HIGH = "HIGH"


# Request Schemas
class AnalysisRequest(BaseModel):
    """Request for new analysis"""
    sample_id: Optional[str] = None
    padj_threshold: float = Field(default=0.05, ge=0, le=1)
    lfc_threshold: float = Field(default=0.5, ge=0)


class DEGeneInput(BaseModel):
    """Single gene from DE results"""
    gene: str
    log2FoldChange: float
    padj: Optional[float] = None
    baseMean: Optional[float] = None
    pvalue: Optional[float] = None


# Response Schemas
class SignatureMatch(BaseModel):
    """Disease signature match result"""
    rank: int
    signature_id: str
    signature_name: str
    category: DiseaseCategory
    stage: DiseaseStage
    score: float
    gene_overlap: int
    overlapping_genes: List[str] = []


class StageAssessment(BaseModel):
    """Stage assessment for a disease category"""
    category: DiseaseCategory
    predicted_stage: DiseaseStage
    confidence: float
    severity_score: float


class BiomarkerAnalysis(BaseModel):
    """Biomarker detection results"""
    detected: List[str]
    missing: List[str]
    top_de_genes: List[str]


class AnalysisResult(BaseModel):
    """Complete analysis result"""
    id: str
    sample_id: str
    created_at: datetime
    status: str

    # Summary
    primary_category: DiseaseCategory
    predicted_stage: DiseaseStage
    severity_score: float
    confidence: float
    risk_level: RiskLevel

    # DE Summary
    total_genes: int
    significant_genes: int
    upregulated: int
    downregulated: int

    # Detailed Results
    top_signatures: List[SignatureMatch]
    stage_assessments: Dict[str, StageAssessment]
    biomarkers: BiomarkerAnalysis
    recommendations: List[str]
    supporting_evidence: List[str]

    # Report
    report_url: Optional[str] = None


class AnalysisListItem(BaseModel):
    """Analysis item for list view"""
    id: str
    sample_id: str
    created_at: datetime
    status: str
    primary_category: Optional[DiseaseCategory] = None
    predicted_stage: Optional[DiseaseStage] = None
    severity_score: Optional[float] = None


class AnalysisListResponse(BaseModel):
    """Paginated list of analyses"""
    items: List[AnalysisListItem]
    total: int
    page: int
    page_size: int


class SignatureListItem(BaseModel):
    """Disease signature for listing"""
    id: str
    name: str
    category: DiseaseCategory
    stage: DiseaseStage
    gene_count: int
    description: Optional[str] = None


class SignatureListResponse(BaseModel):
    """List of disease signatures"""
    items: List[SignatureListItem]
    total: int
    categories: Dict[str, int]


class HealthResponse(BaseModel):
    """Health check response"""
    status: str
    version: str
    database: str
    disease_signatures: int


class ErrorResponse(BaseModel):
    """Error response"""
    detail: str
    error_code: Optional[str] = None
