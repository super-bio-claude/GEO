"""
Analysis API Routes
"""
import io
import pandas as pd
from fastapi import APIRouter, UploadFile, File, HTTPException, Query, BackgroundTasks
from fastapi.responses import FileResponse
from typing import Optional
from pathlib import Path

from app.models.schemas import (
    AnalysisResult, AnalysisRequest, AnalysisListResponse, AnalysisListItem,
    SignatureListResponse, SignatureListItem, ErrorResponse
)
from app.services.analysis_service import analysis_service
from app.services.supabase_service import supabase_service
from app.core.config import settings

router = APIRouter(prefix="/analysis", tags=["Analysis"])


@router.post("/upload", response_model=AnalysisResult)
async def upload_and_analyze(
    file: UploadFile = File(...),
    sample_id: Optional[str] = Query(None, description="Sample identifier"),
    padj_threshold: float = Query(0.05, ge=0, le=1, description="Adjusted p-value threshold"),
    lfc_threshold: float = Query(0.5, ge=0, description="Log2 fold change threshold")
):
    """
    Upload DE results CSV and run disease analysis

    The CSV file must contain columns:
    - gene: Gene symbol
    - log2FoldChange: Log2 fold change value
    - padj: Adjusted p-value (optional but recommended)

    Returns complete analysis results including disease similarity,
    stage prediction, and recommendations.
    """
    # Validate file type
    if not file.filename.endswith('.csv'):
        raise HTTPException(
            status_code=400,
            detail="File must be a CSV file"
        )

    # Read file content
    try:
        content = await file.read()
        df = pd.read_csv(io.StringIO(content.decode('utf-8')))
    except Exception as e:
        raise HTTPException(
            status_code=400,
            detail=f"Failed to read CSV file: {str(e)}"
        )

    # Validate required columns
    if 'gene' not in df.columns:
        raise HTTPException(
            status_code=400,
            detail="CSV must contain 'gene' column"
        )
    if 'log2FoldChange' not in df.columns:
        raise HTTPException(
            status_code=400,
            detail="CSV must contain 'log2FoldChange' column"
        )

    # Generate sample ID if not provided
    if not sample_id:
        sample_id = Path(file.filename).stem

    # Run analysis
    try:
        result = await analysis_service.analyze_de_results(
            de_data=df,
            sample_id=sample_id,
            padj_threshold=padj_threshold,
            lfc_threshold=lfc_threshold
        )

        # Save to database if connected
        if supabase_service.is_connected():
            await supabase_service.save_analysis(result.model_dump())

        return result

    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Analysis failed: {str(e)}"
        )


@router.get("/{analysis_id}", response_model=AnalysisResult)
async def get_analysis(analysis_id: str):
    """Get analysis result by ID"""
    if supabase_service.is_connected():
        result = await supabase_service.get_analysis(analysis_id)
        if result:
            return result

    raise HTTPException(
        status_code=404,
        detail=f"Analysis {analysis_id} not found"
    )


@router.get("/", response_model=AnalysisListResponse)
async def list_analyses(
    page: int = Query(1, ge=1),
    page_size: int = Query(20, ge=1, le=100),
    status: Optional[str] = Query(None)
):
    """List all analyses with pagination"""
    if not supabase_service.is_connected():
        return AnalysisListResponse(items=[], total=0, page=page, page_size=page_size)

    items, total = await supabase_service.list_analyses(
        page=page,
        page_size=page_size,
        status=status
    )

    return AnalysisListResponse(
        items=[AnalysisListItem(**item) for item in items],
        total=total,
        page=page,
        page_size=page_size
    )


@router.delete("/{analysis_id}")
async def delete_analysis(analysis_id: str):
    """Delete analysis by ID"""
    if not supabase_service.is_connected():
        raise HTTPException(status_code=503, detail="Database not connected")

    success = await supabase_service.delete_analysis(analysis_id)
    if not success:
        raise HTTPException(status_code=404, detail="Analysis not found")

    return {"message": "Analysis deleted successfully"}


@router.get("/signatures/list", response_model=SignatureListResponse)
async def list_signatures(
    category: Optional[str] = Query(None, description="Filter by category"),
    limit: int = Query(100, ge=1, le=1000)
):
    """List available disease signatures"""
    signatures, category_counts = analysis_service.get_signatures(
        category=category,
        limit=limit
    )

    return SignatureListResponse(
        items=[SignatureListItem(**s) for s in signatures],
        total=len(signatures),
        categories=category_counts
    )
