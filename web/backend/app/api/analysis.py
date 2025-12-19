"""
Analysis API Routes
"""
import io
import os
import subprocess
import tempfile
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

# Project root path
PROJECT_ROOT = Path(__file__).parent.parent.parent.parent.parent


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


@router.post("/upload-counts", response_model=AnalysisResult)
async def upload_counts_and_analyze(
    file: UploadFile = File(...),
    sample_id: Optional[str] = Query(None, description="Sample identifier"),
    control_pattern: str = Query("control", description="Pattern to identify control samples"),
    treatment_pattern: str = Query("compress", description="Pattern to identify treatment samples"),
    padj_threshold: float = Query(0.05, ge=0, le=1, description="Adjusted p-value threshold"),
    lfc_threshold: float = Query(0.5, ge=0, description="Log2 fold change threshold")
):
    """
    Upload raw counts matrix (TXT/CSV) and run full analysis pipeline:
    1. Parse sample names to identify conditions
    2. Run DESeq2 differential expression analysis
    3. Run disease similarity analysis

    The file must be a tab-separated counts matrix with:
    - First row: sample names
    - First column: gene names
    - Values: raw counts (integers)
    """
    # Validate file type
    if not (file.filename.endswith('.txt') or file.filename.endswith('.csv') or file.filename.endswith('.tsv')):
        raise HTTPException(
            status_code=400,
            detail="File must be a TXT, CSV, or TSV file"
        )

    # Read file content
    try:
        content = await file.read()
        content_str = content.decode('utf-8')

        # Detect separator
        first_line = content_str.split('\n')[0]
        if '\t' in first_line:
            sep = '\t'
        elif ',' in first_line:
            sep = ','
        else:
            sep = '\t'

        df = pd.read_csv(io.StringIO(content_str), sep=sep, index_col=0)
    except Exception as e:
        raise HTTPException(
            status_code=400,
            detail=f"Failed to read counts file: {str(e)}"
        )

    # Generate sample ID if not provided
    if not sample_id:
        sample_id = Path(file.filename).stem

    # Create temporary directory for analysis
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Save counts matrix
        counts_file = temp_path / "counts.txt"
        df.to_csv(counts_file, sep='\t')

        # Generate metadata from sample names
        samples = df.columns.tolist()
        metadata = []

        for sample in samples:
            sample_lower = sample.lower()
            if control_pattern.lower() in sample_lower:
                condition = "control"
            elif treatment_pattern.lower() in sample_lower:
                condition = "treatment"
            else:
                # Default to treatment if no pattern matches
                condition = "treatment"
            metadata.append({"sample": sample, "condition": condition})

        metadata_df = pd.DataFrame(metadata)
        metadata_df.set_index('sample', inplace=True)

        # Check if we have both conditions
        conditions = metadata_df['condition'].unique().tolist()
        if len(conditions) < 2:
            raise HTTPException(
                status_code=400,
                detail=f"Need both control and treatment samples. Found only: {conditions}. "
                       f"Adjust control_pattern ('{control_pattern}') or treatment_pattern ('{treatment_pattern}')"
            )

        metadata_file = temp_path / "metadata.csv"
        metadata_df.to_csv(metadata_file)

        # Output directory
        output_dir = temp_path / "deseq2_output"
        output_dir.mkdir(exist_ok=True)

        # Run DESeq2 R script
        r_script = PROJECT_ROOT / "src" / "R_scripts" / "deseq2_analysis.R"

        try:
            result = subprocess.run(
                [
                    "Rscript",
                    str(r_script),
                    str(counts_file),
                    str(metadata_file),
                    str(output_dir)
                ],
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )

            if result.returncode != 0:
                raise HTTPException(
                    status_code=500,
                    detail=f"DESeq2 analysis failed: {result.stderr}"
                )
        except subprocess.TimeoutExpired:
            raise HTTPException(
                status_code=500,
                detail="DESeq2 analysis timed out (> 5 minutes)"
            )
        except FileNotFoundError:
            raise HTTPException(
                status_code=500,
                detail="Rscript not found. Please ensure R is installed."
            )

        # Read DESeq2 results
        deseq2_results_file = output_dir / "deseq2_results_apeglm.csv"

        if not deseq2_results_file.exists():
            raise HTTPException(
                status_code=500,
                detail="DESeq2 results file not found. Analysis may have failed."
            )

        de_results = pd.read_csv(deseq2_results_file)

        # Run disease analysis
        try:
            analysis_result = await analysis_service.analyze_de_results(
                de_data=de_results,
                sample_id=sample_id,
                padj_threshold=padj_threshold,
                lfc_threshold=lfc_threshold
            )

            # Save to database if connected
            if supabase_service.is_connected():
                await supabase_service.save_analysis(analysis_result.model_dump())

            return analysis_result

        except Exception as e:
            raise HTTPException(
                status_code=500,
                detail=f"Disease analysis failed: {str(e)}"
            )
