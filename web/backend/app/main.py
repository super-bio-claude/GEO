"""
RNA-seq Disease Analysis Platform - FastAPI Backend
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from contextlib import asynccontextmanager
from pathlib import Path

from app.core.config import settings
from app.api import analysis, reports
from app.services.analysis_service import analysis_service
from app.services.supabase_service import supabase_service


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan events"""
    # Startup
    print("Starting RNA-seq Disease Analysis Platform...")
    print(f"Disease signatures loaded: {analysis_service.get_signature_count()}")

    # Create directories
    Path("uploads").mkdir(exist_ok=True)
    Path("reports").mkdir(exist_ok=True)

    yield

    # Shutdown
    print("Shutting down...")


# Create FastAPI app
app = FastAPI(
    title=settings.APP_NAME,
    version=settings.APP_VERSION,
    description="""
    RNA-seq Disease Analysis Platform API

    This API provides endpoints for:
    - Uploading differential expression results
    - Running disease similarity analysis
    - Disease stage prediction
    - Generating clinical reports

    ## Quick Start
    1. Upload a CSV file with DE results using `/api/v1/analysis/upload`
    2. The file must contain `gene` and `log2FoldChange` columns
    3. Optional: Include `padj` column for filtering significant genes
    """,
    lifespan=lifespan,
    docs_url="/docs",
    redoc_url="/redoc"
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(analysis.router, prefix=settings.API_V1_PREFIX)
app.include_router(reports.router, prefix=settings.API_V1_PREFIX)


@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "name": settings.APP_NAME,
        "version": settings.APP_VERSION,
        "docs": "/docs",
        "api": settings.API_V1_PREFIX
    }


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "version": settings.APP_VERSION,
        "database": "connected" if supabase_service.is_connected() else "not configured",
        "disease_signatures": analysis_service.get_signature_count()
    }


@app.get(f"{settings.API_V1_PREFIX}/stats")
async def get_stats():
    """Get platform statistics"""
    return {
        "disease_signatures": analysis_service.get_signature_count(),
        "categories": {
            "cancer": "Cancer-related signatures",
            "neurological": "Neurological disease signatures",
            "metabolic": "Metabolic disease signatures",
            "aging": "Aging-related signatures",
            "inflammatory": "Inflammatory signatures",
            "other": "Other perturbation signatures"
        }
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "app.main:app",
        host="0.0.0.0",
        port=8000,
        reload=True
    )
