"""
Application Configuration
"""
from pydantic_settings import BaseSettings
from functools import lru_cache
from pathlib import Path
from typing import Optional


class Settings(BaseSettings):
    """Application settings loaded from environment variables"""

    # App Settings
    APP_NAME: str = "RNA-seq Disease Analysis Platform"
    APP_VERSION: str = "1.0.0"
    DEBUG: bool = False

    # API Settings
    API_V1_PREFIX: str = "/api/v1"

    # CORS
    CORS_ORIGINS: list = ["http://localhost:3000", "http://127.0.0.1:3000"]

    # Supabase
    SUPABASE_URL: str = ""
    SUPABASE_KEY: str = ""
    SUPABASE_SERVICE_KEY: Optional[str] = None

    # File Storage
    UPLOAD_DIR: Path = Path("uploads")
    REPORTS_DIR: Path = Path("reports")
    MAX_FILE_SIZE: int = 100 * 1024 * 1024  # 100MB

    # Disease Analysis
    DISEASE_DB_PATH: Path = Path("../../data/disease_db")
    CHROMADB_PATH: Path = Path("../../data/chromadb_disease")

    # Background Tasks
    REDIS_URL: str = "redis://localhost:6379/0"

    # Analysis Settings
    DEFAULT_PADJ_THRESHOLD: float = 0.05
    DEFAULT_LFC_THRESHOLD: float = 0.5

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"
        case_sensitive = True


@lru_cache()
def get_settings() -> Settings:
    """Get cached settings instance"""
    return Settings()


settings = get_settings()
