"""
Supabase Database Service
"""
from typing import Dict, List, Optional
from datetime import datetime
import json
from supabase import create_client, Client

from app.core.config import settings


class SupabaseService:
    """Service for Supabase database operations"""

    def __init__(self):
        self.client: Optional[Client] = None

    def connect(self) -> Client:
        """Connect to Supabase"""
        if self.client is None:
            self.client = create_client(
                settings.SUPABASE_URL,
                settings.SUPABASE_KEY
            )
        return self.client

    async def save_analysis(self, analysis_data: Dict) -> Dict:
        """Save analysis result to database"""
        client = self.connect()

        # Prepare data for insertion
        record = {
            "id": analysis_data["id"],
            "sample_id": analysis_data["sample_id"],
            "status": analysis_data["status"],
            "primary_category": analysis_data["primary_category"],
            "predicted_stage": analysis_data["predicted_stage"],
            "severity_score": analysis_data["severity_score"],
            "confidence": analysis_data["confidence"],
            "risk_level": analysis_data["risk_level"],
            "total_genes": analysis_data["total_genes"],
            "significant_genes": analysis_data["significant_genes"],
            "upregulated": analysis_data["upregulated"],
            "downregulated": analysis_data["downregulated"],
            "top_signatures": json.dumps(analysis_data.get("top_signatures", [])),
            "stage_assessments": json.dumps(analysis_data.get("stage_assessments", {})),
            "biomarkers": json.dumps(analysis_data.get("biomarkers", {})),
            "recommendations": analysis_data.get("recommendations", []),
            "supporting_evidence": analysis_data.get("supporting_evidence", []),
            "report_url": analysis_data.get("report_url"),
            "created_at": datetime.utcnow().isoformat()
        }

        result = client.table("analyses").insert(record).execute()
        return result.data[0] if result.data else None

    async def get_analysis(self, analysis_id: str) -> Optional[Dict]:
        """Get analysis by ID"""
        client = self.connect()

        result = client.table("analyses").select("*").eq("id", analysis_id).execute()

        if result.data:
            record = result.data[0]
            # Parse JSON fields
            record["top_signatures"] = json.loads(record.get("top_signatures", "[]"))
            record["stage_assessments"] = json.loads(record.get("stage_assessments", "{}"))
            record["biomarkers"] = json.loads(record.get("biomarkers", "{}"))
            return record

        return None

    async def list_analyses(
        self,
        page: int = 1,
        page_size: int = 20,
        status: Optional[str] = None
    ) -> tuple[List[Dict], int]:
        """List analyses with pagination"""
        client = self.connect()

        query = client.table("analyses").select("*", count="exact")

        if status:
            query = query.eq("status", status)

        # Pagination
        offset = (page - 1) * page_size
        query = query.order("created_at", desc=True).range(offset, offset + page_size - 1)

        result = query.execute()

        return result.data or [], result.count or 0

    async def delete_analysis(self, analysis_id: str) -> bool:
        """Delete analysis by ID"""
        client = self.connect()

        result = client.table("analyses").delete().eq("id", analysis_id).execute()
        return len(result.data) > 0

    async def update_analysis_status(
        self,
        analysis_id: str,
        status: str,
        error_message: Optional[str] = None
    ) -> Optional[Dict]:
        """Update analysis status"""
        client = self.connect()

        update_data = {"status": status}
        if error_message:
            update_data["error_message"] = error_message

        result = client.table("analyses").update(update_data).eq("id", analysis_id).execute()
        return result.data[0] if result.data else None

    def is_connected(self) -> bool:
        """Check if connected to Supabase"""
        try:
            if not settings.SUPABASE_URL or not settings.SUPABASE_KEY:
                return False
            client = self.connect()
            # Simple health check
            client.table("analyses").select("id").limit(1).execute()
            return True
        except Exception:
            return False


# Singleton instance
supabase_service = SupabaseService()
