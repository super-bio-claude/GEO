"""
Reports API Routes
"""
from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, HTMLResponse
from pathlib import Path
import markdown

from app.core.config import settings

router = APIRouter(prefix="/reports", tags=["Reports"])


@router.get("/{filename}")
async def get_report(filename: str, format: str = "md"):
    """
    Get analysis report

    Parameters:
    - filename: Report filename
    - format: Output format (md, html)
    """
    reports_dir = Path(__file__).parent.parent.parent / "reports"
    report_path = reports_dir / filename

    if not report_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Report {filename} not found"
        )

    if format == "html":
        # Convert markdown to HTML
        with open(report_path, 'r') as f:
            md_content = f.read()

        html_content = markdown.markdown(
            md_content,
            extensions=['tables', 'fenced_code']
        )

        # Wrap in basic HTML template
        full_html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <title>Disease Analysis Report</title>
            <style>
                body {{
                    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
                    max-width: 900px;
                    margin: 0 auto;
                    padding: 20px;
                    line-height: 1.6;
                }}
                table {{
                    border-collapse: collapse;
                    width: 100%;
                    margin: 20px 0;
                }}
                th, td {{
                    border: 1px solid #ddd;
                    padding: 12px;
                    text-align: left;
                }}
                th {{
                    background-color: #f4f4f4;
                }}
                h1 {{ color: #2c3e50; }}
                h2 {{ color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
                h3 {{ color: #7f8c8d; }}
                code {{
                    background-color: #f8f8f8;
                    padding: 2px 6px;
                    border-radius: 4px;
                }}
            </style>
        </head>
        <body>
            {html_content}
        </body>
        </html>
        """

        return HTMLResponse(content=full_html)

    # Return raw markdown
    return FileResponse(
        path=report_path,
        media_type="text/markdown",
        filename=filename
    )


@router.get("/")
async def list_reports():
    """List all available reports"""
    reports_dir = Path(__file__).parent.parent.parent / "reports"

    if not reports_dir.exists():
        return {"reports": []}

    reports = []
    for f in reports_dir.glob("*.md"):
        reports.append({
            "filename": f.name,
            "size": f.stat().st_size,
            "modified": f.stat().st_mtime
        })

    return {"reports": sorted(reports, key=lambda x: -x["modified"])}
