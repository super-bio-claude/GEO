# RNA-seq Disease Analysis Web Platform

A web-based platform for analyzing RNA-seq differential expression data to identify disease signatures and predict disease progression.

## Architecture

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│                 │     │                 │     │                 │
│   Next.js       │────▶│   FastAPI       │────▶│   Supabase      │
│   Frontend      │     │   Backend       │     │   Database      │
│   (Port 3000)   │     │   (Port 8000)   │     │                 │
│                 │     │                 │     │                 │
└─────────────────┘     └────────┬────────┘     └─────────────────┘
                                 │
                                 ▼
                        ┌─────────────────┐
                        │                 │
                        │   ChromaDB      │
                        │   Vector DB     │
                        │                 │
                        └─────────────────┘
```

## Features

- **File Upload**: Drag & drop CSV upload with differential expression results
- **Disease Similarity**: Compare against 3,600+ disease signatures using vector search
- **Stage Prediction**: Predict disease stage (preclinical, early, intermediate, advanced)
- **Multi-category Assessment**: Cancer, Neurological, Metabolic, Aging
- **Clinical Reports**: Generate comprehensive Markdown/HTML reports
- **Biomarker Detection**: Identify relevant disease biomarkers

## Quick Start

### 1. Setup Supabase

1. Create a new project at [supabase.com](https://supabase.com)
2. Run the SQL schema in the SQL editor:
   ```bash
   cat backend/supabase_schema.sql
   ```
3. Copy your project URL and anon key

### 2. Start Backend

```bash
cd backend

# Create virtual environment
python -m venv venv
source venv/bin/activate  # or `venv\Scripts\activate` on Windows

# Install dependencies
pip install -r requirements.txt

# Create .env file
cp .env.example .env
# Edit .env with your Supabase credentials

# Start server
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

API docs available at: http://localhost:8000/docs

### 3. Start Frontend

```bash
cd frontend

# Install dependencies
npm install

# Create .env.local
cp .env.example .env.local

# Start development server
npm run dev
```

Open http://localhost:3000

## API Endpoints

### Analysis

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/v1/analysis/upload` | Upload CSV and run analysis |
| GET | `/api/v1/analysis/{id}` | Get analysis result by ID |
| GET | `/api/v1/analysis/` | List all analyses |
| DELETE | `/api/v1/analysis/{id}` | Delete analysis |

### Reports

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/v1/reports/{filename}` | Get report (md or html) |
| GET | `/api/v1/reports/` | List all reports |

### Other

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/health` | Health check |
| GET | `/api/v1/stats` | Platform statistics |
| GET | `/api/v1/analysis/signatures/list` | List disease signatures |

## Input File Format

The CSV file must contain:
- `gene`: Gene symbol (required)
- `log2FoldChange`: Log2 fold change value (required)
- `padj`: Adjusted p-value (optional but recommended)

Example:
```csv
gene,log2FoldChange,padj,baseMean
GFAP,1.23,0.001,1500
FOSB,-2.45,0.0001,800
NPAS4,-3.12,0.00001,500
...
```

## Response Format

```json
{
  "id": "uuid",
  "sample_id": "GSE125583_Alzheimer",
  "primary_category": "neurological",
  "predicted_stage": "advanced",
  "severity_score": 0.70,
  "confidence": 0.47,
  "risk_level": "HIGH",
  "top_signatures": [...],
  "biomarkers": {
    "detected": ["GFAP", "FOSB_down"],
    "missing": ["BDNF"]
  },
  "recommendations": [
    "Neuroinflammation markers present - anti-inflammatory therapy consideration"
  ]
}
```

## Development

### Backend Structure

```
backend/
├── app/
│   ├── api/              # API routes
│   ├── core/             # Configuration
│   ├── models/           # Pydantic schemas
│   └── services/         # Business logic
├── uploads/              # Uploaded files
├── reports/              # Generated reports
└── requirements.txt
```

### Frontend Structure

```
frontend/
├── app/                  # Next.js App Router
├── components/           # React components
├── lib/                  # Utilities
└── types/                # TypeScript types
```

## Deployment

### Docker (Recommended)

```dockerfile
# backend/Dockerfile
FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
```

### Environment Variables

**Backend:**
- `SUPABASE_URL`: Your Supabase project URL
- `SUPABASE_KEY`: Your Supabase anon key
- `DEBUG`: Set to `false` in production

**Frontend:**
- `NEXT_PUBLIC_API_URL`: Backend API URL
- `NEXT_PUBLIC_SUPABASE_URL`: (Optional) Direct Supabase access

## License

MIT
