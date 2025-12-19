-- RNA-seq Disease Analysis Platform - Supabase Schema
-- Run this SQL in your Supabase SQL editor to set up the database

-- Enable UUID extension
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- Analyses table
CREATE TABLE IF NOT EXISTS analyses (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    sample_id VARCHAR(255) NOT NULL,
    status VARCHAR(50) DEFAULT 'pending',
    error_message TEXT,

    -- Disease Assessment
    primary_category VARCHAR(50),
    predicted_stage VARCHAR(50),
    severity_score DECIMAL(5, 3),
    confidence DECIMAL(5, 3),
    risk_level VARCHAR(20),

    -- DE Summary
    total_genes INTEGER,
    significant_genes INTEGER,
    upregulated INTEGER,
    downregulated INTEGER,

    -- Detailed Results (JSON)
    top_signatures JSONB DEFAULT '[]'::JSONB,
    stage_assessments JSONB DEFAULT '{}'::JSONB,
    biomarkers JSONB DEFAULT '{}'::JSONB,

    -- Arrays
    recommendations TEXT[] DEFAULT '{}',
    supporting_evidence TEXT[] DEFAULT '{}',

    -- Report
    report_url VARCHAR(500),

    -- Timestamps
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),

    -- User (optional, for future auth)
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE
);

-- Create indexes
CREATE INDEX IF NOT EXISTS idx_analyses_sample_id ON analyses(sample_id);
CREATE INDEX IF NOT EXISTS idx_analyses_status ON analyses(status);
CREATE INDEX IF NOT EXISTS idx_analyses_created_at ON analyses(created_at DESC);
CREATE INDEX IF NOT EXISTS idx_analyses_primary_category ON analyses(primary_category);
CREATE INDEX IF NOT EXISTS idx_analyses_user_id ON analyses(user_id);

-- Update timestamp trigger
CREATE OR REPLACE FUNCTION update_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER analyses_updated_at
    BEFORE UPDATE ON analyses
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at();

-- Row Level Security (RLS)
ALTER TABLE analyses ENABLE ROW LEVEL SECURITY;

-- Policy: Allow anonymous access for now (adjust for production)
CREATE POLICY "Allow all access to analyses"
    ON analyses
    FOR ALL
    USING (true)
    WITH CHECK (true);

-- Optional: User-specific policy (uncomment for authenticated access)
-- CREATE POLICY "Users can view own analyses"
--     ON analyses
--     FOR SELECT
--     USING (auth.uid() = user_id);

-- CREATE POLICY "Users can insert own analyses"
--     ON analyses
--     FOR INSERT
--     WITH CHECK (auth.uid() = user_id);

-- Storage bucket for uploaded files (run in Supabase dashboard)
-- INSERT INTO storage.buckets (id, name, public)
-- VALUES ('uploads', 'uploads', false);

-- INSERT INTO storage.buckets (id, name, public)
-- VALUES ('reports', 'reports', true);

-- Sample query to verify
-- SELECT * FROM analyses ORDER BY created_at DESC LIMIT 10;

COMMENT ON TABLE analyses IS 'Stores RNA-seq disease analysis results';
COMMENT ON COLUMN analyses.primary_category IS 'Primary disease category: cancer, neurological, metabolic, aging, etc.';
COMMENT ON COLUMN analyses.predicted_stage IS 'Predicted disease stage: preclinical, early, intermediate, advanced';
COMMENT ON COLUMN analyses.severity_score IS 'Severity score from 0 to 1';
COMMENT ON COLUMN analyses.top_signatures IS 'Top matching disease signatures as JSON array';
