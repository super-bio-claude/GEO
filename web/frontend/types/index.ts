export type DiseaseCategory = 'cancer' | 'neurological' | 'metabolic' | 'aging' | 'inflammatory' | 'other'
export type DiseaseStage = 'preclinical' | 'early' | 'intermediate' | 'advanced' | 'unknown'
export type RiskLevel = 'LOW' | 'MODERATE' | 'HIGH'

export interface SignatureMatch {
  rank: number
  signature_id: string
  signature_name: string
  category: DiseaseCategory
  stage: DiseaseStage
  score: number
  gene_overlap: number
  overlapping_genes: string[]
}

export interface StageAssessment {
  category: DiseaseCategory
  predicted_stage: DiseaseStage
  confidence: number
  severity_score: number
}

export interface BiomarkerAnalysis {
  detected: string[]
  missing: string[]
  top_de_genes: string[]
}

export interface AnalysisResultType {
  id: string
  sample_id: string
  created_at: string
  status: string
  primary_category: DiseaseCategory
  predicted_stage: DiseaseStage
  severity_score: number
  confidence: number
  risk_level: RiskLevel
  total_genes: number
  significant_genes: number
  upregulated: number
  downregulated: number
  top_signatures: SignatureMatch[]
  stage_assessments: Record<string, StageAssessment>
  biomarkers: BiomarkerAnalysis
  recommendations: string[]
  supporting_evidence: string[]
  report_url?: string
}
