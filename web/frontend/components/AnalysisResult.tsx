'use client'

import { AnalysisResultType, DiseaseCategory, RiskLevel } from '@/types'

const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'

interface AnalysisResultProps {
  result: AnalysisResultType
}

const categoryColors: Record<DiseaseCategory, string> = {
  cancer: 'bg-red-100 text-red-800',
  neurological: 'bg-purple-100 text-purple-800',
  metabolic: 'bg-amber-100 text-amber-800',
  aging: 'bg-indigo-100 text-indigo-800',
  inflammatory: 'bg-pink-100 text-pink-800',
  other: 'bg-gray-100 text-gray-800',
}

const riskColors: Record<RiskLevel, string> = {
  HIGH: 'bg-red-100 text-red-800 border-red-200',
  MODERATE: 'bg-yellow-100 text-yellow-800 border-yellow-200',
  LOW: 'bg-green-100 text-green-800 border-green-200',
}

export default function AnalysisResult({ result }: AnalysisResultProps) {
  return (
    <div className="space-y-6">
      {/* Summary Card */}
      <div className="card">
        <div className="flex items-start justify-between">
          <div>
            <h2 className="text-2xl font-bold text-gray-900">{result.sample_id}</h2>
            <p className="text-sm text-gray-500 mt-1">
              Analyzed at {new Date(result.created_at).toLocaleString()}
            </p>
          </div>
          <span className={`badge px-4 py-2 text-base font-semibold border ${riskColors[result.risk_level]}`}>
            {result.risk_level} RISK
          </span>
        </div>

        {/* Key Metrics */}
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mt-6">
          <div className="bg-gray-50 rounded-lg p-4">
            <p className="text-sm text-gray-500">Primary Category</p>
            <p className="text-lg font-semibold text-gray-900 capitalize mt-1">
              {result.primary_category}
            </p>
          </div>
          <div className="bg-gray-50 rounded-lg p-4">
            <p className="text-sm text-gray-500">Predicted Stage</p>
            <p className="text-lg font-semibold text-gray-900 capitalize mt-1">
              {result.predicted_stage}
            </p>
          </div>
          <div className="bg-gray-50 rounded-lg p-4">
            <p className="text-sm text-gray-500">Severity Score</p>
            <p className="text-lg font-semibold text-gray-900 mt-1">
              {(result.severity_score * 100).toFixed(0)}%
            </p>
            <div className="w-full bg-gray-200 rounded-full h-2 mt-2">
              <div
                className={`h-2 rounded-full ${
                  result.severity_score > 0.7 ? 'bg-red-500' :
                  result.severity_score > 0.4 ? 'bg-yellow-500' : 'bg-green-500'
                }`}
                style={{ width: `${result.severity_score * 100}%` }}
              />
            </div>
          </div>
          <div className="bg-gray-50 rounded-lg p-4">
            <p className="text-sm text-gray-500">Confidence</p>
            <p className="text-lg font-semibold text-gray-900 mt-1">
              {(result.confidence * 100).toFixed(0)}%
            </p>
          </div>
        </div>
      </div>

      {/* DE Summary */}
      <div className="card">
        <h3 className="text-lg font-semibold text-gray-900 mb-4">Differential Expression Summary</h3>
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
          <div className="text-center p-3 bg-gray-50 rounded-lg">
            <p className="text-2xl font-bold text-gray-900">{result.total_genes.toLocaleString()}</p>
            <p className="text-sm text-gray-500">Total Genes</p>
          </div>
          <div className="text-center p-3 bg-blue-50 rounded-lg">
            <p className="text-2xl font-bold text-blue-600">{result.significant_genes.toLocaleString()}</p>
            <p className="text-sm text-gray-500">Significant</p>
          </div>
          <div className="text-center p-3 bg-green-50 rounded-lg">
            <p className="text-2xl font-bold text-green-600">{result.upregulated.toLocaleString()}</p>
            <p className="text-sm text-gray-500">Upregulated</p>
          </div>
          <div className="text-center p-3 bg-red-50 rounded-lg">
            <p className="text-2xl font-bold text-red-600">{result.downregulated.toLocaleString()}</p>
            <p className="text-sm text-gray-500">Downregulated</p>
          </div>
        </div>
      </div>

      {/* Top Signatures */}
      <div className="card">
        <h3 className="text-lg font-semibold text-gray-900 mb-4">Top Disease Signatures</h3>
        <div className="overflow-x-auto">
          <table className="min-w-full divide-y divide-gray-200">
            <thead className="bg-gray-50">
              <tr>
                <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Rank</th>
                <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Signature</th>
                <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Category</th>
                <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Score</th>
                <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Gene Overlap</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-200">
              {result.top_signatures.slice(0, 10).map((sig) => (
                <tr key={sig.signature_id} className="hover:bg-gray-50">
                  <td className="px-4 py-3 text-sm text-gray-500">{sig.rank}</td>
                  <td className="px-4 py-3 text-sm font-medium text-gray-900">{sig.signature_name}</td>
                  <td className="px-4 py-3">
                    <span className={`badge ${categoryColors[sig.category]}`}>
                      {sig.category}
                    </span>
                  </td>
                  <td className="px-4 py-3 text-sm text-gray-900">{sig.score.toFixed(3)}</td>
                  <td className="px-4 py-3 text-sm text-gray-500">{sig.gene_overlap}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Biomarkers */}
      <div className="grid md:grid-cols-2 gap-6">
        <div className="card">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Biomarkers Detected</h3>
          <div className="flex flex-wrap gap-2">
            {result.biomarkers.detected.length > 0 ? (
              result.biomarkers.detected.map((marker) => (
                <span key={marker} className="badge bg-green-100 text-green-800">
                  {marker}
                </span>
              ))
            ) : (
              <p className="text-gray-500 text-sm">No specific biomarkers detected</p>
            )}
          </div>
        </div>
        <div className="card">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Top DE Genes</h3>
          <div className="flex flex-wrap gap-2">
            {result.biomarkers.top_de_genes.map((gene) => (
              <span key={gene} className="badge bg-blue-100 text-blue-800">
                {gene}
              </span>
            ))}
          </div>
        </div>
      </div>

      {/* Recommendations */}
      <div className="card">
        <h3 className="text-lg font-semibold text-gray-900 mb-4">Recommendations</h3>
        <ul className="space-y-2">
          {result.recommendations.map((rec, index) => (
            <li key={index} className="flex items-start">
              <svg className="w-5 h-5 text-primary-500 mr-2 mt-0.5 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" />
              </svg>
              <span className="text-gray-700">{rec}</span>
            </li>
          ))}
        </ul>
      </div>

      {/* Supporting Evidence */}
      {result.supporting_evidence.length > 0 && (
        <div className="card">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Supporting Evidence</h3>
          <ul className="space-y-1">
            {result.supporting_evidence.map((evidence, index) => (
              <li key={index} className="text-sm text-gray-600">
                • {evidence}
              </li>
            ))}
          </ul>
        </div>
      )}

      {/* Download Report */}
      {result.report_url && (
        <div className="flex justify-center gap-4">
          <a
            href={`${API_URL}${result.report_url}?format=md`}
            download
            className="btn-primary inline-flex items-center"
          >
            <svg className="w-5 h-5 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
            </svg>
            Download Report (MD)
          </a>
          <a
            href={`${API_URL}${result.report_url}?format=html`}
            target="_blank"
            rel="noopener noreferrer"
            className="bg-gray-100 hover:bg-gray-200 text-gray-700 font-medium py-2 px-4 rounded-lg transition-colors inline-flex items-center"
          >
            <svg className="w-5 h-5 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
            </svg>
            View Report (HTML)
          </a>
        </div>
      )}
    </div>
  )
}
