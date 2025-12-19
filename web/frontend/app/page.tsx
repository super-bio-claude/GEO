'use client'

import { useState } from 'react'
import FileUpload from '@/components/FileUpload'
import AnalysisResult from '@/components/AnalysisResult'
import { AnalysisResultType } from '@/types'

export default function Home() {
  const [result, setResult] = useState<AnalysisResultType | null>(null)
  const [isLoading, setIsLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)

  const handleAnalysisComplete = (data: AnalysisResultType) => {
    setResult(data)
    setError(null)
  }

  const handleError = (errorMessage: string) => {
    setError(errorMessage)
    setResult(null)
  }

  return (
    <div className="space-y-8">
      {/* Hero Section */}
      <div className="text-center">
        <h1 className="text-4xl font-bold text-gray-900 mb-4">
          RNA-seq Disease Analysis
        </h1>
        <p className="text-lg text-gray-600 max-w-2xl mx-auto">
          Upload your differential expression results to identify disease signatures,
          predict disease stage, and receive clinical recommendations.
        </p>
      </div>

      {/* Upload Section */}
      <div className="card max-w-3xl mx-auto">
        <h2 className="text-xl font-semibold text-gray-900 mb-4">
          Upload DE Results
        </h2>
        <FileUpload
          onUploadComplete={handleAnalysisComplete}
          onError={handleError}
          setIsLoading={setIsLoading}
        />

        {error && (
          <div className="mt-4 p-4 bg-red-50 border border-red-200 rounded-lg">
            <p className="text-red-800">{error}</p>
          </div>
        )}
      </div>

      {/* Loading State */}
      {isLoading && (
        <div className="card max-w-3xl mx-auto">
          <div className="flex items-center justify-center py-12">
            <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary-600"></div>
            <span className="ml-4 text-gray-600">Analyzing your data...</span>
          </div>
        </div>
      )}

      {/* Results Section */}
      {result && !isLoading && (
        <AnalysisResult result={result} />
      )}

      {/* Features Section */}
      {!result && !isLoading && (
        <div className="grid md:grid-cols-3 gap-6 mt-12">
          <div className="card text-center">
            <div className="w-12 h-12 bg-purple-100 rounded-lg flex items-center justify-center mx-auto mb-4">
              <svg className="w-6 h-6 text-purple-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
              </svg>
            </div>
            <h3 className="text-lg font-semibold text-gray-900 mb-2">Disease Similarity</h3>
            <p className="text-gray-600 text-sm">
              Compare against 3,600+ disease signatures from MSigDB and curated databases
            </p>
          </div>

          <div className="card text-center">
            <div className="w-12 h-12 bg-blue-100 rounded-lg flex items-center justify-center mx-auto mb-4">
              <svg className="w-6 h-6 text-blue-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" />
              </svg>
            </div>
            <h3 className="text-lg font-semibold text-gray-900 mb-2">Stage Prediction</h3>
            <p className="text-gray-600 text-sm">
              Predict disease stage based on biomarker expression patterns
            </p>
          </div>

          <div className="card text-center">
            <div className="w-12 h-12 bg-green-100 rounded-lg flex items-center justify-center mx-auto mb-4">
              <svg className="w-6 h-6 text-green-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
              </svg>
            </div>
            <h3 className="text-lg font-semibold text-gray-900 mb-2">Clinical Report</h3>
            <p className="text-gray-600 text-sm">
              Generate comprehensive reports with recommendations
            </p>
          </div>
        </div>
      )}
    </div>
  )
}
