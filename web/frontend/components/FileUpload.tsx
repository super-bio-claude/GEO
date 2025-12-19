'use client'

import { useCallback, useState } from 'react'
import { useDropzone } from 'react-dropzone'
import { AnalysisResultType } from '@/types'

const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'

type FileType = 'de_results' | 'raw_counts'

interface FileUploadProps {
  onUploadComplete: (result: AnalysisResultType) => void
  onError: (error: string) => void
  setIsLoading: (loading: boolean) => void
}

export default function FileUpload({ onUploadComplete, onError, setIsLoading }: FileUploadProps) {
  const [fileType, setFileType] = useState<FileType>('raw_counts')
  const [sampleId, setSampleId] = useState('')
  const [padjThreshold, setPadjThreshold] = useState(0.05)
  const [lfcThreshold, setLfcThreshold] = useState(0.5)
  const [controlPattern, setControlPattern] = useState('control')
  const [treatmentPattern, setTreatmentPattern] = useState('compress')
  const [selectedFile, setSelectedFile] = useState<File | null>(null)

  const acceptedFileTypes = fileType === 'de_results'
    ? { 'text/csv': ['.csv'] }
    : { 'text/plain': ['.txt', '.tsv'], 'text/csv': ['.csv'] }

  const onDrop = useCallback((acceptedFiles: File[]) => {
    if (acceptedFiles.length > 0) {
      setSelectedFile(acceptedFiles[0])
    }
  }, [])

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: acceptedFileTypes,
    multiple: false
  })

  const handleSubmit = async () => {
    if (!selectedFile) {
      onError('Please select a file first')
      return
    }

    setIsLoading(true)

    try {
      const formData = new FormData()
      formData.append('file', selectedFile)

      let endpoint: string
      const params = new URLSearchParams({
        padj_threshold: padjThreshold.toString(),
        lfc_threshold: lfcThreshold.toString(),
      })

      if (sampleId) {
        params.append('sample_id', sampleId)
      }

      if (fileType === 'raw_counts') {
        endpoint = '/api/v1/analysis/upload-counts'
        params.append('control_pattern', controlPattern)
        params.append('treatment_pattern', treatmentPattern)
      } else {
        endpoint = '/api/v1/analysis/upload'
      }

      const response = await fetch(
        `${API_URL}${endpoint}?${params}`,
        {
          method: 'POST',
          body: formData,
        }
      )

      if (!response.ok) {
        const error = await response.json()
        throw new Error(error.detail || 'Analysis failed')
      }

      const result = await response.json()
      onUploadComplete(result)
    } catch (err) {
      onError(err instanceof Error ? err.message : 'An error occurred')
    } finally {
      setIsLoading(false)
    }
  }

  return (
    <div className="space-y-6">
      {/* File Type Selection */}
      <div className="flex gap-4">
        <button
          onClick={() => { setFileType('raw_counts'); setSelectedFile(null) }}
          className={`flex-1 py-3 px-4 rounded-lg font-medium transition-colors border-2
            ${fileType === 'raw_counts'
              ? 'border-primary-500 bg-primary-50 text-primary-700'
              : 'border-gray-200 bg-white text-gray-600 hover:border-gray-300'}`}
        >
          <div className="text-lg">Raw Counts Matrix</div>
          <div className="text-xs mt-1 opacity-75">TXT/CSV file with gene counts</div>
        </button>
        <button
          onClick={() => { setFileType('de_results'); setSelectedFile(null) }}
          className={`flex-1 py-3 px-4 rounded-lg font-medium transition-colors border-2
            ${fileType === 'de_results'
              ? 'border-primary-500 bg-primary-50 text-primary-700'
              : 'border-gray-200 bg-white text-gray-600 hover:border-gray-300'}`}
        >
          <div className="text-lg">DE Results</div>
          <div className="text-xs mt-1 opacity-75">CSV with gene, log2FC, padj</div>
        </button>
      </div>

      {/* Dropzone */}
      <div
        {...getRootProps()}
        className={`border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-colors
          ${isDragActive ? 'border-primary-500 bg-primary-50' : 'border-gray-300 hover:border-gray-400'}
          ${selectedFile ? 'bg-green-50 border-green-300' : ''}`}
      >
        <input {...getInputProps()} />
        {selectedFile ? (
          <div>
            <svg className="w-12 h-12 text-green-500 mx-auto mb-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" />
            </svg>
            <p className="text-lg font-medium text-green-700">{selectedFile.name}</p>
            <p className="text-sm text-gray-500 mt-1">
              {(selectedFile.size / 1024).toFixed(1)} KB
            </p>
            <p className="text-sm text-primary-600 mt-2">Click or drag to replace</p>
          </div>
        ) : (
          <div>
            <svg className="w-12 h-12 text-gray-400 mx-auto mb-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M15 13l-3-3m0 0l-3 3m3-3v12" />
            </svg>
            <p className="text-lg font-medium text-gray-700">
              {isDragActive ? 'Drop the file here...' : `Drag & drop your ${fileType === 'raw_counts' ? 'counts matrix' : 'CSV'} file here`}
            </p>
            <p className="text-sm text-gray-500 mt-1">or click to browse</p>
            <p className="text-xs text-gray-400 mt-4">
              {fileType === 'raw_counts'
                ? 'Accepts: TXT, TSV, CSV (tab-separated counts matrix)'
                : 'Required columns: gene, log2FoldChange, padj (optional)'}
            </p>
          </div>
        )}
      </div>

      {/* Options */}
      <div className="grid md:grid-cols-2 lg:grid-cols-3 gap-4">
        <div>
          <label className="block text-sm font-medium text-gray-700 mb-1">
            Sample ID (optional)
          </label>
          <input
            type="text"
            value={sampleId}
            onChange={(e) => setSampleId(e.target.value)}
            placeholder="e.g., GSE123456"
            className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-primary-500 focus:border-primary-500"
          />
        </div>

        {fileType === 'raw_counts' && (
          <>
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Control pattern
              </label>
              <input
                type="text"
                value={controlPattern}
                onChange={(e) => setControlPattern(e.target.value)}
                placeholder="e.g., control, ctrl, WT"
                className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-primary-500 focus:border-primary-500"
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Treatment pattern
              </label>
              <input
                type="text"
                value={treatmentPattern}
                onChange={(e) => setTreatmentPattern(e.target.value)}
                placeholder="e.g., treat, disease, KO"
                className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-primary-500 focus:border-primary-500"
              />
            </div>
          </>
        )}

        <div>
          <label className="block text-sm font-medium text-gray-700 mb-1">
            Adjusted p-value threshold
          </label>
          <input
            type="number"
            value={padjThreshold}
            onChange={(e) => setPadjThreshold(parseFloat(e.target.value))}
            min={0}
            max={1}
            step={0.01}
            className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-primary-500 focus:border-primary-500"
          />
        </div>
        <div>
          <label className="block text-sm font-medium text-gray-700 mb-1">
            Log2 fold change threshold
          </label>
          <input
            type="number"
            value={lfcThreshold}
            onChange={(e) => setLfcThreshold(parseFloat(e.target.value))}
            min={0}
            step={0.1}
            className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-primary-500 focus:border-primary-500"
          />
        </div>
      </div>

      {/* Info box for raw counts */}
      {fileType === 'raw_counts' && (
        <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
          <div className="flex items-start">
            <svg className="w-5 h-5 text-blue-500 mt-0.5 mr-3 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
            </svg>
            <div className="text-sm text-blue-800">
              <p className="font-medium mb-1">Full Pipeline Mode</p>
              <p>This will run DESeq2 differential expression analysis followed by disease similarity analysis.
                 Sample conditions are detected from column names using the patterns above.</p>
            </div>
          </div>
        </div>
      )}

      {/* Submit Button */}
      <button
        onClick={handleSubmit}
        disabled={!selectedFile}
        className={`w-full py-3 rounded-lg font-medium transition-colors
          ${selectedFile
            ? 'bg-primary-600 hover:bg-primary-700 text-white'
            : 'bg-gray-200 text-gray-500 cursor-not-allowed'}`}
      >
        {fileType === 'raw_counts' ? 'Run Full Analysis (DESeq2 + Disease Analysis)' : 'Analyze'}
      </button>
    </div>
  )
}
