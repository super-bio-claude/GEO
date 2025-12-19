import type { Metadata } from 'next'
import { Inter } from 'next/font/google'
import './globals.css'

const inter = Inter({ subsets: ['latin'] })

export const metadata: Metadata = {
  title: 'RNA-seq Disease Analysis Platform',
  description: 'Analyze RNA-seq data for disease similarity and stage prediction',
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html lang="en">
      <body className={inter.className}>
        <div className="min-h-screen">
          {/* Navigation */}
          <nav className="bg-white shadow-sm border-b border-gray-100">
            <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
              <div className="flex justify-between h-16">
                <div className="flex items-center">
                  <div className="flex-shrink-0 flex items-center">
                    <svg className="h-8 w-8 text-primary-600" viewBox="0 0 24 24" fill="currentColor">
                      <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-2 15l-5-5 1.41-1.41L10 14.17l7.59-7.59L19 8l-9 9z"/>
                    </svg>
                    <span className="ml-2 text-xl font-semibold text-gray-900">
                      Disease Analysis
                    </span>
                  </div>
                </div>
                <div className="flex items-center space-x-4">
                  <a href="/" className="text-gray-600 hover:text-gray-900 px-3 py-2 text-sm font-medium">
                    Home
                  </a>
                  <a href="/history" className="text-gray-600 hover:text-gray-900 px-3 py-2 text-sm font-medium">
                    History
                  </a>
                  <a href="/signatures" className="text-gray-600 hover:text-gray-900 px-3 py-2 text-sm font-medium">
                    Signatures
                  </a>
                </div>
              </div>
            </div>
          </nav>

          {/* Main content */}
          <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
            {children}
          </main>

          {/* Footer */}
          <footer className="bg-white border-t border-gray-100 mt-auto">
            <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-4">
              <p className="text-center text-gray-500 text-sm">
                RNA-seq Disease Analysis Platform v1.0
              </p>
            </div>
          </footer>
        </div>
      </body>
    </html>
  )
}
