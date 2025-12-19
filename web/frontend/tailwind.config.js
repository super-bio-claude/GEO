/** @type {import('tailwindcss').Config} */
module.exports = {
  content: [
    './app/**/*.{js,ts,jsx,tsx,mdx}',
    './components/**/*.{js,ts,jsx,tsx,mdx}',
  ],
  theme: {
    extend: {
      colors: {
        primary: {
          50: '#f0f9ff',
          100: '#e0f2fe',
          500: '#0ea5e9',
          600: '#0284c7',
          700: '#0369a1',
        },
        disease: {
          cancer: '#ef4444',
          neurological: '#8b5cf6',
          metabolic: '#f59e0b',
          aging: '#6366f1',
          inflammatory: '#ec4899',
          other: '#6b7280',
        }
      }
    },
  },
  plugins: [],
}
