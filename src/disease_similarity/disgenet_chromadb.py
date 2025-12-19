#!/usr/bin/env python3
"""
Disease Similarity Analysis using DisGeNET and ChromaDB
========================================================

This module performs disease similarity analysis by:
1. Downloading disease-gene associations from DisGeNET
2. Creating embeddings for disease gene signatures
3. Storing disease signatures in ChromaDB vector store
4. Querying with user's DE genes to find similar diseases

Usage:
    python disgenet_chromadb.py --de_genes results/deseq2/significant_genes.csv --output results/disease_similarity
"""

import os
import sys
import json
import hashlib
import requests
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# ChromaDB and embeddings
import chromadb
from chromadb.config import Settings
from sentence_transformers import SentenceTransformer


class DisGeNETDownloader:
    """Download and parse DisGeNET disease-gene associations"""

    # DisGeNET public data URLs (curated subset)
    DISGENET_CURATED_URL = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz"
    DISGENET_ALL_URL = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz"

    def __init__(self, data_dir: str = "data/disgenet"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)

    def download_disgenet(self, use_curated: bool = True) -> pd.DataFrame:
        """Download DisGeNET data (or load from cache)"""

        cache_file = self.data_dir / ("curated_gda.parquet" if use_curated else "all_gda.parquet")

        if cache_file.exists():
            print(f"  Loading cached DisGeNET data from {cache_file}")
            return pd.read_parquet(cache_file)

        # DisGeNET requires API key for direct download
        # Use alternative: fetch from GitHub mirror or create from public sources
        print("  DisGeNET requires registration for bulk download.")
        print("  Using alternative approach with curated disease-gene data...")

        return self._create_disease_gene_db()

    def _create_disease_gene_db(self) -> pd.DataFrame:
        """Create disease-gene database from multiple public sources"""

        print("\n  Building disease-gene database from public sources...")

        # 1. OMIM disease genes (from NCBI)
        omim_data = self._fetch_omim_genes()

        # 2. KEGG disease genes
        kegg_data = self._fetch_kegg_disease_genes()

        # 3. Combine all sources
        all_data = pd.concat([omim_data, kegg_data], ignore_index=True)

        # Remove duplicates
        all_data = all_data.drop_duplicates(subset=['disease_id', 'gene_symbol'])

        # Save cache
        cache_file = self.data_dir / "disease_gene_db.parquet"
        all_data.to_parquet(cache_file)
        print(f"    Saved {len(all_data)} disease-gene associations to {cache_file}")

        return all_data

    def _fetch_omim_genes(self) -> pd.DataFrame:
        """Fetch OMIM morbid genes from NCBI"""
        print("    Fetching OMIM morbid genes...")

        # Use NCBI's mim2gene_medgen file
        url = "https://ftp.ncbi.nih.gov/gene/DATA/mim2gene_medgen"

        try:
            df = pd.read_csv(url, sep='\t', comment='#',
                           names=['mim_number', 'gene_id', 'type', 'gene_symbol', 'medgen_cui'])

            # Filter for phenotype entries with gene associations
            df = df[df['type'].isin(['phenotype', 'gene/phenotype'])]
            df = df[df['gene_symbol'].notna() & (df['gene_symbol'] != '-')]

            result = pd.DataFrame({
                'disease_id': 'OMIM:' + df['mim_number'].astype(str),
                'disease_name': 'OMIM_' + df['mim_number'].astype(str),  # Will be updated
                'gene_symbol': df['gene_symbol'],
                'gene_id': df['gene_id'],
                'source': 'OMIM'
            })

            print(f"      Retrieved {len(result)} OMIM associations")
            return result

        except Exception as e:
            print(f"      Warning: Could not fetch OMIM data: {e}")
            return pd.DataFrame(columns=['disease_id', 'disease_name', 'gene_symbol', 'gene_id', 'source'])

    def _fetch_kegg_disease_genes(self) -> pd.DataFrame:
        """Fetch KEGG disease genes"""
        print("    Fetching KEGG disease genes...")

        try:
            # Get list of human diseases
            diseases_url = "https://rest.kegg.jp/list/disease"
            diseases_resp = requests.get(diseases_url)

            results = []
            disease_lines = diseases_resp.text.strip().split('\n')[:200]  # Limit for speed

            for line in disease_lines:
                parts = line.split('\t')
                if len(parts) >= 2:
                    disease_id = parts[0]  # e.g., ds:H00001
                    disease_name = parts[1]

                    # Get genes for this disease
                    try:
                        gene_url = f"https://rest.kegg.jp/link/hsa/{disease_id}"
                        gene_resp = requests.get(gene_url, timeout=5)

                        for gene_line in gene_resp.text.strip().split('\n'):
                            if gene_line:
                                gene_parts = gene_line.split('\t')
                                if len(gene_parts) >= 2:
                                    gene_id = gene_parts[1].replace('hsa:', '')
                                    results.append({
                                        'disease_id': disease_id,
                                        'disease_name': disease_name,
                                        'gene_symbol': gene_id,
                                        'gene_id': gene_id,
                                        'source': 'KEGG'
                                    })
                    except:
                        continue

            df = pd.DataFrame(results)
            print(f"      Retrieved {len(df)} KEGG associations")
            return df

        except Exception as e:
            print(f"      Warning: Could not fetch KEGG data: {e}")
            return pd.DataFrame(columns=['disease_id', 'disease_name', 'gene_symbol', 'gene_id', 'source'])


class DiseaseSignatureEmbedder:
    """Create embeddings for disease gene signatures"""

    def __init__(self, model_name: str = "all-MiniLM-L6-v2"):
        print(f"\n  Loading embedding model: {model_name}")
        self.model = SentenceTransformer(model_name)
        self.embedding_dim = self.model.get_sentence_embedding_dimension()
        print(f"    Embedding dimension: {self.embedding_dim}")

    def create_disease_signatures(self, gda_df: pd.DataFrame,
                                   min_genes: int = 5,
                                   max_genes: int = 500) -> Dict[str, Dict]:
        """Create disease signatures (gene sets) from GDA data"""

        print("\n  Creating disease signatures...")

        # Group genes by disease
        disease_genes = defaultdict(set)
        disease_names = {}

        for _, row in gda_df.iterrows():
            disease_id = row['disease_id']
            gene = row['gene_symbol']

            disease_genes[disease_id].add(gene)
            if 'disease_name' in row:
                disease_names[disease_id] = row['disease_name']

        # Filter by gene count
        signatures = {}
        for disease_id, genes in disease_genes.items():
            if min_genes <= len(genes) <= max_genes:
                signatures[disease_id] = {
                    'genes': sorted(list(genes)),
                    'name': disease_names.get(disease_id, disease_id),
                    'n_genes': len(genes)
                }

        print(f"    Created {len(signatures)} disease signatures (genes: {min_genes}-{max_genes})")
        return signatures

    def embed_disease_signatures(self, signatures: Dict[str, Dict]) -> Dict[str, np.ndarray]:
        """Create embeddings for disease signatures"""

        print("\n  Embedding disease signatures...")

        embeddings = {}

        for disease_id, sig in signatures.items():
            # Create text representation of gene signature
            gene_text = " ".join(sig['genes'])
            disease_text = f"{sig['name']}: {gene_text}"

            # Generate embedding
            embedding = self.model.encode(disease_text, show_progress_bar=False)
            embeddings[disease_id] = embedding

        print(f"    Embedded {len(embeddings)} disease signatures")
        return embeddings

    def embed_query_genes(self, genes: List[str],
                          fold_changes: Optional[Dict[str, float]] = None) -> np.ndarray:
        """Create embedding for query gene set"""

        # Create text representation
        if fold_changes:
            # Include direction information
            up_genes = [g for g in genes if fold_changes.get(g, 0) > 0]
            down_genes = [g for g in genes if fold_changes.get(g, 0) < 0]

            text_parts = []
            if up_genes:
                text_parts.append(f"upregulated: {' '.join(up_genes)}")
            if down_genes:
                text_parts.append(f"downregulated: {' '.join(down_genes)}")

            query_text = " | ".join(text_parts) if text_parts else " ".join(genes)
        else:
            query_text = " ".join(genes)

        return self.model.encode(query_text, show_progress_bar=False)


class DiseaseSimilarityDB:
    """ChromaDB-based disease similarity database"""

    def __init__(self, db_path: str = "data/chromadb"):
        self.db_path = Path(db_path)
        self.db_path.mkdir(parents=True, exist_ok=True)

        print(f"\n  Initializing ChromaDB at {self.db_path}")

        self.client = chromadb.PersistentClient(
            path=str(self.db_path),
            settings=Settings(anonymized_telemetry=False)
        )

        self.collection = None
        self.embedder = None

    def create_collection(self, signatures: Dict[str, Dict],
                          embeddings: Dict[str, np.ndarray],
                          collection_name: str = "disease_signatures"):
        """Create or update ChromaDB collection"""

        print(f"\n  Creating ChromaDB collection: {collection_name}")

        # Delete existing collection if exists
        try:
            self.client.delete_collection(collection_name)
        except:
            pass

        # Create new collection
        self.collection = self.client.create_collection(
            name=collection_name,
            metadata={"hnsw:space": "cosine"}
        )

        # Add disease signatures
        ids = []
        embs = []
        metadatas = []
        documents = []

        for disease_id, sig in signatures.items():
            if disease_id in embeddings:
                ids.append(disease_id)
                embs.append(embeddings[disease_id].tolist())
                metadatas.append({
                    'name': sig['name'],
                    'n_genes': sig['n_genes'],
                    'genes': ','.join(sig['genes'][:50])  # Limit for metadata
                })
                documents.append(f"{sig['name']}: {' '.join(sig['genes'])}")

        # Add in batches
        batch_size = 1000
        for i in range(0, len(ids), batch_size):
            end_idx = min(i + batch_size, len(ids))
            self.collection.add(
                ids=ids[i:end_idx],
                embeddings=embs[i:end_idx],
                metadatas=metadatas[i:end_idx],
                documents=documents[i:end_idx]
            )

        print(f"    Added {len(ids)} disease signatures to ChromaDB")

    def load_collection(self, collection_name: str = "disease_signatures"):
        """Load existing collection"""
        try:
            self.collection = self.client.get_collection(collection_name)
            print(f"    Loaded collection with {self.collection.count()} items")
            return True
        except:
            print(f"    Collection {collection_name} not found")
            return False

    def query_similar_diseases(self, query_embedding: np.ndarray,
                                top_k: int = 20) -> List[Dict]:
        """Query for similar diseases"""

        if self.collection is None:
            raise ValueError("No collection loaded. Call create_collection or load_collection first.")

        results = self.collection.query(
            query_embeddings=[query_embedding.tolist()],
            n_results=top_k,
            include=['metadatas', 'distances', 'documents']
        )

        similar_diseases = []
        for i in range(len(results['ids'][0])):
            similar_diseases.append({
                'disease_id': results['ids'][0][i],
                'disease_name': results['metadatas'][0][i]['name'],
                'n_genes': results['metadatas'][0][i]['n_genes'],
                'similarity': 1 - results['distances'][0][i],  # Convert distance to similarity
                'genes_preview': results['metadatas'][0][i]['genes']
            })

        return similar_diseases


class DiseaseSimilarityAnalyzer:
    """Main class for disease similarity analysis"""

    def __init__(self, data_dir: str = "data", db_dir: str = "data/chromadb"):
        self.data_dir = Path(data_dir)
        self.db_dir = Path(db_dir)

        self.downloader = DisGeNETDownloader(str(self.data_dir / "disgenet"))
        self.embedder = DiseaseSignatureEmbedder()
        self.db = DiseaseSimilarityDB(str(self.db_dir))

        self.signatures = None
        self.gda_df = None

    def build_database(self, use_curated: bool = True):
        """Build the disease similarity database"""

        print("\n=== Building Disease Similarity Database ===")

        # 1. Download/load disease-gene associations
        self.gda_df = self.downloader.download_disgenet(use_curated)

        # 2. Create disease signatures
        self.signatures = self.embedder.create_disease_signatures(
            self.gda_df, min_genes=5, max_genes=500
        )

        # 3. Create embeddings
        embeddings = self.embedder.embed_disease_signatures(self.signatures)

        # 4. Store in ChromaDB
        self.db.create_collection(self.signatures, embeddings)

        # 5. Save signatures for reference
        sig_file = self.data_dir / "disease_signatures.json"
        with open(sig_file, 'w') as f:
            # Convert to JSON-serializable format
            json_sigs = {k: {'name': v['name'], 'n_genes': v['n_genes'], 'genes': v['genes']}
                        for k, v in self.signatures.items()}
            json.dump(json_sigs, f, indent=2)
        print(f"\n  Saved signatures to {sig_file}")

        print("\n=== Database Build Complete ===")

    def find_similar_diseases(self,
                              genes: List[str],
                              fold_changes: Optional[Dict[str, float]] = None,
                              top_k: int = 20) -> pd.DataFrame:
        """Find diseases similar to the query gene signature"""

        print(f"\n=== Finding Similar Diseases ===")
        print(f"  Query genes: {len(genes)}")

        # Load collection if not already loaded
        if self.db.collection is None:
            if not self.db.load_collection():
                raise ValueError("Database not built. Run build_database() first.")

        # Create query embedding
        query_emb = self.embedder.embed_query_genes(genes, fold_changes)

        # Query similar diseases
        results = self.db.query_similar_diseases(query_emb, top_k)

        # Convert to DataFrame
        df = pd.DataFrame(results)

        # Calculate gene overlap
        if self.signatures is None:
            # Load signatures
            sig_file = self.data_dir / "disease_signatures.json"
            if sig_file.exists():
                with open(sig_file) as f:
                    self.signatures = json.load(f)

        if self.signatures:
            query_set = set(genes)
            overlaps = []
            overlap_genes_list = []

            for _, row in df.iterrows():
                disease_id = row['disease_id']
                if disease_id in self.signatures:
                    disease_genes = set(self.signatures[disease_id]['genes'])
                    overlap = query_set.intersection(disease_genes)
                    overlaps.append(len(overlap))
                    overlap_genes_list.append(','.join(sorted(overlap)))
                else:
                    overlaps.append(0)
                    overlap_genes_list.append('')

            df['gene_overlap'] = overlaps
            df['overlapping_genes'] = overlap_genes_list

        return df

    def analyze_de_results(self, de_results_file: str,
                           output_dir: str,
                           padj_threshold: float = 0.05,
                           top_k: int = 30) -> pd.DataFrame:
        """Analyze DE results for disease similarity"""

        print(f"\n=== Analyzing DE Results for Disease Similarity ===")
        print(f"  DE results file: {de_results_file}")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Load DE results
        de_df = pd.read_csv(de_results_file)

        # Get significant genes
        if 'padj' in de_df.columns:
            sig_genes = de_df[de_df['padj'] < padj_threshold]['gene'].tolist()
        else:
            sig_genes = de_df['gene'].tolist()

        print(f"  Significant genes: {len(sig_genes)}")

        # Get fold changes
        fold_changes = {}
        if 'log2FoldChange' in de_df.columns and 'gene' in de_df.columns:
            for _, row in de_df.iterrows():
                fold_changes[row['gene']] = row['log2FoldChange']

        # If very few significant genes, use top genes by fold change
        if len(sig_genes) < 10:
            print("  Few significant genes. Using top 100 genes by |log2FC|...")
            de_df['abs_lfc'] = de_df['log2FoldChange'].abs()
            top_genes = de_df.nlargest(100, 'abs_lfc')['gene'].tolist()
            genes_to_query = list(set(sig_genes + top_genes))
        else:
            genes_to_query = sig_genes

        print(f"  Query genes: {len(genes_to_query)}")

        # Find similar diseases
        results = self.find_similar_diseases(
            genes=genes_to_query,
            fold_changes=fold_changes,
            top_k=top_k
        )

        # Save results
        output_file = output_path / "disease_similarity_results.csv"
        results.to_csv(output_file, index=False)
        print(f"\n  Results saved to: {output_file}")

        # Print top results
        print(f"\n  Top {min(10, len(results))} Similar Diseases:")
        print("  " + "=" * 70)
        for i, row in results.head(10).iterrows():
            print(f"  {i+1}. {row['disease_name'][:50]}")
            print(f"     Similarity: {row['similarity']:.4f} | Gene overlap: {row['gene_overlap']}")
            if row['overlapping_genes']:
                print(f"     Overlapping: {row['overlapping_genes'][:80]}...")
            print()

        return results


def main():
    """Main function for disease similarity analysis"""
    import argparse

    parser = argparse.ArgumentParser(description='Disease Similarity Analysis')
    parser.add_argument('--de_genes', type=str, required=True,
                       help='Path to DE results CSV file')
    parser.add_argument('--output', type=str, default='results/disease_similarity',
                       help='Output directory')
    parser.add_argument('--data_dir', type=str, default='data',
                       help='Data directory for DisGeNET cache')
    parser.add_argument('--rebuild_db', action='store_true',
                       help='Force rebuild of disease database')
    parser.add_argument('--top_k', type=int, default=30,
                       help='Number of similar diseases to return')
    parser.add_argument('--padj', type=float, default=0.05,
                       help='Adjusted p-value threshold for significant genes')

    args = parser.parse_args()

    # Initialize analyzer
    analyzer = DiseaseSimilarityAnalyzer(
        data_dir=args.data_dir,
        db_dir=os.path.join(args.data_dir, "chromadb")
    )

    # Build database if needed
    db_exists = Path(args.data_dir) / "chromadb" / "chroma.sqlite3"
    if args.rebuild_db or not db_exists.exists():
        analyzer.build_database()

    # Analyze DE results
    results = analyzer.analyze_de_results(
        de_results_file=args.de_genes,
        output_dir=args.output,
        padj_threshold=args.padj,
        top_k=args.top_k
    )

    print("\n=== Analysis Complete ===")
    return results


if __name__ == "__main__":
    main()
