#!/usr/bin/env python3
"""
Disease Similarity Analysis using MSigDB and ChromaDB
======================================================

Enhanced version using MSigDB disease-related gene sets for better
gene symbol matching and coverage.

MSigDB Collections used:
- C7: Immunologic signatures (disease-related immune states)
- C2:CGP: Chemical and genetic perturbations (includes disease signatures)
- HALLMARK: Cancer hallmarks and biological states
"""

import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Set
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

import chromadb
from chromadb.config import Settings
from sentence_transformers import SentenceTransformer


class MSigDBDiseaseLoader:
    """Load disease-related gene sets from MSigDB via msigdbr"""

    def __init__(self, data_dir: str = "data/msigdb"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)

    def load_disease_gene_sets(self) -> Dict[str, Dict]:
        """Load disease-related gene sets from MSigDB"""

        print("\n  Loading MSigDB disease-related gene sets...")

        cache_file = self.data_dir / "disease_gene_sets.json"

        if cache_file.exists():
            print(f"    Loading from cache: {cache_file}")
            with open(cache_file) as f:
                raw_data = json.load(f)

            # Convert gene lists if they are nested
            processed = {}
            for gs_name, gs_data in raw_data.items():
                genes = gs_data.get('genes', [])
                # Handle nested lists from R's toJSON
                if genes and isinstance(genes[0], list):
                    genes = [g[0] if isinstance(g, list) else g for g in genes]
                processed[gs_name] = {
                    'name': gs_data.get('name', gs_name),
                    'genes': genes,
                    'category': gs_data.get('category', 'unknown'),
                    'n_genes': len(genes)
                }

            print(f"    Loaded {len(processed)} gene sets")
            return processed

        # Use R to fetch MSigDB data (msigdbr has better coverage)
        gene_sets = self._fetch_via_r()

        # Save cache
        with open(cache_file, 'w') as f:
            json.dump(gene_sets, f)

        print(f"    Saved {len(gene_sets)} gene sets to cache")
        return gene_sets

    def _fetch_via_r(self) -> Dict[str, Dict]:
        """Fetch MSigDB gene sets via R msigdbr package"""
        import subprocess
        import tempfile

        print("    Fetching MSigDB gene sets via R...")

        r_script = '''
        library(msigdbr)
        library(jsonlite)

        # Get disease-related gene sets
        collections <- list(
            # Hallmark gene sets (cancer hallmarks, biological states)
            list(collection = "H", subcollection = NULL, category = "hallmark"),
            # C2:CGP - Chemical and genetic perturbations (disease signatures)
            list(collection = "C2", subcollection = "CGP", category = "perturbation"),
            # C7 - Immunologic signatures
            list(collection = "C7", subcollection = NULL, category = "immunologic")
        )

        all_gene_sets <- list()

        for (coll in collections) {
            cat(sprintf("  Fetching %s...", coll$category))

            tryCatch({
                if (is.null(coll$subcollection)) {
                    msig <- msigdbr(species = "Homo sapiens", collection = coll$collection)
                } else {
                    msig <- msigdbr(species = "Homo sapiens",
                                   collection = coll$collection,
                                   subcollection = coll$subcollection)
                }

                # Group by gene set name
                for (gs_name in unique(msig$gs_name)) {
                    genes <- msig$gene_symbol[msig$gs_name == gs_name]
                    all_gene_sets[[gs_name]] <- list(
                        name = gs_name,
                        genes = unique(genes),
                        category = coll$category,
                        n_genes = length(unique(genes))
                    )
                }
                cat(sprintf(" %d gene sets\\n", length(unique(msig$gs_name))))
            }, error = function(e) {
                cat(sprintf(" Error: %s\\n", e$message))
            })
        }

        # Convert to JSON
        json_output <- toJSON(all_gene_sets, auto_unbox = TRUE)
        cat(json_output)
        '''

        with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
            f.write(r_script)
            r_file = f.name

        try:
            result = subprocess.run(
                ['Rscript', r_file],
                capture_output=True,
                text=True,
                timeout=300
            )

            # Parse JSON from output (after the status messages)
            output = result.stdout
            json_start = output.rfind('{')
            if json_start == -1:
                json_start = output.rfind('[')

            if json_start >= 0:
                json_str = output[json_start:]
                gene_sets = json.loads(json_str)
                print(f"    Retrieved {len(gene_sets)} gene sets")
                return gene_sets
            else:
                print(f"    Warning: Could not parse R output")
                return {}

        except Exception as e:
            print(f"    Error: {e}")
            return {}
        finally:
            os.unlink(r_file)


class EnhancedDiseaseAnalyzer:
    """Enhanced disease similarity analyzer using MSigDB"""

    def __init__(self, data_dir: str = "data", db_dir: str = "data/chromadb_msigdb"):
        self.data_dir = Path(data_dir)
        self.db_dir = Path(db_dir)
        self.db_dir.mkdir(parents=True, exist_ok=True)

        print("\n=== Initializing Enhanced Disease Similarity Analyzer ===")

        # Load embedding model
        print("  Loading embedding model...")
        self.embedder = SentenceTransformer("all-MiniLM-L6-v2")

        # Initialize ChromaDB
        print(f"  Initializing ChromaDB at {self.db_dir}")
        self.client = chromadb.PersistentClient(
            path=str(self.db_dir),
            settings=Settings(anonymized_telemetry=False)
        )

        self.collection = None
        self.gene_sets = None

    def build_database(self):
        """Build the disease gene set database"""

        print("\n=== Building Disease Gene Set Database ===")

        # Load MSigDB gene sets
        loader = MSigDBDiseaseLoader(str(self.data_dir / "msigdb"))
        self.gene_sets = loader.load_disease_gene_sets()

        if not self.gene_sets:
            print("  No gene sets loaded. Cannot build database.")
            return False

        # Filter gene sets by size
        min_genes, max_genes = 10, 500
        filtered_sets = {
            k: v for k, v in self.gene_sets.items()
            if min_genes <= v.get('n_genes', 0) <= max_genes
        }
        print(f"  Filtered to {len(filtered_sets)} gene sets (size: {min_genes}-{max_genes})")

        # Create embeddings
        print("  Creating embeddings...")
        embeddings = {}
        for gs_name, gs_data in filtered_sets.items():
            text = f"{gs_name}: {' '.join(gs_data['genes'][:100])}"
            embeddings[gs_name] = self.embedder.encode(text, show_progress_bar=False)

        # Create ChromaDB collection
        print("  Creating ChromaDB collection...")
        try:
            self.client.delete_collection("msigdb_disease")
        except:
            pass

        self.collection = self.client.create_collection(
            name="msigdb_disease",
            metadata={"hnsw:space": "cosine"}
        )

        # Add to collection
        ids = list(filtered_sets.keys())
        embs = [embeddings[k].tolist() for k in ids]
        metadatas = []
        for v in filtered_sets.values():
            # Ensure all metadata values are strings or numbers, not lists
            name = v['name']
            if isinstance(name, list):
                name = name[0] if name else ''

            category = v.get('category', 'unknown')
            if isinstance(category, list):
                category = category[0] if category else 'unknown'

            genes = v['genes']
            genes_preview = ','.join(str(g) for g in genes[:30])

            metadatas.append({
                'name': str(name),
                'category': str(category),
                'n_genes': int(v['n_genes']),
                'genes_preview': genes_preview
            })
        documents = []
        for v in filtered_sets.values():
            name = v['name']
            if isinstance(name, list):
                name = name[0] if name else ''
            genes = v['genes'][:50]
            genes_str = ' '.join(str(g) for g in genes)
            documents.append(f"{name}: {genes_str}")

        # Add in batches
        batch_size = 500
        for i in range(0, len(ids), batch_size):
            end = min(i + batch_size, len(ids))
            self.collection.add(
                ids=ids[i:end],
                embeddings=embs[i:end],
                metadatas=metadatas[i:end],
                documents=documents[i:end]
            )

        print(f"  Added {len(ids)} gene sets to ChromaDB")

        # Save gene sets for reference
        gs_file = self.data_dir / "msigdb_gene_sets.json"
        with open(gs_file, 'w') as f:
            json.dump(filtered_sets, f)

        self.gene_sets = filtered_sets
        return True

    def load_database(self) -> bool:
        """Load existing database"""
        try:
            self.collection = self.client.get_collection("msigdb_disease")

            # Load gene sets
            gs_file = self.data_dir / "msigdb_gene_sets.json"
            if gs_file.exists():
                with open(gs_file) as f:
                    self.gene_sets = json.load(f)

            print(f"  Loaded collection with {self.collection.count()} gene sets")
            return True
        except:
            return False

    def find_similar_diseases(self,
                              genes: List[str],
                              fold_changes: Optional[Dict[str, float]] = None,
                              top_k: int = 30) -> pd.DataFrame:
        """Find gene sets similar to query"""

        print(f"\n=== Finding Similar Gene Sets ===")
        print(f"  Query genes: {len(genes)}")

        if self.collection is None:
            if not self.load_database():
                print("  Database not found. Building...")
                self.build_database()

        # Create query embedding
        if fold_changes:
            up_genes = [g for g in genes if fold_changes.get(g, 0) > 0]
            down_genes = [g for g in genes if fold_changes.get(g, 0) < 0]
            query_text = f"upregulated: {' '.join(up_genes[:50])} | downregulated: {' '.join(down_genes[:50])}"
        else:
            query_text = " ".join(genes[:100])

        query_emb = self.embedder.encode(query_text, show_progress_bar=False)

        # Query ChromaDB
        results = self.collection.query(
            query_embeddings=[query_emb.tolist()],
            n_results=top_k,
            include=['metadatas', 'distances', 'documents']
        )

        # Process results
        output = []
        query_set = set(genes)

        for i in range(len(results['ids'][0])):
            gs_id = results['ids'][0][i]
            metadata = results['metadatas'][0][i]
            similarity = 1 - results['distances'][0][i]

            # Calculate gene overlap
            if self.gene_sets and gs_id in self.gene_sets:
                gs_genes = set(self.gene_sets[gs_id]['genes'])
                overlap = query_set.intersection(gs_genes)
                overlap_count = len(overlap)
                overlap_genes = ','.join(sorted(overlap)[:20])
                jaccard = len(overlap) / len(query_set.union(gs_genes)) if gs_genes else 0
            else:
                overlap_count = 0
                overlap_genes = ''
                jaccard = 0

            output.append({
                'gene_set_id': gs_id,
                'gene_set_name': metadata['name'],
                'category': metadata['category'],
                'n_genes': metadata['n_genes'],
                'embedding_similarity': similarity,
                'gene_overlap': overlap_count,
                'jaccard_index': jaccard,
                'overlapping_genes': overlap_genes,
                'combined_score': (similarity * 0.5) + (jaccard * 0.5)
            })

        df = pd.DataFrame(output)
        df = df.sort_values('combined_score', ascending=False)

        return df

    def analyze_de_results(self, de_file: str, output_dir: str,
                           padj_threshold: float = 0.05,
                           top_k: int = 30) -> pd.DataFrame:
        """Analyze DE results"""

        print(f"\n=== Analyzing DE Results ===")
        print(f"  File: {de_file}")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Load DE results
        de_df = pd.read_csv(de_file)

        # Get genes
        if 'padj' in de_df.columns:
            sig_genes = de_df[de_df['padj'] < padj_threshold]['gene'].tolist()
        else:
            sig_genes = de_df['gene'].tolist()

        print(f"  Significant genes (padj < {padj_threshold}): {len(sig_genes)}")

        # Get fold changes
        fold_changes = {}
        if 'log2FoldChange' in de_df.columns:
            for _, row in de_df.iterrows():
                fold_changes[row['gene']] = row['log2FoldChange']

        # If few significant genes, use top by fold change
        if len(sig_genes) < 20:
            print("  Using top genes by |log2FC| due to few significant genes...")
            de_df['abs_lfc'] = de_df['log2FoldChange'].abs()
            top_genes = de_df.nlargest(200, 'abs_lfc')['gene'].tolist()
            query_genes = list(set(sig_genes + top_genes))
        else:
            query_genes = sig_genes

        print(f"  Query genes: {len(query_genes)}")

        # Find similar gene sets
        results = self.find_similar_diseases(
            genes=query_genes,
            fold_changes=fold_changes,
            top_k=top_k
        )

        # Save results
        output_file = output_path / "disease_similarity_msigdb.csv"
        results.to_csv(output_file, index=False)

        # Print results
        print(f"\n  Top Similar Gene Sets:")
        print("  " + "=" * 80)

        for i, row in results.head(15).iterrows():
            print(f"\n  {results.index.get_loc(i)+1}. {row['gene_set_name'][:60]}")
            print(f"     Category: {row['category']}")
            print(f"     Embedding similarity: {row['embedding_similarity']:.4f}")
            print(f"     Gene overlap: {row['gene_overlap']} genes (Jaccard: {row['jaccard_index']:.4f})")
            print(f"     Combined score: {row['combined_score']:.4f}")
            if row['overlapping_genes']:
                print(f"     Overlapping: {row['overlapping_genes'][:70]}...")

        print(f"\n  Results saved to: {output_file}")
        return results


def main():
    import argparse

    parser = argparse.ArgumentParser(description='MSigDB Disease Similarity Analysis')
    parser.add_argument('--de_genes', type=str, required=True)
    parser.add_argument('--output', type=str, default='results/disease_similarity')
    parser.add_argument('--data_dir', type=str, default='data')
    parser.add_argument('--rebuild', action='store_true')
    parser.add_argument('--top_k', type=int, default=30)

    args = parser.parse_args()

    analyzer = EnhancedDiseaseAnalyzer(
        data_dir=args.data_dir,
        db_dir=os.path.join(args.data_dir, "chromadb_msigdb")
    )

    if args.rebuild or not analyzer.load_database():
        analyzer.build_database()

    results = analyzer.analyze_de_results(
        de_file=args.de_genes,
        output_dir=args.output,
        top_k=args.top_k
    )

    return results


if __name__ == "__main__":
    main()
