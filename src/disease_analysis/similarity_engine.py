#!/usr/bin/env python3
"""
Disease Similarity Engine
=========================

ChromaDB-based vector similarity search for disease signatures.
Supports multi-modal matching:
1. Embedding similarity (semantic)
2. Gene overlap (Jaccard)
3. Pathway overlap
4. Direction-aware matching (up/down regulated)
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import warnings
warnings.filterwarnings('ignore')

import chromadb
from chromadb.config import Settings
from sentence_transformers import SentenceTransformer


@dataclass
class SimilarityResult:
    """Result of similarity search"""
    signature_id: str
    signature_name: str
    category: str
    stage: str
    severity_score: float
    embedding_similarity: float
    gene_overlap: int
    jaccard_index: float
    direction_score: float  # Agreement in up/down regulation
    combined_score: float
    overlapping_genes: List[str]
    description: str


class DiseaseSimilarityEngine:
    """ChromaDB-based disease similarity search engine"""

    def __init__(self, db_path: str = "data/chromadb_disease",
                 model_name: str = "all-MiniLM-L6-v2"):
        self.db_path = Path(db_path)
        self.db_path.mkdir(parents=True, exist_ok=True)

        print("  Loading embedding model...")
        self.embedder = SentenceTransformer(model_name)
        self.embedding_dim = self.embedder.get_sentence_embedding_dimension()

        print(f"  Initializing ChromaDB at {self.db_path}")
        self.client = chromadb.PersistentClient(
            path=str(self.db_path),
            settings=Settings(anonymized_telemetry=False)
        )

        self.collection = None
        self.signatures_data = None

    def build_index(self, signatures: Dict) -> bool:
        """Build ChromaDB index from signatures"""
        print("\n  Building ChromaDB index...")

        # Delete existing collection
        try:
            self.client.delete_collection("disease_signatures")
        except:
            pass

        self.collection = self.client.create_collection(
            name="disease_signatures",
            metadata={"hnsw:space": "cosine"}
        )

        ids = []
        embeddings = []
        metadatas = []
        documents = []

        for sig_id, sig in signatures.items():
            # Get genes
            genes = sig.genes if hasattr(sig, 'genes') else sig.get('genes', [])
            if not genes:
                continue

            # Handle nested lists
            if genes and isinstance(genes[0], list):
                genes = [g[0] if isinstance(g, list) else g for g in genes]

            # Create text for embedding
            name = sig.name if hasattr(sig, 'name') else sig.get('name', sig_id)
            if isinstance(name, list):
                name = name[0] if name else sig_id

            text = f"{name}: {' '.join(str(g) for g in genes[:100])}"

            # Generate embedding
            embedding = self.embedder.encode(text, show_progress_bar=False)

            # Prepare metadata
            category = sig.category if hasattr(sig, 'category') else sig.get('category', 'unknown')
            stage = sig.stage if hasattr(sig, 'stage') else sig.get('stage', 'unknown')
            severity = sig.severity_score if hasattr(sig, 'severity_score') else sig.get('severity_score', 0.5)
            description = sig.description if hasattr(sig, 'description') else sig.get('description', '')

            if isinstance(category, list):
                category = category[0] if category else 'unknown'
            if isinstance(stage, list):
                stage = stage[0] if stage else 'unknown'

            ids.append(sig_id)
            embeddings.append(embedding.tolist())
            metadatas.append({
                'name': str(name),
                'category': str(category),
                'stage': str(stage),
                'severity_score': float(severity),
                'n_genes': len(genes),
                'genes_preview': ','.join(str(g) for g in genes[:30]),
                'description': str(description)[:500]
            })
            documents.append(text)

        # Add in batches
        batch_size = 500
        for i in range(0, len(ids), batch_size):
            end = min(i + batch_size, len(ids))
            self.collection.add(
                ids=ids[i:end],
                embeddings=embeddings[i:end],
                metadatas=metadatas[i:end],
                documents=documents[i:end]
            )

        print(f"    Indexed {len(ids)} signatures")

        # Store signatures for gene overlap calculation
        self.signatures_data = {}
        for sig_id, sig in signatures.items():
            genes = sig.genes if hasattr(sig, 'genes') else sig.get('genes', [])
            if genes and isinstance(genes[0], list):
                genes = [g[0] if isinstance(g, list) else g for g in genes]

            up_genes = sig.up_genes if hasattr(sig, 'up_genes') else sig.get('up_genes', [])
            down_genes = sig.down_genes if hasattr(sig, 'down_genes') else sig.get('down_genes', [])

            self.signatures_data[sig_id] = {
                'genes': set(genes),
                'up_genes': set(up_genes) if up_genes else set(),
                'down_genes': set(down_genes) if down_genes else set()
            }

        # Save signatures data
        sig_file = self.db_path / "signatures_data.json"
        serializable = {k: {'genes': list(v['genes']),
                           'up_genes': list(v['up_genes']),
                           'down_genes': list(v['down_genes'])}
                       for k, v in self.signatures_data.items()}
        with open(sig_file, 'w') as f:
            json.dump(serializable, f)

        return True

    def load_index(self) -> bool:
        """Load existing index"""
        try:
            self.collection = self.client.get_collection("disease_signatures")

            # Load signatures data
            sig_file = self.db_path / "signatures_data.json"
            if sig_file.exists():
                with open(sig_file) as f:
                    data = json.load(f)
                self.signatures_data = {k: {'genes': set(v['genes']),
                                           'up_genes': set(v['up_genes']),
                                           'down_genes': set(v['down_genes'])}
                                       for k, v in data.items()}

            print(f"    Loaded index with {self.collection.count()} signatures")
            return True
        except Exception as e:
            print(f"    Failed to load index: {e}")
            return False

    def search(self,
               query_genes: List[str],
               query_up_genes: Optional[List[str]] = None,
               query_down_genes: Optional[List[str]] = None,
               fold_changes: Optional[Dict[str, float]] = None,
               top_k: int = 50,
               category_filter: Optional[str] = None) -> List[SimilarityResult]:
        """Search for similar disease signatures"""

        if self.collection is None:
            raise ValueError("Index not loaded. Call build_index or load_index first.")

        # Prepare query
        query_set = set(query_genes)

        # Determine up/down genes from fold changes if provided
        if fold_changes and not query_up_genes:
            query_up_genes = [g for g in query_genes if fold_changes.get(g, 0) > 0]
            query_down_genes = [g for g in query_genes if fold_changes.get(g, 0) < 0]

        query_up_set = set(query_up_genes) if query_up_genes else set()
        query_down_set = set(query_down_genes) if query_down_genes else set()

        # Create query embedding
        if query_up_genes or query_down_genes:
            query_text = f"upregulated: {' '.join(query_up_genes[:50])} | downregulated: {' '.join(query_down_genes[:50])}"
        else:
            query_text = " ".join(query_genes[:100])

        query_emb = self.embedder.encode(query_text, show_progress_bar=False)

        # Query ChromaDB
        where_filter = {"category": category_filter} if category_filter else None

        results = self.collection.query(
            query_embeddings=[query_emb.tolist()],
            n_results=top_k,
            include=['metadatas', 'distances', 'documents'],
            where=where_filter
        )

        # Process results
        output = []

        for i in range(len(results['ids'][0])):
            sig_id = results['ids'][0][i]
            metadata = results['metadatas'][0][i]
            embedding_sim = 1 - results['distances'][0][i]

            # Calculate gene overlap
            if self.signatures_data and sig_id in self.signatures_data:
                sig_genes = self.signatures_data[sig_id]['genes']
                sig_up = self.signatures_data[sig_id]['up_genes']
                sig_down = self.signatures_data[sig_id]['down_genes']

                overlap = query_set.intersection(sig_genes)
                jaccard = len(overlap) / len(query_set.union(sig_genes)) if sig_genes else 0

                # Direction score
                direction_score = 0
                if query_up_set and sig_up:
                    up_agree = len(query_up_set.intersection(sig_up))
                    up_disagree = len(query_up_set.intersection(sig_down))
                    direction_score += (up_agree - up_disagree) / max(len(query_up_set), 1)

                if query_down_set and sig_down:
                    down_agree = len(query_down_set.intersection(sig_down))
                    down_disagree = len(query_down_set.intersection(sig_up))
                    direction_score += (down_agree - down_disagree) / max(len(query_down_set), 1)

                direction_score = max(0, min(1, (direction_score + 1) / 2))  # Normalize to 0-1
            else:
                overlap = set()
                jaccard = 0
                direction_score = 0.5

            # Combined score (weighted average)
            combined = (
                embedding_sim * 0.4 +
                jaccard * 0.3 +
                direction_score * 0.3
            )

            result = SimilarityResult(
                signature_id=sig_id,
                signature_name=metadata['name'],
                category=metadata['category'],
                stage=metadata['stage'],
                severity_score=metadata['severity_score'],
                embedding_similarity=embedding_sim,
                gene_overlap=len(overlap),
                jaccard_index=jaccard,
                direction_score=direction_score,
                combined_score=combined,
                overlapping_genes=sorted(list(overlap))[:20],
                description=metadata.get('description', '')
            )
            output.append(result)

        # Sort by combined score
        output.sort(key=lambda x: x.combined_score, reverse=True)

        return output

    def search_by_category(self,
                           query_genes: List[str],
                           fold_changes: Optional[Dict[str, float]] = None,
                           top_k_per_category: int = 10) -> Dict[str, List[SimilarityResult]]:
        """Search for similar signatures grouped by category"""

        categories = ['cancer', 'neurological', 'metabolic', 'inflammatory',
                     'aging', 'other']

        results_by_category = {}

        for cat in categories:
            try:
                results = self.search(
                    query_genes=query_genes,
                    fold_changes=fold_changes,
                    top_k=top_k_per_category,
                    category_filter=cat
                )
                if results:
                    results_by_category[cat] = results
            except:
                continue

        return results_by_category


if __name__ == "__main__":
    # Test
    engine = DiseaseSimilarityEngine()

    if not engine.load_index():
        from disease_database import DiseaseDatabaseBuilder
        builder = DiseaseDatabaseBuilder()
        sigs = builder.load_database()
        engine.build_index(sigs)

    # Test search
    test_genes = ['BHLHE40', 'DDIT4', 'VEGFA', 'SLC2A3', 'BNIP3']
    results = engine.search(test_genes, top_k=10)

    for r in results:
        print(f"{r.signature_name}: {r.combined_score:.3f} (overlap: {r.gene_overlap})")
