"""
Gene Synonym Service - Efficient lookup for gene name variations
Optimized with batch fetching and comprehensive caching
"""

import json
import os
import requests
import time
from typing import Dict, List, Set, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from src.api.subtiwiki_api import GENE_NAME_TO_ID, get_genes_in_category

class SynonymService:
    def __init__(self):
        self.synonym_cache = {}  # gene_name -> [synonyms]
        self.reverse_synonym_cache = {}  # synonym -> actual_gene_name
        self.cache_file = os.path.join('src', 'data', 'cache', 'gene_synonyms.json')
        self.comprehensive_index_file = os.path.join('src', 'data', 'cache', 'comprehensive_synonym_index.json')
        self.lock = threading.Lock()
        self.load_synonym_cache()
        self.load_comprehensive_index()
    
    def load_synonym_cache(self):
        """Load cached synonym data from file"""
        try:
            if os.path.exists(self.cache_file):
                with open(self.cache_file, 'r', encoding='utf-8') as f:
                    cache_data = json.load(f)
                    self.synonym_cache = cache_data.get('synonyms', {})
                    self.reverse_synonym_cache = cache_data.get('reverse_synonyms', {})
                    print(f"Loaded {len(self.synonym_cache)} genes with synonym data")
        except Exception as e:
            print(f"Error loading synonym cache: {e}")
    
    def load_comprehensive_index(self):
        """Load comprehensive synonym index if available"""
        try:
            if os.path.exists(self.comprehensive_index_file):
                with open(self.comprehensive_index_file, 'r', encoding='utf-8') as f:
                    index_data = json.load(f)
                    forward_index = index_data.get('forward_index', {})
                    reverse_index = index_data.get('reverse_index', {})
                    
                    # Merge with existing cache
                    with self.lock:
                        self.synonym_cache.update(forward_index)
                        self.reverse_synonym_cache.update(reverse_index)
                    
                    print(f"Loaded comprehensive index: {len(forward_index)} genes, {len(reverse_index)} synonyms")
        except Exception as e:
            print(f"Error loading comprehensive index: {e}")
    
    def save_synonym_cache(self):
        """Save synonym data to cache file"""
        try:
            os.makedirs(os.path.dirname(self.cache_file), exist_ok=True)
            cache_data = {
                'synonyms': self.synonym_cache,
                'reverse_synonyms': self.reverse_synonym_cache,
                'total_genes': len(self.synonym_cache)
            }
            with open(self.cache_file, 'w', encoding='utf-8') as f:
                json.dump(cache_data, f, indent=2, ensure_ascii=False)
        except Exception as e:
            print(f"Error saving synonym cache: {e}")
    
    def fetch_gene_synonyms_batch(self, gene_names: List[str], max_workers: int = 5) -> Dict[str, List[str]]:
        """Fetch synonyms for multiple genes in parallel"""
        results = {}
        
        def fetch_single_gene(gene_name: str) -> Tuple[str, List[str]]:
            """Fetch synonyms for a single gene"""
            try:
                base_url = "https://www.subtiwiki.uni-goettingen.de/v5/api"
                response = requests.get(f"{base_url}/gene/{gene_name}", timeout=10)
                
                if response.status_code == 200:
                    data = response.json()
                    synonyms = data.get('data', {}).get('synonyms', [])
                    return gene_name, synonyms
                else:
                    print(f"API Error for {gene_name}: {response.status_code}")
                    return gene_name, []
                    
            except Exception as e:
                print(f"Error fetching synonyms for {gene_name}: {e}")
                return gene_name, []
        
        # Use ThreadPoolExecutor for parallel requests
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all requests
            future_to_gene = {executor.submit(fetch_single_gene, gene): gene for gene in gene_names}
            
            # Collect results
            for future in as_completed(future_to_gene):
                gene_name, synonyms = future.result()
                results[gene_name] = synonyms
                
                # Update caches
                with self.lock:
                    self.synonym_cache[gene_name] = synonyms
                    for synonym in synonyms:
                        self.reverse_synonym_cache[synonym.lower()] = gene_name
        
        return results
    
    def fetch_gene_synonyms(self, gene_name: str) -> List[str]:
        """Fetch synonyms for a specific gene from SubtiWiki API (legacy method)"""
        if gene_name in self.synonym_cache:
            return self.synonym_cache[gene_name]
        
        # Use batch method for single gene
        results = self.fetch_gene_synonyms_batch([gene_name])
        return results.get(gene_name, [])
    
    def build_comprehensive_index(self, all_genes: List[str], max_workers: int = 5):
        """Build comprehensive synonym index for all genes"""
        print(f"Building comprehensive synonym index for {len(all_genes)} genes...")
        
        # Filter out genes already in cache
        uncached_genes = [gene for gene in all_genes if gene not in self.synonym_cache]
        
        if not uncached_genes:
            print("All genes already cached!")
            return
        
        print(f"Fetching synonyms for {len(uncached_genes)} uncached genes...")
        
        # Fetch in batches
        batch_size = 50
        total_batches = (len(uncached_genes) + batch_size - 1) // batch_size
        
        for i in range(0, len(uncached_genes), batch_size):
            batch = uncached_genes[i:i + batch_size]
            batch_num = i // batch_size + 1
            
            print(f"Processing batch {batch_num}/{total_batches} ({len(batch)} genes)...")
            
            # Fetch synonyms for this batch
            batch_results = self.fetch_gene_synonyms_batch(batch, max_workers)
            
            # Save progress after each batch
            self.save_synonym_cache()
            
            # Small delay to be respectful to the API
            time.sleep(0.5)
        
        # Save comprehensive index
        self.save_comprehensive_index()
        print("Comprehensive synonym index built successfully!")
    
    def save_comprehensive_index(self):
        """Save comprehensive synonym index"""
        try:
            os.makedirs(os.path.dirname(self.comprehensive_index_file), exist_ok=True)
            
            index_data = {
                'metadata': {
                    'total_genes': len(self.synonym_cache),
                    'genes_with_synonyms': sum(1 for syns in self.synonym_cache.values() if syns),
                    'total_synonym_entries': len(self.reverse_synonym_cache),
                    'build_timestamp': time.time(),
                    'version': '2.0'
                },
                'forward_index': self.synonym_cache,
                'reverse_index': self.reverse_synonym_cache
            }
            
            with open(self.comprehensive_index_file, 'w', encoding='utf-8') as f:
                json.dump(index_data, f, indent=2, ensure_ascii=False)
                
            print(f"Comprehensive index saved: {len(self.synonym_cache)} genes, {len(self.reverse_synonym_cache)} synonyms")
            
        except Exception as e:
            print(f"Error saving comprehensive index: {e}")
    
    def find_missing_genes_with_synonyms(self, selected_groups: List[int], uploaded_genes: Set[str]) -> Dict:
        """
        Optimized synonym lookup using comprehensive cache - NO API calls during user interaction
        
        Args:
            selected_groups: List of selected group IDs
            uploaded_genes: Set of gene names from uploaded file
        
        Returns:
            Dict with resolved genes and synonym mappings
        """
        print(f"Starting optimized synonym lookup for {len(selected_groups)} groups...")
        
        # Get all genes from selected groups
        expected_genes = set()
        for group_id in selected_groups:
            group_genes = get_genes_in_category(group_id)
            if group_genes:
                expected_genes.update(group_genes)
        
        print(f"Expected genes from groups: {len(expected_genes)}")
        
        # Find missing genes (genes in groups but not in uploaded data)
        uploaded_genes_lower = {gene.lower() for gene in uploaded_genes}
        missing_genes = []
        
        for gene in expected_genes:
            if gene.lower() not in uploaded_genes_lower:
                missing_genes.append(gene)
        
        print(f"Missing genes to check for synonyms: {len(missing_genes)}")
        
        if not missing_genes:
            return {
                'total_expected': len(expected_genes),
                'found_direct': len(expected_genes),
                'found_via_synonyms': 0,
                'still_missing': 0,
                'synonym_mappings': {},
                'missing_genes': []
            }
        
        # Check synonyms using comprehensive cache (NO API calls)
        synonym_mappings = {}
        found_via_synonyms = 0
        still_missing = []
        
        for gene in missing_genes:
            found_synonym = False
            
            # Check if gene has synonyms in cache
            if gene in self.synonym_cache:
                synonyms = self.synonym_cache[gene]
                
                # Check if any synonym matches uploaded genes
                for synonym in synonyms:
                    if synonym.lower() in uploaded_genes_lower:
                        # Find the exact case from uploaded genes
                        matched_uploaded = next(
                            (ug for ug in uploaded_genes if ug.lower() == synonym.lower()), 
                            synonym
                        )
                        synonym_mappings[matched_uploaded] = gene
                        found_via_synonyms += 1
                        found_synonym = True
                        break
            
            # Also check reverse cache for uploaded genes
            if not found_synonym:
                for uploaded_gene in uploaded_genes:
                    if uploaded_gene.lower() in self.reverse_synonym_cache:
                        if self.reverse_synonym_cache[uploaded_gene.lower()] == gene:
                            synonym_mappings[uploaded_gene] = gene
                            found_via_synonyms += 1
                            found_synonym = True
                            break
            
            if not found_synonym:
                still_missing.append(gene)
        
        result = {
            'total_expected': len(expected_genes),
            'found_direct': len(expected_genes) - len(missing_genes),
            'found_via_synonyms': found_via_synonyms,
            'still_missing': len(still_missing),
            'synonym_mappings': synonym_mappings,
            'missing_genes': still_missing
        }
        
        print(f"Synonym lookup complete:")
        print(f"   Direct matches: {result['found_direct']}")
        print(f"   Found via synonyms: {result['found_via_synonyms']}")
        print(f"   Still missing: {result['still_missing']}")
        
        return result
    
    def resolve_gene_name(self, gene_name: str) -> str:
        """
        Resolve a gene name to its canonical form using cache only
        
        Args:
            gene_name: Gene name to resolve
        
        Returns:
            Canonical gene name, or original if not found
        """
        # First check if it's already a known gene
        if gene_name in GENE_NAME_TO_ID:
            return gene_name
        
        # Check reverse synonym cache
        if gene_name.lower() in self.reverse_synonym_cache:
            return self.reverse_synonym_cache[gene_name.lower()]
        
        # If not found, return original
        return gene_name
    
    def get_cache_stats(self) -> Dict:
        """Get statistics about the synonym cache"""
        total_genes = len(self.synonym_cache)
        genes_with_synonyms = sum(1 for syns in self.synonym_cache.values() if syns)
        total_synonyms = len(self.reverse_synonym_cache)
        
        return {
            'total_genes': total_genes,
            'genes_with_synonyms': genes_with_synonyms,
            'total_synonym_entries': total_synonyms,
            'cache_coverage': f"{genes_with_synonyms/total_genes*100:.1f}%" if total_genes > 0 else "0%"
        }

# Global instance
synonym_service = SynonymService()
