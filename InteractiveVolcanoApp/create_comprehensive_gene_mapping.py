#!/usr/bin/env python3
"""
Create comprehensive gene name mapping using SubtiWiki API.
This script fetches all genes from SubtiWiki and creates a complete mapping
between different gene naming systems.
"""

import json
import os
import requests
import time
from typing import Dict, List, Set, Optional
import pandas as pd

# SubtiWiki API base URL
SUBTIWIKI_API_BASE = "https://www.subtiwiki.uni-goettingen.de/v5/api"

def get_all_genes_from_subtiwiki() -> List[Dict]:
    """
    Get all genes from SubtiWiki API.
    
    Returns:
        List of gene information from SubtiWiki
    """
    print("🔍 Fetching all genes from SubtiWiki API...")
    
    try:
        # Get all genes
        genes_url = f"{SUBTIWIKI_API_BASE}/gene/"
        response = requests.get(genes_url, timeout=60)
        
        if response.status_code == 200:
            genes = response.json()
            print(f"✅ Fetched {len(genes)} genes from SubtiWiki")
            return genes
        else:
            print(f"❌ Error fetching genes: HTTP {response.status_code}")
            return []
            
    except Exception as e:
        print(f"❌ Error fetching genes from SubtiWiki: {e}")
        return []

def get_gene_details_batch(gene_ids: List[int], batch_size: int = 10) -> Dict[int, Dict]:
    """
    Get detailed information for a batch of genes.
    
    Args:
        gene_ids: List of gene IDs
        batch_size: Number of genes to process in each batch
        
    Returns:
        Dict mapping gene ID to gene details
    """
    gene_details = {}
    
    for i in range(0, len(gene_ids), batch_size):
        batch = gene_ids[i:i + batch_size]
        print(f"⏳ Processing batch {i//batch_size + 1}/{(len(gene_ids) + batch_size - 1)//batch_size} ({len(batch)} genes)")
        
        for gene_id in batch:
            try:
                gene_url = f"{SUBTIWIKI_API_BASE}/gene/{gene_id}"
                response = requests.get(gene_url, timeout=30)
                
                if response.status_code == 200:
                    gene_details[gene_id] = response.json()
                
                # Rate limiting
                time.sleep(0.1)
                
            except Exception as e:
                print(f"Error getting details for gene {gene_id}: {e}")
                continue
    
    return gene_details

def extract_gene_names_and_synonyms(gene_details: Dict[int, Dict]) -> Dict[str, List[str]]:
    """
    Extract gene names and synonyms from gene details.
    
    Args:
        gene_details: Dict mapping gene ID to gene details
        
    Returns:
        Dict mapping gene name to list of synonyms
    """
    gene_synonyms = {}
    
    for gene_id, details in gene_details.items():
        gene_name = details.get('name', '')
        if not gene_name:
            continue
        
        synonyms = []
        
        # Add synonyms from the gene details
        if 'synonyms' in details:
            synonyms.extend(details['synonyms'])
        
        # Add outlinks which might contain alternative names
        if 'outlinks' in details:
            outlinks = details['outlinks']
            if 'subtilist' in outlinks:
                synonyms.append(outlinks['subtilist'])
            if 'expression_browser' in outlinks:
                # Extract gene name from expression browser ID
                expr_id = outlinks['expression_browser']
                if '_' in expr_id:
                    gene_part = expr_id.split('_')[0]
                    synonyms.append(gene_part)
        
        # Remove duplicates and empty strings
        synonyms = list(set([s for s in synonyms if s and s.strip()]))
        
        if synonyms:
            gene_synonyms[gene_name] = synonyms
    
    return gene_synonyms

def load_expression_cluster_genes() -> Set[str]:
    """Load gene names from the expression cluster mapping."""
    clusters_file = os.path.join('src', 'data', 'gene_expression_clusters.json')
    if os.path.exists(clusters_file):
        with open(clusters_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
            return set(data.keys())
    return set()

def load_volcano_plot_genes() -> Set[str]:
    """Load gene names from the volcano plot data files."""
    plot_genes = set()
    
    # Check both CSV files for gene names
    csv_files = [
        'Ca_D1_3_6_data_with_Regulators.csv',
        'mmc4_CueR_with_Regulators.csv'
    ]
    
    for csv_file in csv_files:
        if os.path.exists(csv_file):
            try:
                df = pd.read_csv(csv_file)
                if 'geneName' in df.columns:
                    plot_genes.update(df['geneName'].dropna().tolist())
                print(f"Loaded {len(df['geneName'].dropna())} genes from {csv_file}")
            except Exception as e:
                print(f"Error reading {csv_file}: {e}")
    
    return plot_genes

def create_comprehensive_mapping(cluster_genes: Set[str], plot_genes: Set[str], gene_synonyms: Dict[str, List[str]]) -> Dict:
    """
    Create a comprehensive gene mapping using SubtiWiki synonyms.
    
    Args:
        cluster_genes: Set of genes from expression clusters
        plot_genes: Set of genes from volcano plot data
        gene_synonyms: Dict mapping gene names to synonyms from SubtiWiki
        
    Returns:
        Comprehensive gene mapping
    """
    print("🔗 Creating comprehensive gene mapping...")
    
    mapping = {
        'cluster_to_plot': {},
        'plot_to_cluster': {},
        'subtiwiki_synonyms': gene_synonyms,
        'unmapped_cluster_genes': [],
        'unmapped_plot_genes': []
    }
    
    # Direct matches
    direct_matches = cluster_genes.intersection(plot_genes)
    print(f"✅ Found {len(direct_matches)} direct matches")
    
    for gene in direct_matches:
        mapping['cluster_to_plot'][gene] = gene
        mapping['plot_to_cluster'][gene] = gene
    
    # Try to find matches through SubtiWiki synonyms
    synonym_matches = 0
    
    for cluster_gene in cluster_genes - direct_matches:
        # Check if cluster gene is in SubtiWiki synonyms
        if cluster_gene in gene_synonyms:
            for synonym in gene_synonyms[cluster_gene]:
                if synonym in plot_genes:
                    mapping['cluster_to_plot'][cluster_gene] = synonym
                    mapping['plot_to_cluster'][synonym] = cluster_gene
                    synonym_matches += 1
                    print(f"🔗 SubtiWiki synonym match: {cluster_gene} -> {synonym}")
                    break
        
        # Check if any plot gene has this cluster gene as a synonym
        for plot_gene in plot_genes - set(mapping['plot_to_cluster'].keys()):
            if plot_gene in gene_synonyms and cluster_gene in gene_synonyms[plot_gene]:
                mapping['cluster_to_plot'][cluster_gene] = plot_gene
                mapping['plot_to_cluster'][plot_gene] = cluster_gene
                synonym_matches += 1
                print(f"🔗 Reverse SubtiWiki synonym match: {cluster_gene} <- {plot_gene}")
                break
    
    print(f"🔗 Found {synonym_matches} synonym matches through SubtiWiki")
    
    # Track unmapped genes
    mapped_cluster_genes = set(mapping['cluster_to_plot'].keys())
    mapped_plot_genes = set(mapping['plot_to_cluster'].keys())
    
    mapping['unmapped_cluster_genes'] = list(cluster_genes - mapped_cluster_genes)
    mapping['unmapped_plot_genes'] = list(plot_genes - mapped_plot_genes)
    
    print(f"❌ {len(mapping['unmapped_cluster_genes'])} cluster genes unmapped")
    print(f"❌ {len(mapping['unmapped_plot_genes'])} plot genes unmapped")
    
    return mapping

def save_mapping(mapping: Dict, output_file: str):
    """Save the gene mapping to a JSON file."""
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(mapping, f, indent=2, ensure_ascii=False)
    print(f"💾 Saved comprehensive gene mapping to {output_file}")

def main():
    """Main function to create comprehensive gene mapping."""
    print("🧬 Creating comprehensive gene name mapping using SubtiWiki API...")
    
    # Get all genes from SubtiWiki
    genes = get_all_genes_from_subtiwiki()
    if not genes:
        print("❌ Failed to fetch genes from SubtiWiki. Exiting.")
        return
    
    # Extract gene IDs
    gene_ids = [gene['id'] for gene in genes if 'id' in gene]
    print(f"📊 Processing {len(gene_ids)} genes from SubtiWiki")
    
    # Get detailed information for all genes
    gene_details = get_gene_details_batch(gene_ids)
    print(f"📊 Retrieved details for {len(gene_details)} genes")
    
    # Extract gene names and synonyms
    gene_synonyms = extract_gene_names_and_synonyms(gene_details)
    print(f"📊 Extracted synonyms for {len(gene_synonyms)} genes")
    
    # Load expression cluster genes
    cluster_genes = load_expression_cluster_genes()
    print(f"📊 Loaded {len(cluster_genes)} expression cluster genes")
    
    # Load volcano plot genes
    plot_genes = load_volcano_plot_genes()
    print(f"📊 Loaded {len(plot_genes)} volcano plot genes")
    
    # Create comprehensive mapping
    mapping = create_comprehensive_mapping(cluster_genes, plot_genes, gene_synonyms)
    
    # Save the mapping
    output_file = os.path.join('src', 'data', 'gene_name_mapping_comprehensive.json')
    save_mapping(mapping, output_file)
    
    # Print summary
    total_mappings = len(mapping['cluster_to_plot'])
    total_cluster_genes = len(mapping['cluster_to_plot']) + len(mapping['unmapped_cluster_genes'])
    coverage = (total_mappings / total_cluster_genes * 100) if total_cluster_genes > 0 else 0
    
    print(f"\n📈 Comprehensive Mapping Summary:")
    print(f"  Total mappings: {total_mappings}")
    print(f"  Coverage: {coverage:.1f}%")
    print(f"  Unmapped cluster genes: {len(mapping['unmapped_cluster_genes'])}")
    print(f"  Unmapped plot genes: {len(mapping['unmapped_plot_genes'])}")
    print(f"  SubtiWiki genes with synonyms: {len(mapping['subtiwiki_synonyms'])}")
    
    # Show some examples
    print(f"\n📋 Sample mappings:")
    for i, (cluster_gene, plot_gene) in enumerate(list(mapping['cluster_to_plot'].items())[:10]):
        print(f"  {cluster_gene} -> {plot_gene}")
    
    # Show some examples of SubtiWiki synonyms
    print(f"\n📋 Sample SubtiWiki synonyms:")
    for i, (gene_name, synonyms) in enumerate(list(mapping['subtiwiki_synonyms'].items())[:5]):
        print(f"  {gene_name}: {synonyms}")

if __name__ == "__main__":
    main()
