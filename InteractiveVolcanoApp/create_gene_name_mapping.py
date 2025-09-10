#!/usr/bin/env python3
"""
Create gene name mapping between expression cluster gene names and volcano plot gene names.
This script will help resolve the gene name mismatch issue.
"""

import json
import os
import pandas as pd
from typing import Dict, List, Set

def load_gene_synonyms() -> Dict[str, List[str]]:
    """Load gene synonyms from the existing cache file."""
    synonyms_file = os.path.join('static', 'api_cache', 'gene_synonyms.json')
    if os.path.exists(synonyms_file):
        with open(synonyms_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
            return data.get('synonyms', {})
    return {}

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

def create_gene_mapping() -> Dict[str, Dict]:
    """Create a comprehensive gene name mapping."""
    print("🔍 Loading gene data...")
    
    # Load all gene name sources
    synonyms = load_gene_synonyms()
    cluster_genes = load_expression_cluster_genes()
    plot_genes = load_volcano_plot_genes()
    
    print(f"📊 Found {len(cluster_genes)} expression cluster genes")
    print(f"📊 Found {len(plot_genes)} volcano plot genes")
    print(f"📊 Found {len(synonyms)} gene synonyms")
    
    # Create mapping structure
    mapping = {
        'cluster_to_plot': {},  # Expression cluster gene -> Volcano plot gene
        'plot_to_cluster': {},  # Volcano plot gene -> Expression cluster gene
        'synonyms': synonyms,
        'unmapped_cluster_genes': [],
        'unmapped_plot_genes': []
    }
    
    # Direct matches
    direct_matches = cluster_genes.intersection(plot_genes)
    print(f"✅ Found {len(direct_matches)} direct matches")
    
    for gene in direct_matches:
        mapping['cluster_to_plot'][gene] = gene
        mapping['plot_to_cluster'][gene] = gene
    
    # Try to find matches through synonyms
    synonym_matches = 0
    for cluster_gene in cluster_genes - direct_matches:
        # Check if cluster gene is in synonyms
        if cluster_gene in synonyms:
            for synonym in synonyms[cluster_gene]:
                if synonym in plot_genes:
                    mapping['cluster_to_plot'][cluster_gene] = synonym
                    mapping['plot_to_cluster'][synonym] = cluster_gene
                    synonym_matches += 1
                    print(f"🔗 Synonym match: {cluster_gene} -> {synonym}")
                    break
        
        # Check if any plot gene has this cluster gene as a synonym
        for plot_gene, plot_synonyms in synonyms.items():
            if plot_gene in plot_genes and cluster_gene in plot_synonyms:
                mapping['cluster_to_plot'][cluster_gene] = plot_gene
                mapping['plot_to_cluster'][plot_gene] = cluster_gene
                synonym_matches += 1
                print(f"🔗 Reverse synonym match: {cluster_gene} <- {plot_gene}")
                break
    
    print(f"🔗 Found {synonym_matches} synonym matches")
    
    # Track unmapped genes
    mapped_cluster_genes = set(mapping['cluster_to_plot'].keys())
    mapped_plot_genes = set(mapping['plot_to_cluster'].keys())
    
    mapping['unmapped_cluster_genes'] = list(cluster_genes - mapped_cluster_genes)
    mapping['unmapped_plot_genes'] = list(plot_genes - mapped_plot_genes)
    
    print(f"❌ {len(mapping['unmapped_cluster_genes'])} cluster genes unmapped")
    print(f"❌ {len(mapping['unmapped_plot_genes'])} plot genes unmapped")
    
    # Show some examples
    print("\n📋 Sample mappings:")
    for i, (cluster_gene, plot_gene) in enumerate(list(mapping['cluster_to_plot'].items())[:10]):
        print(f"  {cluster_gene} -> {plot_gene}")
    
    return mapping

def save_mapping(mapping: Dict[str, Dict], output_file: str):
    """Save the gene mapping to a JSON file."""
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(mapping, f, indent=2, ensure_ascii=False)
    print(f"💾 Saved gene mapping to {output_file}")

def main():
    """Main function to create and save the gene mapping."""
    print("🧬 Creating gene name mapping between expression clusters and volcano plot data...")
    
    # Create the mapping
    mapping = create_gene_mapping()
    
    # Save to the data directory
    output_file = os.path.join('src', 'data', 'gene_name_mapping.json')
    save_mapping(mapping, output_file)
    
    # Print summary
    total_mappings = len(mapping['cluster_to_plot'])
    total_cluster_genes = len(mapping['cluster_to_plot']) + len(mapping['unmapped_cluster_genes'])
    coverage = (total_mappings / total_cluster_genes * 100) if total_cluster_genes > 0 else 0
    
    print(f"\n📈 Mapping Summary:")
    print(f"  Total mappings: {total_mappings}")
    print(f"  Coverage: {coverage:.1f}%")
    print(f"  Unmapped cluster genes: {len(mapping['unmapped_cluster_genes'])}")
    print(f"  Unmapped plot genes: {len(mapping['unmapped_plot_genes'])}")

if __name__ == "__main__":
    main()
