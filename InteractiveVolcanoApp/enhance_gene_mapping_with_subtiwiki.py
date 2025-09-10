#!/usr/bin/env python3
"""
Enhance gene name mapping using SubtiWiki API to resolve unmapped genes.
This script will use the SubtiWiki API to find gene name mappings for genes
that couldn't be matched directly between expression clusters and volcano plot data.
"""

import json
import os
import requests
import time
from typing import Dict, List, Set, Optional
import pandas as pd

# SubtiWiki API base URL
SUBTIWIKI_API_BASE = "https://www.subtiwiki.uni-goettingen.de/v5/api"

def load_existing_mapping() -> Dict:
    """Load the existing gene name mapping."""
    mapping_file = os.path.join('src', 'data', 'gene_name_mapping.json')
    if os.path.exists(mapping_file):
        with open(mapping_file, 'r', encoding='utf-8') as f:
            return json.load(f)
    return {}

def search_gene_in_subtiwiki(gene_name: str) -> Optional[Dict]:
    """
    Search for a gene in SubtiWiki API.
    
    Args:
        gene_name: The gene name to search for
        
    Returns:
        Dict with gene information if found, None otherwise
    """
    try:
        # Search for the gene
        search_url = f"{SUBTIWIKI_API_BASE}/search/"
        params = {
            'q': gene_name,
            'category': 'Gene',
            'mode': 'exact'
        }
        
        response = requests.get(search_url, params=params, timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            
            # Check for exact hits
            exact_hits = data.get('exact_hits', [])
            for hit in exact_hits:
                if hit.get('category') == 'Gene' and hit.get('name', '').lower() == gene_name.lower():
                    # Get detailed gene information
                    gene_id = hit.get('id')
                    if gene_id:
                        return get_gene_details_from_subtiwiki(gene_id)
            
            # Check for partial hits
            partial_hits = data.get('partial_hits', [])
            for hit in partial_hits:
                if hit.get('category') == 'Gene' and hit.get('name', '').lower() == gene_name.lower():
                    gene_id = hit.get('id')
                    if gene_id:
                        return get_gene_details_from_subtiwiki(gene_id)
        
        return None
        
    except Exception as e:
        print(f"Error searching for gene {gene_name}: {e}")
        return None

def get_gene_details_from_subtiwiki(gene_id: int) -> Optional[Dict]:
    """
    Get detailed gene information from SubtiWiki API.
    
    Args:
        gene_id: The gene ID from SubtiWiki
        
    Returns:
        Dict with gene details if found, None otherwise
    """
    try:
        gene_url = f"{SUBTIWIKI_API_BASE}/gene/{gene_id}"
        response = requests.get(gene_url, timeout=30)
        
        if response.status_code == 200:
            return response.json()
        
        return None
        
    except Exception as e:
        print(f"Error getting gene details for ID {gene_id}: {e}")
        return None

def get_gene_synonyms_from_subtiwiki(gene_details: Dict) -> List[str]:
    """
    Extract synonyms from SubtiWiki gene details.
    
    Args:
        gene_details: Gene details from SubtiWiki API
        
    Returns:
        List of gene synonyms
    """
    synonyms = []
    
    # Get synonyms from the gene details
    if 'synonyms' in gene_details:
        synonyms.extend(gene_details['synonyms'])
    
    # Get outlinks which might contain alternative names
    if 'outlinks' in gene_details:
        outlinks = gene_details['outlinks']
        if 'subtilist' in outlinks:
            synonyms.append(outlinks['subtilist'])
        if 'expression_browser' in outlinks:
            # Extract gene name from expression browser ID
            expr_id = outlinks['expression_browser']
            if '_' in expr_id:
                gene_part = expr_id.split('_')[0]
                synonyms.append(gene_part)
    
    return list(set(synonyms))  # Remove duplicates

def check_gene_in_plot_data(gene_name: str, plot_genes: Set[str]) -> bool:
    """
    Check if a gene name exists in the volcano plot data.
    
    Args:
        gene_name: The gene name to check
        plot_genes: Set of genes in the volcano plot data
        
    Returns:
        True if gene exists in plot data, False otherwise
    """
    return gene_name.lower() in {g.lower() for g in plot_genes}

def enhance_mapping_with_subtiwiki(mapping: Dict, plot_genes: Set[str]) -> Dict:
    """
    Enhance the gene mapping using SubtiWiki API.
    
    Args:
        mapping: Existing gene mapping
        plot_genes: Set of genes in volcano plot data
        
    Returns:
        Enhanced mapping with additional matches from SubtiWiki
    """
    print("🔍 Enhancing gene mapping with SubtiWiki API...")
    
    # Get unmapped cluster genes
    unmapped_cluster_genes = mapping.get('unmapped_cluster_genes', [])
    print(f"📊 Found {len(unmapped_cluster_genes)} unmapped cluster genes to process")
    
    # Track new mappings
    new_mappings = 0
    processed_count = 0
    
    for cluster_gene in unmapped_cluster_genes:
        processed_count += 1
        if processed_count % 50 == 0:
            print(f"⏳ Processed {processed_count}/{len(unmapped_cluster_genes)} genes...")
        
        # Search for the gene in SubtiWiki
        gene_details = search_gene_in_subtiwiki(cluster_gene)
        
        if gene_details:
            # Get synonyms from SubtiWiki
            synonyms = get_gene_synonyms_from_subtiwiki(gene_details)
            
            # Check if any synonym matches plot genes
            for synonym in synonyms:
                if check_gene_in_plot_data(synonym, plot_genes):
                    # Find the exact match in plot_genes (case-insensitive)
                    for plot_gene in plot_genes:
                        if plot_gene.lower() == synonym.lower():
                            # Add the mapping
                            mapping['cluster_to_plot'][cluster_gene] = plot_gene
                            mapping['plot_to_cluster'][plot_gene] = cluster_gene
                            new_mappings += 1
                            print(f"🔗 SubtiWiki match: {cluster_gene} -> {plot_gene} (via synonym: {synonym})")
                            break
                    break
        
        # Rate limiting to be respectful to the API
        time.sleep(0.1)
    
    print(f"✅ Enhanced mapping with {new_mappings} new matches from SubtiWiki")
    return mapping

def update_unmapped_lists(mapping: Dict, plot_genes: Set[str]):
    """
    Update the unmapped gene lists after enhancement.
    
    Args:
        mapping: The gene mapping
        plot_genes: Set of genes in volcano plot data
    """
    # Update unmapped cluster genes
    mapped_cluster_genes = set(mapping['cluster_to_plot'].keys())
    all_cluster_genes = set(mapping['cluster_to_plot'].keys()) | set(mapping.get('unmapped_cluster_genes', []))
    mapping['unmapped_cluster_genes'] = list(all_cluster_genes - mapped_cluster_genes)
    
    # Update unmapped plot genes
    mapped_plot_genes = set(mapping['plot_to_cluster'].keys())
    mapping['unmapped_plot_genes'] = list(plot_genes - mapped_plot_genes)

def save_enhanced_mapping(mapping: Dict, output_file: str):
    """
    Save the enhanced gene mapping.
    
    Args:
        mapping: The enhanced gene mapping
        output_file: Output file path
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(mapping, f, indent=2, ensure_ascii=False)
    print(f"💾 Saved enhanced gene mapping to {output_file}")

def load_plot_genes() -> Set[str]:
    """
    Load gene names from volcano plot data files.
    
    Returns:
        Set of gene names in volcano plot data
    """
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

def main():
    """Main function to enhance gene mapping with SubtiWiki API."""
    print("🧬 Enhancing gene name mapping with SubtiWiki API...")
    
    # Load existing mapping
    mapping = load_existing_mapping()
    if not mapping:
        print("❌ No existing gene mapping found. Please run create_gene_name_mapping.py first.")
        return
    
    # Load plot genes
    plot_genes = load_plot_genes()
    print(f"📊 Loaded {len(plot_genes)} genes from volcano plot data")
    
    # Show initial statistics
    initial_mappings = len(mapping['cluster_to_plot'])
    initial_unmapped = len(mapping.get('unmapped_cluster_genes', []))
    print(f"📈 Initial state: {initial_mappings} mappings, {initial_unmapped} unmapped cluster genes")
    
    # Enhance mapping with SubtiWiki API
    enhanced_mapping = enhance_mapping_with_subtiwiki(mapping, plot_genes)
    
    # Update unmapped lists
    update_unmapped_lists(enhanced_mapping, plot_genes)
    
    # Show final statistics
    final_mappings = len(enhanced_mapping['cluster_to_plot'])
    final_unmapped = len(enhanced_mapping.get('unmapped_cluster_genes', []))
    new_mappings = final_mappings - initial_mappings
    
    print(f"\n📈 Final Statistics:")
    print(f"  Total mappings: {final_mappings}")
    print(f"  New mappings from SubtiWiki: {new_mappings}")
    print(f"  Remaining unmapped cluster genes: {final_unmapped}")
    print(f"  Unmapped plot genes: {len(enhanced_mapping.get('unmapped_plot_genes', []))}")
    
    # Calculate coverage
    total_cluster_genes = final_mappings + final_unmapped
    coverage = (final_mappings / total_cluster_genes * 100) if total_cluster_genes > 0 else 0
    print(f"  Coverage: {coverage:.1f}%")
    
    # Save enhanced mapping
    output_file = os.path.join('src', 'data', 'gene_name_mapping_enhanced.json')
    save_enhanced_mapping(enhanced_mapping, output_file)
    
    # Show some examples of new mappings
    if new_mappings > 0:
        print(f"\n📋 Examples of new mappings from SubtiWiki:")
        new_mapping_items = list(enhanced_mapping['cluster_to_plot'].items())[initial_mappings:initial_mappings+10]
        for cluster_gene, plot_gene in new_mapping_items:
            print(f"  {cluster_gene} -> {plot_gene}")

if __name__ == "__main__":
    main()
