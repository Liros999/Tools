#!/usr/bin/env python3
"""
Simple Genomic Visualization Tool for B.subtilis Expression Data
Creates basic heat maps and bar plots showing log2fc values mapped to genome coordinates.
"""

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional

class SimpleGenomicVisualizer:
    def __init__(self, gene_data_file: str):
        """Initialize the visualizer."""
        self.gene_data = self._load_gene_data(gene_data_file)
        
    def _load_gene_data(self, gene_data_file: str) -> Dict:
        """Load gene coordinates from JSON file."""
        try:
            with open(gene_data_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print(f"Warning: Gene data file {gene_data_file} not found.")
            return {}
    
    def _get_gene_coordinates(self, gene_name: str) -> Optional[Tuple[int, int, str]]:
        """Get coordinates for a gene."""
        if gene_name in self.gene_data:
            gene_info = self.gene_data[gene_name]
            return (gene_info['start'], gene_info['end'], gene_info['strand'])
        
        # Try case-insensitive match
        for name, info in self.gene_data.items():
            if name.lower() == gene_name.lower():
                return (info['start'], info['end'], info['strand'])
        
        return None
    
    def create_visualizations(self, data_file: str, output_prefix: str, log2fc_column: str):
        """Create both heat map and bar plot for a dataset."""
        print(f"Processing {log2fc_column}...")
        
        # Load data
        df = pd.read_csv(data_file)
        total_genes_in_file = len(df)
        print(f"  Total genes in file: {total_genes_in_file}")
        
        # Get gene coordinates and log2fc values
        gene_data = []
        genes_with_coords = 0
        genes_with_valid_log2fc = 0
        
        for _, row in df.iterrows():
            gene_name = str(row['geneName'])
            if gene_name != 'nan' and gene_name != 'NA':
                coords = self._get_gene_coordinates(gene_name)
                if coords:
                    genes_with_coords += 1
                    if not pd.isna(row[log2fc_column]):
                        genes_with_valid_log2fc += 1
                        gene_data.append({
                            'gene': gene_name,
                            'start': coords[0],
                            'log2fc': row[log2fc_column]
                        })
        
        # Print statistics
        print(f"  Genes with coordinates found: {genes_with_coords}")
        print(f"  Genes with valid log2fc values: {genes_with_valid_log2fc}")
        print(f"  Genes successfully analyzed: {len(gene_data)}")
        
        if not gene_data:
            print(f"No valid gene data found for {log2fc_column}")
            return 0
        
        # Sort by genome position
        gene_data.sort(key=lambda x: x['start'])
        
        # Create bar plot
        self._create_simple_barplot(gene_data, output_prefix, log2fc_column)
        
        # Create heat map
        self._create_simple_heatmap(gene_data, output_prefix, log2fc_column)
        
        return len(gene_data)
    
    def _create_simple_barplot(self, gene_data: List[Dict], output_prefix: str, title: str):
        """Create a simple bar plot."""
        starts = [g['start'] for g in gene_data]
        log2fc_values = [g['log2fc'] for g in gene_data]
        
        # Convert to Mbp for better readability
        starts_mbp = [s / 1000000 for s in starts]
        
        plt.figure(figsize=(15, 6))
        
        # Color based on log2fc value
        colors = ['red' if x > 0 else 'blue' for x in log2fc_values]
        
        plt.bar(starts_mbp, log2fc_values, width=0.0005, color=colors, alpha=0.7)
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        plt.xlabel('Genome Position (Mbp)')
        plt.ylabel('Log2 Fold Change')
        plt.title(f'Gene Expression - {title}')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_barplot.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Bar plot saved: {output_prefix}_barplot.png")
    
    def _create_simple_heatmap(self, gene_data: List[Dict], output_prefix: str, title: str):
        """Create a simple heat map."""
        # Create genome-wide array
        genome_length = 4215606
        genome_array = np.zeros(genome_length)
        
        for gene in gene_data:
            start = gene['start']
            if 0 <= start < genome_length:
                genome_array[start] = gene['log2fc']
        
        # Find non-zero positions
        non_zero_pos = np.where(genome_array != 0)[0]
        if len(non_zero_pos) == 0:
            print("No gene data for heat map")
            return
        
        # Get values and create heat map
        values = genome_array[non_zero_pos]
        
        # Convert to Mbp for better readability
        non_zero_pos_mbp = non_zero_pos / 1000000
        
        plt.figure(figsize=(15, 4))
        
        # Simple scatter plot as heat map
        plt.scatter(non_zero_pos_mbp, [0] * len(non_zero_pos), 
                   c=values, cmap='RdBu_r', s=10, alpha=0.8)
        
        plt.colorbar(label='Log2 Fold Change')
        plt.xlabel('Genome Position (Mbp)')
        plt.title(f'Gene Expression Heat Map - {title}')
        plt.yticks([])
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_heatmap.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Heat map saved: {output_prefix}_heatmap.png")

def main():
    """Main function to create visualizations for all datasets."""
    
    # Get paths
    current_dir = os.path.dirname(os.path.abspath(__file__))
    gene_data_file = os.path.join(current_dir, "data", "all_gene_data.json")
    volcano_app_dir = os.path.join(os.path.dirname(current_dir), "InteractiveVolcanoApp")
    
    # Check files
    if not os.path.exists(gene_data_file):
        print(f"Error: Gene data file not found at {gene_data_file}")
        return
    
    # Initialize visualizer
    visualizer = SimpleGenomicVisualizer(gene_data_file)
    
    # Create output directory
    output_dir = os.path.join(current_dir, "results", "simple_visualizations")
    os.makedirs(output_dir, exist_ok=True)
    
    print("Creating simple genomic visualizations...")
    print(f"Output directory: {output_dir}")
    
    # Track statistics
    total_genes_analyzed = 0
    
    # Process D1-3 data
    d1_3_file = os.path.join(volcano_app_dir, 'Ca_D1_3_6_data_with_Regulators.csv')
    if os.path.exists(d1_3_file):
        print("\nProcessing D1-3 data...")
        output_prefix = os.path.join(output_dir, 'D1-3')
        gene_count = visualizer.create_visualizations(d1_3_file, output_prefix, 'log2FoldChange_D1-3')
        if gene_count:
            total_genes_analyzed += gene_count
    
    # Process D6 data
    if os.path.exists(d1_3_file):
        print("\nProcessing D6 data...")
        output_prefix = os.path.join(output_dir, 'D6')
        gene_count = visualizer.create_visualizations(d1_3_file, output_prefix, 'log2FoldChangeD-6')
        if gene_count:
            total_genes_analyzed += gene_count
    
    # Process CueR data
    cueR_file = os.path.join(volcano_app_dir, 'mmc4_CueR_with_Regulators.csv')
    if os.path.exists(cueR_file):
        print("\nProcessing CueR data...")
        output_prefix = os.path.join(output_dir, 'CueR')
        gene_count = visualizer.create_visualizations(cueR_file, output_prefix, 'log2FoldChange')
        if gene_count:
            total_genes_analyzed += gene_count
    
    print(f"\n=== SUMMARY ===")
    print(f"Total genes analyzed across all datasets: {total_genes_analyzed}")
    print(f"All visualizations completed! Check: {output_dir}")

if __name__ == "__main__":
    main()
