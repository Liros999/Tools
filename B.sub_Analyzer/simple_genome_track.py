#!/usr/bin/env python3
"""
Simple Genome Track Visualization Tool for B.subtilis Expression Data
Creates a linear genome track showing log2fc values mapped to physical chromosome locations.
"""

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional

class SimpleGenomeTrackVisualizer:
    def __init__(self, gene_data_file: str, genome_length: int = 4215606):
        """
        Initialize the genome track visualizer.
        
        Args:
            gene_data_file: Path to the gene coordinates JSON file
            genome_length: Total length of B.subtilis genome (default: 4215606 bp)
        """
        self.genome_length = genome_length
        self.gene_data = self._load_gene_data(gene_data_file)
        
    def _load_gene_data(self, gene_data_file: str) -> Dict:
        """Load gene coordinates from JSON file."""
        try:
            with open(gene_data_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print(f"Warning: Gene data file {gene_data_file} not found.")
            return {}
    
    def _load_synonyms_data(self) -> Dict:
        """Load gene synonyms data from JSON file."""
        try:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            synonyms_file = os.path.join(os.path.dirname(current_dir), "static", "api_cache", "gene_synonyms.json")
            
            if os.path.exists(synonyms_file):
                with open(synonyms_file, 'r') as f:
                    return json.load(f)
            else:
                print(f"Warning: Synonyms file not found at {synonyms_file}")
                return {}
        except Exception as e:
            print(f"Warning: Could not load synonyms data: {e}")
            return {}
    
    def _get_gene_coordinates(self, gene_name: str, synonyms_data: Dict = None) -> Optional[Tuple[int, int, str]]:
        """
        Get coordinates for a gene, including synonym checking.
        
        Args:
            gene_name: Name of the gene to look up
            synonyms_data: Dictionary containing gene synonyms
            
        Returns:
            Tuple of (start, end, strand) or None if not found
        """
        # Try exact match first
        if gene_name in self.gene_data:
            gene_info = self.gene_data[gene_name]
            return (gene_info['start'], gene_info['end'], gene_info['strand'])
        
        # Try case-insensitive match
        for name, info in self.gene_data.items():
            if name.lower() == gene_name.lower():
                return (info['start'], info['end'], info['strand'])
        
        # Try synonym matching if synonyms data is provided
        if synonyms_data and 'reverse_synonyms' in synonyms_data:
            reverse_synonyms = synonyms_data['reverse_synonyms']
            if gene_name.lower() in reverse_synonyms:
                canonical_name = reverse_synonyms[gene_name.lower()]
                if canonical_name in self.gene_data:
                    gene_info = self.gene_data[canonical_name]
                    return (gene_info['start'], gene_info['end'], gene_info['strand'])
        
        return None
    
    def create_genome_track(self, data_file: str, output_prefix: str, log2fc_column: str):
        """
        Create a simple genome track visualization.
        
        Args:
            data_file: Path to CSV file with expression data
            output_prefix: Prefix for output files
            log2fc_column: Name of the log2fc column to use
        """
        print(f"Processing {log2fc_column}...")
        
        # Load data
        df = pd.read_csv(data_file)
        synonyms_data = self._load_synonyms_data()
        
        total_genes_in_file = len(df)
        print(f"  Total genes in file: {total_genes_in_file}")
        
        # Get gene coordinates and log2fc values
        gene_data = []
        genes_with_coords = 0
        genes_with_valid_log2fc = 0
        genes_with_synonyms = 0
        
        for _, row in df.iterrows():
            gene_name = str(row['geneName'])
            if gene_name != 'nan' and gene_name != 'NA':
                coords = self._get_gene_coordinates(gene_name, synonyms_data)
                if coords:
                    genes_with_coords += 1
                    if not pd.isna(row[log2fc_column]):
                        genes_with_valid_log2fc += 1
                        gene_data.append({
                            'gene': gene_name,
                            'start': coords[0],
                            'end': coords[1],
                            'strand': coords[2],
                            'log2fc': row[log2fc_column]
                        })
                        
                        # Check if this was found via synonym
                        if synonyms_data and 'reverse_synonyms' in synonyms_data:
                            if gene_name.lower() in synonyms_data['reverse_synonyms']:
                                genes_with_synonyms += 1
        
        # Print statistics
        print(f"  Genes with coordinates found: {genes_with_coords}")
        if genes_with_synonyms > 0:
            print(f"  Genes found via synonyms: {genes_with_synonyms}")
        print(f"  Genes with valid log2fc values: {genes_with_valid_log2fc}")
        print(f"  Genes successfully analyzed: {len(gene_data)}")
        
        if not gene_data:
            print(f"No valid gene data found for {log2fc_column}")
            return 0
        
        # Create genome track visualization
        self._create_simple_genome_track(gene_data, output_prefix, log2fc_column)
        
        return len(gene_data)
    
    def _create_simple_genome_track(self, gene_data: List[Dict], output_prefix: str, title: str):
        """
        Create a simple linear genome track visualization.
        
        Args:
            gene_data: List of gene data dictionaries
            output_prefix: Prefix for output file
            title: Title for the plot
        """
        # Sort by genome position
        gene_data.sort(key=lambda x: x['start'])
        
        # Create figure
        fig, ax = plt.subplots(figsize=(20, 8))
        
        # Convert genome positions to Mbp for better readability
        genome_length_mbp = self.genome_length / 1000000
        
        # Draw genome backbone
        ax.plot([0, genome_length_mbp], [0, 0], 'k-', linewidth=2, alpha=0.5)
        
        # Plot genes as bars
        for gene in gene_data:
            start_mbp = gene['start'] / 1000000
            end_mbp = gene['end'] / 1000000
            log2fc = gene['log2fc']
            
            # Color based on log2fc value
            if log2fc > 0:
                color = 'red'
                alpha = min(0.8, abs(log2fc) / 5)  # Scale alpha by magnitude
            else:
                color = 'blue'
                alpha = min(0.8, abs(log2fc) / 5)
            
            # Draw gene bar
            ax.barh(y=0, width=end_mbp - start_mbp, left=start_mbp, 
                   height=log2fc, color=color, alpha=alpha, edgecolor='black', linewidth=0.5)
        
        # Add zero line
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=1)
        
        # Customize plot
        ax.set_xlim(0, genome_length_mbp)
        ax.set_xlabel('Genome Position (Mbp)', fontsize=12)
        ax.set_ylabel('Log2 Fold Change', fontsize=12)
        ax.set_title(f'B. subtilis Genome Track - {title}', fontsize=16, pad=20)
        
        # Add genome position markers every 0.5 Mbp
        tick_positions = np.arange(0, genome_length_mbp + 0.5, 0.5)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([f'{pos:.1f}' for pos in tick_positions], rotation=0)
        
        # Add grid
        ax.grid(True, alpha=0.3, axis='x')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', alpha=0.7, label='Up-regulated'),
            Patch(facecolor='blue', alpha=0.7, label='Down-regulated')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_genome_track.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Genome track saved: {output_prefix}_genome_track.png")

def main():
    """Main function to create genome track visualizations for all datasets."""
    
    # Get paths
    current_dir = os.path.dirname(os.path.abspath(__file__))
    gene_data_file = os.path.join(current_dir, "data", "all_gene_data.json")
    volcano_app_dir = os.path.join(os.path.dirname(current_dir), "InteractiveVolcanoApp")
    
    # Check files
    if not os.path.exists(gene_data_file):
        print(f"Error: Gene data file not found at {gene_data_file}")
        return
    
    # Initialize visualizer
    visualizer = SimpleGenomeTrackVisualizer(gene_data_file)
    
    # Create output directory
    output_dir = os.path.join(current_dir, "results", "genome_tracks")
    os.makedirs(output_dir, exist_ok=True)
    
    print("Creating simple genome track visualizations...")
    print(f"Output directory: {output_dir}")
    
    # Track statistics
    total_genes_analyzed = 0
    
    # Process D1-3 data
    d1_3_file = os.path.join(volcano_app_dir, 'Ca_D1_3_6_data_with_Regulators.csv')
    if os.path.exists(d1_3_file):
        print("\nProcessing D1-3 data...")
        output_prefix = os.path.join(output_dir, 'D1-3')
        gene_count = visualizer.create_genome_track(d1_3_file, output_prefix, 'log2FoldChange_D1-3')
        if gene_count:
            total_genes_analyzed += gene_count
    
    # Process D6 data
    if os.path.exists(d1_3_file):
        print("\nProcessing D6 data...")
        output_prefix = os.path.join(output_dir, 'D6')
        gene_count = visualizer.create_genome_track(d1_3_file, output_prefix, 'log2FoldChangeD-6')
        if gene_count:
            total_genes_analyzed += gene_count
    
    # Process CueR data
    cueR_file = os.path.join(volcano_app_dir, 'mmc4_CueR_with_Regulators.csv')
    if os.path.exists(cueR_file):
        print("\nProcessing CueR data...")
        output_prefix = os.path.join(output_dir, 'CueR')
        gene_count = visualizer.create_genome_track(cueR_file, output_prefix, 'log2FoldChange')
        if gene_count:
            total_genes_analyzed += gene_count
    
    print(f"\n=== SUMMARY ===")
    print(f"Total genes analyzed across all datasets: {total_genes_analyzed}")
    print(f"All genome track visualizations completed! Check: {output_dir}")

if __name__ == "__main__":
    main()



