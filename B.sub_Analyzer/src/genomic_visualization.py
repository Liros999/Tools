#!/usr/bin/env python3
"""
Genomic Visualization Tool for B.subtilis Expression Data
Creates heat maps and bar plots showing log2fc values mapped to genome coordinates.
"""

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

class GenomicVisualizer:
    def __init__(self, gene_data_file: str, genome_length: int = 4215606):
        """
        Initialize the genomic visualizer.
        
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
    
    def _get_gene_coordinates(self, gene_name: str) -> Optional[Tuple[int, int, str]]:
        """
        Get coordinates for a gene.
        
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
        
        return None
    
    def _create_genome_bins(self, bin_size: int = 10000) -> List[Tuple[int, int]]:
        """Create genome bins for heat map visualization."""
        bins = []
        for start in range(1, self.genome_length, bin_size):
            end = min(start + bin_size - 1, self.genome_length)
            bins.append((start, end))
        return bins
    
    def _assign_genes_to_bins(self, gene_coords: Dict[str, Tuple[int, int, str]], 
                            bin_size: int = 10000) -> Dict[int, List[str]]:
        """Assign genes to genome bins."""
        bins = self._create_genome_bins(bin_size)
        bin_assignments = {i: [] for i in range(len(bins))}
        
        for gene_name, (start, end, strand) in gene_coords.items():
            # Find which bin this gene belongs to
            for i, (bin_start, bin_end) in enumerate(bins):
                if start <= bin_end and end >= bin_start:
                    bin_assignments[i].append(gene_name)
        
        return bin_assignments
    
    def create_heatmap_visualization(self, data_file: str, output_prefix: str, 
                                   bin_size: int = 10000, plot_type: str = 'heatmap'):
        """
        Create genomic heat map or bar plot visualization.
        
        Args:
            data_file: Path to CSV file with expression data
            output_prefix: Prefix for output files
            bin_size: Size of genome bins for heat map
            plot_type: 'heatmap' or 'barplot'
        """
        # Load expression data
        df = pd.read_csv(data_file)
        
        # Determine which log2fc column to use based on file content
        log2fc_col = None
        if 'log2FoldChange_D1-3' in df.columns:
            log2fc_col = 'log2FoldChange_D1-3'
        elif 'log2FoldChangeD-6' in df.columns:
            log2fc_col = 'log2FoldChangeD-6'
        elif 'log2FoldChange' in df.columns:
            log2fc_col = 'log2FoldChange'
        else:
            print(f"Error: No log2fc column found in {data_file}")
            return
        
        # Get gene name column
        gene_col = 'geneName' if 'geneName' in df.columns else 'geneName'
        
        # Create gene coordinates mapping
        gene_coords = {}
        for _, row in df.iterrows():
            gene_name = str(row[gene_col])
            if gene_name != 'nan' and gene_name != 'NA':
                coords = self._get_gene_coordinates(gene_name)
                if coords:
                    gene_coords[gene_name] = coords
        
        # Create genome-wide array
        genome_array = np.zeros(self.genome_length)
        gene_positions = []
        
        # Assign log2fc values to genome positions
        for _, row in df.iterrows():
            gene_name = str(row[gene_col])
            if gene_name in gene_coords:
                start, end, strand = gene_coords[gene_name]
                log2fc = row[log2fc_col]
                
                if not pd.isna(log2fc):
                    # Assign the log2fc value to all positions in the gene
                    for pos in range(start, end + 1):
                        if 0 <= pos < len(genome_array):
                            genome_array[pos] = log2fc
                    
                    # Store gene position for bar plot
                    gene_positions.append({
                        'gene': gene_name,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'log2fc': log2fc
                    })
        
        # Create visualization
        if plot_type == 'heatmap':
            self._create_heatmap(genome_array, output_prefix, bin_size)
        elif plot_type == 'barplot':
            self._create_barplot(gene_positions, output_prefix)
        else:
            # Create both
            self._create_heatmap(genome_array, output_prefix, bin_size)
            self._create_barplot(gene_positions, output_prefix)
    
    def _create_heatmap(self, genome_array: np.ndarray, output_prefix: str, bin_size: int):
        """Create heat map visualization."""
        # Create bins and calculate average log2fc per bin
        bins = self._create_genome_bins(bin_size)
        bin_values = []
        
        for start, end in bins:
            # Get values for this bin
            bin_data = genome_array[start-1:end]
            # Calculate mean (excluding zeros)
            non_zero = bin_data[bin_data != 0]
            if len(non_zero) > 0:
                bin_values.append(np.mean(non_zero))
            else:
                bin_values.append(0)
        
        # Create heat map
        plt.figure(figsize=(20, 8))
        
        # Reshape for heat map (assuming circular genome)
        n_bins = len(bin_values)
        # Create a 2D array for heat map
        heatmap_data = np.array(bin_values).reshape(1, -1)
        
        # Create custom colormap
        colors = ['darkblue', 'blue', 'lightblue', 'white', 'lightcoral', 'red', 'darkred']
        n_bins_cmap = 256
        cmap = plt.cm.LinearSegmentedColormap.from_list('custom_diverging', colors, N=n_bins_cmap)
        
        # Plot heat map
        sns.heatmap(heatmap_data, 
                   cmap=cmap, 
                   center=0,
                   cbar_kws={'label': 'Log2 Fold Change'},
                   xticklabels=100,
                   yticklabels=False)
        
        plt.title(f'Genomic Expression Heat Map - {output_prefix}', fontsize=16, pad=20)
        plt.xlabel('Genome Position (kb)', fontsize=12)
        plt.ylabel('Expression Level', fontsize=12)
        
        # Add genome position labels
        tick_positions = np.linspace(0, len(bin_values), 10, dtype=int)
        tick_labels = [f"{i*bin_size//1000}" for i in tick_positions]
        plt.xticks(tick_positions, tick_labels, rotation=0)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_heatmap.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    def _create_barplot(self, gene_positions: List[Dict], output_prefix: str):
        """Create bar plot visualization."""
        if not gene_positions:
            print("No gene positions found for bar plot")
            return
        
        # Sort by genome position
        gene_positions.sort(key=lambda x: x['start'])
        
        # Create figure
        plt.figure(figsize=(20, 10))
        
        # Separate up and down regulated genes
        up_genes = [g for g in gene_positions if g['log2fc'] > 0]
        down_genes = [g for g in gene_positions if g['log2fc'] < 0]
        
        # Plot bars
        if up_genes:
            starts = [g['start'] for g in up_genes]
            log2fc_values = [g['log2fc'] for g in up_genes]
            plt.bar(starts, log2fc_values, width=1000, color='red', alpha=0.7, label='Up-regulated')
        
        if down_genes:
            starts = [g['start'] for g in down_genes]
            log2fc_values = [g['log2fc'] for g in down_genes]
            plt.bar(starts, log2fc_values, width=1000, color='blue', alpha=0.7, label='Down-regulated')
        
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        plt.xlabel('Genome Position (bp)', fontsize=12)
        plt.ylabel('Log2 Fold Change', fontsize=12)
        plt.title(f'Genomic Expression Bar Plot - {output_prefix}', fontsize=16)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Add some gene labels (every 10th gene to avoid overcrowding)
        for i, gene in enumerate(gene_positions[::10]):
            plt.annotate(gene['gene'], 
                        xy=(gene['start'], gene['log2fc']),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8, rotation=45, alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_barplot.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    def create_comparative_visualization(self, data_files: List[str], output_prefix: str):
        """
        Create comparative visualization of multiple datasets.
        
        Args:
            data_files: List of CSV file paths
            output_prefix: Prefix for output files
        """
        # Load all datasets
        datasets = {}
        for file_path in data_files:
            df = pd.read_csv(file_path)
            dataset_name = os.path.basename(file_path).split('.')[0]
            datasets[dataset_name] = df
        
        # Create comparative heat map
        self._create_comparative_heatmap(datasets, output_prefix)
    
    def _create_comparative_heatmap(self, datasets: Dict[str, pd.DataFrame], output_prefix: str):
        """Create comparative heat map of multiple datasets."""
        # Prepare data for comparison
        comparison_data = []
        dataset_names = []
        
        for name, df in datasets.items():
            # Find log2fc column
            log2fc_col = None
            for col in df.columns:
                if 'log2FoldChange' in col:
                    log2fc_col = col
                    break
            
            if log2fc_col:
                # Get gene coordinates and log2fc values
                gene_coords = {}
                for _, row in df.iterrows():
                    gene_name = str(row['geneName'])
                    if gene_name != 'nan' and gene_name != 'NA':
                        coords = self._get_gene_coordinates(gene_name)
                        if coords:
                            gene_coords[gene_name] = coords
                
                # Create genome array for this dataset
                genome_array = np.zeros(self.genome_length)
                for _, row in df.iterrows():
                    gene_name = str(row['geneName'])
                    if gene_name in gene_coords:
                        start, end, strand = gene_coords[gene_name]
                        log2fc = row[log2fc_col]
                        if not pd.isna(log2fc):
                            for pos in range(start, end + 1):
                                if 0 <= pos < len(genome_array):
                                    genome_array[pos] = log2fc
                
                # Bin the data
                bin_size = 10000
                bins = self._create_genome_bins(bin_size)
                bin_values = []
                
                for start, end in bins:
                    bin_data = genome_array[start-1:end]
                    non_zero = bin_data[bin_data != 0]
                    if len(non_zero) > 0:
                        bin_values.append(np.mean(non_zero))
                    else:
                        bin_values.append(0)
                
                comparison_data.append(bin_values)
                dataset_names.append(name)
        
        # Create comparative heat map
        if comparison_data:
            plt.figure(figsize=(20, 8))
            
            # Create 2D array for heat map
            heatmap_data = np.array(comparison_data)
            
            # Create custom colormap
            colors = ['darkblue', 'blue', 'lightblue', 'white', 'lightcoral', 'red', 'darkred']
            cmap = plt.cm.LinearSegmentedColormap.from_list('custom_diverging', colors, N=256)
            
            # Plot heat map
            sns.heatmap(heatmap_data, 
                       cmap=cmap, 
                       center=0,
                       cbar_kws={'label': 'Log2 Fold Change'},
                       xticklabels=100,
                       yticklabels=dataset_names)
            
            plt.title(f'Comparative Genomic Expression Heat Map', fontsize=16, pad=20)
            plt.xlabel('Genome Position (kb)', fontsize=12)
            plt.ylabel('Dataset', fontsize=12)
            
            # Add genome position labels
            tick_positions = np.linspace(0, heatmap_data.shape[1], 10, dtype=int)
            tick_labels = [f"{i*10000//1000}" for i in tick_positions]
            plt.xticks(tick_positions, tick_labels, rotation=0)
            
            plt.tight_layout()
            plt.savefig(f'{output_prefix}_comparative_heatmap.png', dpi=300, bbox_inches='tight')
            plt.show()

def main():
    """Main function to create visualizations for all three datasets."""
    # Initialize visualizer
    gene_data_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "all_gene_data.json")
    visualizer = GenomicVisualizer(gene_data_file)
    
    # Define data files
    data_files = {
        'D1-3': 'Ca_D1_3_6_data_with_Regulators.csv',
        'D6': 'Ca_D1_3_6_data_with_Regulators.csv', 
        'CueR': 'mmc4_CueR_with_Regulators.csv'
    }
    
    # Create individual visualizations
    for dataset_name, file_path in data_files.items():
        print(f"Creating visualizations for {dataset_name}...")
        
        # For D1-3 and D6, we need to specify different columns
        if dataset_name == 'D1-3':
            # Create a temporary file with only D1-3 data
            df = pd.read_csv(file_path)
            temp_file = f'temp_D1-3.csv'
            df_temp = df[['geneName', 'log2FoldChange_D1-3']].copy()
            df_temp.columns = ['geneName', 'log2FoldChange']
            df_temp.to_csv(temp_file, index=False)
            visualizer.create_heatmap_visualization(temp_file, f'D1-3_visualization', plot_type='both')
            os.remove(temp_file)
            
        elif dataset_name == 'D6':
            # Create a temporary file with only D6 data
            df = pd.read_csv(file_path)
            temp_file = f'temp_D6.csv'
            df_temp = df[['geneName', 'log2FoldChangeD-6']].copy()
            df_temp.columns = ['geneName', 'log2FoldChange']
            df_temp.to_csv(temp_file, index=False)
            visualizer.create_heatmap_visualization(temp_file, f'D6_visualization', plot_type='both')
            os.remove(temp_file)
            
        else:  # CueR
            visualizer.create_heatmap_visualization(file_path, f'CueR_visualization', plot_type='both')
    
    # Create comparative visualization
    print("Creating comparative visualization...")
    visualizer.create_comparative_visualization(list(data_files.values()), 'comparative')

if __name__ == "__main__":
    main()



