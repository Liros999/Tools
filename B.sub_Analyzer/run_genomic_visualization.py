#!/usr/bin/env python3
"""
Script to run genomic visualization for B.subtilis expression data.
"""

import os
import sys
import pandas as pd
from src.genomic_visualization import GenomicVisualizer

def main():
    """Main function to create visualizations for all three datasets."""
    
    # Get the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Define paths
    gene_data_file = os.path.join(current_dir, "data", "all_gene_data.json")
    volcano_app_dir = os.path.join(os.path.dirname(current_dir), "InteractiveVolcanoApp")
    
    # Check if gene data file exists
    if not os.path.exists(gene_data_file):
        print(f"Error: Gene data file not found at {gene_data_file}")
        return
    
    # Initialize visualizer
    visualizer = GenomicVisualizer(gene_data_file)
    
    # Define data files with full paths
    data_files = {
        'D1-3': os.path.join(volcano_app_dir, 'Ca_D1_3_6_data_with_Regulators.csv'),
        'D6': os.path.join(volcano_app_dir, 'Ca_D1_3_6_data_with_Regulators.csv'), 
        'CueR': os.path.join(volcano_app_dir, 'mmc4_CueR_with_Regulators.csv')
    }
    
    # Check if data files exist
    for name, file_path in data_files.items():
        if not os.path.exists(file_path):
            print(f"Error: Data file not found for {name}: {file_path}")
            return
    
    # Create output directory
    output_dir = os.path.join(current_dir, "results", "genomic_visualizations")
    os.makedirs(output_dir, exist_ok=True)
    
    print("Creating genomic visualizations...")
    print(f"Output directory: {output_dir}")
    
    # Create individual visualizations
    for dataset_name, file_path in data_files.items():
        print(f"\nCreating visualizations for {dataset_name}...")
        
        # For D1-3 and D6, we need to specify different columns
        if dataset_name == 'D1-3':
            # Create a temporary file with only D1-3 data
            df = pd.read_csv(file_path)
            temp_file = os.path.join(output_dir, f'temp_D1-3.csv')
            df_temp = df[['geneName', 'log2FoldChange_D1-3']].copy()
            df_temp.columns = ['geneName', 'log2FoldChange']
            df_temp.to_csv(temp_file, index=False)
            
            output_prefix = os.path.join(output_dir, 'D1-3_visualization')
            visualizer.create_heatmap_visualization(temp_file, output_prefix, plot_type='both')
            os.remove(temp_file)
            
        elif dataset_name == 'D6':
            # Create a temporary file with only D6 data
            df = pd.read_csv(file_path)
            temp_file = os.path.join(output_dir, f'temp_D6.csv')
            df_temp = df[['geneName', 'log2FoldChangeD-6']].copy()
            df_temp.columns = ['geneName', 'log2FoldChange']
            df_temp.to_csv(temp_file, index=False)
            
            output_prefix = os.path.join(output_dir, 'D6_visualization')
            visualizer.create_heatmap_visualization(temp_file, output_prefix, plot_type='both')
            os.remove(temp_file)
            
        else:  # CueR
            output_prefix = os.path.join(output_dir, 'CueR_visualization')
            visualizer.create_heatmap_visualization(file_path, output_prefix, plot_type='both')
    
    # Create comparative visualization
    print("\nCreating comparative visualization...")
    output_prefix = os.path.join(output_dir, 'comparative')
    visualizer.create_comparative_visualization(list(data_files.values()), output_prefix)
    
    print(f"\nAll visualizations completed! Check the output directory: {output_dir}")

if __name__ == "__main__":
    main()



