#!/usr/bin/env python3
"""
DNAzyme Gene Annotation Analyzer
================================

This script analyzes DNAzyme BLAST results from Homo sapiens and annotates
genomic regions with gene information using the Ensembl REST API.

Scientific Basis:
- Uses Ensembl REST API v1 for genomic region overlap analysis
- Implements proper rate limiting (100ms delay between requests)
- Processes BLAST results with NC_0000XX.X chromosome identifiers
- Annotates gene overlaps with scientific precision

Author: AI Assistant
Date: 2024
"""

import os
import glob
import pandas as pd
import requests
import time
import re
import json
from datetime import datetime
from pathlib import Path

# --- Configuration Constants ---
REQUEST_DELAY = 0.1  # 100ms delay between API requests (scientifically validated rate limiting)
ENSEMBL_BASE_URL = "https://rest.ensembl.org"
CHROMOSOME_PATTERN = r'NC_0000(\d{2})'  # Pattern for standard human chromosome identifiers

# --- Logging Configuration ---
LOG_DIR = Path("logs")
LOG_DIR.mkdir(exist_ok=True)
LOG_FILE = LOG_DIR / f"dnazyme_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

def setup_logging():
    """Initialize logging to file with comprehensive analysis tracking."""
    import logging
    
    # Create formatters
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')
    
    # File handler
    file_handler = logging.FileHandler(LOG_FILE)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(file_formatter)
    
    # Console handler with reduced verbosity
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.WARNING)  # Only show warnings and errors on console
    console_handler.setFormatter(console_formatter)
    
    # Setup logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger

logger = setup_logging()

def cleanup_old_files():
    """Remove old annotated files from previous runs to prevent duplicates."""
    logger.info("Cleaning up old result files...")
    old_files = glob.glob('*_annotated.csv')
    if old_files:
        for f in old_files:
            os.remove(f)
            logger.info(f"Removed old file: {f}")
        logger.info(f"Removed {len(old_files)} old annotated files.")
    else:
        logger.info("No old files to clean up.")

def find_column_name(df, potential_names):
    """
    Find column names in dataframe from a list of potential names.
    
    Args:
        df (pd.DataFrame): Input dataframe
        potential_names (list): List of potential column names to search for
        
    Returns:
        str or None: First matching column name or None if not found
    """
    for name in potential_names:
        if name in df.columns:
            return name
    return None

def get_single_annotation(region_str, logger):
    """
    Fetch gene annotations for a single genomic region using Ensembl REST API.
    
    Args:
        region_str (str): Genomic region in format "chromosome:start-end"
        logger: Logger instance for tracking API calls
        
    Returns:
        str: Gene annotation or error indicator
    """
    endpoint = f"/overlap/region/human/{region_str}"
    url = f"{ENSEMBL_BASE_URL}{endpoint}"
    headers = {"Content-Type": "application/json"}
    params = {"feature": "gene"}
    
    try:
        response = requests.get(url, headers=headers, params=params, timeout=30)
        
        # Log rate limit information
        remaining = response.headers.get('X-RateLimit-Remaining')
        reset_time = response.headers.get('X-RateLimit-Reset')
        logger.debug(f"Rate limit remaining: {remaining}, Resets in: {reset_time}s")

        response.raise_for_status()
        results = response.json()

        if not results:
            return "Intergenic"

        # Extract gene names and handle multiple genes
        gene_names = sorted(list(set(g['external_name'] for g in results if g.get('external_name'))))
        if not gene_names:
            return "Unnamed Gene Region"
        
        if len(gene_names) == 1:
            return f"[{gene_names[0]}]Overlap"
        else:
            return f"[{gene_names[0]}]-------[{gene_names[1]}]"

    except requests.exceptions.HTTPError as e:
        logger.warning(f"HTTP Error for {region_str}: {e}")
        return "API_Error"
    except requests.exceptions.Timeout:
        logger.warning(f"Timeout for {region_str}")
        return "Timeout_Error"
    except Exception as e:
        logger.error(f"Unexpected error for {region_str}: {e}")
        return "Parsing_Error"

def process_csv_file(filepath, logger):
    """
    Process a single CSV file and annotate genomic regions.
    
    Args:
        filepath (str): Path to the CSV file to process
        logger: Logger instance for tracking progress
        
    Returns:
        bool: True if successful, False otherwise
    """
    filename = os.path.basename(filepath)
    logger.info(f"Processing file: {filename}")
    
    try:
        # Read CSV file
        df = pd.read_csv(filepath)
        logger.info(f"Loaded {df.shape[0]} rows from {filename}")
        
        # Find required columns
        col_names = {
            'accession': find_column_name(df, ['subject acc.ver']),
            'start': find_column_name(df, ['s. start']),
            'end': find_column_name(df, ['s. end', 's. End', 'Subject End'])
        }

        # Validate column presence
        if None in col_names.values():
            missing = [key for key, val in col_names.items() if val is None]
            logger.error(f"Missing required columns in {filename}: {', '.join(missing)}")
            return False

        logger.info("All required columns found. Preparing genomic regions...")
        
        # Convert coordinates to numeric and clean data
        df[col_names['start']] = pd.to_numeric(df[col_names['start']], errors='coerce')
        df[col_names['end']] = pd.to_numeric(df[col_names['end']], errors='coerce')
        df.dropna(subset=[col_names['start'], col_names['end']], inplace=True)
        
        # Convert to integer type for genomic coordinates
        df[col_names['start']] = df[col_names['start']].astype('Int64')
        df[col_names['end']] = df[col_names['end']].astype('Int64')
        
        # Initialize annotation column
        df['Gene_Annotation'] = "Queued"
        
        # Process each row for gene annotation
        valid_regions = 0
        total_regions = df.shape[0]
        
        # Create progress bar
        from tqdm import tqdm
        pbar = tqdm(total=total_regions, desc=f"Annotating {filename}", 
                   unit="regions", position=0, leave=True)
        
        for index, row in df.iterrows():
            accession = row[col_names['accession']]
            
            # Check if it's a standard human chromosome
            if re.search(CHROMOSOME_PATTERN, str(accession)):
                match = re.search(CHROMOSOME_PATTERN, str(accession))
                chrom_number = match.group(1).lstrip('0')
                start, end = row[col_names['start']], row[col_names['end']]
                
                # Ensure start < end for proper region formatting
                if start > end:
                    start, end = end, start
                
                region_str = f"{chrom_number}:{start}-{end}"
                annotation = get_single_annotation(region_str, logger)
                df.at[index, 'Gene_Annotation'] = annotation
                valid_regions += 1
                
                # Update progress bar
                pbar.update(1)
                pbar.set_postfix({
                    'Valid': valid_regions,
                    'Current': f"{chrom_number}:{start}-{end}",
                    'Annotation': annotation[:20] + "..." if len(annotation) > 20 else annotation
                })
                
                # Respectful delay between API requests
                time.sleep(REQUEST_DELAY)
            else:
                df.at[index, 'Gene_Annotation'] = "Non-standard Chromosome/Contig"
                pbar.update(1)
        
        # Close progress bar
        pbar.close()

        # Save annotated results
        output_filename = f"{os.path.splitext(filename)[0]}_annotated.csv"
        df.to_csv(output_filename, index=False)
        
        logger.info(f"Successfully annotated {valid_regions} regions")
        logger.info(f"Results saved to: {output_filename}")
        
        # Generate summary statistics
        annotation_counts = df['Gene_Annotation'].value_counts()
        logger.info("Annotation Summary:")
        for annotation, count in annotation_counts.items():
            logger.info(f"  {annotation}: {count}")
        
        return True
        
    except Exception as e:
        logger.error(f"Error processing {filename}: {e}")
        return False

def main():
    """Main execution function for DNAzyme gene annotation analysis."""
    print("=== DNAzyme Gene Annotation Analyzer ===")
    print(f"Analysis started at: {datetime.now()}")
    print(f"Request delay: {REQUEST_DELAY}s")
    print(f"Ensembl API base URL: {ENSEMBL_BASE_URL}")
    
    # Clean up old files
    cleanup_old_files()
    
    # Find CSV files to process
    csv_files = [f for f in glob.glob('*.csv') if '_annotated' not in f and 'Original' not in f]
    
    if not csv_files:
        print("ERROR: No CSV files found to process.")
        print("Expected files: mqsA_cured_Homo_Sapiens_Binding_Substrates_Off-Target_Blast_Results.csv")
        print("                MazE_cured_Homo_Sapiens_Binding_Substrates_Off-Target_Blast_Results.csv")
        return
    
    print(f"Found {len(csv_files)} CSV files to process:")
    for f in csv_files:
        print(f"  - {f}")
    print()  # Empty line before progress bars
    
    # Process each file
    successful_files = 0
    from tqdm import tqdm
    file_pbar = tqdm(csv_files, desc="Processing Files", unit="file", position=0, leave=True)
    
    for filepath in file_pbar:
        file_pbar.set_description(f"Processing {os.path.basename(filepath)}")
        if process_csv_file(filepath, logger):
            successful_files += 1
    
    file_pbar.close()
    
    # Final summary
    print("\n=== Analysis Complete ===")
    print(f"Successfully processed: {successful_files}/{len(csv_files)} files")
    print(f"Log file: {LOG_FILE}")
    
    if successful_files == len(csv_files):
        print("All files processed successfully!")
    else:
        print(f"WARNING: Some files failed to process. Check logs for details.")

if __name__ == "__main__":
    main()
