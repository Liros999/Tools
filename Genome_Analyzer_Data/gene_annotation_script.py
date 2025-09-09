#!/usr/bin/env python3
"""
Gene Annotation Analysis Script
Main script for analyzing BLAST results with gene annotations using Ensembl API
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import List, Dict, Any

# Add src directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from core.gene_annotator import GeneAnnotator
from core.data_processor import DataProcessor
from utils.error_handler import GenomeAnalysisError, handle_exception

def setup_logging(log_dir: str = "logs", log_level: str = "INFO") -> None:
    """
    Set up comprehensive logging system.
    
    Args:
        log_dir: Directory for log files
        log_level: Logging level
    """
    # Create logs directory
    log_path = Path(log_dir)
    log_path.mkdir(exist_ok=True)
    
    # Configure logging format
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    
    # File handler for all logs
    file_handler = logging.FileHandler(log_path / "gene_annotation.log")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(log_format))
    
    # File handler for errors
    error_handler = logging.FileHandler(log_path / "error.log")
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(logging.Formatter(log_format))
    
    # Console handler for important messages
    console_handler = logging.StreamHandler()
    console_handler.setLevel(getattr(logging, log_level))
    console_handler.setFormatter(logging.Formatter(log_format))
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    root_logger.addHandler(file_handler)
    root_logger.addHandler(error_handler)
    root_logger.addHandler(console_handler)
    
    logging.info("Logging system initialized successfully")

def progress_callback(current: int, total: int) -> None:
    """
    Progress callback for annotation progress.
    
    Args:
        current: Current row number
        total: Total number of rows
    """
    if current % 10 == 0 or current == total:
        percentage = (current / total) * 100
        print(f"Progress: {current}/{total} ({percentage:.1f}%)")

@handle_exception
def process_single_file(file_path: str, output_dir: str = ".", 
                       required_columns: List[str] = None) -> Dict[str, Any]:
    """
    Process a single CSV file for gene annotation.
    
    Args:
        file_path: Path to CSV file
        output_dir: Output directory
        required_columns: List of required column types
        
    Returns:
        Processing results dictionary
    """
    if required_columns is None:
        required_columns = ['accession', 'start', 'end']
    
    logging.info(f"Processing file: {file_path}")
    
    # Initialize components
    data_processor = DataProcessor()
    gene_annotator = GeneAnnotator()
    
    # Test API connection
    if not gene_annotator.test_api_connection():
        raise GenomeAnalysisError(
            error_type="API_CONNECTION_ERROR",
            message="Failed to connect to Ensembl API",
            citation="Ensembl REST API Documentation - Connection Requirements",
            recommended_next_step="Check internet connection and API server status"
        )
    
    # Process CSV file
    df, column_mapping = data_processor.process_csv_file(file_path, required_columns)
    
    # Get data summary
    data_summary = data_processor.get_data_summary(df, column_mapping)
    logging.info(f"Data summary: {data_summary}")
    
    # Annotate genes
    logging.info("Starting gene annotation...")
    gene_annotator.annotate_dataframe(
        df, 
        column_mapping['accession'],
        column_mapping['start'],
        column_mapping['end'],
        progress_callback
    )
    
    # Save annotated data
    output_path = data_processor.save_annotated_data(df, file_path, output_dir)
    
    # Get annotation statistics
    annotation_stats = gene_annotator.get_annotation_statistics()
    
    return {
        "file_path": file_path,
        "output_path": output_path,
        "data_summary": data_summary,
        "annotation_stats": annotation_stats,
        "total_rows": len(df)
    }

@handle_exception
def process_multiple_files(file_paths: List[str], output_dir: str = ".", 
                          required_columns: List[str] = None) -> List[Dict[str, Any]]:
    """
    Process multiple CSV files for gene annotation.
    
    Args:
        file_paths: List of CSV file paths
        output_dir: Output directory
        required_columns: List of required column types
        
    Returns:
        List of processing results
    """
    results = []
    
    for i, file_path in enumerate(file_paths, 1):
        try:
            print(f"\n--- Processing file {i}/{len(file_paths)}: {os.path.basename(file_path)} ---")
            
            result = process_single_file(file_path, output_dir, required_columns)
            results.append(result)
            
            print(f"✓ Successfully processed: {os.path.basename(file_path)}")
            print(f"  Output: {os.path.basename(result['output_path'])}")
            print(f"  Rows: {result['total_rows']}")
            
        except Exception as e:
            logging.error(f"Failed to process {file_path}: {e}")
            print(f"✗ Failed to process: {os.path.basename(file_path)}")
            print(f"  Error: {str(e)}")
            
            # Add error result
            results.append({
                "file_path": file_path,
                "error": str(e),
                "status": "failed"
            })
    
    return results

def print_summary(results: List[Dict[str, Any]]) -> None:
    """
    Print processing summary.
    
    Args:
        results: List of processing results
    """
    print("\n" + "="*60)
    print("PROCESSING SUMMARY")
    print("="*60)
    
    successful = [r for r in results if "error" not in r]
    failed = [r for r in results if "error" in r]
    
    print(f"Total files: {len(results)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")
    
    if successful:
        total_rows = sum(r.get('total_rows', 0) for r in successful)
        print(f"Total rows processed: {total_rows}")
        
        print("\nSuccessful files:")
        for result in successful:
            print(f"  ✓ {os.path.basename(result['file_path'])} -> {os.path.basename(result['output_path'])}")
    
    if failed:
        print("\nFailed files:")
        for result in failed:
            print(f"  ✗ {os.path.basename(result['file_path'])}: {result['error']}")

def main():
    """Main application entry point."""
    parser = argparse.ArgumentParser(
        description="Gene Annotation Analysis - BLAST Results with Ensembl Integration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gene_annotation_script.py --file blast_results.csv
  python gene_annotation_script.py --directory . --pattern "*.csv"
  python gene_annotation_script.py --file file1.csv file2.csv --output results/
  python gene_annotation_script.py --file data/*.csv --log-level DEBUG
        """
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--file", 
        nargs="+", 
        help="CSV file(s) to process"
    )
    
    input_group.add_argument(
        "--directory", 
        help="Directory containing CSV files to process"
    )
    
    input_group.add_argument(
        "--pattern", 
        default="*.csv",
        help="File pattern for directory processing (default: *.csv)"
    )
    
    # Output options
    parser.add_argument(
        "--output", 
        default=".",
        help="Output directory for annotated files (default: current directory)"
    )
    
    # Column options
    parser.add_argument(
        "--columns", 
        nargs="+",
        default=["accession", "start", "end"],
        help="Required column types (default: accession start end)"
    )
    
    # Logging options
    parser.add_argument(
        "--log-level", 
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)"
    )
    
    # Processing options
    parser.add_argument(
        "--test-api", 
        action="store_true",
        help="Test Ensembl API connection and exit"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(log_level=args.log_level)
    
    # Test API connection if requested
    if args.test_api:
        print("Testing Ensembl API connection...")
        annotator = GeneAnnotator()
        if annotator.test_api_connection():
            print("✓ API connection successful")
            return 0
        else:
            print("✗ API connection failed")
            return 1
    
    # Determine files to process
    file_paths = []
    
    if args.file:
        # Process specific files
        for file_path in args.file:
            if os.path.exists(file_path):
                file_paths.append(file_path)
            else:
                print(f"Warning: File not found: {file_path}")
    
    elif args.directory:
        # Process directory with pattern
        data_processor = DataProcessor()
        file_paths = data_processor.find_csv_files(args.directory, args.pattern)
    
    if not file_paths:
        print("Error: No valid CSV files found to process")
        return 1
    
    print(f"Found {len(file_paths)} CSV files to process")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    try:
        # Process files
        if len(file_paths) == 1:
            # Single file processing
            result = process_single_file(file_paths[0], args.output, args.columns)
            results = [result]
        else:
            # Multiple file processing
            results = process_multiple_files(file_paths, args.output, args.columns)
        
        # Print summary
        print_summary(results)
        
        return 0
        
    except GenomeAnalysisError as e:
        logging.error(f"Gene annotation error: {e.message}")
        print(f"Error: {e.message}")
        print(f"Recommended next step: {e.recommended_next_step}")
        return 1
        
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        print(f"Unexpected error: {e}")
        return 1

if __name__ == "__main__":
    exit(main())
