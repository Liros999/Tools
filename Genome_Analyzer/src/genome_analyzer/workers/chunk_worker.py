"""
Chunk worker functions for genome analysis parallel processing.
Handles genome chunk search operations and gene data validation.
"""

import os
from typing import List, Dict

from ..analyzer import GenomeAnalyzer


def _search_genome_chunk_worker(analyzer: GenomeAnalyzer, start: int, end: int, 
                               pattern: str, max_mismatches: int, boundary_bp: int) -> List:
    """
    FIXED: Properly bound working intergenic location detection method with gene data validation
    """
    try:
        # Create worker analyzer
        worker_analyzer = GenomeAnalyzer()
        worker_analyzer.sequence = analyzer.sequence[start:end]
        worker_analyzer.current_chromosome = getattr(analyzer, 'current_chromosome', None)
        
        # CRITICAL FIX: Ensure gene data is loaded and valid before passing to worker
        if hasattr(analyzer, 'gene_data') and analyzer.gene_data and len(analyzer.gene_data) > 0:
            worker_analyzer.gene_data = analyzer.gene_data
            print(f"[INFO] Worker {start}-{end}: Using gene data from main analyzer ({len(analyzer.gene_data)} genes)")

        else:
            print(f"[WARNING] Main analyzer has no gene data! Gene count: {len(getattr(analyzer, 'gene_data', {}))}")
            # Try to load gene data directly if main analyzer doesn't have it
            try:
                # Try to detect species and load appropriate gene data
                if hasattr(analyzer, 'current_chromosome') and analyzer.current_chromosome:
                    if 'B.sub' in analyzer.current_chromosome or 'subtilis' in analyzer.current_chromosome:
                        # SIMPLIFIED: Use direct path to known B.Sub gene data file
                        gene_data_path = "data/B.Sub/all_gene_data.json"  # FIXED: Correct relative path
                        if os.path.exists(gene_data_path):
                            import json
                            with open(gene_data_path, 'r') as f:
                                worker_analyzer.gene_data = json.load(f)
                                print(f"[INFO] Worker {start}-{end}: Loaded gene data from file ({len(worker_analyzer.gene_data)} genes)")

                        else:
                            print(f"[ERROR] B.Sub gene data file not found at {gene_data_path}")
                            return []
                    elif 'E.coli' in analyzer.current_chromosome or 'coli' in analyzer.current_chromosome:
                        # Load E.Coli gene data - handle different file format
                        gene_data_path = "data/E.Coli/Gene_Coordiantes.txt"  # FIXED: Correct filename
                        if os.path.exists(gene_data_path):
                            # E.Coli uses a different format - need to parse it
                            try:
                                genes = {}
                                with open(gene_data_path, 'r') as f:
                                    for line in f:
                                        if line.startswith('#'):
                                            continue
                                        parts = line.strip().split('\t')
                                        if len(parts) >= 3:
                                            gene_id = parts[0]
                                            start = int(parts[1])
                                            end = int(parts[2])
                                            genes[gene_id] = {'start': start, 'end': end}
                                worker_analyzer.gene_data = genes

                            except Exception as e:
                                print(f"[ERROR] Failed to parse E.Coli gene file: {e}")
                                return []
                        else:
                            print(f"[ERROR] E.Coli gene data file not found at {gene_data_path}")
                            return []
                    elif 'Homo' in analyzer.current_chromosome or 'sapiens' in analyzer.current_chromosome:
                        # Load Homo Sapiens gene data - handle GTF format
                        gtf_path = "data/Homo_Sapiens/gencode.v48.primary_assembly.annotation.gtf"
                        if os.path.exists(gtf_path):
                            try:
                                # Parse GTF file for genes
                                genes = {}
                                with open(gtf_path, 'r') as f:
                                    for line in f:
                                        if line.startswith('#'):
                                            continue
                                        parts = line.strip().split('\t')
                                        if len(parts) >= 9 and parts[2] == 'gene':
                                            start = int(parts[3])
                                            end = int(parts[4])
                                            # Extract gene ID from attributes
                                            attrs = parts[8]
                                            gene_id = None
                                            for attr in attrs.split(';'):
                                                if 'gene_id' in attr:
                                                    gene_id = attr.split('"')[1]
                                                    break
                                            if gene_id:
                                                genes[gene_id] = {'start': start, 'end': end}
                                worker_analyzer.gene_data = genes

                            except Exception as e:
                                print(f"[ERROR] Failed to parse Homo Sapiens GTF file: {e}")
                                return []
                        else:
                            print(f"[ERROR] Homo Sapiens GTF file not found at {gtf_path}")
                            return []
                    else:
                        print(f"[ERROR] No gene data available for chromosome: {analyzer.current_chromosome}")
                        return []
                else:
                    print(f"[ERROR] Cannot determine species for gene data loading")
                    return []
            except Exception as e:
                print(f"[ERROR] Failed to load gene data: {e}")
                return []
        
        # CRITICAL FIX: Properly bind the working intergenic location method
        if hasattr(analyzer, '_get_intergenic_location_working'):
            # Bind the method to the worker analyzer instance
            worker_analyzer._get_intergenic_location = analyzer._get_intergenic_location_working.__get__(worker_analyzer, GenomeAnalyzer)

        elif hasattr(analyzer, '_get_intergenic_location'):
            # Bind the method to the worker analyzer instance
            worker_analyzer._get_intergenic_location = analyzer._get_intergenic_location.__get__(worker_analyzer, GenomeAnalyzer)

        
        # Validate that worker analyzer has gene data before search
        if not hasattr(worker_analyzer, 'gene_data') or not worker_analyzer.gene_data:
            print(f"[ERROR] Worker analyzer has no gene data after setup!")
            return []
        

        
        # Perform search
        results = worker_analyzer.search_sequence(pattern, max_mismatches, boundary_bp, stream=False)
        
        # Adjust positions for chunk offset
        adjusted_results = []
        for result in results:
            try:
                if len(result) >= 8:  # E.coli/Generic format (8 values)
                    start_pos, end_pos, strand, pat, found, revcomp, mismatches, location = result
                    adjusted_results.append((start_pos + start, end_pos + start, strand, pat, found, revcomp, mismatches, location))
                elif len(result) >= 7:  # B.sub format (7 values) - CONVERT TO 8-VALUE FORMAT
                    start_pos, end_pos, strand, pat, found, location, mismatches = result
                    # Add reverse complement to match 8-value format
                    revcomp = found  # Use found sequence as revcomp for now
                    adjusted_results.append((start_pos + start, end_pos + start, strand, pat, found, revcomp, mismatches, location))
                else:
                    print(f"[WARNING] Skipping malformed result: {result}")
                    continue
            except Exception as e:
                print(f"[ERROR] Failed to process result {result}: {e}")
                continue
        
        return adjusted_results
        
    except Exception as e:
        print(f"[ERROR] Chunk {start}-{end}: {str(e)}")
        return []


def _search_genome_chunk_worker_simple(analyzer: GenomeAnalyzer, start: int, end: int, 
                                      pattern: str, max_mismatches: int, boundary_bp: int) -> List:
    """
    Simplified worker for Homo Sapiens genome search
    """
    try:
        # Use the main analyzer directly for simple search
        results = analyzer.search_sequence(pattern, max_mismatches, boundary_bp)
        
        # Filter results to this chunk
        chunk_results = []
        for result in results:
            if start <= result[0] <= end:  # result[0] is start position
                chunk_results.append(result)
        
        return chunk_results
        
    except Exception as e:
        print(f"[ERROR] Worker error: {e}")
        return []


# Export functions
__all__ = [
    '_search_genome_chunk_worker',
    '_search_genome_chunk_worker_simple'
]
