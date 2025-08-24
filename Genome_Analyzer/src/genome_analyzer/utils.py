"""
Utility functions for genome analysis operations.
Contains helper functions for file operations, memory monitoring, and data processing.
"""

import os
import psutil
from typing import List, Tuple, Optional, Dict, Any
from .analyzer import GenomeAnalysisError

# Constants from main file
MEMORY_MAPPING_THRESHOLD_BYTES = 64 * 4096  # 64 pages = 256KB
PROGRESS_BAR_SEGMENTS = 20  # Number of segments in progress bars
CACHE_MAX_RETRIES = 5  # Maximum retries for cache file creation
CACHE_RETRY_DELAY = 0.5  # Delay between cache creation retries (seconds)

def _monitor_memory_usage(initial_memory: float = None, threshold_mb: int = 512) -> tuple:
    """
    Centralized memory monitoring function to eliminate code duplication.
    Returns (current_memory_mb, memory_increase_mb, should_cleanup).
    
    Scientific basis for 512MB threshold: Based on typical RAM allocation for
    bioinformatics applications (Ref: Bioinformatics algorithms design patterns, 2013).
    """
    try:
        current_memory = psutil.Process().memory_info().rss / (1024 * 1024)
        
        if initial_memory is None:
            return current_memory, 0, False
        
        memory_increase = current_memory - initial_memory
        should_cleanup = memory_increase > threshold_mb
        
        return current_memory, memory_increase, should_cleanup
        
    except ImportError:
        return 0, 0, False

def _handle_file_error(error: Exception, file_path: str, operation: str) -> None:
    """
    Centralized file error handling to eliminate code duplication.
    Raises appropriate GenomeAnalysisError with consistent messaging.
    """
    if isinstance(error, FileNotFoundError):
        raise GenomeAnalysisError(f"{operation} file not found at {file_path}")
    else:
        raise GenomeAnalysisError(f"An error occurred during {operation}: {error}")

def _print_worker_progress(worker_line: int, message: str, progress: int = 100) -> None:
    """
    Centralized worker progress reporting to eliminate code duplication.
    """
    progress_bar = "█" * (progress * PROGRESS_BAR_SEGMENTS // 100)
    print(f"\033[{worker_line+1};0H\033[K[WORKER] {message} | Progress: {progress_bar} {progress:.1f}%")

def _list_subdirectories(path: str) -> List[str]:
    """List subdirectories in the given path."""
    try:
        all_dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
        # FILTER: Only show Homo_Sapiens, remove B.Sub and E.Coli
        filtered_dirs = [d for d in all_dirs if d == 'Homo_Sapiens']
        return filtered_dirs
    except Exception:
        return []

def _list_files_with_ext(path: str, exts: Tuple[str, ...]) -> List[str]:
    """List files with specified extensions in the given path."""
    files = []
    try:
        for name in os.listdir(path):
            full = os.path.join(path, name)
            if os.path.isfile(full) and any(name.lower().endswith(e) for e in exts):
                files.append(name)
    except Exception:
        pass
    return sorted(files)

def _discover_gtf_chromosomes(gtf_path: str, max_to_show: int = 200) -> List[str]:
    """Discover chromosomes from GTF file."""
    chroms = []
    seen = set()
    malformed_line_count = 0
    total_lines_processed = 0
    
    try:
        with open(gtf_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                total_lines_processed += 1
                if not line or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    malformed_line_count += 1
                    continue
                seqid = parts[0]
                if seqid in seen:
                    continue
                seen.add(seqid)
                chroms.append(seqid)
                if len(chroms) >= max_to_show:
                    break

        # Check for data integrity issues
        if malformed_line_count > 0:
            print(f"\n[WARNING] DATA INTEGRITY ISSUE DETECTED:")
            print(f"[WARNING] Ignored {malformed_line_count:,} malformed lines during chromosome discovery")
            print(f"[WARNING] Total lines processed: {total_lines_processed:,}")
            print(f"[WARNING] Malformed line ratio: {(malformed_line_count/total_lines_processed)*100:.2f}%")
            print(f"[WARNING] The discovered chromosome list may be incomplete due to input file corruption")
            print(f"[WARNING] Discovered chromosomes: {len(chroms)}")
            print(f"[WARNING] This could lead to incomplete analysis results!\n")
            
    except Exception as e:
        print(f"[ERROR] Failed to discover chromosomes from {gtf_path}: {e}")
        return []
    
    return chroms

def _discover_fasta_headers(fasta_path: str, max_headers: int = 100000) -> List[str]:
    """Discover FASTA headers from file."""
    headers: List[str] = []
    try:
        with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('>'):
                    headers.append(line[1:].split()[0])
                    if len(headers) >= max_headers:
                        break
    except Exception:
        return []
    return headers

def _filter_human_chromosomes(chroms: List[str], fasta_path: str) -> List[str]:
    """Filter to canonical human chromosomes."""
    fasta_headers = set(_discover_fasta_headers(fasta_path))
    canonical: List[str] = []
    for n in list(map(str, range(1, 23))) + ['X', 'Y', 'MT']:
        candidates = {n, 'chr' + n}
        if n == 'MT':
            candidates.update({'M', 'chrM'})
        for c in candidates:
            if not fasta_headers or c in fasta_headers:
                canonical.append(c)
                break
    allowed = set(canonical) if canonical else set(chroms)
    return [c for c in chroms if c in allowed]

def _get_default_genome_file(species_name: str) -> Optional[str]:
    """Get default genome file for each species."""
    defaults = {
        'Homo_Sapiens': 'GRCh38.primary_assembly.genome.fa'
    }
    return defaults.get(species_name)

def _get_default_annotation_file(species_name: str) -> Optional[str]:
    """Get default annotation file for each species."""
    defaults = {
        'Homo_Sapiens': 'gencode.v48.primary_assembly.annotation.gtf'
    }
    return defaults.get(species_name)

def _get_filtered_chromosomes(gtf_path: str, fasta_path: str) -> List[str]:
    """
    Centralized chromosome discovery and filtering to eliminate code duplication.
    Returns filtered list of chromosomes for parallel processing.
    """
    chroms = _discover_gtf_chromosomes(gtf_path)
    filtered_chroms = _filter_human_chromosomes(chroms, fasta_path)
    print(f"[INFO] Using parallel processing for {len(filtered_chroms)} chromosomes")
    return filtered_chroms

def _perform_cache_cleanup(cache_loader, operation_name: str = "operation"):
    """
    Centralized cache cleanup to eliminate code duplication.
    """
    try:
        cache_loader.cleanup_cache()
        print(f"[INFO] Cache cleanup completed after {operation_name}")
    except Exception as e:
        print(f"[WARNING] Cache cleanup failed after {operation_name}: {e}")

def _process_search_results(matches: List[Tuple], pattern: str, strand: str, 
                          analyzer, boundary_bp: int, 
                          on_row: Optional[callable] = None, stream: bool = False,
                          include_chrom: bool = False) -> List[Tuple]:
    """
    Centralized result processing to eliminate code duplication.
    Returns processed results list.
    """
    processed_results = []
    for start, end, window, _, mismatches, _ in matches:
        location = analyzer.get_location_details(start, end)
        
        # Apply boundary annotations if specified
        if boundary_bp and boundary_bp > 0:
            seq_len = len(analyzer.sequence)
            dist_start = start - 1
            dist_end = seq_len - end
            near_flags = []
            if dist_start <= boundary_bp:
                near_flags.append(f"NearStart:{dist_start}bp")
            if dist_end <= boundary_bp:
                near_flags.append(f"NearEnd:{dist_end}bp")
            if near_flags:
                location = f"{location} [{' & '.join(near_flags)}]"
        
        # Create result row
        row = (start, end, strand, pattern, window, analyzer._reverse_complement(window), mismatches, location)
        processed_results.append(row)
        
        # Handle output based on parameters
        if on_row:
            on_row(row)
        elif stream:
            # Import here to avoid circular imports
            from .ui import _print_row_tuple
            _print_row_tuple(row, include_chrom)
    
    return processed_results

def _handle_sigint(sig, frame):
    """Handle SIGINT signal for graceful exit."""
    try:
        print("\nOperation cancelled by user. Exiting gracefully.")
    finally:
        import sys
        sys.exit(0)
