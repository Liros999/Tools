"""
Chromosome worker functions for genome analysis parallel processing.
Handles individual chromosome search operations with timeout and progress tracking.
"""

import os
import time
import hashlib
import threading
import queue
from typing import Tuple, List, Dict, Optional

from ..analyzer import GenomeAnalyzer
from ..file_processing import load_annotations_from_shared_cache, load_genome_for_chromosome


def enhanced_search_single_chromosome_worker(chromosome: str, fasta_path: str, gtf_path: str, 
                                           pattern: str, max_mismatches: int, boundary_bp: int, 
                                           worker_id: int = 0, pre_loaded_genome: str = None, 
                                           chromosome_boundaries: Dict = None, pre_loaded_annotations: Dict = None) -> Tuple[bool, str, List]:
    """Enhanced worker function with proper progress reporting and result validation."""
    
    try:
        # Create new analyzer instance for this worker (multiprocessing-safe)
        analyzer = GenomeAnalyzer()
        
        # Minimal worker status - dashboard handles detailed progress
        print(f"[WORKER] {chromosome}: Started")
        
        # OPTIMIZATION: Use pre-loaded genome if available, otherwise load individually
        if pre_loaded_genome and chromosome_boundaries and chromosome in chromosome_boundaries:
            # Extract chromosome sequence from pre-loaded genome using boundaries
            # This is much faster than loading from disk again
            boundaries = chromosome_boundaries[chromosome]
            start_pos = boundaries['start'] - 1  # Convert to 0-based indexing
            end_pos = boundaries['end']
            
            # OPTIMIZATION: Use slice operation for efficient sequence extraction
            chromosome_sequence = pre_loaded_genome[start_pos:end_pos]
            analyzer.sequence = chromosome_sequence
            analyzer.current_chromosome = chromosome
            print(f"[WORKER] {chromosome}: Extracted sequence from pre-loaded genome ({len(chromosome_sequence):,} bp)")
            
            # Worker progress tracking (no dashboard communication from worker processes)
            print(f"[WORKER] {chromosome}: Sequence extracted - 25% complete")
        else:
            # Fallback to individual loading (slower but memory-efficient)
            print(f"[WORKER] {chromosome}: Loading individual chromosome sequence...")
            sequence = load_genome_for_chromosome(fasta_path, chromosome, use_memory_mapping=True)
            analyzer.sequence = sequence
            analyzer.current_chromosome = chromosome
            
            # Worker progress tracking (no dashboard communication from worker processes)
            print(f"[WORKER] {chromosome}: Sequence loaded - 25% complete")
        
        # CRITICAL OPTIMIZATION: Use pre-loaded annotations for instant access
        if pre_loaded_annotations and chromosome in pre_loaded_annotations:
            print(f"[WORKER] {chromosome}: Loading pre-loaded annotations...")
            
            # Worker progress tracking (no dashboard communication from worker processes)
            print(f"[WORKER] {chromosome}: Annotations loading - 50% complete")
            
            annotations = pre_loaded_annotations[chromosome]
            genes = annotations['genes']
            cds = annotations['cds']
            
            # Build gene tree instantly (no file I/O)
            analyzer._build_gene_tree_with_intergenic(genes)
            for cds_item in cds:
                analyzer.cds_tree[cds_item['start']:cds_item['end'] + 1] = cds_item
            
            print(f"[WORKER] {chromosome}: Annotations loaded instantly ({len(genes)} genes, {len(cds)} CDS)")
            
            # Worker progress tracking (no dashboard communication from worker processes)
            print(f"[WORKER] {chromosome}: Annotations ready - 75% complete")
        else:
            # Fallback to loading from GTF (slower)
            print(f"[WORKER] {chromosome}: Loading annotations from GTF (fallback)...")
            
            # Worker progress tracking (no dashboard communication from worker processes)
            print(f"[WORKER] {chromosome}: GTF loading - 50% complete")
            
            gtf_stat = os.stat(gtf_path)
            cache_key = f"{gtf_path}_{gtf_stat.st_mtime}_{gtf_stat.st_size}"
            cache_hash = hashlib.md5(cache_key.encode()).hexdigest()
            cache_dir = os.path.join(os.path.dirname(__file__), "cache")
            cache_file = os.path.join(cache_dir, f"gtf_cache_{cache_hash}.pkl")
            
            if os.path.exists(cache_file):
                genes, cds = load_annotations_from_shared_cache(cache_file, chromosome)
                analyzer._build_gene_tree_with_intergenic(genes)
                for cds_item in cds:
                    analyzer.cds_tree[cds_item['start']:cds_item['end'] + 1] = cds_item
            else:
                from ..file_processing import load_annotations_from_gtf
                load_annotations_from_gtf(analyzer, gtf_path, chromosome)
        
        # Execute search with real-time progress reporting
        print(f"[WORKER] {chromosome}: Starting pattern search...")
        
        # Worker progress tracking (no dashboard communication from worker processes)
        print(f"[WORKER] {chromosome}: Search started - 90% complete")
        
        # Start search with progress monitoring and timeout
        start_time = time.time()
        
        # Add timeout for large chromosomes (chr8 is ~145MB, allow 5 minutes max)
        max_search_time = 300  # 5 minutes timeout
        search_timeout = False
        
        try:
            # Use a separate thread for search with timeout
            result_queue = queue.Queue()
            
            def search_with_timeout():
                try:
                    search_results = analyzer.search_sequence(pattern, max_mismatches, boundary_bp, stream=False)
                    result_queue.put(('success', search_results))
                except Exception as e:
                    result_queue.put(('error', str(e)))
            
            search_thread = threading.Thread(target=search_with_timeout)
            search_thread.daemon = True
            search_thread.start()
            
            # Wait for search with timeout
            search_thread.join(timeout=max_search_time)
            
            if search_thread.is_alive():
                print(f"[WARNING] {chromosome}: Search timeout after {max_search_time}s - skipping")
                search_timeout = True
                results = []
            else:
                # Get results from queue
                try:
                    result_type, results = result_queue.get_nowait()
                    if result_type == 'error':
                        print(f"[ERROR] {chromosome}: Search failed: {results}")
                        results = []
                except queue.Empty:
                    print(f"[WARNING] {chromosome}: No results from search thread")
                    results = []
            
        except Exception as e:
            print(f"[ERROR] {chromosome}: Timeout handling failed: {e}")
            results = []
            search_timeout = True
        
        search_duration = time.time() - start_time
        
        # Validate results
        valid_results = []
        for result in results:
            if isinstance(result, (tuple, list)) and len(result) >= 8:
                valid_results.append(result)
        
        if search_timeout:
            print(f"[WORKER] {chromosome}: Search TIMEOUT ({search_duration:.2f}s) - 0 matches returned")
        else:
            print(f"[WORKER] {chromosome}: Search completed ({len(valid_results)} matches) in {search_duration:.2f}s")
        
        # Worker progress tracking (no dashboard communication from worker processes)
        print(f"[WORKER] {chromosome}: Search complete - 100% complete")
        
        return True, chromosome, valid_results
        
    except Exception as e:
        print(f"[ERROR] Worker {chromosome}: {str(e)}")
        return False, chromosome, []


# Export function
__all__ = ['enhanced_search_single_chromosome_worker']
