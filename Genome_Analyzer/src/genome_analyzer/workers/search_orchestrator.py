"""
Search orchestration for genome analysis parallel processing.
Handles the main parallel search coordination and worker distribution.
"""

import os
import threading
import time
from typing import List, Dict, Optional
from concurrent.futures import as_completed

import multiprocessing

from .process_pool import _get_global_process_pool
from .dashboard_communicator import DashboardCommunicator
from .chromosome_worker import enhanced_search_single_chromosome_worker
from .chunk_worker import _search_genome_chunk_worker_simple
from ..analyzer import GenomeAnalyzer
from ..file_processing import create_shared_gtf_cache, load_annotations_from_shared_cache
from .server_monitor import IS_SERVER, SERVER_MAX_WORKERS, _monitor_server_resources


def enhanced_parallel_chromosome_search(chromosomes: List[str], fasta_path: str, gtf_path: str,
                                       pattern: str, max_mismatches: int, boundary_bp: int, 
                                       command_queue: Optional[multiprocessing.Queue] = None,
                                       pre_loaded_genome: str = None, 
                                       pre_loaded_annotations: Dict = None,
                                       pre_loaded_boundaries: Dict = None) -> List:
    """Enhanced parallel search with Linux compatibility and server resource management."""
    
    # Initialize dashboard communicator
    dashboard = DashboardCommunicator(command_queue)
    
    print(f"[INFO] Starting parallel search - {len(chromosomes)} chromosomes")
    print(f"[INFO] Monitor progress in dashboard window")
    
    # Server resource monitoring
    if IS_SERVER:
        print(f"[INFO] Server environment detected - resource limits: {SERVER_MAX_WORKERS} workers, {SERVER_MEMORY_LIMIT_GB}GB memory")
        _monitor_server_resources()
    
    # CREATE SHARED GTF CACHE (eliminates worker reloading overhead)
    print(f"[INFO] Creating shared GTF cache...")
    cache_file = create_shared_gtf_cache(gtf_path)
    print(f"[INFO] Cache ready - workers will load annotations instantly")
    
    # Determine optimal worker count for server environment
    max_workers = min(SERVER_MAX_WORKERS if IS_SERVER else os.cpu_count(), len(chromosomes))
    
    # Send initial dashboard info
    algorithm = "Aho-Corasick" if max_mismatches == 0 else "Myers Bit-Vector"
    dashboard.send_search_info(pattern, chromosomes, algorithm)
    
    # Send sequence file info if available
    if hasattr(enhanced_parallel_chromosome_search, 'sequence_file_path'):
        dashboard.command_queue.put(("SEQUENCE_UPDATE", {
            'filename': enhanced_parallel_chromosome_search.sequence_file_path,
            'total_sequences': len(chromosomes),
            'sequence_number': 0
        }))
    
    # Add sequence tracking for pattern progress
    # NOTE: Don't override the total_sequences count - let the UI handle this
    if command_queue:
        command_queue.put(("SEQUENCE_UPDATE", {
            'current_sequence': f"Pattern: {pattern[:30]}...",
            'sequence_number': 1
            # Don't set total_sequences here - it should come from the UI
        }))
    
    all_results = []
    completed_count = 0
    
    # Get global process pool
    executor = _get_global_process_pool(max_workers)
    
    # Submit tasks
    future_to_chrom = {}
    worker_to_chrom = {}  # Track which worker processes which chromosome
    
    # OPTIMIZATION: Pre-load genome sequence, extract boundaries, AND pre-load ALL annotations
    print(f"[INFO] Pre-loading genome sequence, extracting boundaries, and pre-loading ALL annotations...")
    
    # Update dashboard with initialization stage
    if command_queue:
        command_queue.put(("INIT_STAGE", {
            'stage': 'Genome Loading',
            'details': 'Loading complete genome sequence',
            'progress_percentage': 10
        }))
    
    try:
        # OPTIMIZATION: Use pre-loaded data if available to avoid reloading
        if pre_loaded_genome and pre_loaded_annotations and pre_loaded_boundaries:
            print(f"[INFO] Using pre-loaded genome data (cache hit!)")
            complete_genome = pre_loaded_genome
            chromosome_boundaries = pre_loaded_boundaries
            print(f"[INFO] Pre-loaded genome: {len(complete_genome):,} bp, {len(chromosome_boundaries)} chromosomes")
            
            # Update dashboard with cache hit
            if command_queue:
                command_queue.put(("INIT_STAGE", {
                    'stage': 'Cache Hit - Using Pre-loaded Data',
                    'details': f'{len(complete_genome):,} bp, {len(chromosome_boundaries)} chromosomes',
                    'progress_percentage': 80
                }))
        else:
            # Load the complete genome once (this is the bottleneck)
            print(f"[INFO] Loading complete genome sequence ({len(fasta_path)} bytes)...")
            from ..file_processing import load_genome_for_chromosome
            complete_genome = load_genome_for_chromosome(fasta_path, "ALL", use_memory_mapping=True)
            print(f"[INFO] Complete genome loaded successfully ({len(complete_genome):,} bp)")
            
            # Update dashboard with genome loaded stage
            if command_queue:
                command_queue.put(("INIT_STAGE", {
                    'stage': 'Genome Loaded',
                    'details': f'{len(complete_genome):,} bp loaded',
                    'progress_percentage': 30
                }))
            
            # Extract chromosome boundaries from GTF annotations for accurate sequence extraction
            print(f"[INFO] Extracting chromosome boundaries from GTF annotations...")
            
            # Update dashboard with boundary extraction stage
            if command_queue:
                command_queue.put(("INIT_STAGE", {
                    'stage': 'Extracting Boundaries',
                    'details': 'Parsing GTF for chromosome positions',
                    'progress_percentage': 40
                }))
            
            chromosome_boundaries = {}
            
            # OPTIMIZATION: Use more efficient GTF parsing with early termination
            print(f"[INFO] Loading complete genome sequence from {os.path.basename(fasta_path)}...")
            with open(gtf_path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 9:
                        continue
                    seqid, source, feature_type, start_str, end_str = parts[:5]
                    
                    if feature_type in ('gene', 'CDS') and seqid in chromosomes:
                        try:
                            start = int(start_str)
                            end = int(end_str)
                            
                            if seqid not in chromosome_boundaries:
                                chromosome_boundaries[seqid] = {'start': start, 'end': end}
                            else:
                                chromosome_boundaries[seqid]['start'] = min(chromosome_boundaries[seqid]['start'], start)
                                chromosome_boundaries[seqid]['end'] = max(chromosome_boundaries[seqid]['end'], end)
                        except ValueError:
                            continue
            
            print(f"[INFO] Extracted boundaries for {len(chromosome_boundaries)} chromosomes")
            
            # Update dashboard with boundaries extracted stage
            if command_queue:
                command_queue.put(("INIT_STAGE", {
                    'stage': 'Boundaries Ready',
                    'details': f'{len(chromosome_boundaries)} chromosomes mapped',
                    'progress_percentage': 50
                }))
            
            # CRITICAL OPTIMIZATION: Pre-load ALL chromosome annotations into shared memory
            print(f"[INFO] Pre-loading ALL chromosome annotations for instant worker access...")
            
            # Update dashboard with annotation pre-loading stage
            if command_queue:
                command_queue.put(("INIT_STAGE", {
                    'stage': 'Pre-loading Annotations',
                    'details': 'Loading all chromosome annotations',
                    'progress_percentage': 60
                }))
            
            pre_loaded_annotations = {}
            
            for chrom in chromosomes:
                try:
                    # Load annotations for each chromosome from shared cache
                    genes, cds = load_annotations_from_shared_cache(cache_file, chrom)
                    pre_loaded_annotations[chrom] = {'genes': genes, 'cds': cds}
                    print(f"[INFO] Pre-loaded annotations for {chrom}: {len(genes)} genes, {len(cds)} CDS")
                    
                    # Update dashboard with individual chromosome annotation progress
                    if command_queue:
                        annotation_percentage = 60 + (len(pre_loaded_annotations) / len(chromosomes)) * 20
                        command_queue.put(("ANNOTATION_PROGRESS", {
                            'chromosome': chrom,
                            'genes': len(genes),
                            'cds': len(cds),
                            'total_loaded': len(pre_loaded_annotations),
                            'progress_percentage': int(annotation_percentage)
                        }))
                        
                except Exception as e:
                    print(f"[WARNING] Failed to pre-load annotations for {chrom}: {e}")
                    pre_loaded_annotations[chrom] = {'genes': [], 'cds': []}
            
            print(f"[INFO] Pre-loaded annotations for {len(pre_loaded_annotations)} chromosomes")
            
            # Update dashboard with annotation pre-loading complete stage
            if command_queue:
                command_queue.put(("INIT_STAGE", {
                    'stage': 'Annotations Ready',
                    'details': f'All {len(pre_loaded_annotations)} chromosomes ready',
                    'progress_percentage': 80
                }))
        
    except Exception as e:
        print(f"[ERROR] Failed to pre-load genome: {e}")
        # Fallback to individual loading
        complete_genome = None
        chromosome_boundaries = {}
        pre_loaded_annotations = {}
        
        # Update dashboard with error stage
        if command_queue:
            command_queue.put(("INIT_STAGE", {
                'stage': 'Initialization Failed',
                'details': f'Error: {str(e)[:50]}',
                'progress_percentage': 0
            }))
    
    # Update dashboard with worker distribution stage
    if command_queue:
        command_queue.put(("INIT_STAGE", {
            'stage': 'Distributing Workers',
            'details': f'Assigning {len(chromosomes)} chromosomes to {max_workers} workers'
        }))
    
    for i, chrom in enumerate(chromosomes):
        future = executor.submit(enhanced_search_single_chromosome_worker, 
                               chrom, fasta_path, gtf_path, 
                               pattern, max_mismatches, boundary_bp, i, complete_genome, chromosome_boundaries, pre_loaded_annotations)
        future_to_chrom[future] = chrom
        
        # Update dashboard with worker status
        worker_id = i % max_workers
        worker_to_chrom[chrom] = worker_id  # Store worker assignment
        
        # Calculate assignment percentage
        assignment_percentage = ((i + 1) / len(chromosomes)) * 100
        
        dashboard.update_worker(worker_id, "Starting", int(assignment_percentage), f"Loading {chrom}")
        
        # Send detailed worker initialization info with percentage
        if command_queue:
            command_queue.put(("WORKER_INIT", {
                'worker_id': worker_id,
                'chromosome': chrom,
                'status': 'Initializing',
                'stage': 'Sequence Extraction',
                'progress_percentage': int(assignment_percentage)
            }))
        
        # Update dashboard with worker assignment progress
        if command_queue:
            command_queue.put(("WORKER_DISTRIBUTION", {
                'assigned': i + 1,
                'total': len(chromosomes),
                'current_chromosome': chrom,
                'worker_id': worker_id,
                'assignment_percentage': int(assignment_percentage)
            }))
    
    # Update dashboard with processing stage
    if command_queue:
        command_queue.put(("INIT_STAGE", {
            'stage': 'Processing Results',
            'details': 'Monitoring worker completion and collecting results',
            'progress_percentage': 85
        }))
    
    # Start worker progress simulation for dashboard transparency
    def simulate_worker_progress():
        """Simulate worker progress updates for dashboard transparency."""
        worker_progress = {i: 0 for i in range(max_workers)}
        start_time = time.time()
        
        while completed_count < len(chromosomes):
            elapsed = time.time() - start_time
            estimated_total_time = elapsed * len(chromosomes) / max(completed_count, 1)
            
            for worker_id in range(max_workers):
                if worker_id in worker_progress:
                    # Simulate progress based on time elapsed
                    progress = min(95, int((elapsed / estimated_total_time) * 100))
                    if progress > worker_progress[worker_id]:
                        worker_progress[worker_id] = progress
                        
                        # Update dashboard with simulated worker progress
                        if command_queue:
                            try:
                                command_queue.put(("WORKER_PROGRESS_SIMULATION", {
                                    'worker_id': worker_id,
                                    'progress_percentage': progress,
                                    'stage': 'Processing',
                                    'details': f'Estimated {progress}% complete'
                                }), timeout=0.1)
                            except:
                                pass
            
            time.sleep(2.0)  # Update every 2 seconds to reduce CPU usage
    
    # Start progress simulation in background thread
    progress_thread = threading.Thread(target=simulate_worker_progress, daemon=True)
    progress_thread.start()
    
    # Process results as they complete with enhanced progress tracking
    for future in as_completed(future_to_chrom):
        chrom = future_to_chrom[future]
        completed_count += 1
        
        # Calculate completion percentage
        completion_percentage = (completed_count / len(chromosomes)) * 100
        
        # Update dashboard with individual chromosome completion
        if command_queue:
            command_queue.put(("CHROMOSOME_COMPLETE", {
                'chromosome': chrom,
                'completed': completed_count,
                'total': len(chromosomes),
                'percentage': int(completion_percentage)
            }))
        
        try:
            success, chrom_name, results = future.result()
            
            if success and results:
                # Filter valid results and send to dashboard
                valid_results = []
                for result in results:
                    if isinstance(result, (tuple, list)) and len(result) >= 8:
                        # Add chromosome to result if not present
                        if len(result) == 8:
                            full_result = (chrom_name,) + tuple(result)
                        else:
                            full_result = tuple(result)
                            
                        # Send to dashboard
                        if command_queue:
                            command_queue.put(("NEW_RESULT", full_result))
                        
                        valid_results.append(full_result)
                        all_results.append(full_result)
                
                # Update dashboard with completion status
                worker_id = worker_to_chrom.get(chrom_name, 0)  # Get the correct worker ID
                dashboard.update_worker(worker_id, "Completed", 100, 
                                      f"{chrom_name}: {len(valid_results)} matches")
                
                # Update dashboard with match summary
                if command_queue:
                    command_queue.put(("MATCH_SUMMARY", {
                        'chromosome': chrom_name,
                        'matches': len(valid_results),
                        'total_matches': len(all_results)
                    }))
                
            else:
                worker_id = worker_to_chrom.get(chrom_name, 0)  # Get the correct worker ID
                dashboard.update_worker(worker_id, "Completed", 100, f"{chrom_name}: No matches")
                
                # Update dashboard with no matches
                if command_queue:
                    command_queue.put(("MATCH_SUMMARY", {
                        'chromosome': chrom_name,
                        'matches': 0,
                        'total_matches': len(all_results)
                    }))
            
            # Update sequence progress
            if command_queue:
                command_queue.put(("SEQUENCE_UPDATE", {
                    'current_sequence': f"Pattern: {pattern[:30]}...",
                    'sequence_number': 1,
                    'total_sequences': 1
                }))
                
            # Update overall progress
            dashboard.update_overall_progress(completed_count, len(chromosomes), len(all_results))
            
        except Exception as e:
            print(f"[ERROR] {chrom} failed: {str(e)}")
            worker_id = worker_to_chrom.get(chrom, 0)  # Get the correct worker ID
            dashboard.update_worker(worker_id, "Failed", 0, f"{chrom}: Error")
            
            # Update dashboard with error
            if command_queue:
                command_queue.put(("CHROMOSOME_ERROR", {
                    'chromosome': chrom,
                    'error': str(e)[:100],
                    'completed': completed_count,
                    'total': len(chromosomes)
                }))
    
    print(f"\n[INFO] Search complete: {len(all_results)} matches across {len(chromosomes)} chromosomes")
    print(f"[INFO] Results displayed in dashboard window")
    
    return all_results


def parallel_single_genome_search(analyzer: GenomeAnalyzer, pattern: str, max_mismatches: int, boundary_bp: int, command_queue=None) -> List:
    """
    Perform parallel search on a single genome using multiple CPU cores.
    Simplified for Homo Sapiens only analysis.
    """
    print(f"[INFO] Starting parallel single genome search")
    print(f"[INFO] Detected analyzer type: Homo Sapiens")
    print(f"[DEBUG] Function parameters: pattern='{pattern}', max_mismatches={max_mismatches}, boundary_bp={boundary_bp}")
    print(f"[DEBUG] Analyzer type: {type(analyzer)}")
    print(f"[DEBUG] Analyzer has sequence: {hasattr(analyzer, 'sequence')}")
    if hasattr(analyzer, 'sequence'):
        print(f"[DEBUG] Sequence length: {len(analyzer.sequence):,} bp")
    
    # Send search info to dashboard
    if command_queue:
        print(f"[DEBUG] Sending SEARCH_INFO to dashboard")
        command_queue.put(("SEARCH_INFO", {
            'pattern': pattern,
            'chromosomes': ['single_genome'],
            'algorithm': 'Parallel Single Genome Search'
        }))
        print(f"[DEBUG] SEARCH_INFO sent successfully")
    else:
        print(f"[DEBUG] No command_queue provided")
    
    # Determine optimal worker count
    max_workers = multiprocessing.cpu_count()
    print(f"[INFO] Using {max_workers} workers for parallel processing")
    print(f"[DEBUG] CPU count: {max_workers}")
    
    # Get global process pool
    if command_queue:
        print(f"[DEBUG] Sending UPDATE_OVERALL - Initializing process pool...")
        command_queue.put(("UPDATE_OVERALL", {
            'completed': 0,
            'total': 4,
            'percentage': 0.0,
            'status': 'Initializing process pool...'
        }))
        print(f"[DEBUG] UPDATE_OVERALL sent successfully")
    
    print(f"[DEBUG] About to call _get_global_process_pool({max_workers})")
    print(f"[DEBUG] This is where the function might hang...")
    
    executor = _get_global_process_pool(max_workers)
    
    print(f"[DEBUG] _get_global_process_pool returned: {type(executor)}")
    print(f"[DEBUG] Executor has submit method: {hasattr(executor, 'submit')}")
    
    if command_queue:
        command_queue.put(("UPDATE_OVERALL", {
            'completed': 1,
            'total': 4,
            'percentage': 25.0,
            'status': 'Process pool ready, preparing chunks...'
        }))
    
    print(f"[DEBUG] Process pool ready, preparing chunks...")
    
    # Split the genome into chunks for parallel processing
    sequence = analyzer.sequence
    print(f"[DEBUG] Sequence type: {type(sequence)}")
    print(f"[DEBUG] Sequence length: {len(sequence):,} bp")
    
    chunk_size = len(sequence) // max_workers
    print(f"[DEBUG] Chunk size: {chunk_size:,} bp")
    
    chunks = []
    print(f"[DEBUG] Creating {max_workers} chunks...")
    
    for i in range(max_workers):
        start = i * chunk_size
        end = start + chunk_size if i < max_workers - 1 else len(sequence)
        chunks.append((start, end))
        print(f"[DEBUG] Chunk {i}: {start:,} - {end:,} ({end-start:,} bp)")
        
        # Update worker status for each chunk
        if command_queue:
            print(f"[DEBUG] Sending UPDATE_WORKER for worker {i}")
            command_queue.put(("UPDATE_WORKER", {
                'worker_id': i,
                'status': 'Ready',
                'progress': 0,
                'chromosome': 'single_genome',
                'speed': '0 bp/s',
                'memory': '0MB'
            }))
            print(f"[DEBUG] UPDATE_WORKER sent for worker {i}")
    
    print(f"[DEBUG] Created {len(chunks)} chunks successfully")
    
    # Print header for live output (only if verbose mode)
    if hasattr(analyzer, 'verbose_output') and analyzer.verbose_output:
        from ..ui import _print_results_header
        _print_results_header(include_chrom=False)
    
    if command_queue:
        command_queue.put(("UPDATE_OVERALL", {
            'completed': 2,
            'total': 4,
            'percentage': 50.0,
            'status': 'Submitting search tasks...'
        }))
    
    print(f"[DEBUG] Submitting chunk search tasks...")
    print(f"[DEBUG] Executor type: {type(executor)}")
    print(f"[DEBUG] Number of chunks to submit: {len(chunks)}")
    
    # Submit chunk search tasks
    future_to_chunk = {}
    print(f"[DEBUG] Submitting tasks one by one...")
    
    for i, (start, end) in enumerate(chunks):
        print(f"[DEBUG] Submitting chunk {i}: {start:,} - {end:,}")
        try:
            future = executor.submit(_search_genome_chunk_worker_simple, analyzer, start, end, pattern, max_mismatches, boundary_bp)
            future_to_chunk[future] = i
            print(f"[DEBUG] Chunk {i} submitted successfully, future: {type(future)}")
        except Exception as e:
            print(f"[ERROR] Failed to submit chunk {i}: {e}")
            print(f"[DEBUG] Exception type: {type(e).__name__}")
            print(f"[DEBUG] Exception details: {str(e)}")
    
    print(f"[DEBUG] Submitted {len(future_to_chunk)} tasks successfully")
    
    if command_queue:
        print(f"[DEBUG] Sending UPDATE_OVERALL - Search tasks submitted...")
        command_queue.put(("UPDATE_OVERALL", {
            'completed': 3,
            'total': 4,
            'percentage': 75.0,
            'status': 'Search tasks submitted, processing results...'
        }))
        print(f"[DEBUG] UPDATE_OVERALL sent successfully")
    
    all_results = []
    print(f"[DEBUG] Starting to collect results...")
    
    # Collect results as they complete
    completed = 0
    print(f"[DEBUG] Entering result collection loop with {len(future_to_chunk)} futures")
    
    for future in as_completed(future_to_chunk):
        chunk_id = future_to_chunk[future]
        completed += 1
        
        # Update worker status
        if command_queue:
            command_queue.put(("UPDATE_WORKER", {
                'worker_id': chunk_id,
                'status': 'Processing',
                'progress': 100,  # Chunk completed
                'chromosome': 'single_genome',
                'speed': f'{chunk_id + 1} chunks/sec',
                'memory': '50MB'
            }))
        
        try:
            results = future.result()
            print(f"[INFO] Completed chunk {chunk_id + 1} ({completed}/{len(chunks)})")
            
            # Update overall progress
            if command_queue:
                percentage = (completed / len(chunks)) * 100
                command_queue.put(("UPDATE_OVERALL", {
                    'completed': completed,
                    'total': len(chunks),
                    'percentage': percentage
                }))
            
            # Print results live
            for result in results:
                all_results.append(result)
                
                # Send result to dashboard
                if command_queue and result:
                    command_queue.put(("NEW_RESULT", {
                        'start': result[0],
                        'end': result[1], 
                        'strand': result[2],
                        'pattern': result[3],
                        'sequence': result[4],
                        'chromosome': 'single_genome',
                        'gene': result[7] if len(result) > 7 else 'Unknown'
                    }))
                
                from ..ui import _print_row_tuple
                _print_row_tuple(result, include_chrom=False)
                
        except Exception as e:
            print(f"[INFO] Exception processing chunk {chunk_id}: {e}")
    
    if command_queue:
        command_queue.put(("UPDATE_OVERALL", {
            'completed': 4,
            'total': 4,
            'percentage': 100.0,
            'status': f'Search complete! Found {len(all_results)} matches'
        }))
    
    print(f"[INFO] Parallel genome search complete. Found {len(all_results)} total matches")
    return all_results


# Export functions
__all__ = [
    'enhanced_parallel_chromosome_search',
    'parallel_single_genome_search'
]
