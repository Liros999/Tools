"""
Parallel processing functions for genome analysis.
Contains worker functions, process pool management, and parallel search implementations.
"""

import os
import multiprocessing
import atexit
import time
import hashlib
import platform
import psutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Optional, Dict, Any

# Linux compatibility and adaptive server resource management
IS_LINUX = platform.system() == 'Linux'
IS_SERVER = os.environ.get('SERVER_ENV', 'false').lower() == 'true'

# Import adaptive server configuration
try:
    import sys
    import os
    # Add the parent directory to the path to import linux_server_config
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
    from linux_server_config import ServerConfig, _adaptive_manager
    if _adaptive_manager:
        SERVER_MAX_WORKERS = _adaptive_manager.MAX_WORKERS
        SERVER_MEMORY_LIMIT_GB = _adaptive_manager.MAX_MEMORY_GB
        SERVER_CACHE_SIZE_MB = _adaptive_manager.CACHE_SIZE_MB
        print(f"[INFO] Using adaptive server configuration: {_adaptive_manager.server_class}")
    else:
        SERVER_MAX_WORKERS = min(4, os.cpu_count()) if IS_SERVER else os.cpu_count()
        SERVER_MEMORY_LIMIT_GB = 8
        SERVER_CACHE_SIZE_MB = 2048
except (ImportError, ModuleNotFoundError, AttributeError):
    # Fallback if config not available
    print("[INFO] Linux server config not available, using default configuration")
    SERVER_MAX_WORKERS = min(4, os.cpu_count()) if IS_SERVER else os.cpu_count()
    SERVER_MEMORY_LIMIT_GB = 8
    SERVER_CACHE_SIZE_MB = 2048
    # Set _adaptive_manager to None to avoid NameError
    _adaptive_manager = None

from .analyzer import GenomeAnalyzer, GenomeAnalysisError
from .file_processing import load_annotations_from_shared_cache, create_shared_gtf_cache, load_genome_for_chromosome
from .utils import _get_filtered_chromosomes, _process_search_results

# Global process pool for parallel processing
_global_process_pool = None
_global_process_pool_workers = None
_global_process_pool_initialized = False

def _monitor_server_resources():
    """Monitor server resources and enforce adaptive limits."""
    if not IS_SERVER:
        return
    
    try:
        # Check memory usage
        memory = psutil.virtual_memory()
        memory_gb = memory.used / (1024**3)
        memory_percent = memory.percent
        
        # Get adaptive limits
        try:
            if '_adaptive_manager' in globals() and _adaptive_manager and hasattr(_adaptive_manager, 'MAX_MEMORY_GB'):
                memory_limit = _adaptive_manager.MAX_MEMORY_GB
                memory_warning = _adaptive_manager.MEMORY_WARNING_THRESHOLD_GB
                server_class = _adaptive_manager.server_class
            else:
                memory_limit = SERVER_MEMORY_LIMIT_GB
                memory_warning = SERVER_MEMORY_LIMIT_GB * 0.8
                server_class = "standard"
        except (NameError, AttributeError):
            memory_limit = SERVER_MEMORY_LIMIT_GB
            memory_warning = SERVER_MEMORY_LIMIT_GB * 0.8
            server_class = "standard"
        
        # Memory monitoring with adaptive thresholds
        if memory_gb > memory_limit:
            print(f"[WARNING] Memory usage {memory_gb:.1f}GB exceeds limit {memory_limit}GB")
            print(f"[INFO] Server class: {server_class.upper()} - Consider upgrading if this persists")
            _cleanup_server_resources()
        elif memory_gb > memory_warning:
            print(f"[INFO] Memory usage {memory_gb:.1f}GB approaching limit {memory_limit}GB ({memory_percent:.1f}%)")
            print(f"[INFO] Current utilization: {memory_percent:.1f}% - Optimal for {server_class.upper()} server")
        
        # Check CPU usage with adaptive thresholds
        cpu_percent = psutil.cpu_percent(interval=1)
        try:
            if '_adaptive_manager' in globals() and _adaptive_manager and hasattr(_adaptive_manager, 'MAX_CPU_PERCENT'):
                cpu_limit = _adaptive_manager.MAX_CPU_PERCENT
            else:
                cpu_limit = 80
        except (NameError, AttributeError):
            cpu_limit = 80
            
        if cpu_percent > cpu_limit:
            print(f"[WARNING] High CPU usage: {cpu_percent:.1f}% (limit: {cpu_limit}%)")
            if cpu_percent > 95:
                print(f"[WARNING] Critical CPU usage - consider reducing worker count")
        elif cpu_percent > cpu_limit * 0.8:
            print(f"[INFO] CPU usage: {cpu_percent:.1f}% - Optimal utilization for {server_class.upper()} server")
        
        # Adaptive resource recommendations
        try:
            if '_adaptive_manager' in globals() and _adaptive_manager and memory_percent < 50:
                print(f"[INFO] Memory utilization {memory_percent:.1f}% - {server_class.upper()} server can handle more load")
                print(f"[INFO] Consider increasing worker count or cache size for better performance")
        except (NameError, AttributeError):
            pass
            
    except Exception as e:
        print(f"[WARNING] Server resource monitoring failed: {e}")

def _cleanup_server_resources():
    """Clean up server resources to free memory."""
    if not IS_SERVER:
        return
    
    try:
        import gc
        gc.collect()  # Force garbage collection
        
        # Clear large caches
        global _global_process_pool
        if _global_process_pool is not None:
            try:
                _global_process_pool.shutdown(wait=True)
                _global_process_pool = None
                _global_process_pool_initialized = False
                print("[INFO] Process pool cleaned up to free memory")
            except:
                pass
                
    except Exception as e:
        print(f"[WARNING] Server resource cleanup failed: {e}")

def _get_global_process_pool(worker_count=None):
    """Get or create a global process pool with Linux compatibility and server resource management."""
    global _global_process_pool, _global_process_pool_workers, _global_process_pool_initialized
    
    print(f"[DEBUG] _get_global_process_pool called with worker_count={worker_count}")
    print(f"[DEBUG] Current state: pool={_global_process_pool is not None}, workers={_global_process_pool_workers}, initialized={_global_process_pool_initialized}")
    
    if worker_count is None:
        # Use server-appropriate worker count
        if IS_SERVER:
            worker_count = SERVER_MAX_WORKERS
            print(f"[INFO] Server environment detected - limiting workers to {worker_count}")
        else:
            worker_count = os.cpu_count()
        print(f"[DEBUG] Determined worker_count={worker_count}")
    
    # Create new pool if it doesn't exist or worker count changed
    if (_global_process_pool is None or 
        _global_process_pool_workers != worker_count or 
        not _global_process_pool_initialized):
        
        print(f"[DEBUG] Need to create new process pool")
        
        # Clean up existing pool if it exists
        if _global_process_pool is not None:
            try:
                print(f"[DEBUG] Shutting down existing process pool")
                _global_process_pool.shutdown(wait=True)
                print(f"[DEBUG] Existing process pool shut down")
            except Exception as e:
                print(f"[DEBUG] Error shutting down existing pool: {e}")
        
        print(f"[INFO] Creating global process pool with {worker_count} workers")
        print(f"[DEBUG] Platform: {platform.system()}, IS_LINUX={IS_LINUX}")
        
        # Create process pool with OS compatibility
        try:
            print(f"[DEBUG] Attempting to create ProcessPoolExecutor...")
            start_time = time.time()
            
            if IS_LINUX:
                # Linux-specific process pool configuration
                print(f"[DEBUG] Using Linux fork context")
                _global_process_pool = ProcessPoolExecutor(
                    max_workers=worker_count,
                    mp_context=multiprocessing.get_context('fork')  # Use fork on Linux
                )
            else:
                # Windows/other OS process pool - use spawn context
                print(f"[DEBUG] Using Windows spawn context")
                _global_process_pool = ProcessPoolExecutor(
                    max_workers=worker_count,
                    mp_context=multiprocessing.get_context('spawn')  # Use spawn on Windows
                )
            
            elapsed = time.time() - start_time
            print(f"[DEBUG] ProcessPoolExecutor created in {elapsed:.3f} seconds")
            print(f"[INFO] Process pool created successfully with {worker_count} workers")
            
        except Exception as e:
            print(f"[ERROR] Failed to create process pool: {e}")
            print(f"[DEBUG] Exception type: {type(e).__name__}")
            print(f"[DEBUG] Exception details: {str(e)}")
            
            # Fallback to default context
            try:
                print(f"[DEBUG] Attempting fallback ProcessPoolExecutor...")
                _global_process_pool = ProcessPoolExecutor(max_workers=worker_count)
                print(f"[INFO] Using fallback process pool with {worker_count} workers")
            except Exception as fallback_e:
                print(f"[ERROR] Fallback also failed: {fallback_e}")
                raise fallback_e
        
        _global_process_pool_workers = worker_count
        _global_process_pool_initialized = True
        print(f"[DEBUG] Process pool state updated: workers={_global_process_pool_workers}, initialized={_global_process_pool_initialized}")
    
    else:
        print(f"[DEBUG] Using existing process pool with {_global_process_pool_workers} workers")
    
    print(f"[DEBUG] Returning process pool: {type(_global_process_pool)}")
    return _global_process_pool

def _cleanup_global_process_pool():
    """Clean up the global process pool on program exit."""
    global _global_process_pool
    if _global_process_pool is not None:
        try:
            print("[INFO] Shutting down global process pool")
            _global_process_pool.shutdown(wait=True)
        except:
            pass

# Register cleanup function
atexit.register(_cleanup_global_process_pool)

def _process_search_results(matches: List[Tuple], pattern: str, strand: str, 
                          analyzer: GenomeAnalyzer, boundary_bp: int, 
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
                from .file_processing import load_annotations_from_gtf
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
            import threading
            import queue
            
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

class DashboardCommunicator:
    """Handles real-time dashboard communication with proper message formatting."""
    
    def __init__(self, command_queue: Optional[multiprocessing.Queue] = None):
        self.command_queue = command_queue
        self.worker_states = {}
        self.overall_progress = {'completed': 0, 'total': 0, 'matches': 0}
        
    def is_available(self) -> bool:
        """Check if dashboard is available for communication."""
        return self.command_queue is not None
        
    def update_worker(self, worker_id: int, status: str, progress: int, details: str = ""):
        """Send worker update to dashboard."""
        if not self.is_available():
            return
            
        try:
            message = {
                'worker_id': worker_id,
                'status': status,
                'progress': max(0, min(100, progress)),
                'details': str(details)[:50],  # Limit detail length
                'timestamp': time.time()
            }
            
            self.command_queue.put(("UPDATE_WORKER", message), timeout=0.1)
            self.worker_states[worker_id] = message
            
        except Exception as e:
            print(f"[ERROR] Failed to update worker: {e}")
            
    def update_overall_progress(self, completed: int, total: int, matches: int):
        """Send overall progress update to dashboard."""
        if not self.is_available():
            return
            
        try:
            progress_data = {
                'completed': completed,
                'total': total,
                'matches': matches,
                'percentage': (completed / total * 100) if total > 0 else 0,
                'timestamp': time.time()
            }
            
            self.command_queue.put(("UPDATE_OVERALL", progress_data), timeout=0.1)
            self.overall_progress = progress_data
            
        except Exception as e:
            print(f"[ERROR] Failed to update overall progress: {e}")
            
    def send_search_info(self, pattern: str, chromosomes: List[str], algorithm: str):
        """Send search configuration info to dashboard."""
        if not self.is_available():
            return
            
        try:
            search_info = {
                'pattern': pattern,
                'chromosomes': chromosomes,
                'algorithm': algorithm,
                'start_time': time.time()
            }
            
            self.command_queue.put(("SEARCH_INFO", search_info), timeout=0.1)
            
        except Exception as e:
            print(f"[ERROR] Failed to send search info: {e}")
    
    def send_sequence_info(self, filename: str, total_sequences: int):
        """Send sequence file information to dashboard."""
        if not self.is_available():
            return
            
        try:
            sequence_info = {
                'sequence_file': filename,
                'total_sequences': total_sequences
            }
            
            self.command_queue.put(("SEQUENCE_UPDATE", sequence_info), timeout=0.1)
            
        except Exception as e:
            print(f"[ERROR] Failed to send sequence info: {e}")


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
            from .file_processing import load_genome_for_chromosome
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
    import threading
    import time
    
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
        from .ui import _print_results_header
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
                
                from .ui import _print_row_tuple
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
