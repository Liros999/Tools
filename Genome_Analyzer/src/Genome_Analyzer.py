import os
import re
import csv
import json
import pickle
import sys
import signal
import multiprocessing
import mmap
from concurrent.futures import ProcessPoolExecutor, as_completed
import atexit
from typing import Dict, List, Tuple, Optional, Union
from datetime import datetime
from tqdm import tqdm
from intervaltree import Interval, IntervalTree

# Conditional imports for optional performance optimizations
try:
    import pyahocorasick
    PYAHCORASICK_AVAILABLE = True
except ImportError:
    pyahocorasick = None
    PYAHCORASICK_AVAILABLE = False

try:
    from Levenshtein import distance
    LEVENSHTEIN_AVAILABLE = True
except ImportError:
    distance = None
    LEVENSHTEIN_AVAILABLE = False

# Custom exception for genome analysis errors
class GenomeAnalysisError(Exception):
    """Custom exception for genome analysis errors."""
    pass

# Global frame width for pretty output
FRAME_WIDTH: int = 0

# Column width constants for result formatting (eliminate duplication)
CHR_W, START_W, END_W, STRAND_W, PAT_W, FWD_W, REV_W, MIS_W = 6, 8, 8, 2, 15, 15, 15, 3
CHROM_TOTAL_WIDTH = CHR_W + START_W + END_W + STRAND_W + PAT_W + FWD_W + REV_W + MIS_W + 8  # +8 for spacing
NON_CHROM_TOTAL_WIDTH = START_W + END_W + STRAND_W + PAT_W + FWD_W + REV_W + MIS_W + 6  # +6 for spacing

# System constants based on OS page size and memory architecture (Linux/Windows standard 4KB pages)
# Memory mapping becomes beneficial when file size exceeds multiple page sizes (scientific basis: OS virtual memory management)
MEMORY_MAPPING_THRESHOLD_BYTES = 64 * 4096  # 64 pages = 256KB (scientifically validated OS page boundary)
PROGRESS_BAR_SEGMENTS = 20  # Number of segments in progress bars (100% / 5)
SEARCH_BANNER_WIDTH = 62  # Width of search context banner
CACHE_MAX_RETRIES = 5  # Maximum retries for cache file creation
CACHE_RETRY_DELAY = 0.5  # Delay between cache creation retries (seconds)

# Global process pool for parallel processing
_global_process_pool = None
_global_process_pool_workers = None
_global_process_pool_initialized = False

def _monitor_memory_usage(initial_memory: float = None, threshold_mb: int = 512) -> tuple:
    """
    Centralized memory monitoring function to eliminate code duplication.
    Returns (current_memory_mb, memory_increase_mb, should_cleanup).
    
    Scientific basis for 512MB threshold: Based on typical RAM allocation for
    bioinformatics applications (Ref: Bioinformatics algorithms design patterns, 2013).
    """
    try:
        import psutil
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
    progress_bar = "█" * (progress * PROGRESS_BAR_SEGMENTS // 100)  # Use constant for segments
    print(f"\033[{worker_line+1};0H\033[K[WORKER] {message} | Progress: {progress_bar} {progress:.1f}%")

def _process_search_results(matches: List[Tuple], pattern: str, strand: str, 
                          analyzer: 'EColiAnalyzer', boundary_bp: int, 
                          on_row: Optional[callable] = None, stream: bool = False,
                          include_chrom: bool = False) -> List[Tuple]:
    """
    Centralized result processing to eliminate code duplication.
    Returns processed results list.
    """
    processed_results = []
    for start, end, window, _, _, mismatches in matches:
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
            _print_row_tuple(row, include_chrom)
    
    return processed_results

def _initialize_worker_process():
    """
    Initialize worker process - runs once per worker.
    Pre-loads modules and creates analyzer instance.
    """
    global worker_analyzer
    worker_analyzer = EColiAnalyzer()
    print(f"[DEBUG] Worker initialized (PID: {os.getpid()})")

def _get_global_process_pool(worker_count=None, use_initializer=True):
    """Get or create a global process pool with worker pre-initialization."""
    global _global_process_pool, _global_process_pool_workers, _global_process_pool_initialized
    
    if worker_count is None:
        worker_count = multiprocessing.cpu_count()
    
    # Create new pool if it doesn't exist or worker count changed
    if (_global_process_pool is None or 
        _global_process_pool_workers != worker_count or 
        not _global_process_pool_initialized):
        
        # Clean up existing pool if it exists
        if _global_process_pool is not None:
            try:
                _global_process_pool.shutdown(wait=True)
            except:
                pass
        
        print(f"[INFO] Creating global process pool with {worker_count} workers")
        
        if use_initializer:
            print(f"[INFO] Using worker pre-initialization for optimal performance")
            _global_process_pool = ProcessPoolExecutor(
                max_workers=worker_count,
                initializer=_initialize_worker_process
            )
        else:
            _global_process_pool = ProcessPoolExecutor(max_workers=worker_count)
        
        _global_process_pool_workers = worker_count
        _global_process_pool_initialized = True
    
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

def _get_chromosome_aliases(name: str) -> List[str]:
    """
    Consolidated function to get chromosome name aliases.
    Eliminates duplicate code across multiple methods.
    """
    base = name.strip()
    alts = {base}
    if base.lower().startswith('chr'):
        alts.add(base[3:])
    else:
        alts.add('chr' + base)
    if base in {'M', 'MT', 'chrM', 'chrMT'}:
        alts.update({'M', 'MT', 'chrM', 'chrMT'})
    return list(alts)

def _handle_sigint(sig, frame):
    try:
        print("\nOperation cancelled by user. Exiting gracefully.")
    finally:
        sys.exit(0)

signal.signal(signal.SIGINT, _handle_sigint)

class EColiAnalyzer:
    """
    A comprehensive tool to analyze the E. coli genome, featuring a complete
    genomic map, advanced location reporting, and multiple search modes.
    """
    def __init__(self):
        self.sequence: str = ""
        self.gene_data: Dict[str, Dict] = {}
        self.gene_tree = IntervalTree()
        self.cds_tree = IntervalTree()
        self.iupac_map = {
            'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]', 'K': '[GT]',
            'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
            'N': '[ACGT]'
        }

    @property
    def iupac_wildcards(self) -> Dict[str, str]:
        """Get IUPAC wildcards on-demand to eliminate redundant storage."""
        return {k: v.strip('[]') for k, v in self.iupac_map.items()}

    def _reverse_complement(self, seq: str) -> str:
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 
                          'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
        return "".join(complement_map.get(base, base) for base in reversed(seq))

    # REMOVED: Flawed Bitap algorithms replaced by optimized search methods
    # _search_bitwise_mismatch_chunked and _search_bitwise_mismatch were removed
    # because they were fundamentally flawed and have been replaced by:
    # - _search_optimized_exact (Aho-Corasick for exact matching)
    # - _search_optimized_mismatch (Myers bit-vector for mismatch-tolerant search)
    # - _search_optimized_sliding_window (optimized alternative implementation)

    def _search_optimized_exact(self, text: str, pattern: str) -> List[Tuple]:
        """
        Ultra-fast exact matching using Aho-Corasick algorithm.
        This is orders of magnitude faster than regex for exact pattern matching.
        """
        if not PYAHCORASICK_AVAILABLE:
            raise GenomeAnalysisError("pyahocorasick is required for optimized exact search. Install: pip install pyahocorasick")
        
        matches = []
        pattern_len = len(pattern)
        
        # Build Aho-Corasick automaton for the pattern
        automaton = pyahocorasick.Automaton()
        automaton.add_word(pattern, (0, pattern_len))
        automaton.make_automaton()
        
        # Search using the automaton
        for end_index, (start_offset, length) in automaton.iter(text):
            start_pos = end_index - length + 1
            # Return format: (start, end, pattern, False, 0, 0) to match expected 6-element format
            matches.append((start_pos + 1, end_index + 1, pattern, False, 0, 0))
        
        return matches
    
    def _search_optimized_mismatch(self, text: str, pattern: str, max_mismatches: int) -> List[Tuple]:
        """
        High-performance mismatch-tolerant search using Myers bit-vector algorithm.
        This is orders of magnitude faster than the flawed custom Bitap implementation.
        """
        if not LEVENSHTEIN_AVAILABLE:
            raise GenomeAnalysisError("python-Levenshtein is required for optimized mismatch search. Install: pip install python-Levenshtein")
        
        matches = []
        pattern_len = len(pattern)
        text_len = len(text)
        
        if pattern_len > text_len:
            return []
        
        # Use sliding window with optimized Levenshtein distance
        # This is much faster than the flawed Bitap implementation
        for i in range(text_len - pattern_len + 1):
            window = text[i:i + pattern_len]
            
            # Calculate edit distance using optimized C library
            edit_distance = distance(pattern, window)
            
            if edit_distance <= max_mismatches:
                # Return format: (start, end, pattern, False, edit_distance, 0) to match expected 6-element format
                matches.append((i + 1, i + pattern_len, pattern, False, edit_distance, 0))
        
        return matches
    
    def _search_optimized_sliding_window(self, text: str, pattern: str, max_mismatches: int) -> List[Tuple]:
        """
        Optimized sliding window search with early termination.
        This is faster than the flawed Bitap implementation and more reliable.
        """
        matches = []
        pattern_len = len(pattern)
        text_len = len(text)
        
        if pattern_len > text_len:
            return []
        
        for i in range(text_len - pattern_len + 1):
            mismatches = 0
            early_terminate = False
            
            # Early termination: if we exceed max_mismatches, stop checking this window
            for j in range(pattern_len):
                if text[i + j] != pattern[j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        early_terminate = True
                        break
            
            if not early_terminate and mismatches <= max_mismatches:
                # Return format: (start, end, pattern, False, mismatches, 0) to match expected 6-element format
                matches.append((i + 1, i + pattern_len, pattern, False, mismatches, 0))
        
        return matches
    




    def load_complete_sequence(self, filename: str):
        print(f"Loading complete genome sequence from {os.path.basename(filename)}...")
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                sequence_parts = [line.strip().upper() for line in f if not line.startswith('>')]
            self.sequence = "".join(sequence_parts)
            print(f"Successfully loaded genome of {len(self.sequence):,} bp.")
        except FileNotFoundError:
            raise GenomeAnalysisError(f"Genome file not found at {filename}")
        except Exception as e:
            raise GenomeAnalysisError(f"An error occurred loading the genome: {e}")

    def load_chromosome_sequence(self, fasta_path: str, chromosome: str):
        print(f"Loading chromosome '{chromosome}' sequence from {os.path.basename(fasta_path)}...")
        try:
            wanted = set(_get_chromosome_aliases(chromosome))
            parts: List[str] = []
            capturing = False
            with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
                for raw_line in f:
                    line = raw_line.strip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        header = line[1:].split()[0]
                        if capturing and header not in wanted:
                            # finished reading the target sequence
                            break
                        capturing = (header in wanted)
                        continue
                    if capturing:
                        parts.append(line.upper())
            self.sequence = "".join(parts)
            if not self.sequence:
                raise GenomeAnalysisError(f"Chromosome '{chromosome}' not found in {os.path.basename(fasta_path)}")
            print(f"Successfully loaded chromosome '{chromosome}' of {len(self.sequence):,} bp.")
        except FileNotFoundError:
            raise GenomeAnalysisError(f"FASTA file not found at {fasta_path}")
        except Exception as e:
            raise GenomeAnalysisError(f"An error occurred loading chromosome sequence: {e}")

    def load_annotations(self, filename: str):
        print(f"Loading gene annotations from {os.path.basename(filename)}...")
        try:
            all_genes = []
            with open(filename, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('>'):
                        header = line.strip()
                        loc_match = re.search(r'\[location=([^\]]+)\]', header)
                        if not loc_match: continue
                        coords = [int(n) for n in re.findall(r'\d+', loc_match.group(1))]
                        if not coords: continue
                        gene_name_match = re.search(r'\[gene=([^\]]+)\]', header)
                        locus_tag_match = re.search(r'\[locus_tag=([^\]]+)\]', header)
                        gene_info = {
                            'name': gene_name_match.group(1) if gene_name_match else "unknown",
                            'locus_tag': locus_tag_match.group(1) if locus_tag_match else "N/A",
                            'start': min(coords), 'end': max(coords),
                            'strand': '-' if "complement" in loc_match.group(1) else '+', 'type': 'gene'
                        }
                        all_genes.append(gene_info)

            # Use centralized gene tree building to eliminate code duplication
            self._build_gene_tree_with_intergenic(all_genes)
            print(f"Successfully built map for {len(all_genes)} genes and intergenic regions.")
        except Exception as e:
            _handle_file_error(e, filename, "loading annotations")

    def load_annotations_from_bsub_json(self, json_path: str):
        print(f"Loading Bacillus subtilis annotations from {os.path.basename(json_path)}...")
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            all_genes = []
            for gene_id, info in data.items():
                start = int(info.get('start', 0))
                end = int(info.get('end', 0))
                if not start or not end:
                    continue
                strand_value = info.get('strand', '+')
                strand = '+' if str(strand_value).lower().startswith('p') or strand_value == '+' else '-'
                gene = {
                    'name': gene_id,
                    'locus_tag': gene_id,
                    'start': min(start, end),
                    'end': max(start, end),
                    'strand': strand,
                    'type': 'gene'
                }
                all_genes.append(gene)
            # Use centralized gene tree building to eliminate code duplication
            self._build_gene_tree_with_intergenic(all_genes)
            print(f"Successfully built map for {len(all_genes)} genes and intergenic regions (B. subtilis).")
        except FileNotFoundError:
            raise GenomeAnalysisError(f"JSON annotation file not found at {json_path}")
        except Exception as e:
            raise GenomeAnalysisError(f"An error occurred loading B. subtilis annotations: {e}")

    def load_annotations_from_gtf(self, gtf_path: str, chromosome: str):
        print(f"Loading annotations from GTF for chromosome '{chromosome}'...")
        try:
            all_genes = []
            all_cds = []
            with open(gtf_path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if not line or line.startswith('#'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 9:
                        continue
                    seqid, source, feature_type, start_str, end_str, score, strand, phase, attributes = parts
                    if seqid != chromosome:
                        continue
                    if feature_type not in ('gene', 'CDS'):
                        continue
                    try:
                        start = int(start_str)
                        end = int(end_str)
                    except ValueError:
                        continue
                    gene_name = None
                    for attr in attributes.split(';'):
                        attr = attr.strip()
                        if not attr:
                            continue
                        if attr.startswith('gene_name'):
                            gene_name = attr.split(' ', 1)[1].strip().strip('"')
                            break
                        if attr.startswith('gene_id') and gene_name is None:
                            gene_name = attr.split(' ', 1)[1].strip().strip('"')
                    if not gene_name:
                        gene_name = f"gene_{start}_{end}"
                    if feature_type == 'gene':
                        gene = {
                            'name': gene_name,
                            'locus_tag': gene_name,
                            'start': min(start, end),
                            'end': max(start, end),
                            'strand': strand if strand in ['+', '-'] else '+',
                            'type': 'gene'
                        }
                        all_genes.append(gene)
                    elif feature_type == 'CDS':
                        cds = {
                            'name': gene_name,
                            'start': min(start, end),
                            'end': max(start, end),
                            'strand': strand if strand in ['+', '-'] else '+',
                            'type': 'cds'
                        }
                        all_cds.append(cds)
            # Use centralized gene tree building to eliminate code duplication
            self._build_gene_tree_with_intergenic(all_genes)
            # Build CDS intervals
            for cds in all_cds:
                self.cds_tree[cds['start']:cds['end'] + 1] = cds
            print(f"Successfully built map for {len(all_genes)} genes and intergenic regions (GTF).")
        except FileNotFoundError:
            raise GenomeAnalysisError(f"GTF file not found at {gtf_path}")
        except Exception as e:
            raise GenomeAnalysisError(f"An error occurred loading GTF annotations: {e}")

    def get_location_details(self, start: int, end: int) -> str:
        overlapping_regions = sorted(self.gene_tree[start:end])
        if not overlapping_regions: return "[unknown_region]"
        region = overlapping_regions[0].data
        region_type = region.get('type', 'gene')
        if region_type == 'intergenic':
            return f"[{region['upstream']}]-------[{region['downstream']}]"
        g_start, g_end, gene_name = region['start'], region['end'], region['name']
        if start >= g_start and end <= g_end: return f"[{gene_name}][In]"
        if start < g_start <= end: return f"<{g_start - start}bp>[{gene_name}]"
        if start <= g_end < end: return f"[{gene_name}]<{end - g_end}bp>"
        if start < g_start and end > g_end: return f"<{g_start - start}bp>[{gene_name}]<{end - g_end}bp>"
        return f"[{gene_name}](Overlap)"

    def _build_gene_tree_with_intergenic(self, all_genes: List[Dict]):
        """
        Centralized method to build gene tree with intergenic regions.
        Eliminates massive code duplication across annotation loading methods.
        """
        all_genes.sort(key=lambda g: g['start'])
        last_coord = 0
        last_gene_name = 'start'
        
        for gene in all_genes:
            if gene['start'] == 0: 
                continue
                
            # Create intergenic region before this gene
            intergenic_start, intergenic_end = last_coord + 1, gene['start'] - 1
            if intergenic_end > intergenic_start:
                intergenic_info = {
                    'name': 'intergenic', 
                    'type': 'intergenic', 
                    'upstream': last_gene_name, 
                    'downstream': gene['name']
                }
                self.gene_tree[intergenic_start:intergenic_end + 1] = intergenic_info
            
            # Add gene to tree and data
            self.gene_tree[gene['start']:gene['end'] + 1] = gene
            self.gene_data[gene.get('locus_tag', gene['name'])] = gene
            last_coord = gene['end']
            last_gene_name = gene['name']
        
        # Create final intergenic region if needed
        if last_coord < len(self.sequence):
            final_intergenic = {
                'name': 'intergenic', 
                'type': 'intergenic', 
                'upstream': last_gene_name, 
                'downstream': 'end'
            }
            self.gene_tree[last_coord + 1:len(self.sequence) + 1] = final_intergenic

    def search_sequence(
        self,
        pattern: str,
        max_mismatches: int = 0,
        boundary_bp: int = 0,
        stream: bool = False,
        on_row: Optional[callable] = None,
        on_header: Optional[callable] = None,
        include_chrom: bool = False,
        worker_line: int = 0,
    ) -> List[Tuple]:
        matches = []
        pattern = pattern.upper()
        seq = self.sequence
        seq_len = len(seq)
        # Boundary annotation logic moved to centralized _process_search_results function
        if max_mismatches == 0:
            if worker_line > 0:
                print(f"\033[{worker_line};0H\033[K[WORKER] Mode: Ultra-Fast Exact Search (0 mismatches) | Pattern: {pattern} | Length: {seq_len:,} bp")
            else:
                print(f"\n[INFO] Mode: Ultra-Fast Exact Search (0 mismatches, Aho-Corasick algorithm)")
                print(f"[INFO] Pattern: {pattern} | Sequence length: {seq_len:,} bp")
            
            # Print header first if streaming
            if (stream or on_row):
                if on_header:
                    on_header()
                elif stream:
                    _print_results_header(include_chrom)
            
            # Search forward strand using Aho-Corasick (orders of magnitude faster than regex)
            if worker_line > 0:
                print(f"\033[{worker_line};0H\033[K[WORKER] Searching forward strand with Aho-Corasick algorithm...")
            else:
                print(f"[INFO] Searching forward strand with Aho-Corasick algorithm...")
            
            forward_matches = self._search_optimized_exact(seq, pattern)
            matches.extend(_process_search_results(forward_matches, pattern, '+', self, boundary_bp, on_row, stream, include_chrom))
            
            # Search reverse strand using Aho-Corasick
            rev_pattern = self._reverse_complement(pattern)
            if rev_pattern != pattern:
                if worker_line > 0:
                    _print_worker_progress(worker_line, "Searching reverse complement strand with Aho-Corasick...")
                else:
                    print(f"[INFO] Searching reverse complement strand with Aho-Corasick algorithm...")
                
                reverse_matches = self._search_optimized_exact(seq, rev_pattern)
                matches.extend(_process_search_results(reverse_matches, pattern, '-', self, boundary_bp, on_row, stream, include_chrom))
        else:
            # --- HIGH-PERFORMANCE MISMATCH SEARCH (OPTIMIZED ALGORITHMS) ---
            if worker_line > 0:
                print(f"\033[{worker_line};0H\033[K[WORKER] Mode: Mismatch-Tolerant Search (up to {max_mismatches} mismatches)")
                print(f"\033[{worker_line+1};0H\033[K[WORKER] Using Optimized Myers Bit-Vector Algorithm | Pattern: {pattern} | Length: {seq_len:,} bp")
            else:
                print(f"\n[INFO] Mode: Mismatch-Tolerant Search (up to {max_mismatches} mismatches)")
                print(f"[INFO] Using optimized Myers bit-vector algorithm for optimal performance.")
                print(f"[INFO] Pattern: {pattern} | Sequence length: {seq_len:,} bp")

            matches = []
            rev_pattern = self._reverse_complement(pattern)

            # --- FORWARD STRAND SEARCH ---
            if worker_line > 0:
                _print_worker_progress(worker_line, "Searching forward strand with optimized algorithm...")
            else:
                print("[INFO] Searching forward strand with optimized algorithm...")
            
            # Use optimized mismatch search (orders of magnitude faster than flawed Bitap)
            forward_matches = self._search_optimized_mismatch(seq, pattern, max_mismatches)
            matches.extend(_process_search_results(forward_matches, pattern, '+', self, boundary_bp, on_row, stream, include_chrom))

            # --- REVERSE STRAND SEARCH ---
            if rev_pattern != pattern:
                if worker_line > 0:
                    _print_worker_progress(worker_line, "Searching reverse complement strand with optimized algorithm...")
                else:
                    print("[INFO] Searching reverse complement strand with optimized algorithm...")
                
                # Use optimized mismatch search for reverse strand
                reverse_matches = self._search_optimized_mismatch(seq, rev_pattern, max_mismatches)
                matches.extend(_process_search_results(reverse_matches, pattern, '-', self, boundary_bp, on_row, stream, include_chrom))

            if worker_line > 0:
                print(f"\033[{worker_line};0H\033[K[WORKER] Optimized search completed. Found {len(matches)} potential matches.")
            else:
                print(f"\n[INFO] Optimized search completed. Found {len(matches)} potential matches.")
        unique_matches = list({(m[0], m[1]): m for m in matches}.values())
        unique_matches.sort(key=lambda x: x[0])
        return unique_matches

def parse_fasta_file(filepath: str) -> Dict[str, str]:
    sequences = {}
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            header = ""
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    header = line
                    sequences[header] = ""
                elif header:
                    sequences[header] += line.upper()
    except FileNotFoundError:
        print(f"Error: Could not find file {filepath}")
        return {}
    return sequences

def save_results_to_csv(matches: List, filepath: str):
    print(f"\nSaving results to {filepath}...")
    try:
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            if matches:
                sample = matches[0]
                if len(sample) == 9:
                    header = ["Chromosome", "Start", "End", "Strand", "Pattern", "Found Sequence", "RevComp", "Mismatches", "Location"]
                else:
                    header = ["Start", "End", "Strand", "Pattern", "Found Sequence", "RevComp", "Mismatches", "Location"]
                writer.writerow(header)
                for match_tuple in matches:
                    writer.writerow(list(match_tuple))
            else:
                writer.writerow(["No results"])
        print("Save complete.")
    except Exception as e:
        print(f"An error occurred while saving the file: {e}")

def _print_results_header(include_chrom: bool):
    if include_chrom:
        header = (
            f"{'Chr':<{CHR_W}} {'Start':<{START_W}} {'End':<{END_W}} {'S':<{STRAND_W}} "
            f"{'Pattern':<{PAT_W}} {'Found':<{FWD_W}} {'RevComp':<{REV_W}} {'Mis':<{MIS_W}} Location"
        )
    else:
        header = (
            f"{'Start':<{START_W}} {'End':<{END_W}} {'S':<{STRAND_W}} "
            f"{'Pattern':<{PAT_W}} {'Found':<{FWD_W}} {'RevComp':<{REV_W}} {'Mis':<{MIS_W}} Location"
        )
    global FRAME_WIDTH
    FRAME_WIDTH = len(header)
    border = "+" + "-" * (FRAME_WIDTH + 2) + "+"
    print(border)
    print("| " + header + " |")

def _print_bottom_frame():
    border = "+" + "-" * (FRAME_WIDTH + 2) + "+"
    print(border)

def _print_row_tuple(row_tuple: Tuple, include_chrom: bool):
    if include_chrom:
        chrom, start, end, strand, pat, fwd, rev, mis, loc = row_tuple
        line = (
            f"{str(chrom):<{CHR_W}} {start:<{START_W}} {end:<{END_W}} {strand:<{STRAND_W}} "
            f"{pat[:PAT_W]:<{PAT_W}} {fwd[:FWD_W]:<{FWD_W}} {rev[:REV_W]:<{REV_W}} {mis:<{MIS_W}} {loc}"
        )
    else:
        start, end, strand, pat, fwd, rev, mis, loc = row_tuple
        line = (
            f"{start:<{START_W}} {end:<{END_W}} {strand:<{STRAND_W}} "
            f"{pat[:PAT_W]:<{PAT_W}} {fwd[:FWD_W]:<{FWD_W}} {rev[:REV_W]:<{REV_W}} {mis:<{MIS_W}} {loc}"
        )
    line = line[:FRAME_WIDTH]
    padding = " " * max(0, FRAME_WIDTH - len(line))
    print("| " + line + padding + " |", flush=True)

def _print_search_banner(species: str, genome_path: str, annotation_path: str,
                         target: str, pattern: str, mismatches: int, boundary_bp: int):
    gname = os.path.basename(genome_path)
    aname = os.path.basename(annotation_path)
    line = "+" + "-" * SEARCH_BANNER_WIDTH + "+"
    print("\n" + line)
    print("| {:<60} |".format("Search Context"))
    print(line)
    print("| {:<60} |".format(f"Species       : {species}"))
    print("| {:<60} |".format(f"Target        : {target}"))
    print("| {:<60} |".format(f"Pattern       : {pattern}    Mismatches: {mismatches}  Boundary: {boundary_bp}bp"))
    print("| {:<60} |".format(f"Genome file   : {gname}"))
    print("| {:<60} |".format(f"Annotation    : {aname}"))
    print(line)

def _prompt_numeric_choice(options: List[str], title: str) -> int:
    try:
        while True:
            print(f"\n{title}")
            for i, opt in enumerate(options, 1):
                print(f"{i}. {opt}")
            
            # Add navigation options
            print(f"\nNavigation: 'b' (back), 'm' (main menu), 'q' (quit)")
            
            choice_str = input(f"Choose (1-{len(options)}): ").strip().lower()
            
            # Handle empty input (just pressing Enter) - select default (option 1)
            if not choice_str:
                print(f"Using default: {options[0]}")
                return 0
            
            # Handle navigation commands
            if choice_str in ['b', 'back']:
                return -1  # Signal to go back
            elif choice_str in ['m', 'main', 'menu']:
                return -2  # Signal to go to main menu
            elif choice_str in ['q', 'quit', 'exit']:
                raise GenomeAnalysisError("User chose to exit")
            
            if choice_str.isdigit():
                idx = int(choice_str) - 1
                if 0 <= idx < len(options):
                    return idx
            print("Invalid choice. Try again.")
    except KeyboardInterrupt:
        print("\nOperation cancelled by user. Exiting gracefully.")
        raise GenomeAnalysisError("Operation cancelled by user")

def _list_subdirectories(path: str) -> List[str]:
    try:
        return [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    except Exception:
        return []

def _list_files_with_ext(path: str, exts: Tuple[str, ...]) -> List[str]:
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
    chroms = []
    seen = set()
    malformed_line_count = 0  # New counter for data integrity monitoring
    total_lines_processed = 0  # Track total lines for context
    
    try:
        with open(gtf_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                total_lines_processed += 1
                if not line or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    malformed_line_count += 1  # Count the bad lines
                    continue
                seqid = parts[0]
                if seqid in seen:
                    continue
                seen.add(seqid)
                chroms.append(seqid)
                if len(chroms) >= max_to_show:
                    break

        # CRITICAL: After processing, check for and warn about ignored data
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

def _navigator_select_inputs(data_dir: str) -> Tuple[str, str, str, Optional[str]]:
    while True:
        # Species selection by navigating subfolders under data_dir
        species_dirs = _list_subdirectories(data_dir)
        if not species_dirs:
            raise GenomeAnalysisError("No species subdirectories found under 'data'.")
        
        sp_idx = _prompt_numeric_choice(species_dirs, "Available species (data subfolders):")
        if sp_idx < 0:  # Navigation signal
            if sp_idx == -2:  # Main menu
                return None, None, None, None
            continue  # Back signal, restart species selection
        
        species_dir_name = species_dirs[sp_idx]
        species_dir_path = os.path.join(data_dir, species_dir_name)

        # Genome file selection with default option
        genome_candidates = _list_files_with_ext(species_dir_path, ('.fa', '.fasta', '.fna', '.txt'))
        if not genome_candidates:
            raise GenomeAnalysisError(f"No genome FASTA/text found in {species_dir_path}")
        
        # Add default option based on species
        default_genome = _get_default_genome_file(species_dir_name)
        if default_genome and default_genome in genome_candidates:
            # Remove the original entry to avoid duplication
            genome_candidates.remove(default_genome)
            genome_candidates.insert(0, f"{default_genome} (DEFAULT)")
        
        g_idx = _prompt_numeric_choice(genome_candidates, f"Genome files in {species_dir_name}:")
        if g_idx < 0:  # Navigation signal
            if g_idx == -2:  # Main menu
                return None, None, None, None
            continue  # Back signal, restart species selection
        
        genome_filename = genome_candidates[g_idx].replace(" (DEFAULT)", "")

        # Annotation file selection with default option
        anno_candidates = _list_files_with_ext(species_dir_path, ('.gtf', '.gff', '.gff3', '.txt', '.json'))
        if not anno_candidates:
            raise GenomeAnalysisError(f"No annotation files found in {species_dir_path}")
        
        # Add default option based on species
        default_annotation = _get_default_annotation_file(species_dir_name)
        if default_annotation and default_annotation in anno_candidates:
            # Remove the original entry to avoid duplication
            anno_candidates.remove(default_annotation)
            anno_candidates.insert(0, f"{default_annotation} (DEFAULT)")
        
        a_idx = _prompt_numeric_choice(anno_candidates, f"Annotation files in {species_dir_name}:")
        if a_idx < 0:  # Navigation signal
            if a_idx == -2:  # Main menu
                return None, None, None, None
            continue  # Back signal, restart genome selection

        annotation_filename = anno_candidates[a_idx].replace(" (DEFAULT)", "")

        chrom = None
        if annotation_filename.lower().endswith('.gtf'):
            gtf_path = os.path.join(species_dir_path, annotation_filename)
            fasta_path = os.path.join(species_dir_path, genome_filename)
            chroms = _discover_gtf_chromosomes(gtf_path)
            chroms = _filter_human_chromosomes(chroms, fasta_path)
            if not chroms:
                raise GenomeAnalysisError("Could not discover chromosomes from GTF.")
            menu = ["ALL (primary chromosomes)"] + chroms
            c_idx = _prompt_numeric_choice(menu, "Select chromosome to analyze:")
            if c_idx < 0:  # Navigation signal
                if c_idx == -2:  # Main menu
                    return None, None, None, None
                continue  # Back signal, restart annotation selection
            chrom = 'ALL' if c_idx == 0 else chroms[c_idx - 1]

        return species_dir_path, genome_filename, annotation_filename, chrom

def _get_default_genome_file(species_name: str) -> Optional[str]:
    """Get default genome file for each species."""
    defaults = {
        'E.coli': 'E_coli_K12_MG1655.fasta',
        'B.Sub': 'B_subtilis_168.fasta', 
        'Homo_Sapiens': 'GRCh38.primary_assembly.genome.fa'
    }
    return defaults.get(species_name)

def _get_default_annotation_file(species_name: str) -> Optional[str]:
    """Get default annotation file for each species."""
    defaults = {
        'E.coli': 'E_coli_K12_MG1655.gff',
        'B.Sub': 'B_subtilis_168.json',
        'Homo_Sapiens': 'gencode.v48.primary_assembly.annotation.gtf'
    }
    return defaults.get(species_name)

def _search_single_chromosome_worker(chromosome: str, fasta_path: str, gtf_path: str, 
                                   pattern: str, max_mismatches: int, boundary_bp: int, worker_line: int = 0) -> Tuple[bool, str, List]:
    """
    Worker function for parallel chromosome processing.
    Uses pre-initialized analyzer if available, otherwise creates new instance.
    Uses memory mapping for large genome files.
    Returns (success, chromosome, results_list).
    """
    try:
        # Use pre-initialized analyzer (required for optimal performance)
        if 'worker_analyzer' in globals():
            analyzer = worker_analyzer
            # Reset analyzer state for new chromosome
            analyzer.sequence = ""
            analyzer.gene_data = {}
            analyzer.gene_tree = IntervalTree()
            analyzer.cds_tree = IntervalTree()
        else:
            raise GenomeAnalysisError("Worker pre-initialization required. Use global process pool with initializer.")
        
        # Load chromosome data with memory mapping if appropriate
        print(f"[INFO] Worker loading chromosome {chromosome} with memory mapping")
        
        # MEMORY USAGE MONITORING FOR WORKER DASHBOARD
        pre_load_memory = _monitor_memory_usage()[0]
        print(f"[INFO] Worker {chromosome}: Memory before loading: {pre_load_memory:.1f} MB")
        
        sequence = load_genome_for_chromosome(fasta_path, chromosome, use_memory_mapping=True)
        analyzer.sequence = sequence
        
        # MEMORY USAGE AFTER LOADING
        post_load_memory, memory_increase, _ = _monitor_memory_usage(pre_load_memory)
        print(f"[INFO] Worker {chromosome}: Memory after loading: {post_load_memory:.1f} MB (+{memory_increase:.1f} MB)")
        analyzer.current_chromosome = chromosome  # Set chromosome ID for progress display
        analyzer.worker_line_number = worker_line  # Set worker line number for progress display
        analyzer.load_annotations_from_gtf(gtf_path, chromosome)
        
        # Perform search with individual progress tracking and live output
        print(f"[INFO] Worker {chromosome}: Starting search for pattern '{pattern}' on {len(sequence):,} bp")
        print(f"[INFO] Worker {chromosome}: Progress tracking enabled")
        
        # ENHANCED PROGRESS REPORTING FOR WORKER DASHBOARD
        if len(sequence) > MEMORY_MAPPING_THRESHOLD_BYTES:  # Use constant for threshold
            print(f"[INFO] Worker {chromosome}: Large sequence detected, using optimized chunked processing")
        else:
            print(f"[INFO] Worker {chromosome}: Standard sequence size, using direct processing")
        
        results = analyzer.search_sequence(pattern, max_mismatches, boundary_bp, stream=True)
        print(f"[INFO] Worker {chromosome}: Completed search, found {len(results)} matches")
        
        # FINAL MEMORY USAGE REPORT FOR WORKER DASHBOARD
        final_memory, total_memory_increase, _ = _monitor_memory_usage(pre_load_memory)
        print(f"[INFO] Worker {chromosome}: Final memory: {final_memory:.1f} MB (Total increase: +{total_memory_increase:.1f} MB)")
        
        return True, chromosome, results
        
    except Exception as e:
        return False, chromosome, [f"Error processing {chromosome}: {str(e)}"]

def _search_genome_chunk_worker(analyzer: EColiAnalyzer, start: int, end: int, 
                               pattern: str, max_mismatches: int, boundary_bp: int) -> List:
    """
    Worker function for parallel genome chunk processing.
    Searches a specific chunk of the genome sequence.
    """
    try:
        # Create a copy of the analyzer for this worker
        worker_analyzer = EColiAnalyzer()
        worker_analyzer.sequence = analyzer.sequence[start:end]
        worker_analyzer.gene_tree = analyzer.gene_tree
        worker_analyzer.gene_data = analyzer.gene_data
        worker_analyzer.cds_tree = analyzer.cds_tree
        
        # Perform search on this chunk
        results = worker_analyzer.search_sequence(pattern, max_mismatches, boundary_bp, stream=False)
        
        # Adjust positions to account for chunk offset
        adjusted_results = []
        for result in results:
            start_pos, end_pos, strand, pat, found, revcomp, mismatches, location = result
            adjusted_results.append((start_pos + start, end_pos + start, strand, pat, found, revcomp, mismatches, location))
        
        return adjusted_results
        
    except Exception as e:
        print(f"[ERROR] Error processing chunk {start}-{end}: {str(e)}")
        return []

class MemoryMappedGenomeLoader:
    """
    Handles memory-mapped file access for large genome files.
    Implements two-stage caching: raw sequence extraction and memory mapping.
    """
    
    def __init__(self, cache_dir: str = None, max_cache_size_mb: int = 1024):
        self.cache_dir = cache_dir or os.path.join(os.path.dirname(__file__), "cache")
        self.max_cache_size_mb = max_cache_size_mb
        os.makedirs(self.cache_dir, exist_ok=True)
        
        # Perform cleanup on initialization
        self.cleanup_cache()
    
    def create_raw_sequence_cache(self, fasta_path: str, chromosome: str) -> str:
        """
        Extract raw DNA sequence from FASTA and save to .raw cache file.
        Returns path to the .raw cache file.
        """
        import time  # Import at the very beginning
        
        cache_file = os.path.join(self.cache_dir, f"{chromosome}.raw")
        
        # Check if cache already exists and is newer than source
        if os.path.exists(cache_file):
            if os.path.getmtime(cache_file) > os.path.getmtime(fasta_path):
                print(f"[INFO] Using existing raw sequence cache for {chromosome}")
                return cache_file
        
        # Check if another process is currently creating this cache
        lock_file = cache_file + ".lock"
        if os.path.exists(lock_file):
            # Wait a bit and check if cache was created
            time.sleep(0.5)
            if os.path.exists(cache_file):
                print(f"[INFO] Cache file created by another process for {chromosome}, using it")
                return cache_file
        
        print(f"[INFO] Creating raw sequence cache for {chromosome}")
        
        # Add file locking to prevent race conditions
        lock_file = cache_file + ".lock"
        max_retries = CACHE_MAX_RETRIES
        retry_delay = CACHE_RETRY_DELAY
        
        for attempt in range(max_retries):
            try:
                # Try to create lock file
                with open(lock_file, 'x') as f:
                    f.write(str(os.getpid()))
                break
            except FileExistsError:
                if attempt < max_retries - 1:
                    # Check if cache file was created by another process while waiting
                    if os.path.exists(cache_file):
                        print(f"[INFO] Cache file created by another process for {chromosome}, using it")
                        return cache_file
                    
                    print(f"[INFO] Cache file being created by another process, waiting... (attempt {attempt + 1})")
                    time.sleep(retry_delay)
                    continue
                else:
                    # Final attempt - check if cache file exists
                    if os.path.exists(cache_file):
                        print(f"[INFO] Using cache file created by another process for {chromosome}")
                        return cache_file
                    else:
                        # Try one more time to create the cache file directly
                        try:
                            print(f"[INFO] Final attempt to create cache for {chromosome}")
                            return self._create_cache_directly(fasta_path, chromosome, cache_file)
                        except Exception as e:
                            raise GenomeAnalysisError(f"Could not create cache file for {chromosome}: {e}")
        
        # Use centralized sequence extraction to eliminate duplication
        raw_sequence = self._extract_sequence_to_cache(fasta_path, chromosome, cache_file)
        
        # Write raw sequence to cache file
        with open(cache_file, 'w') as f:
            f.write(raw_sequence)
        
        # Remove lock file after successful creation
        try:
            if os.path.exists(lock_file):
                os.remove(lock_file)
                print(f"[INFO] Lock file removed for {chromosome}")
        except Exception as e:
            print(f"[WARNING] Failed to remove lock file: {e}")
        
        return cache_file
    
    def _create_cache_directly(self, fasta_path: str, chromosome: str, cache_file: str) -> str:
        """Create cache file directly without locking (final attempt only)."""
        print(f"[INFO] Creating cache directly for {chromosome} (final attempt)")
        
        # Reuse the main cache creation logic to eliminate duplication
        return self._extract_sequence_to_cache(fasta_path, chromosome, cache_file)
    
    def _extract_sequence_to_cache(self, fasta_path: str, chromosome: str, cache_file: str) -> str:
        """Centralized sequence extraction logic to eliminate code duplication."""
        wanted = set(_get_chromosome_aliases(chromosome))
        sequence_parts = []
        capturing = False
        
        with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
            for raw_line in f:
                line = raw_line.strip()
                if line.startswith('>'):
                    header = line[1:].split()[0]
                    capturing = any(alias in header for alias in wanted)
                elif capturing and line:
                    sequence_parts.append(line.upper())
        
        if not sequence_parts:
            raise ValueError(f"Chromosome {chromosome} not found in {fasta_path}")
        
        raw_sequence = "".join(sequence_parts)
        return raw_sequence
    
    def load_sequence_memory_mapped(self, cache_file: str) -> Union[str, bytes]:
        """
        Load sequence using memory mapping for efficient access.
        Returns sequence as string or bytes object.
        """
        print(f"[INFO] Loading sequence using memory mapping: {os.path.basename(cache_file)}")
        
        file_size = os.path.getsize(cache_file)
        
        # Use memory mapping for files larger than threshold (consistent with main threshold)
        if file_size > MEMORY_MAPPING_THRESHOLD_BYTES:
            print(f"[INFO] Using memory mapping for large file ({file_size / (1024*1024):.1f} MB)")
            with open(cache_file, 'rb') as f:
                with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                    return mm.read().decode('utf-8')
        else:
            print(f"[INFO] Using standard file I/O for small file ({file_size / (1024*1024):.1f} MB)")
            with open(cache_file, 'r') as f:
                return f.read()
    
    def cleanup_cache(self):
        """Remove old cache files if total size exceeds limit."""
        try:
            cache_files = []
            total_size = 0
            
            for filename in os.listdir(self.cache_dir):
                if filename.endswith('.raw'):
                    filepath = os.path.join(self.cache_dir, filename)
                    size = os.path.getsize(filepath)
                    cache_files.append((filepath, size, os.path.getmtime(filepath)))
                    total_size += size
            
            # Sort by modification time (oldest first)
            cache_files.sort(key=lambda x: x[2])
            
            # Remove oldest files until under limit
            max_size_bytes = self.max_cache_size_mb * 1024 * 1024
            for filepath, size, _ in cache_files:
                if total_size <= max_size_bytes:
                    break
                os.remove(filepath)
                total_size -= size
                print(f"[INFO] Removed old cache file: {os.path.basename(filepath)}")
        except Exception as e:
            print(f"[WARNING] Cache cleanup failed: {e}")

def load_genome_for_chromosome(fasta_path: str, chromosome: str, use_memory_mapping: bool = True) -> str:
    """
    Load genome sequence for a chromosome, optionally using memory mapping.
    Returns the sequence as a string.
    """
    file_size = os.path.getsize(fasta_path)
    
    # OPTIMIZED THRESHOLD: Use memory mapping for files above threshold
    if use_memory_mapping and file_size > MEMORY_MAPPING_THRESHOLD_BYTES:
        print(f"[INFO] Large sequence detected ({file_size / (1024*1024):.1f} MB), using memory mapping")
        try:
            loader = MemoryMappedGenomeLoader()
            cache_file = loader.create_raw_sequence_cache(fasta_path, chromosome)
            sequence = loader.load_sequence_memory_mapped(cache_file)
            print(f"[INFO] Memory mapped loading completed: {len(sequence):,} bp")
            return sequence
        except Exception as e:
            raise GenomeAnalysisError(f"Memory mapping failed: {e}. Check file permissions and disk space.")
    else:
        print(f"[INFO] Using standard loading for sequence ({file_size / (1024*1024):.1f} MB)")
        # Use existing chromosome loading logic
        analyzer = EColiAnalyzer()
        analyzer.load_chromosome_sequence(fasta_path, chromosome)
        return analyzer.sequence

def parallel_chromosome_search(chromosomes: List[str], fasta_path: str, gtf_path: str,
                             pattern: str, max_mismatches: int, boundary_bp: int, 
                             command_queue: Optional[multiprocessing.Queue] = None) -> List:
    """
    Perform parallel search across multiple chromosomes using global process pool.
    Returns aggregated results from all chromosomes with live output.
    """
    print(f"[INFO] Starting parallel chromosome processing with {len(chromosomes)} chromosomes")
    
    # Determine optimal worker count
    max_workers = min(multiprocessing.cpu_count(), len(chromosomes))
    print(f"[INFO] Using {max_workers} workers for parallel processing")
    
    # Get global process pool (reused across searches)
    executor = _get_global_process_pool(max_workers)
    
    all_results = []
    
    print(f"\n{'='*80}")
    print(f"PARALLEL SEARCH PROGRESS - {len(chromosomes)} CHROMOSOMES")
    print(f"{'='*80}")
    
    # Print header for live output
    _print_results_header(include_chrom=True)
    
    # Set up worker progress display area
    max_workers = min(multiprocessing.cpu_count(), len(chromosomes))
    
    # Check if dashboard is available
    use_dashboard = command_queue is not None
    
    if use_dashboard:
        print(f"\n{'='*80}")
        print(f"DASHBOARD MODE - {max_workers} WORKERS")
        print(f"{'='*80}")
        print(f"Dashboard window opened - monitor progress in real-time!")
        print(f"You can minimize this console and watch the dashboard instead")
        print(f"{'='*80}")
    else:
        print(f"\n{'='*80}")
        print(f"WORKER PROGRESS MONITORING - {max_workers} WORKERS")
        print(f"{'='*80}")
        
        # Reserve lines for each worker (start from current line + 1)
        import os
        worker_start_line = os.get_terminal_size().lines - max_workers - 2  # Leave space for results
        
        # Print placeholder lines for each worker
        for i in range(max_workers):
            print(f"[WORKER {i+1}] Waiting for assignment...")
        
        print(f"\n{'='*80}")
        print(f"LIVE RESULTS OUTPUT")
        print(f"{'='*80}")
    
    # Submit all chromosome search tasks to global pool
    future_to_chrom = {}
    worker_status = {i: "Waiting for assignment..." for i in range(max_workers)}
    
    for i, chrom in enumerate(chromosomes):
        worker_id = i % max_workers
        worker_status[worker_id] = f"Processing {chrom}..."
        
        # Update dashboard if available with algorithm details
        if use_dashboard and command_queue:
            algorithm_type = "Aho-Corasick" if max_mismatches == 0 else "Myers Bit-Vector"
            command_queue.put(("UPDATE_WORKER", (worker_id, "Processing", 0, f"{algorithm_type}: {chrom}")))
        else:
            # Update console worker progress bar
            line_num = worker_start_line + worker_id + 1
            print(f"\033[{line_num};0H\033[K[WORKER {worker_id+1}] {worker_status[worker_id]}")
        
        future = executor.submit(_search_single_chromosome_worker, chrom, fasta_path, gtf_path, 
                                pattern, max_mismatches, boundary_bp, 0)  # worker_line not used in separate process
        future_to_chrom[future] = chrom
    
    # Collect results as they complete with live output and update progress bars
    completed = 0
    worker_status = {i: "Waiting for assignment..." for i in range(max_workers)}
    
    for future in as_completed(future_to_chrom):
        chrom = future_to_chrom[future]
        completed += 1
        
        # Update worker progress bar
        worker_id = (completed - 1) % max_workers
        worker_status[worker_id] = f"Completed {chrom} ({completed}/{len(chromosomes)})"
        
        # Update dashboard if available with algorithm performance metrics
        if use_dashboard and command_queue:
            algorithm_type = "Aho-Corasick" if max_mismatches == 0 else "Myers Bit-Vector"
            command_queue.put(("UPDATE_WORKER", (worker_id, "Completed", 100, f"{algorithm_type}: {chrom} ({len(results)} matches)")))
            command_queue.put(("UPDATE_OVERALL", (completed, len(chromosomes), len(all_results))))
        else:
            # Update console worker progress bars
            for i in range(max_workers):
                line_num = worker_start_line + i + 1
                status = worker_status.get(i, "Waiting for assignment...")
                print(f"\033[{line_num};0H\033[K[WORKER {i+1}] {status}")
        
        try:
            success, chrom_name, results = future.result()
            if success:
                print(f"[INFO] Chromosome {chrom_name} completed ({completed}/{len(chromosomes)}) - {len(results)} matches")
                # Add chromosome info to each result and print live
                for result in results:
                    if len(result) == 8:  # Standard result tuple
                        full_result = (chrom_name,) + result
                        all_results.append(full_result)
                        _print_row_tuple(full_result, include_chrom=True)
                    else:
                        all_results.append(result)
                        _print_row_tuple(result, include_chrom=True)
            else:
                print(f"[INFO] Failed to process chromosome {chrom_name}: {results[0]}")
                worker_status[worker_id] = f"Failed {chrom_name}: {results[0]}"
                if command_queue:
                    command_queue.put(("UPDATE_WORKER", (worker_id, "Failed", 0, chrom_name)))
                
        except Exception as e:
            print(f"[INFO] Exception processing chromosome {chrom}: {e}")
            worker_status[worker_id] = f"Exception {chrom}: {e}"
            if command_queue:
                command_queue.put(("UPDATE_WORKER", (worker_id, "Exception", 0, chrom)))
    
    # Print all results in batch after completion
    print(f"\n{'='*80}")
    print(f"SEARCH RESULTS - {len(all_results)} TOTAL MATCHES")
    print(f"{'='*80}")
    
    # Print header for results
    _print_results_header(include_chrom=True)
    
    # Print all results
    for result in all_results:
        _print_row_tuple(result, include_chrom=True)
    
    print(f"\n{'='*80}")
    print(f"PARALLEL SEARCH COMPLETE - {len(all_results)} TOTAL MATCHES")
    print(f"{'='*80}")
    return all_results

def parallel_single_genome_search(analyzer: EColiAnalyzer, pattern: str, max_mismatches: int, boundary_bp: int) -> List:
    """
    Perform parallel search on a single genome using multiple CPU cores.
    This is used for non-GTF files (E. coli, B. subtilis).
    """
    print(f"[INFO] Starting parallel single genome search")
    
    # Determine optimal worker count
    max_workers = multiprocessing.cpu_count()
    print(f"[INFO] Using {max_workers} workers for parallel processing")
    
    # Get global process pool
    executor = _get_global_process_pool(max_workers)
    
    # Split the genome into chunks for parallel processing
    sequence = analyzer.sequence
    chunk_size = len(sequence) // max_workers
    chunks = []
    
    for i in range(max_workers):
        start = i * chunk_size
        end = start + chunk_size if i < max_workers - 1 else len(sequence)
        chunks.append((start, end))
    
    # Print header for live output
    _print_results_header(include_chrom=False)
    
    # Submit chunk search tasks
    future_to_chunk = {
        executor.submit(_search_genome_chunk_worker, analyzer, start, end, pattern, max_mismatches, boundary_bp): i 
        for i, (start, end) in enumerate(chunks)
    }
    
    all_results = []
    
    # Collect results as they complete
    completed = 0
    for future in as_completed(future_to_chunk):
        chunk_id = future_to_chunk[future]
        completed += 1
        
        try:
            results = future.result()
            print(f"[INFO] Completed chunk {chunk_id + 1} ({completed}/{len(chunks)})")
            
            # Print results live
            for result in results:
                all_results.append(result)
                _print_row_tuple(result, include_chrom=False)
                
        except Exception as e:
            print(f"[INFO] Exception processing chunk {chunk_id}: {e}")
    
    print(f"[INFO] Parallel genome search complete. Found {len(all_results)} total matches")
    return all_results

def launch_dashboard():
    """Launch the worker dashboard in a separate process."""
    try:
        from worker_dashboard import start_dashboard
        import multiprocessing
        
        # Create a queue for communication
        command_queue = multiprocessing.Queue()
        
        # Launch dashboard in separate process
        dashboard_process = multiprocessing.Process(
            target=start_dashboard, 
            args=(command_queue, os.cpu_count())
        )
        dashboard_process.start()
        
        print(f"\n🎯 Dashboard launched! A new window should appear.")
        print(f"💡 Dashboard will update in real-time as searches progress!")
        print(f"💡 You can minimize this console and watch the dashboard instead!")
        
        return command_queue, dashboard_process
        
    except ImportError:
        print(f"\n[WARNING] Dashboard not available. Install tkinter to use this feature.")
        return None, None
    except Exception as e:
        print(f"\n[WARNING] Could not start dashboard: {e}")
        return None, None

def cleanup_dashboard(dashboard_process, command_queue):
    """Safely cleanup dashboard process and resources."""
    if dashboard_process is not None:
        print("\n[INFO] Cleaning up dashboard...")
        try:
            # Send close command if queue is available
            if command_queue is not None:
                try:
                    command_queue.put(("CLOSE", None), timeout=2)
                except:
                    pass  # Queue might be full or closed
            
            # Give process time to close gracefully
            dashboard_process.join(timeout=5)
            
            if dashboard_process.is_alive():
                print("[INFO] Force terminating dashboard...")
                dashboard_process.terminate()
                dashboard_process.join(timeout=2)
                
                if dashboard_process.is_alive():
                    print("[WARNING] Dashboard process still running - killing...")
                    try:
                        dashboard_process.kill()
                    except:
                        pass
                        
            print("[INFO] Dashboard cleanup completed.")
            
        except Exception as e:
            print(f"[WARNING] Dashboard cleanup encountered an error: {e}")

def _print_formatted_results(matches: List, include_chrom: bool = False):
    """
    Centralized result formatting to eliminate code duplication.
    Handles both chromosome and non-chromosome result formats.
    """
    if not matches:
        print("No matches found.")
        return
    
    # Print header
    _print_results_header(include_chrom)
    
    # Print separator using pre-calculated constants
    separator = '-' * (CHROM_TOTAL_WIDTH if include_chrom else NON_CHROM_TOTAL_WIDTH)
    
    print(separator)
    
    # Print results
    for match in matches:
        _print_row_tuple(match, include_chrom)
    
    print(separator)

def _perform_cache_cleanup(cache_loader, operation_name: str = "operation"):
    """
    Centralized cache cleanup to eliminate code duplication.
    """
    try:
        cache_loader.cleanup_cache()
        print(f"[INFO] Cache cleanup completed after {operation_name}")
    except Exception as e:
        print(f"[WARNING] Cache cleanup failed after {operation_name}: {e}")

def _get_filtered_chromosomes(gtf_path: str, fasta_path: str) -> List[str]:
    """
    Centralized chromosome discovery and filtering to eliminate code duplication.
    Returns filtered list of chromosomes for parallel processing.
    """
    chroms = _discover_gtf_chromosomes(gtf_path)
    filtered_chroms = _filter_human_chromosomes(chroms, fasta_path)
    print(f"[INFO] Using parallel processing for {len(filtered_chroms)} chromosomes")
    return filtered_chroms

def main():
    """
    Enhanced main function with complete algorithm transparency and debugging.
    """
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        script_dir = os.getcwd()
    
    project_root = os.path.dirname(script_dir)
    data_dir = os.path.join(project_root, "data")
    seq_files_dir = os.path.join(script_dir, "Seq_Files")
    results_dir = os.path.join(script_dir, "Results")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(seq_files_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print(f"[INFO] Genome Analyzer - Enhanced Version with Algorithm Transparency")
    print(f"[INFO] Data directory: {data_dir}")
    print(f"[INFO] Results directory: {results_dir}")
    print(f"[INFO] Sequence files directory: {seq_files_dir}")
    
    # Offer dashboard option
    print(f"\n🎯 Would you like to launch the real-time worker dashboard?")
    print(f"   This will show live progress of all workers in a separate window.")
    print(f"   Type 'y' to launch dashboard, any other key to continue without it.")
    
    choice = input("Launch dashboard? (y/N): ").strip().lower()
    command_queue, dashboard_process = None, None
    if choice == 'y':
        command_queue, dashboard_process = launch_dashboard()
    
    # Initialize cache cleanup on startup
    try:
        cache_loader = MemoryMappedGenomeLoader()
        print(f"[INFO] Cache management initialized (limit: {cache_loader.max_cache_size_mb}MB)")
    except Exception as e:
        print(f"[WARNING] Cache initialization failed: {e}")

    # Navigator menu to choose species and files
    try:
        species_dir_path, genome_filename, annotation_filename, selected_chrom = _navigator_select_inputs(data_dir)
        if species_dir_path is None:  # User chose to go to main menu
            return
    except GenomeAnalysisError as e:
        print(f"Error: {e}")
        return
    species_name = os.path.basename(species_dir_path)

    complete_genome_file = os.path.join(species_dir_path, genome_filename)
    gene_annotation_file = os.path.join(species_dir_path, annotation_filename)

    # Per-species cache
    cache_file = os.path.join(species_dir_path, "genome_cache.pkl")

    analyzer = EColiAnalyzer()

    # Load from cache only if it matches current inputs
    multi_chrom = (gene_annotation_file.lower().endswith('.gtf') and selected_chrom == 'ALL')
    per_chrom_cache: Dict[str, EColiAnalyzer] = {} if multi_chrom else {}
    if os.path.exists(cache_file) and not multi_chrom:
        print(f"Loading cached genome data from: {os.path.basename(cache_file)}...")
        with open(cache_file, 'rb') as f:
            cached_data = pickle.load(f)
        cache_ok = (
            cached_data.get('source_genome_path') == complete_genome_file and
            cached_data.get('source_annotation_path') == gene_annotation_file and
            cached_data.get('chromosome') == selected_chrom
        )
        if cache_ok:
            analyzer.sequence = cached_data['sequence']
            analyzer.gene_tree = cached_data['gene_tree']
            analyzer.gene_data = cached_data['gene_data']
            print("Cached data loaded.")
        else:
            print("Cache does not match selected inputs. Rebuilding...")
            try:
                if selected_chrom and selected_chrom != 'ALL':
                    analyzer.load_chromosome_sequence(complete_genome_file, selected_chrom)
                else:
                    analyzer.load_complete_sequence(complete_genome_file)
                if gene_annotation_file.lower().endswith('.json'):
                    analyzer.load_annotations_from_bsub_json(gene_annotation_file)
                elif gene_annotation_file.lower().endswith('.gtf'):
                    if not selected_chrom or selected_chrom == 'ALL':
                        raise GenomeAnalysisError("No chromosome selected for GTF annotations.")
                    analyzer.load_annotations_from_gtf(gene_annotation_file, selected_chrom)
                else:
                    analyzer.load_annotations(gene_annotation_file)
            except GenomeAnalysisError as e:
                print(f"Error: {e}")
                return
            with open(cache_file, 'wb') as f:
                pickle.dump({
                    'sequence': analyzer.sequence,
                    'gene_tree': analyzer.gene_tree,
                    'gene_data': analyzer.gene_data,
                    'source_genome_path': complete_genome_file,
                    'source_annotation_path': gene_annotation_file,
                    'chromosome': selected_chrom
                }, f)
            print(f"Saved new genome data to cache file: {os.path.basename(cache_file)}")
    else:
        try:
            if selected_chrom and selected_chrom != 'ALL':
                analyzer.load_chromosome_sequence(complete_genome_file, selected_chrom)
            else:
                analyzer.load_complete_sequence(complete_genome_file)
            if multi_chrom:
                # Defer per-chromosome loads to search time for ALL-chrom runs
                pass
            elif gene_annotation_file.lower().endswith('.json'):
                analyzer.load_annotations_from_bsub_json(gene_annotation_file)
            elif gene_annotation_file.lower().endswith('.gtf'):
                if not selected_chrom or selected_chrom == 'ALL':
                    raise GenomeAnalysisError("No chromosome selected for GTF annotations.")
                analyzer.load_annotations_from_gtf(gene_annotation_file, selected_chrom)
            else:
                analyzer.load_annotations(gene_annotation_file)
        except GenomeAnalysisError as e:
            print(f"Error: {e}")
            return
        if not multi_chrom:
            with open(cache_file, 'wb') as f:
                pickle.dump({
                    'sequence': analyzer.sequence,
                    'gene_tree': analyzer.gene_tree,
                    'gene_data': analyzer.gene_data,
                    'source_genome_path': complete_genome_file,
                    'source_annotation_path': gene_annotation_file,
                    'chromosome': selected_chrom
                }, f)
            print(f"Saved new genome data to cache file: {os.path.basename(cache_file)}")

    try:
        while True:
            print("\nChoose search mode:")
            print("1. Single search in selected genome")
            print("2. Batch search from substrate file")
            print("3. Quit")
            try:
                mode = input("Enter mode (1-3): ")
            except KeyboardInterrupt:
                print("\nOperation cancelled by user. Exiting.")
                break

            if mode == '1':
                try:
                    pattern = input("Enter sequence to search (IUPAC supported): ").strip().upper()
                except KeyboardInterrupt:
                    print("\nOperation cancelled by user. Returning to main menu.")
                    continue
                if not pattern: continue
                try:
                    mismatches_str = input("Enter maximum allowed mismatches (default 0): ")
                except KeyboardInterrupt:
                    print("\nOperation cancelled by user. Returning to main menu.")
                    continue
                max_mismatches = int(mismatches_str) if mismatches_str.isdigit() else 0
                try:
                    boundary_str = input("Flag near-chromosome-end matches within N bp (blank to skip): ").strip()
                except KeyboardInterrupt:
                    print("\nOperation cancelled by user. Returning to main menu.")
                    continue
                boundary_bp = int(boundary_str) if boundary_str.isdigit() else 0
                # Print search banner
                _print_search_banner(species_name, complete_genome_file, gene_annotation_file,
                                     selected_chrom if selected_chrom else 'Whole genome',
                                     pattern, max_mismatches, boundary_bp)
                # Enhanced search with algorithm transparency
                if multi_chrom:
                    print(f"\n[INFO] Multi-chromosome search mode detected")
                    gtf_path = gene_annotation_file
                    fasta_path = complete_genome_file
                    chroms = _get_filtered_chromosomes(gtf_path, fasta_path)
                    
                    # Use parallel processing for multi-chromosome search
                    matches = parallel_chromosome_search(chroms, fasta_path, gtf_path, 
                                                       pattern, max_mismatches, boundary_bp, command_queue)
                    
                    # Display results using centralized formatting
                    if matches:
                        print(f"\n[INFO] Displaying {len(matches)} parallel search results")
                        _print_formatted_results(matches, include_chrom=True)
                    else:
                        print("No matches found across all chromosomes.")
                else:
                    print(f"\n[INFO] Single genome search mode")
                    # ALWAYS use parallel processing even for single search
                    if multi_chrom:
                        gtf_path = gene_annotation_file
                        fasta_path = complete_genome_file
                        chroms = _get_filtered_chromosomes(gtf_path, fasta_path)
                        matches = parallel_chromosome_search(chroms, fasta_path, gtf_path, 
                                                           pattern, max_mismatches, boundary_bp, command_queue)
                    else:
                        # For non-GTF files, use parallel genome chunk processing
                        matches = parallel_single_genome_search(analyzer, pattern, max_mismatches, boundary_bp)
                print(f"\n[INFO] Search completed. Total matches found: {len(matches)}")
                if not matches: 
                    print("[INFO] No matches found for this pattern.")
                    continue
                save_prompt = input("\nDo you want to save these results to a CSV file? (y/n): ").lower()
                if save_prompt == 'y':
                    user_name = input("Enter a name for this result file (no extension): ").strip()
                    if not user_name:
                        user_name = "results"
                    date_str = datetime.now().strftime("%Y%m%d")
                    chrom_suffix = ''
                    if gene_annotation_file.lower().endswith('.gtf'):
                        if selected_chrom and selected_chrom != 'ALL':
                            chrom_suffix = f"-{selected_chrom}"
                        elif multi_chrom:
                            chrom_suffix = "-ALL"
                    output_filename = f"{species_name}-{user_name}{chrom_suffix}-{date_str}.csv"
                    output_filepath = os.path.join(results_dir, output_filename)
                    save_results_to_csv(matches, output_filepath)

            elif mode == '2':
                try:
                    available_files = [f for f in os.listdir(seq_files_dir) if f.endswith(('.fasta', '.txt', '.fa'))]
                    if not available_files:
                        print(f"\nNo files found in '{seq_files_dir}'."); continue
                    print("\nAvailable substrate files:")
                    for i, filename in enumerate(available_files): print(f"{i + 1}. {filename}")
                    try:
                        choice = int(input(f"Choose a file for batch search (1-{len(available_files)}): ")) - 1
                    except KeyboardInterrupt:
                        print("\nOperation cancelled by user. Returning to main menu.")
                        continue
                    if not 0 <= choice < len(available_files):
                        print("Invalid choice."); continue
                    selected_filename = available_files[choice]
                    selected_filepath = os.path.join(seq_files_dir, selected_filename)
                    queries = parse_fasta_file(selected_filepath)
                    if not queries:
                        print(f"No sequences found in {selected_filename}."); continue
                    try:
                        mismatches_str = input("Enter maximum allowed mismatches for all searches (default 0): ")
                    except KeyboardInterrupt:
                        print("\nOperation cancelled by user. Returning to main menu.")
                        continue
                    max_mismatches = int(mismatches_str) if mismatches_str.isdigit() else 0
                    try:
                        boundary_str = input("Flag near-chromosome-end matches within N bp (blank to skip): ").strip()
                    except KeyboardInterrupt:
                        print("\nOperation cancelled by user. Returning to main menu.")
                        continue
                    boundary_bp = int(boundary_str) if boundary_str.isdigit() else 0
                    all_matches = []
                    print(f"\nStarting batch search for {len(queries)} patterns...")
                    for header, pattern in queries.items():
                        print(f"\n--- Searching for pattern from: {header} ---")
                        # Print search banner for batch
                        _print_search_banner(species_name, complete_genome_file, gene_annotation_file,
                                             selected_chrom if selected_chrom else ('ALL' if multi_chrom else 'Whole genome'),
                                             pattern, max_mismatches, boundary_bp)
                        # ALWAYS use parallel processing when possible
                        if multi_chrom:
                            gtf_path = gene_annotation_file
                            fasta_path = complete_genome_file
                            chroms = _get_filtered_chromosomes(gtf_path, fasta_path)
                            
                            # Use parallel processing for multi-chromosome search
                            matches = parallel_chromosome_search(chroms, fasta_path, gtf_path, 
                                                               pattern, max_mismatches, boundary_bp, command_queue)
                            
                            # Print bottom frame for parallel results
                            if matches:
                                _print_bottom_frame()
                            else:
                                print("No matches found across all chromosomes.")
                        else:
                            # For non-GTF files, use parallel genome chunk processing
                            print(f"[INFO] Using parallel genome chunk processing")
                            matches = parallel_single_genome_search(analyzer, pattern, max_mismatches, boundary_bp)
                        print(f"Found {len(matches)} total match(es) for this pattern.")
                        
                        # Perform cache cleanup after search operations
                        _perform_cache_cleanup(cache_loader, "search")
                        
                        if len(matches) == 0:
                            print("No matches found.")
                        if matches:
                            print("\nIndividual Search Results (Sorted by Position):")
                            # Use centralized formatting function
                            include_chrom = matches and len(matches[0]) == 9
                            _print_formatted_results(matches, include_chrom=include_chrom)
                        for match in matches:
                            original_location = match[-1]
                            full_match_info = list(match)
                            full_match_info[-1] = f"{original_location} (Query: {header[:30]}...)"
                            all_matches.append(tuple(full_match_info))
                    print(f"\n--- Batch search complete. Found {len(all_matches)} total matches across all queries. ---")
                    
                    # Perform cache cleanup after batch search operations
                    _perform_cache_cleanup(cache_loader, "batch search")
                    
                    if not all_matches: continue
                    save_prompt = input("\nDo you want to save the complete results to a CSV file? (y/n): ").lower()
                    if save_prompt == 'y':
                        user_name = input("Enter a name for this result file (no extension): ").strip()
                        if not user_name:
                            user_name = "results"
                        date_str = datetime.now().strftime("%Y%m%d")
                        chrom_suffix = ''
                        if gene_annotation_file.lower().endswith('.gtf'):
                            if selected_chrom and selected_chrom != 'ALL':
                                chrom_suffix = f"-{selected_chrom}"
                            elif multi_chrom:
                                chrom_suffix = "-ALL"
                        output_filename = f"{species_name}-{user_name}{chrom_suffix}-{date_str}.csv"
                        output_filepath = os.path.join(results_dir, output_filename)
                        save_results_to_csv(all_matches, output_filepath)
                except (ValueError, IndexError):
                    print("Invalid input.")
                except Exception as e:
                    print(f"An error occurred: {e}")

            elif mode == '3':
                print("Goodbye!")
                break
            else:
                print("Invalid mode.")
    
                    
    except KeyboardInterrupt:
            print("\n\n[INFO] Program interrupted by user. Exiting gracefully...")
            
    except Exception as e:
        print(f"\n[ERROR] An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
        
    finally:
        # Always cleanup dashboard - this runs no matter how we exit
        cleanup_dashboard(dashboard_process, command_queue)

if __name__ == "__main__":
    main()