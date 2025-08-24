
import os
import re
import csv
import json
import pickle
import mmap
import time
from typing import Dict, List, Tuple, Union

from .analyzer import GenomeAnalyzer, GenomeAnalysisError

MEMORY_MAPPING_THRESHOLD_BYTES = 64 * 4096  # 64 pages = 256KB (scientifically validated OS page boundary)
CACHE_MAX_RETRIES = 5
CACHE_RETRY_DELAY = 0.5

def create_shared_gtf_cache(gtf_path: str) -> str:
    """
    Create a shared cache file for GTF annotations that can be used across processes.
    This eliminates the need to reload GTF for each worker.
    
    Args:
        gtf_path: Path to the GTF annotation file
        
    Returns:
        Path to the cache file
    """
    import hashlib
    
    # Create cache filename based on GTF file path and modification time
    gtf_stat = os.stat(gtf_path)
    cache_key = f"{gtf_path}_{gtf_stat.st_mtime}_{gtf_stat.st_size}"
    cache_hash = hashlib.md5(cache_key.encode()).hexdigest()
    
    cache_dir = os.path.join(os.path.dirname(__file__), "cache")
    os.makedirs(cache_dir, exist_ok=True)
    
    cache_file = os.path.join(cache_dir, f"gtf_cache_{cache_hash}.pkl")
    
    # Check if cache already exists and is newer than source
    if os.path.exists(cache_file):
        if os.path.getmtime(cache_file) > gtf_stat.st_mtime:
            print(f"[INFO] Using existing GTF cache: {cache_file}")
            return cache_file
    
    print(f"[INFO] Creating shared GTF cache from {gtf_path} (one-time operation)...")
    start_time = time.time()
    
    # Parse GTF file once
    chromosome_annotations = {}
    
    try:
        with open(gtf_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if not line or line.startswith('#'):
                    continue
                
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                
                seqid, source, feature_type, start_str, end_str, score, strand, phase, attributes = parts
                
                # Skip if not a gene or CDS
                if feature_type not in ('gene', 'CDS'):
                    continue
                
                try:
                    start = int(start_str)
                    end = int(end_str)
                except ValueError:
                    continue
                
                # Extract gene name from attributes
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
                
                # Initialize chromosome if not exists
                if seqid not in chromosome_annotations:
                    chromosome_annotations[seqid] = {'genes': [], 'cds': []}
                
                # Add annotation based on type
                if feature_type == 'gene':
                    gene = {
                        'name': gene_name,
                        'locus_tag': gene_name,
                        'start': min(start, end),
                        'end': max(start, end),
                        'strand': strand if strand in ['+', '-'] else '+',
                        'type': 'gene'
                    }
                    chromosome_annotations[seqid]['genes'].append(gene)
                elif feature_type == 'CDS':
                    cds = {
                        'name': gene_name,
                        'start': min(start, end),
                        'end': max(start, end),
                        'strand': strand if strand in ['+', '-'] else '+',
                        'type': 'cds'
                    }
                    chromosome_annotations[seqid]['cds'].append(cds)
    
    except FileNotFoundError:
        raise GenomeAnalysisError(f"GTF file not found at {gtf_path}")
    except Exception as e:
        raise GenomeAnalysisError(f"An error occurred loading GTF annotations: {e}")
    
    # Save to cache file
    try:
        with open(cache_file, 'wb') as f:
            pickle.dump(chromosome_annotations, f)
    except Exception as e:
        raise GenomeAnalysisError(f"Failed to save GTF cache: {e}")
    
    elapsed_time = time.time() - start_time
    total_genes = sum(len(data['genes']) for data in chromosome_annotations.values())
    total_cds = sum(len(data['cds']) for data in chromosome_annotations.values())
    
    print(f"[INFO] Successfully created shared cache: {total_genes} genes and {total_cds} CDS across {len(chromosome_annotations)} chromosomes in {elapsed_time:.2f}s")
    print(f"[INFO] Cache saved to: {cache_file}")
    
    return cache_file

def load_annotations_from_shared_cache(cache_file: str, chromosome: str) -> Tuple[List[Dict], List[Dict]]:
    """
    Load annotations for a specific chromosome from the shared cache file.
    
    Args:
        cache_file: Path to the shared cache file
        chromosome: Chromosome name to retrieve
        
    Returns:
        Tuple of (genes, cds) for the specified chromosome
    """
    try:
        with open(cache_file, 'rb') as f:
            chromosome_annotations = pickle.load(f)
        
        chromosome_data = chromosome_annotations.get(chromosome, {'genes': [], 'cds': []})
        return chromosome_data['genes'], chromosome_data['cds']
        
    except Exception as e:
        print(f"[WARNING] Failed to load from cache {cache_file}: {e}")
        return [], []

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

def load_complete_sequence(analyzer, filename: str):
    print(f"Loading complete genome sequence from {os.path.basename(filename)}...")
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            sequence_parts = [line.strip().upper() for line in f if not line.startswith('>')]
        analyzer.sequence = "".join(sequence_parts)
        print(f"Successfully loaded genome of {len(analyzer.sequence):,} bp.")
    except FileNotFoundError:
        raise GenomeAnalysisError(f"Genome file not found at {filename}")
    except Exception as e:
        raise GenomeAnalysisError(f"An error occurred loading the genome: {e}")

def load_chromosome_sequence(analyzer, fasta_path: str, chromosome: str):
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
        analyzer.sequence = "".join(parts)
        if not analyzer.sequence:
            raise GenomeAnalysisError(f"Chromosome '{chromosome}' not found in {os.path.basename(fasta_path)}")
        print(f"Successfully loaded chromosome '{chromosome}' of {len(analyzer.sequence):,} bp.")
    except FileNotFoundError:
        raise GenomeAnalysisError(f"FASTA file not found at {fasta_path}")
    except Exception as e:
        raise GenomeAnalysisError(f"An error occurred loading chromosome sequence: {e}")

def load_annotations(analyzer, filename: str):
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

        analyzer._build_gene_tree_with_intergenic(all_genes)
        print(f"Successfully built map for {len(all_genes)} genes and intergenic regions.")
    except Exception as e:
        _handle_file_error(e, filename, "loading annotations")

def load_annotations_from_bsub_json(analyzer, json_path: str):
    import os
    print(f"Loading Bacillus subtilis annotations from {os.path.basename(json_path)}...")
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Get genome length for validation
        genome_length = len(analyzer.sequence) if analyzer.sequence else 4215606  # Default B.sub genome length
        
        all_genes = []
        corrupted_genes = []
        suspicious_genes = []
        
        for gene_id, info in data.items():
            start = int(info.get('start', 0))
            end = int(info.get('end', 0))
            
            # Basic validation
            if not start or not end:
                print(f"[WARNING] Skipping {gene_id}: Invalid coordinates (start={start}, end={end})")
                continue
            
            # CRITICAL: Check for coordinate inversions (start > end)
            if start > end:
                print(f"[ERROR] Coordinate inversion detected in {gene_id}: start={start} > end={end}")
                print(f"[INFO] Auto-correcting coordinates for {gene_id}")
                start, end = end, start  # Swap coordinates
            
            # Validate coordinate ranges
            if start < 100:
                print(f"[WARNING] Suspicious early position for {gene_id}: start={start} (position < 100)")
                suspicious_genes.append(gene_id)
            
            if end > genome_length:
                print(f"[ERROR] Gene {gene_id} extends beyond genome: end={end} > genome_length={genome_length}")
                corrupted_genes.append(gene_id)
                continue
            
            # Validate gene length (bacterial genes typically < 50kb)
            gene_length = end - start + 1
            if gene_length > 50000:
                print(f"[WARNING] Suspiciously long gene {gene_id}: {gene_length:,} bp")
                suspicious_genes.append(gene_id)
            
            strand_value = info.get('strand', '+')
            strand = '+' if str(strand_value).lower().startswith('p') or strand_value == '+' else '-'
            
            gene = {
                'name': gene_id,
                'locus_tag': gene_id,
                'start': start,
                'end': end,
                'strand': strand,
                'type': 'gene'
            }
            all_genes.append(gene)
        
        # Report validation results
        if corrupted_genes:
            print(f"[ERROR] Found {len(corrupted_genes)} corrupted genes: {corrupted_genes[:5]}{'...' if len(corrupted_genes) > 5 else ''}")
        
        if suspicious_genes:
            print(f"[WARNING] Found {len(suspicious_genes)} suspicious genes: {suspicious_genes[:5]}{'...' if len(suspicious_genes) > 5 else ''}")
        
        analyzer._build_gene_tree_with_intergenic(all_genes)
        print(f"Successfully built map for {len(all_genes)} genes and intergenic regions (B. subtilis).")
        
        if corrupted_genes or suspicious_genes:
            print(f"[INFO] Data validation complete: {len(all_genes)} valid genes, {len(corrupted_genes)} corrupted, {len(suspicious_genes)} suspicious")
            
    except FileNotFoundError:
        raise GenomeAnalysisError(f"JSON annotation file not found at {json_path}")
    except Exception as e:
        raise GenomeAnalysisError(f"An error occurred loading B. subtilis annotations: {e}")

def load_annotations_from_gtf(analyzer, gtf_path: str, chromosome: str):
    print(f"Loading annotations from GTF for chromosome '{chromosome}'...")
    try:
        # Get genome length for validation
        genome_length = len(analyzer.sequence) if analyzer.sequence else 100000000  # Default large genome length
        
        all_genes = []
        all_cds = []
        corrupted_genes = []
        suspicious_genes = []
        
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
                
                # CRITICAL: Check for coordinate inversions (start > end)
                if start > end:
                    print(f"[ERROR] Coordinate inversion detected in GTF: start={start} > end={end}")
                    print(f"[INFO] Auto-correcting coordinates")
                    start, end = end, start  # Swap coordinates
                
                # Validate coordinate ranges
                if start < 1:
                    print(f"[WARNING] Suspicious position in GTF: start={start} (position < 1)")
                    suspicious_genes.append(f"gene_{start}_{end}")
                
                if end > genome_length:
                    print(f"[ERROR] Gene extends beyond genome: end={end} > genome_length={genome_length}")
                    corrupted_genes.append(f"gene_{start}_{end}")
                    continue
                
                # Validate gene length (human genes can be very long, set higher threshold)
                gene_length = end - start + 1
                if gene_length > 10000000:  # 10Mb threshold for human genes
                    print(f"[WARNING] Extremely long gene: {gene_length:,} bp")
                    suspicious_genes.append(f"gene_{start}_{end}")
                
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
                        'start': start,
                        'end': end,
                        'strand': strand if strand in ['+', '-'] else '+',
                        'type': 'gene'
                    }
                    all_genes.append(gene)
                elif feature_type == 'CDS':
                    cds = {
                        'name': gene_name,
                        'start': start,
                        'end': end,
                        'strand': strand if strand in ['+', '-'] else '+',
                        'type': 'cds'
                    }
                    all_cds.append(cds)
        
        # Report validation results
        if corrupted_genes:
            print(f"[ERROR] Found {len(corrupted_genes)} corrupted genes in GTF")
        
        if suspicious_genes:
            print(f"[WARNING] Found {len(suspicious_genes)} suspicious genes in GTF")
        
        analyzer._build_gene_tree_with_intergenic(all_genes)
        # Build CDS intervals
        for cds in all_cds:
            analyzer.cds_tree[cds['start']:cds['end'] + 1] = cds
        print(f"Successfully built map for {len(all_genes)} genes and intergenic regions (GTF).")
        
        if corrupted_genes or suspicious_genes:
            print(f"[INFO] GTF validation complete: {len(all_genes)} valid genes, {len(corrupted_genes)} corrupted, {len(suspicious_genes)} suspicious")
            
    except FileNotFoundError:
        raise GenomeAnalysisError(f"GTF file not found at {gtf_path}")
    except Exception as e:
        raise GenomeAnalysisError(f"An error occurred loading GTF annotations: {e}")

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
        import time
        
        cache_file = os.path.join(self.cache_dir, f"{chromosome}.raw")
        
        # Check if cache already exists and is newer than source
        if os.path.exists(cache_file):
            if os.path.getmtime(cache_file) > os.path.getmtime(fasta_path):
                return cache_file
        
        lock_file = cache_file + ".lock"
        if os.path.exists(lock_file):
            time.sleep(0.5)
            if os.path.exists(cache_file):
                return cache_file
        
        max_retries = CACHE_MAX_RETRIES
        retry_delay = CACHE_RETRY_DELAY
        
        for attempt in range(max_retries):
            try:
                with open(lock_file, 'x') as f:
                    f.write(str(os.getpid()))
                break
            except FileExistsError:
                if attempt < max_retries - 1:
                    if os.path.exists(cache_file):
                        return cache_file
                    
                    time.sleep(retry_delay)
                    continue
                else:
                    if os.path.exists(cache_file):
                        return cache_file
                    else:
                        try:
                            return self._create_cache_directly(fasta_path, chromosome, cache_file)
                        except Exception as e:
                            raise GenomeAnalysisError(f"Could not create cache file for {chromosome}: {e}")
        
        raw_sequence = self._extract_sequence_to_cache(fasta_path, chromosome, cache_file)
        
        with open(cache_file, 'w') as f:
            f.write(raw_sequence)
        
        try:
            if os.path.exists(lock_file):
                os.remove(lock_file)
        except Exception as e:
            print(f"[WARNING] Failed to remove lock file: {e}")
        
        return cache_file
    
    def _create_cache_directly(self, fasta_path: str, chromosome: str, cache_file: str) -> str:
        raw_sequence = self._extract_sequence_to_cache(fasta_path, chromosome, cache_file)
        
        with open(cache_file, 'w') as f:
            f.write(raw_sequence)
        
        return cache_file
    
    def _extract_sequence_to_cache(self, fasta_path: str, chromosome: str, cache_file: str) -> str:
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
        file_size = os.path.getsize(cache_file)
        
        if file_size > MEMORY_MAPPING_THRESHOLD_BYTES:
            with open(cache_file, 'rb') as f:
                with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                    return mm.read().decode('utf-8')
        else:
            with open(cache_file, 'r') as f:
                return f.read()
    
    def cleanup_cache(self):
        try:
            cache_files = []
            total_size = 0
            
            for filename in os.listdir(self.cache_dir):
                if filename.endswith('.raw'):
                    filepath = os.path.join(self.cache_dir, filename)
                    size = os.path.getsize(filepath)
                    cache_files.append((filepath, size, os.path.getmtime(filepath)))
                    total_size += size
            
            cache_files.sort(key=lambda x: x[2])
            
            max_size_bytes = self.max_cache_size_mb * 1024 * 1024
            for filepath, size, _ in cache_files:
                if total_size <= max_size_bytes:
                    break
                os.remove(filepath)
                total_size -= size
        except Exception as e:
            print(f"[WARNING] Cache cleanup failed: {e}")

def load_genome_for_chromosome(fasta_path: str, chromosome: str, use_memory_mapping: bool = True) -> str:
    file_size = os.path.getsize(fasta_path)
    
    # OPTIMIZATION: Special case for loading complete genome (ALL chromosomes)
    if chromosome == "ALL":
        print(f"[INFO] Loading complete genome sequence ({file_size:,} bytes)...")
        analyzer = GenomeAnalyzer()
        analyzer.load_complete_sequence(fasta_path)
        print(f"[INFO] Complete genome loaded successfully ({len(analyzer.sequence):,} bp)")
        return analyzer.sequence
    
    if use_memory_mapping and file_size > MEMORY_MAPPING_THRESHOLD_BYTES:
        try:
            loader = MemoryMappedGenomeLoader()
            cache_file = loader.create_raw_sequence_cache(fasta_path, chromosome)
            sequence = loader.load_sequence_memory_mapped(cache_file)
            return sequence
        except Exception as e:
            raise GenomeAnalysisError(f"Memory mapping failed: {e}. Check file permissions and disk space.")
    else:
        analyzer = GenomeAnalyzer()
        analyzer.load_chromosome_sequence(fasta_path, chromosome)
        return analyzer.sequence

def _handle_file_error(error: Exception, file_path: str, operation: str) -> None:
    """
    Centralized file error handling to eliminate code duplication.
    Raises appropriate GenomeAnalysisError with consistent messaging.
    """
    if isinstance(error, FileNotFoundError):
        raise GenomeAnalysisError(f"{operation} file not found at {file_path}")
    else:
        raise GenomeAnalysisError(f"An error occurred during {operation}: {error}")
