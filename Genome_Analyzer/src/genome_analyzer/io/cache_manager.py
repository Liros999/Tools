"""
Cache manager for genome analysis.
Handles shared GTF cache creation and loading for efficient annotation access.
"""

import os
import pickle
import hashlib
from typing import Dict, List, Tuple, Optional


def create_shared_gtf_cache(gtf_path: str) -> str:
    """
    Create a shared GTF cache file for efficient annotation loading.
    
    Args:
        gtf_path: Path to GTF annotation file
    
    Returns:
        Path to the created cache file
    """
    try:
        if not os.path.exists(gtf_path):
            raise FileNotFoundError(f"GTF file not found: {gtf_path}")
        
        # Create cache directory if it doesn't exist
        cache_dir = os.path.join(os.path.dirname(gtf_path), "cache")
        os.makedirs(cache_dir, exist_ok=True)
        
        # Generate cache filename based on file properties
        gtf_stat = os.stat(gtf_path)
        cache_key = f"{gtf_path}_{gtf_stat.st_mtime}_{gtf_stat.st_size}"
        cache_hash = hashlib.md5(cache_key.encode()).hexdigest()
        cache_file = os.path.join(cache_dir, f"gtf_cache_{cache_hash}.pkl")
        
        # Check if cache already exists and is valid
        if os.path.exists(cache_file):
            print(f"[INFO] Using existing GTF cache: {os.path.basename(cache_file)}")
            return cache_file
        
        print(f"[INFO] Creating new GTF cache: {os.path.basename(cache_file)}")
        
        # Parse GTF file and create cache
        cache_data = _parse_gtf_for_cache(gtf_path)
        
        # Save cache to file
        with open(cache_file, 'wb') as f:
            pickle.dump(cache_data, f)
        
        print(f"[INFO] GTF cache created successfully: {len(cache_data)} chromosomes")
        return cache_file
        
    except Exception as e:
        print(f"[ERROR] Failed to create GTF cache: {e}")
        raise


def _parse_gtf_for_cache(gtf_path: str) -> Dict[str, Dict]:
    """Parse GTF file and create cache data structure."""
    try:
        cache_data = {}
        current_chromosome = None
        genes = {}
        cds = {}
        
        with open(gtf_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                
                seqid, source, feature_type, start_str, end_str, score, strand, frame, attributes = parts
                
                # Check if we've moved to a new chromosome
                if seqid != current_chromosome:
                    # Save previous chromosome data
                    if current_chromosome and (genes or cds):
                        cache_data[current_chromosome] = {
                            'genes': genes.copy(),
                            'cds': cds.copy()
                        }
                    
                    # Reset for new chromosome
                    current_chromosome = seqid
                    genes = {}
                    cds = {}
                
                try:
                    start = int(start_str)
                    end = int(end_str)
                except ValueError:
                    continue
                
                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key.strip()] = value.strip().strip('"')
                
                if feature_type == 'gene':
                    gene_id = attr_dict.get('gene_id', f'gene_{len(genes)}')
                    genes[gene_id] = {
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'source': source,
                        'chromosome': seqid
                    }
                    
                elif feature_type == 'CDS':
                    cds_id = attr_dict.get('transcript_id', f'cds_{len(cds)}')
                    cds[cds_id] = {
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'gene_id': attr_dict.get('gene_id', ''),
                        'chromosome': seqid
                    }
        
        # Save last chromosome data
        if current_chromosome and (genes or cds):
            cache_data[current_chromosome] = {
                'genes': genes,
                'cds': cds
            }
        
        return cache_data
        
    except Exception as e:
        print(f"[ERROR] Failed to parse GTF for cache: {e}")
        raise


def load_annotations_from_shared_cache(cache_file: str, chromosome: str) -> Tuple[Dict, Dict]:
    """
    Load annotations for a specific chromosome from shared cache.
    
    Args:
        cache_file: Path to cache file
        chromosome: Chromosome to load annotations for
    
    Returns:
        Tuple of (genes_dict, cds_dict)
    """
    try:
        if not os.path.exists(cache_file):
            raise FileNotFoundError(f"Cache file not found: {cache_file}")
        
        # Load cache data
        with open(cache_file, 'rb') as f:
            cache_data = pickle.load(f)
        
        if chromosome not in cache_data:
            print(f"[WARNING] Chromosome {chromosome} not found in cache")
            return {}, {}
        
        chromosome_data = cache_data[chromosome]
        genes = chromosome_data.get('genes', {})
        cds = chromosome_data.get('cds', {})
        
        print(f"[INFO] Loaded {len(genes)} genes and {len(cds)} CDS for {chromosome} from cache")
        return genes, cds
        
    except Exception as e:
        print(f"[ERROR] Failed to load annotations from cache: {e}")
        return {}, {}


# Export functions
__all__ = [
    'create_shared_gtf_cache',
    'load_annotations_from_shared_cache'
]
