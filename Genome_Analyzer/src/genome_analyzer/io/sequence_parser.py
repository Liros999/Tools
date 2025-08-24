"""
Sequence parser for genome analysis.
Handles parsing of FASTA and GTF files for sequence and annotation data.
"""

import os
import re
from typing import Dict, List, Tuple, Optional, Generator


def parse_fasta_file(fasta_path: str, chromosome_filter: Optional[str] = None) -> Generator[Tuple[str, str], None, None]:
    """
    Parse FASTA file and yield (header, sequence) pairs.
    
    Args:
        fasta_path: Path to FASTA file
        chromosome_filter: Optional chromosome name to filter by
    
    Yields:
        Tuple of (header, sequence)
    """
    try:
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
        
        current_header = None
        current_sequence = []
        
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_header and current_sequence:
                        if not chromosome_filter or chromosome_filter in current_header:
                            yield current_header, ''.join(current_sequence)
                    
                    # Start new sequence
                    current_header = line[1:]  # Remove '>'
                    current_sequence = []
                else:
                    # Add sequence line
                    current_sequence.append(line)
        
        # Yield last sequence
        if current_header and current_sequence:
            if not chromosome_filter or chromosome_filter in current_header:
                yield current_header, ''.join(current_sequence)
                
    except Exception as e:
        print(f"[ERROR] Failed to parse FASTA file: {e}")
        raise


def parse_gtf_file(gtf_path: str, chromosome_filter: Optional[str] = None) -> Generator[Dict, None, None]:
    """
    Parse GTF file and yield annotation records.
    
    Args:
        gtf_path: Path to GTF file
        chromosome_filter: Optional chromosome name to filter by
    
    Yields:
        Dictionary containing annotation data
    """
    try:
        if not os.path.exists(gtf_path):
            raise FileNotFoundError(f"GTF file not found: {gtf_path}")
        
        with open(gtf_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                
                seqid, source, feature_type, start_str, end_str, score, strand, frame, attributes = parts
                
                # Filter by chromosome if specified
                if chromosome_filter and seqid != chromosome_filter:
                    continue
                
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
                
                # Create annotation record
                record = {
                    'seqid': seqid,
                    'source': source,
                    'feature_type': feature_type,
                    'start': start,
                    'end': end,
                    'score': score,
                    'strand': strand,
                    'frame': frame,
                    'attributes': attr_dict,
                    'line_number': line_num
                }
                
                yield record
                
    except Exception as e:
        print(f"[ERROR] Failed to parse GTF file: {e}")
        raise


def extract_chromosome_names(gtf_path: str) -> List[str]:
    """
    Extract unique chromosome names from GTF file.
    
    Args:
        gtf_path: Path to GTF file
    
    Returns:
        List of unique chromosome names
    """
    try:
        chromosomes = set()
        
        for record in parse_gtf_file(gtf_path):
            chromosomes.add(record['seqid'])
        
        # Sort chromosomes for consistent ordering
        sorted_chromosomes = sorted(list(chromosomes))
        
        print(f"[INFO] Found {len(sorted_chromosomes)} chromosomes in GTF file")
        return sorted_chromosomes
        
    except Exception as e:
        print(f"[ERROR] Failed to extract chromosome names: {e}")
        return []


def validate_sequence(sequence: str) -> bool:
    """
    Validate DNA sequence for valid characters.
    
    Args:
        sequence: DNA sequence string
    
    Returns:
        True if sequence is valid, False otherwise
    """
    if not sequence:
        return False
    
    # Valid DNA characters (including IUPAC codes)
    valid_chars = set('ACGTNacgtnRYSMKWBDHVrysmkwbdhv')
    
    return all(char in valid_chars for char in sequence)


# Export functions
__all__ = [
    'parse_fasta_file',
    'parse_gtf_file',
    'extract_chromosome_names',
    'validate_sequence'
]
