"""
IO package for genome analysis.
Contains file processing, data loading, and caching functionality.
"""

from .file_loader import load_genome_for_chromosome, load_annotations_from_gtf
from .cache_manager import create_shared_gtf_cache, load_annotations_from_shared_cache
from .sequence_parser import parse_fasta_file, parse_gtf_file

__all__ = [
    'load_genome_for_chromosome',
    'load_annotations_from_gtf',
    'create_shared_gtf_cache',
    'load_annotations_from_shared_cache',
    'parse_fasta_file',
    'parse_gtf_file'
]
