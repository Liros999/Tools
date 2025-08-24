"""
File processing for genome analysis.
Refactored to use modular IO package for better maintainability.
"""

# Import all file processing functions from the IO package
from .io import (
    load_genome_for_chromosome,
    load_annotations_from_gtf,
    create_shared_gtf_cache,
    load_annotations_from_shared_cache,
    parse_fasta_file,
    parse_gtf_file,
    extract_chromosome_names,
    validate_sequence
)

# Re-export for backward compatibility
__all__ = [
    'load_genome_for_chromosome',
    'load_annotations_from_gtf',
    'create_shared_gtf_cache',
    'load_annotations_from_shared_cache',
    'parse_fasta_file',
    'parse_gtf_file',
    'extract_chromosome_names',
    'validate_sequence'
]
