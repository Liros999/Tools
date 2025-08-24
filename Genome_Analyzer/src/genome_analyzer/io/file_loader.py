"""
File loader for genome analysis.
Handles loading of genome sequences and annotation files.
"""

import os
import pickle
import hashlib
from typing import Dict, List, Tuple, Optional
from pathlib import Path


def load_genome_for_chromosome(fasta_path: str, chromosome: str = "ALL", use_memory_mapping: bool = True) -> str:
    """
    Load genome sequence for a specific chromosome or entire genome.
    
    Args:
        fasta_path: Path to FASTA file
        chromosome: Chromosome name or "ALL" for complete genome
        use_memory_mapping: Whether to use memory mapping for large files
    
    Returns:
        Genome sequence as string
    """
    try:
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
        
        print(f"[INFO] Loading genome from: {fasta_path}")
        
        if use_memory_mapping and chromosome == "ALL":
            # Use memory mapping for large complete genome files
            return _load_genome_with_memory_mapping(fasta_path)
        else:
            # Load specific chromosome or use standard loading
            return _load_genome_standard(fasta_path, chromosome)
            
    except Exception as e:
        print(f"[ERROR] Failed to load genome: {e}")
        raise


def _load_genome_with_memory_mapping(fasta_path: str) -> str:
    """Load complete genome using memory mapping for efficiency."""
    try:
        import mmap
        
        with open(fasta_path, 'r') as f:
            # Memory map the file
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                # Convert to string, skipping header lines
                content = mm.read().decode('utf-8')
                
                # Remove FASTA headers and join sequences
                lines = content.split('\n')
                sequence_lines = [line for line in lines if not line.startswith('>')]
                genome_sequence = ''.join(sequence_lines)
                
                print(f"[INFO] Loaded complete genome: {len(genome_sequence):,} bp")
                return genome_sequence
                
    except Exception as e:
        print(f"[WARNING] Memory mapping failed, falling back to standard loading: {e}")
        return _load_genome_standard(fasta_path, "ALL")


def _load_genome_standard(fasta_path: str, chromosome: str) -> str:
    """Load genome using standard file reading."""
    try:
        genome_sequence = ""
        current_chrom = None
        
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Header line - extract chromosome name
                    current_chrom = line.split()[0][1:]  # Remove '>' and get first word
                    if chromosome != "ALL" and current_chrom != chromosome:
                        current_chrom = None
                elif current_chrom and (chromosome == "ALL" or current_chrom == chromosome):
                    # Sequence line for target chromosome
                    genome_sequence += line
        
        if not genome_sequence:
            raise ValueError(f"No sequence found for chromosome: {chromosome}")
        
        print(f"[INFO] Loaded {'complete genome' if chromosome == 'ALL' else chromosome}: {len(genome_sequence):,} bp")
        return genome_sequence
        
    except Exception as e:
        print(f"[ERROR] Standard genome loading failed: {e}")
        raise


def load_annotations_from_gtf(analyzer, gtf_path: str, chromosome: str = "chr1") -> Tuple[Dict, Dict]:
    """
    Load gene annotations from GTF file for a specific chromosome.
    
    Args:
        analyzer: GenomeAnalyzer instance
        gtf_path: Path to GTF annotation file
        chromosome: Chromosome to load annotations for
    
    Returns:
        Tuple of (genes_dict, cds_dict)
    """
    try:
        if not os.path.exists(gtf_path):
            raise FileNotFoundError(f"GTF file not found: {gtf_path}")
        
        print(f"[INFO] Loading GTF annotations for {chromosome} from: {gtf_path}")
        
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
                
                # Filter by chromosome
                if seqid != chromosome:
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
        
        # Validate gene lengths for human genome
        for gene_id, gene_info in genes.items():
            gene_length = gene_info['end'] - gene_info['start']
            if gene_length > 10_000_000:  # 10Mb threshold for human genes
                print(f"[WARNING] Suspiciously long gene {gene_id}: {gene_length:,} bp")
        
        print(f"[INFO] Loaded {len(genes)} genes and {len(cds)} CDS for {chromosome}")
        return genes, cds
        
    except Exception as e:
        print(f"[ERROR] Failed to load GTF annotations: {e}")
        raise


# Export functions
__all__ = [
    'load_genome_for_chromosome',
    'load_annotations_from_gtf'
]
