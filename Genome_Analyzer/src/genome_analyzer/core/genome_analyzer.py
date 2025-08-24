"""
Main GenomeAnalyzer class for genome analysis operations.
Handles sequence analysis, gene location detection, and search operations.
"""

import os
import re
from typing import List, Tuple, Optional, Dict, Any
from intervaltree import IntervalTree
from tqdm import tqdm

from .exceptions import GenomeAnalysisError
from .resource_manager import ResourceManager
from .. import search


class GenomeAnalyzer:
    """A comprehensive tool to analyze genomes with species-specific logic."""
    
    def __init__(self):
        self.sequence: str = ""
        self.current_chromosome: str = ""
        self.gene_data: Dict[str, Dict] = {}
        self.gene_tree = IntervalTree()
        self.cds_tree = IntervalTree()
        self.iupac_map = {
            'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]', 'K': '[GT]',
            'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
            'N': '[ACGT]'
        }
        self.species_type = "generic"  # Homo Sapiens uses generic logic
        self.verbose_output = False
        self.multi_chromosome_mode = False
        self.available_chromosomes = []
        self.genome_file_path = ""
        self.annotation_file_path = ""
        
        # Initialize resource manager
        self.resource_manager = ResourceManager()

    @property
    def iupac_wildcards(self) -> Dict[str, str]:
        """Get IUPAC wildcards on-demand to eliminate redundant storage."""
        return {k: v.strip('[]') for k, v in self.iupac_map.items()}

    def _reverse_complement(self, seq: str) -> str:
        """Generate reverse complement of DNA sequence."""
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                          'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
        return "".join(complement_map.get(base, base) for base in reversed(seq))

    def get_location_details(self, start: int, end: int) -> str:
        """Get precise location details for a match in Homo Sapiens genome."""
        return self._get_location_details_generic(start, end)

    def _get_location_details_generic(self, start: int, end: int) -> str:
        """Generic location detection for Homo Sapiens genome."""
        try:
            # Find genes that overlap with this position
            overlapping_genes = []
            for gene_id, gene_info in self.gene_data.items():
                if (start <= gene_info['end'] and end >= gene_info['start']):
                    overlapping_genes.append((gene_id, gene_info))
            
            if overlapping_genes:
                # Match overlaps with gene(s)
                if len(overlapping_genes) == 1:
                    gene_id, gene_info = overlapping_genes[0]
                    return f"[{gene_id}][Overlap]"
                else:
                    # Multiple gene overlap
                    gene_names = [g[0] for g in overlapping_genes]
                    return f"[MultiGene]({','.join(gene_names)})"
            else:
                # Match is in intergenic region - find flanking genes
                return self._get_intergenic_location_generic(start, end)
                
        except Exception as e:
            return f"[Generic_Error:{str(e)[:20]}]"

    def _get_intergenic_location_generic(self, start: int, end: int) -> str:
        """Generic intergenic location detection for Homo Sapiens."""
        try:
            # Find nearest genes upstream and downstream
            upstream_genes = []
            downstream_genes = []
            
            for gene_id, gene_info in self.gene_data.items():
                if gene_info['end'] < start:  # Upstream
                    distance = start - gene_info['end']
                    upstream_genes.append((gene_id, distance))
                elif gene_info['start'] > end:  # Downstream
                    distance = gene_info['start'] - end
                    downstream_genes.append((gene_id, distance))
            
            # Sort by distance
            upstream_genes.sort(key=lambda x: x[1])
            downstream_genes.sort(key=lambda x: x[1])
            
            # Build location string
            location_parts = []
            
            if upstream_genes:
                nearest_upstream = upstream_genes[0]
                location_parts.append(f"<{nearest_upstream[1]}bp[{nearest_upstream[0]}]")
            
            if downstream_genes:
                nearest_downstream = downstream_genes[0]
                location_parts.append(f"[{nearest_downstream[0]}]<{nearest_downstream[1]}bp")
            
            if location_parts:
                return f"[Intergenic]({', '.join(location_parts)})"
            else:
                return "[Intergenic][Unknown]"
                
        except Exception as e:
            return f"[Intergenic_Error:{str(e)[:20]}]"

    def _build_gene_tree_with_intergenic(self, genes: List[Dict]):
        """Build gene tree with intergenic regions for location detection."""
        try:
            self.gene_tree = IntervalTree()
            
            for gene in genes:
                if 'start' in gene and 'end' in gene:
                    start = gene['start']
                    end = gene['end']
                    self.gene_tree[start:end + 1] = gene
                    
        except Exception as e:
            print(f"[WARNING] Failed to build gene tree: {e}")

    def search_sequence(self, pattern: str, max_mismatches: int = 0, boundary_bp: int = 0, stream: bool = False) -> List[Tuple]:
        """Search for pattern in genome sequence with specified parameters."""
        try:
            if not self.sequence:
                raise GenomeAnalysisError("No sequence loaded for analysis")
            
            # Use search module for pattern matching
            matches = search.search_pattern(self.sequence, pattern, max_mismatches)
            
            if not matches:
                return []
            
            # Process results
            results = []
            for start, end, window, strand, mismatches in matches:
                location = self.get_location_details(start, end)
                
                # Apply boundary annotations if specified
                if boundary_bp and boundary_bp > 0:
                    seq_len = len(self.sequence)
                    dist_start = start - 1
                    dist_end = seq_len - end
                    near_flags = []
                    if dist_start <= boundary_bp:
                        near_flags.append(f"NearStart:{dist_start}bp")
                    if dist_end <= boundary_bp:
                        near_flags.append(f"NearEnd:{dist_end}bp")
                    if near_flags:
                        location = f"{location} [{' & '.join(near_flags)}]"
                
                # Create result tuple
                result = (start, end, strand, pattern, window, self._reverse_complement(window), mismatches, location)
                results.append(result)
                
                # Stream output if requested
                if stream:
                    from ..ui import _print_row_tuple
                    _print_row_tuple(result, include_chrom=False)
            
            return results
            
        except Exception as e:
            raise GenomeAnalysisError(f"Search failed: {str(e)}")

    def load_sequence(self, sequence: str, chromosome: str = ""):
        """Load sequence data into analyzer."""
        self.sequence = sequence
        self.current_chromosome = chromosome
        print(f"[INFO] Loaded sequence: {len(sequence):,} bp")

    def load_gene_data(self, gene_data: Dict[str, Dict]):
        """Load gene annotation data."""
        self.gene_data = gene_data
        print(f"[INFO] Loaded gene data: {len(gene_data)} genes")

    def set_verbose_output(self, verbose: bool):
        """Set verbose output mode."""
        self.verbose_output = verbose

    def get_sequence_stats(self) -> Dict[str, Any]:
        """Get sequence statistics."""
        if not self.sequence:
            return {}
        
        stats = {
            'length': len(self.sequence),
            'gc_content': (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100,
            'chromosome': self.current_chromosome,
            'gene_count': len(self.gene_data),
            'species_type': self.species_type
        }
        return stats


# Export class
__all__ = ['GenomeAnalyzer']
