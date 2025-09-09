import requests
import csv
import os
import json
from typing import Dict, List, Tuple, Optional
import bisect

class GeneData:
    def __init__(self):
        self.genes: Dict[str, Dict] = {}
        self.coordinates: List[Tuple[int, int, str, str]] = []

    def load_from_gff(self, filename: str) -> bool:
        """Load gene data from a GFF file (E. coli format)"""
        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            return False
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue
                    seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
                    if feature_type != 'gene':
                        continue
                    # Extract gene name from attributes
                    gene_name = None
                    for attr in attributes.split(';'):
                        if attr.startswith('Name='):
                            gene_name = attr.split('=', 1)[1]
                            break
                    if not gene_name:
                        continue
                    start = int(start)
                    end = int(end)
                    self.genes[gene_name] = {
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'sequence': ''  # Sequence can be filled if needed
                    }
                    self.coordinates.append((start, end, strand, gene_name))
            self.coordinates.sort(key=lambda x: x[0])
            return True
        except Exception as e:
            print(f"Error loading GFF gene data: {e}")
            return False

    def load_from_fasta_headers(self, filename: str) -> bool:
        """Load gene data from FASTA file headers (E. coli K12 format)"""
        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            return False
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # Example header: >lcl|NC_000913.3_cds_NP_414542.1_1 [gene=thrL] [locus_tag=b0001] ... [location=190..255]
                        gene_name_match = re.search(r'\[gene=([^\]]+)\]', line)
                        location_match = re.search(r'\[location=(\d+)\.\.(\d+)\]', line)
                        strand_match = re.search(r'\[location=complement\((\d+)\.\.(\d+)\)\]', line)

                        gene_name = gene_name_match.group(1) if gene_name_match else None
                        
                        start = None
                        end = None
                        strand = '+'

                        if location_match:
                            start = int(location_match.group(1))
                            end = int(location_match.group(2))
                        elif strand_match:
                            start = int(strand_match.group(1))
                            end = int(strand_match.group(2))
                            strand = '-'

                        if gene_name and start is not None and end is not None:
                            self.genes[gene_name] = {
                                'start': start,
                                'end': end,
                                'strand': strand,
                                'sequence': ''  # Sequence can be filled if needed
                            }
                            self.coordinates.append((start, end, strand, gene_name))
            self.coordinates.sort(key=lambda x: x[0])
            return True
        except Exception as e:
            print(f"Error loading gene data from FASTA headers: {e}")
            return False

    def get_gene_info(self, gene_name: str) -> Optional[Dict]:
        """Get information for a specific gene"""
        return self.genes.get(gene_name)

    def find_genes_for_coordinates(self, start: int, end: int) -> List[str]:
        """Find genes that overlap with given coordinates using binary search."""
        result = []
        # Find the insertion point for the start coordinate
        idx = bisect.bisect_left(self.coordinates, (start, 0, '', ''))

        # Check genes at and before the insertion point
        for i in range(max(0, idx - 1), len(self.coordinates)):
            gene_start, gene_end, _, gene_name = self.coordinates[i]
            # If the gene starts after the end of our search range, we can stop
            if gene_start > end:
                break
            # Check for overlap
            if start <= gene_end and end >= gene_start:
                result.append(gene_name)
        return result