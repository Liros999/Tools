import requests
import csv
import os
import json
from typing import Dict, List, Tuple, Optional

class GeneData:
    def __init__(self):
        self.genes: Dict[str, Dict] = {}  # gene_name -> {start, end, strand, sequence}
        self.coordinates: List[Tuple[int, int, str, str]] = []  # (start, end, strand, gene_name)
        
    def fetch_all_genes(self) -> bool:
        """Fetch all genes from SubtiWiki and cache them"""
        BASE_URL = "https://subtiwiki.uni-goettingen.de/v5/api"
        try:
            # First get list of all genes
            response = requests.get(f"{BASE_URL}/genes")
            response.raise_for_status()
            genes_list = response.json().get("data", [])
            
            for gene in genes_list:
                gene_name = gene.get("name")
                if not gene_name:
                    continue
                    
                # Get detailed gene data
                gene_data = self._fetch_gene_data(gene_name)
                if gene_data:
                    self.genes[gene_name] = gene_data
                    # Add to coordinates list for easy searching
                    self.coordinates.append((
                        gene_data["start"],
                        gene_data["end"],
                        gene_data["strand"],
                        gene_name
                    ))
            
            # Sort coordinates by start position for efficient searching
            self.coordinates.sort(key=lambda x: x[0])
            return True
            
        except Exception as e:
            print(f"Error fetching genes: {e}")
            return False
    
    def _fetch_gene_data(self, gene_name: str) -> Dict:
        """Fetch detailed data for a single gene"""
        BASE_URL = "https://subtiwiki.uni-goettingen.de/v5/api"
        try:
            response = requests.get(f"{BASE_URL}/gene/{gene_name}/all")
            response.raise_for_status()
            data = response.json()
            
            # Extract relevant information
            annotations = data.get("data", {}).get("genomic_annotations", [])
            if not annotations:
                return None
                
            annotation = annotations[0]  # Take first annotation
            return {
                "start": annotation.get("start"),
                "end": annotation.get("end"),
                "strand": annotation.get("orientation"),
                "sequence": annotation.get("dna_sequence", "")
            }
        except Exception as e:
            print(f"Error fetching data for {gene_name}: {e}")
            return None
    
    def save_to_csv(self, filename: str = "gene_data.csv"):
        """Save gene data to CSV file"""
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Gene", "Start", "End", "Strand", "Sequence"])
            for gene_name, data in self.genes.items():
                writer.writerow([
                    gene_name,
                    data["start"],
                    data["end"],
                    data["strand"],
                    data["sequence"]
                ])
    
    def load_from_csv(self, filename: str = "gene_data.csv"):
        """Load gene data from CSV file"""
        if not os.path.exists(filename):
            return False
            
        with open(filename, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                gene_name = row["Gene"]
                self.genes[gene_name] = {
                    "start": int(row["Start"]),
                    "end": int(row["End"]),
                    "strand": row["Strand"],
                    "sequence": row["Sequence"]
                }
                self.coordinates.append((
                    int(row["Start"]),
                    int(row["End"]),
                    row["Strand"],
                    gene_name
                ))
        self.coordinates.sort(key=lambda x: x[0])
        return True
    
    def load_from_json(self, filename: str = "all_gene_data.json") -> bool:
        """Load gene data from JSON file"""
        try:
            if not os.path.exists(filename):
                return False
                
            with open(filename, 'r') as f:
                self.genes = json.load(f)
            return True
        except Exception as e:
            print(f"Error loading gene data: {e}")
            return False
    
    def get_gene_info(self, gene_name: str) -> Optional[Dict]:
        """Get information for a specific gene"""
        return self.genes.get(gene_name)
    
    def find_genes_for_coordinates(self, start: int, end: int) -> List[str]:
        """Find genes that overlap with given coordinates"""
        result = []
        for gene_name, gene_info in self.genes.items():
            gene_start = gene_info.get('start')
            gene_end = gene_info.get('end')
            if gene_start is not None and gene_end is not None:
                # Check if coordinates overlap with gene
                if (start <= gene_end and end >= gene_start):
                    result.append(gene_name)
        return result
    
    def find_nearest_genes(self, start: int, end: int) -> Tuple[Optional[str], Optional[str]]:
        """Find nearest genes on both sides of the given coordinates"""
        left_gene = None
        right_gene = None
        min_left_dist = float('inf')
        min_right_dist = float('inf')
        
        for gene_name, gene_info in self.genes.items():
            gene_start = gene_info.get('start')
            gene_end = gene_info.get('end')
            if gene_start is None or gene_end is None:
                continue
                
            # Gene is to the left
            if gene_end < start:
                dist = start - gene_end
                if dist < min_left_dist:
                    min_left_dist = dist
                    left_gene = gene_name
            # Gene is to the right
            elif gene_start > end:
                dist = gene_start - end
                if dist < min_right_dist:
                    min_right_dist = dist
                    right_gene = gene_name
                    
        return (left_gene, right_gene) 