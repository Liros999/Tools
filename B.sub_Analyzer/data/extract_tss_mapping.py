#!/usr/bin/env python3
"""
Extract TSS coordinates and transcription factor mapping from BSGatlas GFF file.
Creates a clean mapping of TSS coordinates to sigma factors and transcription factors.
"""

import re
import json
import csv

def parse_tss_data():
    """Parse TSS data from the GFF file and extract coordinates and sigma factors."""
    tss_mapping = []
    
    # Read directly from the original GFF file
    with open('TSS_BSGatlas_v1.1.gff', 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            # Parse GFF fields: seqname, source, feature, start, end, score, strand, frame, attributes
            fields = line.split('\t')
            if len(fields) < 9:
                continue
                
            # Only process TSS entries with sigma factor information
            if fields[2] != 'TSS' or 'Sigma:' not in fields[8]:
                continue
                
            seqname = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            # Extract TSS ID
            tss_id_match = re.search(r'ID=([^;]+)', attributes)
            tss_id = tss_id_match.group(1) if tss_id_match else "Unknown"
            
            # Extract sigma factor from comment
            sigma_match = re.search(r'Sigma:\s*([^,;]+)', attributes)
            sigma_factor = sigma_match.group(1).strip() if sigma_match else "Unknown"
            
            # Extract resolution
            resolution_match = re.search(r'Resolution:\s*(\d+)', attributes)
            resolution = int(resolution_match.group(1)) if resolution_match else None
            
            # Extract PubMed references
            pubmed_match = re.search(r'PubMed:\s*([^"]+)', attributes)
            pubmed_refs = pubmed_match.group(1).strip() if pubmed_match else ""
            
            # Extract source information
            source_match = re.search(r'Based on:\s*([^,]+)', attributes)
            data_source = source_match.group(1).strip() if source_match else ""
            
            tss_entry = {
                'tss_id': tss_id,
                'chromosome': seqname,
                'coordinate': start,  # TSS coordinate (start == end for TSS)
                'strand': strand,
                'sigma_factor': sigma_factor,
                'resolution': resolution,
                'data_source': data_source,
                'pubmed_references': pubmed_refs
            }
            
            tss_mapping.append(tss_entry)
    
    return tss_mapping

def analyze_sigma_factors(tss_data):
    """Analyze sigma factor distribution and clean up the data."""
    sigma_counts = {}
    cleaned_data = []
    
    for entry in tss_data:
        sigma = entry['sigma_factor']
        
        # Clean and standardize sigma factor names
        if sigma == "?":
            sigma = "Unknown"
        elif ";" in sigma:
            # Handle multiple sigma factors (e.g., "E;F")
            sigma = sigma.replace(";", "/")
        
        entry['sigma_factor'] = sigma
        cleaned_data.append(entry)
        
        # Count sigma factors
        if sigma in sigma_counts:
            sigma_counts[sigma] += 1
        else:
            sigma_counts[sigma] = 1
    
    return cleaned_data, sigma_counts

def create_coordinate_mapping(tss_data):
    """Create a coordinate-based mapping for quick lookup."""
    coordinate_mapping = {}
    
    for entry in tss_data:
        coord = entry['coordinate']
        if coord not in coordinate_mapping:
            coordinate_mapping[coord] = []
        coordinate_mapping[coord].append(entry)
    
    return coordinate_mapping

def main():
    print("Extracting TSS data from BSGatlas GFF file...")
    
    # Parse TSS data
    tss_data = parse_tss_data()
    print(f"Found {len(tss_data)} TSS entries")
    
    # Clean and analyze sigma factors
    cleaned_data, sigma_counts = analyze_sigma_factors(tss_data)
    
    # Create coordinate mapping
    coordinate_mapping = create_coordinate_mapping(cleaned_data)
    
    # Save results
    print("\nSaving results...")
    
    # Save complete TSS data as JSON
    with open('tss_complete_mapping.json', 'w') as f:
        json.dump(cleaned_data, f, indent=2)
    
    # Save coordinate mapping as JSON
    with open('tss_coordinate_mapping.json', 'w') as f:
        json.dump(coordinate_mapping, f, indent=2)
    
    # Save as CSV for easy viewing
    with open('tss_mapping.csv', 'w', newline='', encoding='utf-8') as f:
        if cleaned_data:
            fieldnames = cleaned_data[0].keys()
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(cleaned_data)
    
    # Save sigma factor analysis
    with open('sigma_factor_analysis.json', 'w') as f:
        json.dump(sigma_counts, f, indent=2)
    
    # Print summary
    print(f"\nSUMMARY:")
    print(f"Total TSS entries: {len(cleaned_data)}")
    print(f"Unique coordinates: {len(coordinate_mapping)}")
    print(f"Sigma factors found: {len(sigma_counts)}")
    
    print(f"\nSigma factor distribution:")
    for sigma, count in sorted(sigma_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {sigma}: {count} TSS")
    
    print(f"\nFiles created:")
    print(f"  - tss_complete_mapping.json (complete data)")
    print(f"  - tss_coordinate_mapping.json (coordinate-indexed)")
    print(f"  - tss_mapping.csv (spreadsheet format)")
    print(f"  - sigma_factor_analysis.json (sigma factor counts)")
    
    # Show a few examples
    print(f"\nExample TSS entries:")
    for i, entry in enumerate(cleaned_data[:5]):
        print(f"  {i+1}. Position {entry['coordinate']} ({entry['strand']}) - Sigma {entry['sigma_factor']} - {entry['tss_id']}")

if __name__ == "__main__":
    main()
