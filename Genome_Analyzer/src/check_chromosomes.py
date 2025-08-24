#!/usr/bin/env python3
"""
Check available chromosomes in GTF file
"""

import os

def check_chromosomes():
    gtf_file = '../data/Homo_Sapiens/gencode.v48.primary_assembly.annotation.gtf'
    
    if not os.path.exists(gtf_file):
        print(f"GTF file not found: {gtf_file}")
        return
    
    chromosomes = set()
    print("Scanning GTF file for chromosomes...")
    
    try:
        with open(gtf_file, 'r') as f:
            for line_num, line in enumerate(f):
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) > 0:
                    chromosomes.add(parts[0])
                
                # Scan entire file to find all chromosomes
                # No line limit - scan everything
        
        print(f"Found chromosomes: {sorted(list(chromosomes))}")
        print(f"Total unique chromosomes: {len(chromosomes)}")
        
        # Filter for main chromosomes only (not contigs)
        main_chromosomes = []
        for c in chromosomes:
            if c.startswith('chr'):
                if c[3:].isdigit() and 1 <= int(c[3:]) <= 22:  # chr1-chr22
                    main_chromosomes.append(c)
                elif c in ['chrX', 'chrY', 'chrM']:  # Sex chromosomes and mitochondrial
                    main_chromosomes.append(c)
        
        main_chromosomes.sort(key=lambda x: (x[3:] if x[3:].isdigit() else x[3:]))
        
        print(f"Main chromosomes (25 total): {main_chromosomes}")
        print(f"Total main chromosomes: {len(main_chromosomes)}")
        print(f"Contigs and other sequences: {len(chromosomes) - len(main_chromosomes)}")
        
    except Exception as e:
        print(f"Error reading GTF file: {e}")

if __name__ == "__main__":
    check_chromosomes()
