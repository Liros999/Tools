#!/usr/bin/env python3
"""
Debug script to test chr1 loading and identify corruption issues.
"""

import os
import sys

def test_chr1_loading():
    """Test chr1 loading from the FASTA file."""
    
    # Get the data directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    data_dir = os.path.join(project_root, "data")
    
    # Check if Homo_sapiens directory exists
    hs_dir = os.path.join(data_dir, "Homo_sapiens")
    if not os.path.exists(hs_dir):
        print(f"❌ Homo_sapiens directory not found at: {hs_dir}")
        return
    
    # Look for genome file
    genome_files = [f for f in os.listdir(hs_dir) if f.endswith('.fa') or f.endswith('.fasta')]
    if not genome_files:
        print(f"❌ No genome files found in: {hs_dir}")
        return
    
    genome_file = genome_files[0]
    genome_path = os.path.join(hs_dir, genome_file)
    print(f"📁 Testing genome file: {genome_file}")
    
    # Test basic file reading
    try:
        file_size = os.path.getsize(genome_path)
        print(f"📊 File size: {file_size / (1024*1024*1024):.2f} GB")
        
        # Check first few lines
        with open(genome_path, 'r', encoding='utf-8', errors='ignore') as f:
            first_lines = []
            for i, line in enumerate(f):
                if i >= 10:  # Read first 10 lines
                    break
                first_lines.append(line.strip())
        
        print(f"\n📝 First 10 lines:")
        for i, line in enumerate(first_lines):
            print(f"  {i+1:2d}: {line[:100]}{'...' if len(line) > 100 else ''}")
        
        # Look for chr1 specifically
        print(f"\n🔍 Searching for chr1 entries...")
        chr1_found = False
        chr1_headers = []
        
        with open(genome_path, 'r', encoding='utf-8', errors='ignore') as f:
            for i, line in enumerate(f):
                if line.startswith('>chr1'):
                    chr1_found = True
                    chr1_headers.append(line.strip())
                    if len(chr1_headers) >= 3:  # Show first 3 chr1 headers
                        break
                if i > 1000000:  # Limit search to first 1M lines
                    break
        
        if chr1_found:
            print(f"✅ chr1 found! First headers:")
            for header in chr1_headers:
                print(f"  {header}")
        else:
            print(f"❌ chr1 not found in first 1M lines")
            
        # Test chromosome loading logic
        print(f"\n🧬 Testing chromosome loading logic...")
        from Genome_Analyzer import EColiAnalyzer
        
        analyzer = EColiAnalyzer()
        try:
            analyzer.load_chromosome_sequence(genome_path, "chr1")
            print(f"✅ chr1 loaded successfully: {len(analyzer.sequence):,} bp")
            
            # Check sequence content
            if analyzer.sequence:
                print(f"📊 Sequence stats:")
                print(f"  - Length: {len(analyzer.sequence):,} bp")
                print(f"  - First 100 bp: {analyzer.sequence[:100]}")
                print(f"  - Last 100 bp: {analyzer.sequence[-100:]}")
                print(f"  - Valid nucleotides: {sum(1 for c in analyzer.sequence if c in 'ACGTN'):,}")
                print(f"  - Invalid characters: {sum(1 for c in analyzer.sequence if c not in 'ACGTN'):,}")
            else:
                print(f"❌ Sequence is empty!")
                
        except Exception as e:
            print(f"❌ Error loading chr1: {e}")
            import traceback
            traceback.print_exc()
            
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    print("🔍 Debugging chr1 loading issues...")
    test_chr1_loading()
    print("\n✅ Debug complete!")
