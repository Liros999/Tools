#!/usr/bin/env python3
"""
Debug script to test live output format and identify display issues.
"""

import os
import sys

def test_live_output_format():
    """Test the live output format functions."""
    
    print("🔍 Testing live output format...")
    
    # Import the functions
    try:
        from Genome_Analyzer import _print_results_header, _print_row_tuple, _print_bottom_frame
        print("✅ Successfully imported output functions")
    except Exception as e:
        print(f"❌ Error importing functions: {e}")
        return
    
    # Test header printing
    print(f"\n📋 Testing header printing:")
    print("="*80)
    
    # Test with chromosome (9 columns)
    print("🧬 Header WITH chromosome column:")
    _print_results_header(include_chrom=True)
    
    # Test without chromosome (8 columns)
    print("\n🧬 Header WITHOUT chromosome column:")
    _print_results_header(include_chrom=False)
    
    # Test row printing
    print(f"\n📊 Testing row printing:")
    print("="*80)
    
    # Test with chromosome data (9-tuple)
    print("🧬 Row WITH chromosome (9-tuple):")
    chrom_row = ("chr1", 1000, 1010, "+", "TATAAA", "TATAAA", "TTTATA", 0, "Gene1")
    _print_row_tuple(chrom_row, include_chrom=True)
    
    # Test without chromosome data (8-tuple)
    print("\n🧬 Row WITHOUT chromosome (8-tuple):")
    normal_row = (1000, 1010, "+", "TATAAA", "TATAAA", "TTTATA", 0, "Gene1")
    _print_row_tuple(normal_row, include_chrom=False)
    
    # Test bottom frame
    print(f"\n📋 Testing bottom frame:")
    _print_bottom_frame()
    
    # Test the exact format we want
    print(f"\n🎯 TESTING THE EXACT FORMAT WE WANT:")
    print("="*80)
    
    # Print header
    _print_results_header(include_chrom=True)
    
    # Print some sample results
    sample_results = [
        ("chr1", 1000, 1010, "+", "TATAAA", "TATAAA", "TTTATA", 0, "Gene1"),
        ("chr1", 5000, 5010, "-", "TATAAA", "TATAAA", "TTTATA", 0, "Gene2"),
        ("chr2", 1000, 1010, "+", "TATAAA", "TATAAA", "TTTATA", 0, "Gene3"),
    ]
    
    for result in sample_results:
        _print_row_tuple(result, include_chrom=True)
    
    # Print bottom frame
    _print_bottom_frame()
    
    print(f"\n✅ Live output format test complete!")

def test_parallel_output_simulation():
    """Simulate the parallel output to see what's happening."""
    
    print(f"\n🔄 Testing parallel output simulation...")
    print("="*80)
    
    # Simulate the parallel search output
    print(f"PARALLEL SEARCH PROGRESS - 3 CHROMOSOMES")
    print("="*80)
    
    # Print header for live output
    from Genome_Analyzer import _print_results_header
    _print_results_header(include_chrom=True)
    
    # Simulate worker progress
    print(f"\n{'='*80}")
    print(f"🎯 DASHBOARD MODE - 4 WORKERS")
    print(f"{'='*80}")
    print(f"📊 Dashboard window opened - monitor progress in real-time!")
    print(f"💡 You can minimize this console and watch the dashboard instead")
    print(f"{'='*80}")
    
    # Simulate results coming in
    print(f"\n{'='*80}")
    print(f"LIVE RESULTS OUTPUT")
    print(f"{'='*80}")
    
    # Simulate live results
    from Genome_Analyzer import _print_row_tuple
    
    print(f"[INFO] ✅ Chromosome chr1 completed (1/3) - 2 matches")
    results = [
        ("chr1", 1000, 1010, "+", "TATAAA", "TATAAA", "TTTATA", 0, "Gene1"),
        ("chr1", 5000, 5010, "-", "TATAAA", "TATAAA", "TTTATA", 0, "Gene2"),
    ]
    
    for result in results:
        _print_row_tuple(result, include_chrom=True)
    
    print(f"[INFO] ✅ Chromosome chr2 completed (2/3) - 1 match")
    results = [
        ("chr2", 1000, 1010, "+", "TATAAA", "TATAAA", "TTTATA", 0, "Gene3"),
    ]
    
    for result in results:
        _print_row_tuple(result, include_chrom=True)
    
    print(f"[INFO] ✅ Chromosome chr3 completed (3/3) - 0 matches")
    
    # Final summary
    print(f"\n{'='*80}")
    print(f"SEARCH RESULTS - 3 TOTAL MATCHES")
    print(f"{'='*80}")
    
    # Print header for results
    _print_results_header(include_chrom=True)
    
    # Print all results
    all_results = [
        ("chr1", 1000, 1010, "+", "TATAAA", "TATAAA", "TTTATA", 0, "Gene1"),
        ("chr1", 5000, 5010, "-", "TATAAA", "TATAAA", "TTTATA", 0, "Gene2"),
        ("chr2", 1000, 1010, "+", "TATAAA", "TATAAA", "TTTATA", 0, "Gene3"),
    ]
    
    for result in all_results:
        _print_row_tuple(result, include_chrom=True)
    
    from Genome_Analyzer import _print_bottom_frame
    _print_bottom_frame()
    
    print(f"\n{'='*80}")
    print(f"PARALLEL SEARCH COMPLETE - 3 TOTAL MATCHES")
    print(f"{'='*80}")
    
    print(f"\n✅ Parallel output simulation complete!")

if __name__ == "__main__":
    print("🔍 Debugging live output format issues...")
    test_live_output_format()
    test_parallel_output_simulation()
    print("\n✅ All debug tests complete!")
