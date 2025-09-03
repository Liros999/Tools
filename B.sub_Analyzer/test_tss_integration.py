#!/usr/bin/env python3
"""
Test script to verify TSS distance integration with B.sub_Analyzer
"""

import os
import sys
sys.path.append('src')

from src.tss_distance_analyzer import TSSDistanceAnalyzer

def test_tss_distance_calculation():
    """Test the TSS distance calculation functionality."""
    print("Testing TSS Distance Analyzer...")
    
    try:
        # Initialize TSS analyzer
        analyzer = TSSDistanceAnalyzer()
        print(f"✓ TSS analyzer initialized successfully")
        print(f"✓ Loaded {len(analyzer.tss_list)} TSS entries")
        
        # Test distance calculation for some known coordinates
        test_cases = [
            # (position, strand, description)
            (3443867, '+', 'copZ TSS-2762 position (forward strand)'),
            (3443867, '-', 'copZ TSS-2762 position (reverse strand)'),
            (3443613, '+', 'copZ gene start (forward strand)'),
            (3443613, '-', 'copZ gene start (reverse strand)'),
            (1000, '+', 'Early genome position (forward strand)'),
            (4000000, '-', 'Late genome position (reverse strand)'),
        ]
        
        print("\nTesting distance calculations:")
        for position, strand, description in test_cases:
            distance_string = analyzer.calculate_tss_distance_string(position, position + 10, strand)
            print(f"  {description}")
            print(f"    Position: {position}, Strand: {strand}")
            print(f"    Result: {distance_string}")
            print()
        
        # Test with a mock search result
        print("Testing with mock search results:")
        mock_results = [
            (3443613, 3443622, '+', 'ATCGATCGAT', 'ATCGATCGAT', '[copZ](in)', 0),
            (3443613, 3443622, '-', 'ATCGATCGAT', 'ATCGATCGAT', '[copZ](in)', 0),
            (100000, 100010, '+', 'TACNNNNGTA', 'TACGGGGGTA', '[gene1][-------][gene2]', 1),
        ]
        
        enhanced_results = analyzer.enhance_live_results(mock_results)
        print(f"✓ Enhanced {len(enhanced_results)} mock results")
        
        for result in enhanced_results:
            start, end, strand, pattern, found_seq, location, tss_distance, mismatches = result
            print(f"  Position: {start}-{end} ({strand})")
            print(f"  Location: {location}")
            print(f"  TSS Distance: {tss_distance}")
            print()
        
        return True
        
    except Exception as e:
        print(f"✗ Error testing TSS analyzer: {e}")
        return False

def test_file_enhancement():
    """Test enhancing existing CSV files."""
    print("Testing file enhancement capabilities...")
    
    # Check if there are any existing result files to enhance
    results_dir = "results"
    if os.path.exists(results_dir):
        csv_files = [f for f in os.listdir(results_dir) if f.endswith('.csv') and not f.endswith('_with_TSS.csv')]
        if csv_files:
            print(f"Found {len(csv_files)} existing CSV files:")
            for csv_file in csv_files[:3]:  # Test first 3 files
                print(f"  - {csv_file}")
            
            print("\nTesting enhancement on first file...")
            try:
                analyzer = TSSDistanceAnalyzer()
                test_file = os.path.join(results_dir, csv_files[0])
                enhanced_file = analyzer.enhance_search_results(test_file)
                print(f"✓ Successfully enhanced {csv_files[0]}")
                print(f"✓ Enhanced file saved as: {os.path.basename(enhanced_file)}")
                return True
            except Exception as e:
                print(f"✗ Error enhancing file: {e}")
                return False
        else:
            print("No existing CSV files found to test enhancement")
            return True
    else:
        print("No results directory found")
        return True

def main():
    """Run all tests."""
    print("=" * 60)
    print("B.sub_Analyzer TSS Integration Test")
    print("=" * 60)
    
    # Test basic TSS functionality
    test1_passed = test_tss_distance_calculation()
    
    print("\n" + "=" * 60)
    
    # Test file enhancement
    test2_passed = test_file_enhancement()
    
    print("\n" + "=" * 60)
    print("Test Summary:")
    print(f"TSS Distance Calculation: {'✓ PASSED' if test1_passed else '✗ FAILED'}")
    print(f"File Enhancement: {'✓ PASSED' if test2_passed else '✗ FAILED'}")
    
    if test1_passed and test2_passed:
        print("\n🎉 All tests passed! TSS integration is working correctly.")
    else:
        print("\n⚠️  Some tests failed. Check the error messages above.")
    
    print("\nUsage:")
    print("1. Run B.sub_Analyzer as usual with: python src/B_Sub2_0.py")
    print("2. The 'Distance from TSS' column will automatically appear in results")
    print("3. CSV files will include TSS distance information")
    print("4. To enhance existing CSV files: python src/tss_distance_analyzer.py [filename]")

if __name__ == "__main__":
    main()

