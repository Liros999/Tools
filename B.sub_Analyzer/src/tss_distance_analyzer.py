#!/usr/bin/env python3
"""
TSS Distance Analyzer for B.sub_Analyzer
Adds distance to closest TSS for each search result with proper strand orientation logic.
"""

import json
import pandas as pd
import os
from typing import List, Tuple, Dict, Optional

class TSSDistanceAnalyzer:
    def __init__(self, tss_mapping_file: str = None):
        """Initialize TSS distance analyzer with TSS coordinate mapping."""
        if tss_mapping_file is None:
            # Default path to the TSS mapping file we created
            current_dir = os.path.dirname(os.path.abspath(__file__))
            tss_mapping_file = os.path.join(os.path.dirname(current_dir), "data", "tss_coordinate_mapping.json")
        
        self.tss_mapping = self._load_tss_mapping(tss_mapping_file)
        self.tss_list = self._create_sorted_tss_list()
        
    def _load_tss_mapping(self, tss_mapping_file: str) -> Dict:
        """Load TSS coordinate mapping from JSON file."""
        try:
            with open(tss_mapping_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print(f"Warning: TSS mapping file {tss_mapping_file} not found.")
            return {}
        except Exception as e:
            print(f"Error loading TSS mapping: {e}")
            return {}
            
    def _create_sorted_tss_list(self) -> List[Tuple[int, str, str]]:
        """Create a sorted list of TSS coordinates for efficient searching."""
        tss_list = []
        for coord_str, tss_entries in self.tss_mapping.items():
            coord = int(coord_str)
            for entry in tss_entries:
                tss_list.append((coord, entry['strand'], entry['sigma_factor']))
        
        # Sort by coordinate for binary search
        tss_list.sort(key=lambda x: x[0])
        return tss_list
    
    def find_closest_tss(self, position: int, strand: str) -> Tuple[Optional[int], Optional[str], Optional[str], Optional[int]]:
        """
        Find the closest TSS to a given position with improved strand consideration.
        Returns (tss_position, tss_strand, sigma_factor, distance) or (None, None, None, None) if no TSS found.
        
        Priority: 1) Same strand TSS, 2) Closest TSS regardless of strand
        """
        if not self.tss_list:
            return None, None, None, None
            
        # Find all TSS within a reasonable distance (e.g., 50kb)
        max_distance = 50000
        same_strand_tss = []
        all_tss = []
        
        for tss_coord, tss_strand, sigma_factor in self.tss_list:
            distance = abs(tss_coord - position)
            if distance <= max_distance:
                all_tss.append((tss_coord, tss_strand, sigma_factor, distance))
                # Prioritize TSS on the same strand
                if tss_strand == strand:
                    same_strand_tss.append((tss_coord, tss_strand, sigma_factor, distance))
        
        if not all_tss:
            return None, None, None, None
        
        # Prefer same strand TSS if available and reasonably close
        if same_strand_tss:
            same_strand_tss.sort(key=lambda x: x[3])
            closest_same_strand = same_strand_tss[0]
            
            # If same-strand TSS is within 10kb, prefer it
            if closest_same_strand[3] <= 10000:
                return closest_same_strand[0], closest_same_strand[1], closest_same_strand[2], closest_same_strand[3]
        
        # Otherwise, return the absolutely closest TSS
        all_tss.sort(key=lambda x: x[3])
        closest_tss = all_tss[0]
        return closest_tss[0], closest_tss[1], closest_tss[2], closest_tss[3]
    
    def calculate_tss_distance_string(self, match_start: int, match_end: int, match_strand: str) -> str:
        """
        Calculate distance string for a match considering strand orientation.
        
        Logic:
        - For forward strand (+): upstream = lower coordinates, downstream = higher coordinates
        - For reverse strand (-): upstream = higher coordinates, downstream = lower coordinates
        - Format: "TSS >>> bp_count" if TSS is upstream, "bp_count <<< TSS" if TSS is downstream
        """
        # Use the start position of the match for distance calculation
        match_position = match_start
        
        tss_position, tss_strand, sigma_factor, distance = self.find_closest_tss(match_position, match_strand)
        
        if tss_position is None:
            return "No TSS within 50kb"
            
        # Calculate actual distance (not absolute)
        bp_distance = abs(tss_position - match_position)
        
        # Determine if TSS is upstream or downstream based on strand
        if match_strand == '+':
            # Forward strand: upstream = lower coordinates, downstream = higher coordinates
            if tss_position < match_position:
                # TSS is upstream (lower coordinate)
                direction_string = f"TSS >>> {bp_distance}"
            else:
                # TSS is downstream (higher coordinate)
                direction_string = f"{bp_distance} <<< TSS"
        else:  # match_strand == '-'
            # Reverse strand: upstream = higher coordinates, downstream = lower coordinates
            if tss_position > match_position:
                # TSS is upstream (higher coordinate for reverse strand)
                direction_string = f"TSS >>> {bp_distance}"
            else:
                # TSS is downstream (lower coordinate for reverse strand)
                direction_string = f"{bp_distance} <<< TSS"
        
        # Add sigma factor and strand information for better debugging
        sigma_info = f" (σ{sigma_factor})" if sigma_factor and sigma_factor != "Unknown" else ""
        tss_strand_info = f"[{tss_strand}]" if tss_strand != match_strand else ""
        
        return f"{direction_string}{sigma_info}{tss_strand_info}"
    
    def enhance_search_results(self, results_file: str, output_file: str = None) -> str:
        """
        Enhance existing search results CSV with TSS distance information.
        
        Args:
            results_file: Path to existing search results CSV
            output_file: Path for enhanced results (if None, appends '_with_TSS' to original name)
            
        Returns:
            Path to the enhanced results file
        """
        if not os.path.exists(results_file):
            raise FileNotFoundError(f"Results file not found: {results_file}")
            
        # Load existing results
        df = pd.read_csv(results_file)
        
        # Verify expected columns exist
        required_columns = ['Start', 'End', 'Strand', 'Location']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        print(f"Enhancing {len(df)} search results with TSS distance information...")
        
        # Calculate TSS distances
        tss_distances = []
        for _, row in df.iterrows():
            start = int(row['Start'])
            end = int(row['End'])
            strand = str(row['Strand'])
            
            distance_string = self.calculate_tss_distance_string(start, end, strand)
            tss_distances.append(distance_string)
        
        # Insert TSS distance column after Location
        location_index = df.columns.get_loc('Location')
        
        # Create new dataframe with TSS distance column inserted
        columns_before = df.columns[:location_index + 1].tolist()
        columns_after = df.columns[location_index + 1:].tolist()
        
        new_df = df[columns_before].copy()
        new_df['Distance_from_TSS'] = tss_distances
        
        for col in columns_after:
            new_df[col] = df[col]
        
        # Determine output file path
        if output_file is None:
            base_name = os.path.splitext(results_file)[0]
            output_file = f"{base_name}_with_TSS.csv"
        
        # Save enhanced results
        new_df.to_csv(output_file, index=False)
        print(f"Enhanced results saved to: {output_file}")
        
        # Print statistics
        tss_found = sum(1 for dist in tss_distances if "No TSS" not in dist)
        print(f"TSS distances calculated for {tss_found}/{len(tss_distances)} results")
        
        return output_file
    
    def enhance_live_results(self, results: List[Tuple]) -> List[Tuple]:
        """
        Enhance live search results from B.sub_Analyzer with TSS distance.
        
        Args:
            results: List of tuples (start, end, strand, pattern, found_sequence, location, mismatches)
            
        Returns:
            Enhanced list with TSS distance inserted after location
        """
        enhanced_results = []
        
        for result in results:
            start, end, strand, pattern, found_sequence, location, mismatches = result
            
            # Calculate TSS distance
            tss_distance = self.calculate_tss_distance_string(start, end, strand)
            
            # Insert TSS distance after location (index 5)
            enhanced_result = (start, end, strand, pattern, found_sequence, location, tss_distance, mismatches)
            enhanced_results.append(enhanced_result)
        
        return enhanced_results
    
    def print_enhanced_results(self, results: List[Tuple], display_count: int = None):
        """Print enhanced results with TSS distance in a formatted table."""
        if not results:
            print("No results to display.")
            return
            
        if display_count is None:
            display_count = len(results)
        
        # Headers for enhanced results
        print("\nEnhanced Search Results with TSS Distance:")
        print("-" * 180)
        print(f"{'Start':<10}{'End':<10}{'Strand':<8}{'Pattern':<18}{'Found Sequence':<20}{'Location':<35}{'Distance from TSS':<25}{'Mismatches':<12}")
        print("-" * 180)
        
        for i, result in enumerate(results[:display_count]):
            start, end, strand, pattern, found_sequence, location, tss_distance, mismatches = result
            
            # Truncate long strings for display
            pattern_display = pattern[:16] + "..." if len(pattern) > 16 else pattern
            sequence_display = found_sequence[:18] + "..." if len(found_sequence) > 18 else found_sequence
            location_display = location[:33] + "..." if len(location) > 33 else location
            tss_display = tss_distance[:23] + "..." if len(tss_distance) > 23 else tss_distance
            
            print(f"{str(start):<10}{str(end):<10}{strand:<8}{pattern_display:<18}{sequence_display:<20}{location_display:<35}{tss_display:<25}{str(mismatches):<12}")
        
        print("-" * 180)
        print(f"Total enhanced results shown: {min(display_count, len(results))}")

def enhance_existing_results(results_directory: str = None):
    """
    Utility function to enhance all CSV files in the results directory.
    """
    if results_directory is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        results_directory = os.path.join(os.path.dirname(current_dir), "results")
    
    if not os.path.exists(results_directory):
        print(f"Results directory not found: {results_directory}")
        return
    
    # Initialize TSS analyzer
    analyzer = TSSDistanceAnalyzer()
    
    # Find all CSV files in results directory
    csv_files = [f for f in os.listdir(results_directory) if f.endswith('.csv') and not f.endswith('_with_TSS.csv')]
    
    if not csv_files:
        print("No CSV files found in results directory.")
        return
    
    print(f"Found {len(csv_files)} CSV files to enhance:")
    for csv_file in csv_files:
        print(f"  - {csv_file}")
    
    # Enhance each file
    for csv_file in csv_files:
        file_path = os.path.join(results_directory, csv_file)
        try:
            enhanced_file = analyzer.enhance_search_results(file_path)
            print(f"✓ Enhanced: {csv_file} -> {os.path.basename(enhanced_file)}")
        except Exception as e:
            print(f"✗ Error enhancing {csv_file}: {e}")

def main():
    """Main function for standalone usage."""
    import sys
    
    if len(sys.argv) > 1:
        # Enhance specific file
        results_file = sys.argv[1]
        analyzer = TSSDistanceAnalyzer()
        try:
            enhanced_file = analyzer.enhance_search_results(results_file)
            print(f"Successfully enhanced {results_file}")
            print(f"Enhanced results saved to: {enhanced_file}")
        except Exception as e:
            print(f"Error: {e}")
    else:
        # Enhance all files in results directory
        enhance_existing_results()

if __name__ == "__main__":
    main()
