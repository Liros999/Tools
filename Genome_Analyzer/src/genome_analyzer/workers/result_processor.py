"""
Result processing functions for genome analysis parallel processing.
Handles result validation, formatting, and output management.
"""

from typing import List, Tuple, Optional, callable

from ..analyzer import GenomeAnalyzer


def _process_search_results(matches: List[Tuple], pattern: str, strand: str, 
                          analyzer: GenomeAnalyzer, boundary_bp: int, 
                          on_row: Optional[callable] = None, stream: bool = False,
                          include_chrom: bool = False) -> List[Tuple]:
    """
    Centralized result processing to eliminate code duplication.
    Returns processed results list.
    """
    processed_results = []
    for start, end, window, _, mismatches, _ in matches:
        location = analyzer.get_location_details(start, end)
        
        # Apply boundary annotations if specified
        if boundary_bp and boundary_bp > 0:
            seq_len = len(analyzer.sequence)
            dist_start = start - 1
            dist_end = seq_len - end
            near_flags = []
            if dist_start <= boundary_bp:
                near_flags.append(f"NearStart:{dist_start}bp")
            if dist_end <= boundary_bp:
                near_flags.append(f"NearEnd:{dist_end}bp")
            if near_flags:
                location = f"{location} [{' & '.join(near_flags)}]"
        
        # Create result row
        row = (start, end, strand, pattern, window, analyzer._reverse_complement(window), mismatches, location)
        processed_results.append(row)
        
        # Handle output based on parameters
        if on_row:
            on_row(row)
        elif stream:
            # Import here to avoid circular imports
            from ..ui import _print_row_tuple
            _print_row_tuple(row, include_chrom)
    
    return processed_results


# Export function
__all__ = ['_process_search_results']
