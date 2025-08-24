"""
Output formatter for genome analysis.
Handles formatting and display of search results and progress information.
"""

from typing import Tuple, Optional


def _print_results_header(include_chrom: bool = False):
    """
    Print the header for search results table.
    
    Args:
        include_chrom: Whether to include chromosome column in header
    """
    try:
        if include_chrom:
            header = "+--------+--------+--------+--------+--------+--------+--------+--------+--------+"
            print(header)
            print("| Chrom  | Start  |  End   |Strand  |Pattern | Found  |RevComp |Mismatch|Location|")
            print(header)
        else:
            header = "+--------+--------+--------+--------+--------+--------+--------+--------+"
            print(header)
            print("| Start  |  End   |Strand  |Pattern | Found  |RevComp |Mismatch|Location|")
            print(header)
            
    except Exception as e:
        print(f"[ERROR] Failed to print results header: {e}")


def _print_row_tuple(row: Tuple, include_chrom: bool = False):
    """
    Print a formatted row of search results.
    
    Args:
        row: Tuple containing search result data
        include_chrom: Whether to include chromosome column
    """
    try:
        if not row or len(row) < 8:
            print(f"[WARNING] Invalid row data: {row}")
            return
        
        if include_chrom and len(row) >= 9:
            # Format with chromosome: (chrom, start, end, strand, pattern, found, revcomp, mismatches, location)
            chrom, start, end, strand, pattern, found, revcomp, mismatches, location = row[:9]
            print(f"| {str(chrom)[:6]:>6} | {str(start):>6} | {str(end):>6} | {strand:>6} | "
                  f"{str(pattern)[:6]:>6} | {str(found)[:6]:>6} | {str(revcomp)[:6]:>6} | "
                  f"{str(mismatches):>6} | {str(location)[:6]:>6} |")
        else:
            # Format without chromosome: (start, end, strand, pattern, found, revcomp, mismatches, location)
            start, end, strand, pattern, found, revcomp, mismatches, location = row[:8]
            print(f"| {str(start):>6} | {str(end):>6} | {strand:>6} | "
                  f"{str(pattern)[:6]:>6} | {str(found)[:6]:>6} | {str(revcomp)[:6]:>6} | "
                  f"{str(mismatches):>6} | {str(location)[:6]:>6} |")
            
    except Exception as e:
        print(f"[ERROR] Failed to print row: {e}")


def format_search_summary(total_matches: int, total_searched: int, search_time: float) -> str:
    """
    Format a summary of search results.
    
    Args:
        total_matches: Total number of matches found
        total_searched: Total number of positions searched
        search_time: Time taken for search in seconds
    
    Returns:
        Formatted summary string
    """
    try:
        summary = f"\n=== SEARCH SUMMARY ===\n"
        summary += f"Total matches found: {total_matches:,}\n"
        summary += f"Total positions searched: {total_searched:,}\n"
        summary += f"Search time: {search_time:.2f} seconds\n"
        
        if total_searched > 0:
            match_rate = (total_matches / total_searched) * 100
            summary += f"Match rate: {match_rate:.4f}%\n"
        
        return summary
        
    except Exception as e:
        return f"[ERROR] Failed to format summary: {e}"


def format_progress_bar(current: int, total: int, width: int = 50) -> str:
    """
    Format a progress bar for display.
    
    Args:
        current: Current progress value
        total: Total value for 100%
        width: Width of progress bar in characters
    
    Returns:
        Formatted progress bar string
    """
    try:
        if total <= 0:
            return "[" + " " * width + "] 0%"
        
        percentage = min(100, (current / total) * 100)
        filled_width = int((percentage / 100) * width)
        
        bar = "[" + "=" * filled_width + " " * (width - filled_width) + "]"
        return f"{bar} {percentage:.1f}%"
        
    except Exception as e:
        return f"[ERROR] Failed to format progress bar: {e}"


# Export functions
__all__ = [
    '_print_results_header',
    '_print_row_tuple',
    'format_search_summary',
    'format_progress_bar'
]
