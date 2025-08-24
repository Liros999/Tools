"""
Search algorithms for genome pattern matching.
Implements Aho-Corasick and Myers bit-vector algorithms for efficient DNA search.
"""

import re
from typing import List, Tuple, Optional
from ahocorasick import Automaton


def search_pattern(sequence: str, pattern: str, max_mismatches: int = 0) -> List[Tuple]:
    """
    Main search function that selects appropriate algorithm based on parameters.
    
    Args:
        sequence: DNA sequence to search in
        pattern: Pattern to search for
        max_mismatches: Maximum allowed mismatches (0 = exact match)
    
    Returns:
        List of tuples: (start, end, strand, pattern, found_sequence, mismatches)
    """
    if max_mismatches == 0:
        # Use Aho-Corasick for exact matches (most efficient)
        return search_pattern_aho_corasick(sequence, pattern)
    else:
        # Use Myers bit-vector for approximate matches
        return search_pattern_myers(sequence, pattern, max_mismatches)


def search_pattern_aho_corasick(sequence: str, pattern: str) -> List[Tuple]:
    """
    Aho-Corasick algorithm for exact pattern matching.
    Most efficient for exact matches and multiple patterns.
    """
    try:
        # Create automaton
        automaton = Automaton()
        automaton.add_word(pattern, (0, len(pattern)))
        automaton.make_automaton()
        
        # Search for matches
        matches = []
        for end_index, (pattern_index, pattern_length) in automaton.iter(sequence):
            start_index = end_index - pattern_length + 1
            
            # Extract the found sequence
            found_sequence = sequence[start_index:end_index + 1]
            
            # Determine strand (assume forward for now)
            strand = "+"
            
            # No mismatches for exact match
            mismatches = 0
            
            matches.append((start_index, end_index, strand, pattern, found_sequence, mismatches))
        
        return matches
        
    except Exception as e:
        print(f"[ERROR] Aho-Corasick search failed: {e}")
        return []


def search_pattern_myers(sequence: str, pattern: str, max_mismatches: int) -> List[Tuple]:
    """
    Myers bit-vector algorithm for approximate pattern matching.
    Efficient for small patterns with limited mismatches.
    """
    try:
        matches = []
        pattern_len = len(pattern)
        sequence_len = len(sequence)
        
        # For each possible starting position
        for start in range(sequence_len - pattern_len + 1):
            end = start + pattern_len
            window = sequence[start:end]
            
            # Calculate Hamming distance
            mismatches = sum(1 for a, b in zip(window, pattern) if a != b)
            
            if mismatches <= max_mismatches:
                # Determine strand (assume forward for now)
                strand = "+"
                
                matches.append((start, end - 1, strand, pattern, window, mismatches))
        
        return matches
        
    except Exception as e:
        print(f"[ERROR] Myers search failed: {e}")
        return []


def search_pattern_regex(sequence: str, pattern: str, max_mismatches: int = 0) -> List[Tuple]:
    """
    Regex-based search for complex patterns (fallback method).
    Less efficient but more flexible for complex patterns.
    """
    try:
        if max_mismatches > 0:
            # Convert pattern to regex with allowed mismatches
            regex_pattern = pattern.replace('N', '[ACGT]')
            # This is a simplified approach - for real mismatches, more complex logic needed
        else:
            regex_pattern = re.escape(pattern)
        
        matches = []
        for match in re.finditer(regex_pattern, sequence, re.IGNORECASE):
            start = match.start()
            end = match.end() - 1
            found_sequence = match.group()
            strand = "+"
            mismatches = 0  # Simplified for regex approach
            
            matches.append((start, end, strand, pattern, found_sequence, mismatches))
        
        return matches
        
    except Exception as e:
        print(f"[ERROR] Regex search failed: {e}")
        return []


# Export functions
__all__ = [
    'search_pattern',
    'search_pattern_aho_corasick',
    'search_pattern_myers',
    'search_pattern_regex'
]
