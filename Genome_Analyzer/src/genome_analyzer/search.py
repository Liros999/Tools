"""
Search module for genome analysis.
Redirects to the new algorithms package for better organization.
"""

# Import from the new algorithms package
from .algorithms import (
    search_pattern,
    search_pattern_aho_corasick,
    search_pattern_myers,
    search_pattern_regex
)

# Re-export for backward compatibility
__all__ = [
    'search_pattern',
    'search_pattern_aho_corasick',
    'search_pattern_myers',
    'search_pattern_regex'
]
