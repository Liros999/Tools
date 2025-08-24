"""
User interface for genome analysis.
Refactored to use modular interface package for better maintainability.
"""

# Import all UI functions from the interface package
from .interface import (
    offer_dashboard,
    launch_dashboard,
    search_mode_prompt,
    start_dashboard,
    _print_results_header,
    _print_row_tuple
)

# Re-export for backward compatibility
__all__ = [
    'offer_dashboard',
    'launch_dashboard',
    'search_mode_prompt',
    'start_dashboard',
    '_print_results_header',
    '_print_row_tuple'
]