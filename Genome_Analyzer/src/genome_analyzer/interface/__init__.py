"""
Interface package for genome analysis.
Contains UI components, user interaction, and dashboard functionality.
"""

from .ui_manager import offer_dashboard, launch_dashboard, search_mode_prompt
from .dashboard import start_dashboard
from .output_formatter import _print_results_header, _print_row_tuple

__all__ = [
    'offer_dashboard',
    'launch_dashboard', 
    'search_mode_prompt',
    'start_dashboard',
    '_print_results_header',
    '_print_row_tuple'
]
