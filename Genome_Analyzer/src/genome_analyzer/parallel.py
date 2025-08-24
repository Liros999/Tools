"""
Parallel processing functions for genome analysis.
Refactored to use modular workers package for better maintainability.
"""

# Import all parallel processing functions from the workers package
from .workers import (
    _get_global_process_pool,
    _cleanup_global_process_pool,
    DashboardCommunicator,
    enhanced_search_single_chromosome_worker,
    _search_genome_chunk_worker,
    _search_genome_chunk_worker_simple,
    _process_search_results,
    enhanced_parallel_chromosome_search,
    parallel_single_genome_search,
    _monitor_server_resources,
    _cleanup_server_resources
)

# Re-export for backward compatibility
__all__ = [
    '_get_global_process_pool',
    '_cleanup_global_process_pool',
    'DashboardCommunicator',
    'enhanced_search_single_chromosome_worker',
    '_search_genome_chunk_worker',
    '_search_genome_chunk_worker_simple',
    '_process_search_results',
    'enhanced_parallel_chromosome_search',
    'parallel_single_genome_search',
    '_monitor_server_resources',
    '_cleanup_server_resources'
]
