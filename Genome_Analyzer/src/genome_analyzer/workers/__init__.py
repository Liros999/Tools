"""
Workers package for genome analysis parallel processing.
Contains process pool management, worker functions, and dashboard communication.
"""

from .process_pool import _get_global_process_pool, _cleanup_global_process_pool
from .dashboard_communicator import DashboardCommunicator
from .chromosome_worker import enhanced_search_single_chromosome_worker
from .chunk_worker import _search_genome_chunk_worker, _search_genome_chunk_worker_simple
from .server_monitor import _monitor_server_resources, _cleanup_server_resources

__all__ = [
    '_get_global_process_pool',
    '_cleanup_global_process_pool', 
    'DashboardCommunicator',
    'enhanced_search_single_chromosome_worker',
    '_search_genome_chunk_worker',
    '_search_genome_chunk_worker_simple',
    '_monitor_server_resources',
    '_cleanup_server_resources'
]
