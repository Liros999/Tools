"""
Dashboard communication for genome analysis parallel processing.
Handles real-time dashboard updates and worker progress tracking.
"""

import time
import multiprocessing
from typing import Optional, Dict, List


class DashboardCommunicator:
    """Handles real-time dashboard communication with proper message formatting."""
    
    def __init__(self, command_queue: Optional[multiprocessing.Queue] = None):
        self.command_queue = command_queue
        self.worker_states = {}
        self.overall_progress = {'completed': 0, 'total': 0, 'matches': 0}
        
    def is_available(self) -> bool:
        """Check if dashboard is available for communication."""
        return self.command_queue is not None
        
    def update_worker(self, worker_id: int, status: str, progress: int, details: str = ""):
        """Send worker update to dashboard."""
        if not self.is_available():
            return
            
        try:
            message = {
                'worker_id': worker_id,
                'status': status,
                'progress': max(0, min(100, progress)),
                'details': str(details)[:50],  # Limit detail length
                'timestamp': time.time()
            }
            
            self.command_queue.put(("UPDATE_WORKER", message), timeout=0.1)
            self.worker_states[worker_id] = message
            
        except Exception as e:
            print(f"[ERROR] Failed to update worker: {e}")
            
    def update_overall_progress(self, completed: int, total: int, matches: int):
        """Send overall progress update to dashboard."""
        if not self.is_available():
            return
            
        try:
            progress_data = {
                'completed': completed,
                'total': total,
                'matches': matches,
                'percentage': (completed / total * 100) if total > 0 else 0,
                'timestamp': time.time()
            }
            
            self.command_queue.put(("UPDATE_OVERALL", progress_data), timeout=0.1)
            self.overall_progress = progress_data
            
        except Exception as e:
            print(f"[ERROR] Failed to update overall progress: {e}")
            
    def send_search_info(self, pattern: str, chromosomes: List[str], algorithm: str):
        """Send search configuration info to dashboard."""
        if not self.is_available():
            return
            
        try:
            search_info = {
                'pattern': pattern,
                'chromosomes': chromosomes,
                'algorithm': algorithm,
                'start_time': time.time()
            }
            
            self.command_queue.put(("SEARCH_INFO", search_info), timeout=0.1)
            
        except Exception as e:
            print(f"[ERROR] Failed to send search info: {e}")
    
    def send_sequence_info(self, filename: str, total_sequences: int):
        """Send sequence file information to dashboard."""
        if not self.is_available():
            return
            
        try:
            sequence_info = {
                'sequence_file': filename,
                'total_sequences': total_sequences
            }
            
            self.command_queue.put(("SEQUENCE_UPDATE", sequence_info), timeout=0.1)
            
        except Exception as e:
            print(f"[ERROR] Failed to send sequence info: {e}")


# Export class
__all__ = ['DashboardCommunicator']
