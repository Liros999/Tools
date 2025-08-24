"""
Process pool management for genome analysis parallel processing.
Handles global process pool creation, cleanup, and OS compatibility.
"""

import os
import multiprocessing
import atexit
import time
import platform
from concurrent.futures import ProcessPoolExecutor

# Import server configuration
from .server_monitor import IS_LINUX, IS_SERVER, SERVER_MAX_WORKERS

# Global process pool for parallel processing
_global_process_pool = None
_global_process_pool_workers = None
_global_process_pool_initialized = False


def _get_global_process_pool(worker_count=None):
    """Get or create a global process pool with Linux compatibility and server resource management."""
    global _global_process_pool, _global_process_pool_workers, _global_process_pool_initialized
    
    print(f"[DEBUG] _get_global_process_pool called with worker_count={worker_count}")
    print(f"[DEBUG] Current state: pool={_global_process_pool is not None}, workers={_global_process_pool_workers}, initialized={_global_process_pool_initialized}")
    
    if worker_count is None:
        # Use server-appropriate worker count
        if IS_SERVER:
            worker_count = SERVER_MAX_WORKERS
            print(f"[INFO] Server environment detected - limiting workers to {worker_count}")
        else:
            worker_count = os.cpu_count()
        print(f"[DEBUG] Determined worker_count={worker_count}")
    
    # Create new pool if it doesn't exist or worker count changed
    if (_global_process_pool is None or 
        _global_process_pool_workers != worker_count or 
        not _global_process_pool_initialized):
        
        print(f"[DEBUG] Need to create new process pool")
        
        # Clean up existing pool if it exists
        if _global_process_pool is not None:
            try:
                print(f"[DEBUG] Shutting down existing process pool")
                _global_process_pool.shutdown(wait=True)
                print(f"[DEBUG] Existing process pool shut down")
            except Exception as e:
                print(f"[DEBUG] Error shutting down existing pool: {e}")
        
        print(f"[INFO] Creating global process pool with {worker_count} workers")
        print(f"[DEBUG] Platform: {platform.system()}, IS_LINUX={IS_LINUX}")
        
        # Create process pool with OS compatibility
        try:
            print(f"[DEBUG] Attempting to create ProcessPoolExecutor...")
            start_time = time.time()
            
            if IS_LINUX:
                # Linux-specific process pool configuration
                print(f"[DEBUG] Using Linux fork context")
                _global_process_pool = ProcessPoolExecutor(
                    max_workers=worker_count,
                    mp_context=multiprocessing.get_context('fork')  # Use fork on Linux
                )
            else:
                # Windows/other OS process pool - use spawn context
                print(f"[DEBUG] Using Windows spawn context")
                _global_process_pool = ProcessPoolExecutor(
                    max_workers=worker_count,
                    mp_context=multiprocessing.get_context('spawn')  # Use spawn on Windows
                )
            
            elapsed = time.time() - start_time
            print(f"[DEBUG] ProcessPoolExecutor created in {elapsed:.3f} seconds")
            print(f"[INFO] Process pool created successfully with {worker_count} workers")
            
        except Exception as e:
            print(f"[ERROR] Failed to create process pool: {e}")
            print(f"[DEBUG] Exception type: {type(e).__name__}")
            print(f"[DEBUG] Exception details: {str(e)}")
            
            # Fallback to default context
            try:
                print(f"[DEBUG] Attempting fallback ProcessPoolExecutor...")
                _global_process_pool = ProcessPoolExecutor(max_workers=worker_count)
                print(f"[INFO] Using fallback process pool with {worker_count} workers")
            except Exception as fallback_e:
                print(f"[ERROR] Fallback also failed: {fallback_e}")
                raise fallback_e
        
        _global_process_pool_workers = worker_count
        _global_process_pool_initialized = True
        print(f"[DEBUG] Process pool state updated: workers={_global_process_pool_workers}, initialized={_global_process_pool_initialized}")
    
    else:
        print(f"[DEBUG] Using existing process pool with {_global_process_pool_workers} workers")
    
    print(f"[DEBUG] Returning process pool: {type(_global_process_pool)}")
    return _global_process_pool


def _cleanup_global_process_pool():
    """Clean up the global process pool on program exit."""
    global _global_process_pool
    if _global_process_pool is not None:
        try:
            print("[INFO] Shutting down global process pool")
            _global_process_pool.shutdown(wait=True)
        except:
            pass


# Register cleanup function
atexit.register(_cleanup_global_process_pool)


# Export functions
__all__ = [
    '_get_global_process_pool',
    '_cleanup_global_process_pool'
]
