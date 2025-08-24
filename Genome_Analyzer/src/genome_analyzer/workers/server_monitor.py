"""
Server resource monitoring and management for genome analysis.
Handles adaptive server configuration and resource cleanup.
"""

import os
import platform
import psutil

# Linux compatibility and adaptive server resource management
IS_LINUX = platform.system() == 'Linux'
IS_SERVER = os.environ.get('SERVER_ENV', 'false').lower() == 'true'

# Import adaptive server configuration
try:
    import sys
    import os
    # Add the parent directory to the path to import linux_server_config
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
    from linux_server_config import ServerConfig, _adaptive_manager
    if _adaptive_manager:
        SERVER_MAX_WORKERS = _adaptive_manager.MAX_WORKERS
        SERVER_MEMORY_LIMIT_GB = _adaptive_manager.MAX_MEMORY_GB
        SERVER_CACHE_SIZE_MB = _adaptive_manager.CACHE_SIZE_MB
        print(f"[INFO] Using adaptive server configuration: {_adaptive_manager.server_class}")
    else:
        SERVER_MAX_WORKERS = min(4, os.cpu_count()) if IS_SERVER else os.cpu_count()
        SERVER_MEMORY_LIMIT_GB = 8
        SERVER_CACHE_SIZE_MB = 2048
except (ImportError, ModuleNotFoundError, AttributeError):
    # Fallback if config not available
    print("[INFO] Linux server config not available, using default configuration")
    SERVER_MAX_WORKERS = min(4, os.cpu_count()) if IS_SERVER else os.cpu_count()
    SERVER_MEMORY_LIMIT_GB = 8
    SERVER_CACHE_SIZE_MB = 2048
    # Set _adaptive_manager to None to avoid NameError
    _adaptive_manager = None


def _monitor_server_resources():
    """Monitor server resources and enforce adaptive limits."""
    if not IS_SERVER:
        return
    
    try:
        # Check memory usage
        memory = psutil.virtual_memory()
        memory_gb = memory.used / (1024**3)
        memory_percent = memory.percent
        
        # Get adaptive limits
        try:
            if '_adaptive_manager' in globals() and _adaptive_manager and hasattr(_adaptive_manager, 'MAX_MEMORY_GB'):
                memory_limit = _adaptive_manager.MAX_MEMORY_GB
                memory_warning = _adaptive_manager.MEMORY_WARNING_THRESHOLD_GB
                server_class = _adaptive_manager.server_class
            else:
                memory_limit = SERVER_MEMORY_LIMIT_GB
                memory_warning = SERVER_MEMORY_LIMIT_GB * 0.8
                server_class = "standard"
        except (NameError, AttributeError):
            memory_limit = SERVER_MEMORY_LIMIT_GB
            memory_warning = SERVER_MEMORY_LIMIT_GB * 0.8
            server_class = "standard"
        
        # Memory monitoring with adaptive thresholds
        if memory_gb > memory_limit:
            print(f"[WARNING] Memory usage {memory_gb:.1f}GB exceeds limit {memory_limit}GB")
            print(f"[INFO] Server class: {server_class.upper()} - Consider upgrading if this persists")
            _cleanup_server_resources()
        elif memory_gb > memory_warning:
            print(f"[INFO] Memory usage {memory_gb:.1f}GB approaching limit {memory_limit}GB ({memory_percent:.1f}%)")
            print(f"[INFO] Current utilization: {memory_percent:.1f}% - Optimal for {server_class.upper()} server")
        
        # Check CPU usage with adaptive thresholds
        cpu_percent = psutil.cpu_percent(interval=1)
        try:
            if '_adaptive_manager' in globals() and _adaptive_manager and hasattr(_adaptive_manager, 'MAX_CPU_PERCENT'):
                cpu_limit = _adaptive_manager.MAX_CPU_PERCENT
            else:
                cpu_limit = 80
        except (NameError, AttributeError):
            cpu_limit = 80
            
        if cpu_percent > cpu_limit:
            print(f"[WARNING] High CPU usage: {cpu_percent:.1f}% (limit: {cpu_limit}%)")
            if cpu_percent > 95:
                print(f"[WARNING] Critical CPU usage - consider reducing worker count")
        elif cpu_percent > cpu_limit * 0.8:
            print(f"[INFO] CPU usage: {cpu_percent:.1f}% - Optimal utilization for {server_class.upper()} server")
        
        # Adaptive resource recommendations
        try:
            if '_adaptive_manager' in globals() and _adaptive_manager and memory_percent < 50:
                print(f"[INFO] Memory utilization {memory_percent:.1f}% - {server_class.upper()} server can handle more load")
                print(f"[INFO] Consider increasing worker count or cache size for better performance")
        except (NameError, AttributeError):
            pass
            
    except Exception as e:
        print(f"[WARNING] Server resource monitoring failed: {e}")


def _cleanup_server_resources():
    """Clean up server resources to free memory."""
    if not IS_SERVER:
        return
    
    try:
        import gc
        gc.collect()  # Force garbage collection
        
        # Clear large caches
        from .process_pool import _global_process_pool, _global_process_pool_initialized
        if _global_process_pool is not None:
            try:
                _global_process_pool.shutdown(wait=True)
                _global_process_pool = None
                _global_process_pool_initialized = False
                print("[INFO] Process pool cleaned up to free memory")
            except:
                pass
                
    except Exception as e:
        print(f"[WARNING] Server resource cleanup failed: {e}")


# Export server configuration
__all__ = [
    'IS_LINUX', 
    'IS_SERVER', 
    'SERVER_MAX_WORKERS', 
    'SERVER_MEMORY_LIMIT_GB', 
    'SERVER_CACHE_SIZE_MB',
    '_monitor_server_resources',
    '_cleanup_server_resources'
]
