"""
Resource management for genome analysis operations.
Handles server resource monitoring, limits, and cleanup operations.
"""

import os
import platform
import psutil
import gc
import sys
import time

# Linux compatibility and resource management
IS_LINUX = platform.system() == 'Linux'
IS_SERVER = os.environ.get('SERVER_ENV', 'false').lower() == 'true'

# Resource limits for server environments
SERVER_MEMORY_LIMIT_GB = 8  # Maximum memory usage on servers
SERVER_CPU_LIMIT_PERCENT = 80  # Maximum CPU usage on servers
SERVER_CACHE_SIZE_MB = 2048  # Maximum cache size on servers (2GB)


class ResourceManager:
    """Manages server resources and prevents excessive consumption."""
    
    def __init__(self):
        self.memory_warning_threshold = SERVER_MEMORY_LIMIT_GB * 0.8  # 80% of limit
        self.cpu_warning_threshold = SERVER_CPU_LIMIT_PERCENT * 0.8  # 80% of limit
    
    def check_memory_usage(self) -> float:
        """Check current memory usage and return percentage."""
        try:
            if IS_LINUX:
                # Linux-specific memory monitoring
                with open('/proc/meminfo', 'r') as f:
                    lines = f.readlines()
                    total_mem = int(lines[0].split()[1]) * 1024  # Convert KB to bytes
                    available_mem = int(lines[2].split()[1]) * 1024
                    used_mem = total_mem - available_mem
                    return (used_mem / total_mem) * 100
            else:
                # Cross-platform memory monitoring
                memory = psutil.virtual_memory()
                return memory.percent
        except Exception:
            return 0.0
    
    def check_cpu_usage(self) -> float:
        """Check current CPU usage and return percentage."""
        try:
            return psutil.cpu_percent(interval=1)
        except Exception:
            return 0.0
    
    def enforce_resource_limits(self):
        """Enforce resource limits and trigger cleanup if needed."""
        memory_usage = self.check_memory_usage()
        cpu_usage = self.check_cpu_usage()
        
        # Memory cleanup if approaching limit
        if memory_usage > self.memory_warning_threshold:
            print(f"[WARNING] High memory usage: {memory_usage:.1f}%")
            self.cleanup_memory()
        
        # CPU throttling if approaching limit
        if cpu_usage > self.cpu_warning_threshold:
            print(f"[WARNING] High CPU usage: {cpu_usage:.1f}%")
            self.throttle_processing()
    
    def cleanup_memory(self):
        """Aggressive memory cleanup for server environments."""
        if IS_SERVER:
            print("[INFO] Performing server memory cleanup...")
            gc.collect()  # Force garbage collection
            
            # Clear Python object caches
            for module in list(sys.modules.keys()):
                if module.startswith('genome_analyzer'):
                    try:
                        del sys.modules[module]
                    except:
                        pass
    
    def throttle_processing(self):
        """Throttle processing to reduce CPU usage."""
        if IS_SERVER:
            time.sleep(0.1)  # Small delay to reduce CPU load


# Export class and constants
__all__ = [
    'ResourceManager',
    'IS_LINUX',
    'IS_SERVER',
    'SERVER_MEMORY_LIMIT_GB',
    'SERVER_CPU_LIMIT_PERCENT',
    'SERVER_CACHE_SIZE_MB'
]
