"""
UI manager for genome analysis.
Handles user interaction, dashboard launching, and search mode prompts.
"""

import os
import sys
import multiprocessing
from typing import Tuple, Optional, List


def offer_dashboard() -> Tuple[Optional[multiprocessing.Queue], Optional[multiprocessing.Queue]]:
    """
    Offer to launch the dashboard for monitoring search progress.
    
    Returns:
        Tuple of (command_queue, dashboard_process) or (None, None) if declined
    """
    try:
        print("\n[INFO] Dashboard available for real-time search monitoring")
        print("[INFO] Launching dashboard automatically...")
        
        # Automatically launch dashboard (user requested to skip prompt)
        return launch_dashboard()
        
    except Exception as e:
        print(f"[ERROR] Failed to offer dashboard: {e}")
        return None, None


def launch_dashboard() -> Tuple[Optional[multiprocessing.Queue], Optional[multiprocessing.Queue]]:
    """
    Launch the worker dashboard for real-time monitoring.
    
    Returns:
        Tuple of (command_queue, dashboard_process) or (None, None) if failed
    """
    try:
        # Create communication queue
        command_queue = multiprocessing.Queue()
        
        # Import and start dashboard
        from .dashboard import start_dashboard
        
        # Start dashboard in separate process
        dashboard_process = multiprocessing.Process(
            target=start_dashboard,
            args=(command_queue,),
            daemon=True
        )
        dashboard_process.start()
        
        print(f"[INFO] Dashboard launched successfully (PID: {dashboard_process.pid})")
        return command_queue, dashboard_process
        
    except Exception as e:
        print(f"[ERROR] Failed to launch dashboard: {e}")
        return None, None


def search_mode_prompt(analyzer, chromosomes: List[str]) -> Tuple[str, str, int, int]:
    """
    Prompt user for search mode and parameters.
    
    Args:
        analyzer: GenomeAnalyzer instance
        chromosomes: List of available chromosomes
    
    Returns:
        Tuple of (mode, pattern, max_mismatches, boundary_bp)
    """
    try:
        print(f"\n[INFO] Available chromosomes: {len(chromosomes)}")
        print(f"[INFO] Chromosomes: {', '.join(chromosomes[:5])}{'...' if len(chromosomes) > 5 else ''}")
        
        print("\n=== SEARCH MODE SELECTION ===")
        print("1. Single chromosome search")
        print("2. Batch search across all chromosomes")
        print("3. Single genome search")
        
        while True:
            try:
                mode = input("\nSelect search mode (1-3): ").strip()
                if mode in ['1', '2', '3']:
                    break
                else:
                    print("[ERROR] Please enter 1, 2, or 3")
            except (EOFError, KeyboardInterrupt):
                print("\n[INFO] Using default mode: Batch search")
                mode = '2'
                break
        
        # Get search pattern
        pattern = input("Enter search pattern (DNA sequence): ").strip().upper()
        if not pattern:
            pattern = "ATCG"  # Default pattern
            print(f"[INFO] Using default pattern: {pattern}")
        
        # Get mismatch tolerance
        try:
            max_mismatches = int(input("Enter maximum mismatches (0 for exact): ").strip() or "0")
        except ValueError:
            max_mismatches = 0
            print("[INFO] Using exact matching (0 mismatches)")
        
        # Get boundary distance
        try:
            boundary_bp = int(input("Enter boundary distance in bp (0 for no boundary): ").strip() or "0")
        except ValueError:
            boundary_bp = 0
            print("[INFO] No boundary restrictions")
        
        return mode, pattern, max_mismatches, boundary_bp
        
    except Exception as e:
        print(f"[ERROR] Search mode prompt failed: {e}")
        # Return defaults
        return '2', 'ATCG', 0, 0


# Export functions
__all__ = [
    'offer_dashboard',
    'launch_dashboard',
    'search_mode_prompt'
]
