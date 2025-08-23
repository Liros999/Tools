#!/usr/bin/env python3
"""
Worker Dashboard for Genome Analyzer
Runs as a separate process to avoid tkinter threading issues on Windows.
"""

# WORKER DASHBOARD FIXES FOR LIVE UPDATES
# These fixes ensure the dashboard properly receives and displays live updates

import tkinter as tk
from tkinter import ttk
import threading
import multiprocessing
import time
from typing import Dict, Any, Optional

class EnhancedWorkerDashboard:
    """Enhanced worker dashboard with proper live update handling."""
    
    def __init__(self, command_queue: multiprocessing.Queue, worker_count: int = 8):
        self.command_queue = command_queue
        self.worker_count = worker_count
        self.running = True
        
        # Initialize data structures
        self.worker_data = {}
        self.overall_progress = {'completed': 0, 'total': 0, 'matches': 0, 'percentage': 0}
        self.search_info = {'pattern': '', 'chromosomes': [], 'algorithm': ''}
        
        # Create GUI
        self.root = tk.Tk()
        self.root.title("Genome Search Progress Dashboard")
        self.root.geometry("1000x700")
        self.root.configure(bg='#2b2b2b')
        
        # Set up GUI components
        self._create_widgets()
        
        # Start update thread
        self.update_thread = threading.Thread(target=self._message_processor, daemon=True)
        self.update_thread.start()
        
        # Start GUI refresh timer
        self._schedule_gui_refresh()
    
    def _create_widgets(self):
        """Create the dashboard GUI components."""
        
        # Title frame
        title_frame = tk.Frame(self.root, bg='#2b2b2b', height=80)
        title_frame.pack(fill='x', padx=10, pady=5)
        title_frame.pack_propagate(False)
        
        title_label = tk.Label(title_frame, text="🧬 Genome Search Dashboard", 
                             font=('Arial', 18, 'bold'), fg='#00ff88', bg='#2b2b2b')
        title_label.pack(pady=10)
        
        # Search info frame
        info_frame = tk.Frame(self.root, bg='#3b3b3b', height=100)
        info_frame.pack(fill='x', padx=10, pady=5)
        info_frame.pack_propagate(False)
        
        self.pattern_label = tk.Label(info_frame, text="Pattern: Loading...", 
                                    font=('Arial', 12), fg='#ffffff', bg='#3b3b3b')
        self.pattern_label.pack(anchor='w', padx=10, pady=2)
        
        self.algorithm_label = tk.Label(info_frame, text="Algorithm: Loading...", 
                                      font=('Arial', 12), fg='#ffffff', bg='#3b3b3b')
        self.algorithm_label.pack(anchor='w', padx=10, pady=2)
        
        self.chromosomes_label = tk.Label(info_frame, text="Chromosomes: Loading...", 
                                        font=('Arial', 12), fg='#ffffff', bg='#3b3b3b')
        self.chromosomes_label.pack(anchor='w', padx=10, pady=2)
        
        # Overall progress frame
        progress_frame = tk.Frame(self.root, bg='#3b3b3b', height=80)
        progress_frame.pack(fill='x', padx=10, pady=5)
        progress_frame.pack_propagate(False)
        
        self.overall_label = tk.Label(progress_frame, text="Overall Progress: 0/0 (0%)", 
                                    font=('Arial', 14, 'bold'), fg='#00ff88', bg='#3b3b3b')
        self.overall_label.pack(pady=5)
        
        self.overall_progress_bar = ttk.Progressbar(progress_frame, length=800, mode='determinate')
        self.overall_progress_bar.pack(pady=5)
        
        self.matches_label = tk.Label(progress_frame, text="Total Matches: 0", 
                                    font=('Arial', 12), fg='#ffaa00', bg='#3b3b3b')
        self.matches_label.pack()
        
        # Worker status frame with scrollable content
        worker_frame = tk.Frame(self.root, bg='#2b2b2b')
        worker_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        # Create scrollable canvas
        canvas = tk.Canvas(worker_frame, bg='#2b2b2b', highlightthickness=0)
        scrollbar = ttk.Scrollbar(worker_frame, orient="vertical", command=canvas.yview)
        self.scrollable_frame = tk.Frame(canvas, bg='#2b2b2b')
        
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Worker status widgets
        self.worker_widgets = {}
        self._create_worker_widgets()
    
    def _create_worker_widgets(self):
        """Create widgets for each worker."""
        for worker_id in range(self.worker_count):
            worker_container = tk.Frame(self.scrollable_frame, bg='#4b4b4b', relief='raised', bd=1)
            worker_container.pack(fill='x', padx=5, pady=2)
            
            # Worker header
            header_frame = tk.Frame(worker_container, bg='#4b4b4b')
            header_frame.pack(fill='x', padx=5, pady=2)
            
            worker_label = tk.Label(header_frame, text=f"Worker {worker_id + 1}", 
                                  font=('Arial', 12, 'bold'), fg='#ffffff', bg='#4b4b4b')
            worker_label.pack(side='left')
            
            status_label = tk.Label(header_frame, text="Waiting...", 
                                  font=('Arial', 10), fg='#888888', bg='#4b4b4b')
            status_label.pack(side='right')
            
            # Progress bar
            progress_bar = ttk.Progressbar(worker_container, length=400, mode='determinate')
            progress_bar.pack(fill='x', padx=5, pady=2)
            
            # Details label
            details_label = tk.Label(worker_container, text="No activity", 
                                   font=('Arial', 9), fg='#cccccc', bg='#4b4b4b')
            details_label.pack(fill='x', padx=5, pady=2)
            
            # Store widget references
            self.worker_widgets[worker_id] = {
                'container': worker_container,
                'status_label': status_label,
                'progress_bar': progress_bar,
                'details_label': details_label
            }
    
    def _message_processor(self):
        """Process messages from the command queue in a separate thread."""
        while self.running:
            try:
                # Get message with timeout
                message_type, data = self.command_queue.get(timeout=0.1)
                
                if message_type == "CLOSE":
                    self.running = False
                    break
                elif message_type == "UPDATE_WORKER":
                    self._update_worker_data(data)
                elif message_type == "UPDATE_OVERALL":
                    self._update_overall_progress(data)
                elif message_type == "SEARCH_INFO":
                    self._update_search_info(data)
                    
            except Exception:
                continue  # Continue if queue is empty or other errors
    
    def _update_worker_data(self, data: Dict[str, Any]):
        """Update worker data from message."""
        if isinstance(data, dict):
            worker_id = data.get('worker_id')
            if worker_id is not None:
                try:
                    worker_id = int(worker_id) if isinstance(worker_id, str) else worker_id
                    if 0 <= worker_id < self.worker_count:
                        self.worker_data[worker_id] = data
                except (ValueError, TypeError):
                    pass
    
    def _update_overall_progress(self, data: Dict[str, Any]):
        """Update overall progress data."""
        if isinstance(data, dict):
            self.overall_progress.update(data)
    
    def _update_search_info(self, data: Dict[str, Any]):
        """Update search information."""
        if isinstance(data, dict):
            self.search_info.update(data)
    
    def _schedule_gui_refresh(self):
        """Schedule GUI refresh at regular intervals."""
        self._refresh_gui()
        if self.running:
            self.root.after(500, self._schedule_gui_refresh)  # Refresh every 500ms
    
    def _refresh_gui(self):
        """Refresh the GUI with current data."""
        try:
            # Update search info labels
            if self.search_info.get('pattern'):
                self.pattern_label.config(text=f"Pattern: {self.search_info['pattern']}")
            
            if self.search_info.get('algorithm'):
                self.algorithm_label.config(text=f"Algorithm: {self.search_info['algorithm']}")
            
            if self.search_info.get('chromosomes'):
                chrom_list = self.search_info['chromosomes']
                chrom_text = f"Chromosomes: {len(chrom_list)} total"
                if len(chrom_list) <= 5:
                    chrom_text += f" ({', '.join(chrom_list)})"
                self.chromosomes_label.config(text=chrom_text)
            
            # Update overall progress
            completed = self.overall_progress.get('completed', 0)
            total = self.overall_progress.get('total', 0)
            matches = self.overall_progress.get('matches', 0)
            percentage = self.overall_progress.get('percentage', 0)
            
            self.overall_label.config(text=f"Overall Progress: {completed}/{total} ({percentage:.1f}%)")
            self.overall_progress_bar['value'] = percentage
            self.matches_label.config(text=f"Total Matches: {matches}")
            
            # Update worker widgets
            for worker_id, widgets in self.worker_widgets.items():
                worker_data = self.worker_data.get(worker_id, {})
                
                # Update status
                status = worker_data.get('status', 'Waiting')
                widgets['status_label'].config(text=status)
                
                # Update progress bar
                progress = worker_data.get('progress', 0)
                widgets['progress_bar']['value'] = progress
                
                # Update details
                details = worker_data.get('details', 'No activity')
                widgets['details_label'].config(text=details[:60] + "..." if len(details) > 60 else details)
                
                # Color coding based on status
                if status == 'Completed':
                    widgets['status_label'].config(fg='#00ff88')
                    widgets['container'].config(bg='#2d5a2d')
                elif status == 'Failed':
                    widgets['status_label'].config(fg='#ff4444')
                    widgets['container'].config(bg='#5a2d2d')
                elif status == 'Processing':
                    widgets['status_label'].config(fg='#ffaa00')
                    widgets['container'].config(bg='#5a4d2d')
                else:
                    widgets['status_label'].config(fg='#888888')
                    widgets['container'].config(bg='#4b4b4b')
                    
        except Exception as e:
            print(f"[DEBUG] GUI refresh error: {e}")
    
    def run(self):
        """Start the dashboard main loop."""
        try:
            self.root.mainloop()
        finally:
            self.running = False

def start_dashboard(command_queue: multiprocessing.Queue, worker_count: int = 8):
    """Start the enhanced dashboard with proper error handling."""
    try:
        dashboard = EnhancedWorkerDashboard(command_queue, worker_count)
        dashboard.run()
    except Exception as e:
        print(f"[ERROR] Dashboard failed to start: {e}")

# Keep backward compatibility
class DashboardApp:
    """Legacy dashboard class for backward compatibility."""
    
    def __init__(self, root: tk.Tk, command_queue: multiprocessing.Queue, worker_count: int):
        self.root = root
        self.dashboard = EnhancedWorkerDashboard(command_queue, worker_count)
        self.dashboard.root = root  # Use the provided root
        self.dashboard._create_widgets()  # Recreate widgets in the provided root
        
        # Start the message processor
        self.dashboard.update_thread = threading.Thread(target=self.dashboard._message_processor, daemon=True)
        self.dashboard.update_thread.start()
        self.dashboard._schedule_gui_refresh()
    
    def process_queue(self):
        """Legacy method - no longer needed."""
        pass
