"""
Dashboard for genome analysis.
Provides real-time monitoring of search progress and worker status.
"""

import tkinter as tk
from tkinter import ttk, scrolledtext
import multiprocessing
import threading
import time
from typing import Optional, Dict, Any


def start_dashboard(command_queue: multiprocessing.Queue):
    """
    Start the worker dashboard for monitoring search progress.
    
    Args:
        command_queue: Queue for receiving commands from workers
    """
    try:
        # Create and run dashboard
        dashboard = WorkerDashboard(command_queue)
        dashboard.run()
        
    except Exception as e:
        print(f"[ERROR] Dashboard failed to start: {e}")


class WorkerDashboard:
    """Real-time dashboard for monitoring genome search progress."""
    
    def __init__(self, command_queue: multiprocessing.Queue):
        self.command_queue = command_queue
        self.root = None
        self.worker_frames = {}
        self.overall_progress = {'completed': 0, 'total': 0, 'matches': 0}
        self.search_info = {}
        self.sequence_info = {}
        
    def run(self):
        """Run the dashboard main loop."""
        try:
            # Create main window
            self.root = tk.Tk()
            self.root.title("Genome Analyzer - Worker Dashboard")
            self.root.geometry("1200x800")
            
            # Setup UI components
            self._setup_ui()
            
            # Start message processing thread
            self._start_message_processor()
            
            # Run main loop
            self.root.mainloop()
            
        except Exception as e:
            print(f"[ERROR] Dashboard run failed: {e}")
    
    def _setup_ui(self):
        """Setup the dashboard user interface."""
        try:
            # Main frame
            main_frame = ttk.Frame(self.root, padding="10")
            main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
            
            # Configure grid weights
            self.root.columnconfigure(0, weight=1)
            self.root.rowconfigure(0, weight=1)
            main_frame.columnconfigure(1, weight=1)
            
            # Title
            title_label = ttk.Label(main_frame, text="Genome Analyzer - Worker Dashboard", 
                                   font=("Arial", 16, "bold"))
            title_label.grid(row=0, column=0, columnspan=2, pady=(0, 20))
            
            # Search info frame
            self._create_search_info_frame(main_frame)
            
            # Overall progress frame
            self._create_overall_progress_frame(main_frame)
            
            # Worker status frame
            self._create_worker_status_frame(main_frame)
            
            # Results log frame
            self._create_results_log_frame(main_frame)
            
        except Exception as e:
            print(f"[ERROR] UI setup failed: {e}")
    
    def _create_search_info_frame(self, parent):
        """Create search information display frame."""
        try:
            frame = ttk.LabelFrame(parent, text="Search Information", padding="10")
            frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
            
            # Pattern info
            self.pattern_label = ttk.Label(frame, text="Pattern: Not set")
            self.pattern_label.grid(row=0, column=0, sticky=tk.W, padx=(0, 20))
            
            # Chromosome info
            self.chromosome_label = ttk.Label(frame, text="Chromosomes: Not set")
            self.chromosome_label.grid(row=0, column=1, sticky=tk.W, padx=(0, 20))
            
            # Algorithm info
            self.algorithm_label = ttk.Label(frame, text="Algorithm: Not set")
            self.algorithm_label.grid(row=0, column=2, sticky=tk.W)
            
        except Exception as e:
            print(f"[ERROR] Search info frame creation failed: {e}")
    
    def _create_overall_progress_frame(self, parent):
        """Create overall progress display frame."""
        try:
            frame = ttk.LabelFrame(parent, text="Overall Progress", padding="10")
            frame.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
            
            # Progress bar
            self.progress_var = tk.DoubleVar()
            self.progress_bar = ttk.Progressbar(frame, variable=self.progress_var, maximum=100)
            self.progress_bar.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
            
            # Progress labels
            self.progress_label = ttk.Label(frame, text="0 / 0 completed (0%)")
            self.progress_label.grid(row=1, column=0, sticky=tk.W)
            
            self.matches_label = ttk.Label(frame, text="Matches: 0")
            self.matches_label.grid(row=1, column=1, sticky=tk.W, padx=(20, 0))
            
        except Exception as e:
            print(f"[ERROR] Overall progress frame creation failed: {e}")
    
    def _create_worker_status_frame(self, parent):
        """Create worker status display frame."""
        try:
            frame = ttk.LabelFrame(parent, text="Worker Status", padding="10")
            frame.grid(row=3, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
            
            # Worker grid
            self.worker_frame = frame
            
        except Exception as e:
            print(f"[ERROR] Worker status frame creation failed: {e}")
    
    def _create_results_log_frame(self, parent):
        """Create results log display frame."""
        try:
            frame = ttk.LabelFrame(parent, text="Results Log", padding="10")
            frame.grid(row=4, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(0, 10))
            
            # Configure grid weights for this frame
            parent.rowconfigure(4, weight=1)
            frame.rowconfigure(0, weight=1)
            frame.columnconfigure(0, weight=1)
            
            # Results text area
            self.results_text = scrolledtext.ScrolledText(frame, height=15, width=100)
            self.results_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
            
        except Exception as e:
            print(f"[ERROR] Results log frame creation failed: {e}")
    
    def _start_message_processor(self):
        """Start the message processing thread."""
        try:
            def process_messages():
                while True:
                    try:
                        # Check for messages with timeout
                        if self.command_queue.poll(timeout=0.1):
                            message_type, data = self.command_queue.get_nowait()
                            self._handle_message(message_type, data)
                    except Exception as e:
                        print(f"[ERROR] Message processing failed: {e}")
                        time.sleep(0.1)
            
            # Start message processor thread
            message_thread = threading.Thread(target=process_messages, daemon=True)
            message_thread.start()
            
        except Exception as e:
            print(f"[ERROR] Message processor start failed: {e}")
    
    def _handle_message(self, message_type: str, data: Any):
        """Handle incoming dashboard messages."""
        try:
            if message_type == "SEARCH_INFO":
                self._update_search_info(data)
            elif message_type == "UPDATE_OVERALL":
                self._update_overall_progress(data)
            elif message_type == "UPDATE_WORKER":
                self._update_worker_status(data)
            elif message_type == "NEW_RESULT":
                self._add_result(data)
            elif message_type == "SEQUENCE_UPDATE":
                self._update_sequence_info(data)
            else:
                # Log unknown message types
                self._log_message(f"Unknown message type: {message_type}")
                
        except Exception as e:
            print(f"[ERROR] Message handling failed: {e}")
    
    def _update_search_info(self, data: Dict):
        """Update search information display."""
        try:
            self.search_info.update(data)
            
            # Update UI elements
            if 'pattern' in data:
                self.pattern_label.config(text=f"Pattern: {data['pattern']}")
            if 'chromosomes' in data:
                chroms = data['chromosomes']
                if len(chroms) > 3:
                    chroms = chroms[:3] + ['...']
                self.chromosome_label.config(text=f"Chromosomes: {', '.join(chroms)}")
            if 'algorithm' in data:
                self.algorithm_label.config(text=f"Algorithm: {data['algorithm']}")
                
        except Exception as e:
            print(f"[ERROR] Search info update failed: {e}")
    
    def _update_overall_progress(self, data: Dict):
        """Update overall progress display."""
        try:
            self.overall_progress.update(data)
            
            # Update progress bar
            if 'percentage' in data:
                self.progress_var.set(data['percentage'])
            
            # Update progress label
            completed = self.overall_progress.get('completed', 0)
            total = self.overall_progress.get('total', 0)
            percentage = self.overall_progress.get('percentage', 0)
            self.progress_label.config(text=f"{completed} / {total} completed ({percentage:.1f}%)")
            
            # Update matches label
            matches = self.overall_progress.get('matches', 0)
            self.matches_label.config(text=f"Matches: {matches}")
            
        except Exception as e:
            print(f"[ERROR] Overall progress update failed: {e}")
    
    def _update_worker_status(self, data: Dict):
        """Update worker status display."""
        try:
            worker_id = data.get('worker_id', 0)
            
            # Create worker frame if it doesn't exist
            if worker_id not in self.worker_frames:
                self._create_worker_frame(worker_id)
            
            # Update worker status
            worker_frame = self.worker_frames[worker_id]
            
            if 'status' in data:
                worker_frame['status_label'].config(text=f"Status: {data['status']}")
            if 'progress' in data:
                worker_frame['progress_var'].set(data['progress'])
                worker_frame['progress_label'].config(text=f"{data['progress']}%")
            if 'details' in data:
                worker_frame['details_label'].config(text=f"Details: {data['details']}")
                
        except Exception as e:
            print(f"[ERROR] Worker status update failed: {e}")
    
    def _create_worker_frame(self, worker_id: int):
        """Create a frame for displaying worker status."""
        try:
            # Create worker frame
            worker_frame = ttk.Frame(self.worker_frame)
            worker_frame.grid(row=worker_id // 4, column=worker_id % 4, padx=5, pady=5, sticky=(tk.W, tk.E))
            
            # Worker ID label
            id_label = ttk.Label(worker_frame, text=f"Worker {worker_id}")
            id_label.grid(row=0, column=0, columnspan=2, sticky=tk.W)
            
            # Status label
            status_label = ttk.Label(worker_frame, text="Status: Idle")
            status_label.grid(row=1, column=0, columnspan=2, sticky=tk.W)
            
            # Progress bar
            progress_var = tk.DoubleVar()
            progress_bar = ttk.Progressbar(worker_frame, variable=progress_var, maximum=100, width=15)
            progress_bar.grid(row=2, column=0, sticky=tk.W, pady=(5, 0))
            
            # Progress label
            progress_label = ttk.Label(worker_frame, text="0%")
            progress_label.grid(row=2, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
            
            # Details label
            details_label = ttk.Label(worker_frame, text="Details: Ready")
            details_label.grid(row=3, column=0, columnspan=2, sticky=tk.W, pady=(5, 0))
            
            # Store references
            self.worker_frames[worker_id] = {
                'frame': worker_frame,
                'status_label': status_label,
                'progress_var': progress_var,
                'progress_label': progress_label,
                'details_label': details_label
            }
            
        except Exception as e:
            print(f"[ERROR] Worker frame creation failed: {e}")
    
    def _add_result(self, data: Any):
        """Add a new result to the log."""
        try:
            if isinstance(data, (tuple, list)) and len(data) >= 8:
                # Format result for display
                if len(data) == 9:  # With chromosome
                    chrom, start, end, strand, pattern, found, revcomp, mismatches, location = data
                    result_text = f"[{chrom}] {start}-{end} ({strand}) | {pattern} | {found} | {location}"
                else:  # Without chromosome
                    start, end, strand, pattern, found, revcomp, mismatches, location = data
                    result_text = f"{start}-{end} ({strand}) | {pattern} | {found} | {location}"
                
                # Add to results log
                self.results_text.insert(tk.END, result_text + "\n")
                self.results_text.see(tk.END)
                
                # Update matches count
                current_matches = self.overall_progress.get('matches', 0)
                self.overall_progress['matches'] = current_matches + 1
                self.matches_label.config(text=f"Matches: {self.overall_progress['matches']}")
                
        except Exception as e:
            print(f"[ERROR] Result addition failed: {e}")
    
    def _update_sequence_info(self, data: Dict):
        """Update sequence information display."""
        try:
            self.sequence_info.update(data)
            
            # Log sequence update
            if 'current_sequence' in data:
                self._log_message(f"Processing: {data['current_sequence']}")
            if 'total_sequences' in data:
                self._log_message(f"Total sequences: {data['total_sequences']}")
                
        except Exception as e:
            print(f"[ERROR] Sequence info update failed: {e}")
    
    def _log_message(self, message: str):
        """Add a log message to the results area."""
        try:
            timestamp = time.strftime("%H:%M:%S")
            log_text = f"[{timestamp}] {message}\n"
            self.results_text.insert(tk.END, log_text)
            self.results_text.see(tk.END)
            
        except Exception as e:
            print(f"[ERROR] Log message failed: {e}")


# Export function
__all__ = ['start_dashboard']
