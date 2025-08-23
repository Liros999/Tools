#!/usr/bin/env python3
"""
Worker Dashboard for Genome Analyzer
Runs as a separate process to avoid tkinter threading issues on Windows.
"""

import tkinter as tk
from tkinter import ttk
import queue
import sys

class DashboardApp:
    """A tkinter-based dashboard to monitor worker progress."""
    
    def __init__(self, root: tk.Tk, command_queue: queue.Queue, worker_count: int):
        self.root = root
        self.command_queue = command_queue
        self.worker_count = worker_count
        
        # Configure the main window
        self.root.title("🧬 Genome Analyzer - Worker Dashboard")
        self.root.geometry("900x600")
        self.root.configure(bg='#2b2b2b')
        
        # Create the UI
        self.create_ui()
        
        # Start processing messages
        self.process_queue()
    
    def create_ui(self):
        """Create the dashboard UI components."""
        # Title
        title_label = tk.Label(
            self.root, 
            text="🧬 GENOME ANALYZER WORKER DASHBOARD", 
            font=('Arial', 16, 'bold'),
            bg='#2b2b2b', fg='#ffffff'
        )
        title_label.pack(pady=20)
        
        # Overall Progress Section
        overall_frame = tk.Frame(self.root, bg='#2b2b2b')
        overall_frame.pack(pady=10, padx=20, fill='x')
        
        tk.Label(
            overall_frame, 
            text="Overall Progress:", 
            font=('Arial', 12, 'bold'),
            bg='#2b2b2b', fg='#ffffff'
        ).pack(anchor='w')
        
        self.overall_progress = ttk.Progressbar(
            overall_frame, 
            length=700, 
            mode='determinate',
            style='Custom.Horizontal.TProgressbar'
        )
        self.overall_progress.pack(pady=5, fill='x')
        
        self.overall_label = tk.Label(
            overall_frame, 
            text="0 / 0 chromosomes completed (0%)",
            font=('Arial', 10),
            bg='#2b2b2b', fg='#cccccc'
        )
        self.overall_label.pack(anchor='w')
        
        # Results Counter
        results_frame = tk.Frame(self.root, bg='#2b2b2b')
        results_frame.pack(pady=10, padx=20, fill='x')
        
        tk.Label(
            results_frame, 
            text="Results Found:", 
            font=('Arial', 12, 'bold'),
            bg='#2b2b2b', fg='#ffffff'
        ).pack(anchor='w')
        
        self.results_label = tk.Label(
            results_frame, 
            text="0 results found",
            font=('Arial', 14, 'bold'),
            bg='#2b2b2b', fg='#4CAF50'
        )
        self.results_label.pack(anchor='w')
        
        # Worker Details Section
        worker_frame = tk.Frame(self.root, bg='#2b2b2b')
        worker_frame.pack(pady=20, padx=20, fill='both', expand=True)
        
        tk.Label(
            worker_frame, 
            text="Worker Status:", 
            font=('Arial', 12, 'bold'),
            bg='#2b2b2b', fg='#ffffff'
        ).pack(anchor='w')
        
        # Create worker widgets
        self.worker_widgets = []
        for i in range(self.worker_count):
            worker_container = tk.Frame(worker_frame, bg='#2b2b2b')
            worker_container.pack(fill='x', pady=4)
            
            # Worker label
            label = tk.Label(
                worker_container, 
                text=f"Worker {i+1}: [Idle]", 
                width=60, 
                anchor="w",
                font=('Consolas', 9),
                bg='#404040', fg='#ffffff',
                relief='solid', bd=1
            )
            label.pack(side='left', padx=5, pady=2)
            
            # Progress bar
            progress = ttk.Progressbar(
                worker_container, 
                length=200, 
                mode='determinate'
            )
            progress.pack(side='left', padx=5, pady=2, fill='x', expand=True)
            
            self.worker_widgets.append({
                'label': label, 
                'progress': progress
            })
        
        # Status messages area
        status_frame = tk.Frame(self.root, bg='#2b2b2b')
        status_frame.pack(pady=10, padx=20, fill='x')
        
        tk.Label(
            status_frame, 
            text="Status Messages:", 
            font=('Arial', 12, 'bold'),
            bg='#2b2b2b', fg='#ffffff'
        ).pack(anchor='w')
        
        self.status_text = tk.Text(
            status_frame, 
            height=8, 
            width=80,
            bg='#404040', 
            fg='#ffffff',
            font=('Consolas', 9)
        )
        self.status_text.pack(pady=5, fill='x')
        
        # Scrollbar for status text
        scrollbar = tk.Scrollbar(self.status_text)
        scrollbar.pack(side='right', fill='y')
        self.status_text.config(yscrollcommand=scrollbar.set)
        scrollbar.config(command=self.status_text.yview)
    
    def process_queue(self):
        """Check the command queue for new messages and update the GUI."""
        try:
            while True:
                command, data = self.command_queue.get_nowait()
                
                if command == "UPDATE_WORKER":
                    worker_id, status, progress, task = data
                    if 0 <= worker_id < len(self.worker_widgets):
                        widget = self.worker_widgets[worker_id]
                        widget['label'].config(text=f"Worker {worker_id+1}: [{status}] - {task}")
                        widget['progress']['value'] = progress
                        
                        # Color coding based on status
                        if status == "Completed":
                            widget['label'].config(bg='#4CAF50', fg='#ffffff')
                        elif status == "Processing":
                            widget['label'].config(bg='#2196F3', fg='#ffffff')
                        elif status == "Error":
                            widget['label'].config(bg='#f44336', fg='#ffffff')
                        else:
                            widget['label'].config(bg='#404040', fg='#ffffff')
                
                elif command == "UPDATE_OVERALL":
                    completed, total, found = data
                    if total > 0:
                        percentage = (completed / total) * 100
                        self.overall_progress['value'] = percentage
                        self.overall_label.config(
                            text=f"{completed} / {total} chromosomes completed ({percentage:.1f}%)"
                        )
                    
                    self.results_label.config(text=f"{found:,} results found")
                
                elif command == "STATUS_MESSAGE":
                    message = data
                    self.status_text.insert(tk.END, f"{message}\n")
                    self.status_text.see(tk.END)
                
                elif command == "CLOSE":
                    self.root.quit()
                    return
                    
        except queue.Empty:
            pass  # No new commands
        
        finally:
            # Schedule next update
            self.root.after(100, self.process_queue)
    
    def add_status_message(self, message: str):
        """Add a status message to the display."""
        self.status_text.insert(tk.END, f"{message}\n")
        self.status_text.see(tk.END)

def start_dashboard(command_queue: queue.Queue, worker_count: int):
    """Entry point for the subprocess to run the dashboard."""
    try:
        # Create the main window
        root = tk.Tk()
        
        # Create the dashboard app
        dashboard = DashboardApp(root, command_queue, worker_count)
        
        # Add initial status message
        dashboard.add_status_message(f"[INFO] Dashboard started with {worker_count} workers")
        dashboard.add_status_message("[INFO] Waiting for search to begin...")
        
        # Start the main loop
        root.mainloop()
        
    except Exception as e:
        print(f"Dashboard Error: {e}")
        # Try to send error back to main process
        try:
            command_queue.put(("STATUS_MESSAGE", f"[ERROR] Dashboard failed: {e}"))
        except:
            pass

if __name__ == "__main__":
    # Test the dashboard standalone
    print("Dashboard module loaded successfully!")
    print("This module is designed to be imported and run as a subprocess.")
