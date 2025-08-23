#!/usr/bin/env python3
"""
ULTIMATE Genome Analyzer Dashboard
Professional, feature-rich dashboard for genome search results and analysis
Enhanced with comprehensive styling, filtering, and data management
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import multiprocessing
import threading
import time
import csv
import os
import queue
from typing import Dict, Any, Optional, List, Tuple
from datetime import datetime

class UltimateGenomeDashboard:
    """Professional genome analysis dashboard with comprehensive features."""
    
    def __init__(self, command_queue: multiprocessing.Queue, worker_count: int = 8):
        self.command_queue = command_queue
        self.worker_count = worker_count
        self.running = True
        
        # Thread safety for dashboard updates
        self.update_lock = threading.Lock()
        
        # Enhanced data storage
        self.search_results = []  # All accumulated results
        self.filtered_results = []  # Currently filtered/displayed results
        self.worker_data = {}
        self.overall_progress = {'completed': 0, 'total': 0, 'matches': 0, 'percentage': 0}
        self.search_info = {
            'pattern': '', 'chromosomes': [], 'algorithm': '', 
            'sequence_file': '', 'current_sequence': '', 'sequence_number': 0, 'total_sequences': 0
        }
        self.search_statistics = {
            'start_time': time.time(),
            'end_time': None,
            'duration': 0,
            'patterns_processed': 0,
            'total_matches': 0,
            'matches_by_chromosome': {},
            'matches_by_pattern': {}
        }
        
        # Create main window
        self.root = tk.Tk()
        self.root.title("🧬 ULTIMATE Genome Analyzer Dashboard")
        self.root.geometry("1400x900")
        self.root.configure(bg='#1e1e1e')
        self.root.minsize(1200, 700)
        
        # Configure professional styles
        self._configure_styles()
        
        # Create enhanced GUI
        self._create_main_interface()
        
        # Start message processor
        self.update_thread = threading.Thread(target=self._message_processor, daemon=True)
        self.update_thread.start()
        
        # Schedule regular updates
        self._schedule_updates()
    
    def _configure_styles(self):
        """Configure professional styling and themes."""
        self.style = ttk.Style()
        self.style.theme_use('clam')
        
        # Professional color scheme
        self.colors = {
            'bg_primary': '#1e1e1e',
            'bg_secondary': '#2d2d2d', 
            'bg_tertiary': '#3d3d3d',
            'accent_primary': '#00ff88',
            'accent_secondary': '#00ccff',
            'accent_warning': '#ffaa00',
            'accent_error': '#ff4444',
            'text_primary': '#ffffff',
            'text_secondary': '#cccccc',
            'text_muted': '#888888'
        }
        
        # Configure ttk styles
        self.style.configure('Header.TLabel', 
                           background=self.colors['bg_primary'], 
                           foreground=self.colors['accent_primary'],
                           font=('Segoe UI', 14, 'bold'))
        
        self.style.configure('Info.TLabel',
                           background=self.colors['bg_secondary'],
                           foreground=self.colors['text_primary'],
                           font=('Segoe UI', 10))
        
        self.style.configure('Results.Treeview',
                           background=self.colors['bg_tertiary'],
                           foreground=self.colors['text_primary'],
                           fieldbackground=self.colors['bg_tertiary'],
                           font=('Consolas', 9))
        
        # Custom progress bar styles for worker status
        self.style.configure('Success.Horizontal.TProgressbar',
                           troughcolor=self.colors['bg_tertiary'],
                           background=self.colors['accent_primary'])
        
        self.style.configure('Error.Horizontal.TProgressbar',
                           troughcolor=self.colors['bg_tertiary'],
                           background=self.colors['accent_error'])
        
        self.style.configure('Active.Horizontal.TProgressbar',
                           troughcolor=self.colors['bg_tertiary'],
                           background=self.colors['accent_warning'])
        
        self.style.configure('Results.Treeview.Heading',
                           background=self.colors['bg_secondary'],
                           foreground=self.colors['accent_primary'],
                           font=('Segoe UI', 10, 'bold'))
    
    def _create_main_interface(self):
        """Create the enhanced main dashboard interface."""
        
        # Main container with professional layout
        main_container = tk.Frame(self.root, bg=self.colors['bg_primary'])
        main_container.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Enhanced sections
        self._create_header_section(main_container)
        self._create_progress_section(main_container)
        self._create_results_section(main_container)
        self._create_status_section(main_container)
    
    def _create_header_section(self, parent):
        """Create enhanced header with search information and real-time clock."""
        header_frame = tk.Frame(parent, bg=self.colors['bg_secondary'], relief='raised', bd=1)
        header_frame.pack(fill='x', pady=(0, 5))
        
        # Title and clock
        title_frame = tk.Frame(header_frame, bg=self.colors['bg_secondary'])
        title_frame.pack(fill='x', padx=10, pady=5)
        
        title_label = tk.Label(title_frame, text="🧬 ULTIMATE Genome Analyzer Dashboard",
                             font=('Segoe UI', 16, 'bold'), 
                             fg=self.colors['accent_primary'], 
                             bg=self.colors['bg_secondary'])
        title_label.pack(side='left')
        
        # Real-time clock
        self.clock_label = tk.Label(title_frame, text="", 
                                  font=('Segoe UI', 10),
                                  fg=self.colors['text_secondary'],
                                  bg=self.colors['bg_secondary'])
        self.clock_label.pack(side='right')
        
        # Search information grid
        info_frame = tk.Frame(header_frame, bg=self.colors['bg_secondary'])
        info_frame.pack(fill='x', padx=10, pady=5)
        
        # Left column - Current search info
        left_info = tk.Frame(info_frame, bg=self.colors['bg_secondary'])
        left_info.pack(side='left', fill='both', expand=True)
        
        self.pattern_label = tk.Label(left_info, text="Pattern: Loading...",
                                    font=('Consolas', 11, 'bold'),
                                    fg=self.colors['text_primary'],
                                    bg=self.colors['bg_secondary'])
        self.pattern_label.pack(anchor='w')
        
        self.algorithm_label = tk.Label(left_info, text="Algorithm: Loading...",
                                      font=('Segoe UI', 10),
                                      fg=self.colors['text_secondary'],
                                      bg=self.colors['bg_secondary'])
        self.algorithm_label.pack(anchor='w')
        
        self.sequence_label = tk.Label(left_info, text="Sequence: Loading...",
                                     font=('Segoe UI', 10, 'bold'),
                                     fg=self.colors['accent_secondary'],
                                     bg=self.colors['bg_secondary'])
        self.sequence_label.pack(anchor='w')
        
        # Right column - File and batch info
        right_info = tk.Frame(info_frame, bg=self.colors['bg_secondary'])
        right_info.pack(side='right', fill='both', expand=True)
        
        self.file_label = tk.Label(right_info, text="File: No file loaded",
                                 font=('Segoe UI', 10),
                                 fg=self.colors['text_secondary'],
                                 bg=self.colors['bg_secondary'])
        self.file_label.pack(anchor='e')
        
        self.batch_progress_label = tk.Label(right_info, text="Batch: 0/0 patterns",
                                           font=('Segoe UI', 10, 'bold'),
                                           fg=self.colors['accent_warning'],
                                           bg=self.colors['bg_secondary'])
        self.batch_progress_label.pack(anchor='e')
        
        self.chromosomes_label = tk.Label(right_info, text="Chromosomes: Loading...",
                                        font=('Segoe UI', 10),
                                        fg=self.colors['text_secondary'],
                                        bg=self.colors['bg_secondary'])
        self.chromosomes_label.pack(anchor='e')
    
    def _create_progress_section(self, parent):
        """Create enhanced progress monitoring section."""
        progress_frame = tk.Frame(parent, bg=self.colors['bg_secondary'], relief='raised', bd=1)
        progress_frame.pack(fill='x', pady=(0, 5))
        
        # Overall progress
        overall_frame = tk.Frame(progress_frame, bg=self.colors['bg_secondary'])
        overall_frame.pack(fill='x', padx=10, pady=5)
        
        tk.Label(overall_frame, text="Overall Progress:",
               font=('Segoe UI', 12, 'bold'),
               fg=self.colors['accent_primary'],
               bg=self.colors['bg_secondary']).pack(side='left')
        
        self.overall_label = tk.Label(overall_frame, text="0/0 chromosomes (0%)",
                                    font=('Segoe UI', 12),
                                    fg=self.colors['text_primary'],
                                    bg=self.colors['bg_secondary'])
        self.overall_label.pack(side='left', padx=(10, 0))
        
        self.matches_label = tk.Label(overall_frame, text="📊 0 matches found",
                                    font=('Segoe UI', 12, 'bold'),
                                    fg=self.colors['accent_secondary'],
                                    bg=self.colors['bg_secondary'])
        self.matches_label.pack(side='right')
        
        # Enhanced progress bar
        self.overall_progress_bar = ttk.Progressbar(progress_frame, length=1000, mode='determinate')
        self.overall_progress_bar.pack(fill='x', padx=10, pady=(0, 10))
        
        # Enhanced Worker Status Panel
        worker_frame = tk.Frame(progress_frame, bg=self.colors['bg_secondary'])
        worker_frame.pack(fill='x', padx=10, pady=(0, 10))
        
        # Worker Status Header
        tk.Label(worker_frame, text="🔧 WORKER STATUS PANEL",
               font=('Segoe UI', 12, 'bold'),
               fg=self.colors['accent_primary'],
               bg=self.colors['bg_secondary']).pack(side='left')
        
        # Create detailed worker status table
        self.worker_status_frame = tk.Frame(worker_frame, bg=self.colors['bg_secondary'])
        self.worker_status_frame.pack(side='left', padx=(20, 0), fill='x', expand=True)
        
        # Worker Status Table Headers
        headers_frame = tk.Frame(self.worker_status_frame, bg=self.colors['bg_tertiary'])
        headers_frame.pack(fill='x', pady=(0, 2))
        
        header_labels = ['Worker ID', 'Chromosome', 'Progress', 'Speed', 'Status', 'Memory']
        header_widths = [8, 12, 10, 10, 10, 8]
        
        for i, (header, width) in enumerate(zip(header_labels, header_widths)):
            tk.Label(headers_frame, text=header, width=width, height=1,
                   font=('Segoe UI', 9, 'bold'),
                   fg=self.colors['accent_primary'],
                   bg=self.colors['bg_tertiary'],
                   relief='flat').pack(side='left', padx=1)
        
        # Worker Status Rows
        self.worker_rows = {}
        for i in range(self.worker_count):
            row_frame = tk.Frame(self.worker_status_frame, bg=self.colors['bg_secondary'])
            row_frame.pack(fill='x', pady=1)
            
            # Worker ID
            worker_id_label = tk.Label(row_frame, text=f"W-{i+1:03d}", width=8, height=1,
                                     font=('Consolas', 9),
                                     fg=self.colors['text_primary'],
                                     bg=self.colors['bg_secondary'],
                                     relief='flat')
            worker_id_label.pack(side='left', padx=1)
            
            # Chromosome
            chr_label = tk.Label(row_frame, text="--", width=12, height=1,
                               font=('Consolas', 9),
                               fg=self.colors['text_muted'],
                               bg=self.colors['bg_secondary'],
                               relief='flat')
            chr_label.pack(side='left', padx=1)
            
            # Progress Bar
            progress_frame = tk.Frame(row_frame, bg=self.colors['bg_secondary'])
            progress_frame.pack(side='left', padx=1)
            progress_bar = ttk.Progressbar(progress_frame, length=80, mode='determinate')
            progress_bar.pack(side='left')
            
            # Speed
            speed_label = tk.Label(row_frame, text="--", width=10, height=1,
                                 font=('Consolas', 9),
                                 fg=self.colors['text_muted'],
                                 bg=self.colors['bg_secondary'],
                                 relief='flat')
            speed_label.pack(side='left', padx=1)
            
            # Status with color coding
            status_label = tk.Label(row_frame, text="Idle", width=10, height=1,
                                  font=('Consolas', 9, 'bold'),
                                  fg=self.colors['text_primary'],
                                  bg=self.colors['bg_tertiary'],
                                  relief='flat')
            status_label.pack(side='left', padx=1)
            
            # Memory
            memory_label = tk.Label(row_frame, text="--", width=8, height=1,
                                  font=('Consolas', 9),
                                  fg=self.colors['text_muted'],
                                  bg=self.colors['bg_secondary'],
                                  relief='flat')
            memory_label.pack(side='left', padx=1)
            
            # Store references for updates
            self.worker_rows[i] = {
                'id': worker_id_label,
                'chromosome': chr_label,
                'progress': progress_bar,
                'speed': speed_label,
                'status': status_label,
                'memory': memory_label
            }
    
    def _create_results_section(self, parent):
        """Create the enhanced results display section with filtering and controls."""
        results_frame = tk.Frame(parent, bg=self.colors['bg_primary'])
        results_frame.pack(fill='both', expand=True, pady=(0, 5))
        
        # Enhanced results controls toolbar
        controls_frame = tk.Frame(results_frame, bg=self.colors['bg_secondary'], relief='raised', bd=1)
        controls_frame.pack(fill='x', pady=(0, 5))
        
        # Left side - Search and filter controls
        left_controls = tk.Frame(controls_frame, bg=self.colors['bg_secondary'])
        left_controls.pack(side='left', padx=10, pady=5)
        
        tk.Label(left_controls, text="🔍 Filter:",
               font=('Segoe UI', 10, 'bold'),
               fg=self.colors['text_primary'],
               bg=self.colors['bg_secondary']).pack(side='left')
        
        self.filter_var = tk.StringVar()
        self.filter_entry = tk.Entry(left_controls, textvariable=self.filter_var,
                                   font=('Segoe UI', 10), width=20)
        self.filter_entry.pack(side='left', padx=(5, 10))
        self.filter_entry.bind('<KeyRelease>', self._on_filter_change)
        
        # Filter by chromosome dropdown
        tk.Label(left_controls, text="Chr:",
               font=('Segoe UI', 10),
               fg=self.colors['text_primary'],
               bg=self.colors['bg_secondary']).pack(side='left')
        
        self.chromosome_filter = ttk.Combobox(left_controls, width=8, state='readonly')
        self.chromosome_filter.pack(side='left', padx=(5, 10))
        self.chromosome_filter.bind('<<ComboboxSelected>>', self._on_filter_change)
        
        # Start point filter
        tk.Label(left_controls, text="Start:",
               font=('Segoe UI', 10),
               fg=self.colors['text_primary'],
               bg=self.colors['bg_secondary']).pack(side='left')
        
        self.start_filter = tk.Entry(left_controls, font=('Segoe UI', 9), width=8)
        self.start_filter.pack(side='left', padx=(5, 10))
        self.start_filter.bind('<KeyRelease>', self._on_filter_change)
        
        # End point filter
        tk.Label(left_controls, text="End:",
               font=('Segoe UI', 10),
               fg=self.colors['text_primary'],
               bg=self.colors['bg_secondary']).pack(side='left')
        
        self.end_filter = tk.Entry(left_controls, font=('Segoe UI', 9), width=8)
        self.end_filter.pack(side='left', padx=(5, 10))
        self.end_filter.bind('<KeyRelease>', self._on_filter_change)
        
        # Strand filter
        tk.Label(left_controls, text="Strand:",
               font=('Segoe UI', 10),
               fg=self.colors['text_primary'],
               bg=self.colors['bg_secondary']).pack(side='left')
        
        self.strand_filter = ttk.Combobox(left_controls, width=6, state='readonly',
                                         values=['All', '+', '-'])
        self.strand_filter.set('All')
        self.strand_filter.pack(side='left', padx=(5, 10))
        self.strand_filter.bind('<<ComboboxSelected>>', self._on_filter_change)
        
        # Mismatches filter
        tk.Label(left_controls, text="Mis:",
               font=('Segoe UI', 10),
               fg=self.colors['text_primary'],
               bg=self.colors['bg_secondary']).pack(side='left')
        
        self.mismatches_filter = ttk.Combobox(left_controls, width=6, state='readonly',
                                             values=['All', '0', '1', '2', '3', '4', '5+'])
        self.mismatches_filter.set('All')
        self.mismatches_filter.pack(side='left', padx=(5, 10))
        self.mismatches_filter.bind('<<ComboboxSelected>>', self._on_filter_change)
        
        # Right side - Export and utility controls
        right_controls = tk.Frame(controls_frame, bg=self.colors['bg_secondary'])
        right_controls.pack(side='right', padx=10, pady=5)
        
        self.export_btn = tk.Button(right_controls, text="📄 Export Results",
                                  command=self._export_results,
                                  font=('Segoe UI', 9),
                                  bg=self.colors['bg_tertiary'],
                                  fg=self.colors['text_primary'])
        self.export_btn.pack(side='right', padx=5)
        
        self.clear_btn = tk.Button(right_controls, text="🗑️ Clear Results",
                                 command=self._clear_results,
                                 font=('Segoe UI', 9),
                                 bg=self.colors['accent_error'],
                                 fg=self.colors['text_primary'])
        self.clear_btn.pack(side='right', padx=5)
        
        self.clear_filters_btn = tk.Button(right_controls, text="🔍 Clear Filters",
                                         command=self._clear_filters,
                                         font=('Segoe UI', 9),
                                         bg=self.colors['accent_warning'],
                                         fg=self.colors['text_primary'])
        self.clear_filters_btn.pack(side='right', padx=5)
        
        self.results_count_label = tk.Label(right_controls, text="Showing: 0 / 0 results",
                                          font=('Segoe UI', 10, 'bold'),
                                          fg=self.colors['accent_primary'],
                                          bg=self.colors['bg_secondary'])
        self.results_count_label.pack(side='right', padx=10)
        
        # Enhanced results table with scrollbars
        table_container = tk.Frame(results_frame, bg=self.colors['bg_primary'])
        table_container.pack(fill='both', expand=True)
        
        # Create enhanced Treeview for results
        self.results_tree = ttk.Treeview(table_container, style='Results.Treeview')
        
        # Define enhanced columns
        columns = ('seq_num', 'pattern', 'chr', 'start', 'end', 'strand', 
                  'found_sequence', 'rev_comp', 'mismatches', 'location', 'timestamp')
        self.results_tree['columns'] = columns
        self.results_tree['show'] = 'headings'
        
        # Configure column headers and widths
        column_config = {
            'seq_num': ('Seq#', 50, 'center'),
            'pattern': ('Pattern', 120, 'center'),
            'chr': ('Chr', 60, 'center'),
            'start': ('Start', 80, 'center'),
            'end': ('End', 80, 'center'),
            'strand': ('Strand', 50, 'center'),
            'found_sequence': ('Found Sequence', 120, 'center'),
            'rev_comp': ('Rev Complement', 120, 'center'),
            'mismatches': ('Mis', 40, 'center'),
            'location': ('Gene/Location', 150, 'w'),
            'timestamp': ('Time', 80, 'center')
        }
        
        for col, (heading, width, anchor) in column_config.items():
            self.results_tree.heading(col, text=heading, anchor='center')
            self.results_tree.column(col, width=width, anchor=anchor, minwidth=30)
        
        # Add scrollbars
        v_scrollbar = ttk.Scrollbar(table_container, orient='vertical', command=self.results_tree.yview)
        h_scrollbar = ttk.Scrollbar(table_container, orient='horizontal', command=self.results_tree.xview)
        self.results_tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
        
        # Pack table and scrollbars
        self.results_tree.grid(row=0, column=0, sticky='nsew')
        v_scrollbar.grid(row=0, column=1, sticky='ns')
        h_scrollbar.grid(row=1, column=0, sticky='ew')
        
        table_container.grid_rowconfigure(0, weight=1)
        table_container.grid_columnconfigure(0, weight=1)
        
        # Bind double-click for detailed view
        self.results_tree.bind('<Double-1>', self._show_result_details)
        
        # Add enhanced context menu
        self._create_context_menu()
    
    def _create_status_section(self, parent):
        """Create enhanced bottom status section with performance metrics."""
        status_frame = tk.Frame(parent, bg=self.colors['bg_secondary'], relief='sunken', bd=1)
        status_frame.pack(fill='x')
        
        # Statistics display
        stats_frame = tk.Frame(status_frame, bg=self.colors['bg_secondary'])
        stats_frame.pack(fill='x', padx=10, pady=5)
        
        # Left side - Performance stats
        left_stats = tk.Frame(stats_frame, bg=self.colors['bg_secondary'])
        left_stats.pack(side='left')
        
        self.duration_label = tk.Label(left_stats, text="⏱️ Duration: 00:00:00",
                                     font=('Segoe UI', 10),
                                     fg=self.colors['text_secondary'],
                                     bg=self.colors['bg_secondary'])
        self.duration_label.pack(side='left', padx=(0, 20))
        
        self.speed_label = tk.Label(left_stats, text="⚡ Speed: 0 matches/sec",
                                  font=('Segoe UI', 10),
                                  fg=self.colors['text_secondary'],
                                  bg=self.colors['bg_secondary'])
        self.speed_label.pack(side='left', padx=(0, 20))
        
        # Right side - Current status
        self.status_label = tk.Label(stats_frame, text="🔄 Ready for search...",
                                   font=('Segoe UI', 10),
                                   fg=self.colors['accent_primary'],
                                   bg=self.colors['bg_secondary'])
        self.status_label.pack(side='right')
    
    def _create_context_menu(self):
        """Create enhanced context menu for results table."""
        self.context_menu = tk.Menu(self.root, tearoff=0)
        self.context_menu.add_command(label="📋 Copy Row", command=self._copy_selected_row)
        self.context_menu.add_command(label="📊 Show Details", command=self._show_result_details)
        self.context_menu.add_separator()
        self.context_menu.add_command(label="🔍 Filter by Chromosome", command=self._filter_by_selected_chr)
        self.context_menu.add_command(label="🧬 Filter by Pattern", command=self._filter_by_selected_pattern)
        
        self.results_tree.bind('<Button-3>', self._show_context_menu)
    
    def _show_context_menu(self, event):
        """Show context menu at cursor position."""
        try:
            item = self.results_tree.selection()[0]
            self.context_menu.post(event.x_root, event.y_root)
        except IndexError:
            pass
    
    def _copy_selected_row(self):
        """Copy selected row to clipboard."""
        try:
            item = self.results_tree.selection()[0]
            values = self.results_tree.item(item)['values']
            row_text = '\t'.join(str(v) for v in values)
            self.root.clipboard_clear()
            self.root.clipboard_append(row_text)
            self.status_label.config(text="📋 Row copied to clipboard")
        except IndexError:
            pass
    
    def _show_result_details(self, event=None):
        """Show detailed view of selected result."""
        try:
            item = self.results_tree.selection()[0]
            values = self.results_tree.item(item)['values']
            
            # Create detail window
            detail_window = tk.Toplevel(self.root)
            detail_window.title("Result Details")
            detail_window.geometry("500x400")
            detail_window.configure(bg=self.colors['bg_primary'])
            
            # Display all details
            detail_text = tk.Text(detail_window, wrap='word', font=('Consolas', 10),
                                bg=self.colors['bg_tertiary'], fg=self.colors['text_primary'])
            detail_text.pack(fill='both', expand=True, padx=10, pady=10)
            
            # Format details
            details = f"""
GENOME SEARCH RESULT DETAILS
{'='*40}

Sequence Number: {values[0]}
Pattern Searched: {values[1]}
Chromosome: {values[2]}
Start Position: {values[3]:,}
End Position: {values[4]:,}
Strand: {values[5]}
Found Sequence: {values[6]}
Reverse Complement: {values[7]}
Mismatches: {values[8]}
Gene/Location: {values[9]}
Discovery Time: {values[10]}

COORDINATES
Position Range: {values[3]:,} - {values[4]:,}
Length: {int(values[4]) - int(values[3]) + 1 if str(values[4]).isdigit() and str(values[3]).isdigit() else 'N/A'} bp

SEQUENCE ANALYSIS
Pattern: {values[1]}
Found:   {values[6]}
Match Quality: {'Perfect' if values[8] == '0' else f'{values[8]} mismatches'}
"""
            detail_text.insert('1.0', details)
            detail_text.config(state='disabled')
            
        except (IndexError, tk.TclError):
            pass
    
    def _filter_by_selected_chr(self):
        """Filter by selected chromosome."""
        try:
            item = self.results_tree.selection()[0]
            chr_value = self.results_tree.item(item)['values'][2]
            self.chromosome_filter.set(chr_value)
            self._apply_filters()
        except IndexError:
            pass
    
    def _filter_by_selected_pattern(self):
        """Filter by selected pattern."""
        try:
            item = self.results_tree.selection()[0]
            pattern_value = self.results_tree.item(item)['values'][1]
            self.filter_var.set(pattern_value)
            self._apply_filters()
        except IndexError:
            pass
    
    def _on_filter_change(self, event=None):
        """Handle filter changes."""
        self._apply_filters()
    
    def _apply_filters(self):
        """Apply current filters to results with enhanced filtering options."""
        filter_text = self.filter_var.get().lower()
        chr_filter = self.chromosome_filter.get()
        start_filter = self.start_filter.get().strip()
        end_filter = self.end_filter.get().strip()
        strand_filter = self.strand_filter.get()
        mismatches_filter = self.mismatches_filter.get()
        
        # Filter results
        self.filtered_results = []
        for result in self.search_results:
            # Apply text filter (searches pattern, sequence, location)
            if filter_text:
                searchable_text = f"{result.get('pattern', '')} {result.get('found_sequence', '')} {result.get('location', '')}".lower()
                if filter_text not in searchable_text:
                    continue
            
            # Apply chromosome filter
            if chr_filter and chr_filter != 'All':
                if result.get('chromosome', '') != chr_filter:
                    continue
            
            # Apply start point filter
            if start_filter:
                try:
                    start_value = int(start_filter)
                    result_start = int(result.get('start', 0))
                    if result_start < start_value:
                        continue
                except ValueError:
                    # If not a number, skip this filter
                    pass
            
            # Apply end point filter
            if end_filter:
                try:
                    end_value = int(end_filter)
                    result_end = int(result.get('end', 0))
                    if result_end > end_value:
                        continue
                except ValueError:
                    # If not a number, skip this filter
                    pass
            
            # Apply strand filter
            if strand_filter and strand_filter != 'All':
                if result.get('strand', '') != strand_filter:
                    continue
            
            # Apply mismatches filter
            if mismatches_filter and mismatches_filter != 'All':
                try:
                    result_mismatches = int(result.get('mismatches', 0))
                    if mismatches_filter == '5+':
                        if result_mismatches < 5:
                            continue
                    else:
                        filter_mismatches = int(mismatches_filter)
                        if result_mismatches != filter_mismatches:
                            continue
                except ValueError:
                    # If not a number, skip this filter
                    pass
                    
            self.filtered_results.append(result)
        
        # Update display
        self._update_results_display()
    
    def _update_results_display(self):
        """Update the results table display with enhanced formatting."""
        # Clear current display
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)
        
        # Add filtered results
        for result in self.filtered_results:
            values = (
                result.get('sequence_number', ''),
                result.get('pattern', ''),
                result.get('chromosome', ''),
                result.get('start', ''),
                result.get('end', ''),
                result.get('strand', ''),
                result.get('found_sequence', ''),
                result.get('rev_complement', ''),
                result.get('mismatches', ''),
                result.get('location', ''),
                result.get('timestamp', '')
            )
            
            # Add row with color coding based on mismatches
            item_id = self.results_tree.insert('', 'end', values=values)
            
            # Enhanced color coding
            mismatches = result.get('mismatches', 0)
            try:
                mismatch_count = int(mismatches)
                if mismatch_count == 0:
                    self.results_tree.set(item_id, 'mismatches', '✓')
                elif mismatch_count <= 2:
                    pass  # Default color
                else:
                    pass  # Could add tags for different colors
            except (ValueError, TypeError):
                pass
        
        # Update count label
        total_results = len(self.search_results)
        filtered_results = len(self.filtered_results)
        self.results_count_label.config(text=f"Showing: {filtered_results:,} / {total_results:,} results")
        
        # Update chromosome filter options
        chromosomes = sorted(set(r.get('chromosome', '') for r in self.search_results if r.get('chromosome', '')))
        self.chromosome_filter['values'] = ['All'] + chromosomes
        if not self.chromosome_filter.get():
            self.chromosome_filter.set('All')
    
    def _export_results(self):
        """Export results to CSV file with enhanced formatting."""
        if not self.filtered_results:
            messagebox.showwarning("No Data", "No results to export.")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension='.csv',
            filetypes=[('CSV files', '*.csv'), ('All files', '*.*')],
            title='Export Results'
        )
        
        if filename:
            try:
                with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
                    writer = csv.writer(csvfile)
                    
                    # Write header
                    header = ['Sequence_Number', 'Pattern', 'Chromosome', 'Start', 'End', 
                             'Strand', 'Found_Sequence', 'Reverse_Complement', 'Mismatches', 
                             'Location', 'Timestamp']
                    writer.writerow(header)
                    
                    # Write data
                    for result in self.filtered_results:
                        row = [
                            result.get('sequence_number', ''),
                            result.get('pattern', ''),
                            result.get('chromosome', ''),
                            result.get('start', ''),
                            result.get('end', ''),
                            result.get('strand', ''),
                            result.get('found_sequence', ''),
                            result.get('rev_complement', ''),
                            result.get('mismatches', ''),
                            result.get('location', ''),
                            result.get('timestamp', '')
                        ]
                        writer.writerow(row)
                
                self.status_label.config(text=f"📄 Exported {len(self.filtered_results)} results to {os.path.basename(filename)}")
                
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export results:\n{e}")
    
    def _clear_results(self):
        """Clear all results after confirmation."""
        if messagebox.askyesno("Clear Results", "Are you sure you want to clear all results?"):
            self.search_results.clear()
            self.filtered_results.clear()
            self._update_results_display()
            self.status_label.config(text="🗑️ Results cleared")
    
    def _clear_filters(self):
        """Clear all active filters and show all results."""
        # Reset filter values
        self.filter_var.set('')
        self.chromosome_filter.set('All')
        self.start_filter.delete(0, tk.END)
        self.end_filter.delete(0, tk.END)
        self.strand_filter.set('All')
        self.mismatches_filter.set('All')
        
        # Apply filters (which will show all results)
        self._apply_filters()
        
        print("[INFO] All filters cleared - showing all results")
    
    def _message_processor(self):
        """Process incoming messages from the search process."""
        while self.running:
            try:
                message_type, data = self.command_queue.get(timeout=0.1)
                
                if message_type == "CLOSE":
                    self.running = False
                    break
                elif message_type == "SEARCH_INFO":
                    self._update_search_info(data)
                elif message_type == "UPDATE_WORKER":
                    self._update_worker_status(data)
                elif message_type == "UPDATE_OVERALL":
                    self._update_overall_progress(data)
                elif message_type == "NEW_RESULT":
                    self._add_new_result(data)
                elif message_type == "SEQUENCE_UPDATE":
                    self._update_sequence_tracking(data)
                elif message_type == "RESULTS_BATCH":
                    self._add_results_batch(data)
                    
            except queue.Empty:
                continue
            except Exception as e:
                # Log dashboard message processing errors for debugging
                print(f"[DASHBOARD ERROR] Message processing failed: {e}")
                continue
    
    def _update_search_info(self, data):
        """Update search information."""
        self.search_info.update(data)
    
    def send_search_info(self, pattern: str, chromosomes: List[str], algorithm: str):
        """
        Send search information to dashboard - compatibility method for Genome_Analyzer.py calls.
        
        Args:
            pattern: Search pattern
            chromosomes: List of chromosomes to search
            algorithm: Search algorithm being used
        """
        # Update search info
        self.search_info.update({
            'pattern': pattern,
            'chromosomes': chromosomes,
            'algorithm': algorithm
        })
        
        # Update the search info display if it exists
        if hasattr(self, 'search_info_label'):
            self.search_info_label.config(
                text=f"Pattern: {pattern} | Algorithm: {algorithm} | Chromosomes: {len(chromosomes)}"
            )
    
    def _update_worker_status(self, data):
        """Update detailed worker status panel with real-time metrics."""
        with self.update_lock:  # Thread-safe updates
            worker_id = data.get('worker_id', 0)
            status = data.get('status', 'Idle')
            progress = data.get('progress', 0)
            details = data.get('details', '')
            chromosome = data.get('chromosome', '--')
            speed = data.get('speed', '--')
            memory = data.get('memory', '--')
            
            if 0 <= worker_id < len(self.worker_rows):
                row = self.worker_rows[worker_id]
                
                # Update chromosome
                row['chromosome'].config(text=chromosome)
                
                # Update progress bar
                row['progress']['value'] = progress
                
                # Update speed
                row['speed'].config(text=speed)
                
                # Update memory
                row['memory'].config(text=memory)
                
                # Update status with color coding
                status_colors = {
                    'Idle': (self.colors['bg_tertiary'], self.colors['text_muted']),
                    'Starting': (self.colors['accent_warning'], self.colors['bg_primary']),
                    'Processing': (self.colors['accent_primary'], self.colors['bg_primary']),
                    'Active': (self.colors['accent_primary'], self.colors['bg_primary']),
                    'Completed': (self.colors['accent_secondary'], self.colors['bg_primary']),
                    'Failed': (self.colors['accent_error'], self.colors['text_primary']),
                    'Waiting': (self.colors['bg_tertiary'], self.colors['text_primary'])
                }
                
                bg_color, fg_color = status_colors.get(status, status_colors['Idle'])
                row['status'].config(text=status, bg=bg_color, fg=fg_color)
                
                # Update progress bar color based on status
                if status == 'Completed':
                    row['progress']['style'] = 'Success.Horizontal.TProgressbar'
                elif status == 'Failed':
                    row['progress']['style'] = 'Error.Horizontal.TProgressbar'
                elif status in ['Processing', 'Active']:
                    row['progress']['style'] = 'Active.Horizontal.TProgressbar'
                else:
                    row['progress']['style'] = 'Horizontal.TProgressbar'
    
    def _update_overall_progress(self, data):
        """Update overall progress information."""
        self.overall_progress.update(data)
    
    def update_overall_progress(self, completed: int, total: int, matches: int):
        """
        Update overall progress - compatibility method for Genome_Analyzer.py calls.
        
        Args:
            completed: Number of completed chromosomes
            total: Total number of chromosomes
            matches: Total number of matches found
        """
        # Calculate percentage
        percentage = int((completed / total) * 100) if total > 0 else 0
        
        # Update progress data
        self.overall_progress.update({
            'completed': completed,
            'total': total,
            'matches': matches,
            'percentage': percentage
        })
        
        # Update the overall progress bar
        if hasattr(self, 'overall_progress_bar'):
            self.overall_progress_bar['value'] = percentage
        
        # Update the progress label
        if hasattr(self, 'overall_progress_label'):
            self.overall_progress_label.config(
                text=f"Progress: {completed}/{total} chromosomes ({percentage}%) - {matches} matches found"
            )
    
    def update_worker(self, worker_id: int, status: str, progress: int, details: str = ""):
        """
        Update worker status - compatibility method for Genome_Analyzer.py calls.
        
        Args:
            worker_id: Worker ID (0-based)
            status: Worker status (Idle, Starting, Processing, Active, Completed, Failed)
            progress: Progress percentage (0-100)
            details: Additional details
        """
        if 0 <= worker_id < len(self.worker_rows):
            # Extract chromosome from details if possible
            chromosome = '--'
            if ':' in details:
                chromosome = details.split(':')[0].strip()
            
            # Prepare data for the update method
            update_data = {
                'worker_id': worker_id,
                'status': status,
                'progress': progress,
                'chromosome': chromosome,
                'speed': '--',
                'memory': '--',
                'details': details
            }
            
            # Update the worker status
            self._update_worker_status(update_data)
    
    def update_worker_metrics(self, worker_id: int, **metrics):
        """
        Update detailed worker metrics in real-time.
        
        Args:
            worker_id: Worker ID (0-based)
            **metrics: Keyword arguments for various metrics
                - status: Worker status (Idle, Starting, Processing, Active, Completed, Failed)
                - progress: Progress percentage (0-100)
                - chromosome: Current chromosome being processed
                - speed: Processing speed (e.g., "2.3M/s")
                - memory: Memory usage (e.g., "145MB")
                - details: Additional details
        """
        if 0 <= worker_id < len(self.worker_rows):
            # Prepare data for the update method
            update_data = {
                'worker_id': worker_id,
                'status': metrics.get('status', 'Idle'),
                'progress': metrics.get('progress', 0),
                'chromosome': metrics.get('chromosome', '--'),
                'speed': metrics.get('speed', '--'),
                'memory': metrics.get('memory', '--'),
                'details': metrics.get('details', '')
            }
            
            # Update the worker status
            self._update_worker_status(update_data)
    
    def _add_new_result(self, result_data):
        """Add a new search result to the dashboard with enhanced formatting."""
        # Convert tuple to dict if needed
        if isinstance(result_data, (tuple, list)):
            # Handle the case where result_data is a tuple from the old format
            if len(result_data) >= 8:
                # Convert to dict format
                result_dict = {
                    'chromosome': result_data[0] if len(result_data) > 0 else 'N/A',
                    'start': result_data[1] if len(result_data) > 1 else '',
                    'end': result_data[2] if len(result_data) > 2 else '',
                    'strand': result_data[3] if len(result_data) > 3 else '',
                    'pattern': result_data[4] if len(result_data) > 4 else '',
                    'found_sequence': result_data[5] if len(result_data) > 5 else '',
                    'rev_complement': result_data[6] if len(result_data) > 6 else '',
                    'mismatches': result_data[7] if len(result_data) > 7 else 0,
                    'location': result_data[8] if len(result_data) > 8 else '',
                    'sequence_number': 0  # Will be set by caller
                }
                result_data = result_dict
            else:
                # Silent invalid result data format
                return
        
        # Ensure result_data is a dict
        if not isinstance(result_data, dict):
            # Silent result data processing error
            return
            
        # Add timestamp
        result_data['timestamp'] = datetime.now().strftime('%H:%M:%S')
        
        # Add to results
        self.search_results.append(result_data)
        
        # Update statistics
        self.search_statistics['total_matches'] += 1
        
        # Update chromosome statistics
        chromosome = result_data.get('chromosome', 'unknown')
        if chromosome not in self.search_statistics['matches_by_chromosome']:
            self.search_statistics['matches_by_chromosome'][chromosome] = 0
        self.search_statistics['matches_by_chromosome'][chromosome] += 1
        
        # Update pattern statistics
        pattern = result_data.get('pattern', 'unknown')
        if pattern not in self.search_statistics['matches_by_pattern']:
            self.search_statistics['matches_by_pattern'][pattern] = 0
        self.search_statistics['matches_by_pattern'][pattern] += 1
        
        # If no filters are active, add to filtered results
        if not self.filter_var.get() and (not self.chromosome_filter.get() or self.chromosome_filter.get() == 'All'):
            self.filtered_results.append(result_data)
            
            # Add to display
            values = (
                result_data.get('sequence_number', ''),
                result_data.get('pattern', ''),
                result_data.get('chromosome', ''),
                result_data.get('start', ''),
                result_data.get('end', ''),
                result_data.get('strand', ''),
                result_data.get('found_sequence', ''),
                result_data.get('rev_complement', ''),
                result_data.get('mismatches', ''),
                result_data.get('location', ''),
                result_data.get('timestamp', '')
            )
            
            item_id = self.results_tree.insert('', 'end', values=values)
            
            # Auto-scroll to new result
            self.results_tree.see(item_id)
        
        # Update chromosome filter options if new chromosome
        current_chromosomes = set(self.chromosome_filter['values']) - {'All'}
        if chromosome not in current_chromosomes:
            new_chromosomes = sorted(current_chromosomes | {chromosome})
            self.chromosome_filter['values'] = ['All'] + new_chromosomes
    
    def _update_sequence_tracking(self, data):
        """Update sequence tracking information."""
        if isinstance(data, dict):
            self.search_info.update(data)
    
    def _add_results_batch(self, results_data):
        """Add multiple results at once."""
        if isinstance(results_data, list):
            for result in results_data:
                self._add_new_result(result)
    
    def _schedule_updates(self):
        """Schedule regular GUI updates."""
        self._update_display()
        if self.running:
            self.root.after(1000, self._schedule_updates)  # Update every second
    
    def _update_display(self):
        """Update all display elements with enhanced information."""
        try:
            # Update clock
            current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            self.clock_label.config(text=current_time)
            
            # Update search info
            if self.search_info.get('pattern'):
                self.pattern_label.config(text=f"Pattern: {self.search_info['pattern']}")
            
            if self.search_info.get('algorithm'):
                self.algorithm_label.config(text=f"Algorithm: {self.search_info['algorithm']}")
            
            # Update sequence info with number
            current_seq = self.search_info.get('current_sequence', '')
            seq_num = self.search_info.get('sequence_number', 0)
            total_seqs = self.search_info.get('total_sequences', 0)
            
            if current_seq:
                seq_text = f"Sequence {seq_num}/{total_seqs}: {current_seq[:30]}{'...' if len(current_seq) > 30 else ''}"
                self.sequence_label.config(text=seq_text)
            
            # Update file info
            sequence_file = self.search_info.get('sequence_file', '')
            if sequence_file:
                file_name = os.path.basename(sequence_file)
                self.file_label.config(text=f"File: {file_name}")
            
            # Update batch progress
            if total_seqs > 0:
                self.batch_progress_label.config(text=f"Batch: {seq_num}/{total_seqs} patterns")
            
            # Update chromosomes
            chromosomes = self.search_info.get('chromosomes', [])
            if chromosomes:
                chrom_text = f"Chromosomes: {len(chromosomes)} total"
                if len(chromosomes) <= 5:
                    chrom_text += f" ({', '.join(chromosomes)})"
                self.chromosomes_label.config(text=chrom_text)
            
            # Update overall progress
            completed = self.overall_progress.get('completed', 0)
            total = self.overall_progress.get('total', 0)
            percentage = self.overall_progress.get('percentage', 0)
            
            self.overall_label.config(text=f"{completed}/{total} chromosomes ({percentage:.1f}%)")
            self.overall_progress_bar['value'] = percentage
            
            # Update matches count
            total_matches = len(self.search_results)
            self.matches_label.config(text=f"📊 {total_matches:,} matches found")
            
            # Update results count
            filtered_count = len(self.filtered_results)
            self.results_count_label.config(text=f"Showing: {filtered_count:,} / {total_matches:,} results")
            
            # Update statistics
            current_time_val = time.time()
            duration = current_time_val - self.search_statistics['start_time']
            
            hours = int(duration // 3600)
            minutes = int((duration % 3600) // 60)
            seconds = int(duration % 60)
            duration_str = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
            
            self.duration_label.config(text=f"⏱️ Duration: {duration_str}")
            
            # Calculate speed
            if duration > 0:
                speed = total_matches / duration
                self.speed_label.config(text=f"⚡ Speed: {speed:.1f} matches/sec")
            
            # Update status based on overall progress
            if percentage == 100 and total > 0:
                self.status_label.config(text="✅ Search completed!", fg=self.colors['accent_primary'])
            elif percentage > 0:
                self.status_label.config(text=f"🔄 Searching... {percentage:.1f}% complete", fg=self.colors['accent_warning'])
            else:
                self.status_label.config(text="🔄 Ready for search...", fg=self.colors['text_secondary'])
                
        except Exception as e:
            print(f"[DEBUG] Display update error: {e}")
    
    def run(self):
        """Start the enhanced dashboard."""
        try:
            self.root.mainloop()
        finally:
            self.running = False


# Legacy class name for backward compatibility
class EnhancedWorkerDashboard(UltimateGenomeDashboard):
    """Legacy class name - use UltimateGenomeDashboard instead."""
    pass


# Function required by Genome_Analyzer.py
def start_dashboard(command_queue: multiprocessing.Queue, worker_count: int = 8):
    """Start the dashboard - required by the main genome analyzer."""
    try:
        dashboard = UltimateGenomeDashboard(command_queue, worker_count)
        dashboard.run()
    except Exception as e:
        print(f"[ERROR] Dashboard failed to start: {e}")
        import traceback
        traceback.print_exc()
