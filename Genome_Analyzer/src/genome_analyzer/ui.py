import os
import sys
import multiprocessing
from datetime import datetime
from typing import Dict, List, Tuple, Optional

from .analyzer import GenomeAnalyzer, GenomeAnalysisError
from . import file_processing
from .utils import _get_filtered_chromosomes, _perform_cache_cleanup
from . import file_processing

# Global frame width for pretty output
FRAME_WIDTH: int = 0

# Column width constants for result formatting (eliminate duplication)
CHR_W, START_W, END_W, STRAND_W, PAT_W, FWD_W, REV_W, MIS_W = 6, 8, 8, 2, 25, 25, 25, 3
CHROM_TOTAL_WIDTH = CHR_W + START_W + END_W + STRAND_W + PAT_W + FWD_W + REV_W + MIS_W + 8  # +8 for spacing
NON_CHROM_TOTAL_WIDTH = START_W + END_W + STRAND_W + PAT_W + FWD_W + REV_W + MIS_W + 6  # +6 for spacing
SEARCH_BANNER_WIDTH = 62

def _print_results_header(include_chrom: bool):
    if include_chrom:
        header = (
            f"{ 'Chr':<{CHR_W}} {'Start':<{START_W}} {'End':<{END_W}} {'S':<{STRAND_W}} "
            f"{ 'Pattern':<{PAT_W}} {'Found':<{FWD_W}} {'RevComp':<{REV_W}} {'Mis':<{MIS_W}} Location"
        )
    else:
        header = (
            f"{ 'Start':<{START_W}} {'End':<{END_W}} {'S':<{STRAND_W}} "
            f"{ 'Pattern':<{PAT_W}} {'Found':<{FWD_W}} {'RevComp':<{REV_W}} {'Mis':<{MIS_W}} Location"
        )
    global FRAME_WIDTH
    FRAME_WIDTH = len(header)
    border = "+" + "-" * (FRAME_WIDTH + 2) + "+"
    print(border)
    print("| " + header + " |")

def _print_bottom_frame():
    border = "+" + "-" * (FRAME_WIDTH + 2) + "+"
    print(border)

def _print_row_tuple(row_tuple: Tuple, include_chrom: bool):
    if include_chrom:
        chrom, start, end, strand, pat, fwd, rev, mis, loc = row_tuple
        line = (
            f"{str(chrom):<{CHR_W}} {start:<{START_W}} {end:<{END_W}} {strand:<{STRAND_W}} "
            f"{pat[:PAT_W]:<{PAT_W}} {fwd[:FWD_W]:<{FWD_W}} {rev[:REV_W]:<{REV_W}} {mis:<{MIS_W}} {loc}"
        )
    else:
        start, end, strand, pat, fwd, rev, mis, loc = row_tuple
        line = (
            f"{start:<{START_W}} {end:<{END_W}} {strand:<{STRAND_W}} "
            f"{pat[:PAT_W]:<{PAT_W}} {fwd[:FWD_W]:<{FWD_W}} {rev[:REV_W]:<{REV_W}} {mis:<{MIS_W}} {loc}"
        )
    # FIXED: Remove brutal line truncation that was cutting off location information
    # line = line[:FRAME_WIDTH]  # This was chopping off gene names!
    padding = " " * max(0, FRAME_WIDTH - len(line))
    print("| " + line + padding + " |", flush=True)

def _print_search_banner(species: str, genome_path: str, annotation_path: str,
                         target: str, pattern: str, mismatches: int, boundary_bp: int):
    gname = os.path.basename(genome_path)
    aname = os.path.basename(annotation_path)
    line = "+" + "-" * SEARCH_BANNER_WIDTH + "+"
    print("\n" + line)
    print("| {:<60} |".format("Search Context"))
    print(line)
    print("| {:<60} |".format(f"Species       : {species}"))
    print("| {:<60} |".format(f"Target        : {target}"))
    print("| {:<60} |".format(f"Pattern       : {pattern}    Mismatches: {mismatches}  Boundary: {boundary_bp}bp"))
    print("| {:<60} |".format(f"Genome file   : {gname}"))
    print("| {:<60} |".format(f"Annotation    : {aname}"))
    print(line)

def _prompt_numeric_choice(options: List[str], title: str) -> int:
    try:
        while True:
            print(f"\n{title}")
            for i, opt in enumerate(options, 1):
                print(f"{i}. {opt}")
            
            print(f"\nNavigation: 'b' (back), 'm' (main menu), 'q' (quit)")
            
            choice_str = input(f"Choose (1-{len(options)}): ").strip().lower()
            
            if not choice_str:
                print(f"Using default: {options[0]}")
                return 0
            
            if choice_str in ['b', 'back']:
                return -1
            elif choice_str in ['m', 'main', 'menu']:
                return -2
            elif choice_str in ['q', 'quit', 'exit']:
                raise GenomeAnalysisError("User chose to exit")
            
            if choice_str.isdigit():
                idx = int(choice_str) - 1
                if 0 <= idx < len(options):
                    return idx
            print("Invalid choice. Try again.")
    except KeyboardInterrupt:
        print("\nOperation cancelled by user. Exiting gracefully.")
        raise GenomeAnalysisError("Operation cancelled by user")

def launch_dashboard():
    try:
        import sys
        import os
        # Add parent directory to path to import worker_dashboard
        parent_dir = os.path.dirname(os.path.dirname(__file__))
        sys.path.insert(0, parent_dir)
        from worker_dashboard import start_dashboard
        
        command_queue = multiprocessing.Queue()
        
        dashboard_process = multiprocessing.Process(
            target=start_dashboard, 
            args=(command_queue, os.cpu_count())
        )
        dashboard_process.start()
        
        print(f"\n🎯 Dashboard launched! A new window should appear.")
        print(f"💡 Dashboard will update in real-time as searches progress!")
        print(f"💡 You can minimize this console and watch the dashboard instead!")
        
        return command_queue, dashboard_process
        
    except ImportError:
        print(f"\n[WARNING] Dashboard not available. Install tkinter to use this feature.")
        return None, None
    except Exception as e:
        print(f"\n[WARNING] Could not start dashboard: {e}")
        return None, None

def cleanup_dashboard(dashboard_process, command_queue):
    if dashboard_process is not None:
        print("\n[INFO] Cleaning up dashboard...")
        try:
            if command_queue is not None:
                try:
                    command_queue.put(("CLOSE", None), timeout=2)
                except:
                    pass
            
            dashboard_process.join(timeout=5)
            
            if dashboard_process.is_alive():
                print("[INFO] Force terminating dashboard...")
                dashboard_process.terminate()
                dashboard_process.join(timeout=2)
                
                if dashboard_process.is_alive():
                    print("[WARNING] Dashboard process still running - killing...")
                    try:
                        dashboard_process.kill()
                    except:
                        pass
                        
            print("[INFO] Dashboard cleanup completed.")
            
        except Exception as e:
            print(f"[WARNING] Dashboard cleanup encountered an error: {e}")

def offer_dashboard():
    print(f"\nWould you like to launch the real-time worker dashboard?")
    print(f"   This will show live progress of all workers in a separate window.")
    print(f"   Type 'y' to launch dashboard, any other key to continue without it.")
    
    choice = input("Launch dashboard? (y/N): ").strip().lower()
    if choice == 'y':
        return launch_dashboard()
    return None, None

def search_mode_prompt(analyzer, species_name, complete_genome_file, gene_annotation_file, selected_chrom, multi_chrom, results_dir, seq_files_dir, command_queue, cache_loader):
    print("\nChoose search mode:")
    print("1. Single search in selected genome")
    print("2. Batch search from substrate file")
    print("3. Quit")
    try:
        mode = input("Enter mode (1-3): ")
    except KeyboardInterrupt:
        print("\nOperation cancelled by user. Exiting.")
        return

    if mode == '1':
        try:
            pattern = input("Enter sequence to search (IUPAC supported): ").strip().upper()
        except KeyboardInterrupt:
            print("\nOperation cancelled by user. Returning to main menu.")
            return
        if not pattern: return
        try:
            mismatches_str = input("Enter maximum allowed mismatches (default 0): ")
        except KeyboardInterrupt:
            print("\nOperation cancelled by user. Returning to main menu.")
            return
        max_mismatches = int(mismatches_str) if mismatches_str.isdigit() else 0
        # Removed boundary prompt as requested
        boundary_bp = 0
        _print_search_banner(species_name, complete_genome_file, gene_annotation_file,
                             selected_chrom if selected_chrom else 'Whole genome',
                             pattern, max_mismatches, boundary_bp)
        # Simplified: Always use single genome search for Homo Sapiens
        print(f"\n[INFO] Single genome search mode")
        from .parallel import parallel_single_genome_search
        matches = parallel_single_genome_search(analyzer, pattern, max_mismatches, boundary_bp, command_queue)
        print(f"\n[INFO] Search completed. Total matches found: {len(matches)}")
        if not matches: 
            print("[INFO] No matches found for this pattern.")
            return
        save_prompt = input("\nDo you want to save these results to a CSV file? (y/n): ").lower()
        if save_prompt == 'y':
            user_name = input("Enter a name for this result file (no extension): ").strip()
            if not user_name:
                user_name = "results"
            date_str = datetime.now().strftime("%Y%m%d")
            chrom_suffix = ''
            if gene_annotation_file.lower().endswith('.gtf'):
                if selected_chrom and selected_chrom != 'ALL':
                    chrom_suffix = f"-{selected_chrom}"
                elif multi_chrom:
                    chrom_suffix = "-ALL"
            output_filename = f"{species_name}-{user_name}{chrom_suffix}-{date_str}.csv"
            output_filepath = os.path.join(results_dir, output_filename)
            file_processing.save_results_to_csv(matches, output_filepath)

    elif mode == '2':
        try:
            available_files = [f for f in os.listdir(seq_files_dir) if f.endswith(('.fasta', '.txt', '.fa'))]
            if not available_files:
                print(f"\nNo files found in '{seq_files_dir}'."); return
            print("\nAvailable substrate files:")
            for i, filename in enumerate(available_files): print(f"{i + 1}. {filename}")
            try:
                choice = int(input(f"Choose a file for batch search (1-{len(available_files)}): ")) - 1
            except KeyboardInterrupt:
                print("\nOperation cancelled by user. Returning to main menu.")
                return
            if not 0 <= choice < len(available_files):
                print("Invalid choice."); return
            selected_filename = available_files[choice]
            selected_filepath = os.path.join(seq_files_dir, selected_filename)
            queries = file_processing.parse_fasta_file(selected_filepath)
            if not queries:
                print(f"No sequences found in {selected_filename}."); return
            try:
                mismatches_str = input("Enter maximum allowed mismatches for all searches (default 0): ")
            except KeyboardInterrupt:
                print("\nOperation cancelled by user. Returning to main menu.")
                return
            max_mismatches = int(mismatches_str) if mismatches_str.isdigit() else 0
            try:
                boundary_str = input("Flag near-chromosome-end matches within N bp (blank to skip): ").strip()
            except KeyboardInterrupt:
                print("\nOperation cancelled by user. Returning to main menu.")
                return
            boundary_bp = int(boundary_str) if boundary_str.isdigit() else 0
            all_matches = []
            total_sequences = len(queries)
            current_sequence = 0
            
            # CACHE: Store genome data to avoid reloading for each pattern
            cached_genome_data = None
            cached_annotations = None
            cached_boundaries = None
            
            print(f"\n[INFO] Starting batch search for {total_sequences} patterns")
            print(f"[INFO] Results displayed in dashboard window")
            
            if command_queue:
                command_queue.put(("SEQUENCE_UPDATE", {
                    'filename': selected_filename,
                    'total_sequences': total_sequences,
                    'sequence_number': 0
                }))
                
                command_queue.put(("SEARCH_INFO", {
                    'sequence_file': selected_filepath,
                    'total_sequences': total_sequences
                }))
            
            for header, pattern in queries.items():
                current_sequence += 1
                
                if command_queue:
                    command_queue.put(("SEQUENCE_UPDATE", {
                        'sequence_number': current_sequence,
                        'total_sequences': total_sequences,
                        'current_sequence': header[:50]
                    }))
                
                print(f"[INFO] Pattern {current_sequence}/{total_sequences}: {header[:50]}{'...' if len(header) > 50 else ''}")
                
                if multi_chrom:
                    gtf_path = gene_annotation_file
                    fasta_path = complete_genome_file
                    
                    # Use the main chromosomes list directly (25 chromosomes: chr1-22, chrX, chrY, chrM)
                    if hasattr(analyzer, 'available_chromosomes') and analyzer.available_chromosomes:
                        chroms = analyzer.available_chromosomes
                        print(f"[DEBUG] Using analyzer's chromosome list: {len(chroms)} chromosomes")
                    else:
                        # Fallback to filtered chromosomes
                        chroms = _get_filtered_chromosomes(gtf_path, fasta_path)
                        print(f"[DEBUG] Using filtered chromosomes: {len(chroms)} chromosomes")
                    
                    print(f"[DEBUG] Chromosomes for search: {chroms}")
                    
                    # Import the function here to avoid scope issues
                    from .parallel import enhanced_parallel_chromosome_search
                    
                    # CACHE OPTIMIZATION: Use cached data to avoid reloading for each pattern
                    matches = enhanced_parallel_chromosome_search(chroms, fasta_path, gtf_path, 
                                                                   pattern, max_mismatches, boundary_bp, command_queue,
                                                                   cached_genome_data, cached_annotations, cached_boundaries)
                    
                else:
                    # Import the function here to avoid scope issues
                    from .parallel import parallel_single_genome_search
                    matches = parallel_single_genome_search(analyzer, pattern, max_mismatches, boundary_bp, command_queue)
                    
                    if command_queue and matches:
                        for match in matches:
                            if len(match) == 8:
                                full_match = ("N/A",) + tuple(match)
                            else:
                                full_match = tuple(match)
                            command_queue.put(("NEW_RESULT", full_match))
                        
                for match in matches:
                    original_location = match[-1]
                    full_match_info = list(match)
                    full_match_info[-1] = f"{original_location} (Query: {header[:30]}...)"
                    all_matches.append(tuple(full_match_info))
            
            print(f"\n[INFO] Batch search complete: {len(all_matches)} total matches found")
            print(f"[INFO] All results displayed in dashboard window")
            
            _perform_cache_cleanup(cache_loader, "batch search")
            
            if not all_matches: return
            save_prompt = input("\nDo you want to save the complete results to a CSV file? (y/n): ").lower()
            if save_prompt == 'y':
                user_name = input("Enter a name for this result file (no extension): ").strip()
                if not user_name:
                    user_name = "results"
                date_str = datetime.now().strftime("%Y%m%d")
                chrom_suffix = ''
                if gene_annotation_file.lower().endswith('.gtf'):
                    if selected_chrom and selected_chrom != 'ALL':
                        chrom_suffix = f"-{selected_chrom}"
                    elif multi_chrom:
                        chrom_suffix = "-ALL"
                output_filename = f"{species_name}-{user_name}{chrom_suffix}-{date_str}.csv"
                output_filepath = os.path.join(results_dir, output_filename)
                file_processing.save_results_to_csv(all_matches, output_filepath)
        except (ValueError, IndexError):
            print("Invalid input.")
        except Exception as e:
            print(f"An error occurred: {e}")

    elif mode == '3':
        print("Goodbye!")
        sys.exit(0)
    else:
        print("Invalid mode.")

def _navigator_select_inputs(data_dir: str) -> Tuple[str, str, str, Optional[str]]:
    """Navigator to select species, genome, and annotation files."""
    from .utils import (_list_subdirectories, _list_files_with_ext, _discover_gtf_chromosomes, 
                       _filter_human_chromosomes, _get_default_genome_file, _get_default_annotation_file)
    
    while True:
        # Species selection by navigating subfolders under data_dir
        species_dirs = _list_subdirectories(data_dir)
        if not species_dirs:
            raise GenomeAnalysisError("No species subdirectories found under 'data'.")
        
        sp_idx = _prompt_numeric_choice(species_dirs, "Available species (data subfolders):")
        if sp_idx < 0:  # Navigation signal
            if sp_idx == -2:  # Main menu
                return None, None, None, None
            continue  # Back signal, restart species selection
        
        species_dir_name = species_dirs[sp_idx]
        species_dir_path = os.path.join(data_dir, species_dir_name)

        # Genome file selection with default option
        genome_candidates = _list_files_with_ext(species_dir_path, ('.fa', '.fasta', '.fna', '.txt'))
        if not genome_candidates:
            raise GenomeAnalysisError(f"No genome FASTA/text found in {species_dir_path}")
        
        # Add default option based on species
        default_genome = _get_default_genome_file(species_dir_name)
        if default_genome and default_genome in genome_candidates:
            # Remove the original entry to avoid duplication
            genome_candidates.remove(default_genome)
            genome_candidates.insert(0, f"{default_genome} (DEFAULT)")
        
        g_idx = _prompt_numeric_choice(genome_candidates, f"Genome files in {species_dir_name}:")
        if g_idx < 0:  # Navigation signal
            if g_idx == -2:  # Main menu
                return None, None, None, None
            continue  # Back signal, restart species selection
        
        genome_filename = genome_candidates[g_idx].replace(" (DEFAULT)", "")

        # Annotation file selection with default option
        anno_candidates = _list_files_with_ext(species_dir_path, ('.gtf', '.gff', '.gff3', '.txt', '.json'))
        if not anno_candidates:
            raise GenomeAnalysisError(f"No annotation files found in {species_dir_path}")
        
        # Add default option based on species
        default_annotation = _get_default_annotation_file(species_dir_name)
        if default_annotation and default_annotation in anno_candidates:
            # Remove the original entry to avoid duplication
            anno_candidates.remove(default_annotation)
            anno_candidates.insert(0, f"{default_annotation} (DEFAULT)")
        
        a_idx = _prompt_numeric_choice(anno_candidates, f"Annotation files in {species_dir_name}:")
        if a_idx < 0:  # Navigation signal
            if a_idx == -2:  # Main menu
                return None, None, None, None
            continue  # Back signal, restart genome selection

        annotation_filename = anno_candidates[a_idx].replace(" (DEFAULT)", "")

        chrom = None
        if annotation_filename.lower().endswith('.gtf'):
            gtf_path = os.path.join(species_dir_path, annotation_filename)
            fasta_path = os.path.join(species_dir_path, genome_filename)
            chroms = _discover_gtf_chromosomes(gtf_path)
            chroms = _filter_human_chromosomes(chroms, fasta_path)
            if not chroms:
                raise GenomeAnalysisError("Could not discover chromosomes from GTF.")
            menu = ["ALL (primary chromosomes)"] + chroms
            c_idx = _prompt_numeric_choice(menu, "Select chromosome to analyze:")
            if c_idx < 0:  # Navigation signal
                if c_idx == -2:  # Main menu
                    return None, None, None, None
                continue  # Back signal, restart annotation selection
            chrom = 'ALL' if c_idx == 0 else chroms[c_idx - 1]

        return species_dir_path, genome_filename, annotation_filename, chrom

def _print_formatted_results(matches: List, include_chrom: bool = False):
    """
    Centralized result formatting to eliminate code duplication.
    Handles both chromosome and non-chromosome result formats.
    """
    if not matches:
        print("No matches found.")
        return
    
    # Print header
    _print_results_header(include_chrom)
    
    # Print separator using pre-calculated constants
    separator = '-' * (CHROM_TOTAL_WIDTH if include_chrom else NON_CHROM_TOTAL_WIDTH)
    
    print(separator)
    
    # Print results
    for match in matches:
        _print_row_tuple(match, include_chrom)
    
    print(separator)