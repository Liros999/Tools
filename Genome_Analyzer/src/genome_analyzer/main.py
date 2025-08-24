import os
import sys
import pickle
import multiprocessing
import signal
from datetime import datetime
from typing import Dict

from .analyzer import GenomeAnalyzer, GenomeAnalysisError
from . import file_processing
from . import ui
from .utils import _handle_sigint

def main():
    """
    SIMPLIFIED: Direct genome analysis for Homo Sapiens only.
    """
    # Set up signal handling for graceful exit
    signal.signal(signal.SIGINT, _handle_sigint)
    
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        script_dir = os.getcwd()
    
    # SIMPLIFIED: Use direct paths to known genome files for Homo Sapiens only
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # FIXED: Correct path logic - go up two levels from src/genome_analyzer to Genome_Analyzer
    src_dir = os.path.dirname(script_dir)  # Go up from src/genome_analyzer to src
    project_root = os.path.dirname(src_dir)  # Go up from src to Genome_Analyzer
    data_dir = os.path.join(project_root, "data")  # Genome_Analyzer/data
    seq_files_dir = os.path.join(src_dir, "Seq_Files")  # src/Seq_Files
    results_dir = os.path.join(src_dir, "Results")  # src/Results
    
    # Define paths for Homo Sapiens only
    species_paths = {
        "Homo_Sapiens": {
            "data_dir": os.path.join(data_dir, "Homo_Sapiens"),
            "genome_file": "GRCh38.primary_assembly.genome.fa",
            "annotation_file": "gencode.v48.primary_assembly.annotation.gtf"
        }
    }
    
    # Create directories if they don't exist
    os.makedirs(seq_files_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print(f"[INFO] Genome Analyzer - Homo Sapiens Only Version")
    print(f"[INFO] Available species: Homo Sapiens")
    print(f"[INFO] Results directory: {results_dir}")
    print(f"[INFO] Sequence files directory: {seq_files_dir}")
    
    # Homo Sapiens is the only option
    selected_species = "Homo_Sapiens"
    print(f"[INFO] Selected species: {selected_species}")
    
    # Use selected species paths
    species_dir_path = species_paths[selected_species]['data_dir']
    genome_filename = species_paths[selected_species]['genome_file']
    annotation_filename = species_paths[selected_species]['annotation_file']
    
    print(f"[INFO] Using {selected_species} genome: {genome_filename}")
    print(f"[INFO] Using {selected_species} annotations: {annotation_filename}")
    
    # Launch dashboard automatically (necessary for monitoring search progress)
    print(f"\n[INFO] Launching real-time worker dashboard...")
    command_queue, dashboard_process = ui.launch_dashboard()
    
    # Define species_name for compatibility
    species_name = "homo sapiens"
    selected_chrom = None  # Single genome file

    complete_genome_file = os.path.join(species_dir_path, genome_filename)
    gene_annotation_file = os.path.join(species_dir_path, annotation_filename)
    
    # Define these variables for compatibility
    genome_file_path = os.path.join(species_dir_path, genome_filename)
    annotation_file_path = os.path.join(species_dir_path, annotation_filename)

    # Per-species cache
    cache_file = os.path.join(species_dir_path, "genome_cache.pkl")

    # Initialize the analyzer for Homo Sapiens only
    analyzer = GenomeAnalyzer()
    analyzer.species_type = "generic"  # Homo Sapiens uses generic logic
    
    # For compatibility, set species_analyzer to the same analyzer
    species_analyzer = analyzer

    # Initialize cache management
    from . import file_processing
    try:
        cache_loader = file_processing.MemoryMappedGenomeLoader()
        print(f"[INFO] Cache management initialized (limit: {cache_loader.max_cache_size_mb}MB)")
    except Exception as e:
        print(f"[WARNING] Cache initialization failed: {e}")
        cache_loader = None

    # Load from cache only if it matches current inputs
    multi_chrom = True  # Homo Sapiens has 25 main chromosomes (chr1-22, chrX, chrY, chrM)
    per_chrom_cache: Dict[str, GenomeAnalyzer] = {}
    cache_ok = False # Initialize cache_ok
    
    # Get list of main chromosomes to analyze
    main_chromosomes = []
    for i in range(1, 23):  # chr1-chr22
        main_chromosomes.append(f"chr{i}")
    main_chromosomes.extend(['chrX', 'chrY', 'chrM'])  # Sex chromosomes and mitochondrial
    
    print(f"[INFO] Homo Sapiens has {len(main_chromosomes)} main chromosomes: {main_chromosomes}")
    
    # For multi-chromosome analysis, we don't use per-species cache
    # Instead, the enhanced_parallel_chromosome_search function will handle loading each chromosome as needed
    if multi_chrom:
        print(f"[INFO] Multi-chromosome mode enabled for {len(main_chromosomes)} chromosomes")
        print(f"[INFO] Chromosomes will be loaded on-demand during search")
        print(f"[INFO] Main analyzer will be used for coordination only")
        
        # Set up analyzer for multi-chromosome coordination
        analyzer.species_type = "generic"
        analyzer.multi_chromosome_mode = True
        analyzer.available_chromosomes = main_chromosomes
        analyzer.genome_file_path = complete_genome_file
        analyzer.annotation_file_path = annotation_file_path
        
        print(f"[INFO] Analyzer configured for multi-chromosome processing")
    else:
        # Single chromosome mode (fallback)
        print(f"[INFO] Single chromosome mode - loading chr1 only")
        
        if os.path.exists(cache_file):
            print(f"Loading cached genome data from: {os.path.basename(cache_file)}...")
            with open(cache_file, 'rb') as f:
                cached_data = pickle.load(f)
            cache_ok = (
                cached_data.get('source_genome_path') == complete_genome_file and
                cached_data.get('source_annotation_path') == annotation_file_path and
                cached_data.get('chromosome') == None # No chromosome for single genome
            )
            if cache_ok:
                # Load data into analyzer if available
                analyzer.sequence = cached_data['sequence']
                analyzer.gene_tree = cached_data['gene_tree']
                analyzer.gene_data = cached_data['gene_data']

                # Set species type for the analyzer
                analyzer.species_type = "generic"

                print("Cached data loaded successfully.")
            else:
                print("Cache does not match selected inputs. Rebuilding...")

                # Load genome sequence using standard file processing
                file_processing.load_complete_sequence(analyzer, complete_genome_file)

                # Load annotations using GTF format for Homo Sapiens
                file_processing.load_annotations_from_gtf(analyzer, annotation_file_path, "chr1")

                # Save to cache
                with open(cache_file, 'wb') as f:
                    pickle.dump({
                        'sequence': analyzer.sequence,
                        'gene_tree': analyzer.gene_tree,
                        'gene_data': analyzer.gene_data,
                        'source_genome_path': complete_genome_file,
                        'source_annotation_path': annotation_file_path,
                        'chromosome': None
                    }, f)
                print(f"Saved new genome data to cache file: {os.path.basename(cache_file)}")
        else:
            # Load genome sequence using standard file processing
            file_processing.load_complete_sequence(analyzer, complete_genome_file)

            # Load annotations using GTF format for Homo Sapiens  
            file_processing.load_annotations_from_gtf(analyzer, annotation_file_path, "chr1")

            # Save to cache
            with open(cache_file, 'wb') as f:
                pickle.dump({
                    'sequence': analyzer.sequence,
                    'gene_tree': analyzer.gene_tree,
                    'gene_data': analyzer.gene_data,
                    'source_genome_path': complete_genome_file,
                    'source_annotation_path': annotation_file_path,
                    'chromosome': None
                }, f)
            print(f"Saved new genome data to cache file: {os.path.basename(cache_file)}")

    # CRITICAL FIX: Always show search interface after loading data (cache or fresh)
    print(f"[INFO] Files found and ready for analysis")
    
    try:
        while True:
            # Use analyzer for Homo Sapiens analysis with chromosome list
            ui.search_mode_prompt(analyzer, species_name, complete_genome_file, gene_annotation_file, main_chromosomes, multi_chrom, results_dir, seq_files_dir, command_queue, cache_loader)
    except KeyboardInterrupt:
        print("\n\n[INFO] Program interrupted by user. Exiting gracefully...")
    except Exception as e:
        print(f"\n[ERROR] An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Always cleanup dashboard - this runs no matter how we exit
        ui.cleanup_dashboard(dashboard_process, command_queue)

if __name__ == "__main__":
    main()