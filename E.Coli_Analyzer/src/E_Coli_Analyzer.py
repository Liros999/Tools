# C:\Users\לאה\Desktop\Tools\E.Coli_Analyzer\src\E_Coli_Analyzer.py
# --- FINAL, COMPLETE AND WORKING VERSION ---

import os
import re
import csv
import pickle
import sys
from typing import Dict, List, Tuple
from tqdm import tqdm
from intervaltree import Interval, IntervalTree

class EColiAnalyzer:
    """
    A comprehensive tool to analyze the E. coli genome, featuring a complete
    genomic map, advanced location reporting, and multiple search modes.
    """
    def __init__(self):
        self.sequence: str = ""
        self.gene_data: Dict[str, Dict] = {}
        self.gene_tree = IntervalTree()
        self.iupac_map = {
            'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]', 'K': '[GT]',
            'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
            'N': '[ACGT]'
        }
        self.iupac_wildcards = {k: v.strip('[]') for k, v in self.iupac_map.items()}

    def _reverse_complement(self, seq: str) -> str:
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 
                          'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
        return "".join(complement_map.get(base, base) for base in reversed(seq))

    def _count_mismatches(self, pattern: str, window: str) -> int:
        mismatches = 0
        for p_char, w_char in zip(pattern, window):
            if p_char in self.iupac_wildcards:
                if w_char not in self.iupac_wildcards[p_char]: mismatches += 1
            elif p_char != w_char: mismatches += 1
        return mismatches

    def _extract_flanked_window(self, start: int, end: int, strand: str, left_flank: int, right_flank: int) -> Tuple[str, str]:
        """
        Extract a flanked sequence window around a 1-based inclusive match [start,end].
        Applies strand-aware flanks and reverse-complements for '-' strand.
        Returns (sequence, coordinates_str_0_based_half_open)
        """
        if strand == '-':
            flank_up = right_flank
            flank_down = left_flank
        else:
            flank_up = left_flank
            flank_down = right_flank
        slice_start = max(0, (start - 1) - flank_up)
        slice_end = min(len(self.sequence), end + flank_down)
        window = self.sequence[slice_start:slice_end]
        if strand == '-':
            window = self._reverse_complement(window)
        coord = f"{slice_start}-{slice_end}"
        return window, coord

    def display_alignment_pairwise(self, reference_seqs: List[str], targets: List[Dict], top_n: int = 100):
        """
        Display best-window pairwise alignment: each target compared to each reference with star line.
        Targets are sorted by best mismatch score; best printed last.
        """
        if not reference_seqs or not targets:
            print("No references or targets to display."); return
        name_width = max(len(t['name']) for t in targets) if targets else 10
        info_width = 28
        print(f"\nPairwise alignment relative to references (top {top_n} shown; best scores last):")
        print("=" * 80)
        for r_index, ref_raw in enumerate(reference_seqs, start=1):
            ref_name = f"Reference_{r_index}"
            # Compute best score per target
            scored = []
            win_len = len(ref_raw)
            for t in targets:
                seq = t['sequence']
                best_mm, best_idx = None, 0
                for i in range(0, max(1, len(seq) - win_len + 1)):
                    mm = self._count_mismatches(ref_raw, seq[i:i+win_len])
                    if best_mm is None or mm < best_mm:
                        best_mm, best_idx = mm, i
                scored.append({
                    'name': t['name'],
                    'coords': t.get('coordinates', 'N/A'),
                    'strand': t.get('strand', 'N/A'),
                    'best_mm': 0 if best_mm is None else best_mm,
                    'window': seq[best_idx:best_idx+win_len]
                })
            scored.sort(key=lambda x: (x['best_mm'], x['name']))
            shown = scored[:top_n]
            for item in reversed(shown):
                tgt_name = item['name']
                tgt_window = item['window']
                info = f"({item['coords']}) [{item['strand']}] score={item['best_mm']}"
                print(f"\nReference: {ref_name}    Target: {tgt_name} {info}")
                stars = ''.join('*' if a == b else ' ' for a, b in zip(ref_raw, tgt_window))
                print(f"{ref_name.ljust(name_width)}  {ref_raw}")
                print(f"{tgt_name.ljust(name_width)}  {tgt_window}  {info.ljust(info_width)}")
                print(f"{'':<{name_width}}  {stars}")

    def perform_alignment_analysis(self):
        """
        Mode: search target pattern with mismatches, extract flanks around hits,
        then display pairwise star-line alignments versus provided reference(s).
        """
        print("\nSequence Alignment with Flanking Regions (E. coli)")
        pattern = input("Enter the target sequence pattern: ").strip().upper()
        if not pattern:
            print("Please enter a non-empty sequence pattern."); return
        valid_chars = set('ATCGRYSWKMBDHVN')
        if set(pattern) - valid_chars:
            print("Invalid characters in pattern."); return
        max_mismatches = input("Enter maximum allowed mismatches for search (default 0): ").strip()
        max_mismatches = int(max_mismatches) if max_mismatches.isdigit() else 0
        left_flank = input("Enter number of nucleotides to include from the left: ").strip()
        right_flank = input("Enter number of nucleotides to include from the right: ").strip()
        left_flank = int(left_flank) if left_flank.isdigit() else 50
        right_flank = int(right_flank) if right_flank.isdigit() else 50
        # Reference sequences
        references = []
        print("\nEnter reference sequences (one per line, empty line to finish):")
        while True:
            ref_seq = input("Reference sequence: ").strip().upper()
            if not ref_seq: break
            if set(ref_seq) - valid_chars:
                print("Invalid characters in reference."); continue
            references.append(ref_seq)
        if not references:
            print("No reference sequences provided."); return
        # Search genome
        matches = self.search_sequence(pattern, max_mismatches)
        if matches:
            filter_option = input("\nWould you like to filter results? (y/n): ").lower().strip()
            if filter_option == 'y':
                print("\nFilter options:")
                print("1. Show only matches inside genes")
                print("2. Show only matches between genes (intergenic)")
                print("3. Show only matches with specific number of mismatches")
                print("4. Show only matches in specific region")
                choice = input("Enter filter choice (1-4): ").strip()
                if choice == '1':
                    matches = [m for m in matches if '[In]' in m[-1]] or matches
                elif choice == '2':
                    matches = [m for m in matches if '[In]' not in m[-1]] or matches
                elif choice == '3':
                    try:
                        k = int(input("Enter maximum number of mismatches: "))
                        matches = [m for m in matches if m[6] <= k] or matches
                    except ValueError:
                        print("Invalid input. Skipping filter.")
                elif choice == '4':
                    try:
                        region_input = input("Enter region (start,end): ")
                        rstart, rend = map(int, region_input.split(','))
                        filtered = []
                        for m in matches:
                            if (rstart <= m[0] <= rend or rstart <= m[1] <= rend or m[0] <= rstart <= m[1] or m[0] <= rend <= m[1]):
                                filtered.append(m)
                        matches = filtered or matches
                    except ValueError:
                        print("Invalid region format. Skipping filter.")
        if not matches:
            print("No matches found for the pattern."); return

        # Ask how many top results to align and pre-rank to avoid heavy processing
        try:
            top_n_input = input("\nHow many top results do you want to align? (default 100): ").strip()
            top_n = int(top_n_input) if top_n_input else 100
        except ValueError:
            top_n = 100

        # Rank by raw mismatches asc, then by start; slice to top_n
        matches.sort(key=lambda m: (m[6], m[0]))
        if len(matches) > top_n:
            matches = matches[:top_n]

        print(f"\nPreparing {len(matches)} top matches. Extracting flanking sequences...")
        # Build target sequences with flanks
        targets = []
        for i, (start, end, strand, pat, found_seq, rev_comp, mis, loc) in enumerate(matches):
            seq, coords = self._extract_flanked_window(start, end, strand, left_flank, right_flank)
            targets.append({
                'sequence': seq,
                'name': f"Match_{i+1}_{loc.replace('[','').replace(']','').replace('(','_').replace(')','')}",
                'coordinates': coords,
                'strand': strand,
            })
        # Display alignments (best last)
        self.display_alignment_pairwise(references, targets, top_n=top_n)

    def load_complete_sequence(self, filename: str):
        print(f"Loading complete genome sequence from {os.path.basename(filename)}...")
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                sequence_parts = [line.strip().upper() for line in f if not line.startswith('>')]
            self.sequence = "".join(sequence_parts)
            print(f"Successfully loaded genome of {len(self.sequence):,} bp.")
        except FileNotFoundError:
            print(f"Error: Genome file not found at {filename}"); sys.exit(1)
        except Exception as e:
            print(f"An error occurred loading the genome: {e}"); sys.exit(1)

    def load_annotations(self, filename: str):
        print(f"Loading gene annotations from {os.path.basename(filename)}...")
        try:
            all_genes = []
            with open(filename, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('>'):
                        header = line.strip()
                        loc_match = re.search(r'\[location=([^\]]+)\]', header)
                        if not loc_match: continue
                        coords = [int(n) for n in re.findall(r'\d+', loc_match.group(1))]
                        if not coords: continue
                        gene_name_match = re.search(r'\[gene=([^\]]+)\]', header)
                        locus_tag_match = re.search(r'\[locus_tag=([^\]]+)\]', header)
                        gene_info = {
                            'name': gene_name_match.group(1) if gene_name_match else "unknown",
                            'locus_tag': locus_tag_match.group(1) if locus_tag_match else "N/A",
                            'start': min(coords), 'end': max(coords),
                            'strand': '-' if "complement" in loc_match.group(1) else '+', 'type': 'gene'
                        }
                        all_genes.append(gene_info)

            all_genes.sort(key=lambda g: g['start'])
            last_coord = 0
            last_gene_name = 'start'
            for gene in all_genes:
                if gene['start'] == 0: continue
                intergenic_start, intergenic_end = last_coord + 1, gene['start'] - 1
                if intergenic_end > intergenic_start:
                    intergenic_info = {'name': f"intergenic", 'type': 'intergenic', 'upstream': last_gene_name, 'downstream': gene['name']}
                    self.gene_tree[intergenic_start:intergenic_end + 1] = intergenic_info
                self.gene_tree[gene['start']:gene['end'] + 1] = gene
                self.gene_data[gene.get('locus_tag', gene['name'])] = gene
                last_coord = gene['end']
                last_gene_name = gene['name']
            if last_coord < len(self.sequence):
                 self.gene_tree[last_coord + 1:len(self.sequence) + 1] = {'name': f"intergenic", 'type': 'intergenic', 'upstream': last_gene_name, 'downstream': 'end'}
            print(f"Successfully built map for {len(all_genes)} genes and intergenic regions.")
        except FileNotFoundError:
             print(f"Error: Annotation file not found at {filename}"); sys.exit(1)
        except Exception as e:
            print(f"An error occurred loading annotations: {e}"); sys.exit(1)

    def get_location_details(self, start: int, end: int) -> str:
        """
        CRITICAL FIX: Get precise location details for a match, properly identifying:
        - Exact gene matches
        - Gene overlaps with boundary distances  
        - Intergenic regions between genes
        """
        overlapping_regions = sorted(self.gene_tree[start:end])
        
        if not overlapping_regions:
            return "[unknown_region]"
        
        # Get the primary region (most overlap)
        primary_region = overlapping_regions[0].data
        region_type = primary_region.get('type', 'gene')
        
        if region_type == 'intergenic':
            # CRITICAL FIX: This is an intergenic region - report the genes it's between
            upstream = primary_region.get('upstream', 'start')
            downstream = primary_region.get('downstream', 'end')
            return f"[{upstream}][-------][{downstream}]"
        
        # Handle gene regions
        gene_name = primary_region['name']
        g_start = primary_region['start']
        g_end = primary_region['end']
        
        # Check if match is completely within the gene
        if start >= g_start and end <= g_end:
            return f"[{gene_name}][In]"
        
        # Check for overlaps with gene boundaries
        overlap_info = []
        
        if start < g_start:
            overlap_bp = g_start - start
            overlap_info.append(f"<{overlap_bp}bp>[{gene_name}]")
        
        if end > g_end:
            overlap_bp = end - g_end
            overlap_info.append(f"[{gene_name}]<{overlap_bp}bp>")
        
        if overlap_info:
            return "".join(overlap_info)
        
        # If match spans the entire gene
        if start < g_start and end > g_end:
            return f"<{g_start - start}bp>[{gene_name}]<{end - g_end}bp>"
        
        # Default case
        return f"[{gene_name}](Overlap)"
    
    def _get_intergenic_context(self, start: int, end: int) -> str:
        """
        CRITICAL FIX: Get precise intergenic context between genes.
        This method specifically handles cases like between CopZ and CsoR.
        """
        # Find the closest upstream and downstream genes
        upstream_gene = None
        downstream_gene = None
        
        # Get all genes sorted by position
        all_genes = []
        for gene_name, gene_info in self.gene_data.items():
            if 'start' in gene_info and 'end' in gene_info:
                all_genes.append((gene_name, gene_info['start'], gene_info['end']))
        
        if not all_genes:
            return "[start][-------][end]"
        
        # Sort genes by start position
        all_genes.sort(key=lambda x: x[1])
        
        # Find upstream gene (last gene that ends before our match starts)
        for gene_name, gene_start, gene_end in all_genes:
            if gene_end < start:
                upstream_gene = gene_name
            else:
                break
        
        # Find downstream gene (first gene that starts after our match ends)
        for gene_name, gene_start, gene_end in all_genes:
            if gene_start > end:
                downstream_gene = gene_name
                break
        
        # Handle edge cases
        if upstream_gene is None:
            upstream_gene = 'start'
        if downstream_gene is None:
            downstream_gene = 'end'
        
        # Return intergenic context
        return f"[{upstream_gene}][-------][{downstream_gene}]"

    def search_sequence(self, pattern: str, max_mismatches: int = 0) -> List[Tuple]:
        """
        CRITICAL FIX: Search for pattern in genome sequence with improved accuracy.
        Now properly handles IUPAC codes and intergenic region detection.
        """
        matches = []
        pattern = pattern.upper()
        seq = self.sequence
        
        if max_mismatches == 0:
            print("Using fast search mode (0 mismatches)...")
            
            # CRITICAL FIX: Improved IUPAC to regex conversion
            def iupac_to_regex(p): 
                return "".join(self.iupac_map.get(char, re.escape(char)) for char in p)
            
            # Search forward strand
            fwd_regex = iupac_to_regex(pattern)
            print(f"Searching forward strand for pattern: {pattern}")
            
            for match in tqdm(re.finditer(fwd_regex, seq), desc=f"Forward search", leave=False):
                start, end, window = match.start() + 1, match.end(), match.group(0)
                # CRITICAL FIX: Use improved location detection
                location = self.get_location_details(start, end)
                matches.append((start, end, '+', pattern, window, self._reverse_complement(window), 0, location))
            
            # Search reverse strand
            rev_pattern = self._reverse_complement(pattern)
            if rev_pattern != pattern:
                print(f"Searching reverse strand for pattern: {rev_pattern}")
                rev_regex = iupac_to_regex(rev_pattern)
                
                for match in tqdm(re.finditer(rev_regex, seq), desc=f"Reverse search", leave=False):
                    start, end, window = match.start() + 1, match.end(), match.group(0)
                    # CRITICAL FIX: Use improved location detection and report correct sequence
                    location = self.get_location_details(start, end)
                    # Report the reverse complement of what was found in the genome
                    rev_found_seq = self._reverse_complement(window)
                    matches.append((start, end, '-', pattern, rev_found_seq, self._reverse_complement(rev_found_seq), 0, location))
        else:
            print(f"Using mismatch search mode ({max_mismatches} mismatches allowed)...")
            pattern_len = len(pattern)
            rev_pattern = self._reverse_complement(pattern)
            matches_found_counter = 0
            
            # CRITICAL FIX: Improved progress reporting
            progress_bar = tqdm(range(len(seq) - pattern_len + 1), desc=f"Searching {pattern[:15]}...", leave=False)
            
            for i in progress_bar:
                window = seq[i:i+pattern_len]
                
                # Search forward strand
                mismatches_fwd = self._count_mismatches(pattern, window)
                if mismatches_fwd <= max_mismatches:
                    start, end = i + 1, i + pattern_len
                    # CRITICAL FIX: Use improved location detection
                    location = self.get_location_details(start, end)
                    matches.append((start, end, '+', pattern, window, self._reverse_complement(window), mismatches_fwd, location))
                    matches_found_counter += 1
                    progress_bar.set_postfix_str(f"Matches: {matches_found_counter}")
                
                # Search reverse strand
                mismatches_rev = self._count_mismatches(rev_pattern, window)
                if mismatches_rev <= max_mismatches and rev_pattern != pattern:
                    start, end = i + 1, i + pattern_len
                    # CRITICAL FIX: Use improved location detection and report correct sequence
                    location = self.get_location_details(start, end)
                    # Report the reverse complement of what was found in the genome
                    rev_found_seq = self._reverse_complement(window)
                    matches.append((start, end, '-', pattern, rev_found_seq, self._reverse_complement(rev_found_seq), mismatches_rev, location))
                    matches_found_counter += 1
                    progress_bar.set_postfix_str(f"Matches: {matches_found_counter}")
        
        # CRITICAL FIX: Improved duplicate removal and sorting
        print(f"Processing {len(matches)} raw matches...")
        
        # Remove duplicates based on position, strand, and pattern
        unique_matches = []
        seen_positions = set()
        
        for match in matches:
            # Create unique key for each match
            position_key = (match[0], match[1], match[2], match[3])  # start, end, strand, pattern
            if position_key not in seen_positions:
                unique_matches.append(match)
                seen_positions.add(position_key)
        
        # Sort by start position
        unique_matches.sort(key=lambda x: x[0])
        
        print(f"Found {len(unique_matches)} unique matches after deduplication")
        return unique_matches

def parse_fasta_file(filepath: str) -> Dict[str, str]:
    sequences = {}
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            header = ""
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    header = line
                    sequences[header] = ""
                elif header:
                    sequences[header] += line.upper()
    except FileNotFoundError:
        print(f"Error: Could not find file {filepath}")
        return {}
    return sequences

def save_results_to_csv(matches: List, filepath: str):
    print(f"\nSaving results to {filepath}...")
    try:
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            header = ["Start", "End", "Strand", "Pattern", "Found Sequence", "RevComp", "Mismatches", "Location"]
            writer.writerow(header)
            for match_tuple in matches:
                writer.writerow(list(match_tuple))
        print("Save complete.")
    except Exception as e:
        print(f"An error occurred while saving the file: {e}")

def main():
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        script_dir = os.getcwd()
    
    project_root = os.path.dirname(script_dir)
    data_dir = os.path.join(project_root, "data")
    seq_files_dir = os.path.join(script_dir, "Seq_Files")
    results_dir = os.path.join(script_dir, "Results")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(seq_files_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    complete_genome_file = os.path.join(data_dir, "E.Coli_K12_Seq.fasta")
    gene_annotation_file = os.path.join(data_dir, "Gene_Coordiantes.txt")
    cache_file = os.path.join(data_dir, "genome_cache.pkl")
    
    analyzer = EColiAnalyzer()

    if os.path.exists(cache_file):
        print(f"Loading cached genome data from: {os.path.basename(cache_file)}...")
        with open(cache_file, 'rb') as f:
            cached_data = pickle.load(f)
        analyzer.sequence = cached_data['sequence']
        analyzer.gene_tree = cached_data['gene_tree']
        analyzer.gene_data = cached_data['gene_data']
        print("Cached data loaded.")
    else:
        analyzer.load_complete_sequence(complete_genome_file)
        analyzer.load_annotations(gene_annotation_file)
        with open(cache_file, 'wb') as f:
            pickle.dump({
                'sequence': analyzer.sequence,
                'gene_tree': analyzer.gene_tree,
                'gene_data': analyzer.gene_data
            }, f)
        print(f"Saved new genome data to cache file: {os.path.basename(cache_file)}")

    try:
        while True:
            print("\nChoose search mode:")
            print("1. Single search in E. coli genome")
            print("2. Batch search from substrate file")
            print("3. Quit")
            mode = input("Enter mode (1-3): ")

            if mode == '1':
                pattern = input("Enter sequence to search (IUPAC supported): ").strip().upper()
                if not pattern: continue
                mismatches_str = input("Enter maximum allowed mismatches (default 0): ")
                max_mismatches = int(mismatches_str) if mismatches_str.isdigit() else 0
                matches = analyzer.search_sequence(pattern, max_mismatches)
                print(f"\nTotal matches found: {len(matches)}")
                if not matches: continue
                print("\nSearch Results:")
                print("-" * 180)
                print(f"{'Start':<9}{'End':<9}{'Strand':<9}{'Pattern':<35}{'Found Sequence':<35}{'RevComp':<35}{'Mismatches':<12}{'Location'}")
                print("-" * 180)
                for start, end, strand, pat, fwd, rev, mis, loc in matches:
                    print(f"{start:<9}{end:<9}{strand:<9}{pat[:32]:<35}{fwd[:32]:<35}{rev[:32]:<35}{mis:<12}{loc}")
                print("-" * 180)

            elif mode == '2':
                try:
                    available_files = [f for f in os.listdir(seq_files_dir) if f.endswith(('.fasta', '.txt', '.fa'))]
                    if not available_files:
                        print(f"\nNo files found in '{seq_files_dir}'."); continue
                    print("\nAvailable substrate files:")
                    for i, filename in enumerate(available_files): print(f"{i + 1}. {filename}")
                    choice = int(input(f"Choose a file for batch search (1-{len(available_files)}): ")) - 1
                    if not 0 <= choice < len(available_files):
                        print("Invalid choice."); continue
                    selected_filename = available_files[choice]
                    selected_filepath = os.path.join(seq_files_dir, selected_filename)
                    queries = parse_fasta_file(selected_filepath)
                    if not queries:
                        print(f"No sequences found in {selected_filename}."); continue
                    mismatches_str = input("Enter maximum allowed mismatches for all searches (default 0): ")
                    max_mismatches = int(mismatches_str) if mismatches_str.isdigit() else 0
                    all_matches = []
                    print(f"\nStarting batch search for {len(queries)} patterns...")
                    for header, pattern in queries.items():
                        print(f"\n--- Searching for pattern from: {header} ---")
                        matches = analyzer.search_sequence(pattern, max_mismatches)
                        print(f"Found {len(matches)} total match(es) for this pattern.")
                        if matches:
                            print("\nIndividual Search Results (Sorted by Position):")
                            print("-" * 180)
                            print(f"{'Start':<9}{'End':<9}{'Strand':<9}{'Pattern':<35}{'Found Sequence':<35}{'RevComp':<35}{'Mismatches':<12}{'Location'}")
                            print("-" * 180)
                            for match_tuple in matches:
                                start, end, strand, pat, fwd, rev, mis, loc = match_tuple
                                print(f"{start:<9}{end:<9}{strand:<9}{pat[:32]:<35}{fwd[:32]:<35}{rev[:32]:<35}{mis:<12}{loc}")
                            print("-" * 180)
                        for match in matches:
                            original_location = match[-1]
                            full_match_info = list(match)
                            full_match_info[-1] = f"{original_location} (Query: {header[:30]}...)"
                            all_matches.append(tuple(full_match_info))
                    print(f"\n--- Batch search complete. Found {len(all_matches)} total matches across all queries. ---")
                    if not all_matches: continue
                    save_prompt = input("\nDo you want to save the complete results to a CSV file? (y/n): ").lower()
                    if save_prompt == 'y':
                        base_name, _ = os.path.splitext(selected_filename)
                        output_filename = f"{base_name}_Results.csv"
                        output_filepath = os.path.join(results_dir, output_filename)
                        save_results_to_csv(all_matches, output_filepath)
                except (ValueError, IndexError):
                    print("Invalid input.")
                except Exception as e:
                    print(f"An error occurred: {e}")

            elif mode == '3':
                print("Goodbye!")
                break
            else:
                print("Invalid mode.")
    except KeyboardInterrupt:
        print("\n\nProgram interrupted by user. Exiting gracefully.")
        sys.exit(0)

if __name__ == "__main__":
    main()