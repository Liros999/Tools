import re
from typing import List, Tuple, Optional
from gene_coordinates import GeneData
import os
from tqdm import tqdm
from tss_distance_analyzer import TSSDistanceAnalyzer

class BSubAnalyzer:
    def __init__(self):
        self.sequence = ""
        self.gene_data = GeneData()
        # Initialize TSS distance analyzer
        try:
            self.tss_analyzer = TSSDistanceAnalyzer()
            self.tss_available = True
            print("TSS distance analysis enabled.")
        except Exception as e:
            print(f"Warning: TSS distance analysis not available: {e}")
            self.tss_analyzer = None
            self.tss_available = False
        
    def load_sequence(self, filename: str) -> bool:
        """Load B. subtilis sequence from file"""
        try:
            with open(filename, 'r') as f:
                # Skip header if it exists
                first_line = f.readline()
                if not first_line.startswith('>'):
                    f.seek(0)  # Go back to start if no header
                self.sequence = f.read().replace('\n', '').upper()
            return True
        except Exception as e:
            print(f"Error loading sequence: {e}")
            return False
    
    def _iupac_to_regex(self, pattern: str) -> str:
        """Convert IUPAC codes to regex pattern"""
        iupac_map = {
            'R': '[AG]',    # A or G
            'Y': '[CT]',    # C or T
            'S': '[GC]',    # G or C
            'W': '[AT]',    # A or T
            'K': '[GT]',    # G or T
            'M': '[AC]',    # A or C
            'B': '[CGT]',   # C or G or T
            'D': '[AGT]',   # A or G or T
            'H': '[ACT]',   # A or C or T
            'V': '[ACG]',   # A or C or G
            'N': '[ACGT]'   # Any base
        }
        
        # Convert pattern to uppercase
        pattern = pattern.upper()
        
        # Replace IUPAC codes with regex patterns
        for code, regex in iupac_map.items():
            pattern = pattern.replace(code, regex)
            
        return pattern

    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of a DNA sequence"""
        complement = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
            'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N'
        }
        return ''.join(complement[base] for base in reversed(seq))

    def _count_mismatches(self, pattern_seq: str, genome_seq: str) -> int:
        """
        Count mismatches between a pattern (may include IUPAC codes) and a genome window.
        - If the pattern base is a concrete nucleotide (A/T/C/G), it must equal the genome base.
        - If the pattern base is an IUPAC code, any base in its allowed set is accepted without penalty.
        """
        iupac_allowed = {
            'A': {'A'}, 'T': {'T'}, 'C': {'C'}, 'G': {'G'},
            'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
            'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
            'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}
        }
        mismatches = 0
        pattern_seq = pattern_seq.upper()
        genome_seq = genome_seq.upper()
        for pattern_base, genome_base in zip(pattern_seq, genome_seq):
            allowed = iupac_allowed.get(pattern_base, {pattern_base})
            if genome_base not in allowed:
                mismatches += 1
        return mismatches

    def _count_palindromic_mismatches(self, seq: str) -> int:
        """Count palindromic mismatches: number of positions where base is not the complement of the symmetric base."""
        complement = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
            'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N'
        }
        seq = seq.upper()
        mismatches = 0
        for i in range(len(seq) // 2):
            left = seq[i]
            right = seq[-(i+1)]
            if complement.get(left, 'N') != right:
                mismatches += 1
        return mismatches

    def _is_palindrome(self, seq: str, max_mismatches: int = 0) -> bool:
        """Check if sequence is a palindrome with allowed palindromic mismatches"""
        return self._count_palindromic_mismatches(seq) <= max_mismatches

    def search_sequence(self, pattern: str, max_mismatches: int = 0) -> List[Tuple[int, int, str, str, str, str, int]]:
        """
        CRITICAL FIX: Search for pattern and its reverse complement in the genome, allowing for mismatches.
        Returns list of (start, end, strand, pattern, found_sequence, location, mismatches)
        """
        matches = []
        pattern = pattern.upper()
        pattern_len = len(pattern)
        seq = self.sequence.upper()
        rev_pattern = self._reverse_complement(pattern)
        
        if max_mismatches == 0:
            # CRITICAL FIX: Use improved regex for exact matches with proper IUPAC handling
            print("\nSearching forward strand (exact match)...")
            regex_pattern = self._iupac_to_regex(pattern)
            
            for match in tqdm(re.finditer(regex_pattern, seq), desc="Progress"):
                start = match.start() + 1  # Convert to 1-based coordinates
                end = match.end()
                found_seq = match.group()
                location = self._get_match_location(start, end, '+')
                matches.append((start, end, '+', pattern, found_seq, location, 0))
            
            # CRITICAL FIX: Always search reverse strand, even if rev_pattern == pattern
            print("\nSearching reverse strand (exact match)...")
            rev_regex_pattern = self._iupac_to_regex(rev_pattern)
            
            for match in tqdm(re.finditer(rev_regex_pattern, seq), desc="Progress"):
                start = match.start() + 1  # Convert to 1-based coordinates
                end = match.end()
                found_seq = match.group()
                # CRITICAL FIX: Report the reverse complement of what was found
                rev_found_seq = self._reverse_complement(found_seq)
                location = self._get_match_location(start, end, '-')
                matches.append((start, end, '-', pattern, rev_found_seq, location, 0))
        else:
            # CRITICAL FIX: Improved mismatch search with better progress reporting
            print(f"\nSearching with {max_mismatches} mismatches allowed...")
            matches_found = 0
            
            # Search forward strand
            print("Searching forward strand...")
            for i in tqdm(range(len(seq) - pattern_len + 1), desc="Forward strand"):
                window = seq[i:i + pattern_len]
                mismatches = self._count_mismatches(pattern, window)
                
                if mismatches <= max_mismatches:
                    start = i + 1  # Convert to 1-based coordinates
                    end = i + pattern_len
                    location = self._get_match_location(start, end, '+')
                    matches.append((start, end, '+', pattern, window, location, mismatches))
                    matches_found += 1
            
            # Search reverse strand (always scan)
            print("Searching reverse strand...")
            for i in tqdm(range(len(seq) - pattern_len + 1), desc="Reverse strand"):
                window = seq[i:i + pattern_len]
                mismatches = self._count_mismatches(rev_pattern, window)
                
                if mismatches <= max_mismatches:
                    start = i + 1  # Convert to 1-based coordinates
                    end = i + pattern_len
                    # CRITICAL FIX: Report the reverse complement of what was found
                    rev_found_seq = self._reverse_complement(window)
                    location = self._get_match_location(start, end, '-')
                    matches.append((start, end, '-', pattern, rev_found_seq, location, mismatches))
                    matches_found += 1
            
            print(f"Total matches found: {matches_found}")
        
        # CRITICAL FIX: Remove duplicates and sort by position
        unique_matches = []
        seen_positions = set()
        
        for match in matches:
            position_key = (match[0], match[1], match[2])  # start, end, strand
            if position_key not in seen_positions:
                unique_matches.append(match)
                seen_positions.add(position_key)
        
        # Sort by start position
        unique_matches.sort(key=lambda x: x[0])
        
        return unique_matches

    def search_palindrome(self, pattern: str, max_mismatches: int = 0) -> List[Tuple[int, int, str, str, str, str, int]]:
        """
        Search for palindromic sequences with allowed mismatches
        Returns list of (start, end, strand, pattern, found_sequence, location, mismatches)
        """
        regex_pattern = self._iupac_to_regex(pattern)
        matches = []
        # Forward strand
        print("\nSearching forward strand...")
        for match in tqdm(re.finditer(regex_pattern, self.sequence), desc="Progress"):
            start = match.start()
            end = match.end()
            found_seq = match.group()
            mismatches = self._count_mismatches(found_seq, pattern)
            if self._is_palindrome(found_seq, max_mismatches):
                # Check if sequence is within any gene
                location = self._get_match_location(start, end, '+')
                
                # Determine strand based on gene context
                strand = '+'
                genes = self.gene_data.find_genes_for_coordinates(start, end)
                if genes:
                    gene_info = self.gene_data.get_gene_info(genes[0])
                    if gene_info and gene_info.get('strand') == '-':
                        strand = '-'
                
                matches.append((start, end, strand, pattern, found_seq, location, mismatches))
        
        # Search reverse strand
        print("\nSearching reverse strand...")
        rev_pattern = self._reverse_complement(pattern)
        rev_regex = self._iupac_to_regex(rev_pattern)
        rev_sequence = self._reverse_complement(self.sequence)
        
        for match in tqdm(re.finditer(rev_regex, rev_sequence), desc="Progress"):
            # Convert coordinates back to original sequence
            start = len(self.sequence) - match.end()
            end = len(self.sequence) - match.start()
            found_seq = self._reverse_complement(match.group())
            mismatches = self._count_mismatches(found_seq, pattern)
            
            # Check if it's a palindrome with allowed mismatches
            if self._is_palindrome(found_seq, max_mismatches):
                # Check if sequence is within any gene
                location = self._get_match_location(start, end, '-')
                
                # Determine strand based on gene context
                strand = '-'
                genes = self.gene_data.find_genes_for_coordinates(start, end)
                if genes:
                    gene_info = self.gene_data.get_gene_info(genes[0])
                    if gene_info and gene_info.get('strand') == '+':
                        strand = '+'
                
                matches.append((start, end, strand, pattern, found_seq, location, mismatches))
        # Sort matches by start position
        matches.sort(key=lambda x: x[0])
        return matches

    def search_region_palindromes(self, start_pos: int, end_pos: int, 
                                min_length: int = 6, max_length: int = 30,
                                max_mismatches: int = 0) -> List[Tuple[int, int, str, str, str, str, int]]:
        """
        Search for palindromes in a specific region
        Returns list of (start, end, strand, pattern, found_sequence, location, mismatches)
        """
        matches = []
        # Use 1-based inclusive coordinates for biological convention
        region = self.sequence[start_pos-1:end_pos]
        
        # Search for palindromes of different lengths
        print("\nSearching for palindromes...")
        for length in tqdm(range(min_length, min(max_length + 1, len(region) + 1)), desc="Length progress"):
            # Search forward strand
            print(f"\nChecking length {length} (forward strand)...")
            for i in tqdm(range(len(region) - length + 1), desc="Position progress"):
                seq = region[i:i + length]
                mismatches = self._count_palindromic_mismatches(seq)
                if self._is_palindrome(seq, max_mismatches):
                    # Check if sequence is within any gene
                    abs_start = start_pos - 1 + i + 1  # Convert back to 1-based
                    abs_end = abs_start + length - 1   # Inclusive end
                    location = self._get_match_location(abs_start, abs_end, '+')
                    
                    # Determine strand
                    strand = '+'
                    genes = self.gene_data.find_genes_for_coordinates(abs_start, abs_end)
                    if genes:
                        gene_info = self.gene_data.get_gene_info(genes[0])
                        if gene_info and gene_info.get('strand') == '-':
                            strand = '-'
                    
                    matches.append((abs_start, abs_end, strand, "Palindrome", seq, location, mismatches))
            
            # Search reverse strand
            print(f"\nChecking length {length} (reverse strand)...")
            rev_region = self._reverse_complement(region)
            for i in tqdm(range(len(rev_region) - length + 1), desc="Position progress"):
                seq = rev_region[i:i + length]
                mismatches = self._count_palindromic_mismatches(seq)
                if self._is_palindrome(seq, max_mismatches):
                    # Convert coordinates back to original sequence
                    abs_start = end_pos - (i + length) + 1  # 1-based
                    abs_end = abs_start + length - 1        # Inclusive end
                    location = self._get_match_location(abs_start, abs_end, '-')
                    
                    # Determine strand
                    strand = '-'
                    genes = self.gene_data.find_genes_for_coordinates(abs_start, abs_end)
                    if genes:
                        gene_info = self.gene_data.get_gene_info(genes[0])
                        if gene_info and gene_info.get('strand') == '+':
                            strand = '+'
                    
                    matches.append((abs_start, abs_end, strand, "Palindrome", seq, location, mismatches))
        
        return matches

    def save_results(self, results: List[Tuple[int, int, str, str, str, str, int]], filename: str):
        """Save search results to CSV file in the results folder with TSS distance analysis"""
        results_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results')
        os.makedirs(results_dir, exist_ok=True)
        filepath = os.path.join(results_dir, filename)
        
        # Check if TSS analysis is available
        if self.tss_available and self.tss_analyzer:
            # Enhanced results with TSS distance
            with open(filepath, 'w', newline='') as f:
                f.write("Start,End,Strand,Pattern,Found_Sequence,Location,Distance_from_TSS,Mismatches\n")
                for start, end, strand, pattern, found_seq, location, mismatches in results:
                    # Calculate TSS distance
                    tss_distance = self.tss_analyzer.calculate_tss_distance_string(start, end, strand)
                    f.write(f"{start},{end},{strand},{pattern},{found_seq},{location},{tss_distance},{mismatches}\n")
            print(f"Results with TSS distance analysis saved to {filepath}")
        else:
            # Standard results without TSS analysis
            with open(filepath, 'w', newline='') as f:
                f.write("Start,End,Strand,Pattern,Found_Sequence,Location,Mismatches\n")
                for start, end, strand, pattern, found_seq, location, mismatches in results:
                    f.write(f"{start},{end},{strand},{pattern},{found_seq},{location},{mismatches}\n")
            print(f"Results saved to {filepath} (TSS analysis not available)")

    def _is_match_in_gene(self, start: int, end: int) -> bool:
        """
        Check if a match is completely inside a gene
        Returns True if the match is inside a gene, False if it's between genes
        """
        genes = self.gene_data.find_genes_for_coordinates(start, end)
        return len(genes) > 0

    def _get_match_location(self, start: int, end: int, strand: str) -> str:
        """
        Get the location string for a match, properly identifying if it's in a gene or between genes.
        This is the CRITICAL FIX for proper intergenic region detection.
        """
        # First, check if match is completely within any gene
        genes = self.gene_data.find_genes_for_coordinates(start, end)
        
        # Filter genes to only those on the same strand as the match
        same_strand_genes = []
        for gene in genes:
            gene_info = self.gene_data.get_gene_info(gene)
            if gene_info and normalize_strand(gene_info.get('strand')) == normalize_strand(strand):
                same_strand_genes.append(gene)
        
        if same_strand_genes:
            # Match is inside gene(s) on the same strand
            if len(same_strand_genes) > 1:
                return f"[Overlap]({same_strand_genes[0]})"
            else:
                gene = same_strand_genes[0]
                return f"[{gene}](in)"
        else:
            # Match is between genes (intergenic) - this is the CRITICAL FIX
            return self._get_intergenic_location(start, end, strand)
    
    def _get_intergenic_location(self, start: int, end: int, strand: str) -> str:
        """
        CRITICAL FIX: Properly identify intergenic regions between genes.
        This method specifically handles the case between CopZ and CsoR.
        """
        # Normalize strand
        strand = normalize_strand(strand)
        
        # Get all genes on the same strand, sorted by position
        genes = [
            (g, d["start"], d["end"]) 
            for g, d in self.gene_data.genes.items() 
            if d.get("start") is not None and d.get("end") is not None and normalize_strand(d.get("strand")) == strand
        ]
        
        if not genes:
            return "[start]---[end]"
        
        # Sort genes by start position
        genes.sort(key=lambda x: x[1])
        
        # Find the specific intergenic region where this match is located
        upstream_gene = None
        downstream_gene = None
        
        for i, (gene_name, gene_start, gene_end) in enumerate(genes):
            # Check if match is in the intergenic region after this gene
            if gene_end < start:
                # This gene is upstream of the match
                upstream_gene = gene_name
                
                # Check if there's a downstream gene
                if i + 1 < len(genes):
                    next_gene_name, next_gene_start, next_gene_end = genes[i + 1]
                    if next_gene_start > end:
                        # Match is in intergenic region between these two genes
                        downstream_gene = next_gene_name
                        break
                else:
                    # No more genes downstream
                    downstream_gene = 'end'
                    break
        
        # If we didn't find a downstream gene, check if we're before the first gene
        if upstream_gene is None and genes:
            first_gene_name, first_gene_start, first_gene_end = genes[0]
            if start < first_gene_start:
                upstream_gene = 'start'
                downstream_gene = first_gene_name
        
        # If still no upstream gene, check if we're after the last gene
        if upstream_gene is None and genes:
            last_gene_name, last_gene_start, last_gene_end = genes[-1]
            if end > last_gene_end:
                upstream_gene = last_gene_name
                downstream_gene = 'end'
        
        # Calculate distances for precise intergenic reporting
        if upstream_gene and downstream_gene:
            if upstream_gene == 'start':
                # Match is at the beginning of the genome
                return f"[{upstream_gene}]---[{downstream_gene}]"
            elif downstream_gene == 'end':
                # Match is at the end of the genome
                return f"[{upstream_gene}]---[{downstream_gene}]"
            else:
                # Match is between two genes - this is the CRITICAL CASE for CopZ-CsoR
                return f"[{upstream_gene}]---[{downstream_gene}]"
        else:
            # Fallback to generic intergenic reporting
            return "[intergenic_region]"

    def get_gene_context(self, start: int, end: int, strand: str) -> str:
        """
        Return context string: [upstream_gene]---[downstream_gene] based on match strand.
        Only consider genes on the same strand as the match.
        """
        # Normalize strand
        strand = normalize_strand(strand)
        
        # Prepare gene list (only same strand)
        genes = [
            (g, d["start"], d["end"]) 
            for g, d in self.gene_data.genes.items() 
            if d.get("start") is not None and d.get("end") is not None and normalize_strand(d.get("strand")) == strand
        ]
        
        if not genes:
            return "[start]---[end]"
        
        # Sort genes by start position
        genes.sort(key=lambda x: x[1])
        
        # Find the closest upstream and downstream genes
        upstream = None
        downstream = None
        
        for gene_name, gene_start, gene_end in genes:
            # Check if the match is completely between this gene and the next one
            if gene_end < start:  # Gene is upstream
                upstream = gene_name
            elif gene_start > end:  # Gene is downstream
                downstream = gene_name
                break  # Found the first downstream gene, no need to continue
        
        # If no upstream gene found, check if there are genes before this position
        if upstream is None:
            for gene_name, gene_start, gene_end in reversed(genes):
                if gene_end < start:
                    upstream = gene_name
                    break
        
        # If no downstream gene found, check if there are genes after this position
        if downstream is None:
            for gene_name, gene_start, gene_end in genes:
                if gene_start > end:
                    downstream = gene_name
                    break
        
        up_name = upstream if upstream else 'start'
        down_name = downstream if downstream else 'end'
        
        return f"[{up_name}]---[{down_name}]"

    def _get_gene_context_any_strand(self, start: int, end: int) -> str:
        """
        Return context string without filtering by strand:
        [upstream_gene]---[downstream_gene]
        """
        genes = [
            (g, d["start"], d["end"]) 
            for g, d in self.gene_data.genes.items() 
            if d.get("start") is not None and d.get("end") is not None
        ]
        if not genes:
            return "[start]---[end]"
        genes.sort(key=lambda x: x[1])
        upstream = None
        downstream = None
        for gene_name, gene_start, gene_end in genes:
            if gene_end < start:
                upstream = gene_name
            elif gene_start > end:
                downstream = gene_name
                break
        if upstream is None:
            for gene_name, gene_start, gene_end in reversed(genes):
                if gene_end < start:
                    upstream = gene_name
                    break
        if downstream is None:
            for gene_name, gene_start, gene_end in genes:
                if gene_start > end:
                    downstream = gene_name
                    break
        up_name = upstream if upstream else 'start'
        down_name = downstream if downstream else 'end'
        return f"[{up_name}]---[{down_name}]"

    def _get_match_location_any_strand(self, start: int, end: int) -> str:
        """
        Get location without strand filtering. If inside gene(s), report the first gene with (in),
        else report unstranded intergenic context.
        """
        genes = self.gene_data.find_genes_for_coordinates(start, end)
        if genes:
            gene = genes[0]
            return f"[{gene}](in)"
        return self._get_gene_context_any_strand(start, end)

    def analyze_results(self, results: List[Tuple[int, int, str, str, str, str, int]], 
                       filter_type: str, filter_value: str) -> List[Tuple[int, int, str, str, str, str, int]]:
        """Filter results based on user criteria"""
        filtered_results = []
        
        if filter_type == "1":  # Filter by regions
            try:
                start, end = map(int, filter_value.split(','))
                filtered_results = [r for r in results if start <= r[0] <= end or start <= r[1] <= end]
            except ValueError:
                print("Invalid region format. Please use format: start,end")
                return results
                
        elif filter_type == "2":  # Filter by mismatches
            try:
                max_mismatches = int(filter_value)
                filtered_results = [r for r in results if r[6] <= max_mismatches]
            except ValueError:
                print("Invalid mismatch number")
                return results
                
        elif filter_type == "3":  # Filter by location type
            if filter_value.lower() == "in":
                filtered_results = [r for r in results if "(in)" in r[5]]
            elif filter_value.lower() == "between":
                filtered_results = [r for r in results if "(in)" not in r[5]]
            else:
                print("Invalid location type. Use 'in' or 'between'")
                return results
                
        return filtered_results if filtered_results else results

    def analyze_palindrome_context(self, results: List[Tuple[int, int, str, str, str, str, int]], 
                                 context_size: int, max_mismatches: int = 0) -> List[Tuple[int, int, str, str, str, str, int, bool, str]]:
        """Analyze if sequences are part of palindromes in their context"""
        analyzed_results = []
        
        for result in results:
            start, end, strand, pattern, found_seq, location, mismatches = result
            # Get context sequence
            context_start = max(0, start - context_size)
            context_end = min(len(self.sequence), end + context_size)
            context_seq = self.sequence[context_start:context_end]
            
            # Check if sequence is part of a palindrome
            is_palindrome = self._is_palindrome(context_seq, max_mismatches)
            palindrome_seq = context_seq if is_palindrome else ""
            
            analyzed_results.append((*result, is_palindrome, palindrome_seq))
            
        return analyzed_results

    def perform_alignment_analysis(self):
        """Perform ClustalW-style alignment with flanking sequences"""
        print("\nSequence Alignment with Flanking Regions")
        
        # 1) Target sequence to find in genome
        pattern = input("Enter the target sequence pattern: ").strip().upper()
        if not pattern:
            print("Please enter a non-empty sequence pattern.")
            return
        
        # Check for invalid characters
        valid_chars = set('ATCGRYSWKMBDHVN')
        invalid_chars = set(pattern) - valid_chars
        
        if invalid_chars:
            print(f"Invalid characters found: {', '.join(sorted(invalid_chars))}")
            print("Valid characters are: A, T, C, G, R, Y, S, W, K, M, B, D, H, V, N")
            print("IUPAC codes: R=AG, Y=CT, S=GC, W=AT, K=GT, M=AC, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT")
            return

        # 2) Allowed mismatches for target search
        max_mismatches = get_int_input("Enter maximum allowed mismatches for search (default 0): ", default=0)

        # Search genome using the target pattern with allowed mismatches (not the refs)
        matches = self.search_sequence(pattern, max_mismatches)

        # Optional filtering (same UX as other features) BEFORE any flanks or alignment
        if matches:
            filter_option = input("\nWould you like to filter results? (y/n): ").lower().strip()
            if filter_option == 'y':
                print("\nFilter options:")
                print("1. Show only matches inside genes")
                print("2. Show only matches between genes (intergenic)")
                print("3. Show only matches with specific number of mismatches")
                print("4. Show only matches in specific region")
                filter_choice = input("Enter filter choice (1-4): ").strip()
                if filter_choice == "1":
                    print("\nApplying filter: inside genes...")
                    filtered = []
                    for m in tqdm(matches, desc="Filtering", leave=False):
                        if "(in)" in m[5]:
                            filtered.append(m)
                    matches = filtered or matches
                elif filter_choice == "2":
                    print("\nApplying filter: intergenic matches...")
                    filtered = []
                    for m in tqdm(matches, desc="Filtering", leave=False):
                        # Trust the location string from search: intergenic lacks '(in)'
                        if "(in)" not in m[5]:
                            filtered.append(m)
                    matches = filtered or matches
                elif filter_choice == "3":
                    try:
                        k = int(input("Enter maximum number of mismatches: "))
                        print("\nApplying filter: mismatches ≤ ", k)
                        filtered = []
                        for m in tqdm(matches, desc="Filtering", leave=False):
                            if m[6] <= k:
                                filtered.append(m)
                        matches = filtered or matches
                    except ValueError:
                        print("Invalid input. Skipping filter.")
                elif filter_choice == "4":
                    try:
                        region_input = input("Enter region (start,end): ")
                        rstart, rend = map(int, region_input.split(','))
                        print(f"\nApplying filter: region {rstart}-{rend}...")
                        filtered = []
                        for m in tqdm(matches, desc="Filtering", leave=False):
                            if (rstart <= m[0] <= rend or rstart <= m[1] <= rend or m[0] <= rstart <= m[1] or m[0] <= rend <= m[1]):
                                filtered.append(m)
                        matches = filtered or matches
                    except ValueError:
                        print("Invalid region format. Skipping filter.")

        if not matches:
            print("No matches found for the pattern.")
            return

        # 3) Flanking sizes to extract around filtered hits
        left_flank = get_int_input("Enter number of nucleotides to include from the left: ", default=50)
        right_flank = get_int_input("Enter number of nucleotides to include from the right: ", default=50)

        # 4) Ask how many top results to align (limit workload)
        try:
            top_n_input = input("\nHow many top results do you want to align? (default 100): ").strip()
            top_n = int(top_n_input) if top_n_input else 100
        except ValueError:
            top_n = 100

        print(f"\nPreparing matches. Extracting flanking sequences...")
        
        # Extract flanking sequences
        target_sequences = []
        for i, (start, end, strand, pat, found_seq, location, mismatches) in enumerate(matches):
            # Calculate flanking region coordinates (start/end are 1-based inclusive)
            # For reverse strand, swap left/right flanks in genomic coordinates
            if strand == '-':
                flank_upstream = right_flank  # "left" relative to motif becomes downstream in genome
                flank_downstream = left_flank
            else:
                flank_upstream = left_flank
                flank_downstream = right_flank

            # Convert to 0-based half-open slice: [slice_start, slice_end)
            slice_start = max(0, (start - 1) - flank_upstream)
            slice_end = min(len(self.sequence), end + flank_downstream)
            
            # Extract sequence with flanks
            flanked_seq = self.sequence[slice_start:slice_end]
            
            # If on reverse strand, take reverse complement
            if strand == '-':
                flanked_seq = self._reverse_complement(flanked_seq)
            
            # Use the ORIGINAL search-time location string for naming and filtering traceability
            pretty_loc = location.replace('[', '').replace(']', '').replace('(', '_').replace(')', '')
            target_sequences.append({
                'sequence': flanked_seq,
                'name': f"Match_{i+1}_{pretty_loc}",
                'coordinates': f"{slice_start}-{slice_end}",
                'strand': strand,
                'mismatches': mismatches,
                'orig_location': location
            })
        
        # 5) Reference sequences (unchanged, no flanks)
        ref_sequences = []
        print("\nEnter reference sequences (one per line, empty line to finish):")
        while True:
            ref_seq = input("Reference sequence: ").strip().upper()
            if not ref_seq:
                break
            valid_chars = set('ATCGRYSWKMBDHVN')
            invalid_chars = set(ref_seq) - valid_chars
            if invalid_chars:
                print(f"Invalid characters found: {', '.join(sorted(invalid_chars))}")
                print("Valid characters are: A, T, C, G, R, Y, S, W, K, M, B, D, H, V, N")
                print("IUPAC codes: R=AG, Y=CT, S=GC, W=AT, K=GT, M=AC, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT")
                continue
            ref_sequences.append(ref_seq)
        if not ref_sequences:
            print("No reference sequences provided.")
            return

        # 6) Rank targets by best-window score vs the first reference; keep top_n
        primary_ref = ref_sequences[0]
        scored_targets = []
        win_len = len(primary_ref)
        for t in target_sequences:
            seq = t['sequence']
            best_mm = None
            for i in range(0, max(1, len(seq) - win_len + 1)):
                mm = self._count_mismatches(primary_ref, seq[i:i+win_len])
                if best_mm is None or mm < best_mm:
                    best_mm = mm
            t['best_mm'] = 0 if best_mm is None else best_mm
            scored_targets.append(t)
        scored_targets.sort(key=lambda x: (x.get('best_mm', 0), x.get('name', '')))
        if len(scored_targets) > top_n:
            scored_targets = scored_targets[:top_n]

        # 7) Perform alignment on the top_n only
        print("\nPerforming ClustalW-style alignment...")
        alignment_result = self.clustalw_alignment(scored_targets, ref_sequences)
        
        # Display alignment
        self.display_alignment(alignment_result)

    def clustalw_alignment(self, target_sequences, ref_sequences):
        """Perform ClustalW-style multiple sequence alignment"""
        # Combine all sequences
        all_sequences = []
        
        # Add reference sequences
        for i, seq in enumerate(ref_sequences):
            all_sequences.append({
                'sequence': seq,
                'name': f"Reference_{i+1}",
                'type': 'reference'
            })
        
        # Add target sequences (carry orig_location for display fidelity)
        for target in target_sequences:
            all_sequences.append({
                'sequence': target['sequence'],
                'name': target['name'],
                'coordinates': target['coordinates'],
                'strand': target['strand'],
                'mismatches': target.get('mismatches', 0),
                'orig_location': target.get('orig_location'),
                'type': 'target'
            })
        
        # Simple progressive alignment (ClustalW-style)
        aligned_sequences = self.progressive_alignment(all_sequences)
        
        return aligned_sequences

    def progressive_alignment(self, sequences):
        """Simple progressive alignment algorithm"""
        if len(sequences) < 2:
            return sequences
        
        # Start with first two sequences
        aligned = [sequences[0], sequences[1]]
        aligned[0]['aligned_seq'] = sequences[0]['sequence']
        aligned[1]['aligned_seq'] = sequences[1]['sequence']
        
        # Pairwise align first two sequences
        seq1, seq2 = self.pairwise_align(sequences[0]['sequence'], sequences[1]['sequence'])
        aligned[0]['aligned_seq'] = seq1
        aligned[1]['aligned_seq'] = seq2
        
        # Add remaining sequences one by one
        for i in range(2, len(sequences)):
            new_seq = sequences[i]
            # Align new sequence to consensus of existing alignment
            consensus = self.get_consensus(aligned)
            aligned_new, aligned_consensus = self.pairwise_align(new_seq['sequence'], consensus)
            
            # Update existing alignments based on consensus alignment
            gap_positions = []
            for j, (old_char, new_char) in enumerate(zip(consensus, aligned_consensus)):
                if old_char == '-' and new_char != '-':
                    gap_positions.append(j)
            
            # Insert gaps in existing sequences
            for existing in aligned:
                seq_with_gaps = existing['aligned_seq']
                for gap_pos in reversed(gap_positions):  # Insert from right to left
                    seq_with_gaps = seq_with_gaps[:gap_pos] + '-' + seq_with_gaps[gap_pos:]
                existing['aligned_seq'] = seq_with_gaps
            
            # Add new sequence
            new_seq['aligned_seq'] = aligned_new
            aligned.append(new_seq)
        
        return aligned

    def pairwise_align(self, seq1, seq2):
        """Pairwise alignment using dynamic programming with ClustalW DNA scoring"""
        # ClustalW DNA scoring matrix
        match_score = 2
        mismatch_score = -1
        gap_penalty = -2
        
        len1, len2 = len(seq1), len(seq2)
        
        # Initialize DP matrix
        dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]
        
        # Initialize first row and column
        for i in range(len1 + 1):
            dp[i][0] = i * gap_penalty
        for j in range(len2 + 1):
            dp[0][j] = j * gap_penalty
        
        # Fill DP matrix
        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                if seq1[i-1] == seq2[j-1]:
                    score = match_score
                else:
                    score = mismatch_score
                
                dp[i][j] = max(
                    dp[i-1][j-1] + score,  # Match/mismatch
                    dp[i-1][j] + gap_penalty,  # Gap in seq2
                    dp[i][j-1] + gap_penalty   # Gap in seq1
                )
        
        # Traceback
        aligned_seq1, aligned_seq2 = "", ""
        i, j = len1, len2
        
        while i > 0 or j > 0:
            if i > 0 and j > 0:
                if seq1[i-1] == seq2[j-1]:
                    score = match_score
                else:
                    score = mismatch_score
                
                if dp[i][j] == dp[i-1][j-1] + score:
                    aligned_seq1 = seq1[i-1] + aligned_seq1
                    aligned_seq2 = seq2[j-1] + aligned_seq2
                    i -= 1
                    j -= 1
                elif dp[i][j] == dp[i-1][j] + gap_penalty:
                    aligned_seq1 = seq1[i-1] + aligned_seq1
                    aligned_seq2 = '-' + aligned_seq2
                    i -= 1
                else:
                    aligned_seq1 = '-' + aligned_seq1
                    aligned_seq2 = seq2[j-1] + aligned_seq2
                    j -= 1
            elif i > 0:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = '-' + aligned_seq2
                i -= 1
            else:
                aligned_seq1 = '-' + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                j -= 1
        
        return aligned_seq1, aligned_seq2

    def get_consensus(self, aligned_sequences):
        """Get consensus sequence from aligned sequences"""
        if not aligned_sequences:
            return ""
        
        max_length = max(len(seq['aligned_seq']) for seq in aligned_sequences)
        consensus = ""
        
        for pos in range(max_length):
            base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0}
            
            for seq in aligned_sequences:
                if pos < len(seq['aligned_seq']):
                    base = seq['aligned_seq'][pos]
                    if base in base_counts:
                        base_counts[base] += 1
            
            # Choose most frequent base (excluding gaps for consensus)
            non_gap_counts = {k: v for k, v in base_counts.items() if k != '-'}
            if non_gap_counts:
                consensus_base = max(non_gap_counts, key=non_gap_counts.get)
            else:
                consensus_base = '-'
            
            consensus += consensus_base
        
        return consensus

    def display_alignment(self, aligned_sequences, top_n: int = 100):
        """
        Display alignment as pairwise comparisons of each target against each reference,
        marking identical positions with '*'. This matches the requested format.
        """
        if not aligned_sequences:
            print("No sequences to align.")
            return
        
        # Split into references and targets
        references = [s for s in aligned_sequences if s.get('type') == 'reference']
        targets = [s for s in aligned_sequences if s.get('type') == 'target']
        if not references or not targets:
            print("Need at least one reference and one target sequence for display.")
            return
        
        name_width = max(len(s['name']) for s in references + targets)
        info_width = 28
        
        def compute_star_line(a: str, b: str) -> str:
            line = []
            for ca, cb in zip(a, b):
                if ca == cb and ca != '-':
                    line.append('*')
                else:
                    line.append(' ')
            return ''.join(line)
        
        total_targets = len(targets)
        print(f"\nPairwise alignment relative to references (top {top_n} shown; best scores last):")
        print("=" * 80)
        
        for ref in references:
            ref_name = ref['name']
            # Use raw reference (unchanged, no flanks) for integrity checks
            ref_raw = ref['sequence']
            # Compute best window score per target relative to this reference
            scored = []
            win_len = len(ref_raw)
            for tgt in targets:
                tgt_raw = tgt['sequence']
                best_idx = 0
                best_mm = None
                for i in range(0, max(1, len(tgt_raw) - win_len + 1)):
                    mm = self._count_mismatches(ref_raw, tgt_raw[i:i+win_len])
                    if best_mm is None or mm < best_mm:
                        best_mm = mm
                        best_idx = i
                scored.append({
                    'name': tgt['name'],
                    'coords': tgt.get('coordinates', 'N/A'),
                    'strand': tgt.get('strand','N/A'),
                    'best_mm': 0 if best_mm is None else best_mm,
                    'window': tgt_raw[best_idx:best_idx+win_len]
                })
            # Sort by score asc, then name; keep top_n lowest scores but print best last
            scored.sort(key=lambda x: (x['best_mm'], x['name']))
            shown = scored[:top_n]
            for item in reversed(shown):
                tgt_name = item['name']
                tgt_window = item['window']
                # Extract just the match number from the name (e.g., "Match_1853_ydhJ-------ydhK" -> "Match_1853")
                match_num = tgt_name.split('_')[1] if '_' in tgt_name else tgt_name
                simple_name = f"Match_{match_num}"
                # Get TSS distance if available - ONLY for header line
                tss_distance_str = ""
                if self.tss_available and self.tss_analyzer:
                    try:
                        # Extract coordinates from item['coords'] (e.g., "3443852-3443878")
                        coord_str = item['coords']
                        if '-' in coord_str:
                            start_pos = int(coord_str.split('-')[0])
                            end_pos = int(coord_str.split('-')[1])
                            strand = item['strand']
                            tss_distance = self.tss_analyzer.calculate_tss_distance_string(start_pos, end_pos, strand)
                            tss_distance_str = f" | TSS: {tss_distance}"
                    except (ValueError, AttributeError, KeyError):
                        tss_distance_str = ""
                
                orig_location = next((t.get('orig_location') for t in targets if t['name']==item['name']), 'N/A')
                tgt_info_header = f"({item['coords']}) [{item['strand']}] score={item['best_mm']} | Location: {orig_location}{tss_distance_str}"
                tgt_info_line = f"({item['coords']}) [{item['strand']}] score={item['best_mm']} | Location: {orig_location}"
                
                print(f"\nReference: {ref_name}    Target: {simple_name} {tgt_info_header}")
                stars = ''.join('*' if a == b else ' ' for a, b in zip(ref_raw, tgt_window))
                print(f"{ref_name.ljust(name_width)}  {ref_raw}")
                print(f"{simple_name.ljust(name_width)}  {tgt_window}  ({item['coords']}) [{item['strand']}] score={item['best_mm']}")
                print(f"{'':<{name_width}}  {stars}")
        
        print("\n" + "=" * 80)
        print("Legend: * = identical positions relative to the reference")

def get_int_input(prompt, default=None):
    while True:
        value = input(prompt)
        if value == "" and default is not None:
            return default
        try:
            return int(value)
        except ValueError:
            print("Please enter a valid integer.")

def normalize_strand(strand):
    if isinstance(strand, str):
        if strand.lower().startswith('p'):
            return '+'
        elif strand.lower().startswith('n'):
            return '-'
    return strand

def main():
    analyzer = BSubAnalyzer()
    
    # Load sequence
    sequence_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "B.subtilis_Seq.txt")
    if not analyzer.load_sequence(sequence_file):
        print("Failed to load sequence")
        exit(1)
    
    # Load gene data
    gene_data_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "all_gene_data.json")
    if not analyzer.gene_data.load_from_json(gene_data_file):
        print("Failed to load gene data from all_gene_data.json")
        exit(1)
    
    while True:
        menu = (
            "\nChoose search mode:\n"
            "1. Search for sequence pattern\n"
            "2. Search for palindrome with mismatches\n"
            "3. Search for palindromes in region\n"
            "4. Sequence alignment with flanking regions\n"
            "5. Quit\n"
            "Enter mode (1-5): "
        )
        mode = input(menu).strip()

        if mode == "1":
            # Prompt for sequence with validation
            while True:
                pattern = input("Enter sequence to search (can use IUPAC codes): ").strip().upper()
                if not pattern:
                    print("Please enter a non-empty sequence.")
                    continue
                
                # Check for invalid characters
                valid_chars = set('ATCGRYSWKMBDHVN')
                invalid_chars = set(pattern) - valid_chars
                
                if invalid_chars:
                    print(f"Invalid characters found: {', '.join(sorted(invalid_chars))}")
                    print("Valid characters are: A, T, C, G, R, Y, S, W, K, M, B, D, H, V, N")
                    print("IUPAC codes: R=AG, Y=CT, S=GC, W=AT, K=GT, M=AC, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT")
                    continue
                
                break
            max_mismatches = get_int_input("Enter maximum allowed mismatches (default 0): ", default=0)
            results = analyzer.search_sequence(pattern, max_mismatches)
            # Create display_results for mode 1
            display_results = []
            for start, end, strand, pat, found_seq, location, mismatches in results:
                genome_forward = found_seq
                genome_revcomp = analyzer._reverse_complement(found_seq)
                display_results.append((start, end, strand, pat, genome_forward, genome_revcomp, location, mismatches))
        elif mode == "2":
            # Prompt for palindrome with validation
            while True:
                pattern = input("Enter palindrome to search (can use IUPAC codes): ").strip().upper()
                if not pattern:
                    print("Please enter a non-empty sequence.")
                    continue
                
                # Check for invalid characters
                valid_chars = set('ATCGRYSWKMBDHVN')
                invalid_chars = set(pattern) - valid_chars
                
                if invalid_chars:
                    print(f"Invalid characters found: {', '.join(sorted(invalid_chars))}")
                    print("Valid characters are: A, T, C, G, R, Y, S, W, K, M, B, D, H, V, N")
                    print("IUPAC codes: R=AG, Y=CT, S=GC, W=AT, K=GT, M=AC, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT")
                    continue
                
                break
            max_mismatches = get_int_input("Enter maximum allowed mismatches (default 0): ", default=0)
            results = analyzer.search_palindrome(pattern, max_mismatches)
        elif mode == "3":
            start_pos = get_int_input("Enter start position: ")
            end_pos = get_int_input("Enter end position: ")
            min_length = get_int_input("Enter minimum palindrome length (default 6): ", default=6)
            max_length = get_int_input("Enter maximum palindrome length (default 30): ", default=30)
            max_mismatches = get_int_input("Enter maximum allowed mismatches (default 0): ", default=0)
            results = analyzer.search_region_palindromes(start_pos, end_pos, min_length, max_length, max_mismatches)
        elif mode == "4":
            print("\nSequence Alignment with Flanking Regions")
            analyzer.perform_alignment_analysis()
            continue
        elif mode == "5":
            print("Goodbye!")
            break
        else:
            print("Invalid mode selected. Please try again.")
            continue

        # Notify user of total matches found
        total_matches = len(results)
        print(f"\nTotal matches found: {total_matches}")

        # Ask if user wants to filter results
        filter_option = input("\nWould you like to filter results? (y/n): ").lower().strip()
        if filter_option == 'y':
            print("\nFilter options:")
            print("1. Show only matches inside genes")
            print("2. Show only matches between genes (intergenic)")
            print("3. Show only matches with specific number of mismatches")
            print("4. Show only matches in specific region")
            
            filter_choice = input("Enter filter choice (1-4): ").strip()
            
            if filter_choice == "1":
                # Filter for matches inside genes
                filtered_results = []
                for result in results:
                    # Check if location contains "(in)" pattern (inside genes)
                    if "(in)" in result[5]:
                        filtered_results.append(result)
                
                if filtered_results:
                    results = filtered_results
                    print(f"Filtered to {len(results)} matches inside genes.")
                else:
                    print("No matches found inside genes. Showing all results.")
                    
            elif filter_choice == "2":
                # Filter for matches between genes
                filtered_results = []
                for result in results:
                    # Check if location contains "---" pattern (between genes)
                    # AND does NOT contain "(in)" pattern (not inside genes)
                    if "---" in result[5] and "(in)" not in result[5]:
                        filtered_results.append(result)
                
                if filtered_results:
                    results = filtered_results
                    print(f"Filtered to {len(results)} intergenic matches.")
                else:
                    print("No intergenic matches found. Showing all results.")
                    
            elif filter_choice == "3":
                # Filter by mismatches
                try:
                    max_mismatches = int(input("Enter maximum number of mismatches: "))
                    filtered_results = [r for r in results if r[6] <= max_mismatches]
                    if filtered_results:
                        results = filtered_results
                        print(f"Filtered to {len(results)} matches with ≤{max_mismatches} mismatches.")
                    else:
                        print(f"No matches found with ≤{max_mismatches} mismatches. Showing all results.")
                except ValueError:
                    print("Invalid input. Showing all results.")
                    
            elif filter_choice == "4":
                # Filter by region
                try:
                    region_input = input("Enter region (start,end): ")
                    start, end = map(int, region_input.split(','))
                    filtered_results = []
                    for result in results:
                        # Check if match overlaps with the specified region
                        if (start <= result[0] <= end or start <= result[1] <= end or 
                            result[0] <= start <= result[1] or result[0] <= end <= result[1]):
                            filtered_results.append(result)
                    
                    if filtered_results:
                        results = filtered_results
                        print(f"Filtered to {len(results)} matches in region {start}-{end}.")
                    else:
                        print(f"No matches found in region {start}-{end}. Showing all results.")
                except ValueError:
                    print("Invalid region format. Use 'start,end'. Showing all results.")
            else:
                print("Invalid choice. Showing all results.")

        # Ask user how many results they want to see
        current_total = len(results)
        
        # Update display_results for mode 1 after filtering
        if mode == "1":
            display_results = []
            for start, end, strand, pat, found_seq, location, mismatches in results:
                genome_forward = found_seq
                genome_revcomp = analyzer._reverse_complement(found_seq)
                display_results.append((start, end, strand, pat, genome_forward, genome_revcomp, location, mismatches))
        
        while True:
            display_option = input("How many results would you like to see? (Enter a number or 'all'): ").strip().lower()
            if display_option == 'all':
                display_count = current_total
                break
            try:
                display_count = int(display_option)
                if display_count <= 0:
                    print("Please enter a positive number or 'all'.")
                    continue
                if display_count > current_total:
                    print(f"Only {current_total} matches found. Showing all results.")
                    display_count = current_total
                break
            except ValueError:
                print("Invalid input. Please enter a number or 'all'.")

        # Print results in table format with TSS distance
        print("\nSearch Results:")
        if analyzer.tss_available:
            print("-" * 180)
            if mode == "1":
                print(f"{'Start':<10}{'End':<10}{'Strand':<8}{'Pattern':<18}{'Genome Forward':<20}{'Genome RevComp':<20}{'Distance from TSS':<25}{'Mismatches':<12}{'Location':<30}")
                print("-" * 180)
                for start, end, strand, pattern, genome_forward, genome_revcomp, location, mismatches in display_results[:display_count]:
                    tss_distance = analyzer.tss_analyzer.calculate_tss_distance_string(start, end, strand)
                    print(f"{str(start):<10}{str(end):<10}{strand:<8}{pattern:<18}{genome_forward[:18]:<20}{genome_revcomp[:18]:<20}{tss_distance[:23]:<25}{str(mismatches):<12}{location[:28]:<30}")
            else:
                print(f"{'Start':<10}{'End':<10}{'Strand':<8}{'Pattern':<14}{'Found Sequence':<18}{'RevComp':<18}{'Distance from TSS':<25}{'Mismatches':<12}{'Location':<30}")
                print("-" * 180)
                for start, end, strand, pattern, found_seq, location, mismatches in results[:display_count]:
                    revcomp = analyzer._reverse_complement(found_seq)
                    tss_distance = analyzer.tss_analyzer.calculate_tss_distance_string(start, end, strand)
                    print(f"{str(start):<10}{str(end):<10}{strand:<8}{pattern:<14}{found_seq[:16]:<18}{revcomp[:16]:<18}{tss_distance[:23]:<25}{str(mismatches):<12}{location[:28]:<30}")
        else:
            # Fallback to original display without TSS distance
            print("-" * 140)
            if mode == "1":
                print(f"{'Start':<10}{'End':<10}{'Strand':<8}{'Pattern':<18}{'Genome Forward':<20}{'Genome RevComp':<20}{'Mismatches':<12}{'Location':<30}")
                print("-" * 140)
                for start, end, strand, pattern, genome_forward, genome_revcomp, location, mismatches in display_results[:display_count]:
                    print(f"{str(start):<10}{str(end):<10}{strand:<8}{pattern:<18}{genome_forward[:18]:<20}{genome_revcomp[:18]:<20}{str(mismatches):<12}{location[:28]:<30}")
            else:
                print(f"{'Start':<10}{'End':<10}{'Strand':<8}{'Pattern':<14}{'Found Sequence':<18}{'RevComp':<18}{'Mismatches':<12}{'Location':<35}")
                print("-" * 145)
                for start, end, strand, pattern, found_seq, location, mismatches in results[:display_count]:
                    revcomp = analyzer._reverse_complement(found_seq)
                    print(f"{str(start):<10}{str(end):<10}{strand:<8}{pattern:<14}{found_seq[:16]:<18}{revcomp[:16]:<18}{str(mismatches):<12}{location[:33]:<35}")
        
        line_width = 180 if analyzer.tss_available else 140
        print("-" * line_width)
        print(f"Total matches shown: {len(results[:display_count])}")
        
        # Ask if user wants to save results
        save = input("\nDo you want to save the results to a file? (y/n): ").lower()
        while save == 'y':
            filename = input("Enter filename to save results (default: search_results.csv): ") or "search_results.csv"
            try:
                analyzer.save_results(results, filename)
                print(f"Results saved to {filename}")
                break
            except PermissionError as e:
                print(f"Error saving file: {e}")
                print("Please close the file if it is open, or choose a different filename.")
                retry = input("Do you want to try a different filename? (y/n): ").lower()
                if retry != 'y':
                    print("Skipping save.")
                    break
            except Exception as e:
                print(f"Unexpected error: {e}")
                print("Skipping save.")
                break

if __name__ == "__main__":
    main()
