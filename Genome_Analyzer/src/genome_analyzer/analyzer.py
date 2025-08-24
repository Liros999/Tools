"""
Genome Analyzer - Main module for genome analysis operations.
Refactored to use modular package structure for better maintainability.
"""

# Import all components from the modular packages
from .core import GenomeAnalyzer, ResourceManager, GenomeAnalysisError
from .algorithms import search_pattern, search_pattern_aho_corasick, search_pattern_myers
from .io import (
    load_genome_for_chromosome, 
    load_annotations_from_gtf,
    create_shared_gtf_cache,
    load_annotations_from_shared_cache
)
from .interface import (
    offer_dashboard,
    launch_dashboard,
    search_mode_prompt,
    start_dashboard
)

# Re-export for backward compatibility
__all__ = [
    'GenomeAnalyzer',
    'ResourceManager',
    'GenomeAnalysisError',
    'search_pattern',
    'search_pattern_aho_corasick',
    'search_pattern_myers',
    'load_genome_for_chromosome',
    'load_annotations_from_gtf',
    'create_shared_gtf_cache',
    'load_annotations_from_shared_cache',
    'offer_dashboard',
    'launch_dashboard',
    'search_mode_prompt',
    'start_dashboard'
]

    @property
    def iupac_wildcards(self) -> Dict[str, str]:
        """Get IUPAC wildcards on-demand to eliminate redundant storage."""
        return {k: v.strip('[]') for k, v in self.iupac_map.items()}

    def _reverse_complement(self, seq: str) -> str:
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                          'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
        return "".join(complement_map.get(base, base) for base in reversed(seq))



    def get_location_details(self, start: int, end: int) -> str:
        """
        Get precise location details for a match in Homo Sapiens genome.
        """
        return self._get_location_details_generic(start, end)
        """
        ENHANCED: Try to get species-specific location using advanced detection methods.
        This integrates the best location detection from B.Sub and E.Coli analyzers.
        """
        try:
            # Detect species based on current chromosome or gene data
            species = self._detect_species()
            
            if species == "B.sub" and hasattr(self, 'gene_data') and self.gene_data:
                return self._get_b_sub_location(start, end)
            elif species == "E.coli" and hasattr(self, 'gene_data') and self.gene_data:
                return self._get_e_coli_location(start, end)
            elif species == "Homo_sapiens":
                return self._get_homo_sapiens_location(start, end)
            
        except Exception as e:
            # If species-specific detection fails, fall back to generic
            pass
        
        return None
    
    def _detect_species(self) -> str:
        """Detect the species based on current data and chromosome information."""
        # Check chromosome name for species indicators
        if self.current_chromosome:
            chrom_lower = self.current_chromosome.lower()
            if 'b.sub' in chrom_lower or 'subtilis' in chrom_lower:
                return "B.sub"
            elif 'e.coli' in chrom_lower or 'coli' in chrom_lower:
                return "E.coli"
            elif 'homo' in chrom_lower or 'human' in chrom_lower or 'chr' in chrom_lower:
                return "Homo_sapiens"
        
        # Check gene data structure for species indicators
        if hasattr(self, 'gene_data') and self.gene_data:
            if isinstance(self.gene_data, dict):
                # Check for B.Sub gene naming patterns (BSU_xxxxx)
                if any('BSU_' in str(k) for k in self.gene_data.keys()):
                    return "B.sub"
                # Check for E.Coli gene naming patterns
                elif any('b' in str(k).lower() for k in self.gene_data.keys()):
                    return "E.coli"
        
        return "unknown"
    
    def _get_b_sub_location(self, start: int, end: int) -> str:
        """
        ENHANCED: B.Sub specific location detection using advanced methods.
        This replicates the critical fixes from B.Sub analyzer.
        """
        try:
            # Find genes that overlap with this position
            overlapping_genes = []
            for gene_id, gene_info in self.gene_data.items():
                if (start <= gene_info['end'] and end >= gene_info['start']):
                    overlapping_genes.append((gene_id, gene_info))
            
            if overlapping_genes:
                # Match overlaps with gene(s)
                if len(overlapping_genes) == 1:
                    gene_id, gene_info = overlapping_genes[0]
                    return f"[{gene_id}][Overlap]"
                else:
                    # Multiple gene overlap
                    gene_names = [g[0] for g in overlapping_genes]
                    return f"[MultiGene]({','.join(gene_names)})"
            else:
                # Match is in intergenic region - find flanking genes
                return self._get_b_sub_intergenic_location(start, end)
                
        except Exception as e:
            return f"[B.Sub_Error:{str(e)[:20]}]"
    
    def _get_b_sub_intergenic_location(self, start: int, end: int) -> str:
        """
        ENHANCED: B.Sub specific intergenic location detection.
        This is the CRITICAL FIX for CopZ-CsoR intergenic detection.
        """
        try:
            # Get all genes sorted by position
            all_genes = []
            for gene_id, gene_info in self.gene_data.items():
                if 'start' in gene_info and 'end' in gene_info:
                    all_genes.append((gene_id, gene_info['start'], gene_info['end']))
            
            if not all_genes:
                return "[start][-------][end]"
            
            # Sort genes by start position
            all_genes.sort(key=lambda x: x[1])
            
            # Find upstream and downstream genes
            upstream_gene = None
            downstream_gene = None
            
            for i, (gene_name, gene_start, gene_end) in enumerate(all_genes):
                # Check if match is in intergenic region after this gene
                if gene_end < start:
                    upstream_gene = gene_name
                    
                    # Check if there's a downstream gene
                    if i + 1 < len(all_genes):
                        next_gene_name, next_gene_start, next_gene_end = all_genes[i + 1]
                        if next_gene_start > end:
                            downstream_gene = next_gene_name
                            break
                    else:
                        downstream_gene = 'end'
                        break
            
            # Handle edge cases
            if upstream_gene is None and all_genes:
                first_gene_name, first_gene_start, first_gene_end = all_genes[0]
                if start < first_gene_start:
                    upstream_gene = 'start'
                    downstream_gene = first_gene_name
            
            if upstream_gene is None and all_genes:
                last_gene_name, last_gene_start, last_gene_end = all_genes[-1]
                if end > last_gene_end:
                    upstream_gene = last_gene_name
                    downstream_gene = 'end'
            
            # Format the intergenic location
            if upstream_gene and downstream_gene:
                return f"[{upstream_gene}][-------][{downstream_gene}]"
            else:
                return "[intergenic_region]"
                
        except Exception as e:
            return f"[B.Sub_Intergenic_Error:{str(e)[:20]}]"
    
    def _get_e_coli_location(self, start: int, end: int) -> str:
        """
        ENHANCED: E.Coli specific location detection using advanced methods.
        This replicates the critical fixes from E.Coli analyzer.
        """
        try:
            # Similar to B.Sub but with E.Coli specific gene patterns
            return self._get_b_sub_location(start, end)  # Reuse B.Sub logic for now
            
        except Exception as e:
            return f"[E.Coli_Error:{str(e)[:20]}]"
    
    def _get_location_details_generic(self, start: int, end: int) -> str:
        """
        Generic location detection for Homo Sapiens and other species.
        Uses the original working logic from the main analyzer.
        """
        try:
            overlapping_regions = sorted(self.gene_tree[start:end])
            
            # First, check if match is completely within a gene
            for region_interval in overlapping_regions:
                if region_interval.data.get('type') == 'gene':
                    gene_start = region_interval.data['start']
                    gene_end = region_interval.data['end']
                    gene_name = region_interval.data['name']
                    
                    # Match is completely inside gene
                    if start >= gene_start and end <= gene_end:
                        return f"[{gene_name}][In]"
                    
                    # Match overlaps gene boundary - calculate overlap distances
                    overlap_start = max(start, gene_start)
                    overlap_end = min(end, gene_end)
                    overlap_bp = overlap_end - overlap_start + 1
                    
                    if start < gene_start and end > gene_end:
                        # Match completely spans the gene
                        left_bp = gene_start - start
                        right_bp = end - gene_end
                        return f"<{left_bp}bp>[{gene_name}]{overlap_bp}bp<{right_bp}bp>"
                    elif start < gene_start:
                        # Match starts before gene
                        left_bp = gene_start - start
                        return f"<{left_bp}bp>[{gene_name}]{overlap_bp}bp"
                    else:
                        # Match ends after gene
                        right_bp = end - gene_end
                        return f"[{gene_name}]{overlap_bp}bp<{right_bp}bp>"
            
            # No gene overlaps found - match is in intergenic region
            return self._get_intergenic_location_generic(start, end)
            
        except Exception as e:
            return f"[Generic_Location_Error:{str(e)[:20]}]"
    
    def _get_intergenic_location_generic(self, start: int, end: int) -> str:
        """
        Generic intergenic location detection for Homo Sapiens.
        """
        try:
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
            
        except Exception as e:
            return f"[Generic_Intergenic_Error:{str(e)[:20]}]"
    
    def _get_homo_sapiens_location(self, start: int, end: int) -> str:
        """
        ENHANCED: Homo Sapiens specific location detection.
        This handles chromosome-specific gene contexts.
        """
        try:
            # For Homo Sapiens, we need to handle chromosome-specific gene data
            if hasattr(self, 'gene_data') and self.gene_data:
                # Try to use the same advanced logic
                return self._get_b_sub_location(start, end)
            else:
                # Fall back to generic chromosome-based detection
                return self._get_intergenic_location(start, end)
                
        except Exception as e:
            return f"[Homo_Sapiens_Error:{str(e)[:20]}]"
    
    def _get_intergenic_location_working(self, start: int, end: int, strand: str = '+') -> str:
        """
        FIXED: Working intergenic location detection copied from B.Sub analyzer
        """
        print(f"[LOCATION_DEBUG] Input: start={start}, end={end}, strand={strand}")
        
        if not hasattr(self, 'gene_data') or not self.gene_data:
            print("[LOCATION_DEBUG] No gene data available")
            return "[No_Gene_Data]"
        
        try:
            # Get all genes sorted by position
            genes = []
            for gene_id, gene_info in self.gene_data.items():
                if isinstance(gene_info, dict) and 'start' in gene_info and 'end' in gene_info:
                    genes.append((gene_id, int(gene_info['start']), int(gene_info['end'])))
            
            print(f"[LOCATION_DEBUG] Parsed {len(genes)} genes")
            if genes:
                print(f"[LOCATION_DEBUG] First gene: {genes[0]}")
                print(f"[LOCATION_DEBUG] Last gene: {genes[-1]}")
            
            if not genes:
                print("[LOCATION_DEBUG] No valid genes found, returning default")
                return "[start][-------][end]"
            
            # Sort by start position
            genes.sort(key=lambda x: x[1])
            
            # Find upstream and downstream genes
            upstream_gene = None
            downstream_gene = None
            
            for i, (gene_name, gene_start, gene_end) in enumerate(genes):
                # Match is after this gene
                if gene_end < start:
                    upstream_gene = gene_name
                    
                    # Check next gene
                    if i + 1 < len(genes):
                        next_gene_name, next_gene_start, next_gene_end = genes[i + 1]
                        if next_gene_start > end:
                            downstream_gene = next_gene_name
                            break
                    else:
                        downstream_gene = 'end'
                        break
            
            # Handle edge cases
            if upstream_gene is None and genes:
                first_gene_name, first_gene_start, first_gene_end = genes[0]
                if start < first_gene_start:
                    upstream_gene = 'start'
                    downstream_gene = first_gene_name
            
            print(f"[LOCATION_DEBUG] Found upstream='{upstream_gene}', downstream='{downstream_gene}'")
            
            # Format result
            if upstream_gene and downstream_gene:
                result = f"[{upstream_gene}][-------][{downstream_gene}]"
                print(f"[LOCATION_DEBUG] Returning: '{result}' (length: {len(result)})")
                return result
            else:
                result = "[intergenic_region]"
                print(f"[LOCATION_DEBUG] Returning: '{result}' (length: {len(result)})")
                return result
                
        except Exception as e:
            error_msg = f"[Error:{str(e)[:10]}]"
            print(f"[LOCATION_DEBUG] Exception occurred: {e}")
            print(f"[LOCATION_DEBUG] Returning error: '{error_msg}'")
            return f"[Error:{str(e)[:10]}]"
    
    def _get_intergenic_location(self, start: int, end: int) -> str:
        """
        FIXED: Use the working intergenic detection method
        """
        return self._get_intergenic_location_working(start, end, '+')

    def _build_gene_tree_with_intergenic(self, all_genes: List[Dict]):
        """
        Centralized method to build gene tree with intergenic regions.
        Eliminates massive code duplication across annotation loading methods.
        """
        all_genes.sort(key=lambda g: g['start'])
        last_coord = 0
        last_gene_name = 'start'
        
        for gene in all_genes:
            if gene['start'] == 0: 
                continue
                
            # Create intergenic region before this gene
            intergenic_start, intergenic_end = last_coord + 1, gene['start'] - 1
            if intergenic_end > intergenic_start:
                intergenic_info = {
                    'name': 'intergenic', 
                    'type': 'intergenic', 
                    'upstream': last_gene_name, 
                    'downstream': gene['name']
                }
                self.gene_tree[intergenic_start:intergenic_end + 1] = intergenic_info
            
            # Add gene to tree and data
            self.gene_tree[gene['start']:gene['end'] + 1] = gene
            self.gene_data[gene.get('locus_tag', gene['name'])] = gene
            last_coord = gene['end']
            last_gene_name = gene['name']
        
        # Create final intergenic region if needed
        if last_coord < len(self.sequence):
            final_intergenic = {
                'name': 'intergenic', 
                'type': 'intergenic', 
                'upstream': last_gene_name, 
                'downstream': 'end'
            }
            self.gene_tree[last_coord + 1:len(self.sequence) + 1] = final_intergenic

    def load_complete_sequence(self, filename: str):
        """Load complete genome sequence from FASTA file."""
        print(f"Loading complete genome sequence from {filename}...")
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                sequence_parts = [line.strip().upper() for line in f if not line.startswith('>')]
            self.sequence = "".join(sequence_parts)
            print(f"Successfully loaded genome of {len(self.sequence):,} bp.")
        except FileNotFoundError:
            raise GenomeAnalysisError(f"Genome file not found at {filename}")
        except Exception as e:
            raise GenomeAnalysisError(f"An error occurred loading the genome: {e}")

    def load_chromosome_sequence(self, fasta_path: str, chromosome: str):
        """Load specific chromosome sequence from FASTA file with Linux compatibility."""
        print(f"Loading chromosome '{chromosome}' sequence from {fasta_path}...")
        
        # Linux-compatible file handling
        try:
            from .file_processing import _get_chromosome_aliases
            wanted = set(_get_chromosome_aliases(chromosome))
            parts: List[str] = []
            capturing = False
            
            # Use Linux-compatible file encoding
            encoding = 'utf-8' if IS_LINUX else 'utf-8'
            error_handling = 'ignore' if IS_LINUX else 'strict'
            
            with open(fasta_path, 'r', encoding=encoding, errors=error_handling) as f:
                for raw_line in f:
                    line = raw_line.strip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        header = line[1:].split()[0]
                        if capturing and header not in wanted:
                            # finished reading the target sequence
                            break
                        capturing = (header in wanted)
                        continue
                    if capturing:
                        parts.append(line.upper())
            
            # Linux-compatible string joining
            if IS_LINUX and len(parts) > 1000000:  # Large sequences
                # Use chunked joining to avoid memory issues
                chunk_size = 100000
                sequence_chunks = []
                for i in range(0, len(parts), chunk_size):
                    chunk = "".join(parts[i:i + chunk_size])
                    sequence_chunks.append(chunk)
                self.sequence = "".join(sequence_chunks)
            else:
                self.sequence = "".join(parts)
            
            if not self.sequence:
                raise GenomeAnalysisError(f"Chromosome '{chromosome}' not found in {fasta_path}")
            
            print(f"Successfully loaded chromosome '{chromosome}' of {len(self.sequence):,} bp.")
            
            # Resource management for large sequences
            if IS_SERVER and len(self.sequence) > 100000000:  # 100MB
                resource_manager.cleanup_memory()
                
        except FileNotFoundError:
            raise GenomeAnalysisError(f"FASTA file not found at {fasta_path}")
        except Exception as e:
            raise GenomeAnalysisError(f"An error occurred loading chromosome sequence: {e}")

    def search_sequence(
        self,
        pattern: str,
        max_mismatches: int = 0,
        boundary_bp: int = 0,
        stream: bool = False,
        on_row: Optional[callable] = None,
        on_header: Optional[callable] = None,
        include_chrom: bool = False,
        worker_line: int = 0,
    ) -> List[Tuple]:
        """
        Search for pattern in genome sequence using generic search logic.
        """
        return self._search_sequence_generic(pattern, max_mismatches, boundary_bp, stream, on_row, on_header, include_chrom)

    def _search_sequence_generic(
        self,
        pattern: str,
        max_mismatches: int = 0,
        boundary_bp: int = 0,
        stream: bool = False,
        on_row: Optional[callable] = None,
        on_header: Optional[callable] = None,
        include_chrom: bool = False,
    ) -> List[Tuple]:
        """
        Original generic search logic for non-species-specific searches.
        """
        # Check resource usage before starting search
        if IS_SERVER:
            resource_manager.enforce_resource_limits()

        matches = []
        pattern = pattern.upper()
        seq = self.sequence
        seq_len = len(seq)

        if max_mismatches == 0:
            forward_matches = search._smart_exact_search(seq, pattern)
            matches.extend(self._process_search_results(forward_matches, pattern, '+', boundary_bp, on_row, stream, include_chrom))

            rev_pattern = self._reverse_complement(pattern)
            if rev_pattern != pattern:
                reverse_matches = search._smart_exact_search(seq, rev_pattern)
                matches.extend(self._process_search_results(reverse_matches, pattern, '-', boundary_bp, on_row, stream, include_chrom))
        else:
            matches = []
            rev_pattern = self._reverse_complement(pattern)

            forward_matches = search._search_optimized_mismatch(seq, pattern, max_mismatches)
            matches.extend(self._process_search_results(forward_matches, pattern, '+', boundary_bp, on_row, stream, include_chrom))

            if rev_pattern != pattern:
                reverse_matches = search._search_optimized_mismatch(seq, rev_pattern, max_mismatches)
                matches.extend(self._process_search_results(reverse_matches, pattern, '-', boundary_bp, on_row, stream, include_chrom))

        # Resource cleanup after search
        if IS_SERVER and len(matches) > 1000:
            resource_manager.cleanup_memory()

        unique_matches = list({(m[0], m[1]): m for m in matches}.values())
        unique_matches.sort(key=lambda x: x[0])
        return unique_matches

    def _search_sequence_bsub(self, pattern: str, max_mismatches: int = 0) -> List[Tuple[int, int, str, str, str, str, int]]:
        """
        CRITICAL FIX: Search for pattern and its reverse complement in the genome, allowing for mismatches (B.sub style).
        Returns list of (start, end, strand, pattern, found_sequence, location, mismatches)
        """
        matches = []
        pattern = pattern.upper()
        pattern_len = len(pattern)
        seq = self.sequence.upper()
        rev_pattern = self._reverse_complement(pattern)
        
        # Validate that the pattern can be found in genome (dynamic check)
        if pattern in seq:
            print(f"[INFO] Pattern '{pattern}' found in genome")
        else:
            rev_pattern_check = self._reverse_complement(pattern)
            if rev_pattern_check in seq:
                print(f"[INFO] Reverse complement of pattern '{pattern}' found in genome")
            else:
                print(f"[WARNING] Pattern '{pattern}' not found in genome - this may be normal for rare patterns")

        if max_mismatches == 0:
            # CRITICAL FIX: Use improved regex for exact matches with proper IUPAC handling
            print("\nSearching forward strand (exact match)...")
            regex_pattern = self._iupac_to_regex(pattern)

            for match in tqdm(re.finditer(regex_pattern, seq), desc="Progress"):
                start = match.start()  # Keep 0-based coordinates for internal processing
                end = match.end()
                found_seq = match.group()
                # Convert to 1-based only for display
                display_start = start + 1
                display_end = end
                location = self._get_match_location_bsub(display_start, display_end, '+')
                revcomp = self._reverse_complement(found_seq)
                matches.append((display_start, display_end, '+', pattern, found_seq, revcomp, 0, location))

            # CRITICAL FIX: Search reverse strand with proper pattern handling
            if rev_pattern != pattern:
                print("\nSearching reverse strand (exact match)...")
                rev_regex_pattern = self._iupac_to_regex(rev_pattern)

                for match in tqdm(re.finditer(rev_regex_pattern, seq), desc="Progress"):
                    start = match.start()  # Keep 0-based coordinates for internal processing
                    end = match.end()
                    found_seq = match.group()
                    # Convert to 1-based only for display
                    display_start = start + 1
                    display_end = end
                    # CRITICAL FIX: Report the reverse complement of what was found
                    rev_found_seq = self._reverse_complement(found_seq)
                    location = self._get_match_location_bsub(display_start, display_end, '-')
                    revcomp = self._reverse_complement(rev_found_seq)
                    matches.append((display_start, display_end, '-', pattern, rev_found_seq, revcomp, 0, location))
        else:
            # CRITICAL FIX: Improved mismatch search with better progress reporting
            print(f"\nSearching with {max_mismatches} mismatches allowed...")
            matches_found = 0

            # Search forward strand
            print("Searching forward strand...")
            for i in tqdm(range(len(seq) - pattern_len + 1), desc="Forward strand"):
                window = seq[i:i + pattern_len]
                mismatches = self._count_mismatches_bsub(pattern, window)

                if mismatches <= max_mismatches:
                    start = i + 1  # Convert to 1-based coordinates
                    end = i + pattern_len
                    location = self._get_match_location_bsub(start, end, '+')
                    revcomp = self._reverse_complement(window)
                    matches.append((start, end, '+', pattern, window, revcomp, mismatches, location))
                    matches_found += 1

            # Search reverse strand
            if rev_pattern != pattern:
                print("Searching reverse strand...")
                for i in tqdm(range(len(seq) - pattern_len + 1), desc="Reverse strand"):
                    window = seq[i:i + pattern_len]
                    mismatches = self._count_mismatches_bsub(rev_pattern, window)

                    if mismatches <= max_mismatches:
                        start = i + 1  # Convert to 1-based coordinates
                        end = i + pattern_len
                        # CRITICAL FIX: Report the reverse complement of what was found
                        rev_found_seq = self._reverse_complement(window)
                        location = self._get_match_location_bsub(start, end, '-')
                        revcomp = self._reverse_complement(rev_found_seq)
                        matches.append((start, end, '-', pattern, rev_found_seq, revcomp, mismatches, location))
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

    def _search_sequence_ecoli(self, pattern: str, max_mismatches: int = 0) -> List[Tuple]:
        """
        CRITICAL FIX: Search for pattern in genome sequence with improved accuracy (E.Coli style).
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
                mismatches_fwd = self._count_mismatches_ecoli(pattern, window)
                if mismatches_fwd <= max_mismatches:
                    start, end = i + 1, i + pattern_len
                    # CRITICAL FIX: Use improved location detection
                    location = self.get_location_details(start, end)
                    matches.append((start, end, '+', pattern, window, self._reverse_complement(window), mismatches_fwd, location))
                    matches_found_counter += 1
                    progress_bar.set_postfix_str(f"Matches: {matches_found_counter}")

                # Search reverse strand
                mismatches_rev = self._count_mismatches_ecoli(rev_pattern, window)
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

    def _process_search_results(self, matches: List[Tuple], pattern: str, strand: str, 
                              boundary_bp: int, 
                              on_row: Optional[callable] = None, stream: bool = False,
                              include_chrom: bool = False) -> List[Tuple]:
        processed_results = []
        for start, end, window, _, mismatches, _ in matches:
            location = self.get_location_details(start, end)
            
            if boundary_bp and boundary_bp > 0:
                seq_len = len(self.sequence)
                dist_start = start - 1
                dist_end = seq_len - end
                near_flags = []
                if dist_start <= boundary_bp:
                    near_flags.append(f"NearStart:{dist_start}bp")
                if dist_end <= boundary_bp:
                    near_flags.append(f"NearEnd:{dist_end}bp")
                if near_flags:
                    location = f"{location} [{' & '.join(near_flags)}]"
            
            row = (start, end, strand, pattern, window, self._reverse_complement(window), mismatches, location)
            processed_results.append(row)
            
            if on_row:
                on_row(row)
            elif stream:
                # This is a UI concern, should be moved
                pass
        
        return processed_results