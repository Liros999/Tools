"""Core scanning logic and scheduling."""
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Iterable
from datetime import datetime

from .pairing import revcomp, is_wc, is_gu
from .io_utils import load_text_fasta_or_plain, load_gene_names, load_json, print_error
from .progress import ProgressBar, progress_print_line

# Optional tqdm for nicer progress bars
try:
    from tqdm import tqdm  # type: ignore
except Exception:
    tqdm = None

def extend_match(srna: str, genome: str, srna_start: int, genome_anchor: int, seed_len: int, max_gu: int) -> Tuple[int, int, int, int]:
    gu_count = 0
    for k in range(seed_len):
        a = srna[srna_start + k]
        b = genome[genome_anchor - k]
        if is_wc(a, b):
            continue
        if is_gu(a, b):
            gu_count += 1
            if gu_count > max_gu:
                return (-1, -1, -1, -1)
        else:
            return (-1, -1, -1, -1)

    sr = srna_start + seed_len
    gr = genome_anchor - seed_len
    while sr < len(srna) and gr >= 0:
        a = srna[sr]
        b = genome[gr]
        if is_wc(a, b):
            sr += 1
            gr -= 1
            continue
        if is_gu(a, b):
            if gu_count + 1 > max_gu:
                break
            gu_count += 1
            sr += 1
            gr -= 1
            continue
        break

    sl = srna_start - 1
    gl = genome_anchor + 1
    while sl >= 0 and gl < len(genome):
        a = srna[sl]
        b = genome[gl]
        if is_wc(a, b):
            sl -= 1
            gl += 1
            continue
        if is_gu(a, b):
            if gu_count + 1 > max_gu:
                break
            gu_count += 1
            sl -= 1
            gl += 1
            continue
        break

    srna_left = sl + 1
    srna_right = sr - 1
    genome_left = gr + 1
    genome_right = gl - 1
    return (srna_left, srna_right, genome_left, genome_right)

def build_seed_index(srna: str, seed_len: int) -> Dict[str, List[int]]:
    index: Dict[str, List[int]] = {}
    if seed_len <= 0 or seed_len > len(srna):
        return index
    for i in range(0, len(srna) - seed_len + 1):
        seed = srna[i:i + seed_len]
        if any(ch not in 'ACGTU' for ch in seed):
            continue
        seed_rc = revcomp(seed)
        index.setdefault(seed_rc, []).append(i)
    return index

def scan_chunk(args: Tuple[str, int, int, str, int, int, bool, int]) -> List[Dict[str, object]]:
    genome, chunk_start, chunk_end, srna, seed_len, max_gu, is_rc_genome, genome_len = args
    
    results: List[Dict[str, object]] = []
    right = min(len(genome), chunk_end)
    left = max(0, chunk_start)
    if seed_len <= 0 or right - left < seed_len:
        return []

    # Scan all possible starting positions in the genome chunk
    for pos in range(left, right - seed_len + 1):
        # Get genome window at this position
        genome_window = genome[pos:pos + seed_len]
        if any(ch not in 'ACGTU' for ch in genome_window):
            continue
            
        # For each possible sRNA starting position
        for srna_start in range(0, len(srna) - seed_len + 1):
            # Check if this seed window is complementary (allowing GU wobbles)
            gu_count = 0
            valid_seed = True
            
            for k in range(seed_len):
                srna_base = srna[srna_start + k]
                genome_base = genome_window[k]
                
                # Check Watson-Crick pairing
                if is_wc(srna_base, genome_base):
                    continue
                # Check G-U wobble pairing
                elif is_gu(srna_base, genome_base):
                    gu_count += 1
                    if gu_count > max_gu:
                        valid_seed = False
                        break
                else:
                    valid_seed = False
                    break
            
            if not valid_seed:
                continue
                
            # Fixed window: use exactly seed_len positions
            srna_left = srna_start
            srna_right = srna_start + seed_len - 1
            genome_left = pos
            genome_right = pos + seed_len - 1
            
            # Calculate target coordinates and sequence
            if not is_rc_genome:
                target_start = genome_left + 1
                target_end = genome_right + 1
                target_seq = genome[genome_left:genome_right + 1]
            else:
                # For reverse complement genome, convert coordinates
                f_start = (genome_len - 1) - genome_right
                f_end = (genome_len - 1) - genome_left
                target_start = f_start + 1
                target_end = f_end + 1
                target_seq = revcomp(genome[genome_left:genome_right + 1])
            
            smallrna_seq = srna[srna_left:srna_right + 1]
            
            # Count GU pairs in the final alignment
            final_gu_count = 0
            for k in range(seed_len):
                srna_base = smallrna_seq[k]
                target_base = target_seq[k]
                if is_gu(srna_base, target_base):
                    final_gu_count += 1
            
            # Enforce wobble threshold
            if final_gu_count > max_gu:
                continue
            
            results.append({
                "SmallRNA_region": [srna_left + 1, srna_right + 1],
                "Number_of_GU": final_gu_count,
                "SmallRNA_seq": smallrna_seq,
                "Target_region": [target_start, target_end],
                "Target_seq": target_seq,
                "Strand": "-" if is_rc_genome else "+",
                "Genome_Forward": revcomp(target_seq) if is_rc_genome else target_seq,
                "Genome_RevComp": target_seq if is_rc_genome else revcomp(target_seq),
                "Window": seed_len,
            })
    
    return results

def format_location(target_start, target_end, gene_data):
    # Find all genes that overlap with the match
    overlapping_genes = []
    for gene_name, gene_info in gene_data.items():
        gene_start = gene_info['start']
        gene_end = gene_info['end']
        if max(target_start, gene_start) <= min(target_end, gene_end):
            overlapping_genes.append({'name': gene_name, 'start': gene_start, 'end': gene_end})

    if not overlapping_genes:
        # No overlapping genes, find genes between which the match lies
        genes_before = [g for g in gene_data.values() if g['end'] < target_start]
        genes_after = [g for g in gene_data.values() if g['start'] > target_end]

        if not genes_before or not genes_after:
            return "Intergenic"

        # Find the closest gene before and after
        closest_gene_before = max(genes_before, key=lambda g: g['end'])
        closest_gene_after = min(genes_after, key=lambda g: g['start'])
        
        # Need to get the name of the gene
        gene_name_before = [name for name, info in gene_data.items() if info == closest_gene_before][0]
        gene_name_after = [name for name, info in gene_data.items() if info == closest_gene_after][0]

        return f"[{gene_name_before}]------[{gene_name_after}]"

    # Handle single overlapping gene
    if len(overlapping_genes) == 1:
        gene = overlapping_genes[0]
        gene_name = gene['name']
        gene_start = gene['start']
        gene_end = gene['end']

        # Case 1: Match is completely inside the gene
        if gene_start <= target_start and target_end <= gene_end:
            return f"[{gene_name}](in)"

        # Case 2: Match contains the gene
        if target_start <= gene_start and gene_end <= target_end:
            left_overhang = gene_start - target_start
            right_overhang = target_end - gene_end
            return f"<{left_overhang}>[{gene_name}]<{right_overhang}>"

        # Case 3: Partial overlap
        if target_start < gene_start: # Overlap on the left of the gene
            return f"<{gene_start - target_start}>[{gene_name}]"
        else: # Overlap on the right of the gene
            return f"[{gene_name}]<{target_end - gene_end}>"

    # Handle multiple overlapping genes
    return ",".join([f"[{g['name']}]" for g in overlapping_genes])

def run_search(srna_seq: str, chunk_size: int, max_gu: int, min_len: int, genome_path: str, synonyms_path: str, gene_data_path: str, save_csv_path: str = None, results_dir: str = None, window_sizes: Iterable[int] = None, prompt_save: bool = True) -> None:
    if not srna_seq or any(ch not in 'ACGTUacgtu' for ch in srna_seq):
        print_error(
            err_type="InvalidInput",
            message="sRNA must be a non-empty DNA sequence consisting of A/C/G/T (U allowed)",
            citation="See Tools/Subtiwiki_API.txt",
            recommended_next_step="Provide the nucleotide sequence of the sRNA.",
        )

    bar = ProgressBar(2, prefix="Initializing", file=open(os.devnull, 'w', encoding='utf-8'))
    bar.update()

    genome_seq = load_text_fasta_or_plain(genome_path)
    if not genome_seq:
        print_error(
            err_type="DataError",
            message="Loaded genome sequence is empty after sanitization",
            citation="Genome file content validation",
            recommended_next_step="Validate the genome file content and format (FASTA/plain).",
        )
    bar.update()
    bar.finish()

    bar = ProgressBar(1, prefix="Loading synonyms", file=open(os.devnull, 'w', encoding='utf-8'))
    _ = load_gene_names(synonyms_path)
    bar.update()
    bar.finish()

    bar = ProgressBar(1, prefix="Loading gene data", file=open(os.devnull, 'w', encoding='utf-8'))
    gene_data = load_json(gene_data_path)
    bar.update()
    bar.finish()

    srna = srna_seq.upper().replace('U', 'T')
    genome_fwd = genome_seq.upper().replace('U', 'T')
    genome_rc = revcomp(genome_fwd)

    n = len(genome_fwd)
    if chunk_size is None or chunk_size <= 0:
        cpu = os.cpu_count() or 4
        base_chunk = max(25000, n // (cpu * 4))
        srna_factor = min(2.0, max(0.5, len(srna) / 100))
        chunk_size = int(base_chunk * srna_factor)
        chunk_size = max(10000, min(100000, chunk_size))
        sys.stderr.write(f"Using adaptive chunk size: {chunk_size}\n")
    step = max(1, chunk_size)

    # Determine windows to scan
    windows: List[int] = []
    if window_sizes is not None:
        for w in window_sizes:
            if isinstance(w, int) and w > 0:
                windows.append(w)
    if not windows:
        windows = [min_len]

    bar = ProgressBar(n * 2 * len(windows), prefix="Preparing chunks", file=open(os.devnull, 'w'))
    jobs: List[Tuple[str, int, int, str, int, int, bool, int]] = []
    for w in windows:
        for start in range(0, n, step):
            end = min(n, start + step)
            jobs.append((genome_fwd, start, end, srna, w, max_gu, False, n))
            jobs.append((genome_rc, start, end, srna, w, max_gu, True, n))
            bar.update(step * 2)
    bar.finish()

    
    seen = set()
    lines: List[str] = []
    total_jobs = len(jobs)

    executor = None
    try:
        executor = ProcessPoolExecutor()
        futures = [executor.submit(scan_chunk, job) for job in jobs]

        completed = 0
        # Use tqdm if available to keep a single progress bar at the bottom
        if tqdm is not None:
            with tqdm(total=total_jobs, desc="Scanning", unit="job", dynamic_ncols=True, leave=True) as pbar:
                for fut in as_completed(futures):
                    res = fut.result()
                    if res:
                        for r in res:
                            key = (
                                tuple(r.get("SmallRNA_region") or []),
                                r.get("Number_of_GU"),
                                r.get("SmallRNA_seq"),
                                tuple(r.get("Target_region") or []),
                                r.get("Target_seq"),
                            )
                            if key in seen:
                                continue
                            seen.add(key)
                            sr_start, sr_end = r["SmallRNA_region"]
                            tg_start, tg_end = r["Target_region"]
                            location = f"{format_location(tg_start, tg_end, gene_data)} [{tg_start},{tg_end}]"
                            lines.append([
                                str(r.get('Window', 0)),
                                str(sr_start),
                                str(sr_end),
                                str(r['Number_of_GU']),
                                r['Strand'],
                                r['SmallRNA_seq'],
                                r['Target_seq'],
                                location
                            ])
                    completed += 1
                    # Debug info inside the progress bar (does not spam the console)
                    if tqdm is not None:
                        pbar.set_postfix({"done": completed, "unique": len(seen)})
                        pbar.update(1)
        else:
            # Fallback minimal progress when tqdm is unavailable
            bar = ProgressBar(total_jobs, prefix="Scanning", file=open(os.devnull, 'w', encoding='utf-8'))
            for fut in as_completed(futures):
                res = fut.result()
                if res:
                    for r in res:
                        key = (
                            tuple(r.get("SmallRNA_region") or []),
                            r.get("Number_of_GU"),
                            r.get("SmallRNA_seq"),
                            tuple(r.get("Target_region") or []),
                            r.get("Target_seq"),
                        )
                        if key in seen:
                            continue
                        seen.add(key)
                        sr_start, sr_end = r["SmallRNA_region"]
                        tg_start, tg_end = r["Target_region"]
                        location = f"{format_location(tg_start, tg_end, gene_data)} [{tg_start},{tg_end}]"
                        lines.append([
                            str(r.get('Window', 0)),
                            str(sr_start),
                            str(sr_end),
                            str(r['Number_of_GU']),
                            r['Strand'],
                            r['SmallRNA_seq'],
                            r['Target_seq'],
                            location
                        ])
                bar.update()
            bar.finish()

    except KeyboardInterrupt:
        if executor is not None:
            try:
                executor.shutdown(wait=False, cancel_futures=True)
            except Exception:
                pass
        sys.stderr.write("\n")
        print_error(
            err_type="Interrupted",
            message="Operation cancelled by user during parallel scan",
            citation="User interruption (Ctrl+C)",
            recommended_next_step="Re-run the command when ready; consider smaller window size to reduce runtime.",
            exit_code=130,
        )

    # Table formatting and printing
    if lines:
        # Render each match as a boxed block for clarity.
        # Box contains a compact metadata line (with # counter) and two sequence lines.
        for idx, row in enumerate(lines, start=1):
            # row indices: 0=Window,1=Start,2=End,3=Wobble,4=Strand,5=sRNA,6=Target,7=Location
            ext_len = len(row[5])
            meta_line = f"# {idx} | Window: {row[0]} | Length: {ext_len} | Start: {row[1]} | End: {row[2]} | Wobble: {row[3]} | Strand: {row[4]} | Location: {row[7]}"
            srna_line = f"sRNA:   {row[5]}"
            target_line = f"Target: {row[6]}"
            box_width = max(len(meta_line), len(srna_line), len(target_line)) + 2
            def box_print(content: str) -> None:
                print("|" + content.ljust(box_width) + "|")
            box_print(meta_line)
            box_print(srna_line)
            box_print(target_line)
            print("|" + ("_" * box_width) + "|")

    # Optionally prompt to save if no output location was provided
    if (save_csv_path or results_dir) or (prompt_save and lines and not (save_csv_path or results_dir)):
        # Default directory per user: Noa_Tools/Results
        base_noa_tools = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
        default_results_dir = os.path.join(base_noa_tools, "Results")
        out_dir = results_dir or default_results_dir
        # If prompting, confirm save and ask for a filename to construct default name with date and user-provided Name
        chosen_path = None
        if prompt_save and not (save_csv_path or results_dir):
            try:
                yn = input("Save results? [Y/N]: ").strip().lower()
            except EOFError:
                yn = "n"
            if yn.startswith('y'):
                # Ask for a human-friendly name to include in the filename
                try:
                    name_input = input("Enter file name (without extension): ").strip()
                except EOFError:
                    name_input = "analysis"
                date_str = datetime.now().strftime("%Y-%m-%d")
                safe_name = "".join(ch for ch in name_input if ch.isalnum() or ch in ("-", "_")) or "analysis"
                filename = f"{date_str}_{safe_name}.csv"
                chosen_path = os.path.join(out_dir, filename)
            else:
                # user chose not to save
                chosen_path = None
                out_dir = None
        if chosen_path is None and not (save_csv_path or results_dir):
            return
        out_path = chosen_path or save_csv_path or os.path.join(out_dir, "srna_genome_matches.csv")
        try:
            os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        except Exception as e:
            print_error(
                err_type="FilesystemError",
                message=f"Could not create directory for: {out_path}. Error: {e}",
                citation="Local filesystem permissions",
                recommended_next_step="Create the directory manually or choose a writable path with --results-dir.",
            )
        try:
            bar = ProgressBar(len(lines), prefix="Saving results")
            with open(out_path, 'w', encoding='utf-8') as f:
                f.write("Window,Start,End,Wobble,Strand,sRNA Sequence,Target Sequence,Location\n")
                for ln in lines:
                    f.write(','.join(ln) + '\n')
                    bar.update()
            bar.finish()
        except Exception as e:
            print_error(
                err_type="FilesystemError",
                message=f"Failed to write CSV: {out_path}. Error: {e}",
                citation="Local filesystem permissions",
                recommended_next_step="Choose a different path or adjust permissions.",
            )