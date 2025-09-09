"""
Interactive/CLI sRNA–Genome complementary finder with configurable G–U wobble pairs.

- Interactive prompts or CLI flags for sRNA, divisions, max-GU, min-len
- Parallel scanning over B. subtilis genome, both strands
- Considers synonyms file availability (loaded for readiness)
- Structured JSON errors to stderr
- CSV output columns:
  SmallRNA:[start,end],NumberOfGU,SmallRNASeq,Target:[start,end],TargetSeq

Traceability (per project requirements):
- See Tools/Subtiwiki_API.txt and Tools/AlphaFold_API.txt for API references and validation.
"""

import argparse
import sys
import os
from typing import Iterable, Tuple, List

from nfgn_lib.core import run_search
from nfgn_lib.io_utils import print_error, resolve_paths

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Find complementary regions between sRNA and B.subtilis genome with G–U wobble support.")
    parser.add_argument("--srna", type=str, help="sRNA DNA sequence (A/C/G/T; U accepted)")
    parser.add_argument("--max-gu", type=int, help="Maximum number of allowed G–U wobble pairs")
    parser.add_argument("--min-len", type=int, default=None, help="Seed/window size for exact seeding (default 14 if omitted)")
    parser.add_argument("--win-range", type=str, default=None, help="Window size range or list, e.g. '8-14' or '8,10,12'")
    parser.add_argument("--genome", type=str, default=None, help="Path to B.subtilis genome file (FASTA/plain)")
    parser.add_argument("--synonyms", type=str, default=None, help="Path to gene_synonyms.json")
    parser.add_argument("--gene-data", type=str, default=None, help="Path to all_gene_data.json")
    parser.add_argument("--interactive", action="store_true", help="Prompt for any missing parameters interactively")
    parser.add_argument("--save-csv", type=str, default=None, help="Optional: exact CSV output path")
    parser.add_argument("--results-dir", type=str, default=None, help="Optional: directory to save default CSV (results folder)")
    parser.add_argument("--legacy", action="store_true", help="Use legacy non-optimized scheduling (default is optimized)")
    return parser.parse_args(list(argv))


def parse_window_sizes(min_len: int, win_range: str) -> List[int]:
    windows: List[int] = []
    if win_range:
        txt = win_range.strip()
        if '-' in txt:
            a, b = txt.split('-', 1)
            try:
                start = int(a)
                end = int(b)
                if start <= end:
                    windows = list(range(start, end + 1))
                else:
                    windows = list(range(end, start + 1))
            except ValueError:
                windows = []
        else:
            parts = [p.strip() for p in txt.split(',') if p.strip()]
            try:
                windows = [int(p) for p in parts if int(p) > 0]
            except ValueError:
                windows = []
    if not windows and min_len:
        windows = [min_len]
    if not windows:
        windows = [14]
    return sorted(set(windows))


def main() -> None:
    args = parse_args(sys.argv[1:])

    srna = args.srna
    max_gu = args.max_gu
    min_len = args.min_len
    win_range = args.win_range
    save_csv = args.save_csv
    results_dir = args.results_dir
    gene_data_path = args.gene_data

    if not srna:
        srna = input("Enter sRNA DNA sequence (A/C/G/T; U allowed): ").strip()
    if max_gu is None:
        try:
            max_gu = int(input("Enter maximum number of G–U wobble pairs (>=0): ").strip())
        except Exception:
            print_error(
                err_type="InvalidInput",
                message="max_gu must be an integer",
                citation="Interactive input validation",
                recommended_next_step="Re-run and provide a non-negative integer for max_gu.",
            )
    # Upfront prompt for window size/range when neither is provided
    if min_len is None and not win_range:
        wtxt = input("Enter window size or range/list (e.g., 14 or 8-14 or 8,10,12) [default 14]: ").strip()
        if wtxt:
            if any(ch in wtxt for ch in ('-', ',')):
                win_range = wtxt
                min_len = None
            else:
                try:
                    min_len = int(wtxt)
                except Exception:
                    min_len = 14
        else:
            min_len = 14

    if srna is None or max_gu is None or (min_len is None and not win_range):
        print_error(
            err_type="MissingParameters",
            message="Missing required parameters. Provide --srna, --max-gu, and window size (--min-len or --win-range) or follow prompts.",
            citation="CLI usage",
            recommended_next_step="Pass all flags or respond to prompts.",
        )

    genome_path, synonyms_path, gene_data_path = resolve_paths(args.genome, args.synonyms, args.gene_data)

    while True:
        windows = parse_window_sizes(min_len, win_range)
        run_search(
            srna_seq=srna,
            chunk_size=None,
            max_gu=max_gu,
            min_len=min_len or (windows[0] if windows else 14),
            genome_path=genome_path,
            synonyms_path=synonyms_path,
            gene_data_path=gene_data_path,
            save_csv_path=save_csv,
            results_dir=results_dir,
            window_sizes=windows,
            prompt_save=True,
        )

        # Post-run menu: ask if run another analysis
        again = input("Run another analysis? [Y/N]: ").strip().lower()
        if not again.startswith('y'):
            break

        # Ask whether to modify current parameters or start a new one
        reuse = input("Reuse current parameters? [Y] reuse / [M] modify / [N] new: ").strip().lower()
        if reuse.startswith('y'):
            continue
        if reuse.startswith('n'):
            # Start new: prompt essentials again
            srna = input("sRNA sequence: ").strip()
            try:
                max_gu = int(input("max-GU (>=0): ").strip())
            except Exception:
                pass
            wtxt = input("window (single) or range/list (e.g., 14 or 8-14 or 8,10,12): ").strip()
            if wtxt:
                if any(ch in wtxt for ch in ('-', ',')):
                    win_range = wtxt
                    min_len = None
                else:
                    try:
                        min_len = int(wtxt)
                        win_range = None
                    except Exception:
                        pass
            gp = input("genome path (leave empty to keep current): ").strip()
            if gp:
                genome_path = gp
            continue

        # Modify current
        if True:
            # Simple menu to adjust high-level choices
            print("Modify: 1) sRNA  2) max-GU  3) window single  4) window range/list  5) genome path  6) save csv path  7) results dir  8) synonyms  9) gene data")
            try:
                sel = int(input("Select item to modify (1-9): ").strip())
            except Exception:
                sel = 0
            if sel == 1:
                srna = input("sRNA sequence: ").strip()
            elif sel == 2:
                try:
                    max_gu = int(input("max-GU (>=0): ").strip())
                except Exception:
                    pass
            elif sel == 3:
                try:
                    min_len = int(input("single window size: ").strip())
                    win_range = None
                except Exception:
                    pass
            elif sel == 4:
                win_range = input("window range/list (e.g., 8-14 or 8,10,12): ").strip()
            elif sel == 5:
                genome_path = input("genome path: ").strip()
            elif sel == 6:
                p = input("save CSV path (empty to disable): ").strip()
                save_csv = p or None
            elif sel == 7:
                p = input("results dir (empty to default): ").strip()
                results_dir = p or None
            elif sel == 8:
                synonyms_path = input("synonyms path: ").strip()
            elif sel == 9:
                gene_data_path = input("gene data (all_gene_data.json) path: ").strip()
            continue
        # Fallback full modify flow
            # Quick modify all in one flow
            srna = input(f"sRNA sequence [{srna}]: ").strip() or srna
            try:
                mg = input(f"max-GU [{max_gu}]: ").strip()
                if mg:
                    max_gu = int(mg)
            except Exception:
                pass
            wtxt = input(f"window (single) or range/list [{win_range or min_len}]: ").strip()
            if wtxt:
                if any(ch in wtxt for ch in ('-', ',')):
                    win_range = wtxt
                    min_len = None
                else:
                    try:
                        min_len = int(wtxt)
                        win_range = None
                    except Exception:
                        pass
            gp = input(f"genome path [{genome_path}]: ").strip()
            if gp:
                genome_path = gp
            sc = input(f"save CSV path [{save_csv or ''}]: ").strip()
            if sc:
                save_csv = sc
            rd = input(f"results dir [{results_dir or ''}]: ").strip()
            if rd:
                results_dir = rd
            sp = input(f"synonyms path [{synonyms_path}]: ").strip()
            if sp:
                synonyms_path = sp
            gd = input(f"gene data path [{gene_data_path}]: ").strip()
            if gd:
                gene_data_path = gd
            continue
        # default: loop again
        continue


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print_error(
            err_type="Interrupted",
            message="Operation cancelled by user",
            citation="User interruption",
            recommended_next_step="Re-run the command when ready.",
            exit_code=130,
        )