"""I/O utilities and error handling."""
import json
import os
import sys

ESSENTIAL_BASES = set('ACGTU')

def load_text_fasta_or_plain(path: str) -> str:
    with open(path, 'r', encoding='utf-8') as f:
        data = f.read()
    lines = [ln.strip() for ln in data.splitlines() if ln and not ln.startswith('>')]
    seq = ''.join(lines).upper()
    return ''.join(ch for ch in seq if ch in ESSENTIAL_BASES)


def load_gene_names(path: str):
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f]
    except FileNotFoundError:
        return []

def load_json(path: str):
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except FileNotFoundError:
        return {}


def print_error(err_type: str, message: str, citation: str, recommended_next_step: str, exit_code: int = 1) -> None:
    payload = {
        "type": err_type,
        "message": message,
        "citation": citation,
        "recommended_next_step": recommended_next_step,
    }
    sys.stderr.write(json.dumps(payload) + "\n")
    sys.exit(exit_code)


def resolve_paths(genome_arg: str, synonyms_arg: str, gene_data_arg: str):
    base_tools = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    default_genome = os.path.join(base_tools, "data", "B.subtilis_Seq.txt")
    default_synonyms = os.path.join(base_tools, "data", "all_gene_names.txt")
    default_gene_data = os.path.join(base_tools, "data", "all_gene_data.json")
    return (genome_arg or default_genome, synonyms_arg or default_synonyms, gene_data_arg or default_gene_data)


