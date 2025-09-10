import json
import os
from typing import Dict, Any, List


_CACHE_DIR = os.path.join('src', 'data', 'cache')


def _read_json_safe(path: str):
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return None


def get_gene_to_metabolites() -> Dict[str, List[Dict[str, Any]]]:
    path = os.path.join(_CACHE_DIR, 'gene_to_metabolites.json')
    data = _read_json_safe(path)
    return data or {}


def get_metabolite_to_genes() -> Dict[str, List[str]]:
    path = os.path.join(_CACHE_DIR, 'metabolite_to_genes.json')
    data = _read_json_safe(path)
    return data or {}


def get_metabolite_categories() -> List[Dict[str, Any]]:
    path = os.path.join(_CACHE_DIR, 'metabolite_categories.json')
    data = _read_json_safe(path)
    return data or []


def get_metabolite_to_categories() -> Dict[str, List[Dict[str, Any]]]:
    path = os.path.join(_CACHE_DIR, 'metabolite_to_categories.json')
    data = _read_json_safe(path)
    return data or {}


