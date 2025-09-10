#!/usr/bin/env python3
"""
Build SubtiWiki metabolite mappings BEFORE app start.

Outputs (written under src/data/cache/):
  - metabolite_categories.json            # hierarchical category list (with children ids)
  - metabolite_to_categories.json         # { metabolite_id: [ {id,name,dot_notation,parent_id} ] }
  - metabolite_to_genes.json              # { metabolite_id: [ gene_name, ... ] }
  - gene_to_metabolites.json              # { gene_name: [ {metabolite_id, metabolite_name, type, transport_type?} ] }

Data sources (SubtiWiki API; see Tools/Subtiwiki_API.txt):
  - GET /interaction/metabolite-protein
  - GET /protein/{id} → gene_id
  - GET /gene/{gene_id}?representation=minimal → name
  - GET /metabolite-category/ and /metabolite-category/{id}

This script uses only real API responses, parallelizes requests, and writes
deterministic JSON caches for the app to load at runtime.
"""

import json
import os
import time
from typing import Dict, Any, List, Tuple, Set, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests


BASE_URL = "https://www.subtiwiki.uni-goettingen.de/v5/api"


def _get_data(payload: Any) -> Any:
    if isinstance(payload, dict) and "data" in payload:
        return payload.get("data")
    return payload


def fetch_json(path: str, *, timeout: int = 60) -> Any:
    url = f"{BASE_URL}{path}"
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return _get_data(r.json())


def fetch_interactions() -> List[Dict[str, Any]]:
    """Fetch all metabolite-protein interactions.

    Returns list of interaction objects with fields including:
      - metabolite {id,name}? or metabolite_name
      - protein {id,name,gene_id}? or protein_id
      - type in [Binding, Cofactor, Effector, Transport]
      - transport_type when type == Transport
    """
    print("Fetching /interaction/metabolite-protein ...")
    data = fetch_json("/interaction/metabolite-protein")
    if not isinstance(data, list):
        return []
    print(f"Fetched {len(data)} metabolite-protein interactions")
    return data


def resolve_protein_to_gene_names_minimal(protein_ids: Set[int], gene_id_to_name: Dict[int, str], *, fallback_limit: int = 200, max_workers: int = 16) -> Dict[int, str]:
    """Minimal resolver:
    - Prefer to avoid expensive lookups.
    - We only call /protein/{id}/gene for at most fallback_limit proteins without inline gene ids.
    - Then map gene_id → gene_name using the pre-fetched /gene/ list.
    """
    protein_to_gene: Dict[int, str] = {}

    def protein_to_geneid(pid: int) -> Tuple[int, Optional[int]]:
        try:
            data = fetch_json(f"/protein/{pid}/gene")
            if isinstance(data, dict) and 'gene_id' in data:
                return pid, data['gene_id']
        except Exception as e:
            print(f"[warn] /protein/{pid}/gene failed: {e}")
        return pid, None

    limited = list(protein_ids)[:max(0, fallback_limit)]
    if not limited:
        return {}

    print(f"Resolving up to {len(limited)} proteins via /protein/{{id}}/gene (cap {fallback_limit})...")
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(protein_to_geneid, pid) for pid in limited]
        for fut in as_completed(futures):
            pid, gid = fut.result()
            if gid is not None:
                name = gene_id_to_name.get(int(gid))
                if name:
                    protein_to_gene[int(pid)] = name

    print(f"Resolved {len(protein_to_gene)} protein ids via minimal fallback")
    return protein_to_gene


def build_gene_metabolite_maps(interactions: List[Dict[str, Any]]) -> Tuple[Dict[str, List[Dict[str, Any]]], Dict[int, Set[str]]]:
    """Build gene→metabolites and metabolite_id→genes using a simple, robust strategy:
    1) Fetch all genes once via /gene/ to map id→name.
    2) For each interaction, prefer protein.gene_id if present; otherwise, defer resolution for a capped set via /protein/{id}/gene.
    3) Produce linear JSON-friendly mappings.
    """
    # 1) Prefetch all genes (minimal list)
    print("Fetching /gene/ for id→name map ...")
    all_genes = fetch_json("/gene/")
    gene_id_to_name: Dict[int, str] = {}
    if isinstance(all_genes, list):
        for g in all_genes:
            gid = g.get('id')
            nm = g.get('name')
            if isinstance(gid, int) and nm:
                gene_id_to_name[gid] = nm
    print(f"Mapped {len(gene_id_to_name)} genes (id→name)")

    # 2) Prepare optional minimal protein→gene fallback set
    protein_ids_needing_fallback: Set[int] = set()

    def extract_gene_name(it: Dict[str, Any]) -> Optional[str]:
        protein = it.get('protein') or {}
        if isinstance(protein, dict):
            gid = protein.get('gene_id') or (protein.get('gene') or {}).get('id')
            if isinstance(gid, int):
                return gene_id_to_name.get(gid)
        pid = None
        if isinstance(protein, dict):
            pid = protein.get('id') or it.get('protein_id')
        else:
            pid = it.get('protein_id')
        if isinstance(pid, int):
            protein_ids_needing_fallback.add(pid)
        return None

    # First pass to assemble entries and collect fallbacks
    tmp_entries: List[Tuple[Optional[str], Dict[str, Any]]] = []
    for it in interactions:
        meta = it.get('metabolite') or {}
        metabolite_id = meta.get('id') if isinstance(meta, dict) else None
        metabolite_name = (meta.get('name') if isinstance(meta, dict) else None) or it.get('metabolite_name')
        if not metabolite_name:
            continue
        gname = extract_gene_name(it)
        tmp_entries.append((gname, {
            'metabolite_id': metabolite_id,
            'metabolite_name': metabolite_name,
            'type': it.get('type'),
            'transport_type': it.get('transport_type') if it.get('type') == 'Transport' else None,
            'protein_id': (it.get('protein') or {}).get('id') or it.get('protein_id')
        }))

    # Minimal fallback resolution (capped)
    fallback_map = resolve_protein_to_gene_names_minimal(protein_ids_needing_fallback, gene_id_to_name)

    # 3) Build mappings
    gene_to_metabolites: Dict[str, List[Dict[str, Any]]] = {}
    metabolite_to_genes: Dict[int, Set[str]] = {}

    for gname, info in tmp_entries:
        if not gname:
            pid = info.get('protein_id')
            if isinstance(pid, int):
                gname = fallback_map.get(pid)
        if not gname:
            continue
        entry = {
            'metabolite_id': info.get('metabolite_id'),
            'metabolite_name': info.get('metabolite_name'),
            'type': info.get('type'),
        }
        if info.get('transport_type'):
            entry['transport_type'] = info.get('transport_type')
        gene_to_metabolites.setdefault(gname, []).append(entry)
        mid = info.get('metabolite_id')
        if isinstance(mid, int):
            metabolite_to_genes.setdefault(mid, set()).add(gname)

    return gene_to_metabolites, metabolite_to_genes


def fetch_metabolite_categories_tree() -> List[Dict[str, Any]]:
    """Fetch metabolite categories and hydrate children via /metabolite-category/{id}."""
    print("Fetching /metabolite-category/ ...")
    cats = fetch_json("/metabolite-category/")
    if not isinstance(cats, list):
        cats = []

    detailed: List[Dict[str, Any]] = []
    for idx, c in enumerate(cats):
        cid = c.get("id")
        try:
            d = fetch_json(f"/metabolite-category/{cid}")
            if not isinstance(d, dict):
                continue
            detailed.append({
                "id": d.get("id"),
                "name": d.get("name"),
                "number": d.get("number"),
                "parent_id": d.get("parent_id"),
                "dot_notation": d.get("dot_notation"),
                "children": [child.get("id") for child in (d.get("children") or []) if isinstance(child, dict)],
            })
            if (idx + 1) % 25 == 0:
                print(f"  hydrated {idx+1}/{len(cats)} categories...")
            time.sleep(0.01)
        except Exception as e:
            print(f"[warn] category {cid} failed: {e}")

    print(f"Hydrated {len(detailed)} metabolite categories")
    return detailed


def build_metabolite_to_categories(detailed_categories: List[Dict[str, Any]]) -> Dict[int, List[Dict[str, Any]]]:
    """Invert category→metabolites using /metabolite-category/{id} details.
    Some deployments include a 'metabolites' list in the category payload.
    We re-fetch per category to collect member metabolite ids/names safely.
    """
    result: Dict[int, List[Dict[str, Any]]] = {}
    for dc in detailed_categories:
        cid = dc.get("id")
        try:
            d = fetch_json(f"/metabolite-category/{cid}")
            metabolites = d.get("metabolites") or []
            for m in metabolites:
                mid = m.get("id")
                if not isinstance(mid, int):
                    continue
                result.setdefault(mid, []).append({
                    "id": cid,
                    "name": d.get("name"),
                    "dot_notation": d.get("dot_notation"),
                    "parent_id": d.get("parent_id"),
                })
            time.sleep(0.02)
        except Exception as e:
            print(f"[warn] invert category {cid} failed: {e}")

    print(f"Built metabolite→categories for {len(result)} metabolites")
    return result


def ensure_cache_dir() -> str:
    cache_dir = os.path.join("src", "data", "cache")
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


def write_json(path: str, obj: Any):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
    print(f"Wrote {path}")


def main():
    print("Building metabolite caches from SubtiWiki API ...")

    # 1) Interactions → gene/metabolite mappings
    interactions = fetch_interactions()
    gene_to_metabolites, metabolite_to_genes_sets = build_gene_metabolite_maps(interactions)

    metabolite_to_genes = {mid: sorted(list(genes)) for mid, genes in metabolite_to_genes_sets.items()}

    # 2) Categories tree and inversion
    metabolite_categories = fetch_metabolite_categories_tree()
    metabolite_to_categories = build_metabolite_to_categories(metabolite_categories)

    # 3) Persist caches
    cache_dir = ensure_cache_dir()
    write_json(os.path.join(cache_dir, "gene_to_metabolites.json"), gene_to_metabolites)
    write_json(os.path.join(cache_dir, "metabolite_to_genes.json"), metabolite_to_genes)
    write_json(os.path.join(cache_dir, "metabolite_categories.json"), metabolite_categories)
    write_json(os.path.join(cache_dir, "metabolite_to_categories.json"), metabolite_to_categories)

    print("All metabolite caches built successfully.")


if __name__ == "__main__":
    main()


