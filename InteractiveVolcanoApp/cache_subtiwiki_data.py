import json
import time
from typing import Dict, Any, List, Set
import requests

BASE_URL = "https://www.subtiwiki.uni-goettingen.de/v5/api"

def fetch_json(path: str, *, retries: int = 3, timeout: int = 30) -> Any:
    url = f"{BASE_URL}{path}"
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            data = r.json()
            return data.get('data', data)
        except Exception as e:
            last_err = e
            print(f"[warn] fetch {url} failed (attempt {attempt}/{retries}): {e}")
            time.sleep(2 * attempt)
    print(f"[error] giving up on {url}: {last_err}")
    return None

def build_protein_complex_map() -> Dict[str, List[str]]:
    """Fetch interaction graph and derive a protein complex map.

    Strategy:
    - Fetch /interaction/graph once (contains molecules and interactions)
    - Build protein_id -> gene_name map from molecule entries
      (use gene_name or fallback to /protein/{id} then /gene/{gene_id}?representation=minimal)
    - For each interaction that involves >=2 proteins, treat as a complex edge set and
      aggregate member genes.
    - Return { "Complex: A & B": ["A","B"], ... }
    """
    print("Fetching interaction graph...")
    graph = fetch_json('/interaction/graph')
    if not graph or not isinstance(graph, dict):
        print("[warn] No interaction graph available; skipping complex map build")
        return {}

    molecules = graph.get('molecules', []) or []
    interactions = graph.get('interactions', []) or []

    print(f"Graph summary: molecules={len(molecules)} interactions={len(interactions)}")

    protein_to_gene: Dict[str, str] = {}
    # First pass: try to resolve from molecule records
    for m in molecules:
        if not isinstance(m, dict):
            continue
        if m.get('type') != 'Protein':
            continue
        pid = m.get('id')
        gname = m.get('gene_name') or m.get('name') or (m.get('gene') or {}).get('name')
        if pid is not None and gname:
            protein_to_gene[str(pid)] = gname

    # Second pass: resolve missing proteins by fetching /protein/{id} then /gene/{gene_id}
    missing = [m.get('id') for m in molecules if isinstance(m, dict) and m.get('type') == 'Protein' and str(m.get('id')) not in protein_to_gene]
    if missing:
        print(f"Resolving {len(missing)} protein->gene gaps via /protein/{{id}}")
    for pid in missing[:2000]:  # safety cap
        pdata = fetch_json(f"/protein/{pid}")
        if not pdata or not isinstance(pdata, dict):
            continue
        gene_id = pdata.get('gene_id') or (pdata.get('gene') or {}).get('id')
        gname = None
        if gene_id:
            gdata = fetch_json(f"/gene/{gene_id}?representation=minimal")
            if isinstance(gdata, dict):
                gname = gdata.get('name')
        if not gname:
            gname = pdata.get('name')  # last resort
        if gname:
            protein_to_gene[str(pid)] = gname
        time.sleep(0.05)  # be kind

    # Build complexes by aggregating interactions with >=2 protein members
    complexes: Dict[str, Set[str]] = {}
    for inter in interactions:
        mids = inter.get('molecule_ids', []) if isinstance(inter, dict) else []
        genes = sorted({protein_to_gene.get(str(mid)) for mid in mids if protein_to_gene.get(str(mid))})
        if len(genes) >= 2:
            key = ' & '.join(genes)
            complexes.setdefault(key, set()).update(genes)

    # Format output
    complex_map: Dict[str, List[str]] = {f"Complex: {k}": sorted(list(v)) for k, v in complexes.items()}
    print(f"Built {len(complex_map)} complexes from interactions")
    return complex_map

def fetch_and_cache_data():
    """
    Fetches all genes, categories, and gene-category relationships from the SubtiWiki API
    and caches the data in local JSON files.
    """
    # Fetch all genes
    print("Fetching all genes...")
    all_genes_url = f"{BASE_URL}/gene/"
    response = requests.get(all_genes_url)
    if response.status_code == 200:
        all_genes = response.json()["data"]
        with open("static/api_cache/all_genes.json", "w", encoding="utf-8") as f:
            json.dump(all_genes, f, indent=2, ensure_ascii=False)
        print("Successfully fetched and cached all genes.")
    else:
        print(f"Failed to fetch all genes. Status code: {response.status_code}")

    # Fetch all categories
    print("Fetching all categories...")
    all_categories_url = f"{BASE_URL}/gene-category/"
    response = requests.get(all_categories_url)
    if response.status_code == 200:
        all_categories = response.json()["data"]
        with open("static/api_cache/all_categories.json", "w", encoding="utf-8") as f:
            json.dump(all_categories, f, indent=2, ensure_ascii=False)
        print("Successfully fetched and cached all categories.")

        # Fetch detailed category information including hierarchical structure
        print("Fetching detailed category information and genes...")
        gene_category_map = {}
        detailed_categories = []
        
        for category in all_categories:
            category_id = category["id"]
            category_name = category["name"]
            print(f"  Processing category {category_id}: {category_name}")
            
            genes_in_category_url = f"{BASE_URL}/gene-category/{category_id}"
            response = requests.get(genes_in_category_url)
            if response.status_code == 200:
                detailed_data = response.json()["data"]
                
                # Extract genes for this category
                genes = [gene["name"] for gene in detailed_data.get("genes", [])]
                gene_category_map[category_name] = genes
                
                # Create detailed category info with hierarchical data
                detailed_category = {
                    "id": category_id,
                    "name": category_name,
                    "type": detailed_data.get("type", "Gene"),
                    "number": detailed_data.get("number"),
                    "parent_id": detailed_data.get("parent_id"),
                    "dot_notation": detailed_data.get("dot_notation"),
                    "gene_count": len(genes),
                    "children_ids": [child["id"] for child in detailed_data.get("children", [])]
                }
                detailed_categories.append(detailed_category)
            else:
                print(f"    Failed to fetch details for category {category_name}. Status code: {response.status_code}")
        
        # Save the detailed categories with hierarchical information
        with open("static/api_cache/detailed_categories.json", "w", encoding="utf-8") as f:
            json.dump(detailed_categories, f, indent=2, ensure_ascii=False)
        print("Successfully fetched and cached detailed category hierarchy.")
        
        with open("static/api_cache/gene_category_map.json", "w", encoding="utf-8") as f:
            json.dump(gene_category_map, f, indent=2, ensure_ascii=False)
        print("Successfully fetched and cached gene-category relationships.")

    else:
        print(f"Failed to fetch all categories. Status code: {response.status_code}")

    # Build and cache protein complex map (optional, but recommended)
    try:
        complex_map = build_protein_complex_map()
        with open("static/api_cache/protein_complex_map.json", "w", encoding="utf-8") as f:
            json.dump(complex_map, f, indent=2, ensure_ascii=False)
        print(f"Successfully built and cached protein_complex_map.json with {len(complex_map)} complexes.")
    except Exception as e:
        print(f"[warn] Failed to build protein_complex_map.json: {e}")

if __name__ == "__main__":
    fetch_and_cache_data()
