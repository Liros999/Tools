"""
regulatory_network_analysis.py

Purpose:
- Connects to SubtiWiki API and fetches the regulatory graph.
- Checks which of your genes are in the network.
- For each gene pair, finds the shortest path and direct connection.
- Outputs a summary table and network statistics.
- Provides a simple network visualization.

How to use:
    from regulatory_network_analysis import analyze_gene_connectivity
    my_genes = ["sigA", "spo0A", "abrB"]  # Replace with your gene list
    analyze_gene_connectivity(my_genes)

Requirements:
    pip install requests networkx matplotlib tqdm

Author: [Your Name]
Date: [Today's Date]
"""

import sys
import requests
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm
import time
from gene_synonym_utils import load_synonym_index, resolve_gene_name

BASE_URL = "https://subtiwiki.uni-goettingen.de/v5/api"

def search_gene_id(gene_name):
    """Search for a gene by name and return its SubtiWiki gene ID (or None if not found)."""
    endpoint = f"/search/"
    params = {"q": gene_name, "category": "Gene", "mode": "exact"}
    try:
        response = requests.get(BASE_URL + endpoint, params=params, timeout=20)
        response.raise_for_status()
        data = response.json()
        hits = data.get("data", {}).get("exact_hits", [])
        for hit in hits:
            if hit.get("category") == "Gene" and hit.get("name", "").lower() == gene_name.lower():
                return hit.get("id")
        hits = data.get("data", {}).get("partial_hits", [])
        for hit in hits:
            if hit.get("category") == "Gene" and hit.get("name", "").lower() == gene_name.lower():
                return hit.get("id")
        return None
    except Exception as e:
        print(f"Error searching for gene ID of {gene_name}: {e}")
        return None

def fetch_regulation_graph():
    """Fetch the regulatory graph from SubtiWiki and return genes and regulations."""
    print("Fetching regulatory graph from SubtiWiki", end="", flush=True)
    try:
        for _ in range(6):
            print(".", end="", flush=True)
            time.sleep(0.3)
        response = requests.get(f"{BASE_URL}/regulation/graph", timeout=120)
        print(" done.")
        response.raise_for_status()
        data = response.json().get("data", {})
        genes = data.get("genes", [])
        regulations = data.get("regulations", [])
        if not genes or not regulations:
            print("Error: 'genes' or 'regulations' missing or empty in API response.")
            return None
        return {"genes": genes, "regulations": regulations}
    except requests.exceptions.Timeout:
        print("\nError: The request to SubtiWiki /regulation/graph timed out after 120 seconds.")
        return None
    except Exception as e:
        print(f"\nError fetching regulation graph: {e}")
        return None

class GeneMapper:
    def __init__(self, genes):
        self.id_to_name = {g['id']: g['name'] for g in genes}
        self.name_to_id = {g['name']: g['id'] for g in genes}
    def get_name(self, gene_id):
        return self.id_to_name.get(gene_id, str(gene_id))
    def get_id(self, gene_name):
        if gene_name in self.name_to_id:
            return self.name_to_id[gene_name]
        gene_id = search_gene_id(gene_name)
        if gene_id is not None:
            return gene_id
        return None
    def all_names(self):
        return set(self.name_to_id.keys())

def build_networkx_graph(graph_data):
    """Build a NetworkX DiGraph from SubtiWiki regulation graph data (genes/regulations)."""
    G = nx.DiGraph()
    genes = graph_data["genes"]
    regulations = graph_data["regulations"]
    mapper = GeneMapper(genes)
    for g in genes:
        G.add_node(g['id'], name=g['name'])
    for reg in regulations:
        src = reg['regulator_gene_id']
        tgt = reg['regulated_gene_id']
        if src is None or tgt is None:
            continue  # Skip invalid edges
        G.add_edge(src, tgt, **reg)
    return G, mapper

def prompt_gene_list(mapper):
    while True:
        print("Enter gene names separated by commas (e.g., sigA,spo0A,abrB):")
        genes = [g.strip() for g in input().split(',') if g.strip()]
        if not genes:
            print("No gene names entered. Please try again.")
            continue
        gene_ids = []
        for g in genes:
            gid = mapper.get_id(g)
            if gid is not None:
                gene_ids.append(gid)
            else:
                print(f"Warning: Gene '{g}' not found in SubtiWiki.")
        if gene_ids:
            return gene_ids
        else:
            print("None of the entered genes were found. Please try again.")

def prompt_max_path_length():
    while True:
        try:
            depth = int(input("Enter maximum path length (neighborhood depth, e.g., 2-5): ") or 3)
            if depth > 0:
                return depth
            print("Please enter a positive integer.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def estimate_time(num_genes, depth, pessimistic=False):
    if pessimistic:
        # O(n^3) for all-pairs shortest path
        est = min(5 + 0.05 * num_genes**3, 180)
    else:
        est = min(2 + 0.01 * num_genes**2 * depth, 60)
    return est

def analyze_gene_connectivity(gene_list, interactive=True):
    CACHE_DIR = 'Final_Project/API_Genes/results/subtiwiki_cache'
    name_to_id, synonym_to_id, synonym_to_canonical = load_synonym_index(CACHE_DIR)
    input_ids = []
    for g in gene_list:
        gid, canonical, synonym = resolve_gene_name(g, name_to_id, synonym_to_id, synonym_to_canonical)
        if gid is not None:
            input_ids.append(gid)
            if synonym:
                print(f"Using canonical name '{canonical}' for synonym '{g}'")
        else:
            print(f"Gene/molecule '{g}' not found in synonym index.")
    if not input_ids:
        print("No valid input genes found.")
        return
    graph_data = fetch_regulation_graph()
    if not graph_data:
        print("Failed to fetch or parse regulatory graph. Aborting analysis.")
        return
    G, mapper = build_networkx_graph(graph_data)

    # Map gene names to IDs
    gene_ids = []
    missing_names = []
    for g in input_ids:
        gid = mapper.get_id(g)
        if gid is not None:
            gene_ids.append(gid)
        else:
            print(f"Warning: Gene '{g}' not found in SubtiWiki.")
            missing_names.append(g)

    present_genes = [gid for gid in gene_ids if gid in G.nodes]
    missing_genes = [g for g, gid in zip(input_ids, gene_ids) if gid not in G.nodes]
    missing_genes += missing_names

    if missing_genes:
        print(f"Warning: These gene names are not in the regulatory graph: {missing_genes}")
        if interactive and not present_genes:
            print("No valid genes to analyze. Exiting.")
            return

    print(f"\nGenes to analyze: {[mapper.get_name(gid) for gid in present_genes]}")
    print(f"Missing genes: {missing_genes}")

    if interactive:
        max_path_length = prompt_max_path_length()
    else:
        max_path_length = 3

    est = estimate_time(len(present_genes), max_path_length, pessimistic=True)
    print(f"\nPessimistic estimated time for all-pairs paths: {est:.1f} seconds.")
    if interactive:
        cont = input("Continue? (y/n): ").strip().lower()
        if cont != 'y':
            print("Aborted by user.")
            return

    print("\nGene Connectivity Analysis (all paths up to length {}):".format(max_path_length))
    print(f"{'Gene A':<15}{'Gene B':<15}{'Direct?':<10}{'Num Paths':<10}{'Example Path'}")
    all_paths = []
    direct_edges = set()
    path_counter = Counter()
    edge_labels = {}
    # Collect direct regulatory edges and their labels
    for u, v, data in G.edges(data=True):
        if u in present_genes and v in present_genes:
            direct_edges.add((u, v))
            label = f"{data.get('custom_effect_specifier', data.get('mode', ''))} ({data.get('mechanism', '')})"
            edge_labels[(u, v)] = label
    # Find all shortest paths (direct and indirect)
    for i, gene_a in enumerate(tqdm(present_genes, desc="Gene pairs")):
        for gene_b in present_genes[i+1:]:
            try:
                paths = list(nx.all_shortest_paths(G, gene_a, gene_b)) if nx.has_path(G, gene_a, gene_b) else []
                if not paths:
                    print(f"{mapper.get_name(gene_a):<15}{mapper.get_name(gene_b):<15}{'No':<10}{'0':<10}{'No path found'}")
                else:
                    is_direct = any(len(path) == 2 for path in paths)
                    direct_str = 'Yes' if is_direct else 'No'
                    print(f"{mapper.get_name(gene_a):<15}{mapper.get_name(gene_b):<15}{direct_str:<10}{len(paths):<10}{' -> '.join(mapper.get_name(gid) for gid in paths[0])}")
                for path in paths:
                    all_paths.append(path)
                    for node in path[1:-1]:
                        path_counter[node] += 1
            except Exception as e:
                print(f"{mapper.get_name(gene_a):<15}{mapper.get_name(gene_b):<15}{'ERR':<10}{str(e)}")

    if not all_paths:
        print("No paths found between the selected genes. Only plotting the selected nodes (unconnected).")
        all_paths = [[gid] for gid in present_genes]

    # Summary table of all regulatory relationships
    print("\nSummary Table of Regulatory Relationships:")
    print(f"{'Source':<15}{'Target':<15}{'Mode':<12}{'Effect':<25}{'Mechanism':<15}")
    for (u, v) in direct_edges:
        data = G[u][v]
        mode = data.get('mode', '') or ''
        effect = data.get('custom_effect_specifier', '') or ''
        mechanism = data.get('mechanism', '') or ''
        print(f"{mapper.get_name(u):<15}{mapper.get_name(v):<15}{mode:<12}{effect:<25}{mechanism:<15}")

    print("\nDetecting feedback loops (cycles) involving your genes...")
    cycles = [cycle for cycle in nx.simple_cycles(G) if any(g in cycle for g in present_genes)]
    if cycles:
        print(f"Found {len(cycles)} feedback loops involving your genes.")
        for cycle in cycles[:5]:
            print(f"  Loop: {' -> '.join(mapper.get_name(gid) for gid in cycle)}")
        if len(cycles) > 5:
            print("  ... (more loops not shown)")
    else:
        print("No feedback loops found involving your genes.")

    print("\nGenerating network plot...")
    plot_gene_network(G, present_genes, all_paths, direct_edges, edge_labels, mapper, path_counter)
    print("Done.")

def plot_gene_network(G, genes, paths, direct_edges, edge_labels, mapper, path_counter):
    nodes_in_paths = set()
    edges_in_paths = set()
    for path in paths:
        if path:
            nodes_in_paths.update(path)
            edges_in_paths.update(zip(path[:-1], path[1:]))
    if not nodes_in_paths:
        print("No nodes to plot. The selected genes are not connected and have no reachable subgraphs.")
        return
    subG = G.subgraph(nodes_in_paths)
    pos = nx.spring_layout(subG, seed=42)
    plt.figure(figsize=(10, 8))
    node_sizes = [300 + 700 * (path_counter.get(n, 0) / max(path_counter.values(), default=1)) for n in subG.nodes()]
    nx.draw_networkx_nodes(subG, pos, node_color='skyblue', node_size=node_sizes)
    # Draw direct edges with color and label
    direct_edge_list = [e for e in subG.edges() if e in direct_edges]
    indirect_edge_list = [e for e in subG.edges() if e not in direct_edges]
    direct_colors = []
    for u, v in direct_edge_list:
        mode = (G[u][v].get('mode') or '').lower()
        if mode == 'positive':
            direct_colors.append('green')
        elif mode == 'negative':
            direct_colors.append('red')
        else:
            direct_colors.append('gray')
    nx.draw_networkx_edges(subG, pos, edgelist=direct_edge_list, edge_color=direct_colors, arrows=True, width=2)
    # Draw indirect edges as dashed gray
    nx.draw_networkx_edges(subG, pos, edgelist=indirect_edge_list, edge_color='gray', arrows=True, style='dashed', width=1)
    # Edge labels
    all_edge_labels = {}
    for e in direct_edge_list:
        if e in edge_labels:
            all_edge_labels[e] = edge_labels[e]
    for e in indirect_edge_list:
        all_edge_labels[e] = 'indirect'
    nx.draw_networkx_edge_labels(subG, pos, edge_labels=all_edge_labels, font_size=8)
    labels = {n: mapper.get_name(n) for n in subG.nodes()}
    nx.draw_networkx_labels(subG, pos, labels=labels, font_size=10, font_weight='bold')
    plt.title('Regulatory Network Paths Between Genes')
    plt.axis('off')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    print("Regulatory Network Analysis Tool\n")
    graph_data = fetch_regulation_graph()
    if not graph_data:
        print("Failed to fetch or parse regulatory graph. Aborting analysis.")
        sys.exit(1)
    G, mapper = build_networkx_graph(graph_data)
    if len(sys.argv) > 1:
        genes = sys.argv[1].split(',')
        gene_ids = [mapper.get_id(g) for g in genes if mapper.get_id(g) is not None]
        analyze_gene_connectivity(gene_ids, interactive=False)
    else:
        gene_ids = prompt_gene_list(mapper)
        analyze_gene_connectivity(gene_ids, interactive=True) 