#!/usr/bin/env python3
"""
Script to fetch interaction data from SubtiWiki API and build adjacency map
"""

import requests
import json
import os
import sys
from collections import defaultdict

def fetch_interaction_graph():
    """Fetch the complete interaction graph from SubtiWiki"""
    print("Fetching interaction graph from SubtiWiki...")
    
    url = "https://www.subtiwiki.uni-goettingen.de/v5/api/interaction/graph"
    
    try:
        response = requests.get(url, timeout=60)
        if response.status_code == 200:
            data = response.json()
            if data.get('isSuccess'):
                return data.get('data', {})
            else:
                print(f"API returned success=false: {data.get('message')}")
                return None
        else:
            print(f"HTTP error: {response.status_code}")
            return None
    except Exception as e:
        print(f"Request failed: {e}")
        return None

def save_raw_data(interaction_data):
    """Save the raw interaction data"""
    cache_dir = os.path.join('src', 'data', 'cache')
    os.makedirs(cache_dir, exist_ok=True)
    
    # Save raw interaction data
    raw_file = os.path.join(cache_dir, 'raw_interaction_graph.json')
    with open(raw_file, 'w', encoding='utf-8') as f:
        json.dump(interaction_data, f, indent=2, ensure_ascii=False)
    
    print(f"Saved raw interaction data to {raw_file}")
    
    molecules = interaction_data.get('molecules', [])
    interactions = interaction_data.get('interactions', [])
    
    print(f"Total molecules: {len(molecules)}")
    print(f"Total interactions: {len(interactions)}")
    
    # Show some sample molecules
    proteins = [m for m in molecules if m.get('type') == 'Protein']
    print(f"Proteins: {len(proteins)}")
    
    for mol in proteins[:10]:
        print(f"  {mol.get('name')} (ID: {mol.get('id')})")
    
    # Show some sample interactions
    print(f"\nSample interactions:")
    for interaction in interactions[:5]:
        print(f"  {interaction}")

def build_simple_adjacency(interaction_data):
    """Build a simple adjacency map using molecule names as gene names"""
    print("\nBuilding simple adjacency map...")
    
    molecules = interaction_data.get('molecules', [])
    interactions = interaction_data.get('interactions', [])
    
    # Build molecule ID to name map
    molecule_names = {mol.get('id'): mol.get('name') for mol in molecules}
    
    # Build adjacency map
    adjacency = defaultdict(set)
    
    for interaction in interactions:
        molecule_ids = interaction.get('molecule_ids', [])
        
        # Only process interactions with 2+ molecules
        if len(molecule_ids) >= 2:
            # Get names for all molecules in this interaction
            names_in_interaction = []
            for mol_id in molecule_ids:
                if mol_id in molecule_names:
                    names_in_interaction.append(molecule_names[mol_id])
            
            # Add edges between all molecules in this interaction
            for i, name1 in enumerate(names_in_interaction):
                for name2 in names_in_interaction[i+1:]:
                    if name1 != name2:
                        adjacency[name1].add(name2)
                        adjacency[name2].add(name1)
    
    print(f"Built adjacency map with {len(adjacency)} molecules")
    
    # Save adjacency map
    cache_dir = os.path.join('src', 'data', 'cache')
    adjacency_file = os.path.join(cache_dir, 'protein_interactions.json')
    
    # Convert sets to lists for JSON serialization
    adjacency_serializable = {name: list(neighbors) for name, neighbors in adjacency.items()}
    
    with open(adjacency_file, 'w', encoding='utf-8') as f:
        json.dump(adjacency_serializable, f, indent=2, ensure_ascii=False)
    
    print(f"Saved adjacency map to {adjacency_file}")
    
    # Show some sample data
    sample_names = list(adjacency.keys())[:10]
    for name in sample_names:
        print(f"  {name}: {len(adjacency[name])} neighbors")
        if adjacency[name]:
            print(f"    Neighbors: {list(adjacency[name])[:5]}...")
    
    return adjacency

def main():
    print("Building interaction cache from SubtiWiki API...")
    
    # Fetch interaction data
    interaction_data = fetch_interaction_graph()
    if not interaction_data:
        print("Failed to fetch interaction data")
        return
    
    # Save raw data
    save_raw_data(interaction_data)
    
    # Build simple adjacency map
    adjacency_map = build_simple_adjacency(interaction_data)
    
    print("\nCache building complete!")

if __name__ == "__main__":
    main()
